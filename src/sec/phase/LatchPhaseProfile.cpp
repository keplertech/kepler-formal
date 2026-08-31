// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "phase/LatchPhaseProfile.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numeric>

namespace KEPLER_FORMAL::SEC::PHASE {

namespace {

struct SpecializationScratch {
  std::vector<ExprID> memo;
  std::vector<uint32_t> epochs;
  uint32_t epoch = 0;

  void begin(size_t nodeCount) {
    if (memo.size() < nodeCount) {
      memo.resize(nodeCount, InvalidExpr);
      epochs.resize(nodeCount, 0);
    }
    ++epoch;
    if (epoch == 0) {
      std::fill(epochs.begin(), epochs.end(), 0);
      epoch = 1;
    }
  }
};

size_t effectiveConcurrency(size_t requested) {
  return std::getenv("KEPLER_NO_MT") != nullptr ? 1 : requested;
}

template <typename Function>
void parallelIndexes(size_t count, size_t concurrency, Function&& function) {
  if (count == 0) {
    return;
  }
  const size_t effective = effectiveConcurrency(concurrency);
  if (effective == 1 || count == 1) {
    function(0, count);
    return;
  }
  const int arenaConcurrency =
      effective == 0
          ? tbb::task_arena::automatic
          : static_cast<int>(std::min<size_t>(
                effective,
                static_cast<size_t>(std::numeric_limits<int>::max())));
  tbb::task_arena arena(arenaConcurrency);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, count),
                      [&](const tbb::blocked_range<size_t>& range) {
                        function(range.begin(), range.end());
                      });
  });
}

ExprID specializeExpression(const NormalizedTransitionSystem& system,
                            ExprID root,
                            const std::vector<TernaryValue>& assignments,
                            ExpressionArena& destination,
                            SpecializationScratch& scratch) {
  const auto& source = system.expressions();
  scratch.begin(source.nodes().size());
  std::function<ExprID(ExprID)> copy = [&](ExprID id) -> ExprID {
    if (!source.isValid(id)) {
      return InvalidExpr;
    }
    if (scratch.epochs[id] == scratch.epoch) {
      return scratch.memo[id];
    }
    const auto& node = source.node(id);
    ExprID result = InvalidExpr;
    switch (node.op) {
      case ExpressionOp::Constant:
        result = destination.constant(node.constant);
        break;
      case ExpressionOp::Variable:
        if (node.variable >= assignments.size()) {
          return InvalidExpr;
        }
        if (assignments[node.variable] == TernaryValue::Zero) {
          result = destination.constant(false);
        } else if (assignments[node.variable] == TernaryValue::One) {
          result = destination.constant(true);
        } else {
          result = destination.variable(node.variable);
        }
        break;
      case ExpressionOp::Not:
        result = destination.makeNot(copy(node.left));
        break;
      case ExpressionOp::And:
        result = destination.makeAnd(copy(node.left), copy(node.right));
        break;
      case ExpressionOp::Or:
        result = destination.makeOr(copy(node.left), copy(node.right));
        break;
      case ExpressionOp::Xor:
        result = destination.makeXor(copy(node.left), copy(node.right));
        break;
    }
    scratch.memo[id] = result;
    scratch.epochs[id] = scratch.epoch;
    return result;
  };
  return copy(root);
}

}  // namespace

const char* toString(TransparencyKind kind) {
  switch (kind) {
    case TransparencyKind::Closed:
      return "closed";
    case TransparencyKind::Transparent:
      return "transparent";
    case TransparencyKind::Conditional:
      return "conditional";
  }
  return "conditional";
}

LatchPhaseProfileResult buildLatchPhaseProfiles(
    const NormalizedTransitionSystem& system,
    const PeriodicSignalAnalysisResult& periodic,
    const LatchPhaseProfileOptions& options) {
  LatchPhaseProfileResult result;
  if (!periodic.complete || periodic.phaseCount == 0) {
    result.complete = false;
    result.diagnostic =
        "latch profiling requires a complete periodic-signal analysis";
    return result;
  }

  std::vector<std::vector<TernaryValue>> phaseAssignments(
      periodic.phaseCount, std::vector<TernaryValue>(system.variables().size(),
                                                     TernaryValue::Unknown));
  for (const auto& generator : periodic.generators) {
    if (!generator.usableAtSelectedPhaseCount) {
      continue;
    }
    if (generator.stateIndex >= system.states().size() ||
        generator.expandedWord.size() != periodic.phaseCount) {
      result.complete = false;
      result.diagnostic = "periodic generator is inconsistent with the model";
      return result;
    }
    const size_t variable = system.states()[generator.stateIndex].variable;
    for (size_t phase = 0; phase < periodic.phaseCount; ++phase) {
      auto& assignment = phaseAssignments[phase][variable];
      const auto value = generator.expandedWord[phase];
      if (assignment != TernaryValue::Unknown && assignment != value) {
        result.complete = false;
        result.diagnostic = "periodic generators assign conflicting values";
        return result;
      }
      assignment = value;
    }
  }

  if (!system.latches().empty() &&
      periodic.phaseCount >
          std::numeric_limits<size_t>::max() / system.latches().size()) {
    result.complete = false;
    result.diagnostic = "latch phase-profile size overflows size_t";
    return result;
  }
  const size_t slotCount = system.latches().size() * periodic.phaseCount;
  std::vector<TernaryValue> evaluated(slotCount, TernaryValue::Unknown);
  parallelIndexes(slotCount, options.maxConcurrency,
                  [&](size_t begin, size_t end) {
                    for (size_t slot = begin; slot < end; ++slot) {
                      const size_t latchIndex = slot / periodic.phaseCount;
                      const size_t phase = slot % periodic.phaseCount;
                      evaluated[slot] = evaluateTernaryExpression(
                          system, system.latches()[latchIndex].transparency,
                          phaseAssignments[phase]);
                    }
                  });

  std::vector<size_t> stableLatchOrder(system.latches().size());
  std::iota(stableLatchOrder.begin(), stableLatchOrder.end(), 0);
  std::sort(stableLatchOrder.begin(), stableLatchOrder.end(),
            [&](size_t lhs, size_t rhs) {
              const auto& left = system.latches()[lhs];
              const auto& right = system.latches()[rhs];
              if (left.instancePath != right.instancePath) {
                return left.instancePath < right.instancePath;
              }
              if (left.modelName != right.modelName) {
                return left.modelName < right.modelName;
              }
              return left.stateIndex < right.stateIndex;
            });

  result.profiles.reserve(stableLatchOrder.size());
  SpecializationScratch specializationScratch;
  for (const size_t latchIndex : stableLatchOrder) {
    LatchPhaseProfile profile;
    profile.latchIndex = latchIndex;
    profile.phases.reserve(periodic.phaseCount);
    for (size_t phase = 0; phase < periodic.phaseCount; ++phase) {
      const auto value = evaluated[latchIndex * periodic.phaseCount + phase];
      if (value == TernaryValue::Zero) {
        profile.phases.push_back({TransparencyKind::Closed,
                                  result.residualExpressions.constant(false)});
      } else if (value == TernaryValue::One) {
        profile.phases.push_back({TransparencyKind::Transparent,
                                  result.residualExpressions.constant(true)});
      } else {
        const ExprID residual = specializeExpression(
            system, system.latches()[latchIndex].transparency,
            phaseAssignments[phase], result.residualExpressions,
            specializationScratch);
        if (residual == InvalidExpr) {
          result.complete = false;
          result.diagnostic = "failed to specialize latch transparency";
          return result;
        }
        profile.phases.push_back({TransparencyKind::Conditional, residual});
      }
    }
    result.profiles.push_back(std::move(profile));
  }
  return result;
}

}  // namespace KEPLER_FORMAL::SEC::PHASE
