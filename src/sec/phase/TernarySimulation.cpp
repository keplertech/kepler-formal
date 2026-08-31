// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "phase/TernarySimulation.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/info.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task_arena.h>

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <unordered_map>

namespace KEPLER_FORMAL::SEC::PHASE {

namespace {

struct EvaluationScratch {
  std::vector<TernaryValue> values;
  std::vector<uint32_t> epochs;
  uint32_t epoch = 0;
  struct Frame {
    ExprID expression = InvalidExpr;
    bool visited = false;
  };
  std::vector<Frame> stack;

  void begin(size_t nodeCount) {
    if (values.size() < nodeCount) {
      values.resize(nodeCount, TernaryValue::Unknown);
      epochs.resize(nodeCount, 0);
    }
    ++epoch;
    if (epoch == 0) {
      std::fill(epochs.begin(), epochs.end(), 0);
      epoch = 1;
    }
    stack.clear();
  }
};

TernaryValue ternaryNot(TernaryValue value) {
  if (value == TernaryValue::Zero) {
    return TernaryValue::One;
  }
  if (value == TernaryValue::One) {
    return TernaryValue::Zero;
  }
  return TernaryValue::Unknown;
}

TernaryValue ternaryAnd(TernaryValue lhs, TernaryValue rhs) {
  if (lhs == TernaryValue::Zero || rhs == TernaryValue::Zero) {
    return TernaryValue::Zero;
  }
  if (lhs == TernaryValue::One && rhs == TernaryValue::One) {
    return TernaryValue::One;
  }
  return TernaryValue::Unknown;
}

TernaryValue ternaryOr(TernaryValue lhs, TernaryValue rhs) {
  if (lhs == TernaryValue::One || rhs == TernaryValue::One) {
    return TernaryValue::One;
  }
  if (lhs == TernaryValue::Zero && rhs == TernaryValue::Zero) {
    return TernaryValue::Zero;
  }
  return TernaryValue::Unknown;
}

TernaryValue ternaryXor(TernaryValue lhs, TernaryValue rhs) {
  if (lhs == TernaryValue::Unknown || rhs == TernaryValue::Unknown) {
    return TernaryValue::Unknown;
  }
  return lhs == rhs ? TernaryValue::Zero : TernaryValue::One;
}

bool isValidTernaryValue(TernaryValue value) {
  switch (value) {
    case TernaryValue::Zero:
    case TernaryValue::One:
    case TernaryValue::Unknown:
      return true;
  }
  return false;
}

TernaryValue evaluateWithScratch(
    const NormalizedTransitionSystem& system, ExprID expression,
    const std::vector<TernaryValue>& variableValues,
    EvaluationScratch& scratch) {
  const auto& arena = system.expressions();
  if (!arena.isValid(expression)) {
    throw std::invalid_argument("invalid ternary expression root");
  }
  if (variableValues.size() != system.variables().size()) {
    throw std::invalid_argument(
        "ternary variable environment has the wrong width");
  }
  scratch.begin(arena.nodes().size());
  scratch.stack.push_back({expression, false});
  while (!scratch.stack.empty()) {
    const auto frame = scratch.stack.back();
    scratch.stack.pop_back();
    if (scratch.epochs[frame.expression] == scratch.epoch) {
      continue;
    }
    const auto& node = arena.node(frame.expression);
    if (!frame.visited) {
      switch (node.op) {
        case ExpressionOp::Constant:
          scratch.values[frame.expression] =
              node.constant ? TernaryValue::One : TernaryValue::Zero;
          scratch.epochs[frame.expression] = scratch.epoch;
          continue;
        case ExpressionOp::Variable:
          if (node.variable >= variableValues.size()) {
            throw std::invalid_argument(
                "ternary expression references an invalid variable");
          }
          if (!isValidTernaryValue(variableValues[node.variable])) {
            throw std::invalid_argument(
                "ternary expression reads an invalid variable value");
          }
          scratch.values[frame.expression] = variableValues[node.variable];
          scratch.epochs[frame.expression] = scratch.epoch;
          continue;
        case ExpressionOp::Not:
          if (!arena.isValid(node.left) || node.left >= frame.expression) {
            throw std::invalid_argument("malformed ternary NOT expression");
          }
          scratch.stack.push_back({frame.expression, true});
          scratch.stack.push_back({node.left, false});
          break;
        case ExpressionOp::And:
        case ExpressionOp::Or:
        case ExpressionOp::Xor:
          if (!arena.isValid(node.left) || !arena.isValid(node.right) ||
              node.left >= frame.expression || node.right >= frame.expression) {
            throw std::invalid_argument("malformed ternary binary expression");
          }
          scratch.stack.push_back({frame.expression, true});
          scratch.stack.push_back({node.right, false});
          scratch.stack.push_back({node.left, false});
          break;
        default:
          throw std::invalid_argument(
              "unsupported ternary expression operator");
      }
      continue;
    }

    const auto lhs = scratch.values[node.left];
    switch (node.op) {
      case ExpressionOp::Not:
        scratch.values[frame.expression] = ternaryNot(lhs);
        break;
      case ExpressionOp::And:
        scratch.values[frame.expression] =
            ternaryAnd(lhs, scratch.values[node.right]);
        break;
      case ExpressionOp::Or:
        scratch.values[frame.expression] =
            ternaryOr(lhs, scratch.values[node.right]);
        break;
      case ExpressionOp::Xor:
        scratch.values[frame.expression] =
            ternaryXor(lhs, scratch.values[node.right]);
        break;
      case ExpressionOp::Constant:
      case ExpressionOp::Variable:
        throw std::invalid_argument("malformed ternary expression DAG");
      default:
        throw std::invalid_argument("unsupported ternary expression operator");
    }
    scratch.epochs[frame.expression] = scratch.epoch;
  }
  return scratch.values[expression];
}

struct TernaryStateHash {
  size_t operator()(const TernaryState& state) const noexcept {
    size_t hash = 1469598103934665603ULL;
    for (const auto value : state) {
      hash ^= static_cast<size_t>(value) + 1;
      hash *= 1099511628211ULL;
    }
    return hash;
  }
};

size_t effectiveConcurrency(size_t requested) {
  if (std::getenv("KEPLER_NO_MT") != nullptr || requested == 1) {
    return 1;
  }
  const size_t available =
      static_cast<size_t>(std::max(1, tbb::info::default_concurrency()));
  return requested == 0 ? available : std::min(requested, available);
}

bool canStoreState(size_t currentCount, size_t stateWidth, size_t byteLimit) {
  if (byteLimit == 0 || stateWidth == 0) {
    return true;
  }
  if (stateWidth > byteLimit / sizeof(TernaryValue)) {
    return false;
  }
  const size_t bytesPerState = stateWidth * sizeof(TernaryValue);
  return currentCount < byteLimit / bytesPerState;
}

}  // namespace

TernaryValue evaluateTernaryExpression(
    const NormalizedTransitionSystem& system, ExprID expression,
    const std::vector<TernaryValue>& variableValues) {
  // Latch profiling evaluates many independent roots on the same workers.
  // Reusing scratch per worker avoids allocating an arena-sized memo table for
  // every latch/phase slot while remaining reentrant across threads.
  thread_local EvaluationScratch scratch;
  return evaluateWithScratch(system, expression, variableValues, scratch);
}

const char* toString(TernarySimulationStatus status) {
  switch (status) {
    case TernarySimulationStatus::Complete:
      return "complete";
    case TernarySimulationStatus::InvalidSystem:
      return "invalid-system";
    case TernarySimulationStatus::StepLimitExceeded:
      return "step-limit-exceeded";
    case TernarySimulationStatus::MemoryLimitExceeded:
      return "memory-limit-exceeded";
    case TernarySimulationStatus::EvaluationFailed:
      return "evaluation-failed";
  }
  return "invalid-system";
}

TernarySimulationResult simulateTernary(
    const NormalizedTransitionSystem& system,
    const TernarySimulationOptions& options) {
  TernarySimulationResult result;
  result.statistics.stateWidth = system.states().size();
  result.statistics.expressionNodes = system.expressions().nodes().size();
  result.statistics.effectiveMaxConcurrency =
      effectiveConcurrency(options.maxConcurrency);

  std::string validationReason;
  if (system.hasFatalDiagnostics() || !system.validate(&validationReason)) {
    result.status = TernarySimulationStatus::InvalidSystem;
    result.diagnostic =
        system.hasFatalDiagnostics()
            ? "normalized transition system is incomplete"
            : "invalid normalized transition system: " + validationReason;
    return result;
  }
  for (size_t stateIndex = 0; stateIndex < system.states().size();
       ++stateIndex) {
    if (!isValidTernaryValue(system.states()[stateIndex].initialValue)) {
      result.status = TernarySimulationStatus::InvalidSystem;
      result.diagnostic = "invalid normalized transition system: state " +
                          std::to_string(stateIndex) +
                          " has an invalid initial value";
      return result;
    }
  }
  if (!canStoreState(0, system.states().size(), options.maxStoredStateBytes)) {
    result.status = TernarySimulationStatus::MemoryLimitExceeded;
    result.diagnostic = "ternary trace memory limit cannot hold initial state";
    return result;
  }

  TernaryState initial;
  initial.reserve(system.states().size());
  for (const auto& state : system.states()) {
    initial.push_back(state.initialValue);
  }
  result.trace.push_back(initial);
  result.statistics.storedStates = 1;
  std::unordered_map<TernaryState, size_t, TernaryStateHash> firstOccurrence;
  firstOccurrence.emplace(initial, 0);

  std::vector<TernaryValue> variableValues(system.variables().size(),
                                           TernaryValue::Unknown);
  const size_t effective = result.statistics.effectiveMaxConcurrency;
  const int arenaConcurrency = static_cast<int>(effective);
  tbb::task_arena arena(arenaConcurrency);
  tbb::enumerable_thread_specific<EvaluationScratch> threadScratch;

  for (size_t step = 0; step < options.maxSteps; ++step) {
    std::fill(variableValues.begin(), variableValues.end(),
              TernaryValue::Unknown);
    const auto& current = result.trace.back();
    for (size_t stateIndex = 0; stateIndex < system.states().size();
         ++stateIndex) {
      variableValues[system.states()[stateIndex].variable] =
          current[stateIndex];
    }

    TernaryState next(system.states().size(), TernaryValue::Unknown);
    // One failure slot per state makes diagnostics independent of scheduler
    // order: after the barrier we always report the lowest failing state index.
    std::vector<std::string> failures(system.states().size());
    auto evaluateRange = [&](size_t begin, size_t end) {
      auto& scratch = threadScratch.local();
      for (size_t stateIndex = begin; stateIndex < end; ++stateIndex) {
        try {
          next[stateIndex] =
              evaluateWithScratch(system, system.states()[stateIndex].nextState,
                                  variableValues, scratch);
        } catch (const std::exception& exception) {
          failures[stateIndex] = exception.what();
        } catch (...) {
          failures[stateIndex] = "unknown expression evaluation failure";
        }
      }
    };
    try {
      arena.execute([&] {
        if (effective == 1 || system.states().size() <= 1) {
          evaluateRange(0, system.states().size());
        } else {
          tbb::parallel_for(
              tbb::blocked_range<size_t>(0, system.states().size()),
              [&](const tbb::blocked_range<size_t>& range) {
                evaluateRange(range.begin(), range.end());
              },
              tbb::static_partitioner{});
        }
      });
    } catch (const std::exception& exception) {
      result.status = TernarySimulationStatus::EvaluationFailed;
      result.diagnostic = "ternary expression evaluation failed: " +
                          std::string(exception.what());
      return result;
    } catch (...) {
      result.status = TernarySimulationStatus::EvaluationFailed;
      result.diagnostic = "ternary expression evaluation failed";
      return result;
    }
    const auto failed = std::find_if(
        failures.begin(), failures.end(),
        [](const std::string& message) { return !message.empty(); });
    if (failed != failures.end()) {
      result.status = TernarySimulationStatus::EvaluationFailed;
      result.diagnostic =
          "ternary expression evaluation failed for state " +
          std::to_string(static_cast<size_t>(failed - failures.begin())) +
          ": " + *failed;
      return result;
    }

    ++result.statistics.simulatedSteps;
    if (const auto repeated = firstOccurrence.find(next);
        repeated != firstOccurrence.end()) {
      result.stemLength = repeated->second;
      result.cycleLength = result.trace.size() - repeated->second;
      result.status = TernarySimulationStatus::Complete;
      return result;
    }
    if (!canStoreState(result.trace.size(), system.states().size(),
                       options.maxStoredStateBytes)) {
      result.status = TernarySimulationStatus::MemoryLimitExceeded;
      result.diagnostic = "ternary trace memory limit exceeded";
      return result;
    }
    const size_t index = result.trace.size();
    result.trace.push_back(std::move(next));
    firstOccurrence.emplace(result.trace.back(), index);
    result.statistics.storedStates = result.trace.size();
  }

  result.status = TernarySimulationStatus::StepLimitExceeded;
  result.diagnostic = "ternary reset simulation step limit exceeded";
  return result;
}

}  // namespace KEPLER_FORMAL::SEC::PHASE
