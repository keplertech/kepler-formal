// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "kinduction/OutputBatching.h"

#include <algorithm>
#include <unordered_set>

#include "common/BoolExprUtils.h"

namespace KEPLER_FORMAL::SEC {

namespace {

// A single monolithic OR-of-all-bads query can be too wide for ASIC netlists,
// but proving one output per solver call repeats setup work hundreds of times.
// These limits keep nearby outputs together while preventing one batch from
// dragging most of the design into one SAT cone.
constexpr OutputBatchingLimits kDefaultOutputBatchingLimits;
constexpr OutputBatchingLimits kHugeDualRailOutputBatchingLimits{1, 1};
constexpr OutputBatchingLimits kDualRailOutputBatchingLimits{8, 256};
constexpr OutputBatchingLimits kSmallDualRailOutputBatchingLimits{32, 4096};
constexpr size_t kMaxSmallDualRailBatchOutputs = 32;
constexpr size_t kMaxSmallDualRailBatchStateSymbols = 512;
// Swerv has 45,096 same-design rail pairs, i.e. 90,192 rail symbols after
// dual-rail expansion.  Treat that surface like BP-sized residuals so strict KI
// starts directly at one-output leaves instead of rebuilding broad UNKNOWN
// batches only to split them.
constexpr size_t kMinHugeDualRailBatchStateSymbols = 65536;

void appendOutputSupport(const KInductionProblem& problem,
                         size_t outputIndex,
                         std::unordered_set<size_t>& support) {
  for (const auto symbol : problem.observedOutputExprs0[outputIndex]->getSupportVars()) {
    support.insert(symbol);
  }
  for (const auto symbol : problem.observedOutputExprs1[outputIndex]->getSupportVars()) {
    support.insert(symbol);
  }
}

size_t dualRailStateSymbolCount(const KInductionProblem& problem) {
  if (!problem.dualRailStatePairs.empty()) {
    return problem.dualRailStatePairs.size() * 2;
  }
  return problem.state0Symbols.size() + problem.state1Symbols.size();
}

bool isSmallDualRailBatchingSurface(const KInductionProblem& problem) {
  return problem.observedOutputExprs0.size() <= kMaxSmallDualRailBatchOutputs &&
         dualRailStateSymbolCount(problem) <=
             kMaxSmallDualRailBatchStateSymbols;
}

bool isHugeDualRailBatchingSurface(const KInductionProblem& problem) {
  return dualRailStateSymbolCount(problem) >=
         kMinHugeDualRailBatchStateSymbols;
}

bool hasSelectedOutputSkips(const KInductionProblem& problem) {
  if (problem.dualRailOutputSkipReasons.size() !=
      problem.observedOutputExprs0.size()) {
    return false;
  }
  for (const auto& reason : problem.dualRailOutputSkipReasons) {
    if (!reason.empty()) {
      return true; // LCOV_EXCL_LINE
    }
  }
  return false;
}

void rememberDualRailResidualPublicHypothesis(const KInductionProblem& problem) {
  if (!problem.usesDualRailStateEncoding || problem.lazyTransitions == nullptr ||
      problem.property == nullptr || problem.bad == nullptr ||
      problem.observedOutputExprs0.size() <= 1 ||
      hasSelectedOutputSkips(problem)) {
    return;
  }

  LazyTransitionStore& store = *problem.lazyTransitions;
  if (store.dualRailResidualPublicOutputCount >=
      problem.observedOutputExprs0.size()) {
    return;
  }

  // Store only the public residual output conjunction.  Later KI slices use it
  // as their induction hypothesis; no internal cross-design relation is added.
  store.dualRailResidualPublicProperty = problem.property;
  store.dualRailResidualPublicBad = problem.bad;
  store.dualRailResidualPublicOutputCount = problem.observedOutputExprs0.size();
}

}  // namespace

std::vector<std::pair<size_t, size_t>> buildSupportBoundedOutputBatches(
    const KInductionProblem& problem) {
  return buildSupportBoundedOutputBatches(
      problem, defaultOutputBatchingLimitsForProblem(problem));
}

std::vector<std::pair<size_t, size_t>> buildSupportBoundedOutputBatches(
    const KInductionProblem& problem,
    const OutputBatchingLimits& limits) {
  std::vector<std::pair<size_t, size_t>> batches;
  size_t firstOutput = 0;
  std::unordered_set<size_t> batchSupport;
  batchSupport.reserve(limits.outputBatchSupportLimit);

  while (firstOutput < problem.observedOutputExprs0.size()) {
    size_t endOutput = firstOutput;
    batchSupport.clear();
    while (endOutput < problem.observedOutputExprs0.size()) {
      std::unordered_set<size_t> candidateSupport = batchSupport;
      appendOutputSupport(problem, endOutput, candidateSupport);

      const bool batchAlreadyHasOutput = endOutput > firstOutput;
      const bool exceedsCount =
          endOutput - firstOutput + 1 > limits.maxOutputBatchSize;
      const bool exceedsSupport =
          candidateSupport.size() > limits.outputBatchSupportLimit;
      if (batchAlreadyHasOutput && (exceedsCount || exceedsSupport)) {
        break;
      }

      batchSupport = std::move(candidateSupport);
      ++endOutput;
    }
    batches.emplace_back(firstOutput, endOutput);
    firstOutput = endOutput;
  }

  return batches;
}

OutputBatchingLimits defaultOutputBatchingLimitsForProblem(
    const KInductionProblem& problem) {
  if (problem.usesDualRailStateEncoding) {
    rememberDualRailResidualPublicHypothesis(problem);
    if (isSmallDualRailBatchingSurface(problem)) {
      // GCD-sized dual-rail designs often need the public output conjunction as
      // the strict KI property.  Keep those small surfaces together; large
      // rail-state ASICs still use tiny batches to avoid wide failed probes.
      return kSmallDualRailOutputBatchingLimits;
    }
    if (isHugeDualRailBatchingSurface(problem)) {
      // BP-sized residuals repeatedly hit UNKNOWN on 8/4/2-output probes before
      // reaching the same one-output strict KI leaf. Start at that leaf size so
      // no proof strength changes and no UNKNOWN batch is rebuilt only to split.
      return kHugeDualRailOutputBatchingLimits;
    }
    // Dual-rail output obligations already carry both may-one/may-zero rails.
    // Start with small exact OR batches.  Wide dual-rail OR batches often hit
    // the resource-limited step path only after building a large transition
    // CNF, and then immediately split; tiny batches keep useful grouping
    // without paying for the failed wide proof first.
    return kDualRailOutputBatchingLimits;
  }
  return kDefaultOutputBatchingLimits;
}

void configureOutputBatchProblem(KInductionProblem& batch,
                                 const KInductionProblem& source,
                                 size_t firstOutput,
                                 size_t endOutput) {
  batch.originalObservedOutputCount =
      source.originalObservedOutputCount == 0
          ? source.observedOutputExprs0.size()
          : source.originalObservedOutputCount;
  if (source.observedOutputs.size() == source.observedOutputExprs0.size()) {
    batch.observedOutputs.assign(
        source.observedOutputs.begin() + firstOutput,
        source.observedOutputs.begin() + endOutput);
  } else {
    batch.observedOutputs.clear();
  }
  batch.observedOutputNames.assign(
      source.observedOutputNames.begin() + firstOutput,
      source.observedOutputNames.begin() + endOutput);
  batch.observedOutputExprs0.assign(
      source.observedOutputExprs0.begin() + firstOutput,
      source.observedOutputExprs0.begin() + endOutput);
  batch.observedOutputExprs1.assign(
      source.observedOutputExprs1.begin() + firstOutput,
      source.observedOutputExprs1.begin() + endOutput);
  if (source.dualRailOutputStrictEqualityExprs.size() ==
      source.observedOutputExprs0.size()) {
    batch.dualRailOutputStrictEqualityExprs.assign(
        source.dualRailOutputStrictEqualityExprs.begin() + firstOutput,
        source.dualRailOutputStrictEqualityExprs.begin() + endOutput);
  } else {
    batch.dualRailOutputStrictEqualityExprs.clear();
  }
  if (source.dualRailOutputSkipReasons.size() ==
      source.observedOutputExprs0.size()) {
    batch.dualRailOutputSkipReasons.assign(
        source.dualRailOutputSkipReasons.begin() + firstOutput,
        source.dualRailOutputSkipReasons.begin() + endOutput);
  } else {
    batch.dualRailOutputSkipReasons.clear();
  }
  // Every caller creates the reusable batch from `source` before selecting an
  // output range.  Same-design state relations are immutable model data, so
  // leave those vectors in place instead of copying the whole design per range.

  // SEC output equality is a conjunction. Proving smaller conjunctions and
  // combining the results is logically equivalent to one monolithic property,
  // while allowing the base/induction encoders to run COI on much smaller
  // output cones.
  BoolExpr* property = BoolExpr::createTrue();
  for (size_t i = 0; i < batch.observedOutputExprs0.size(); ++i) {
    property = BoolExpr::And(
        property,
        makeEqualityExpr(batch.observedOutputExprs0[i], batch.observedOutputExprs1[i]));
  }
  batch.property = BoolExpr::simplify(property);
  batch.bad = BoolExpr::simplify(BoolExpr::Not(batch.property));

  // Prove only this output slice in the induction query.  The caller combines
  // successful slices as the SEC output conjunction, while each slice remains a
  // normal strict k-induction obligation over its own base and step cases.
  batch.inductionProperty = batch.property;
  batch.inductionBad = batch.bad;
  batch.description = source.description + " output batch";
}

}  // namespace KEPLER_FORMAL::SEC
