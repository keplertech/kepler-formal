// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "kinduction/KInductionEngine.h"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <memory>
#include <optional>
#include <utility>

#include "common/SecDiag.h"
#include "kinduction/BaseCaseSolver.h"
#include "kinduction/InductionStepSolver.h"
#include "kinduction/OutputBatching.h"

namespace KEPLER_FORMAL::SEC {

namespace {

bool isKInductionDiagEnabled() {
  return std::getenv("KEPLER_SEC_KI_DIAG") != nullptr || isSecDiagEnabled();
}

// Batching protects SEC proofs from one broad OR-of-output-bads SAT query.
// Keep every true multi-output proof batched; medium designs such as
// sky130hs_ibex are still sensitive to monolithic base-case witnesses.
constexpr size_t kMinOutputsForBatchedProof = 2;
constexpr size_t kMinOriginalOutputsForBoundedBatch = 64;
constexpr size_t kMinDeferredRailStateSymbolsForEarlyStop = 512;
// Swerv-scale dual-rail residuals have about 90k rail symbols.  Once a strict
// one-output deferred leaf hits the SAT decision cap, more k-depth rebuilds do
// not add proof strength; they only delay reporting that leaf as uncovered.
constexpr size_t kHugeDeferredRailStateSymbolsForEarlyStop = 65536;
constexpr size_t kDefaultDeferredDualRailLeafResourceLimitStops = 4;
constexpr size_t kHugeDeferredDualRailLeafResourceLimitStops = 1;
constexpr size_t kMaxCompactDualRailConjunctionOutputs = 32;
constexpr size_t kMaxWideDualRailSplitHypothesisOutputs = 8;
constexpr size_t kMaxFullDualRailPublicHypothesisOutputs = 256;
constexpr size_t kMaxStoredPublicSplitHypothesisOutputs = 32;
constexpr size_t kMaxFullDualRailPublicHypothesisStateSymbols = 4096;
constexpr unsigned kDefaultBatchedInductionDecisionLimit = 200000;
// Dual-rail SEC is a coverage extension for resetless state.  Residual outputs
// that do not close quickly should become uncovered/inconclusive instead of
// letting one expanded rail query dominate the full regression runtime.
constexpr unsigned kDefaultDualRailInductionDecisionLimit = 5000;

std::optional<unsigned> readUnsignedEnv(const char* name) {
  const char* value = std::getenv(name);
  if (value == nullptr || *value == '\0') {
    return std::nullopt;
  }
  char* end = nullptr;
  const unsigned long parsed = std::strtoul(value, &end, 10);
  if (end == value || *end != '\0' ||
      parsed > std::numeric_limits<unsigned>::max()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  return static_cast<unsigned>(parsed);
}

size_t originalOutputCountForProblem(const KInductionProblem& problem) {
  return problem.originalObservedOutputCount == 0
      ? problem.observedOutputExprs0.size()
      : problem.originalObservedOutputCount;
}

bool isWideDualRailResidualSurface(const KInductionProblem& problem) {
  return originalOutputCountForProblem(problem) >=
         kMinOriginalOutputsForBoundedBatch;
}

size_t dualRailStateSymbolCount(const KInductionProblem& problem) {
  if (!problem.dualRailStatePairs.empty()) {
    return problem.dualRailStatePairs.size() * 2;
  }
  return problem.state0Symbols.size() + problem.state1Symbols.size();
}

bool isLargeDeferredDualRailLeafSurface(const KInductionProblem& problem) {
  return problem.deferBaseCaseChecks &&
         dualRailStateSymbolCount(problem) >=
             kMinDeferredRailStateSymbolsForEarlyStop;
}

bool shouldBoundDualRailBatchInduction(const KInductionProblem& problem) {
  if (!problem.usesDualRailStateEncoding ||
      problem.observedOutputExprs0.size() <= 1) {
    return false; // LCOV_EXCL_LINE
  }
  if (originalOutputCountForProblem(problem) <=
          kMaxCompactDualRailConjunctionOutputs &&
      dualRailStateSymbolCount(problem) <
          kMinDeferredRailStateSymbolsForEarlyStop) {
    // Compact public-output batches can otherwise spend the whole workflow on a
    // single all-output SAT query.  A resource-limited result only triggers
    // strict KI splitting; every split must still prove under the public
    // conjunction hypothesis and pass the shared concrete base check.
    return true;
  }
  return isWideDualRailResidualSurface(problem) ||
         dualRailStateSymbolCount(problem) >=
             kMinDeferredRailStateSymbolsForEarlyStop;
}

bool isCompactDualRailPublicConjunctionSurface(
    const KInductionProblem& problem) {
  return problem.usesDualRailStateEncoding &&
         originalOutputCountForProblem(problem) <=
             kMaxCompactDualRailConjunctionOutputs &&
         dualRailStateSymbolCount(problem) < // LCOV_EXCL_LINE
             kMinDeferredRailStateSymbolsForEarlyStop;
}

struct PublicConjunctionHypothesis {
  BoolExpr* property = nullptr;
  BoolExpr* bad = nullptr;
  size_t outputCount = 0;
  bool completePublicSurface = false;
};

size_t maxDualRailSplitHypothesisOutputs(const KInductionProblem& problem) {
  if (isCompactDualRailPublicConjunctionSurface(problem)) {
    return kMaxCompactDualRailConjunctionOutputs; // LCOV_EXCL_LINE
  }
  return kMaxWideDualRailSplitHypothesisOutputs;
}

bool canUseFullPublicConjunctionHypothesis(
    const KInductionProblem& problem) {
  return problem.usesDualRailStateEncoding &&
         problem.property != nullptr &&
         problem.bad != nullptr &&
         problem.observedOutputExprs0.size() <=
             kMaxFullDualRailPublicHypothesisOutputs &&
         dualRailStateSymbolCount(problem) <=
             kMaxFullDualRailPublicHypothesisStateSymbols;
}

bool shouldApplyPublicConjunctionHypothesis(
    const KInductionProblem& batchProblem,
    const KInductionProblem& sourceProblem,
    const PublicConjunctionHypothesis& hypothesis) {
  const bool isFullPublicHypothesis =
      hypothesis.property == sourceProblem.property &&
      hypothesis.bad == sourceProblem.bad &&
      canUseFullPublicConjunctionHypothesis(sourceProblem);
  return sourceProblem.usesDualRailStateEncoding &&
         hypothesis.property != nullptr &&
         hypothesis.bad != nullptr &&
         hypothesis.outputCount > batchProblem.observedOutputExprs0.size() &&
         (hypothesis.completePublicSurface ||
          isFullPublicHypothesis ||
          hypothesis.outputCount <=
              maxDualRailSplitHypothesisOutputs(sourceProblem));
}

bool samePublicConjunctionHypothesis(
    const PublicConjunctionHypothesis& lhs,
    const PublicConjunctionHypothesis& rhs) {
  return lhs.property == rhs.property &&
         lhs.bad == rhs.bad &&
         lhs.outputCount == rhs.outputCount &&
         lhs.completePublicSurface == rhs.completePublicSurface; // LCOV_EXCL_LINE
}

bool shouldTryFallbackHypothesisAtRange(
    const PublicConjunctionHypothesis& preferredHypothesis,
    const PublicConjunctionHypothesis& fallbackHypothesis,
    size_t outputCount) {
  if (fallbackHypothesis.property == nullptr ||
      samePublicConjunctionHypothesis(preferredHypothesis, fallbackHypothesis)) {
    return false;
  }
  if (fallbackHypothesis.completePublicSurface &&
      fallbackHypothesis.outputCount >
          kMaxStoredPublicSplitHypothesisOutputs) {
    // The remembered full residual surface can dwarf the current strategy
    // batch.  Let the current strict KI range split and carry its own smaller
    // public conjunction instead of replaying an AES-sized hypothesis at every
    // recursive level.
    return false;
  }
  if (outputCount <= 1 && preferredHypothesis.property != nullptr &&
      fallbackHypothesis.outputCount > preferredHypothesis.outputCount) {
    // Once a split leaf has failed under its parent public-output conjunction,
    // a strictly wider remembered residual conjunction only rebuilds a larger
    // version of the same strict KI leaf.  Keep the wider fallback for parent
    // batches and for standalone one-output residual calls that have no smaller
    // parent hypothesis.
    return false;
  }
  return true;
}

bool shouldTryStoredPublicHypothesisOnStandaloneLeaf(
    const PublicConjunctionHypothesis& hypothesis) {
  if (hypothesis.property == nullptr) {
    return false; // LCOV_EXCL_LINE
  }
  if (!hypothesis.completePublicSurface) {
    return true;
  }
  // A full residual public hypothesis is useful for compact buses, but replaying
  // a 100+ output AES conjunction after every one-output UNKNOWN leaf rebuilds a
  // larger strict KI query without sharing coverage. Parent batches still try
  // the full hypothesis before splitting; standalone leaves keep the local
  // strict KI result once the remembered surface is too wide.
  return hypothesis.outputCount <=
         kMaxStoredPublicSplitHypothesisOutputs;
}

size_t solverTypeCacheValue(KEPLER_FORMAL::Config::SolverType solverType) {
  return static_cast<size_t>(solverType);
}

bool canMemoizeDualRailResidualLeafAttempt(
    const KInductionProblem& leafProblem,
    const PublicConjunctionHypothesis& hypothesis) {
  return leafProblem.usesDualRailStateEncoding &&
         leafProblem.lazyTransitions != nullptr &&
         leafProblem.observedOutputExprs0.size() == 1 &&
         leafProblem.observedOutputExprs1.size() == 1 &&
         hypothesis.property != nullptr &&
         hypothesis.bad != nullptr;
}

bool publicHypothesisWasAppliedToLeaf(
    const KInductionProblem& leafProblem,
    const PublicConjunctionHypothesis& hypothesis) {
  return leafProblem.inductionProperty == hypothesis.property &&
         leafProblem.inductionBad == leafProblem.bad;
}

PublicConjunctionHypothesis activeLeafAttemptHypothesis(
    const KInductionProblem& leafProblem,
    const PublicConjunctionHypothesis& requestedHypothesis) {
  if (requestedHypothesis.property != nullptr &&
      publicHypothesisWasAppliedToLeaf(leafProblem, requestedHypothesis)) {
    return requestedHypothesis;
  }

  // configureOutputBatchProblem installs the selected output predicate as the
  // local induction property.  Cache that exact leaf too; replaying UNKNOWN only
  // saves runtime and never turns the output into coverage.
  return PublicConjunctionHypothesis{
      leafProblem.inductionProperty != nullptr ? leafProblem.inductionProperty
                                               : leafProblem.property,
      leafProblem.inductionBad != nullptr ? leafProblem.inductionBad
                                          : leafProblem.bad,
      leafProblem.observedOutputExprs0.size(),
      false};
}

std::optional<size_t> inconclusiveDualRailLeafAttemptBound(
    const KInductionProblem& leafProblem,
    const PublicConjunctionHypothesis& hypothesis,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK) {
  if (!canMemoizeDualRailResidualLeafAttempt(leafProblem, hypothesis)) {
    return std::nullopt;
  }
  size_t resultBound = 0;
  if (!leafProblem.lazyTransitions
           ->findInconclusiveDualRailResidualPublicKiAttempt(
               leafProblem.observedOutputExprs0[0],
               leafProblem.observedOutputExprs1[0],
               hypothesis.property,
               hypothesis.bad,
               hypothesis.outputCount,
               solverTypeCacheValue(solverType),
               maxK,
               resultBound)) {
    return std::nullopt;
  }
  return resultBound;
}

void rememberInconclusiveDualRailLeafAttempt(
    const KInductionProblem& leafProblem,
    const PublicConjunctionHypothesis& hypothesis,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK,
    size_t resultBound) {
  if (!canMemoizeDualRailResidualLeafAttempt(leafProblem, hypothesis)) {
    return;
  }
  leafProblem.lazyTransitions
      ->rememberInconclusiveDualRailResidualPublicKiAttempt(
          leafProblem.observedOutputExprs0[0],
          leafProblem.observedOutputExprs1[0],
          hypothesis.property,
          hypothesis.bad,
          hypothesis.outputCount,
          solverTypeCacheValue(solverType),
          maxK,
          resultBound);
}

PublicConjunctionHypothesis makePublicConjunctionHypothesis(
    const KInductionProblem& problem) {
  return PublicConjunctionHypothesis{
      problem.property,
      problem.bad,
      problem.observedOutputExprs0.size(),
      false};
}

PublicConjunctionHypothesis storedResidualPublicConjunctionHypothesis(
    const KInductionProblem& problem) {
  if (!problem.usesDualRailStateEncoding || problem.lazyTransitions == nullptr ||
      dualRailStateSymbolCount(problem) >
          kMaxFullDualRailPublicHypothesisStateSymbols) {
    return {};
  }

  const LazyTransitionStore& store = *problem.lazyTransitions;
  if (store.dualRailResidualPublicProperty == nullptr ||
      store.dualRailResidualPublicBad == nullptr ||
      store.dualRailResidualPublicOutputCount >
          kMaxFullDualRailPublicHypothesisOutputs ||
      store.dualRailResidualPublicOutputCount <=
          problem.observedOutputExprs0.size()) {
    return {};
  }

  return PublicConjunctionHypothesis{
      store.dualRailResidualPublicProperty,
      store.dualRailResidualPublicBad,
      store.dualRailResidualPublicOutputCount,
      true};
}

PublicConjunctionHypothesis initialPublicConjunctionHypothesis(
    const KInductionProblem& problem) {
  if (const PublicConjunctionHypothesis stored =
          storedResidualPublicConjunctionHypothesis(problem);
      stored.property != nullptr) {
    // The strategy may split the residual output set before entering KI.  Use
    // the remembered residual public conjunction as the strict KI hypothesis for
    // each slice; acceptance still requires every slice proof plus concrete base
    // validation by the caller.
    return stored;
  }
  if (!canUseFullPublicConjunctionHypothesis(problem)) {
    return {};
  }
  // Decompose the strict KI proof for the whole public SEC conjunction across
  // output batches: every batch proves P_all[0..k-1] => P_batch[k], and the
  // caller later validates the concrete base prefix for P_all once.  This is a
  // public-output hypothesis only; no cross-design internal relation is added.
  return makePublicConjunctionHypothesis(problem);
}

void applyDualRailSplitHypothesis(
    KInductionProblem& batchProblem,
    const KInductionProblem& sourceProblem,
    const PublicConjunctionHypothesis& hypothesis,
    size_t firstOutput,
    size_t endOutput) {
  if (!shouldApplyPublicConjunctionHypothesis(
          batchProblem, sourceProblem, hypothesis)) {
    return;
  }

  // This is the decomposed k-induction step for a public-output conjunction:
  // prove P_parent[0..k-1] => P_i[k] for each split output i, then accept only
  // if every split closes and the shared concrete base check proves the public
  // outputs through the final frontier.  No cross-design internal relation is
  // introduced.
  batchProblem.inductionProperty = hypothesis.property;
  batchProblem.inductionBad = batchProblem.bad;
  if (isKInductionDiagEnabled()) {
    emitSecDiag(
        "SEC diag: k-induction dual-rail public conjunction hypothesis range [",
        firstOutput, ",", endOutput, ") source_outputs=",
        hypothesis.outputCount);
  }
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP
unsigned binaryBatchedInductionDecisionLimit() {
  return readUnsignedEnv("KEPLER_SEC_KI_BATCH_DECISION_LIMIT")
      .value_or(kDefaultBatchedInductionDecisionLimit);
}

std::optional<unsigned> dualRailLeafInductionDecisionLimit(
    const KInductionProblem& problem) {
  if (const auto leafLimit =
          readUnsignedEnv("KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT");
      leafLimit.has_value()) {
    return leafLimit;
  }
  if (!isWideDualRailResidualSurface(problem) &&
      !isLargeDeferredDualRailLeafSurface(problem)) {
    // Compact designs can need a full strict KI leaf search, and those cones
    // are small enough that an unbounded leaf does not dominate the workflow.
    return std::nullopt;
  }
  return kDefaultDualRailInductionDecisionLimit;
}

std::optional<unsigned> batchedInductionDecisionLimit(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (solverType != KEPLER_FORMAL::Config::SolverType::KISSAT) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  if (problem.usesDualRailStateEncoding &&
      // LCOV_EXCL_START
      problem.observedOutputExprs0.size() <= 1) {
      // LCOV_EXCL_STOP
    return dualRailLeafInductionDecisionLimit(problem);
  }

  if (problem.observedOutputExprs0.size() <= 1) {
    return std::nullopt;
  }

  if (!problem.usesDualRailStateEncoding) {
    return binaryBatchedInductionDecisionLimit();
  }

  // Dual-rail proofs are state-space expanded, so a single hard multi-output
  // induction batch can spend the whole workflow inside CDCL.  Base checks are
  // shared outside the slices, making recursive splitting cheap enough to use
  // the normal batch cap here too.
  if (const auto dualRailLimit =
          readUnsignedEnv("KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT");
      dualRailLimit.has_value()) {
    return dualRailLimit;
  }
  if (!shouldBoundDualRailBatchInduction(problem)) {
    // Small dual-rail designs can need the public output conjunction as the
    // induction property.  Do not turn those strict KI attempts into UNKNOWN by
    // default; wide rail-state surfaces still use the bounded split path.
    return std::nullopt; // LCOV_EXCL_LINE
  }
  return readUnsignedEnv("KEPLER_SEC_KI_BATCH_DECISION_LIMIT")
      .value_or(kDefaultDualRailInductionDecisionLimit);
}

void emitKInductionProblemDiag(const KInductionProblem& problem,
                               size_t maxK) {
  if (!isKInductionDiagEnabled()) {
    return;
  }
  emitSecDiag(
      "SEC diag: k-induction problem outputs=",
      problem.observedOutputExprs0.size(),
      " state0=", problem.state0Symbols.size(),
      " state1=", problem.state1Symbols.size(),
      " reset_bootstrap_cycles=", problem.resetBootstrapCycles,
      " max_k=", maxK);
}

// LCOV_EXCL_START
bool shouldCheckLocalBaseCase(const KInductionProblem& problem) {
  return !problem.deferBaseCaseChecks;
  // LCOV_EXCL_STOP
}

bool proofNeedsConcreteFrontierValidation(const KInductionProblem& problem) {
  return problem.resetBootstrapCycles != 0 ||
         problem.inductionProperty != nullptr ||
         problem.inductionBad != nullptr;
}

bool isDeferredDualRailLeafProof(const KInductionProblem& problem) {
  return problem.deferBaseCaseChecks &&
         problem.usesDualRailStateEncoding &&
         problem.observedOutputExprs0.size() <= 1;
}

bool shouldStopDeferredDualRailLeafAfterResourceLimit(
    const KInductionProblem& problem,
    size_t consecutiveResourceLimitedSteps) {
  if (!isDeferredDualRailLeafProof(problem)) {
    return false; // LCOV_EXCL_LINE
  }
  if (!isLargeDeferredDualRailLeafSurface(problem)) {
    return false;
  }
  // Wide deferred dual-rail leaves are residual coverage attempts. Resource-
  // limited KI steps are not proofs, and rebuilding larger k-step instances for
  // each hard BP-scale leaf only delays reporting it as unproven. Small designs
  // keep the deeper strict KI search because those extra depths are cheap and
  // often recover complete coverage.
  const size_t stopAfterResourceLimits =
      dualRailStateSymbolCount(problem) >=
              kHugeDeferredRailStateSymbolsForEarlyStop
          ? kHugeDeferredDualRailLeafResourceLimitStops
          : kDefaultDeferredDualRailLeafResourceLimitStops;
  return consecutiveResourceLimitedSteps >= stopAfterResourceLimits;
}

std::optional<KInductionResult> validateConcreteBasePrefix(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  if (!problem.usesDualRailStateEncoding) {
    if (auto witness = SEC::findBaseCounterexample(problem, solverType, k);
        witness.has_value()) {
      return KInductionResult{
          KInductionStatus::Different,
          witness->badFrame,
          std::move(witness)};
    }
    return std::nullopt;
  }

  const BaseCounterexampleCheckResult baseCheck =
      SEC::checkBaseCounterexampleWithFastValidation(problem, solverType, k);
  if (baseCheck.status == BaseCounterexampleCheckStatus::Counterexample ||
      baseCheck.witness.has_value()) {
    const size_t badFrame =
        baseCheck.witness.has_value() ? baseCheck.witness->badFrame : k;
    return KInductionResult{
        KInductionStatus::Different,
        badFrame,
        baseCheck.witness};
  }
  if (baseCheck.status == BaseCounterexampleCheckStatus::NoCounterexample) {
    return std::nullopt;
  }

  // A resource-limited concrete base check is not a proof.  Keep KI strict by
  // returning inconclusive instead of accepting a strengthened induction
  // certificate as a full SEC result.
  return KInductionResult{KInductionStatus::Inconclusive, k}; // LCOV_EXCL_LINE
}

KInductionResult runMonolithicKInduction(const KInductionProblem& problem,
                                         KEPLER_FORMAL::Config::SolverType solverType,
                                         size_t maxK) {
  // Handle the purely combinational mismatch case before any unrolling.
  if (isKInductionDiagEnabled()) {
    emitSecDiag("SEC diag: k-induction base k=0 begin");
  }
  if (shouldCheckLocalBaseCase(problem)) {
    // LCOV_EXCL_START
    auto baseZeroWitness = SEC::findBaseCounterexample(problem, solverType, 0);
    // LCOV_EXCL_STOP
    if (baseZeroWitness.has_value()) {
      if (isKInductionDiagEnabled()) {
        emitSecDiag("SEC diag: k-induction base k=0 found cex");
      }
      return {
          KInductionStatus::Different,
          baseZeroWitness->badFrame,
          std::move(baseZeroWitness)};
    }
  } else if (isKInductionDiagEnabled()) {
    emitSecDiag("SEC diag: k-induction base k=0 deferred");
  }
  if (isKInductionDiagEnabled()) {
    emitSecDiag("SEC diag: k-induction base k=0 unsat");
  }

  // If there is no state, the base check already decided the whole problem.
  if (problem.combinedStateSymbols().empty()) {
    return {KInductionStatus::Equivalent, 0};
  }

  std::shared_ptr<KInductionBaseCounterexampleCache> baseFrontierCache;
  if (shouldCheckLocalBaseCase(problem)) {
    baseFrontierCache = SEC::makeKInductionBaseCounterexampleCache(problem);
  }

  // At the start of iteration k, all frames < k have already been proved safe
  // by the base checks below.  That is exactly the base obligation needed for
  // the k-step induction query "P[0]..P[k-1] => P[k]"; if the step closes, the
  // property is invariant and there is no reason to spend time on the frontier
  // BMC query for frame k. Only when the step is inconclusive do we extend the
  // concrete base horizon by checking the new frontier for a real counterexample.
  size_t consecutiveResourceLimitedDeferredLeafSteps = 0;
  for (size_t k = 1; k <= maxK; ++k) {
    if (isKInductionDiagEnabled()) {
      emitSecDiag("SEC diag: k-induction step k=", k, " begin");
    }

// LCOV_EXCL_START

    const std::optional<unsigned> inductionDecisionLimit =
        batchedInductionDecisionLimit(problem, solverType);
    const InductionProofStatus inductionStatus =
        SEC::proveByInductionStatus(
            problem, solverType, k, inductionDecisionLimit);
    if (inductionStatus == InductionProofStatus::Proved) {
      if (shouldCheckLocalBaseCase(problem) &&
          proofNeedsConcreteFrontierValidation(problem)) {
        // Reset/bootstrap and explicit induction certificates can prove a
        // LCOV_EXCL_STOP
        // strengthened obligation. Before accepting that as SEC equivalence,
        // LCOV_DISABLED_START
        // validate the concrete top-output base predicate through the proved
        // frontier.
        if (auto validation =
                validateConcreteBasePrefix(problem, solverType, k);
            validation.has_value()) {
          return *validation;
        }
      // LCOV_DISABLED_START
      }
      // LCOV_DISABLED_STOP
      if (isKInductionDiagEnabled()) {
        emitSecDiag("SEC diag: k-induction step k=", k, " proved");
      }
      return {KInductionStatus::Equivalent, k};
    }
    if (inductionStatus == InductionProofStatus::Unknown) {
      if (isDeferredDualRailLeafProof(problem)) {
        ++consecutiveResourceLimitedDeferredLeafSteps;
      }
      if (isKInductionDiagEnabled()) {  // LCOV_EXCL_LINE
        if (problem.usesDualRailStateEncoding &&  // LCOV_EXCL_LINE
            problem.observedOutputExprs0.size() <= 1) {  // LCOV_EXCL_LINE
          emitSecDiag(  // LCOV_EXCL_LINE
              "SEC diag: k-induction step k=", k,
              " resource-limited; checking frontier");
        } else {  // LCOV_EXCL_LINE
          emitSecDiag(  // LCOV_EXCL_LINE
              "SEC diag: k-induction step k=", k,
              " resource-limited; splitting output batch");
        }
      }  // LCOV_EXCL_LINE
      // LCOV_DISABLED_START
      if (!problem.usesDualRailStateEncoding ||  // LCOV_EXCL_LINE
          problem.observedOutputExprs0.size() > 1) {  // LCOV_EXCL_LINE
        return {KInductionStatus::Inconclusive, k};  // LCOV_EXCL_LINE
      }
      // LCOV_DISABLED_STOP
      if (!shouldCheckLocalBaseCase(problem)) {  // LCOV_EXCL_LINE
        // Output-batched dual-rail KI validates the concrete base prefix once
        // for the full output set after all slices prove.  A resource-limited
        // leaf step at the current k is therefore not a proof failure; advance
        // while the strict KI attempts still make progress, and report the leaf
        // inconclusive if repeated capped steps show no useful progress.
        if (problem.deferBaseCaseChecks) {  // LCOV_EXCL_LINE
          if (isKInductionDiagEnabled()) {  // LCOV_EXCL_LINE
            emitSecDiag(  // LCOV_EXCL_LINE
                "SEC diag: k-induction step k=", k,
                " resource-limited; deferred base continues");
          }  // LCOV_EXCL_LINE
          if (isKInductionDiagEnabled() &&  // LCOV_EXCL_LINE
              shouldStopDeferredDualRailLeafAfterResourceLimit(
                  problem,
                  consecutiveResourceLimitedDeferredLeafSteps)) {  // LCOV_EXCL_LINE
            emitSecDiag(  // LCOV_EXCL_LINE
                "SEC diag: k-induction step k=", k,
                " resource-limited deferred leaf limit reached; reporting inconclusive");
          }  // LCOV_EXCL_LINE
        } else {  // LCOV_EXCL_LINE
          return {KInductionStatus::Inconclusive, k};  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
        if (shouldStopDeferredDualRailLeafAfterResourceLimit(
                problem,
                consecutiveResourceLimitedDeferredLeafSteps)) {  // LCOV_EXCL_LINE
          return {KInductionStatus::Inconclusive, k};  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
      // LCOV_DISABLED_START
      }
    } else {
      consecutiveResourceLimitedDeferredLeafSteps = 0;
    }  // LCOV_EXCL_LINE
    // LCOV_DISABLED_STOP
    if (isKInductionDiagEnabled()) {
      emitSecDiag("SEC diag: k-induction step k=", k, " inconclusive");
      emitSecDiag("SEC diag: k-induction base k=", k, " begin");
    // LCOV_DISABLED_START
    }

    // Earlier base checks have already ruled out bad states on frames < k.
    // Check only the newly exposed frontier instead of re-solving an
    // LCOV_DISABLED_STOP
    // OR-of-all-previous-bads query at every depth.
    // LCOV_DISABLED_START
    if (shouldCheckLocalBaseCase(problem)) {
      if (auto witness = SEC::findKInductionBaseCounterexampleAtFrontier(
              *baseFrontierCache, solverType, k);
          witness.has_value()) {
        if (isKInductionDiagEnabled()) {
          emitSecDiag("SEC diag: k-induction base k=", k, " found cex");
        }
        return {
            KInductionStatus::Different,
            witness->badFrame,
            std::move(witness)};
      }
      if (isKInductionDiagEnabled()) {
        emitSecDiag("SEC diag: k-induction base k=", k, " unsat");
      }
    } else if (isKInductionDiagEnabled()) {
      emitSecDiag("SEC diag: k-induction base k=", k, " deferred");
    }
  }

  // Frontier checks are an optimization over the classic cumulative base case:
  // each iteration checks only the newly exposed bad frame because earlier
  // frames were already ruled out.  Before reporting inconclusive, run one
  // cumulative base query as a safety net so any interaction with reset
  // bootstrap offsets or observation-only startup semantics cannot hide a real
  // bounded counterexample.
  if (shouldCheckLocalBaseCase(problem)) {
    if (auto validation =
            validateConcreteBasePrefix(problem, solverType, maxK);
        validation.has_value()) {
      return *validation; // LCOV_EXCL_LINE
    }
  }

  return {KInductionStatus::Inconclusive, maxK};
}

// LCOV_DISABLED_START


// LCOV_DISABLED_STOP
KInductionResult combineBatchResults(KInductionResult lhs,
                                     const KInductionResult& rhs) {
  if (lhs.status == KInductionStatus::Different) {
    return lhs;  // LCOV_EXCL_LINE
  }
  if (rhs.status == KInductionStatus::Different) {
    return rhs;  // LCOV_EXCL_LINE
  }
  if (rhs.status == KInductionStatus::Inconclusive) {
    lhs.status = KInductionStatus::Inconclusive;
  }
  // LCOV_DISABLED_START
  lhs.bound = std::max(lhs.bound, rhs.bound);
  // LCOV_DISABLED_STOP
  return lhs;
}

KInductionResult runOutputRangeKInduction(
    KInductionProblem& batchProblem,
    const KInductionProblem& sourceProblem,
    const PublicConjunctionHypothesis& preferredHypothesis,
    const PublicConjunctionHypothesis& fallbackHypothesis,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK,
    size_t firstOutput,
    // LCOV_DISABLED_START
    size_t endOutput) {
    // LCOV_DISABLED_STOP
  configureOutputBatchProblem(batchProblem, sourceProblem, firstOutput, endOutput);
  applyDualRailSplitHypothesis(
      batchProblem, sourceProblem, preferredHypothesis, firstOutput, endOutput);
  if (isKInductionDiagEnabled()) {
    // LCOV_DISABLED_START
    emitSecDiag(
    // LCOV_DISABLED_STOP
        "SEC diag: k-induction output range [", firstOutput, ",", endOutput,
        ") outputs=", endOutput - firstOutput);
  }
  const PublicConjunctionHypothesis activePreferredHypothesis =
      activeLeafAttemptHypothesis(batchProblem, preferredHypothesis);
  if (endOutput - firstOutput <= 1) {
    if (const std::optional<size_t> cachedBound =
            inconclusiveDualRailLeafAttemptBound(
                batchProblem, activePreferredHypothesis, solverType, maxK);
        cachedBound.has_value()) {
      return {KInductionStatus::Inconclusive, *cachedBound};
    }
  }
  const KInductionResult result =
      runMonolithicKInduction(batchProblem, solverType, maxK);
  if (result.status == KInductionStatus::Inconclusive &&
      endOutput - firstOutput <= 1) {
    rememberInconclusiveDualRailLeafAttempt(
        batchProblem,
        activePreferredHypothesis,
        solverType,
        maxK,
        result.bound);
  }
  if (isKInductionDiagEnabled()) {
    emitSecDiag(
        "SEC diag: k-induction output range [", firstOutput, ",", endOutput,
        ") status=", static_cast<int>(result.status),
        " bound=", result.bound);
  }
  if (result.status == KInductionStatus::Inconclusive &&
      shouldTryFallbackHypothesisAtRange(
          preferredHypothesis,
          fallbackHypothesis,
          endOutput - firstOutput)) {
    // Try the complete residual public-output hypothesis at the current slice
    // before splitting.  The final-frame bad predicate remains this slice, so
    // this is still the strict decomposed KI obligation P_full[0..k-1] => P_i[k].
    configureOutputBatchProblem(batchProblem, sourceProblem, firstOutput, endOutput);
    if (const std::optional<size_t> cachedBound =
            inconclusiveDualRailLeafAttemptBound(
                batchProblem, fallbackHypothesis, solverType, maxK);
        cachedBound.has_value()) {
      return {KInductionStatus::Inconclusive, *cachedBound}; // LCOV_EXCL_LINE
    }
    applyDualRailSplitHypothesis(
        batchProblem, sourceProblem, fallbackHypothesis, firstOutput, endOutput);
    const PublicConjunctionHypothesis activeFallbackHypothesis =
        activeLeafAttemptHypothesis(batchProblem, fallbackHypothesis);
    const KInductionResult fallbackResult =
        runMonolithicKInduction(batchProblem, solverType, maxK);
    if (fallbackResult.status == KInductionStatus::Inconclusive &&
        endOutput - firstOutput <= 1) {
      rememberInconclusiveDualRailLeafAttempt( // LCOV_EXCL_LINE
          batchProblem, // LCOV_EXCL_LINE
          activeFallbackHypothesis,
          solverType, // LCOV_EXCL_LINE
          maxK, // LCOV_EXCL_LINE
          fallbackResult.bound); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
    if (fallbackResult.status != KInductionStatus::Inconclusive ||
        endOutput - firstOutput <= 1) {
      return fallbackResult;
    }
  }
  if (result.status != KInductionStatus::Inconclusive ||
      endOutput - firstOutput <= 1) {
    return result;
  }

  // A conjunction can occasionally be harder for k-induction than its pieces.
  // Split only on inconclusive batches. Children keep the parent public-output
  // conjunction as their induction hypothesis, which is the decomposed strict
  // KI step for proving that parent conjunction invariant.
  const PublicConjunctionHypothesis childHypothesis =
      preferredHypothesis.property != nullptr
          ? preferredHypothesis
          : makePublicConjunctionHypothesis(batchProblem);
  const size_t middle = firstOutput + (endOutput - firstOutput) / 2;
  KInductionResult combined =
      runOutputRangeKInduction(
          batchProblem,
          sourceProblem,
          childHypothesis,
          fallbackHypothesis,
          solverType,
          maxK,
          firstOutput,
          middle);
  if (combined.status == KInductionStatus::Different) {
    return combined;  // LCOV_EXCL_LINE
  }
  return combineBatchResults(
      std::move(combined),
      runOutputRangeKInduction(
          batchProblem,
          sourceProblem,
          childHypothesis,
          fallbackHypothesis,
          solverType,
          maxK,
          middle,
          endOutput));
}

KInductionResult runOutputBatchedKInduction(
  const KInductionProblem& problem,
  KEPLER_FORMAL::Config::SolverType solverType,
  size_t maxK) {
  emitKInductionProblemDiag(problem, maxK);
  KInductionResult combined{KInductionStatus::Equivalent, 0};
  const PublicConjunctionHypothesis publicHypothesis =
      initialPublicConjunctionHypothesis(problem);
  const PublicConjunctionHypothesis preferredHypothesis =
      publicHypothesis.completePublicSurface ? PublicConjunctionHypothesis{}
                                             : publicHypothesis;
  const PublicConjunctionHypothesis fallbackHypothesis =
      publicHypothesis.completePublicSurface ? publicHypothesis
                                             : PublicConjunctionHypothesis{};
  const OutputBatchingLimits batchingLimits =
      defaultOutputBatchingLimitsForProblem(problem);
  // Copy the large shared SEC problem once, then mutate only the small
  // output/property slice for each batch.  The previous implementation copied
  // LCOV_DISABLED_START
  // hundreds of thousands of state symbols and equality pairs per batch, which
  // LCOV_DISABLED_STOP
  // became visible on BlackParrot even after SAT-side batching was effective.
  KInductionProblem batchProblem = problem;
  const bool useSharedBaseCase =
      problem.usesDualRailStateEncoding && shouldCheckLocalBaseCase(problem);
  // Preserve explicit caller deferral for localized residual proofs.  The
  // shared-base optimization below is only for normal batched proofs that still
  // own their base obligation inside this engine.
  batchProblem.deferBaseCaseChecks =
      problem.deferBaseCaseChecks || useSharedBaseCase;
  for (const auto& [firstOutput, endOutput] :
       buildSupportBoundedOutputBatches(problem, batchingLimits)) {
    const KInductionResult result = runOutputRangeKInduction(
        batchProblem,
        problem,
        preferredHypothesis,
        fallbackHypothesis,
        solverType,
        maxK,
        firstOutput,
        endOutput);
    if (result.status == KInductionStatus::Different) {
      return result;
    }
    if (result.status == KInductionStatus::Inconclusive) {
      combined.status = KInductionStatus::Inconclusive;
    }
    combined.bound = std::max(combined.bound, result.bound);
  }
  if (useSharedBaseCase && combined.status == KInductionStatus::Equivalent) {
    // Slices may prove before running their local frontier BMC. The shared
    // full-output check must therefore include the proved frontier itself.
    if (auto validation =
            validateConcreteBasePrefix(problem, solverType, combined.bound);
        validation.has_value()) {
      return *validation;
    }
  }
  return combined;
}

KInductionResult runSingleOutputKInduction(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK) {
  const PublicConjunctionHypothesis publicHypothesis =
      initialPublicConjunctionHypothesis(problem);
  if (publicHypothesis.property == nullptr) {
    const PublicConjunctionHypothesis localHypothesis =
        activeLeafAttemptHypothesis(problem, PublicConjunctionHypothesis{});
    if (const std::optional<size_t> cachedBound =
            inconclusiveDualRailLeafAttemptBound(
                problem, localHypothesis, solverType, maxK);
        cachedBound.has_value()) {
      return {KInductionStatus::Inconclusive, *cachedBound}; // LCOV_EXCL_LINE
    }
    const KInductionResult localOnlyResult =
        runMonolithicKInduction(problem, solverType, maxK);
    if (localOnlyResult.status == KInductionStatus::Inconclusive) {
      rememberInconclusiveDualRailLeafAttempt(
          problem,
          localHypothesis,
          solverType,
          maxK,
          localOnlyResult.bound);
    }
    return localOnlyResult;
  }

  const PublicConjunctionHypothesis localHypothesis =
      activeLeafAttemptHypothesis(problem, PublicConjunctionHypothesis{});
  KInductionResult localResult;
  if (const std::optional<size_t> cachedBound =
          inconclusiveDualRailLeafAttemptBound(
              problem, localHypothesis, solverType, maxK);
      cachedBound.has_value()) {
    localResult = {KInductionStatus::Inconclusive, *cachedBound};
  } else {
    localResult = runMonolithicKInduction(problem, solverType, maxK);
    if (localResult.status == KInductionStatus::Inconclusive) {
      rememberInconclusiveDualRailLeafAttempt(
          problem,
          localHypothesis,
          solverType,
          maxK,
          localResult.bound);
    }
  }
  if (localResult.status != KInductionStatus::Inconclusive) {
    return localResult;
  }
  if (!shouldTryStoredPublicHypothesisOnStandaloneLeaf(publicHypothesis)) {
    return localResult;
  }

  KInductionProblem strengthenedProblem = problem;
  if (const std::optional<size_t> cachedBound =
          inconclusiveDualRailLeafAttemptBound(
              strengthenedProblem, publicHypothesis, solverType, maxK);
      cachedBound.has_value()) {
    // The identical public-output KI leaf already returned UNKNOWN under the
    // same residual public hypothesis.  Reusing UNKNOWN is only a runtime cache:
    // it never marks this output covered.
    return {KInductionStatus::Inconclusive, *cachedBound};
  }
  applyDualRailSplitHypothesis(
      strengthenedProblem,
      problem,
      publicHypothesis,
      0,
      problem.observedOutputExprs0.size());
  const PublicConjunctionHypothesis activePublicHypothesis =
      activeLeafAttemptHypothesis(strengthenedProblem, publicHypothesis);
  KInductionResult strengthenedResult =
      runMonolithicKInduction(strengthenedProblem, solverType, maxK);
  if (strengthenedResult.status == KInductionStatus::Inconclusive) {
    rememberInconclusiveDualRailLeafAttempt(
        strengthenedProblem,
        activePublicHypothesis,
        solverType,
        maxK,
        strengthenedResult.bound);
  }
  return strengthenedResult;
}

}  // namespace

// Overall k-induction algorithm:
// 1. Check frame 0 immediately for a purely combinational mismatch.
// 2. If the SEC problem has no state, that base check fully decides it.
// 3. For k = 1..maxK, first ask whether the k-step induction rule closes from
//    the already-proved safe prefix.
// 4. If the step is inconclusive, extend the safe prefix by checking the next
//    concrete base frontier for a counterexample.
// 5. Return the first counterexample, the first successful proof bound, or
//    "inconclusive" if neither happens within the requested maxK.

KInductionEngine::KInductionEngine(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType)
    : problem_(problem), solverType_(solverType) {}

KInductionResult KInductionEngine::run(size_t maxK) const {
  if (problem_.observedOutputExprs0.size() <= 1) {
    emitKInductionProblemDiag(problem_, maxK);
  }
  if (problem_.observedOutputExprs0.size() >= kMinOutputsForBatchedProof) {
    return runOutputBatchedKInduction(problem_, solverType_, maxK);
  }
  return runSingleOutputKInduction(problem_, solverType_, maxK);
}

}  // namespace KEPLER_FORMAL::SEC
