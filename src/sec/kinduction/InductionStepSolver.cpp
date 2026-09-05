// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "kinduction/InductionStepSolver.h"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>
#include <unordered_set>

#include "../../config/Config.h"
#include "common/SecDiag.h"
#include "kinduction/SatEncoding.h"
#include "proof/TransitionExprResolver.h"

namespace KEPLER_FORMAL::SEC {

struct TransitionSupportCacheEntry {
  std::vector<size_t> targets;
  std::vector<size_t> stateSupport;
  std::vector<size_t> nonStateSupport;

  bool matches(const std::vector<size_t>& candidateTargets) const {
    return targets == candidateTargets;
  }

  size_t symbolCount() const {
    return targets.size() + stateSupport.size() + nonStateSupport.size();
  }

  void appendTo(std::unordered_set<size_t>& stateOutput,
                std::unordered_set<size_t>& nonStateOutput) const {
    stateOutput.insert(stateSupport.begin(), stateSupport.end());
    nonStateOutput.insert(nonStateSupport.begin(), nonStateSupport.end());
  }
};

struct InductionTransitionSupportCache {
  static constexpr size_t kMaxEntries = 8192;
  static constexpr size_t kMaxCachedSymbols = 4 * 1024 * 1024;

  size_t signature = 0;
  size_t entryCount = 0;
  size_t cachedSymbolCount = 0;
  std::unordered_map<size_t, std::vector<TransitionSupportCacheEntry>>
      entriesByHash;

  void resetForSignature(size_t newSignature) {
    if (signature == newSignature) {
      return;
    }
    signature = newSignature;
    entryCount = 0;
    cachedSymbolCount = 0;
    entriesByHash.clear();
  }

  const TransitionSupportCacheEntry* find(
      size_t targetHash,
      const std::vector<size_t>& targets) const {
    const auto bucketIt = entriesByHash.find(targetHash);
    if (bucketIt == entriesByHash.end()) {
      return nullptr;
    }
    for (const auto& entry : bucketIt->second) {
      if (entry.matches(targets)) {
        return &entry;
      }
    }
    return nullptr; // LCOV_EXCL_LINE
  }

  void store(size_t targetHash,
             const std::vector<size_t>& targets,
             const std::unordered_set<size_t>& stateSupport,
             const std::unordered_set<size_t>& nonStateSupport) {
    TransitionSupportCacheEntry entry{
        targets,
        sortedVector(stateSupport),
        sortedVector(nonStateSupport)};
    const size_t symbolCount = entry.symbolCount();
    if (entryCount >= kMaxEntries ||
        cachedSymbolCount + symbolCount > kMaxCachedSymbols) {
      return; // LCOV_EXCL_LINE
    }
    entriesByHash[targetHash].push_back(std::move(entry));
    ++entryCount;
    cachedSymbolCount += symbolCount;
  }

 private:
  static std::vector<size_t> sortedVector(
      const std::unordered_set<size_t>& symbols) {
    std::vector<size_t> sorted(symbols.begin(), symbols.end());
    std::sort(sorted.begin(), sorted.end());
    return sorted;
  }
};

namespace {

constexpr size_t kMaxSimplePathStateSymbols = 4096;
constexpr size_t kMinOriginalOutputsForCompactDualRailProfile = 64;
constexpr size_t kMinDeferredRailStateSymbolsForDirectProfile = 512;
constexpr size_t kMaxExactTransitionNodeCountHintTargets = 2048;
constexpr size_t kMaxTransitionNodeCountHint = 262144;
constexpr unsigned kDefaultDualRailLeafInductionDecisionLimit = 5000;
// Matches SATSolverWrapper's dual-rail no-preprocess cutoff.  Batched
// resource-limited probes still need the UNSAT-oriented proof profile, but
// samples show Kissat preprocessing can dominate these small-looking rail cones
// before the decision cap is observed.
constexpr size_t kDirectDualRailProofProfileSymbolFloor = 4096;
struct InductionCoi {
  std::vector<std::vector<size_t>> transitionTargetsByFrame;
  std::vector<std::vector<size_t>> transitionSupportByFrame;
  std::vector<size_t> relevantStateSymbols;
  std::vector<size_t> solverSymbols;
  std::unordered_set<size_t> solverSymbolSet;
};

KEPLER_FORMAL::Config::SolverType inductionStepSolverType(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType) {
  (void)problem;
  return solverType;
}

std::optional<unsigned> readUnsignedEnv(const char* name) {
  const char* value = std::getenv(name);
  if (value == nullptr || *value == '\0') {
    return std::nullopt; // LCOV_EXCL_LINE
  }
  char* end = nullptr;
  const unsigned long parsed = std::strtoul(value, &end, 10);
  if (end == value || *end != '\0' ||
      parsed > std::numeric_limits<unsigned>::max()) {
    // LCOV_EXCL_START
    return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  return static_cast<unsigned>(parsed);
}

bool isInductionStepCoiDiagEnabled() {
  return std::getenv("KEPLER_SEC_KI_COI_DIAG") != nullptr || isSecDiagEnabled();
}

size_t dualRailStateSymbolCount(const KInductionProblem& problem) {
  if (!problem.dualRailStatePairs.empty()) {
    return problem.dualRailStatePairs.size() * 2;
  }
  return problem.state0Symbols.size() + problem.state1Symbols.size();
}

size_t originalOutputCountForProblem(const KInductionProblem& problem) { // LCOV_EXCL_LINE
  return problem.originalObservedOutputCount == 0 // LCOV_EXCL_LINE
      ? problem.observedOutputExprs0.size() // LCOV_EXCL_LINE
      : problem.originalObservedOutputCount; // LCOV_EXCL_LINE
}

bool isLargeDeferredDualRailLeafSurface(const KInductionProblem& problem) {
  return problem.deferBaseCaseChecks &&
         dualRailStateSymbolCount(problem) >=
             kMinDeferredRailStateSymbolsForDirectProfile;
}

std::optional<unsigned> directInductionDecisionLimit(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (!problem.usesDualRailStateEncoding ||
      problem.observedOutputExprs0.size() > 1 ||
      solverType != KEPLER_FORMAL::Config::SolverType::KISSAT) {
    return std::nullopt;
  }
  if (const auto limit =
          readUnsignedEnv("KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT");
      limit.has_value()) {
    return limit;
  }
  if (originalOutputCountForProblem(problem) < // LCOV_EXCL_LINE
          kMinOriginalOutputsForCompactDualRailProfile && // LCOV_EXCL_LINE
      !isLargeDeferredDualRailLeafSurface(problem)) { // LCOV_EXCL_LINE
    // Compact dual-rail leaves can need the full strict KI search, and their
    // cones are small enough that an unbounded leaf is not the workflow wall.
    return std::nullopt; // LCOV_EXCL_LINE
  }
  return kDefaultDualRailLeafInductionDecisionLimit; // LCOV_EXCL_LINE
}

bool shouldUseDirectDualRailLimitedProofProfile(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::optional<unsigned>& kissatDecisionLimit) {
  return problem.usesDualRailStateEncoding &&
         solverType == KEPLER_FORMAL::Config::SolverType::KISSAT &&
         kissatDecisionLimit.has_value() && *kissatDecisionLimit > 0;
}

bool isWideDualRailResidualSurface(const KInductionProblem& problem) { // LCOV_EXCL_LINE
  return originalOutputCountForProblem(problem) >= // LCOV_EXCL_LINE
         kMinOriginalOutputsForCompactDualRailProfile;
}

size_t directDualRailProofProfileSymbols(const KInductionProblem& problem,
                                         size_t solverSymbols) {
  if (problem.deferBaseCaseChecks) {
    return std::max(solverSymbols, kDirectDualRailProofProfileSymbolFloor);
  }
  if (!isWideDualRailResidualSurface(problem) && // LCOV_EXCL_LINE
      !isLargeDeferredDualRailLeafSurface(problem)) { // LCOV_EXCL_LINE
    return solverSymbols; // LCOV_EXCL_LINE
  }
  return std::max(solverSymbols, kDirectDualRailProofProfileSymbolFloor); // LCOV_EXCL_LINE
}

bool shouldUseDirectCdclProfileForLimitedDualRailLeaf(
    const KInductionProblem& problem,
    const InductionCoi& coi) {
  (void)coi;
  // Deferred leaves are residual coverage attempts created after output
  // splitting.  Use the direct CDCL option mix only when the parent rail-state
  // surface is large enough to get a default decision cap; UNKNOWN still means
  // "not proved" and never becomes coverage.
  return isLargeDeferredDualRailLeafSurface(problem);
}

bool shouldUseDirectCdclProfileForCompactDualRailQuery(
    const KInductionProblem& problem) {
  // Small-to-medium dual-rail batches and their split leaves are real strict KI
  // obligations.  When the query is not using simple-path strengthening,
  // KISSAT's rephase/local-walk profile can dominate these rebuilt one-shot
  // proofs, so keep these rail surfaces on the direct CDCL option mix.
  return problem.usesDualRailStateEncoding &&
         dualRailStateSymbolCount(problem) <
             kMinDeferredRailStateSymbolsForDirectProfile;
}

size_t countTransitionTargets(
    const std::vector<std::vector<size_t>>& transitionTargetsByFrame) {
  size_t count = 0;
  for (const auto& targets : transitionTargetsByFrame) {
    count += targets.size();
  }
  return count;
}

size_t countFrameSupportSymbols(
    const std::vector<std::vector<size_t>>& supportByFrame) {
  size_t count = 0;
  for (const auto& support : supportByFrame) {
    count += support.size();
  }
  return count;
}

size_t boundedAddForReserveHint(size_t lhs, size_t rhs, size_t limit) {
  if (lhs >= limit || rhs >= limit - lhs) {
    return limit; // LCOV_EXCL_LINE
  }
  return lhs + rhs;
}

size_t transitionRelationNodeReserveHint(
    const TransitionExprResolver& transitionByState,
    const std::vector<size_t>& targets,
    const std::vector<size_t>& supportSymbols) {
  if (targets.size() > kMaxExactTransitionNodeCountHintTargets) {
    return supportSymbols.size(); // LCOV_EXCL_LINE
  }

  size_t exactNodeHint = 0;
  for (const auto stateSymbol : targets) {
    exactNodeHint = boundedAddForReserveHint(
        exactNodeHint,
        transitionByState.nodeCount(stateSymbol),
        kMaxTransitionNodeCountHint);
  }
  return std::max(supportSymbols.size(), exactNodeHint);
}

bool shouldAddSimplePathConstraint(const KInductionProblem& problem,
                                   const InductionCoi& coi) {
  if (!problem.usesDualRailStateEncoding) {
    return coi.relevantStateSymbols.size() <= kMaxSimplePathStateSymbols;
  }
  // Simple-path constraints are the standard loop-free strengthening used by
  // practical k-induction.  Apply them whenever the reduced dual-rail cone is
  // still below the existing state threshold: AES residual leaves need this
  // strict-KI strengthening, while very wide ASIC cones remain protected from
  // quadratic frame-pair clauses by the same bound.
  return coi.relevantStateSymbols.size() <= kMaxSimplePathStateSymbols;
}

std::unordered_set<size_t> buildStateSymbolSet(const KInductionProblem& problem) {
  std::unordered_set<size_t> stateSymbols;
  stateSymbols.reserve(problem.state0Symbols.size() + problem.state1Symbols.size());
  stateSymbols.insert(problem.state0Symbols.begin(), problem.state0Symbols.end());
  stateSymbols.insert(problem.state1Symbols.begin(), problem.state1Symbols.end());
  return stateSymbols;
}

std::unordered_map<size_t, size_t> buildPrimaryByComplementSymbol(
    const KInductionProblem& problem) {
  std::unordered_map<size_t, size_t> primaryByComplement;
  primaryByComplement.reserve(
      problem.complementedStatePairs0.size() +
      problem.complementedStatePairs1.size());
  for (const auto& [primarySymbol, complementedSymbol] :
       problem.complementedStatePairs0) {
    primaryByComplement.emplace(complementedSymbol, primarySymbol);
  }
  for (const auto& [primarySymbol, complementedSymbol] :
       problem.complementedStatePairs1) {
    primaryByComplement.emplace(complementedSymbol, primarySymbol);
  }
  return primaryByComplement;
}

std::vector<size_t> sortedSymbols(const std::unordered_set<size_t>& symbols) {
  std::vector<size_t> sorted(symbols.begin(), symbols.end());
  std::sort(sorted.begin(), sorted.end());
  return sorted;
}

size_t mixHashValue(size_t seed, size_t value) {
  seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
  return seed;
}

size_t hashSymbols(const std::vector<size_t>& symbols) {
  size_t seed = symbols.size();
  for (const auto symbol : symbols) {
    seed = mixHashValue(seed, symbol);
  }
  return seed;
}

size_t transitionSupportCacheSignature(const KInductionProblem& problem) {
  size_t seed = reinterpret_cast<size_t>(problem.lazyTransitions.get());
  seed = mixHashValue(seed, problem.transitions0.size());
  seed = mixHashValue(
      seed, reinterpret_cast<size_t>(problem.transitions0.data()));
  seed = mixHashValue(seed, problem.transitions1.size());
  seed = mixHashValue(
      seed, reinterpret_cast<size_t>(problem.transitions1.data()));
  seed = mixHashValue(seed, problem.state0Symbols.size());
  seed = mixHashValue(
      seed, reinterpret_cast<size_t>(problem.state0Symbols.data()));
  seed = mixHashValue(seed, problem.state1Symbols.size());
  seed = mixHashValue(
      seed, reinterpret_cast<size_t>(problem.state1Symbols.data()));
  return seed;
}

InductionTransitionSupportCache& transitionSupportCacheForProblem(
    const KInductionProblem& problem) {
  if (!problem.inductionTransitionSupportCache) {
    problem.inductionTransitionSupportCache =
        std::make_shared<InductionTransitionSupportCache>();
  }
  auto& cache = *problem.inductionTransitionSupportCache;
  cache.resetForSignature(transitionSupportCacheSignature(problem));
  return cache;
}

void collectTransitionSupportWithCache(
    const KInductionProblem& problem,
    const TransitionExprResolver& transitionByState,
    const std::vector<size_t>& targets,
    const std::unordered_set<size_t>& stateSymbols,
    std::unordered_set<size_t>& stateOutput,
    std::unordered_set<size_t>& nonStateOutput) {
  if (targets.empty()) {
    return;
  }

  auto& cache = transitionSupportCacheForProblem(problem);
  const size_t targetHash = hashSymbols(targets);
  if (const auto* cached = cache.find(targetHash, targets);
      cached != nullptr) {
    cached->appendTo(stateOutput, nonStateOutput);
    return;
  }

  std::unordered_set<size_t> stateSupport;
  std::unordered_set<size_t> nonStateSupport;
  stateSupport.reserve(targets.size());
  nonStateSupport.reserve(targets.size());
  transitionByState.collectSupportForTargets(
      targets,
      stateSymbols,
      stateSupport,
      nonStateSupport);
  stateOutput.insert(stateSupport.begin(), stateSupport.end());
  nonStateOutput.insert(nonStateSupport.begin(), nonStateSupport.end());
  cache.store(targetHash, targets, stateSupport, nonStateSupport);
}

std::vector<size_t> mergeTransitionSupportForEncoding(
    const std::unordered_set<size_t>& stateSupport,
    const std::unordered_set<size_t>& nonStateSupport) {
  std::unordered_set<size_t> merged;
  merged.reserve(stateSupport.size() + nonStateSupport.size());
  merged.insert(stateSupport.begin(), stateSupport.end());
  merged.insert(nonStateSupport.begin(), nonStateSupport.end());
  return sortedSymbols(merged);
}

void addFormulaStateSupport(BoolExpr* formula,
                            const std::unordered_set<size_t>& stateSymbols,
                            std::unordered_set<size_t>& output) {
  if (formula == nullptr) {
    // LCOV_EXCL_START
    return;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  for (const auto symbol : formula->getSupportVars()) {
    if (stateSymbols.find(symbol) != stateSymbols.end()) {
      output.insert(symbol);
    }
  }
}

void addFormulaSupport(BoolExpr* formula, std::unordered_set<size_t>& output) {
  if (formula == nullptr) {
    // LCOV_EXCL_START
    return;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  for (const auto symbol : formula->getSupportVars()) {
    if (symbol >= 2) {
      output.insert(symbol);
    }
  }
}

std::vector<size_t> expandTransitionTargets(
    const std::unordered_set<size_t>& requestedTargets,
    const TransitionExprResolver& transitionByState,
    const std::unordered_map<size_t, size_t>& primaryByComplement) {
  std::unordered_set<size_t> expanded;
  expanded.reserve(requestedTargets.size());
  for (const auto symbol : requestedTargets) {
    if (transitionByState.contains(symbol)) {
      expanded.insert(symbol);
      continue;
    }
    // LCOV_EXCL_START
    if (const auto primaryIt = primaryByComplement.find(symbol);  // LCOV_EXCL_LINE
        primaryIt != primaryByComplement.end() &&  // LCOV_EXCL_LINE
        transitionByState.contains(primaryIt->second)) {  // LCOV_EXCL_LINE
      expanded.insert(primaryIt->second);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  return sortedSymbols(expanded);
}

void addRelevantComplementPartners(
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    std::unordered_set<size_t>& solverSymbols) {
  for (const auto& [primarySymbol, complementedSymbol] : complementedStatePairs) {
    if (solverSymbols.find(primarySymbol) != solverSymbols.end() ||
        // LCOV_EXCL_START
        solverSymbols.find(complementedSymbol) != solverSymbols.end()) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      solverSymbols.insert(primarySymbol);
      solverSymbols.insert(complementedSymbol);
    }
  }
}

void addRelevantDualRailPartners(
    const std::vector<DualRailSymbolPair>& railPairs,
    std::unordered_set<size_t>& solverSymbols) {
  for (const auto& rails : railPairs) {
    if (solverSymbols.find(rails.mayBeOne) != solverSymbols.end() ||
        solverSymbols.find(rails.mayBeZero) != solverSymbols.end()) {  // LCOV_EXCL_LINE
      solverSymbols.insert(rails.mayBeOne);
      solverSymbols.insert(rails.mayBeZero);
    }
  }
}

void addPostBootstrapResetInputSymbols(
    const KInductionProblem& problem,
    std::unordered_set<size_t>& solverSymbols) {
  if (problem.resetBootstrapCycles == 0) {
    return;
  }
  for (const auto& [symbol, _] : problem.resetBootstrapInputs) {
    solverSymbols.insert(symbol);
  }
}

void closeStateEqualityDependencies(
    const std::vector<std::pair<size_t, size_t>>& equalityPairs,
    std::unordered_set<size_t>& stateSymbols) {
  bool changed = true;
  while (changed) {
    changed = false;
    for (const auto& [lhsSymbol, rhsSymbol] : equalityPairs) {
      const bool lhsNeeded = stateSymbols.find(lhsSymbol) != stateSymbols.end();
      const bool rhsNeeded = stateSymbols.find(rhsSymbol) != stateSymbols.end();
      if (!lhsNeeded && !rhsNeeded) {
        continue; // LCOV_EXCL_LINE
      }
      changed |= stateSymbols.insert(lhsSymbol).second;
      changed |= stateSymbols.insert(rhsSymbol).second;
    }
  }
}

void closeSameFrameStateEqualityDependencies(
    const KInductionProblem& problem,
    std::unordered_set<size_t>& stateSymbols) {
  closeStateEqualityDependencies(problem.sameFrameStateEqualityPairs0, stateSymbols);
  closeStateEqualityDependencies(problem.sameFrameStateEqualityPairs1, stateSymbols);
}

void addRelevantSameFrameStateEqualityPartners(
    const KInductionProblem& problem,
    std::unordered_set<size_t>& solverSymbols) {
  closeSameFrameStateEqualityDependencies(problem, solverSymbols);
}

InductionCoi buildInductionCoi(const KInductionProblem& problem,
                               BoolExpr* inductionProperty,
                               BoolExpr* inductionBad,
                               size_t k) {
  // Cone-of-influence reduction for the actual k-induction SAT problem:
  // start from the formulas that are asserted at each frame, then walk
  // backwards through only the transition equations needed to define those
  // state bits. Large ASICs can have hundreds of thousands of modeled memory
  // and flop bits; leaving unconstrained, irrelevant bits in the solver makes
  // Kissat spend time deciding variables that cannot affect the proof.
  const auto stateSymbols = buildStateSymbolSet(problem);
  const TransitionExprResolver transitionByState(problem);
  const auto primaryByComplement = buildPrimaryByComplementSymbol(problem);

  std::vector<std::unordered_set<size_t>> requiredStates(k + 1);
  std::unordered_set<size_t> transitionSupportSymbols;
  transitionSupportSymbols.reserve(1024);
  for (size_t frame = 0; frame < k; ++frame) {
    addFormulaStateSupport(inductionProperty, stateSymbols, requiredStates[frame]);
  }
  addFormulaStateSupport(inductionBad, stateSymbols, requiredStates[k]);

  std::vector<std::vector<size_t>> transitionTargetsByFrame(k);
  std::vector<std::vector<size_t>> transitionSupportByFrame(k);
  for (size_t frame = k; frame > 0; --frame) {
    closeSameFrameStateEqualityDependencies(problem, requiredStates[frame]);
    auto targets = expandTransitionTargets(
        requiredStates[frame],
        transitionByState,
        primaryByComplement);
    transitionTargetsByFrame[frame - 1] = targets;
    // Wide dual-rail buses share large transition cones.  Ask the resolver for
    // the whole frame target set at once so KI/IMC reuse the same DAG walk
    // while still collecting exactly the symbols needed by the proof.
    std::unordered_set<size_t> frameStateSupport;
    std::unordered_set<size_t> frameNonStateSupport;
    collectTransitionSupportWithCache(
        problem,
        transitionByState,
        targets,
        stateSymbols,
        frameStateSupport,
        frameNonStateSupport);
    requiredStates[frame - 1].insert(
        frameStateSupport.begin(), frameStateSupport.end());
    transitionSupportSymbols.insert(
        frameNonStateSupport.begin(), frameNonStateSupport.end());
    transitionSupportByFrame[frame - 1] =
        mergeTransitionSupportForEncoding(
            frameStateSupport, frameNonStateSupport);
  }
  closeSameFrameStateEqualityDependencies(problem, requiredStates[0]);

  std::unordered_set<size_t> solverSymbols;
  solverSymbols.reserve(1024);
  addFormulaSupport(inductionProperty, solverSymbols);
  addFormulaSupport(inductionBad, solverSymbols);

  std::unordered_set<size_t> relevantStateSymbols;
  for (const auto& frameStates : requiredStates) {
    relevantStateSymbols.insert(frameStates.begin(), frameStates.end());
    solverSymbols.insert(frameStates.begin(), frameStates.end());
  }
  for (const auto& targets : transitionTargetsByFrame) {
    for (const auto target : targets) {
      relevantStateSymbols.insert(target);
      solverSymbols.insert(target);
    }
  }
  solverSymbols.insert(
      transitionSupportSymbols.begin(), transitionSupportSymbols.end());
  addRelevantComplementPartners(problem.complementedStatePairs0, solverSymbols);
  addRelevantComplementPartners(problem.complementedStatePairs1, solverSymbols);
  addRelevantSameFrameStateEqualityPartners(problem, solverSymbols);
  addRelevantDualRailPartners(problem.dualRailStatePairs, solverSymbols);
  addPostBootstrapResetInputSymbols(problem, solverSymbols);

  InductionCoi coi;
  coi.transitionTargetsByFrame = std::move(transitionTargetsByFrame);
  coi.transitionSupportByFrame = std::move(transitionSupportByFrame);
  coi.relevantStateSymbols = sortedSymbols(relevantStateSymbols);
  coi.solverSymbols = sortedSymbols(solverSymbols);
  coi.solverSymbolSet.insert(coi.solverSymbols.begin(), coi.solverSymbols.end());
  return coi;
}

void emitInductionStepCoiDiag(const KInductionProblem& problem,
                              const InductionCoi& coi,
                              size_t k) {
  if (!isInductionStepCoiDiagEnabled()) {
    return;
  }
  emitSecDiag(
      "SEC diag: k-induction step coi k=", k,
      " solver_symbols=", coi.solverSymbols.size(),
      " transition_targets=", countTransitionTargets(coi.transitionTargetsByFrame),
      " transition_support=", countFrameSupportSymbols(coi.transitionSupportByFrame),
      " relevant_states=", coi.relevantStateSymbols.size());
}

void addComplementedStateRelations(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    const std::unordered_set<size_t>& solverSymbols,
    size_t numFrames) {
  for (size_t frame = 0; frame < numFrames; ++frame) {
    for (const auto& [primarySymbol, complementedSymbol] : complementedStatePairs) {
      if (solverSymbols.find(primarySymbol) == solverSymbols.end() ||
          solverSymbols.find(complementedSymbol) == solverSymbols.end()) {
        // LCOV_EXCL_START
        continue;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      addLiteralEquivalence(
          solver,
          variables.getLiteral(complementedSymbol, frame),
          -variables.getLiteral(primarySymbol, frame));
    }
  }
}

void addSameFrameStateEqualities(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<std::pair<size_t, size_t>>& equalityPairs,
    const std::unordered_set<size_t>& solverSymbols,
    size_t numFrames) {
  for (size_t frame = 0; frame < numFrames; ++frame) {
    for (const auto& [lhsSymbol, rhsSymbol] : equalityPairs) {
      if (solverSymbols.find(lhsSymbol) == solverSymbols.end() ||
          solverSymbols.find(rhsSymbol) == solverSymbols.end()) {
        continue;  // LCOV_EXCL_LINE
      }
      addLiteralEquivalence(
          solver,
          variables.getLiteral(lhsSymbol, frame),
          variables.getLiteral(rhsSymbol, frame));
    }
  }
}

void addSameFrameStateEqualities(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& solverSymbols,
    size_t numFrames) {
  addSameFrameStateEqualities(
      solver,
      variables,
      problem.sameFrameStateEqualityPairs0,
      solverSymbols,
      numFrames);
  addSameFrameStateEqualities(
      solver,
      variables,
      problem.sameFrameStateEqualityPairs1,
      solverSymbols,
      numFrames);
}

void addTransitionRelation(SATSolverWrapper& solver,
                           const FrameVariableStore& variables,
                           const TransitionExprResolver& transitionByState,
                           const std::vector<size_t>& targets,
                           const std::vector<size_t>& supportSymbols,
                           size_t frame) {
  // The targets of a single frame often share most of their input cone. Keep a
  // frame-local encoder so common BoolExpr nodes are Tseitin-encoded once and
  // reused by every next-state equality in this frame.  Use the exact
  // transition read-set computed by the strict KI COI pass instead of every
  // solver symbol; this emits the same transition clauses while avoiding a wide
  // per-frame leaf table for unrelated state and output cones.
  const size_t reserveHint =
      transitionRelationNodeReserveHint(transitionByState, targets, supportSymbols);
  if (isInductionStepCoiDiagEnabled() && reserveHint > supportSymbols.size()) {
    emitSecDiag(
        "SEC diag: k-induction transition node reserve targets=",
        targets.size(),
        " support=", supportSymbols.size(),
        " expected_nodes=", reserveHint);
  }
  FrameFormulaEncoder encoder(
      solver,
      variables.makeLeafLits(frame, supportSymbols),
      reserveHint);
  for (const auto stateSymbol : targets) {
    BoolExpr* expr = transitionByState.at(stateSymbol);
    addLiteralEquivalence(
        solver,
        variables.getLiteral(stateSymbol, frame + 1),
        encoder.encode(expr));
  }
}

void addDualRailStateValidity(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<DualRailSymbolPair>& railPairs,
    const std::unordered_set<size_t>& solverSymbols,
    size_t numFrames) {
  for (size_t frame = 0; frame < numFrames; ++frame) {
    for (const auto& rails : railPairs) {
      if (solverSymbols.find(rails.mayBeOne) == solverSymbols.end() ||
          solverSymbols.find(rails.mayBeZero) == solverSymbols.end()) {
        continue;  // LCOV_EXCL_LINE
      }
      // Dual-rail encodes a non-empty possible-value set.  The empty set
      // (may1=0, may0=0) is not a legal ternary value and must not be available
      // to induction as a synthetic predecessor.
      solver.addClause({
          variables.getLiteral(rails.mayBeOne, frame),
          variables.getLiteral(rails.mayBeZero, frame)});
    }
  }
}

void addPostBootstrapResetInputConstraints(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    size_t numFrames) {
  if (problem.resetBootstrapCycles == 0) {
    return;
  }

  for (const auto& [symbol, assertedValue] : problem.resetBootstrapInputs) {
    if (!variables.hasSymbol(symbol)) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    for (size_t frame = 0; frame < numFrames; ++frame) {
      // The induction query starts at the post-bootstrap frontier.  Match the
      // base-case/PDR environment by keeping reset controls deasserted
      // throughout this window instead of proving across arbitrary reset
      // reassertion.
      solver.addClause(
          {assertedValue ? -variables.getLiteral(symbol, frame)
                         // LCOV_EXCL_START
                         : variables.getLiteral(symbol, frame)});  // LCOV_EXCL_LINE
                         // LCOV_EXCL_STOP
    }
  }
}

}  // namespace

InductionProofStatus proveByInductionStatus(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k,
    std::optional<unsigned> kissatDecisionLimit) {
  const bool hasExplicitInductionInvariant = problem.inductionProperty != nullptr;
  BoolExpr* inductionProperty =
      hasExplicitInductionInvariant ? problem.inductionProperty : problem.property;
  BoolExpr* inductionBad =
      problem.inductionBad != nullptr ? problem.inductionBad : problem.bad;
  const InductionCoi coi = buildInductionCoi(
      problem,
      inductionProperty,
      inductionBad,
      k);
  emitInductionStepCoiDiag(problem, coi, k);
  const auto inductionPropertySupport = inductionProperty->getSupportVars();
  const auto inductionBadSupport = inductionBad->getSupportVars();
  const TransitionExprResolver transitionByState(problem);
  const bool addSimplePath = shouldAddSimplePathConstraint(problem, coi);

  const auto localSolverType = inductionStepSolverType(problem, solverType);
  SATSolverWrapper solver(localSolverType);
  if (shouldUseDirectDualRailLimitedProofProfile(
          problem, localSolverType, kissatDecisionLimit)) {
    // Resource-limited dual-rail induction queries are speculative KI attempts:
    // UNKNOWN just asks the caller to split or try the next k. Keep the UNSAT
    // proof profile, and force the existing no-preprocess branch only for wide
    // residual surfaces where Kissat preprocessing can dominate before the
    // decision cap is seen.
    const size_t profileSymbols =
        directDualRailProofProfileSymbols(problem, coi.solverSymbols.size());
    solver.configureForSecDualRailConeProof(
        profileSymbols);
    const bool useDirectCdcl =
        shouldUseDirectCdclProfileForLimitedDualRailLeaf(problem, coi);
    if (useDirectCdcl) {
      solver.configureForSecLocalBooleanCheck(profileSymbols);
    }
    if (isInductionStepCoiDiagEnabled()) {
      emitSecDiag(
          "SEC diag: k-induction direct dual-rail capped proof profile outputs=",
          problem.observedOutputExprs0.size(),
          " solver_symbols=", coi.solverSymbols.size(),
          " profile_symbols=", profileSymbols,
          " direct_cdcl=",
          useDirectCdcl ? 1 : 0);
    }
  } else if (problem.usesDualRailStateEncoding) {
    // Deferred residual leaves are rebuilt one output at a time.  Even when a
    // small leaf is unbounded so it can genuinely prove (GCD), Kissat's
    // standalone preprocessing/sweeping can dominate these one-shot strict KI
    // formulas.  Use the dual-rail direct profile without adding any decision
    // cap; SAT/UNSAT/UNKNOWN semantics are unchanged.
    const size_t profileSymbols =
        problem.deferBaseCaseChecks
            ? directDualRailProofProfileSymbols(problem, coi.solverSymbols.size())
            : coi.solverSymbols.size();
    solver.configureForSecDualRailConeProof(profileSymbols);
    if (shouldUseDirectCdclProfileForCompactDualRailQuery(problem) &&
        !addSimplePath) {
      solver.configureForSecLocalBooleanCheck(profileSymbols); // LCOV_EXCL_LINE
      if (isInductionStepCoiDiagEnabled()) { // LCOV_EXCL_LINE
        emitSecDiag( // LCOV_EXCL_LINE
            "SEC diag: k-induction compact dual-rail direct profile outputs=",
            problem.observedOutputExprs0.size(), // LCOV_EXCL_LINE
            " solver_symbols=", coi.solverSymbols.size(), // LCOV_EXCL_LINE
            " profile_symbols=", profileSymbols,
            " direct_cdcl=1");
      } // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
  } else {
    solver.configureForSecConeProof(coi.solverSymbols.size());
  }
  FrameVariableStore variables(solver, coi.solverSymbols, k + 1);
  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs0, coi.solverSymbolSet, k + 1);
  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs1, coi.solverSymbolSet, k + 1);
  addSameFrameStateEqualities(
      solver, variables, problem, coi.solverSymbolSet, k + 1);
  addDualRailStateValidity(
      solver,
      variables,
      problem.dualRailStatePairs,
      coi.solverSymbolSet,
      k + 1);
  addPostBootstrapResetInputConstraints(solver, variables, problem, k + 1);

  for (size_t frame = 0; frame < k; ++frame) {
    addTransitionRelation(
        solver,
        variables,
        transitionByState,
        coi.transitionTargetsByFrame[frame],
        coi.transitionSupportByFrame[frame],
        frame);
  }

  for (size_t frame = 0; frame < k; ++frame) {
    FrameFormulaEncoder encoder(
        solver, variables.makeLeafLits(frame, inductionPropertySupport));
    solver.addClause({encoder.encode(inductionProperty)});
  }
  if (addSimplePath) {
    // The simple-path refinement is a completeness aid, not a soundness
    // requirement for classic k-induction. On large gate-level SEC problems it
    // creates one XOR per state bit per frame-pair, which can dominate the SAT
    // instance and drown the actual proof. Keep it for small dual-rail cones
    // where strict KI needs loop-free paths to avoid unreachable rail states.
    if (isInductionStepCoiDiagEnabled()) {
      emitSecDiag(
          "SEC diag: k-induction simple path states=",
          coi.relevantStateSymbols.size());
    }
    addSimplePathConstraint(solver, variables, coi.relevantStateSymbols, k + 1);
  }

  FrameFormulaEncoder lastFrameEncoder(
      solver, variables.makeLeafLits(k, inductionBadSupport));
  solver.addClause({lastFrameEncoder.encode(inductionBad)});
  SATSolverWrapper::SolveStatus solveStatus;
  if (localSolverType == KEPLER_FORMAL::Config::SolverType::KISSAT &&
      kissatDecisionLimit.has_value()) {
    solveStatus = solver.solveWithKissatResourceLimits(
        std::numeric_limits<unsigned>::max(), *kissatDecisionLimit);
  } else if (localSolverType == KEPLER_FORMAL::Config::SolverType::CADICAL &&
             // LCOV_EXCL_START
             kissatDecisionLimit.has_value()) {  // LCOV_EXCL_LINE
    solveStatus = solver.solveWithAssumptionsStatus(  // LCOV_EXCL_LINE
        {},  // LCOV_EXCL_LINE
        *kissatDecisionLimit,  // LCOV_EXCL_LINE
        *kissatDecisionLimit);  // LCOV_EXCL_LINE
  } else {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    solveStatus = solver.solveStatus();
  }
  switch (solveStatus) {
    case SATSolverWrapper::SolveStatus::Unsat:
      return InductionProofStatus::Proved;
    case SATSolverWrapper::SolveStatus::Sat:
      return InductionProofStatus::NotProved;
    case SATSolverWrapper::SolveStatus::Unknown:
      // LCOV_EXCL_START
      return InductionProofStatus::Unknown;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return InductionProofStatus::Unknown;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

bool provesByInduction(const KInductionProblem& problem,
                       KEPLER_FORMAL::Config::SolverType solverType,
                       size_t k) {
  return proveByInductionStatus(
             problem,
             solverType,
             k,
             directInductionDecisionLimit(problem, solverType)) ==
         InductionProofStatus::Proved;
}

}  // namespace KEPLER_FORMAL::SEC
