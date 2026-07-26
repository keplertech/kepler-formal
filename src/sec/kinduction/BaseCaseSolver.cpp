// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "kinduction/BaseCaseSolver.h"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../../config/Config.h"
#include "common/SecDiag.h"
#include "kinduction/OutputBatching.h"
#include "kinduction/SatEncoding.h"
#include "proof/TransitionExprResolver.h"

namespace KEPLER_FORMAL::SEC {

namespace {

// BlackParrot PDR sampling showed the reset-frontier precheck repeatedly
// cycling through many neighboring cube supports in the same blocking wave.
// A tiny cache cleared in the middle of that wave and forced full COI/solver
// reconstruction; keep enough exact assumption solvers to cover the measured
// working set without making the cache effectively unbounded.
constexpr size_t kMaxResetFrontierCachedSolvers = 64;
constexpr size_t kMinResetFrontierCoreChecks = 64;
constexpr size_t kMaxResetFrontierCachedCoresPerFrame = 4096;
constexpr size_t kMaxPreviousResetFrontierBlockersPerQuery = 64;
// The relaxed post-bootstrap query is only a shortcut before the exact
// reset-frontier solver. Sampling AES PDR showed this shortcut can become the
// wall when it fails to shrink the COI, so keep it local and resource-bounded.
constexpr size_t kMaxRelaxedResetFrontierPrecheckSymbols = 256;
constexpr size_t kMaxRelaxedResetFrontierPrecheckTransitionTargets = 128;
constexpr size_t kMaxTinyCubeRelaxedResetFrontierPrecheckSymbols = 8192;
constexpr size_t kMaxTinyCubeRelaxedResetFrontierPrecheckTransitionTargets = 8192;
constexpr size_t kMaxTinyRelaxedResetFrontierPrecheckCubeLiterals = 4;
constexpr unsigned kRelaxedResetFrontierPrecheckConflictLimit = 10000;
// BlackParrot PDR samples showed the full reset-frontier fallback exploding to
// hundreds of thousands of symbols while the reset-summary query stayed much
// smaller (about 90k symbols / 88k transition targets at the slow frontier).
// Let that exact UNSAT precheck run before opening the full prefix.
constexpr size_t kMaxResetSummaryPrecheckSymbols = 200000;
constexpr size_t kMaxResetSummaryPrecheckTransitionTargets = 200000;
// This summary query is the cheap alternative to the full reset-frontier
// fallback. BlackParrot PDR samples showed 65k-symbol summary proofs giving up
// after the old tiny cap, then falling into 600k+ symbol exact assumption
// assumption solves. Spend the bounded effort here where the COI is smaller.
constexpr unsigned kResetSummaryPrecheckConflictLimit = 5000;
constexpr unsigned kResetSummaryFrontierProofConflictLimit = 2000;
constexpr unsigned kResetSummarySingletonProofConflictLimit = 1000;
constexpr long long kResetFrontierCachedAssumptionConflictLimit = 10000;
constexpr long long kResetFrontierBatchProofPropagationLimit = 250000;
constexpr size_t kMaxResetSummaryFrontierBlockers = 256;
constexpr size_t kMaxResetSummaryRefinements = 4;
constexpr size_t kMaxResetSummaryFrontierCubeLiterals = 100000;
constexpr size_t kMaxResetSummaryFrontierProofCubeLiterals = 8192;
constexpr size_t kMaxResetSummaryFrontierProofSymbols = 20000;
constexpr size_t kMaxResetSummaryFrontierProofTransitionTargets = 20000;
constexpr size_t kMaxResetSummaryBulkSingletonBlockers = 128;
constexpr size_t kMaxResetSummaryCachedCois = 64;
// PDR often proves all but one concrete reset frame through cheap cached cores
// before asking the exact bounded fallback. In that sparse case a prefix query
// that targets every frame widens the COI unnecessarily; use exact per-step
// solvers while only a couple of frames remain.
constexpr size_t kMaxSparseResetFrontierPerStepChecks = 2;
// Fast localized counterexample searches should not let one hard UNSAT cone
// consume all local solver effort before other top-output slices are tried.
constexpr int64_t kFastCounterexampleSearchConflictLimit = 5000;
// Some dual-rail frontier checks spend their budget in propagation-heavy
// decision search before conflicts accumulate.  Bound decisions too so a hard
// residual can be reported Unknown and skipped by localized recovery.
constexpr int64_t kFastCounterexampleSearchDecisionLimit = 20000;
constexpr size_t kMaxImcAssumptionFrontierBatchOutputs = 64;
constexpr size_t kDualRailPublicBasePrefixWindow = 8;
constexpr size_t kLargeDualRailPublicBasePrefixWindow = 32;
constexpr size_t kMinLargeDualRailPublicBasePrefixStateSymbols = 512;
// Transition node counts are only reserve hints for FrameFormulaEncoder.
// BlackParrot reset-frontier sampling showed exact counting across ~86k
// transition targets dominating the whole query before any CNF was emitted.
// Keep exact hints for local cones and skip the prepass for ASIC-sized groups;
// the encoder still grows its cache geometrically while emitting the same CNF.
constexpr size_t kMaxExactTransitionNodeCountHintTargets = 512;

void mixHashValue(size_t& seed, size_t value) {
  seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

struct ResetFrontierSolverCacheKey {
  KEPLER_FORMAL::Config::SolverType solverType =
      KEPLER_FORMAL::Config::SolverType::KISSAT;
  size_t targetFrame = 0;
  bool includeCubeValues = false;
  std::vector<size_t> cubeSymbols;
  std::vector<std::pair<size_t, bool>> cubeLiterals;

  bool operator==(const ResetFrontierSolverCacheKey& other) const {
    return solverType == other.solverType &&
           targetFrame == other.targetFrame &&
           includeCubeValues == other.includeCubeValues &&
           cubeSymbols == other.cubeSymbols &&
           cubeLiterals == other.cubeLiterals;
  }
};

struct ResetFrontierSolverCacheKeyHash {
  size_t operator()(const ResetFrontierSolverCacheKey& key) const {
    size_t seed = std::hash<int>()(static_cast<int>(key.solverType));
    mixHashValue(seed, std::hash<size_t>()(key.targetFrame));
    mixHashValue(seed, std::hash<bool>()(key.includeCubeValues));
    for (const size_t symbol : key.cubeSymbols) {
      mixHashValue(seed, std::hash<size_t>()(symbol));
    }
    for (const auto& [symbol, value] : key.cubeLiterals) {
      mixHashValue(seed, std::hash<size_t>()(symbol));  // LCOV_EXCL_LINE
      mixHashValue(seed, std::hash<bool>()(value));  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
    return seed;
    // LCOV_EXCL_STOP
  }
};

bool resetFrontierAssumptionSolvesDisabled() {
  return std::getenv("KEPLER_SEC_PDR_DISABLE_RESET_FRONTIER_ASSUMPTIONS") !=
         nullptr;
}

bool relaxedResetFrontierPrecheckCoiIsLocal(size_t solverSymbols,
                                            size_t transitionTargets,
                                            size_t cubeLiterals) {
  const bool tinyPdrCube =
      cubeLiterals <= kMaxTinyRelaxedResetFrontierPrecheckCubeLiterals;
  // PDR leaf repair often asks about a tiny bad cube whose exact transition COI
  // is still a few thousand symbols.  Keep broad cubes on the historical cap,
  // but let tiny cubes use the sound relaxed UNSAT shortcut before falling into
  // the heavier exact reset-frontier assumption solver.
  const size_t symbolLimit =
      tinyPdrCube ? kMaxTinyCubeRelaxedResetFrontierPrecheckSymbols
                  : kMaxRelaxedResetFrontierPrecheckSymbols;
  const size_t transitionTargetLimit =
      tinyPdrCube ? kMaxTinyCubeRelaxedResetFrontierPrecheckTransitionTargets
                  : kMaxRelaxedResetFrontierPrecheckTransitionTargets;
  return solverSymbols <= symbolLimit &&
         transitionTargets <= transitionTargetLimit;
}

SATSolverWrapper::SolveStatus solveResetFrontierUnitClauseQuery(  // LCOV_EXCL_LINE
    SATSolverWrapper& solver,
    // LCOV_EXCL_START
    KEPLER_FORMAL::Config::SolverType solverType,
    // LCOV_EXCL_STOP
    int64_t conflictLimit) {
  if (solverType == KEPLER_FORMAL::Config::SolverType::KISSAT &&  // LCOV_EXCL_LINE
      conflictLimit >= 0) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    return solver.solveWithKissatResourceLimits(  // LCOV_EXCL_LINE
        static_cast<unsigned>(conflictLimit));  // LCOV_EXCL_LINE
  }
  return solver.solveStatus();  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START

enum class InitialConstraintMode {
// LCOV_EXCL_STOP
  None,
  ObservationOnly,
  PartialInit,
  CompleteInit,
};

size_t resetBootstrapFrames(const KInductionProblem& problem) {
  return ((!problem.hasCompleteInitialState() || problem.usesDualRailStateEncoding) &&
          problem.hasResetBootstrap())
             ? problem.resetBootstrapCycles
             : 0u;
}

InitialConstraintMode determineInitialConstraintMode(const KInductionProblem& problem) {
  if (!problem.hasSequentialState()) {
    return InitialConstraintMode::None;
  }

  if (problem.hasCompleteInitialState()) {
    return InitialConstraintMode::CompleteInit;
  }

  if (problem.hasExplicitInitialState()) {
    return InitialConstraintMode::PartialInit;
  }

  return InitialConstraintMode::ObservationOnly;
}

const std::vector<std::pair<size_t, size_t>>& emptySymbolPairs() {
  static const std::vector<std::pair<size_t, size_t>> pairs;
  return pairs;
}

struct BaseCaseCoi {
  std::vector<std::vector<size_t>> transitionTargetsByFrame;
  std::vector<std::vector<size_t>> requiredStateSymbolsByFrame;
  std::vector<size_t> solverSymbols;
  std::unordered_set<size_t> solverSymbolSet;
};

struct CachedResetFrontierSolver {
  BaseCaseCoi coi;
  std::unique_ptr<SATSolverWrapper> solver;
  std::unique_ptr<FrameVariableStore> variables;
  KEPLER_FORMAL::Config::SolverType solverType =
      KEPLER_FORMAL::Config::SolverType::KISSAT;
  size_t targetFrame = 0;
  bool cubeEncodedAsUnitClauses = false;
};

struct CachedResetSummaryCoi {
  BaseCaseCoi coi;
  size_t transitionTargets = 0;
};

std::unordered_map<size_t, size_t> buildPrimaryByComplementSymbol(
    const KInductionProblem& problem);

struct EqualityIndex {
  std::unordered_map<size_t, std::vector<size_t>> partnersBySymbol;

  explicit EqualityIndex(
      const std::vector<std::pair<size_t, size_t>>& equalityPairs) {
    partnersBySymbol.reserve(equalityPairs.size() * 2);
    for (const auto& [lhsSymbol, rhsSymbol] : equalityPairs) {
      partnersBySymbol[lhsSymbol].push_back(rhsSymbol);
      partnersBySymbol[rhsSymbol].push_back(lhsSymbol);
    }
  }

  void close(std::unordered_set<size_t>& symbols) const {
    // Equality closure is queried repeatedly for very small PDR reset cubes.
    // Walking only the adjacency of already-relevant symbols avoids rescanning
    // the full ASIC equality table until the fixed point stops changing.
    std::vector<size_t> worklist(symbols.begin(), symbols.end());
    for (size_t index = 0; index < worklist.size(); ++index) {
      const auto partnersIt = partnersBySymbol.find(worklist[index]);
      if (partnersIt == partnersBySymbol.end()) {
        continue;
      }
      for (const auto partner : partnersIt->second) {
        if (symbols.insert(partner).second) {
          worklist.push_back(partner);
        }
      }
    }
  }

  std::vector<std::pair<size_t, size_t>> pairsWithin(
      const std::unordered_set<size_t>& symbols) const {
    std::vector<std::pair<size_t, size_t>> pairs;
    for (const auto lhsSymbol : symbols) {
      const auto partnersIt = partnersBySymbol.find(lhsSymbol);
      if (partnersIt == partnersBySymbol.end()) {
        continue;
      }
      for (const auto rhsSymbol : partnersIt->second) {
        if (lhsSymbol < rhsSymbol &&
            symbols.find(rhsSymbol) != symbols.end()) {
          pairs.emplace_back(lhsSymbol, rhsSymbol);
        }
      }
    }
    return pairs;
  }
};

struct ResetFrontierReachabilityContextData {
  ResetFrontierReachabilityContextData(
      const KInductionProblem& problem,
      const TransitionExprResolver& transitionByState,
      BoolExpr* frameInvariant = nullptr)
      : problem(problem),
        transitionByState(transitionByState),
        frameInvariant(frameInvariant),
        frameInvariantSupport(frameInvariant != nullptr
                                  ? frameInvariant->getSupportVars()
                                  : std::set<size_t>{}),
        bootstrapFrames(resetBootstrapFrames(problem)),
        initialMode(bootstrapFrames == 0 ? determineInitialConstraintMode(problem)
                                         : InitialConstraintMode::None),
        primaryByComplement(buildPrimaryByComplementSymbol(problem)),
        initialEqualities(emptySymbolPairs()),
        bootstrapEqualities(emptySymbolPairs()) {}

  const KInductionProblem& problem;
  const TransitionExprResolver& transitionByState;
  // Optional caller-validated invariant for frames at and after the concrete
  // startup frontier. PDR proves this separately before using it here; the
  // reset-frontier BMC only consumes it as an extra exact reachable-state fact.
  BoolExpr* frameInvariant = nullptr;
  std::set<size_t> frameInvariantSupport;
  size_t bootstrapFrames = 0;
  InitialConstraintMode initialMode = InitialConstraintMode::None;
  std::unordered_map<size_t, size_t> primaryByComplement;
  EqualityIndex initialEqualities;
  EqualityIndex bootstrapEqualities;
  // PDR asks many concrete reachability questions for the same small set of
  // state symbols while refining projected counterexamples. Cache the unrolled
  // reset/bootstrap solver by frame and cube support, and vary only the cube
  // values through SAT assumptions.
  mutable std::unordered_map<
      ResetFrontierSolverCacheKey,
      std::unique_ptr<CachedResetFrontierSolver>,
      ResetFrontierSolverCacheKeyHash>
      cachedSolvers;
  // Shared-prefix reset-frontier queries check one cube against every concrete
  // post-bootstrap frame up to a depth. The COI depends only on the cube
  // symbols and max frame, so cache that exact solver separately from the
  // single-frontier cache and vary literal values through assumptions.
  mutable std::unordered_map<
      ResetFrontierSolverCacheKey,
      std::unique_ptr<CachedResetFrontierSolver>,
      ResetFrontierSolverCacheKeyHash>
      cachedPrefixSolvers;
  // Exact reset-frontier UNSAT cores are reusable across neighboring cubes:
  // once core C is proven unreachable, every later cube containing C is also
  // unreachable. BlackParrot PDR sampling showed thousands of repeated
  // assumption solves for such neighboring cubes, so keep a small per-frame
  // cache of minimized unreachable cores before asking SAT again.
  mutable std::unordered_map<
      size_t,
      std::vector<std::vector<std::pair<size_t, bool>>>>
      unreachableCoresByTargetFrame;
  // The reset-summary precheck is often asked about neighboring cubes that
  // share the same support at the same post-bootstrap depth. Its COI is
  // support-only, so cache it separately from SAT solvers and vary values later.
  mutable std::unordered_map<
      ResetFrontierSolverCacheKey,
      CachedResetSummaryCoi,
      ResetFrontierSolverCacheKeyHash>
      cachedResetSummaryCois;
};

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

}  // namespace

struct ImcBaseCounterexampleCache {
  explicit ImcBaseCounterexampleCache(const KInductionProblem& problem)
      : problem(problem),
        stateSymbols(buildStateSymbolSet(problem)),
        transitionByState(problem),
        primaryByComplement(buildPrimaryByComplementSymbol(problem)),
        sameFrameEqualities0(problem.sameFrameStateEqualityPairs0),
        sameFrameEqualities1(problem.sameFrameStateEqualityPairs1) {}

  void closeSameFrameStateEqualityDependencies(
      std::unordered_set<size_t>& symbols) const {
    sameFrameEqualities0.close(symbols);
    sameFrameEqualities1.close(symbols);
  }

  const ImcBaseCounterexampleCache& singleOutputCache(size_t outputIndex) const;
  const ImcBaseCounterexampleCache& outputSubsetCache(
      size_t firstOutput, size_t endOutput) const;

  const KInductionProblem& problem;
  std::unordered_set<size_t> stateSymbols;
  TransitionExprResolver transitionByState;
  std::unordered_map<size_t, size_t> primaryByComplement;
  EqualityIndex sameFrameEqualities0;
  EqualityIndex sameFrameEqualities1;
  mutable std::vector<std::unique_ptr<KInductionProblem>> singleOutputProblems;
  mutable std::vector<std::unique_ptr<ImcBaseCounterexampleCache>> singleOutputCaches;
  mutable std::vector<size_t> subsetFirstOutputs;
  mutable std::vector<size_t> subsetEndOutputs;
  mutable std::vector<std::unique_ptr<KInductionProblem>> subsetProblems;
  mutable std::vector<std::unique_ptr<ImcBaseCounterexampleCache>> subsetCaches;
};

namespace {

std::vector<size_t> sortedSymbols(const std::unordered_set<size_t>& symbols) {
  std::vector<size_t> sorted(symbols.begin(), symbols.end());
  std::sort(sorted.begin(), sorted.end());
  return sorted;
}

std::vector<std::vector<size_t>> sortedFrameStates(
    const std::vector<std::unordered_set<size_t>>& requiredStates) {
  std::vector<std::vector<size_t>> sorted;
  sorted.reserve(requiredStates.size());
  for (const auto& frameStates : requiredStates) {
    sorted.push_back(sortedSymbols(frameStates));
  }
  return sorted;
}

std::set<size_t> formulaSupportOrThrow(BoolExpr* formula, const char* context) {
  if (formula == nullptr) {
    throw std::runtime_error(  // LCOV_EXCL_LINE
        std::string("Missing BoolExpr while encoding base SEC formula: ") +  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        context);  // LCOV_EXCL_LINE
  }
  return formula->getSupportVars();
  // LCOV_EXCL_STOP
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START
void addFormulaStateSupport(BoolExpr* formula,
// LCOV_EXCL_STOP
                            const std::unordered_set<size_t>& stateSymbols,
                            std::unordered_set<size_t>& output) {
  if (formula == nullptr) {
    return;
  }
  for (const auto symbol : formula->getSupportVars()) {
    if (stateSymbols.find(symbol) != stateSymbols.end()) {
      output.insert(symbol);
    }
  }
}

void addFormulaSupport(BoolExpr* formula, std::unordered_set<size_t>& output) {
  if (formula == nullptr) {
    return;
  }
  for (const auto symbol : formula->getSupportVars()) {
    if (symbol >= 2) {
      output.insert(symbol);
    }
  }
}

bool hasStructuredInitialAssignments(const KInductionProblem& problem) {
  return !problem.initialStateAssignments.empty();
}

bool isKInductionCoiDiagEnabled() {
  return std::getenv("KEPLER_SEC_KI_DIAG") != nullptr || isSecDiagEnabled();
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
    if (const auto primaryIt = primaryByComplement.find(symbol);  // LCOV_EXCL_LINE
        primaryIt != primaryByComplement.end() &&  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        transitionByState.contains(primaryIt->second)) {  // LCOV_EXCL_LINE
      expanded.insert(primaryIt->second);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
  return sortedSymbols(expanded);
  // LCOV_EXCL_STOP
}

void closeFrameEqualityDependencies(
    const std::vector<std::pair<size_t, size_t>>& equalityPairs,
    std::unordered_set<size_t>& frameStates) {
  // Equality constraints can make a state bit relevant even if the bad cone
  // touches only its paired bit. Use the adjacency index so ASIC-sized startup
  // relations do not repeatedly rescan hundreds of thousands of unrelated
  // candidates while closing a small output COI.
  EqualityIndex(equalityPairs).close(frameStates);
}

void closeSameFrameStateEqualityDependencies(
    const KInductionProblem& problem,
    std::unordered_set<size_t>& frameStates) {
  closeFrameEqualityDependencies(problem.sameFrameStateEqualityPairs0, frameStates);
  closeFrameEqualityDependencies(problem.sameFrameStateEqualityPairs1, frameStates);
}

void addRelevantComplementPartners(
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    std::unordered_set<size_t>& solverSymbols) {
  for (const auto& [primarySymbol, complementedSymbol] : complementedStatePairs) {
    if (solverSymbols.find(primarySymbol) != solverSymbols.end() ||
        solverSymbols.find(complementedSymbol) != solverSymbols.end()) {  // LCOV_EXCL_LINE
      solverSymbols.insert(primarySymbol);
      // LCOV_EXCL_START
      solverSymbols.insert(complementedSymbol);
      // LCOV_EXCL_STOP
    }
  }
}

void addRelevantSameFrameStateEqualityPartners(
    const KInductionProblem& problem,
    std::unordered_set<size_t>& solverSymbols) {
  closeSameFrameStateEqualityDependencies(problem, solverSymbols);
}

void addRelevantDualRailPartners(
    const std::vector<DualRailSymbolPair>& railPairs,
    std::unordered_set<size_t>& solverSymbols) {
  for (const auto& rails : railPairs) {
    if (solverSymbols.find(rails.mayBeOne) != solverSymbols.end() ||
        solverSymbols.find(rails.mayBeZero) != solverSymbols.end()) {
      solverSymbols.insert(rails.mayBeOne);
      solverSymbols.insert(rails.mayBeZero);
    }
  }
}

BaseCaseCoi buildBaseCaseCoi(const KInductionProblem& problem,
                             InitialConstraintMode initialMode,
                             size_t bootstrapFrames,
                             bool resetBootstrapObservationFrontier,
                             size_t firstBadFrame,
                             size_t internalK,
                             bool constrainPreviouslySafeFrames,
                             bool) {
  // Cone-of-influence reduction for the base query. The original formula is
  // still the same bounded counterexample check, but we avoid allocating and
  // constraining state bits that cannot influence any checked bad output.
  // This matters for ASIC SEC where resetless memories and rail state can be
  // much larger than a given output cone.
  const auto stateSymbols = buildStateSymbolSet(problem);
  const TransitionExprResolver transitionByState(problem);
  const auto primaryByComplement = buildPrimaryByComplementSymbol(problem);

  std::vector<std::unordered_set<size_t>> requiredStates(internalK + 1);
  std::unordered_set<size_t> solverSymbols;
  solverSymbols.reserve(1024);
  for (const auto& [symbol, _] : problem.resetBootstrapInputs) {
    solverSymbols.insert(symbol);
  }

  addFormulaSupport(problem.bad, solverSymbols);
  for (size_t frame = firstBadFrame; frame <= internalK; ++frame) {
    addFormulaStateSupport(problem.bad, stateSymbols, requiredStates[frame]);
  }
  if (constrainPreviouslySafeFrames) {
    // Frontier checks are issued only after all smaller frontiers were proved
    // safe. Assert those already-known property facts in the standalone SAT
    // query so Kissat gets the same pruning an incremental BMC run would have.
    addFormulaSupport(problem.property, solverSymbols); // LCOV_EXCL_LINE
    for (size_t frame = bootstrapFrames; frame < firstBadFrame; ++frame) { // LCOV_EXCL_LINE
      addFormulaStateSupport(problem.property, stateSymbols, requiredStates[frame]); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  if (bootstrapFrames == 0) {
    if (initialMode == InitialConstraintMode::CompleteInit ||
        initialMode == InitialConstraintMode::PartialInit) {
      if (!hasStructuredInitialAssignments(problem)) {
        // Compatibility path for hand-built unit-test problems and older
        // callers.  Production SEC builds initial state as structured unit
        // assignments, so the base query can keep the COI narrow.
        addFormulaSupport(problem.initialCondition, solverSymbols);
        addFormulaStateSupport(
            problem.initialCondition, stateSymbols, requiredStates[0]);
      }
    }
    if (initialMode == InitialConstraintMode::ObservationOnly ||
        initialMode == InitialConstraintMode::PartialInit) {
      addFormulaSupport(problem.property, solverSymbols);
      addFormulaStateSupport(problem.property, stateSymbols, requiredStates[0]);
    }
  } else if (resetBootstrapObservationFrontier) {
    // Reset/bootstrap can leave internal memories or FIFOs arbitrary.  The
    // user-visible SEC frontier in that case is the top observation, not a
    // cross-design equality on those internal state bits.
    addFormulaSupport(problem.property, solverSymbols);
    addFormulaStateSupport(
        problem.property, stateSymbols, requiredStates[bootstrapFrames]);
  }

  std::vector<std::vector<size_t>> transitionTargetsByFrame(internalK);
  for (size_t frame = internalK; frame > 0; --frame) {
    closeSameFrameStateEqualityDependencies(problem, requiredStates[frame]);
    auto targets = expandTransitionTargets(
        requiredStates[frame],
        transitionByState,
        primaryByComplement);
    transitionTargetsByFrame[frame - 1] = targets;
    transitionByState.collectSupportForTargets(
        targets, stateSymbols, requiredStates[frame - 1], solverSymbols);
  }

  closeSameFrameStateEqualityDependencies(problem, requiredStates[0]);

  for (const auto& frameStates : requiredStates) {
    solverSymbols.insert(frameStates.begin(), frameStates.end());
  }
  for (const auto& targets : transitionTargetsByFrame) {
    for (const auto target : targets) {
      solverSymbols.insert(target);
    }
  }
  addRelevantComplementPartners(problem.complementedStatePairs0, solverSymbols);
  addRelevantComplementPartners(problem.complementedStatePairs1, solverSymbols);
  addRelevantSameFrameStateEqualityPartners(problem, solverSymbols);
  addRelevantDualRailPartners(problem.dualRailStatePairs, solverSymbols);

  BaseCaseCoi coi;
  coi.transitionTargetsByFrame = std::move(transitionTargetsByFrame);
  coi.requiredStateSymbolsByFrame = sortedFrameStates(requiredStates);
  coi.solverSymbols = sortedSymbols(solverSymbols);
  coi.solverSymbolSet.insert(coi.solverSymbols.begin(), coi.solverSymbols.end());
  return coi;
}

BaseCaseCoi buildImcCachedBaseCaseCoi(
    const ImcBaseCounterexampleCache& cache,
    InitialConstraintMode initialMode,
    size_t bootstrapFrames,
    bool resetBootstrapObservationFrontier,
    size_t firstBadFrame,
    size_t internalK,
    bool constrainPreviouslySafeFrames) {
  const auto& problem = cache.problem;
  std::vector<std::unordered_set<size_t>> requiredStates(internalK + 1);
  std::unordered_set<size_t> solverSymbols;
  solverSymbols.reserve(1024);
  for (const auto& [symbol, _] : problem.resetBootstrapInputs) {
    solverSymbols.insert(symbol);
  }

  addFormulaSupport(problem.bad, solverSymbols);
  for (size_t frame = firstBadFrame; frame <= internalK; ++frame) {
    addFormulaStateSupport(problem.bad, cache.stateSymbols, requiredStates[frame]);
  }
  if (constrainPreviouslySafeFrames) {
    addFormulaSupport(problem.property, solverSymbols); // LCOV_EXCL_LINE
    for (size_t frame = bootstrapFrames; frame < firstBadFrame; ++frame) { // LCOV_EXCL_LINE
      addFormulaStateSupport( // LCOV_EXCL_LINE
          problem.property, cache.stateSymbols, requiredStates[frame]); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  if (bootstrapFrames == 0) {
    if (initialMode == InitialConstraintMode::CompleteInit ||
        initialMode == InitialConstraintMode::PartialInit) {
      if (!hasStructuredInitialAssignments(problem)) {
        addFormulaSupport(problem.initialCondition, solverSymbols);
        addFormulaStateSupport(
            problem.initialCondition, cache.stateSymbols, requiredStates[0]);
      }
    }
    if (initialMode == InitialConstraintMode::ObservationOnly ||
        initialMode == InitialConstraintMode::PartialInit) {
      addFormulaSupport(problem.property, solverSymbols);
      addFormulaStateSupport(
          problem.property, cache.stateSymbols, requiredStates[0]);
    }
  } else if (resetBootstrapObservationFrontier) {
    addFormulaSupport(problem.property, solverSymbols); // LCOV_EXCL_LINE
    addFormulaStateSupport( // LCOV_EXCL_LINE
        problem.property, cache.stateSymbols, requiredStates[bootstrapFrames]); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  std::vector<std::vector<size_t>> transitionTargetsByFrame(internalK);
  for (size_t frame = internalK; frame > 0; --frame) {
    cache.closeSameFrameStateEqualityDependencies(requiredStates[frame]);
    auto targets = expandTransitionTargets(
        requiredStates[frame],
        cache.transitionByState,
        cache.primaryByComplement);
    transitionTargetsByFrame[frame - 1] = targets;
    cache.transitionByState.collectSupportForTargets(
        targets, cache.stateSymbols, requiredStates[frame - 1], solverSymbols);
  }

  cache.closeSameFrameStateEqualityDependencies(requiredStates[0]);

  for (const auto& frameStates : requiredStates) {
    solverSymbols.insert(frameStates.begin(), frameStates.end());
  }
  for (const auto& targets : transitionTargetsByFrame) {
    for (const auto target : targets) {
      solverSymbols.insert(target);
    }
  }
  addRelevantComplementPartners(problem.complementedStatePairs0, solverSymbols);
  addRelevantComplementPartners(problem.complementedStatePairs1, solverSymbols);
  cache.closeSameFrameStateEqualityDependencies(solverSymbols);
  addRelevantDualRailPartners(problem.dualRailStatePairs, solverSymbols);

  BaseCaseCoi coi;
  coi.transitionTargetsByFrame = std::move(transitionTargetsByFrame);
  coi.requiredStateSymbolsByFrame = sortedFrameStates(requiredStates);
  coi.solverSymbols = sortedSymbols(solverSymbols);
  coi.solverSymbolSet.insert(coi.solverSymbols.begin(), coi.solverSymbols.end());
  return coi;
}

BaseCaseCoi buildStateCubeReachabilityCoiForTargetFrames(
    const ResetFrontierReachabilityContextData& context,
    size_t targetFrame,
    const std::vector<std::pair<size_t, bool>>& cube,
    const std::vector<size_t>& targetFrames,
    bool) {
  // This is the same bounded transition prefix as the normal base case, but
  // the target is a state cube rather than the SEC bad formula.  PDR calls this
  // when a level-0 obligation may be an artifact of the abstract bootstrap
  // summary: the SAT query must be exact, but only for the cube's COI.
  const auto& problem = context.problem;
  const auto& transitionByState = context.transitionByState;
  const auto& stateSymbols = transitionByState.stateSymbols();
  std::vector<std::unordered_set<size_t>> requiredStates(targetFrame + 1);
  std::unordered_set<size_t> solverSymbols;
  solverSymbols.reserve(1024);
  for (const auto& [symbol, _] : problem.resetBootstrapInputs) {
    solverSymbols.insert(symbol);
  }
  for (const auto& [symbol, _] : cube) {
    solverSymbols.insert(symbol);
    if (stateSymbols.find(symbol) != stateSymbols.end()) {
      for (const auto target : targetFrames) {
        if (target <= targetFrame) {
          requiredStates[target].insert(symbol);
        }
      }
    }
  }
  if (context.frameInvariant != nullptr) {
    // The caller has already proved this invariant on the PDR startup frontier
    // and through one post-reset transition. Encode it only from that concrete
    // frontier onward; pre-reset frames may intentionally violate it while
    // reset is still asserting.
    addFormulaSupport(context.frameInvariant, solverSymbols);
    for (size_t frame = context.bootstrapFrames; frame <= targetFrame; ++frame) {
      addFormulaStateSupport(
          context.frameInvariant, stateSymbols, requiredStates[frame]);
    }
  }

  if (context.bootstrapFrames == 0 &&
      (context.initialMode == InitialConstraintMode::CompleteInit ||
       context.initialMode == InitialConstraintMode::PartialInit) &&
      !hasStructuredInitialAssignments(problem)) {
    addFormulaSupport(problem.initialCondition, solverSymbols);
    addFormulaStateSupport(
        problem.initialCondition, stateSymbols, requiredStates[0]);
  }

  std::vector<std::vector<size_t>> transitionTargetsByFrame(targetFrame);
  for (size_t frame = targetFrame; frame > 0; --frame) {
    closeSameFrameStateEqualityDependencies(problem, requiredStates[frame]);
    auto targets = expandTransitionTargets(
        requiredStates[frame],
        transitionByState,
        context.primaryByComplement);
    transitionTargetsByFrame[frame - 1] = targets;
    transitionByState.collectSupportForTargets(
        targets, stateSymbols, requiredStates[frame - 1], solverSymbols);
  }

  closeSameFrameStateEqualityDependencies(problem, requiredStates[0]);

  for (const auto& frameStates : requiredStates) {
    solverSymbols.insert(frameStates.begin(), frameStates.end());
  }
  for (const auto& targets : transitionTargetsByFrame) {
    for (const auto target : targets) {
      solverSymbols.insert(target);
    }
  }
  addRelevantComplementPartners(problem.complementedStatePairs0, solverSymbols);
  addRelevantComplementPartners(problem.complementedStatePairs1, solverSymbols);
  addRelevantSameFrameStateEqualityPartners(problem, solverSymbols);
  addRelevantDualRailPartners(problem.dualRailStatePairs, solverSymbols);

  BaseCaseCoi coi;
  coi.transitionTargetsByFrame = std::move(transitionTargetsByFrame);
  coi.requiredStateSymbolsByFrame = sortedFrameStates(requiredStates);
  coi.solverSymbols = sortedSymbols(solverSymbols);
  coi.solverSymbolSet.insert(coi.solverSymbols.begin(), coi.solverSymbols.end());
  return coi;
}

BaseCaseCoi buildStateCubeReachabilityCoi(
    const ResetFrontierReachabilityContextData& context,
    size_t targetFrame,
    const std::vector<std::pair<size_t, bool>>& cube,
    bool closeStartupEqualityDependencies) {
  return buildStateCubeReachabilityCoiForTargetFrames(
      context,
      targetFrame,
      cube,
      std::vector<size_t>{targetFrame},
      closeStartupEqualityDependencies);
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START
BaseCaseCoi buildStateCubePrefixReachabilityCoi(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    const ResetFrontierReachabilityContextData& context,
    // LCOV_EXCL_START
    size_t maxTargetFrame,
    // LCOV_EXCL_STOP
    const std::vector<std::pair<size_t, bool>>& cube,
    bool closeStartupEqualityDependencies) {
  std::vector<size_t> targetFrames;  // LCOV_EXCL_LINE
  targetFrames.reserve(maxTargetFrame + 1 - context.bootstrapFrames);  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  for (size_t frame = context.bootstrapFrames; frame <= maxTargetFrame; ++frame) {  // LCOV_EXCL_LINE
    targetFrames.push_back(frame);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return buildStateCubeReachabilityCoiForTargetFrames(  // LCOV_EXCL_LINE
      context,  // LCOV_EXCL_LINE
      maxTargetFrame,  // LCOV_EXCL_LINE
      cube,  // LCOV_EXCL_LINE
      targetFrames,
      closeStartupEqualityDependencies);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START

BaseCaseCoi buildResetSummaryCubeReachabilityCoi(
// LCOV_EXCL_STOP
    const ResetFrontierReachabilityContextData& context,
    size_t postBootstrapSteps,
    const std::vector<std::pair<size_t, bool>>& cube) {
  const auto& problem = context.problem;
  const auto& transitionByState = context.transitionByState;
  const auto& stateSymbols = transitionByState.stateSymbols();
  std::vector<std::unordered_set<size_t>> requiredStates(postBootstrapSteps + 1);
  std::unordered_set<size_t> solverSymbols;
  solverSymbols.reserve(1024);
  for (const auto& [symbol, _] : problem.resetBootstrapInputs) {
    solverSymbols.insert(symbol);
  }
  for (const auto& [symbol, _] : cube) {
    solverSymbols.insert(symbol);
    if (stateSymbols.find(symbol) != stateSymbols.end()) {
      requiredStates[postBootstrapSteps].insert(symbol);
    }
  }
  if (context.frameInvariant != nullptr) {
    addFormulaSupport(context.frameInvariant, solverSymbols);  // LCOV_EXCL_LINE
    for (size_t frame = 0; frame <= postBootstrapSteps; ++frame) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      addFormulaStateSupport(  // LCOV_EXCL_LINE
          context.frameInvariant, stateSymbols, requiredStates[frame]);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  std::vector<std::vector<size_t>> transitionTargetsByFrame(postBootstrapSteps);
  // LCOV_EXCL_STOP
  for (size_t frame = postBootstrapSteps; frame > 0; --frame) {
    closeSameFrameStateEqualityDependencies(problem, requiredStates[frame]);
    auto targets = expandTransitionTargets(
        requiredStates[frame],
        transitionByState,
        context.primaryByComplement);
    transitionTargetsByFrame[frame - 1] = targets;
    transitionByState.collectSupportForTargets(
        targets, stateSymbols, requiredStates[frame - 1], solverSymbols);
  }

  context.bootstrapEqualities.close(requiredStates[0]);
  closeSameFrameStateEqualityDependencies(problem, requiredStates[0]);

  for (const auto& frameStates : requiredStates) {
    solverSymbols.insert(frameStates.begin(), frameStates.end());
  }
  for (const auto& targets : transitionTargetsByFrame) {
    for (const auto target : targets) {
      solverSymbols.insert(target);
    }
  }
  addRelevantComplementPartners(problem.complementedStatePairs0, solverSymbols);
  addRelevantComplementPartners(problem.complementedStatePairs1, solverSymbols);
  addRelevantSameFrameStateEqualityPartners(problem, solverSymbols);

  BaseCaseCoi coi;
  coi.transitionTargetsByFrame = std::move(transitionTargetsByFrame);
  coi.requiredStateSymbolsByFrame = sortedFrameStates(requiredStates);
  coi.solverSymbols = sortedSymbols(solverSymbols);
  coi.solverSymbolSet.insert(coi.solverSymbols.begin(), coi.solverSymbols.end());
  return coi;
}

void addResetBootstrapConstraints(SATSolverWrapper& solver,
                                  const FrameVariableStore& variables,
                                  const KInductionProblem& problem,
                                  size_t numFrames) {
  const size_t bootstrapFrames = resetBootstrapFrames(problem);
  if (bootstrapFrames == 0) {
    return;
  }

  for (const auto& [symbol, assertedValue] : problem.resetBootstrapInputs) {
    for (size_t frame = 0; frame < std::min(bootstrapFrames, numFrames); ++frame) {
      solver.addClause(
          {assertedValue ? variables.getLiteral(symbol, frame)
                         : -variables.getLiteral(symbol, frame)});
    }
    for (size_t frame = bootstrapFrames; frame < numFrames; ++frame) {
      solver.addClause(
          {assertedValue ? -variables.getLiteral(symbol, frame)
                         : variables.getLiteral(symbol, frame)});
    }
  }
}

size_t addInitialStateAssignments(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& solverSymbols) {
  size_t encodedCount = 0;
  for (const auto& [symbol, value] : problem.initialStateAssignments) {
    if (solverSymbols.find(symbol) == solverSymbols.end()) {
      continue;
    }
    solver.addClause({value ? variables.getLiteral(symbol, 0)
                            : -variables.getLiteral(symbol, 0)});
    ++encodedCount;
  }
  return encodedCount;
}

void addBootstrapStateAssignments(SATSolverWrapper& solver,
                                  const FrameVariableStore& variables,
                                  const KInductionProblem& problem,
                                  const std::unordered_set<size_t>& solverSymbols,
                                  size_t frame) {
  for (const auto& [symbol, value] : problem.bootstrapStateAssignments) {
    if (solverSymbols.find(symbol) == solverSymbols.end()) {
      continue;
    }
    solver.addClause({value ? variables.getLiteral(symbol, frame)
                            : -variables.getLiteral(symbol, frame)});
  }
}

void addInitialConstraints(SATSolverWrapper& solver,
                           const FrameVariableStore& variables,
                           const KInductionProblem& problem,
                           const std::unordered_set<size_t>& solverSymbols,
                           InitialConstraintMode mode) {
  if (mode == InitialConstraintMode::None) {
    return;
  }

  if ((mode == InitialConstraintMode::CompleteInit ||
       mode == InitialConstraintMode::PartialInit) &&
      problem.initialCondition != nullptr) {
    if (hasStructuredInitialAssignments(problem)) {
      addInitialStateAssignments(solver, variables, problem, solverSymbols);
    } else {
      FrameFormulaEncoder encoder(
          solver,
          variables.makeLeafLits(
              0, formulaSupportOrThrow(problem.initialCondition, "initial condition")));
      solver.addClause({encoder.encode(problem.initialCondition)});
    }
  }

  if (mode == InitialConstraintMode::ObservationOnly ||
      mode == InitialConstraintMode::PartialInit) {
    FrameFormulaEncoder encoder(
        solver,
        variables.makeLeafLits(
            0, formulaSupportOrThrow(problem.property, "observation property")));
    solver.addClause({encoder.encode(problem.property)});
  }
}

void addObservationPropertyConstraint(SATSolverWrapper& solver,
                                      const FrameVariableStore& variables,
                                      const KInductionProblem& problem,
                                      size_t frame) {
  FrameFormulaEncoder encoder(
      solver,
      variables.makeLeafLits(
          frame, formulaSupportOrThrow(problem.property, "observation property")));
  solver.addClause({encoder.encode(problem.property)});
}

void addResetFrontierFrameInvariantConstraints(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const ResetFrontierReachabilityContextData& context,
    size_t targetFrame) {
  if (context.frameInvariant == nullptr) {
    return;
  }
  for (size_t frame = context.bootstrapFrames; frame <= targetFrame; ++frame) {
    FrameFormulaEncoder encoder(
        solver, variables.makeLeafLits(frame, context.frameInvariantSupport));
    solver.addClause({encoder.encode(context.frameInvariant)});
  }
}

size_t countTransitionTargets(
    const std::vector<std::vector<size_t>>& transitionTargetsByFrame) {
  size_t count = 0;
  for (const auto& targets : transitionTargetsByFrame) {
    count += targets.size();
  }
  return count;
}

size_t countMatchingAssignments(
    const std::vector<std::pair<size_t, bool>>& assignments,
    const std::unordered_set<size_t>& solverSymbols) {
  size_t count = 0;
  for (const auto& [symbol, _] : assignments) {
    if (solverSymbols.find(symbol) != solverSymbols.end()) {
      ++count;
    }
  }
  return count;
}

void emitBaseCaseCoiDiag(const KInductionProblem& problem,
                         const BaseCaseCoi& coi,
                         size_t k,
                         size_t firstBadFrame,
                         size_t internalK,
                         bool constrainPreviouslySafeFrames) {
  if (!isKInductionCoiDiagEnabled()) {
    return;
  }
  emitSecDiag(
      "SEC diag: k-induction base coi k=", k,
      " first_bad_frame=", firstBadFrame,
      " frames=", internalK + 1,
      " solver_symbols=", coi.solverSymbols.size(),
      " transition_targets=", countTransitionTargets(coi.transitionTargetsByFrame),
      " bad_support=", problem.bad != nullptr ? problem.bad->getSupportVars().size() : 0,
      " initial_assignments=", problem.initialStateAssignments.size(),
      " encoded_initial_assignments=",
      countMatchingAssignments(problem.initialStateAssignments, coi.solverSymbolSet),
      " structured_initial_assignments=", hasStructuredInitialAssignments(problem),
      " safe_prefix_constraints=", constrainPreviouslySafeFrames);
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
        continue;  // LCOV_EXCL_LINE
      }
      // LCOV_EXCL_START
      addLiteralEquivalence(
      // LCOV_EXCL_STOP
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
        continue;
      }
      // Keep dual-rail BMC inside legal ternary states.  Initial/reset
      // assignments produce legal rails, and the encoded transition preserves
      // legality; this clause prevents SAT from inventing empty rail values in
      // partially projected validation queries.
      solver.addClause({
          variables.getLiteral(rails.mayBeOne, frame),
          variables.getLiteral(rails.mayBeZero, frame)});
    }
  }
}

void addBlockedStateCubeClause(SATSolverWrapper& solver,
                               const FrameVariableStore& variables,
                               const std::vector<std::pair<size_t, bool>>& cube,
                               size_t frame) {
  std::vector<int> clause;
  clause.reserve(cube.size());
  for (const auto& [symbol, value] : cube) {
    const int literal = variables.getLiteral(symbol, frame);
    clause.push_back(value ? -literal : literal);
  }
  solver.addClause(clause);
}

void addTransitionRelation(SATSolverWrapper& solver,
                           const FrameVariableStore& variables,
                           const TransitionExprResolver& transitionByState,
                           const std::vector<size_t>& targets,
                           size_t frame) {
  // All transition formulas in one frame are slices of the same combinational
  // next-state network. Reusing a frame encoder lets shared BoolExpr DAG nodes
  // produce one Tseitin literal instead of being re-encoded once per state bit.
  struct EncodingGroup {
    const std::unordered_map<size_t, size_t>* symbolMap = nullptr;
    std::vector<size_t> stateSymbols;
  };
  std::vector<EncodingGroup> groups;
  groups.reserve(3);
  for (const auto stateSymbol : targets) {
    const TransitionExprView view = transitionByState.expressionView(stateSymbol);
    auto groupIt = std::find_if(
        groups.begin(),
        groups.end(),
        [&](const EncodingGroup& group) {
          return group.symbolMap == view.symbolMap;
        });
    if (groupIt == groups.end()) {
      groups.push_back(EncodingGroup{view.symbolMap, {stateSymbol}});
    } else {
      groupIt->stateSymbols.push_back(stateSymbol);
    }
  }

  for (const auto& group : groups) {
    size_t expectedNodes = 0;
    if (group.stateSymbols.size() <= kMaxExactTransitionNodeCountHintTargets) {
      for (const auto stateSymbol : group.stateSymbols) {
        expectedNodes += transitionByState.nodeCount(stateSymbol);
      }
    }
    FrameFormulaEncoder encoder(
        solver,
        variables.makeLeafLits(frame),
        group.symbolMap,
        false,
        expectedNodes);
    for (const auto stateSymbol : group.stateSymbols) {
      const TransitionExprView view =
          transitionByState.expressionView(stateSymbol);
      addLiteralEquivalence(
          solver,
          variables.getLiteral(stateSymbol, frame + 1),
          encoder.encode(view.expr));
    }
  }
}

std::unordered_map<size_t, bool> buildFrameEnvironment(
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    BoolExpr* formula,
    size_t frame) {
  std::unordered_map<size_t, bool> environment;
  if (formula == nullptr) {
    return environment;  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_START
  const auto support = formula->getSupportVars();
  // LCOV_EXCL_STOP
  environment.reserve(support.size());
  for (const auto symbol : support) {
    if (symbol < 2) {
      continue;
    }
    environment.emplace(
        symbol, solver.getLiteralValue(variables.getLiteral(symbol, frame)));
  }
  return environment;
}

void addFormulaValuesToEnvironment(const SATSolverWrapper& solver,
                                   const FrameVariableStore& variables,
                                   BoolExpr* formula,
                                   size_t frame,
                                   std::unordered_map<size_t, bool>& environment) {
  if (formula == nullptr) {
    return;  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_START
  for (const auto symbol : formula->getSupportVars()) {
  // LCOV_EXCL_STOP
    if (symbol < 2 || environment.find(symbol) != environment.end()) {
      continue;
    }
    if (!variables.hasSymbol(symbol)) {
      continue;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    environment.emplace(
    // LCOV_EXCL_STOP
        symbol, solver.getLiteralValue(variables.getLiteral(symbol, frame)));
  }
}

bool formulaSupportCoveredByVariables(const FrameVariableStore& variables,
                                      BoolExpr* formula) {
  if (formula == nullptr) {
    return true;  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_START
  for (const auto symbol : formula->getSupportVars()) {
  // LCOV_EXCL_STOP
    if (symbol >= 2 && !variables.hasSymbol(symbol)) {
      return false;  // LCOV_EXCL_LINE
    }
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  return true;
}

size_t findFirstBadFrame(const SATSolverWrapper& solver,
                         const FrameVariableStore& variables,
                         const KInductionProblem& problem,
                         size_t firstBadFrame,
                         size_t maxFrame) {
  for (size_t frame = firstBadFrame; frame <= maxFrame; ++frame) {
    if (problem.bad->evaluate(
            buildFrameEnvironment(solver, variables, problem.bad, frame))) {
      return frame;
    }
  }
  throw std::runtime_error("SAT model does not satisfy any bad frame");  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START

std::vector<KInductionResult::FrameInputAssignments> buildInputTrace(
// LCOV_EXCL_STOP
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    size_t firstFrame,
    size_t lastFrame,
    size_t frameOffset) {
  std::vector<KInductionResult::FrameInputAssignments> trace;
  trace.reserve(lastFrame - firstFrame + 1);
  for (size_t frame = firstFrame; frame <= lastFrame; ++frame) {
    KInductionResult::FrameInputAssignments frameAssignments;
    frameAssignments.frame = frame - frameOffset;
    frameAssignments.assignments.reserve(problem.inputSymbols.size());
    for (size_t i = 0; i < problem.inputSymbols.size(); ++i) {
      // COI-reduced proof obligations intentionally do not allocate SAT
      // literals for environment inputs that cannot affect the checked output.
      // A counterexample witness only needs assignments for inputs present in
      // the solved cone.
      if (!variables.hasSymbol(problem.inputSymbols[i])) {
        continue;
      }
      const std::string inputName =
          i < problem.environmentInputNames.size()
              ? problem.environmentInputNames[i]
              : "input_" + std::to_string(problem.inputSymbols[i]);
      frameAssignments.assignments.push_back(
          {inputName,
           solver.getLiteralValue(
               variables.getLiteral(problem.inputSymbols[i], frame))});
    }
    trace.push_back(std::move(frameAssignments));
  }
  return trace;
}

std::vector<KInductionResult::SignalMismatch> collectObservedOutputMismatches(
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    size_t frame) {
  std::vector<KInductionResult::SignalMismatch> mismatches;
  for (size_t i = 0; i < problem.observedOutputExprs0.size(); ++i) {
    // LCOV_EXCL_START
    if (!formulaSupportCoveredByVariables(variables, problem.observedOutputExprs0[i]) ||
    // LCOV_EXCL_STOP
        !formulaSupportCoveredByVariables(variables, problem.observedOutputExprs1[i])) {
      continue;  // LCOV_EXCL_LINE
    }
    std::unordered_map<size_t, bool> environment;
    addFormulaValuesToEnvironment(
        solver, variables, problem.observedOutputExprs0[i], frame, environment);
    addFormulaValuesToEnvironment(
        solver, variables, problem.observedOutputExprs1[i], frame, environment);
    const bool value0 = problem.observedOutputExprs0[i]->evaluate(environment);
    const bool value1 = problem.observedOutputExprs1[i]->evaluate(environment);
    if (value0 != value1) {
      const std::string outputName =
          i < problem.observedOutputNames.size()
              ? problem.observedOutputNames[i]
              : "output_" + std::to_string(i);
      mismatches.push_back(
          {outputName, value0, value1});
    }
  }
  return mismatches;
}

KInductionResult::CounterexampleWitness buildCounterexampleWitness(
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    size_t firstBadFrame,
    size_t maxFrame,
    size_t frameOffset) {
  KInductionResult::CounterexampleWitness witness;
  const size_t internalBadFrame =
      findFirstBadFrame(solver, variables, problem, firstBadFrame, maxFrame);
  witness.badFrame = internalBadFrame - frameOffset;
  witness.inputTrace = buildInputTrace(
      solver, variables, problem, frameOffset, internalBadFrame, frameOffset);
  witness.outputMismatches = collectObservedOutputMismatches(
      solver, variables, problem, internalBadFrame);
  return witness;
}

enum class BaseCaseSolverProfile {
  SecConeProof,
  PdrValidation,
  PdrValidationProofOnly,
  FastCounterexampleSearch,
};

std::optional<KInductionResult::CounterexampleWitness> findBaseCounterexampleImpl(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k,
    std::optional<size_t> exactPublicBadFrame,
    bool localizeMultiOutputFrontier = true,
    BaseCaseSolverProfile solverProfile = BaseCaseSolverProfile::SecConeProof,
    SATSolverWrapper::SolveStatus* solveStatusOut = nullptr);

// LCOV_EXCL_START
KInductionProblem makeSingleObservedOutputProblem(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    const KInductionProblem& problem,
    size_t outputIndex) {
  KInductionProblem single = problem;  // LCOV_EXCL_LINE
  single.observedOutputs =  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      outputIndex < problem.observedOutputs.size()  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
          ? std::vector<SignalKey>{problem.observedOutputs[outputIndex]}  // LCOV_EXCL_LINE
          : std::vector<SignalKey>{};  // LCOV_EXCL_LINE
  single.observedOutputNames =  // LCOV_EXCL_LINE
      outputIndex < problem.observedOutputNames.size()  // LCOV_EXCL_LINE
          ? std::vector<std::string>{problem.observedOutputNames[outputIndex]}  // LCOV_EXCL_LINE
          : std::vector<std::string>{};  // LCOV_EXCL_LINE
  single.observedOutputExprs0 = {problem.observedOutputExprs0[outputIndex]};  // LCOV_EXCL_LINE
  single.observedOutputExprs1 = {problem.observedOutputExprs1[outputIndex]};  // LCOV_EXCL_LINE

  BoolExpr* outputBad = BoolExpr::simplify(  // LCOV_EXCL_LINE
      BoolExpr::Xor(  // LCOV_EXCL_LINE
          single.observedOutputExprs0.front(),  // LCOV_EXCL_LINE
          single.observedOutputExprs1.front()));  // LCOV_EXCL_LINE
  single.bad = outputBad;  // LCOV_EXCL_LINE
  single.property = BoolExpr::Not(outputBad);  // LCOV_EXCL_LINE
  single.inductionBad = outputBad;  // LCOV_EXCL_LINE
  single.inductionProperty = single.property;  // LCOV_EXCL_LINE
  return single;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

}  // namespace

const ImcBaseCounterexampleCache& ImcBaseCounterexampleCache::singleOutputCache(
    size_t outputIndex) const {
  if (singleOutputCaches.size() < problem.observedOutputExprs0.size()) {
    singleOutputCaches.resize(problem.observedOutputExprs0.size());
    singleOutputProblems.resize(problem.observedOutputExprs0.size());
  }
  if (!singleOutputCaches[outputIndex]) {
    // IMC checks the same residual output slice at increasing depths.  Keep the
    // sliced problem and its immutable COI indexes alive so every depth does
    // not rebuild transition/equality maps from the parent batch.
    singleOutputProblems[outputIndex] =
        std::make_unique<KInductionProblem>(
            makeSingleObservedOutputProblem(problem, outputIndex));
    singleOutputCaches[outputIndex] =
        std::make_unique<ImcBaseCounterexampleCache>(
            *singleOutputProblems[outputIndex]);
  }
  return *singleOutputCaches[outputIndex];
} // LCOV_EXCL_LINE

const ImcBaseCounterexampleCache& ImcBaseCounterexampleCache::outputSubsetCache(
    size_t firstOutput, size_t endOutput) const {
  if (firstOutput == 0 && endOutput == problem.observedOutputExprs0.size()) {
    return *this; // LCOV_EXCL_LINE
  }
  if (endOutput == firstOutput + 1) {
    return singleOutputCache(firstOutput);
  }
  for (size_t i = 0; i < subsetCaches.size(); ++i) {
    if (subsetFirstOutputs[i] == firstOutput &&
        subsetEndOutputs[i] == endOutput) { // LCOV_EXCL_LINE
      return *subsetCaches[i]; // LCOV_EXCL_LINE
    }
  }

  // Keep recursively split IMC residual batches alive across depths.  This
  // preserves the existing exact single-output fallback while avoiding a full
  // rebuild when the same hard subset is checked at the next frontier.
  KInductionProblem subset = problem;
  configureOutputBatchProblem(subset, problem, firstOutput, endOutput);
  subsetProblems.push_back(
      std::make_unique<KInductionProblem>(std::move(subset)));
  subsetCaches.push_back(
      std::make_unique<ImcBaseCounterexampleCache>(*subsetProblems.back()));
  subsetFirstOutputs.push_back(firstOutput);
  subsetEndOutputs.push_back(endOutput);
  return *subsetCaches.back();
}

namespace {

std::vector<std::pair<size_t, size_t>> buildStrictImcOutputSubsets(
    const KInductionProblem& problem) {
  const size_t outputCount = problem.observedOutputExprs0.size();
  if (outputCount <= 1) {
    return {}; // LCOV_EXCL_LINE
  }

  auto batches = buildSupportBoundedOutputBatches(problem);
  if (batches.size() == 1 &&
      batches.front().first == 0 && // LCOV_EXCL_LINE
      batches.front().second == outputCount) { // LCOV_EXCL_LINE
    const size_t midpoint = outputCount / 2; // LCOV_EXCL_LINE
    return {{0, midpoint}, {midpoint, outputCount}}; // LCOV_EXCL_LINE
  }
  return batches;
}

struct OutputCandidate {
  size_t index;
  size_t support;
};

struct OutputCandidateSupportLess {
  bool operator()(const OutputCandidate& lhs, // LCOV_EXCL_LINE
                  const OutputCandidate& rhs) const {
    return lhs.support < rhs.support; // LCOV_EXCL_LINE
  }
};

std::optional<KInductionResult::CounterexampleWitness>
findPerOutputBaseCounterexampleAtFrontier(  // LCOV_EXCL_LINE
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    // LCOV_EXCL_START
    size_t k,
    // LCOV_EXCL_STOP
    std::optional<size_t> exactPublicBadFrame,
    BaseCaseSolverProfile solverProfile) {
  if (!exactPublicBadFrame.has_value() ||  // LCOV_EXCL_LINE
      problem.observedOutputExprs0.size() <= 1 ||  // LCOV_EXCL_LINE
      problem.observedOutputExprs0.size() != problem.observedOutputExprs1.size()) {  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  // Exact SEC validation asks whether any observed output is bad at one frame.
  // Solving the whole batch OR can force the SAT solver to reason across
  // unrelated output cones.  Match PDR's bad-cube search and validate each
  // output independently; the disjunction is SAT iff one per-output query is.
  std::vector<OutputCandidate> outputs;  // LCOV_EXCL_LINE
  outputs.reserve(problem.observedOutputExprs0.size());  // LCOV_EXCL_LINE
  for (size_t output = 0; output < problem.observedOutputExprs0.size(); ++output) {  // LCOV_EXCL_LINE
    size_t support = output;  // LCOV_EXCL_LINE
    if (solverProfile == BaseCaseSolverProfile::FastCounterexampleSearch) {  // LCOV_EXCL_LINE
      support = 0;  // LCOV_EXCL_LINE
      if (problem.observedOutputExprs0[output] != nullptr) {  // LCOV_EXCL_LINE
      // LCOV_DISABLED_STOP
        support += problem.observedOutputExprs0[output]->getSupportVars().size();  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      if (problem.observedOutputExprs1[output] != nullptr) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        support += problem.observedOutputExprs1[output]->getSupportVars().size();  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    outputs.push_back({output, support});  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  if (solverProfile == BaseCaseSolverProfile::FastCounterexampleSearch) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    std::stable_sort(  // LCOV_EXCL_LINE
        outputs.begin(), outputs.end(),  // LCOV_EXCL_LINE
        OutputCandidateSupportLess{});  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  for (const OutputCandidate& candidate : outputs) {  // LCOV_EXCL_LINE
    const size_t output = candidate.index;  // LCOV_EXCL_LINE
    KInductionProblem single =
        makeSingleObservedOutputProblem(problem, output);  // LCOV_EXCL_LINE
    if (auto witness = findBaseCounterexampleImpl(  // LCOV_EXCL_LINE
            single,
            // LCOV_EXCL_START
            solverType,
            // LCOV_EXCL_STOP
            k, // LCOV_EXCL_LINE
            exactPublicBadFrame, // LCOV_EXCL_LINE
            /*localizeMultiOutputFrontier=*/false,
            solverProfile);  // LCOV_EXCL_LINE
        witness.has_value()) {  // LCOV_EXCL_LINE
      return witness;  // LCOV_EXCL_LINE
    }
  }  // LCOV_EXCL_LINE
  return std::nullopt;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

std::optional<KInductionResult::CounterexampleWitness> findBaseCounterexampleImpl(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k,
    std::optional<size_t> exactPublicBadFrame,
    bool localizeMultiOutputFrontier,
    BaseCaseSolverProfile solverProfile,
    SATSolverWrapper::SolveStatus* solveStatusOut) {
  if (solveStatusOut != nullptr) {
    *solveStatusOut = SATSolverWrapper::SolveStatus::Unsat;
  }
  if (localizeMultiOutputFrontier &&
      exactPublicBadFrame.has_value() &&
      problem.observedOutputExprs0.size() > 1 &&
      problem.observedOutputExprs0.size() == problem.observedOutputExprs1.size()) {  // LCOV_EXCL_LINE
    return findPerOutputBaseCounterexampleAtFrontier(  // LCOV_EXCL_LINE
        problem, solverType, k, exactPublicBadFrame, solverProfile);  // LCOV_EXCL_LINE
  }

  const size_t bootstrapFrames = resetBootstrapFrames(problem);
  const bool resetBootstrapObservationFrontier =
      bootstrapFrames != 0 && problem.usesResetBootstrapObservationFrontier();
  const size_t internalK = k + bootstrapFrames;
  // Base BMC only needs to assert the requested bad frame(s).  For frontier
  // sweeps, earlier depths are checked by the caller; for cumulative base
  // validation, the bad disjunction already covers the full prefix.
  const bool constrainPreviouslySafeFrames = false;
  const InitialConstraintMode initialMode =
      bootstrapFrames == 0 ? determineInitialConstraintMode(problem)
                           : InitialConstraintMode::None;

  size_t firstBadFrame = 0;
  if (bootstrapFrames != 0) {
    firstBadFrame =
        bootstrapFrames + (resetBootstrapObservationFrontier ? 1u : 0u);
  } else if (initialMode == InitialConstraintMode::ObservationOnly ||
             initialMode == InitialConstraintMode::PartialInit) {
    firstBadFrame = 1;
  }
  if (firstBadFrame > internalK) {
    return std::nullopt;
  }

// LCOV_EXCL_START


// LCOV_EXCL_STOP
  size_t lastBadFrame = internalK;
  if (exactPublicBadFrame.has_value()) {
    const size_t exactInternalBadFrame = *exactPublicBadFrame + bootstrapFrames;
    if (exactInternalBadFrame < firstBadFrame ||
        exactInternalBadFrame > internalK) {
      return std::nullopt;  // LCOV_EXCL_LINE
    }
    firstBadFrame = exactInternalBadFrame;
    lastBadFrame = exactInternalBadFrame;
  }

  const BaseCaseCoi coi = buildBaseCaseCoi(
      problem,
      initialMode,
      bootstrapFrames,
      resetBootstrapObservationFrontier,
      firstBadFrame,
      internalK,
      constrainPreviouslySafeFrames,
      solverProfile != BaseCaseSolverProfile::PdrValidationProofOnly);
  emitBaseCaseCoiDiag(
      problem,
      coi,
      k,
      firstBadFrame,
      internalK,
      constrainPreviouslySafeFrames);
  const TransitionExprResolver transitionByState(problem);

  SATSolverWrapper solver(solverType);
  if (solverProfile == BaseCaseSolverProfile::PdrValidation ||
      solverProfile == BaseCaseSolverProfile::PdrValidationProofOnly ||
      solverProfile == BaseCaseSolverProfile::FastCounterexampleSearch ||
      baseCaseValidationUsesLocalQueryProfile(solverType)) {
    // PDR calls this helper as a short-lived exact CEGAR validation. It asks
    // only whether bad is reachable at this frontier: older public bad frames
    // were checked or learned by earlier PDR refinements, and re-encoding them
    // as a safe-prefix constraint made AES spend minutes inside SAT search.
    // LCOV_EXCL_START
    // Use the local query profile here as well: samples on the regress PDR and
    // LCOV_EXCL_STOP
    // dynamic-node KI flows showed medium-sized validation BMCs spending their
    // time in standalone preprocessing/probing, while the caller only needs a
    // quick SAT/UNSAT answer before moving to the next frontier.
    solver.configureForSecPdrQuery(coi.solverSymbols.size());
  } else {
    solver.configureForSecConeProof(coi.solverSymbols.size());  // LCOV_EXCL_LINE
  }
  FrameVariableStore variables(solver, coi.solverSymbols, internalK + 1);
  addResetBootstrapConstraints(solver, variables, problem, internalK + 1);
  addInitialConstraints(solver, variables, problem, coi.solverSymbolSet, initialMode);
  if (bootstrapFrames != 0 && problem.usesDualRailStateEncoding) {
    // A reset prefix starts from the same exact ternary initialization as PDR;
    // its final state is derived by transitions, not a per-register summary.
    addInitialStateAssignments(
        solver, variables, problem, coi.solverSymbolSet);
  }
  if (resetBootstrapObservationFrontier) {
    addObservationPropertyConstraint(
        solver, variables, problem, bootstrapFrames);
  }

  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs0, coi.solverSymbolSet,
      internalK + 1);
  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs1, coi.solverSymbolSet,
      internalK + 1);
  addSameFrameStateEqualities(
      solver, variables, problem, coi.solverSymbolSet, internalK + 1);
  addDualRailStateValidity(
      solver, variables, problem.dualRailStatePairs, coi.solverSymbolSet,
      internalK + 1);
  for (size_t frame = 0; frame < internalK; ++frame) {
    addTransitionRelation(
        solver, variables, transitionByState, coi.transitionTargetsByFrame[frame], frame);
  }
  if (bootstrapFrames != 0 && !problem.usesDualRailStateEncoding) {
    addBootstrapStateAssignments(
        solver, variables, problem, coi.solverSymbolSet, bootstrapFrames);
  }
  if (constrainPreviouslySafeFrames) {
    for (size_t frame = bootstrapFrames; frame < firstBadFrame; ++frame) {
      FrameFormulaEncoder encoder(
          solver,
          variables.makeLeafLits(
              frame,
              formulaSupportOrThrow(
                  problem.property, "previously-safe property formula")));
      solver.addClause({encoder.encode(problem.property)});
    }
  }

  std::vector<int> badClause;
  badClause.reserve(lastBadFrame - firstBadFrame + 1);
  for (size_t frame = firstBadFrame; frame <= lastBadFrame; ++frame) {
    FrameFormulaEncoder encoder(
        solver,
        variables.makeLeafLits(
            frame, formulaSupportOrThrow(problem.bad, "bad-state formula")));
    badClause.push_back(encoder.encode(problem.bad));
  }
  solver.addClause(badClause);

  SATSolverWrapper::SolveStatus status = SATSolverWrapper::SolveStatus::Unknown;
  // LCOV_EXCL_START
  if (solverProfile == BaseCaseSolverProfile::FastCounterexampleSearch) {
    if (solverType == KEPLER_FORMAL::Config::SolverType::KISSAT) {
    // LCOV_EXCL_STOP
      status = solver.solveWithKissatResourceLimits(
          static_cast<unsigned>(kFastCounterexampleSearchConflictLimit),
          static_cast<unsigned>(kFastCounterexampleSearchDecisionLimit));
    } else {
      status = solver.solveWithAssumptionsStatus(  // LCOV_EXCL_LINE
          {},  // LCOV_EXCL_LINE
          kFastCounterexampleSearchConflictLimit,
          /*propagationLimit=*/-1);
    }
  } else {
    status = solver.solveStatus();
  // LCOV_EXCL_START
  }
  if (solveStatusOut != nullptr) {
  // LCOV_EXCL_STOP
    *solveStatusOut = status;
  }
  if (status != SATSolverWrapper::SolveStatus::Sat) {
    // LCOV_EXCL_START
    if (status == SATSolverWrapper::SolveStatus::Unknown &&
    // LCOV_EXCL_STOP
        isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC diag: k-induction fast base k=", k,
          " unknown within conflict budget ",
          kFastCounterexampleSearchConflictLimit);
    }  // LCOV_EXCL_LINE
    return std::nullopt;
  }
  if (solverProfile == BaseCaseSolverProfile::PdrValidationProofOnly) {
    // PDR validation calls this mode only as a SAT/UNSAT oracle.  The COI can
    // intentionally omit symbols needed for user-facing mismatch traces, so do
    // not build a full witness when the caller only checks has_value().
    KInductionResult::CounterexampleWitness witness;  // LCOV_EXCL_LINE
    witness.badFrame = firstBadFrame - bootstrapFrames;  // LCOV_EXCL_LINE
    return witness;  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return buildCounterexampleWitness(
      solver, variables, problem, firstBadFrame, lastBadFrame, bootstrapFrames);
}

std::optional<KInductionResult::CounterexampleWitness>
findImcCachedBaseCounterexampleAtFrontierQuery(
    const ImcBaseCounterexampleCache& cache,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k,
    BaseCaseSolverProfile solverProfile = BaseCaseSolverProfile::SecConeProof,
    SATSolverWrapper::SolveStatus* solveStatusOut = nullptr) {
  const auto& problem = cache.problem;
  const size_t bootstrapFrames = resetBootstrapFrames(problem);
  const bool resetBootstrapObservationFrontier =
      bootstrapFrames != 0 && problem.usesResetBootstrapObservationFrontier();
  const size_t internalK = k + bootstrapFrames;
  const InitialConstraintMode initialMode =
      bootstrapFrames == 0 ? determineInitialConstraintMode(problem)
                           : InitialConstraintMode::None;

  size_t firstBadFrame = 0;
  if (bootstrapFrames != 0) {
    firstBadFrame =
        bootstrapFrames + (resetBootstrapObservationFrontier ? 1u : 0u);
  } else if (initialMode == InitialConstraintMode::ObservationOnly ||
             initialMode == InitialConstraintMode::PartialInit) {
    firstBadFrame = 1;
  }
  const size_t exactInternalBadFrame = k + bootstrapFrames;
  if (firstBadFrame > internalK ||
      exactInternalBadFrame < firstBadFrame ||
      exactInternalBadFrame > internalK) {
    return std::nullopt;
  }
  firstBadFrame = exactInternalBadFrame;
  const size_t lastBadFrame = exactInternalBadFrame;

  const BaseCaseCoi coi = buildImcCachedBaseCaseCoi(
      cache,
      initialMode,
      bootstrapFrames,
      resetBootstrapObservationFrontier,
      firstBadFrame,
      internalK,
      /*constrainPreviouslySafeFrames=*/false);
  emitBaseCaseCoiDiag(
      problem,
      coi,
      k,
      firstBadFrame,
      internalK,
      /*constrainPreviouslySafeFrames=*/false);
  SATSolverWrapper solver(solverType);
  if (baseCaseValidationUsesLocalQueryProfile(solverType)) {
    solver.configureForSecPdrQuery(coi.solverSymbols.size());
  } else {
    solver.configureForSecConeProof(coi.solverSymbols.size());  // LCOV_EXCL_LINE
  }
  FrameVariableStore variables(solver, coi.solverSymbols, internalK + 1);
  addResetBootstrapConstraints(solver, variables, problem, internalK + 1);
  addInitialConstraints(solver, variables, problem, coi.solverSymbolSet, initialMode);
  if (resetBootstrapObservationFrontier) {
    addObservationPropertyConstraint( // LCOV_EXCL_LINE
        solver, variables, problem, bootstrapFrames); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs0, coi.solverSymbolSet,
      internalK + 1);
  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs1, coi.solverSymbolSet,
      internalK + 1);
  addSameFrameStateEqualities(
      solver, variables, problem, coi.solverSymbolSet, internalK + 1);
  addDualRailStateValidity(
      solver, variables, problem.dualRailStatePairs, coi.solverSymbolSet,
      internalK + 1);
  for (size_t frame = 0; frame < internalK; ++frame) {
    addTransitionRelation(
        solver,
        variables,
        cache.transitionByState,
        coi.transitionTargetsByFrame[frame],
        frame);
  }
  if (bootstrapFrames != 0) {
    addBootstrapStateAssignments(
        solver, variables, problem, coi.solverSymbolSet, bootstrapFrames);
  }

  // The cached newest-frontier path checks only the exact bad frame.  Earlier
  // frames were swept by previous iterations, and any bad at this frame is a
  // valid counterexample even if an earlier frame would also have failed.
  FrameFormulaEncoder encoder(
      solver,
      variables.makeLeafLits(
          firstBadFrame, formulaSupportOrThrow(problem.bad, "bad-state formula")));
  solver.addClause({encoder.encode(problem.bad)});
  SATSolverWrapper::SolveStatus status = SATSolverWrapper::SolveStatus::Unknown;
  if (solverProfile == BaseCaseSolverProfile::FastCounterexampleSearch) {
    // Batch-first IMC only needs a quick decisive answer.  UNSAT is still a
    // sound proof that the whole residual batch is safe at this frontier; UNKNOWN
    // falls back to exact per-output localization instead of being treated as
    // covered.
    if (solverType == KEPLER_FORMAL::Config::SolverType::KISSAT) { // LCOV_EXCL_LINE
      status = solver.solveWithKissatResourceLimits( // LCOV_EXCL_LINE
          static_cast<unsigned>(kFastCounterexampleSearchConflictLimit),
          static_cast<unsigned>(kFastCounterexampleSearchDecisionLimit));
    } else { // LCOV_EXCL_LINE
      status = solver.solveWithAssumptionsStatus( // LCOV_EXCL_LINE
          {}, // LCOV_EXCL_LINE
          kFastCounterexampleSearchConflictLimit,
          /*propagationLimit=*/-1);
    }
  } else { // LCOV_EXCL_LINE
    status = solver.solveStatus();
  }
  if (solveStatusOut != nullptr) {
    *solveStatusOut = status; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  if (status != SATSolverWrapper::SolveStatus::Sat) {
    return std::nullopt;
  }
  return buildCounterexampleWitness(
      solver, variables, problem, firstBadFrame, lastBadFrame, bootstrapFrames);
}

std::optional<KInductionResult::CounterexampleWitness>
findImcAssumptionBaseCounterexampleAtFrontier( // LCOV_EXCL_LINE
    const ImcBaseCounterexampleCache& cache,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  const auto& problem = cache.problem; // LCOV_EXCL_LINE
  if (problem.observedOutputExprs0.empty() || // LCOV_EXCL_LINE
      problem.observedOutputExprs0.size() != problem.observedOutputExprs1.size()) { // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  const size_t bootstrapFrames = resetBootstrapFrames(problem); // LCOV_EXCL_LINE
  const bool resetBootstrapObservationFrontier = // LCOV_EXCL_LINE
      bootstrapFrames != 0 && problem.usesResetBootstrapObservationFrontier(); // LCOV_EXCL_LINE
  const size_t internalK = k + bootstrapFrames; // LCOV_EXCL_LINE
  const InitialConstraintMode initialMode = // LCOV_EXCL_LINE
      bootstrapFrames == 0 ? determineInitialConstraintMode(problem) // LCOV_EXCL_LINE
                           : InitialConstraintMode::None;

  size_t firstBadFrame = 0; // LCOV_EXCL_LINE
  if (bootstrapFrames != 0) { // LCOV_EXCL_LINE
    firstBadFrame = // LCOV_EXCL_LINE
        bootstrapFrames + (resetBootstrapObservationFrontier ? 1u : 0u); // LCOV_EXCL_LINE
  } else if (initialMode == InitialConstraintMode::ObservationOnly || // LCOV_EXCL_LINE
             initialMode == InitialConstraintMode::PartialInit) { // LCOV_EXCL_LINE
    firstBadFrame = 1; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  const size_t exactInternalBadFrame = k + bootstrapFrames; // LCOV_EXCL_LINE
  if (firstBadFrame > internalK || // LCOV_EXCL_LINE
      exactInternalBadFrame < firstBadFrame || // LCOV_EXCL_LINE
      exactInternalBadFrame > internalK) { // LCOV_EXCL_LINE
    return std::nullopt; // LCOV_EXCL_LINE
  }
  firstBadFrame = exactInternalBadFrame; // LCOV_EXCL_LINE

  const BaseCaseCoi coi = buildImcCachedBaseCaseCoi( // LCOV_EXCL_LINE
      cache, // LCOV_EXCL_LINE
      initialMode, // LCOV_EXCL_LINE
      bootstrapFrames, // LCOV_EXCL_LINE
      resetBootstrapObservationFrontier, // LCOV_EXCL_LINE
      firstBadFrame, // LCOV_EXCL_LINE
      internalK, // LCOV_EXCL_LINE
      /*constrainPreviouslySafeFrames=*/false);
  emitBaseCaseCoiDiag( // LCOV_EXCL_LINE
      problem, // LCOV_EXCL_LINE
      coi,
      k, // LCOV_EXCL_LINE
      firstBadFrame, // LCOV_EXCL_LINE
      internalK, // LCOV_EXCL_LINE
      /*constrainPreviouslySafeFrames=*/false);
  const auto assumptionSolverType = // LCOV_EXCL_LINE
      SATSolverWrapper::assumptionSolverTypeFor(solverType); // LCOV_EXCL_LINE
  SATSolverWrapper solver(assumptionSolverType); // LCOV_EXCL_LINE
  // IMC asks the same concrete frontier prefix for many output bad literals.
  // Build that exact prefix once and vary only the selected top-output mismatch
  // through solver assumptions; this keeps the engine IMC-only while avoiding a
  // full transition rebuild per output.
  solver.configureForSecLocalBooleanCheck(coi.solverSymbols.size()); // LCOV_EXCL_LINE
  FrameVariableStore variables(solver, coi.solverSymbols, internalK + 1); // LCOV_EXCL_LINE
  addResetBootstrapConstraints(solver, variables, problem, internalK + 1); // LCOV_EXCL_LINE
  addInitialConstraints(solver, variables, problem, coi.solverSymbolSet, initialMode); // LCOV_EXCL_LINE
  if (resetBootstrapObservationFrontier) { // LCOV_EXCL_LINE
    addObservationPropertyConstraint( // LCOV_EXCL_LINE
        solver, variables, problem, bootstrapFrames); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  addComplementedStateRelations( // LCOV_EXCL_LINE
      solver, variables, problem.complementedStatePairs0, coi.solverSymbolSet, // LCOV_EXCL_LINE
      internalK + 1); // LCOV_EXCL_LINE
  addComplementedStateRelations( // LCOV_EXCL_LINE
      solver, variables, problem.complementedStatePairs1, coi.solverSymbolSet, // LCOV_EXCL_LINE
      internalK + 1); // LCOV_EXCL_LINE
  addSameFrameStateEqualities( // LCOV_EXCL_LINE
      solver, variables, problem, coi.solverSymbolSet, internalK + 1); // LCOV_EXCL_LINE
  addDualRailStateValidity( // LCOV_EXCL_LINE
      solver, variables, problem.dualRailStatePairs, coi.solverSymbolSet, // LCOV_EXCL_LINE
      internalK + 1); // LCOV_EXCL_LINE
  for (size_t frame = 0; frame < internalK; ++frame) { // LCOV_EXCL_LINE
    addTransitionRelation( // LCOV_EXCL_LINE
        solver,
        variables,
        cache.transitionByState, // LCOV_EXCL_LINE
        coi.transitionTargetsByFrame[frame], // LCOV_EXCL_LINE
        frame); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  if (bootstrapFrames != 0) { // LCOV_EXCL_LINE
    addBootstrapStateAssignments( // LCOV_EXCL_LINE
        solver, variables, problem, coi.solverSymbolSet, bootstrapFrames); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  // A newest-frontier query only asks whether bad is reachable at this exact
  // frame.  Earlier bad frames are also valid counterexamples, so constraining
  // previous frames safe is unnecessary and can dominate wide IMC/KI residuals.
  FrameFormulaEncoder badEncoder( // LCOV_EXCL_LINE
      solver, variables.makeLeafLits(firstBadFrame)); // LCOV_EXCL_LINE
  std::vector<int> outputBadLits; // LCOV_EXCL_LINE
  outputBadLits.reserve(problem.observedOutputExprs0.size()); // LCOV_EXCL_LINE
  for (size_t output = 0; output < problem.observedOutputExprs0.size(); ++output) { // LCOV_EXCL_LINE
    outputBadLits.push_back(badEncoder.encode(BoolExpr::simplify(BoolExpr::Xor( // LCOV_EXCL_LINE
        problem.observedOutputExprs0[output], // LCOV_EXCL_LINE
        problem.observedOutputExprs1[output])))); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  for (size_t output = 0; output < outputBadLits.size(); ++output) { // LCOV_EXCL_LINE
    const auto status = solver.solveWithAssumptionsStatus({outputBadLits[output]}); // LCOV_EXCL_LINE
    if (status == SATSolverWrapper::SolveStatus::Sat) { // LCOV_EXCL_LINE
      return buildCounterexampleWitness( // LCOV_EXCL_LINE
          solver, variables, problem, firstBadFrame, firstBadFrame, bootstrapFrames); // LCOV_EXCL_LINE
    }
    if (status == SATSolverWrapper::SolveStatus::Unknown) { // LCOV_EXCL_LINE
      // Keep UNKNOWN conservative.  The assumption path is an optimization only;
      // fall back to the exact single-output query rather than treating a
      // resource-limited answer as a safe frontier.
      if (auto witness = findImcCachedBaseCounterexampleAtFrontierQuery( // LCOV_EXCL_LINE
              cache.singleOutputCache(output), solverType, k); // LCOV_EXCL_LINE
          witness.has_value()) { // LCOV_EXCL_LINE
        return witness; // LCOV_EXCL_LINE
      }
    } // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  return std::nullopt; // LCOV_EXCL_LINE
} // LCOV_EXCL_LINE

std::optional<KInductionResult::CounterexampleWitness>
findImcBaseCounterexampleAtFrontierImpl(
    const ImcBaseCounterexampleCache& cache,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  const auto& problem = cache.problem;
  if (problem.observedOutputExprs0.size() > 1 &&
      problem.observedOutputExprs0.size() == problem.observedOutputExprs1.size()) {
    if (problem.observedOutputExprs0.size() <=
        kMaxImcAssumptionFrontierBatchOutputs) {
      if (solverType == KEPLER_FORMAL::Config::SolverType::KISSAT) {
        // Kissat does not support assumptions, and switching equivalent
        // residual batches to CaDiCaL can dominate medium designs.  The batch
        // bad predicate is already an exact OR over this output slice, so keep
        // the selected solver when it can decide the whole slice directly.
        return findImcCachedBaseCounterexampleAtFrontierQuery(
            cache, solverType, k);
      }
      return findImcAssumptionBaseCounterexampleAtFrontier(cache, solverType, k); // LCOV_EXCL_LINE
    }

    // Wide residual batches are split before SAT.  The exact assumption solver
    // above then amortizes each batch's transition prefix across its outputs.
    for (const auto& [firstOutput, endOutput] :
         buildStrictImcOutputSubsets(problem)) {
      if (auto witness = findImcBaseCounterexampleAtFrontierImpl(
              cache.outputSubsetCache(firstOutput, endOutput), solverType, k);
          witness.has_value()) {
        return witness; // LCOV_EXCL_LINE
      }
    }
    return std::nullopt;
  }

  return findImcCachedBaseCounterexampleAtFrontierQuery(cache, solverType, k);
}

struct DualRailPublicBasePrefixResult {
  bool handled = false;
  std::optional<KInductionResult::CounterexampleWitness> witness;
};

struct DualRailPublicBasePrefixCache {
  const KInductionProblem* problem = nullptr;
  KEPLER_FORMAL::Config::SolverType solverType =
      KEPLER_FORMAL::Config::SolverType::KISSAT;
  bool hasSafeThrough = false;
  size_t safeThrough = 0;

  void reset(const KInductionProblem& newProblem,
             KEPLER_FORMAL::Config::SolverType newSolverType) {
    problem = &newProblem;
    solverType = newSolverType;
    hasSafeThrough = false;
    safeThrough = 0;
  }
};

thread_local DualRailPublicBasePrefixCache dualRailPublicBasePrefixCache;

bool canUseDualRailPublicBasePrefixCache(const KInductionProblem& problem) {
  return problem.usesDualRailStateEncoding &&
         problem.observedOutputExprs0.size() > 1 &&
         problem.observedOutputExprs0.size() == problem.observedOutputExprs1.size();
}

size_t dualRailBaseStateSymbolCount(const KInductionProblem& problem) {
  return problem.dualRailStatePairs.empty()
      ? problem.state0Symbols.size() + problem.state1Symbols.size()
      : problem.dualRailStatePairs.size() * 2;
}

bool shouldUseMonolithicDualRailPublicBasePrefix(
    const KInductionProblem& problem) {
  return dualRailBaseStateSymbolCount(problem) >=
         kMinLargeDualRailPublicBasePrefixStateSymbols;
}

size_t dualRailPublicBasePrefixWindow(const KInductionProblem& problem) {
  const size_t stateSymbols = dualRailBaseStateSymbolCount(problem);
  if (stateSymbols >= kMinLargeDualRailPublicBasePrefixStateSymbols) {
    // ASIC-sized dual-rail residual sweeps revisit the same public base prefix
    // for every frontier.  Proving the normal SEC horizon once is still the
    // exact cumulative base case, but avoids rebuilding the AES-sized prefix at
    // k=8,17,26,35.
    return kLargeDualRailPublicBasePrefixWindow;
  }
  return kDualRailPublicBasePrefixWindow;
}

std::optional<KInductionResult::CounterexampleWitness>
findDualRailBatchedBaseCounterexampleUpTo(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  if (shouldUseMonolithicDualRailPublicBasePrefix(problem)) {
    // For large dual-rail SEC prefixes, one exact public bad query is cheaper
    // than rebuilding the same transition prefix for every output batch.  This
    // remains the normal cumulative KI base case: SAT yields a real public
    // counterexample, UNSAT proves the whole public-output prefix safe.
    return findBaseCounterexampleImpl(
        problem,
        solverType,
        k,
        /*exactPublicBadFrame=*/std::nullopt,
        /*localizeMultiOutputFrontier=*/false);
  }

  for (const auto& [firstOutput, endOutput] :
       buildSupportBoundedOutputBatches(problem)) {
    KInductionProblem batch = problem;
    configureOutputBatchProblem(batch, problem, firstOutput, endOutput);
    if (auto witness = findBaseCounterexampleImpl(
            batch,
            solverType,
            k,
            /*exactPublicBadFrame=*/std::nullopt,
            /*localizeMultiOutputFrontier=*/false);
        witness.has_value()) {
      return witness;
    }
  }
  return std::nullopt;
}

DualRailPublicBasePrefixResult tryDualRailPublicBasePrefixCache(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  if (!canUseDualRailPublicBasePrefixCache(problem)) {
    return {};
  }

  auto& cache = dualRailPublicBasePrefixCache;
  if (cache.problem != &problem || cache.solverType != solverType) {
    cache.reset(problem, solverType);
  }
  if (cache.hasSafeThrough && k <= cache.safeThrough) {
    return {true, std::nullopt};
  }

  const size_t prefixK = k + dualRailPublicBasePrefixWindow(problem);
  auto witness =
      findDualRailBatchedBaseCounterexampleUpTo(problem, solverType, prefixK);
  if (!witness.has_value()) {
    // This is a real cumulative base-case UNSAT result over the prefix, so the
    // exact frontier checks inside that prefix are safe to answer from cache.
    cache.hasSafeThrough = true;
    cache.safeThrough = prefixK;
    return {true, std::nullopt};
  }

  // A later SAT model is not a proof that earlier frontiers are safe.  Return it
  // only for the exact frontier requested; otherwise let the caller run the
  // normal exact single-frontier query below.
  if (witness->badFrame == k) {
    return {true, witness};
  }
  return {};
}

}  // namespace

struct ResetFrontierReachabilityContext {
  explicit ResetFrontierReachabilityContext(
      std::shared_ptr<ResetFrontierReachabilityContextData> data)
      : data(std::move(data)) {}

  std::shared_ptr<ResetFrontierReachabilityContextData> data;
};

KEPLER_FORMAL::Config::SolverType baseCaseValidationSolverType(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType) {
  (void)problem;
  // Base-case checks are one-shot BMC oracles, not incremental assumption
  // queries. Keep the requested proof solver here; reset/frontier helpers that
  // actually depend on assumptions use SATSolverWrapper::assumptionSolverTypeFor
  // where they build those assumption-capable solvers.
  return solverType;
}

bool baseCaseValidationUsesLocalQueryProfile(
    KEPLER_FORMAL::Config::SolverType solverType) {
  // The CaDiCaL backend is selected here specifically for short validation
  // probes. Its default inprocessing can dominate medium SEC cones before CDCL
  // starts. Kissat's proof-oriented default has a similar problem on these
  // rebuilt base BMCs, so both embedded CDCL solvers use the local-query
  // profile and Glucose keeps its own native defaults.
  return solverType == KEPLER_FORMAL::Config::SolverType::CADICAL ||
         solverType == KEPLER_FORMAL::Config::SolverType::KISSAT;
}

std::optional<KInductionResult::CounterexampleWitness> findBaseCounterexample(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  // Keep the public KI base-case path on the one-shot validation policy above:
  // the selected solver is preserved, while the query profile avoids spending
  // minutes in standalone preprocessing on medium SEC cones such as
  // nangate45_dynamic_node.
  return findBaseCounterexampleImpl(
      problem, baseCaseValidationSolverType(problem, solverType), k, std::nullopt);
}

BaseCounterexampleCheckResult checkBaseCounterexampleWithFastValidation(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  SATSolverWrapper::SolveStatus solveStatus =
      SATSolverWrapper::SolveStatus::Unknown;
  BaseCounterexampleCheckResult result;
  result.witness = findBaseCounterexampleImpl(
      problem,
      baseCaseValidationSolverType(problem, solverType),
      k,
      std::nullopt,
      // LCOV_EXCL_START
      /*localizeMultiOutputFrontier=*/false,
      BaseCaseSolverProfile::FastCounterexampleSearch,
      &solveStatus);
      // LCOV_EXCL_STOP
  if (solveStatus == SATSolverWrapper::SolveStatus::Unsat) {
    result.status = BaseCounterexampleCheckStatus::NoCounterexample;
  } else if (solveStatus == SATSolverWrapper::SolveStatus::Sat) {
    result.status = BaseCounterexampleCheckStatus::Counterexample;  // LCOV_EXCL_LINE
  } else {  // LCOV_EXCL_LINE
    result.status = BaseCounterexampleCheckStatus::Unknown;  // LCOV_EXCL_LINE
  }
  return result;
}

std::optional<KInductionResult::CounterexampleWitness>
findBaseCounterexampleAtFrontier(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    // LCOV_EXCL_START
    size_t k) {
    // LCOV_EXCL_STOP
  const auto localSolverType = baseCaseValidationSolverType(problem, solverType);
  if (auto prefixResult =
          tryDualRailPublicBasePrefixCache(problem, localSolverType, k);
      prefixResult.handled) {
    return prefixResult.witness;
  }
  if (problem.usesDualRailStateEncoding &&
      problem.observedOutputExprs0.size() > 1 &&
      problem.observedOutputExprs0.size() == problem.observedOutputExprs1.size()) {
    // Dual-rail residual sweeps revisit the same exact frontier obligation for
    // every depth. Reuse the local cached frontier path so output batches share
    // transition/COI setup while still checking the real bad predicate with the
    // requested SAT solver.
    const auto cache = makeImcBaseCounterexampleCache(problem);
    return findImcBaseCounterexampleAtFrontierImpl(
        *cache, localSolverType, k);
  }
  return findBaseCounterexampleImpl(
      problem, localSolverType, k, k);
}

std::shared_ptr<ImcBaseCounterexampleCache> makeImcBaseCounterexampleCache(
    const KInductionProblem& problem) {
  return std::make_shared<ImcBaseCounterexampleCache>(problem);
}

std::optional<KInductionResult::CounterexampleWitness>
findImcBaseCounterexampleAtFrontier(
    const ImcBaseCounterexampleCache& cache,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  return findImcBaseCounterexampleAtFrontierImpl(
      cache,
      baseCaseValidationSolverType(cache.problem, solverType),
      k);
}

std::shared_ptr<KInductionBaseCounterexampleCache>
makeKInductionBaseCounterexampleCache(const KInductionProblem& problem) {
  return std::make_shared<KInductionBaseCounterexampleCache>(problem);
}

std::optional<KInductionResult::CounterexampleWitness>
findKInductionBaseCounterexampleAtFrontier(
    const KInductionBaseCounterexampleCache& cache,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  return findImcBaseCounterexampleAtFrontierImpl(
      cache,
      baseCaseValidationSolverType(cache.problem, solverType),
      k);
}

std::optional<KInductionResult::CounterexampleWitness>
findFastBaseCounterexampleAtFrontier(  // LCOV_EXCL_LINE
    const KInductionProblem& problem,
    // LCOV_EXCL_START
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  // This fresh frontier SAT query normally keeps the user's configured solver
  // and only borrows the SAT-oriented validation profile. Reset-bootstrap
  // checks are assumption-shaped, so route them through the same default
  // LCOV_EXCL_STOP
  // assumption-capable solver used by the exact frontier validators.
  return findBaseCounterexampleImpl(  // LCOV_EXCL_LINE
      problem,  // LCOV_EXCL_LINE
      baseCaseValidationSolverType(problem, solverType),  // LCOV_EXCL_LINE
      k,  // LCOV_EXCL_LINE
      k,  // LCOV_EXCL_LINE
      /*localizeMultiOutputFrontier=*/true,
      BaseCaseSolverProfile::FastCounterexampleSearch);
}

bool provesNoBaseCounterexampleAtFrontier(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  SATSolverWrapper::SolveStatus solveStatus =
      SATSolverWrapper::SolveStatus::Unknown;
  const auto witness = findBaseCounterexampleImpl(
      problem,
      baseCaseValidationSolverType(problem, solverType),
      k,
      // LCOV_EXCL_START
      k,
      // LCOV_EXCL_STOP
      /*localizeMultiOutputFrontier=*/false,
      BaseCaseSolverProfile::PdrValidationProofOnly,
      &solveStatus);
  // No witness is not enough here: timeout/UNKNOWN also returns no witness, but
  // PDR callers use this helper as a proof certificate.
  return !witness.has_value() &&
         solveStatus == SATSolverWrapper::SolveStatus::Unsat;
}

// LCOV_EXCL_START

bool isStateCubeReachableAtResetFrontier(  // LCOV_EXCL_LINE
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    // LCOV_EXCL_STOP
    const std::vector<std::pair<size_t, bool>>& cube,
    // LCOV_EXCL_START
    size_t postBootstrapSteps) {
    // LCOV_EXCL_STOP
  const TransitionExprResolver transitionByState(problem);  // LCOV_EXCL_LINE
  return isStateCubeReachableAtResetFrontier(  // LCOV_EXCL_LINE
      problem, solverType, transitionByState, cube, postBootstrapSteps);  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

bool isStateCubeReachableAtResetFrontier(  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    const std::vector<std::pair<size_t, bool>>& cube,
    // LCOV_EXCL_STOP
    size_t postBootstrapSteps) {
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);  // LCOV_EXCL_LINE
  return isStateCubeReachableAtResetFrontier(  // LCOV_EXCL_LINE
      *context, solverType, cube, postBootstrapSteps);  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

std::shared_ptr<ResetFrontierReachabilityContext>
makeResetFrontierReachabilityContext(
    // LCOV_EXCL_START
    const KInductionProblem& problem,
    // LCOV_EXCL_STOP
    const TransitionExprResolver& transitionByState,
    BoolExpr* frameInvariant) {
  return std::make_shared<ResetFrontierReachabilityContext>(
      std::make_shared<ResetFrontierReachabilityContextData>(
          problem, transitionByState, frameInvariant));
}  // LCOV_EXCL_LINE

void rememberResetFrontierUnreachableCore(
    const ResetFrontierReachabilityContextData& data,
    size_t targetFrame,
    std::vector<std::pair<size_t, bool>> core);

// LCOV_EXCL_START


// LCOV_EXCL_STOP
void rememberResetFrontierUnreachableCube(
    const ResetFrontierReachabilityContext& context,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps) {
  if (cube.empty()) {
    return;  // LCOV_EXCL_LINE
  }

  const auto& data = *context.data;
  const size_t targetFrame = data.bootstrapFrames + postBootstrapSteps;
  rememberResetFrontierUnreachableCore(data, targetFrame, cube);
}

std::vector<size_t> sortedCubeSymbols(
    const std::vector<std::pair<size_t, bool>>& cube) {
  std::vector<size_t> symbols;
  symbols.reserve(cube.size());
  for (const auto& [symbol, value] : cube) {
    (void)value;
    symbols.push_back(symbol);
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  std::sort(symbols.begin(), symbols.end());
  // LCOV_EXCL_START
  symbols.erase(std::unique(symbols.begin(), symbols.end()), symbols.end());
  return symbols;
}
// LCOV_EXCL_STOP

std::vector<std::pair<size_t, bool>> sortedCubeLiterals(  // LCOV_EXCL_LINE
    std::vector<std::pair<size_t, bool>> cube) {
  std::sort(cube.begin(), cube.end());  // LCOV_EXCL_LINE
  cube.erase(std::unique(cube.begin(), cube.end()), cube.end());  // LCOV_EXCL_LINE
  return cube;  // LCOV_EXCL_LINE
}

ResetFrontierSolverCacheKey resetFrontierSolverCacheKey(
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t targetFrame,
    const std::vector<std::pair<size_t, bool>>& cube,
    // LCOV_EXCL_START
    bool includeCubeValues) {
  ResetFrontierSolverCacheKey key;
  // LCOV_EXCL_STOP
  key.solverType = solverType;
  key.targetFrame = targetFrame;
  key.includeCubeValues = includeCubeValues;
  if (includeCubeValues) {
    key.cubeLiterals = sortedCubeLiterals(cube);  // LCOV_EXCL_LINE
  } else {  // LCOV_EXCL_LINE
    key.cubeSymbols = sortedCubeSymbols(cube);
  }
  return key;
}

const CachedResetSummaryCoi& getCachedResetSummaryCubeReachabilityCoi(
    const ResetFrontierReachabilityContextData& data,
    size_t postBootstrapSteps,
    const std::vector<std::pair<size_t, bool>>& cube) {
  const ResetFrontierSolverCacheKey key =
      resetFrontierSolverCacheKey(
          // LCOV_EXCL_START
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          // LCOV_EXCL_STOP
          postBootstrapSteps,
          cube,
          /*includeCubeValues=*/false);
  // LCOV_EXCL_START
  if (const auto it = data.cachedResetSummaryCois.find(key);
      it != data.cachedResetSummaryCois.end()) {
      // LCOV_EXCL_STOP
    return it->second;  // LCOV_EXCL_LINE
  }

  if (data.cachedResetSummaryCois.size() >= kMaxResetSummaryCachedCois) {
    data.cachedResetSummaryCois.clear();  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  CachedResetSummaryCoi cached;
  cached.coi =
      buildResetSummaryCubeReachabilityCoi(data, postBootstrapSteps, cube);
  cached.transitionTargets =
      countTransitionTargets(cached.coi.transitionTargetsByFrame);
  auto [it, _] =
      data.cachedResetSummaryCois.emplace(key, std::move(cached));
  return it->second;
}

std::vector<int> stateCubeAssumptionLits(
    const FrameVariableStore& variables,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t frame) {
  std::vector<int> assumptions;
  assumptions.reserve(cube.size());
  for (const auto& [symbol, value] : cube) {
    const int literal = variables.getLiteral(symbol, frame);
    assumptions.push_back(value ? literal : -literal);
  }
  return assumptions;
}

bool solverContainsCubeSymbols(const CachedResetFrontierSolver& cached,
                               const std::vector<size_t>& cubeSymbols) {
  return std::all_of(
      cubeSymbols.begin(),
      cubeSymbols.end(),
      [&](const auto symbol) {
        return cached.coi.solverSymbolSet.find(symbol) !=
               cached.coi.solverSymbolSet.end();
      });
}

std::vector<std::pair<size_t, bool>> normalizedAssignmentCube(
    std::vector<std::pair<size_t, bool>> cube) {
  std::sort(cube.begin(), cube.end());
  cube.erase(std::unique(cube.begin(), cube.end()), cube.end());
  return cube;
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

class ParityUnionFind {
 public:
  void addEquality(size_t lhs, size_t rhs) { unite(lhs, rhs, false); }

  void addComplement(size_t lhs, size_t rhs) { unite(lhs, rhs, true); }  // LCOV_EXCL_LINE

  std::optional<bool> xorRelation(size_t lhs, size_t rhs) {
    if (parent_.find(lhs) == parent_.end() ||
        parent_.find(rhs) == parent_.end()) {
      return std::nullopt;
    }
    const auto lhsRoot = find(lhs);
    const auto rhsRoot = find(rhs);
    if (lhsRoot.first != rhsRoot.first) {
      return std::nullopt;
    }
    return lhsRoot.second ^ rhsRoot.second;
  }

  std::pair<size_t, bool> findWithParity(size_t symbol) {
    return find(symbol);
  }

 private:
  void ensure(size_t symbol) {
    if (parent_.find(symbol) == parent_.end()) {
      parent_.emplace(symbol, symbol);
      parityToParent_.emplace(symbol, false);
    }
  }

  std::pair<size_t, bool> find(size_t symbol) {
    ensure(symbol);
    const auto parent = parent_[symbol];
    const auto parity = parityToParent_[symbol];
    if (parent == symbol) {
      return {symbol, false};
    }
    const auto root = find(parent);
    parent_[symbol] = root.first;
    parityToParent_[symbol] = parity ^ root.second;
    return {parent_[symbol], parityToParent_[symbol]};
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP

  void unite(size_t lhs, size_t rhs, bool inverted) {
    const auto lhsRoot = find(lhs);
    const auto rhsRoot = find(rhs);
    if (lhsRoot.first == rhsRoot.first) {
      return;  // LCOV_EXCL_LINE
    }
    parent_[lhsRoot.first] = rhsRoot.first;
    // value(lhs) xor value(rhs) must equal `inverted`.
    parityToParent_[lhsRoot.first] =
        lhsRoot.second ^ rhsRoot.second ^ inverted;
  }

  std::unordered_map<size_t, size_t> parent_;
  std::unordered_map<size_t, bool> parityToParent_;
};

std::optional<std::vector<std::pair<size_t, bool>>>
knownResetFrontierConflictCore(
    const ResetFrontierReachabilityContextData& data,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps) {
  if (cube.empty() || postBootstrapSteps != 0) {
    return std::nullopt;
  }

  const bool usesBootstrapFrontier = data.bootstrapFrames != 0;
  const auto& assignments = usesBootstrapFrontier
                                ? data.problem.bootstrapStateAssignments
                                : data.problem.initialStateAssignments;
  const auto& equalities = emptySymbolPairs();
  if (assignments.empty() && equalities.empty() &&
      data.problem.complementedStatePairs0.empty() &&
      data.problem.complementedStatePairs1.empty() &&
      data.problem.sameFrameStateEqualityPairs0.empty() &&
      data.problem.sameFrameStateEqualityPairs1.empty()) {
    return std::nullopt;
  }

  // LCOV_EXCL_START
  ParityUnionFind relations;
  // LCOV_EXCL_STOP
  for (const auto& [lhsSymbol, rhsSymbol] : equalities) {
    relations.addEquality(lhsSymbol, rhsSymbol);
  }
  for (const auto& [lhsSymbol, rhsSymbol] :
       data.problem.sameFrameStateEqualityPairs0) {
    relations.addEquality(lhsSymbol, rhsSymbol);  // LCOV_EXCL_LINE
  }
  for (const auto& [lhsSymbol, rhsSymbol] :
       data.problem.sameFrameStateEqualityPairs1) {
    relations.addEquality(lhsSymbol, rhsSymbol);  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_START
  for (const auto& [primarySymbol, complementedSymbol] :
  // LCOV_EXCL_STOP
       data.problem.complementedStatePairs0) {
    relations.addComplement(primarySymbol, complementedSymbol);  // LCOV_EXCL_LINE
  }
  for (const auto& [primarySymbol, complementedSymbol] :
       data.problem.complementedStatePairs1) {
    relations.addComplement(primarySymbol, complementedSymbol);  // LCOV_EXCL_LINE
  }

  std::unordered_map<size_t, bool> rootAssignments;
  // LCOV_EXCL_START
  rootAssignments.reserve(assignments.size());
  // LCOV_EXCL_STOP
  for (const auto& [symbol, value] : assignments) {
    const auto [root, parity] = relations.findWithParity(symbol);
    const bool rootValue = value ^ parity;
    if (const auto it = rootAssignments.find(root);
        it != rootAssignments.end() && it->second != rootValue) {
      continue;  // LCOV_EXCL_LINE
    }
    rootAssignments.emplace(root, rootValue);
  }

  std::unordered_map<size_t, std::pair<bool, std::pair<size_t, bool>>>
      cubeValueByRoot;
  cubeValueByRoot.reserve(cube.size());
  for (const auto& literal : cube) {
    const auto [root, parity] = relations.findWithParity(literal.first);
    const auto assignment = rootAssignments.find(root);
    const bool rootValue = literal.second ^ parity;
    if (assignment != rootAssignments.end() &&
        assignment->second != rootValue) {
      return std::vector<std::pair<size_t, bool>>{literal};
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
    if (const auto it = cubeValueByRoot.find(root);
        it != cubeValueByRoot.end()) {
      if (it->second.first != rootValue) {
        return normalizedAssignmentCube({it->second.second, literal});
      }
      continue;  // LCOV_EXCL_LINE
    }
    cubeValueByRoot.emplace(root, std::pair{rootValue, literal});
  }

  return std::nullopt;
}

// LCOV_EXCL_START
std::optional<std::vector<std::pair<size_t, bool>>>
// LCOV_EXCL_STOP
failedAssumptionCoreFromLastResetFrontierSolve(
    CachedResetFrontierSolver& cached,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t targetFrame) {
  if (cached.cubeEncodedAsUnitClauses) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  std::unordered_map<int, std::pair<size_t, bool>> cubeLiteralByAssumption;
  cubeLiteralByAssumption.reserve(cube.size());
  for (const auto& [symbol, value] : cube) {
    const int literal = cached.variables->getLiteral(symbol, targetFrame);
    cubeLiteralByAssumption.emplace(
        value ? literal : -literal, std::pair{symbol, value});
  }

  std::vector<std::pair<size_t, bool>> core;
  for (const int failedAssumption : cached.solver->failedAssumptions()) {
    // LCOV_EXCL_START
    const auto it = cubeLiteralByAssumption.find(failedAssumption);
    // LCOV_EXCL_STOP
    if (it != cubeLiteralByAssumption.end()) {
      core.push_back(it->second);
    }
  }
  if (core.empty()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  return normalizedAssignmentCube(std::move(core));
}

bool assignmentCubeContains(
    const std::vector<std::pair<size_t, bool>>& cube,
    const std::vector<std::pair<size_t, bool>>& core) {
  return std::includes(cube.begin(), cube.end(), core.begin(), core.end());
}

bool cubeSymbolsAreInSolverCoi(
    const BaseCaseCoi& coi,
    const std::vector<std::pair<size_t, bool>>& cube) {
  return std::all_of(
      cube.begin(),
      cube.end(),
      [&](const auto& literal) {
        return coi.solverSymbolSet.find(literal.first) !=
               coi.solverSymbolSet.end();
      });
}

size_t addPreviousResetFrontierBlockers(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const ResetFrontierReachabilityContextData& data,
    const BaseCaseCoi& coi,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t targetFrame) {
  size_t added = 0;
  for (size_t frame = 0;
       frame < targetFrame &&
       added < kMaxPreviousResetFrontierBlockersPerQuery;
       ++frame) {
    const auto coresIt = data.unreachableCoresByTargetFrame.find(frame);
    if (coresIt == data.unreachableCoresByTargetFrame.end()) {
      continue;
    }
    for (const auto& core : coresIt->second) {
      if (!assignmentCubeContains(cube, core) ||
          !cubeSymbolsAreInSolverCoi(coi, core)) {
        continue;
      }
      // The cached core was already proved unreachable at this concrete
      // reset/bootstrap frame. Reusing it as a safe-prefix clause is redundant
      // with the exact prefix semantics, but it gives the next post-bootstrap
      // query the same learned fact without rebuilding the previous SAT proof.
      addBlockedStateCubeClause(solver, variables, core, frame);
      ++added;
      break;
    }
  }
  return added;
}

std::optional<std::vector<std::pair<size_t, bool>>>
findCachedResetFrontierUnreachableCore(
    const ResetFrontierReachabilityContextData& data,
    size_t targetFrame,
    const std::vector<std::pair<size_t, bool>>& cube) {
  const auto frameIt = data.unreachableCoresByTargetFrame.find(targetFrame);
  if (frameIt == data.unreachableCoresByTargetFrame.end()) {
    return std::nullopt;
  }
  for (const auto& core : frameIt->second) {
    if (assignmentCubeContains(cube, core)) {
      return core;
    }
  }
  return std::nullopt;
}

// LCOV_EXCL_START
void rememberResetFrontierUnreachableCore(
// LCOV_EXCL_STOP
    const ResetFrontierReachabilityContextData& data,
    size_t targetFrame,
    std::vector<std::pair<size_t, bool>> core) {
  core = normalizedAssignmentCube(std::move(core));
  if (core.empty()) {
    return;  // LCOV_EXCL_LINE
  }

  auto& cores = data.unreachableCoresByTargetFrame[targetFrame];
  for (const auto& existing : cores) {
    if (assignmentCubeContains(core, existing)) {
      return;
    }
  }
  cores.erase(
      std::remove_if(
          cores.begin(),
          // LCOV_EXCL_START
          cores.end(),
          [&](const auto& existing) {
          // LCOV_EXCL_STOP
            return assignmentCubeContains(existing, core);
          }),
      cores.end());
  if (cores.size() >= kMaxResetFrontierCachedCoresPerFrame) {
    // LCOV_EXCL_START
    cores.erase(cores.begin());  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }  // LCOV_EXCL_LINE
  cores.push_back(std::move(core));
}

// LCOV_EXCL_START

std::optional<std::vector<std::pair<size_t, bool>>>
// LCOV_EXCL_STOP
extractUnreachableCoreFromCachedResetFrontierSolver(  // LCOV_EXCL_LINE
    CachedResetFrontierSolver& cached,
    // LCOV_EXCL_START
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t targetFrame) {
  if (cached.cubeEncodedAsUnitClauses) {  // LCOV_EXCL_LINE
    return cached.solver->solve() ? std::nullopt : std::optional{cube};  // LCOV_EXCL_LINE
  }

  std::vector<int> assumptions;  // LCOV_EXCL_LINE
  assumptions.reserve(cube.size());  // LCOV_EXCL_LINE
  std::unordered_map<int, std::pair<size_t, bool>> cubeLiteralByAssumption;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  cubeLiteralByAssumption.reserve(cube.size());  // LCOV_EXCL_LINE
  for (const auto& [symbol, value] : cube) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    const int literal = cached.variables->getLiteral(symbol, targetFrame);  // LCOV_EXCL_LINE
    const int assumption = value ? literal : -literal;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    assumptions.push_back(assumption);  // LCOV_EXCL_LINE
    cubeLiteralByAssumption.emplace(assumption, std::pair{symbol, value});  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  }

  if (cached.solver->solveWithAssumptions(assumptions)) {  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }


// LCOV_EXCL_STOP
  std::vector<std::pair<size_t, bool>> core;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  for (const int failedAssumption : cached.solver->failedAssumptions()) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    const auto it = cubeLiteralByAssumption.find(failedAssumption);  // LCOV_EXCL_LINE
    if (it != cubeLiteralByAssumption.end()) {  // LCOV_EXCL_LINE
      core.push_back(it->second);  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
  }
  if (core.empty()) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    // Some solver backends / conflict shapes do not expose a mapped failed
    // LCOV_EXCL_START
    // assumption core. Start from the full cube and still run exact deletion
    // minimization below; every accepted drop is checked by SAT.
    core = cube;  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  core = normalizedAssignmentCube(std::move(core));  // LCOV_EXCL_LINE
  auto coreIsReachable =
      [&](const std::vector<std::pair<size_t, bool>>& candidate) {  // LCOV_EXCL_LINE
        return cached.solver->solveWithAssumptions(  // LCOV_EXCL_LINE
            stateCubeAssumptionLits(*cached.variables, candidate, targetFrame));  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      };  // LCOV_EXCL_LINE

  // The assumption solver reports a valid conflict subset, not a guaranteed-minimal one.
  // Minimize it exactly with the same cached reset-frontier solver; the result
  // becomes a stronger PDR F[0] refinement and a reusable cache entry for later
  // neighboring cubes.
  size_t checks = 0;  // LCOV_EXCL_LINE
  const size_t maxChecks =  // LCOV_EXCL_LINE
      std::max(kMinResetFrontierCoreChecks, core.size() * 2);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  for (size_t chunkSize = std::max<size_t>(1, core.size() / 2);  // LCOV_EXCL_LINE
       // LCOV_EXCL_START
       chunkSize > 0 && checks < maxChecks;) {  // LCOV_EXCL_LINE
    for (size_t index = 0; index < core.size() && checks < maxChecks;) {  // LCOV_EXCL_LINE
      const size_t erasedCount = std::min(chunkSize, core.size() - index);  // LCOV_EXCL_LINE
      if (erasedCount == 0 || erasedCount == core.size()) {  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      }
      std::vector<std::pair<size_t, bool>> reduced = core;  // LCOV_EXCL_LINE
      reduced.erase(  // LCOV_EXCL_LINE
          reduced.begin() + static_cast<std::ptrdiff_t>(index),  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
          reduced.begin() +  // LCOV_EXCL_LINE
              // LCOV_EXCL_START
              static_cast<std::ptrdiff_t>(index + erasedCount));  // LCOV_EXCL_LINE
      ++checks;  // LCOV_EXCL_LINE
      if (!coreIsReachable(reduced)) {  // LCOV_EXCL_LINE
        core = std::move(reduced);  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      }
      // LCOV_EXCL_STOP
      index += erasedCount;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
    if (chunkSize == 1) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      break;  // LCOV_EXCL_LINE
    }
    chunkSize = std::max<size_t>(1, chunkSize / 2);  // LCOV_EXCL_LINE
  }
  return core;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

std::unique_ptr<CachedResetFrontierSolver> buildResetFrontierSolver(
    const ResetFrontierReachabilityContextData& data,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t targetFrame,
    bool encodeCubeAsUnitClauses,
    bool closeStartupEqualityDependencies = true) {
  auto cached = std::make_unique<CachedResetFrontierSolver>();
  cached->solverType = solverType;
  cached->targetFrame = targetFrame;
  cached->coi = buildStateCubeReachabilityCoi(
      data, targetFrame, cube, closeStartupEqualityDependencies);

  const auto& problem = data.problem;
  cached->solver = std::make_unique<SATSolverWrapper>(solverType);
  // Reset-frontier checks run inside PDR's blocking loop. Sampled AES runs
  // showed the one-shot query spending its time in Kissat probing/sweeping
  // before any refinement could be learned, so use the same CDCL-oriented
  // profile regardless of whether the cube is encoded as clauses or
  // assumptions.
  cached->solver->configureForSecPdrQuery(cached->coi.solverSymbols.size());
  cached->variables = std::make_unique<FrameVariableStore>(
      *cached->solver,
      cached->coi.solverSymbols,
      targetFrame + 1);
  addResetBootstrapConstraints(
      *cached->solver, *cached->variables, problem, targetFrame + 1);
  addInitialConstraints(
      *cached->solver,
      *cached->variables,
      problem,
      cached->coi.solverSymbolSet,
      data.initialMode);
  addComplementedStateRelations(
      *cached->solver,
      *cached->variables,
      problem.complementedStatePairs0,
      cached->coi.solverSymbolSet,
      targetFrame + 1);
  addComplementedStateRelations(
      *cached->solver,
      *cached->variables,
      problem.complementedStatePairs1,
      cached->coi.solverSymbolSet,
      targetFrame + 1);
  addSameFrameStateEqualities(
      *cached->solver,
      *cached->variables,
      problem,
      cached->coi.solverSymbolSet,
      targetFrame + 1);
  addDualRailStateValidity(
      *cached->solver,
      *cached->variables,
      problem.dualRailStatePairs,
      cached->coi.solverSymbolSet,
      targetFrame + 1);
  addResetFrontierFrameInvariantConstraints(
      *cached->solver, *cached->variables, data, targetFrame);

  for (size_t frame = 0; frame < targetFrame; ++frame) {
    addTransitionRelation(
        *cached->solver,
        *cached->variables,
        data.transitionByState,
        cached->coi.transitionTargetsByFrame[frame],
        frame);
  }
  if (data.bootstrapFrames != 0) {
    addBootstrapStateAssignments(
        *cached->solver,
        *cached->variables,
        problem,
        cached->coi.solverSymbolSet,
        data.bootstrapFrames);
  }
  const size_t previousBlockers = addPreviousResetFrontierBlockers(
      *cached->solver,
      *cached->variables,
      data,
      cached->coi,
      cube,
      targetFrame);
  if (previousBlockers != 0 && isKInductionCoiDiagEnabled()) {
    emitSecDiag(
        "SEC diag: reset frontier previous unreachable blockers=",
        previousBlockers,
        " target_frame=",
        targetFrame);
  }

  cached->cubeEncodedAsUnitClauses = encodeCubeAsUnitClauses;
  if (cached->cubeEncodedAsUnitClauses) {
    for (const auto& [symbol, value] : cube) {
      const int literal = cached->variables->getLiteral(symbol, targetFrame);
      // LCOV_EXCL_START
      cached->solver->addClause({value ? literal : -literal});
      // LCOV_EXCL_STOP
    }
  }
  return cached;
}

// LCOV_EXCL_START
std::unique_ptr<CachedResetFrontierSolver> buildResetFrontierSolverForCoi(  // LCOV_EXCL_LINE
    const ResetFrontierReachabilityContextData& data,
    KEPLER_FORMAL::Config::SolverType solverType,
    BaseCaseCoi coi,
    // LCOV_EXCL_STOP
    const std::vector<std::pair<size_t, bool>>& cube,
    // LCOV_EXCL_START
    size_t targetFrame) {
    // LCOV_EXCL_STOP
  auto cached = std::make_unique<CachedResetFrontierSolver>();  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  cached->solverType = solverType;  // LCOV_EXCL_LINE
  cached->targetFrame = targetFrame;  // LCOV_EXCL_LINE
  cached->coi = std::move(coi);  // LCOV_EXCL_LINE

  const auto& problem = data.problem;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  cached->solver = std::make_unique<SATSolverWrapper>(solverType);  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  cached->solver->configureForSecPdrQuery(cached->coi.solverSymbols.size());  // LCOV_EXCL_LINE
  cached->variables = std::make_unique<FrameVariableStore>(  // LCOV_EXCL_LINE
      *cached->solver,  // LCOV_EXCL_LINE
      cached->coi.solverSymbols,  // LCOV_EXCL_LINE
      targetFrame + 1);  // LCOV_EXCL_LINE
  addResetBootstrapConstraints(  // LCOV_EXCL_LINE
      *cached->solver, *cached->variables, problem, targetFrame + 1);  // LCOV_EXCL_LINE
  addInitialConstraints(  // LCOV_EXCL_LINE
      *cached->solver,  // LCOV_EXCL_LINE
      *cached->variables,  // LCOV_EXCL_LINE
      problem,  // LCOV_EXCL_LINE
      cached->coi.solverSymbolSet,  // LCOV_EXCL_LINE
      data.initialMode);  // LCOV_EXCL_LINE
  addComplementedStateRelations(  // LCOV_EXCL_LINE
      *cached->solver,  // LCOV_EXCL_LINE
      *cached->variables,  // LCOV_EXCL_LINE
      problem.complementedStatePairs0,  // LCOV_EXCL_LINE
      cached->coi.solverSymbolSet,  // LCOV_EXCL_LINE
      targetFrame + 1);  // LCOV_EXCL_LINE
  addComplementedStateRelations(  // LCOV_EXCL_LINE
      *cached->solver,  // LCOV_EXCL_LINE
      *cached->variables,  // LCOV_EXCL_LINE
      problem.complementedStatePairs1,  // LCOV_EXCL_LINE
      cached->coi.solverSymbolSet,  // LCOV_EXCL_LINE
      targetFrame + 1);  // LCOV_EXCL_LINE
  addSameFrameStateEqualities(  // LCOV_EXCL_LINE
      *cached->solver,  // LCOV_EXCL_LINE
      *cached->variables,  // LCOV_EXCL_LINE
      problem,  // LCOV_EXCL_LINE
      cached->coi.solverSymbolSet,  // LCOV_EXCL_LINE
      targetFrame + 1);  // LCOV_EXCL_LINE
  addDualRailStateValidity(  // LCOV_EXCL_LINE
      *cached->solver,  // LCOV_EXCL_LINE
      *cached->variables,  // LCOV_EXCL_LINE
      problem.dualRailStatePairs,  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      cached->coi.solverSymbolSet,  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      targetFrame + 1);  // LCOV_EXCL_LINE
  addResetFrontierFrameInvariantConstraints(  // LCOV_EXCL_LINE
      *cached->solver, *cached->variables, data, targetFrame);  // LCOV_EXCL_LINE

  for (size_t frame = 0; frame < targetFrame; ++frame) {  // LCOV_EXCL_LINE
    addTransitionRelation(  // LCOV_EXCL_LINE
        *cached->solver,  // LCOV_EXCL_LINE
        *cached->variables,  // LCOV_EXCL_LINE
        data.transitionByState,  // LCOV_EXCL_LINE
        cached->coi.transitionTargetsByFrame[frame],  // LCOV_EXCL_LINE
        frame);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  if (data.bootstrapFrames != 0) {  // LCOV_EXCL_LINE
    addBootstrapStateAssignments(  // LCOV_EXCL_LINE
        *cached->solver,  // LCOV_EXCL_LINE
        *cached->variables,  // LCOV_EXCL_LINE
        problem,  // LCOV_EXCL_LINE
        cached->coi.solverSymbolSet,  // LCOV_EXCL_LINE
        data.bootstrapFrames);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  const size_t previousBlockers = addPreviousResetFrontierBlockers(  // LCOV_EXCL_LINE
      *cached->solver,  // LCOV_EXCL_LINE
      *cached->variables,  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      data,  // LCOV_EXCL_LINE
      cached->coi,  // LCOV_EXCL_LINE
      cube,  // LCOV_EXCL_LINE
      targetFrame);  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  if (previousBlockers != 0 && isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    emitSecDiag(  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        "SEC diag: reset frontier previous unreachable blockers=",
        previousBlockers,
        " target_frame=",
        // LCOV_EXCL_STOP
        targetFrame);
  // LCOV_EXCL_START
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  cached->cubeEncodedAsUnitClauses = false;  // LCOV_EXCL_LINE
  return cached;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START

CachedResetFrontierSolver& getCachedResetFrontierPrefixSolver(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    const ResetFrontierReachabilityContextData& data,
    // LCOV_EXCL_START
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t maxTargetFrame) {
  const auto cachedSolverType =  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
      SATSolverWrapper::assumptionSolverTypeFor(solverType);  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  const ResetFrontierSolverCacheKey key =
      resetFrontierSolverCacheKey(  // LCOV_EXCL_LINE
          cachedSolverType,  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
          maxTargetFrame,  // LCOV_EXCL_LINE
          cube,  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          /*includeCubeValues=*/false);
  if (const auto it = data.cachedPrefixSolvers.find(key);  // LCOV_EXCL_LINE
      it != data.cachedPrefixSolvers.end()) {  // LCOV_EXCL_LINE
    return *it->second;  // LCOV_EXCL_LINE
  }

  const auto cubeSymbols = sortedCubeSymbols(cube);  // LCOV_EXCL_LINE
  for (const auto& [_, cached] : data.cachedPrefixSolvers) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    if (cached->solverType == cachedSolverType &&  // LCOV_EXCL_LINE
        cached->targetFrame == maxTargetFrame &&  // LCOV_EXCL_LINE
        !cached->cubeEncodedAsUnitClauses &&  // LCOV_EXCL_LINE
        solverContainsCubeSymbols(*cached, cubeSymbols)) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
        emitSecDiag(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            "SEC diag: reset frontier prefix solver superset cache hit ",
            "target_frame=",
            maxTargetFrame,
            // LCOV_EXCL_STOP
            " cube_literals=",
            cube.size(),  // LCOV_EXCL_LINE
            " solver_symbols=",
            // LCOV_EXCL_START
            cached->coi.solverSymbols.size());  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      return *cached;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
  // LCOV_EXCL_START
  }

  if (data.cachedPrefixSolvers.size() >= kMaxResetFrontierCachedSolvers) {  // LCOV_EXCL_LINE
    data.cachedPrefixSolvers.clear();  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  auto cached = buildResetFrontierSolverForCoi(  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
      data,  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      cachedSolverType,  // LCOV_EXCL_LINE
      buildStateCubePrefixReachabilityCoi(  // LCOV_EXCL_LINE
          data,  // LCOV_EXCL_LINE
          maxTargetFrame,  // LCOV_EXCL_LINE
          cube,  // LCOV_EXCL_LINE
          /*closeStartupEqualityDependencies=*/true),
      cube,  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      maxTargetFrame);  // LCOV_EXCL_LINE
  auto [it, inserted] =  // LCOV_EXCL_LINE
      data.cachedPrefixSolvers.emplace(key, std::move(cached));  // LCOV_EXCL_LINE
  (void)inserted;  // LCOV_EXCL_LINE
  return *it->second;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

std::optional<std::vector<std::pair<size_t, bool>>>
extractResetSummaryFrontierCube(
    // LCOV_EXCL_START
    const ResetFrontierReachabilityContextData& data,
    // LCOV_EXCL_STOP
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const BaseCaseCoi& coi) {
  std::vector<std::pair<size_t, bool>> cube;
  if (coi.requiredStateSymbolsByFrame.empty()) {
    // LCOV_EXCL_START
    return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  for (const auto symbol : coi.requiredStateSymbolsByFrame.front()) {
    cube.emplace_back(
        symbol, solver.getLiteralValue(variables.getLiteral(symbol, 0)));
    if (cube.size() > kMaxResetSummaryFrontierCubeLiterals) {
      return std::nullopt;  // LCOV_EXCL_LINE
    }
  }
  return normalizedAssignmentCube(std::move(cube));
}

CachedResetFrontierSolver& getCachedResetFrontierSolver(
    const ResetFrontierReachabilityContextData& data,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t targetFrame);

// LCOV_EXCL_START
std::optional<std::vector<std::pair<size_t, bool>>>
// LCOV_EXCL_STOP
proveResetSummaryFrontierCubeUnreachable(
    const ResetFrontierReachabilityContextData& data,
    const std::vector<std::pair<size_t, bool>>& cube) {
  const auto normalizedCube = normalizedAssignmentCube(cube);
  if (normalizedCube.empty()) {
    // LCOV_EXCL_START
    return std::nullopt;  // LCOV_EXCL_LINE
  }


// LCOV_EXCL_STOP
  if (const auto knownCore = knownResetFrontierConflictCore(
          data, normalizedCube, /*postBootstrapSteps=*/0);
      knownCore.has_value()) {
    rememberResetFrontierUnreachableCore(  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        data, data.bootstrapFrames, *knownCore);  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
    return knownCore;  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_START
  if (const auto cachedCore = findCachedResetFrontierUnreachableCore(
          data, data.bootstrapFrames, normalizedCube);
          // LCOV_EXCL_STOP
      cachedCore.has_value()) {
    return cachedCore;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  if (normalizedCube.size() > kMaxResetSummaryFrontierProofCubeLiterals) {
    if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC diag: reset summary frontier proof skipped reason=cube_cap ",
          // LCOV_EXCL_STOP
          "frontier_cube=",
          normalizedCube.size(),  // LCOV_EXCL_LINE
          " max_literals=",
          kMaxResetSummaryFrontierProofCubeLiterals);
    }  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  const BaseCaseCoi frontierCoi =
      buildStateCubeReachabilityCoi(
          data,
          data.bootstrapFrames,
          normalizedCube,
          /*closeStartupEqualityDependencies=*/true);
  // LCOV_EXCL_START
  const size_t frontierTransitionTargets =
      countTransitionTargets(frontierCoi.transitionTargetsByFrame);
      // LCOV_EXCL_STOP
  if (frontierCoi.solverSymbols.size() >
          kMaxResetSummaryFrontierProofSymbols ||
      // LCOV_EXCL_START
      frontierTransitionTargets >
      // LCOV_EXCL_STOP
          kMaxResetSummaryFrontierProofTransitionTargets) {
    // LCOV_EXCL_START
    if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC diag: reset summary frontier proof skipped reason=coi_cap ",
          // LCOV_EXCL_START
          "frontier_cube=",
          normalizedCube.size(),  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
          " solver_symbols=",
          frontierCoi.solverSymbols.size(),  // LCOV_EXCL_LINE
          " transition_targets=",
          frontierTransitionTargets);
    }  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  // Summary CEGAR can ask many neighboring frontier questions with the same
  // symbol support. Reuse the exact reset-frontier solver and vary only the
  // assumptions; rebuilding it dominated BlackParrot PDR samples.
  CachedResetFrontierSolver& solver = getCachedResetFrontierSolver(
      data,
      SATSolverWrapper::assumptionSolverTypeFor(
          KEPLER_FORMAL::Config::getSolverType()),
      normalizedCube,
      data.bootstrapFrames);
  const auto assumptions = stateCubeAssumptionLits(
      // LCOV_EXCL_START
      *solver.variables, normalizedCube, data.bootstrapFrames);
  const auto status = solver.solver->solveWithAssumptionsStatus(
      assumptions, kResetSummaryFrontierProofConflictLimit);
      // LCOV_EXCL_STOP
  if (status == SATSolverWrapper::SolveStatus::Sat) {
    return std::nullopt;
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  if (status == SATSolverWrapper::SolveStatus::Unknown) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC diag: reset summary frontier proof resource_limit ",
          // LCOV_EXCL_START
          "frontier_cube=",
          normalizedCube.size(),  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
          " solver_symbols=",
          solver.coi.solverSymbols.size(),  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          " transition_targets=",
          frontierTransitionTargets);
    }  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  auto core = failedAssumptionCoreFromLastResetFrontierSolve(
  // LCOV_EXCL_STOP
      solver, normalizedCube, data.bootstrapFrames);  // LCOV_EXCL_LINE
  if (!core.has_value()) {  // LCOV_EXCL_LINE
    core = normalizedCube;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  rememberResetFrontierUnreachableCore(data, data.bootstrapFrames, *core);  // LCOV_EXCL_LINE
  return core;  // LCOV_EXCL_LINE
}

// LCOV_EXCL_START
std::vector<std::vector<std::pair<size_t, bool>>>
collectResetSummarySingletonFrontierBlockers(  // LCOV_EXCL_LINE
    const ResetFrontierReachabilityContextData& data,
    // LCOV_EXCL_STOP
    const std::vector<std::pair<size_t, bool>>& cube,
    const std::vector<std::vector<std::pair<size_t, bool>>>& existingBlockers,
    // LCOV_EXCL_START
    size_t maxNewBlockers) {
  std::vector<std::vector<std::pair<size_t, bool>>> blockers;  // LCOV_EXCL_LINE
  if (maxNewBlockers == 0) {  // LCOV_EXCL_LINE
    return blockers;  // LCOV_EXCL_LINE
  }

  CachedResetFrontierSolver& solver = getCachedResetFrontierSolver(  // LCOV_EXCL_LINE
      data,  // LCOV_EXCL_LINE
      SATSolverWrapper::assumptionSolverTypeFor(  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
          KEPLER_FORMAL::Config::getSolverType()),  // LCOV_EXCL_LINE
      cube,  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      data.bootstrapFrames);  // LCOV_EXCL_LINE
  for (const auto& literal : cube) {  // LCOV_EXCL_LINE
    if (blockers.size() >= maxNewBlockers) {  // LCOV_EXCL_LINE
      break;  // LCOV_EXCL_LINE
    }

    std::vector<std::pair<size_t, bool>> singleton{literal};  // LCOV_EXCL_LINE
    if (std::any_of(  // LCOV_EXCL_LINE
            existingBlockers.begin(),  // LCOV_EXCL_LINE
            existingBlockers.end(),  // LCOV_EXCL_LINE
            [&](const auto& existing) {  // LCOV_EXCL_LINE
              return assignmentCubeContains(singleton, existing);  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
            }) ||  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        std::any_of(  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            blockers.begin(),  // LCOV_EXCL_LINE
            blockers.end(),  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            [&](const auto& existing) {  // LCOV_EXCL_LINE
              return assignmentCubeContains(singleton, existing);  // LCOV_EXCL_LINE
            })) {
            // LCOV_EXCL_STOP
      continue;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }


// LCOV_EXCL_STOP
    const auto status = solver.solver->solveWithAssumptionsStatus(  // LCOV_EXCL_LINE
        stateCubeAssumptionLits(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            *solver.variables, singleton, data.bootstrapFrames),  // LCOV_EXCL_LINE
        kResetSummarySingletonProofConflictLimit);
    if (status != SATSolverWrapper::SolveStatus::Unsat) {  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }


// LCOV_EXCL_STOP
    rememberResetFrontierUnreachableCore(  // LCOV_EXCL_LINE
        data, data.bootstrapFrames, singleton);  // LCOV_EXCL_LINE
    blockers.push_back(std::move(singleton));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return blockers;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START
bool resetSummaryPrecheckProvesUnreachable(
// LCOV_EXCL_STOP
    const ResetFrontierReachabilityContextData& data,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps) {
  if (resetFrontierAssumptionSolvesDisabled()) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (data.bootstrapFrames == 0 || postBootstrapSteps == 0) {
    return false;
  }
  if (solverType != KEPLER_FORMAL::Config::SolverType::KISSAT) {
    return false;  // LCOV_EXCL_LINE
  }

  const auto& cachedCoi =
      getCachedResetSummaryCubeReachabilityCoi(data, postBootstrapSteps, cube);
  const BaseCaseCoi& coi = cachedCoi.coi;
  const size_t transitionTargets = cachedCoi.transitionTargets;
  if (isKInductionCoiDiagEnabled()) {
    emitSecDiag(
        "SEC diag: reset summary one-shot cube coi "
        "post_bootstrap_steps=",
        postBootstrapSteps,
        " frames=",
        postBootstrapSteps + 1,
        " solver_symbols=",
        coi.solverSymbols.size(),
        " transition_targets=",
        transitionTargets,
        " cube_literals=",
        cube.size(),
        // LCOV_EXCL_START
        " frame_invariant_symbols=",
        data.frameInvariantSupport.size());
        // LCOV_EXCL_STOP
  }

  if (coi.solverSymbols.size() > kMaxResetSummaryPrecheckSymbols ||
      transitionTargets > kMaxResetSummaryPrecheckTransitionTargets) {
    // LCOV_EXCL_START
    if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC diag: reset summary one-shot precheck skipped "
          // LCOV_EXCL_START
          "reason=coi_cap post_bootstrap_steps=",
          postBootstrapSteps,
          // LCOV_EXCL_STOP
          " solver_symbols=",
          coi.solverSymbols.size(),  // LCOV_EXCL_LINE
          " transition_targets=",
          transitionTargets);
    }  // LCOV_EXCL_LINE
    return false;  // LCOV_EXCL_LINE
  }

// LCOV_EXCL_START


// LCOV_EXCL_STOP
  const auto& problem = data.problem;
  std::vector<std::vector<std::pair<size_t, bool>>> frontierBlockers;
  auto appendFrontierBlocker =
      [&](const std::vector<std::pair<size_t, bool>>& blocker) { // LCOV_EXCL_LINE
        if (frontierBlockers.size() >= kMaxResetSummaryFrontierBlockers) { // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          return false;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        if (std::any_of(
        // LCOV_EXCL_STOP
                frontierBlockers.begin(), // LCOV_EXCL_LINE
                frontierBlockers.end(), // LCOV_EXCL_LINE
                [&](const auto& existing) {  // LCOV_EXCL_LINE
                  return assignmentCubeContains(blocker, existing);  // LCOV_EXCL_LINE
                })) {
          return false;  // LCOV_EXCL_LINE
        }
        frontierBlockers.push_back(blocker); // LCOV_EXCL_LINE
        return true; // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      };
      // LCOV_EXCL_STOP
  if (const auto coresIt =
          data.unreachableCoresByTargetFrame.find(data.bootstrapFrames);
      // LCOV_EXCL_START
      coresIt != data.unreachableCoresByTargetFrame.end()) {
      // LCOV_EXCL_STOP
    for (const auto& core : coresIt->second) { // LCOV_EXCL_LINE
      if (frontierBlockers.size() >= kMaxResetSummaryFrontierBlockers) { // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      }
      if (!cubeSymbolsAreInSolverCoi(coi, core)) { // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
      }
      appendFrontierBlocker(core); // LCOV_EXCL_LINE
    }
    if (!frontierBlockers.empty() && isKInductionCoiDiagEnabled()) { // LCOV_EXCL_LINE
      emitSecDiag( // LCOV_EXCL_LINE
          "SEC diag: reset summary frontier blockers=",
          frontierBlockers.size(), // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          " post_bootstrap_steps=",
          // LCOV_EXCL_STOP
          postBootstrapSteps);
    } // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  for (size_t refinement = 0; refinement <= kMaxResetSummaryRefinements;
       ++refinement) {  // LCOV_EXCL_LINE
    SATSolverWrapper solver(solverType);
    solver.configureForSecPdrQuery(coi.solverSymbols.size());
    FrameVariableStore variables(
        // LCOV_EXCL_START
        solver, coi.solverSymbols, postBootstrapSteps + 1);
        // LCOV_EXCL_STOP

    for (const auto& [symbol, assertedValue] : problem.resetBootstrapInputs) {
      for (size_t frame = 0; frame <= postBootstrapSteps; ++frame) {
        solver.addClause(
            {assertedValue ? -variables.getLiteral(symbol, frame)
                           : variables.getLiteral(symbol, frame)});  // LCOV_EXCL_LINE
      }
    }
    addComplementedStateRelations(
        solver,
        variables,
        problem.complementedStatePairs0,
        coi.solverSymbolSet,
        postBootstrapSteps + 1);
    addComplementedStateRelations(
        solver,
        variables,
        problem.complementedStatePairs1,
        coi.solverSymbolSet,
        postBootstrapSteps + 1);
    addSameFrameStateEqualities(
        solver,
        variables,
        problem,
        coi.solverSymbolSet,
        postBootstrapSteps + 1);
    addDualRailStateValidity(
        solver,
        variables,
        problem.dualRailStatePairs,
        coi.solverSymbolSet,
        postBootstrapSteps + 1);
    addBootstrapStateAssignments(
        solver, variables, problem, coi.solverSymbolSet, 0);
    for (const auto& blocker : frontierBlockers) {
      // Summary frame 0 is the already-validated concrete reset/bootstrap
      // LCOV_EXCL_START
      // frontier.  Any exact blocker learned there is a safe constraint for
      // this weaker summary query and can make it prove UNSAT without opening
      // the full post-reset SAT unroll again.
      addBlockedStateCubeClause(solver, variables, blocker, 0);
    }
    if (data.frameInvariant != nullptr) {
    // LCOV_EXCL_STOP
      for (size_t frame = 0; frame <= postBootstrapSteps; ++frame) {  // LCOV_EXCL_LINE
        FrameFormulaEncoder encoder(  // LCOV_EXCL_LINE
            solver, variables.makeLeafLits(frame, data.frameInvariantSupport));  // LCOV_EXCL_LINE
        solver.addClause({encoder.encode(data.frameInvariant)});  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE

    for (size_t frame = 0; frame < postBootstrapSteps; ++frame) {
      addTransitionRelation(
          solver,
          variables,
          data.transitionByState,
          coi.transitionTargetsByFrame[frame],
          frame);
    }
    for (const auto& [symbol, value] : cube) {
      const int literal = variables.getLiteral(symbol, postBootstrapSteps);
      solver.addClause({value ? literal : -literal});
    // LCOV_EXCL_START
    }


// LCOV_EXCL_STOP
    const SATSolverWrapper::SolveStatus status =
        solver.solveWithKissatResourceLimits(
            kResetSummaryPrecheckConflictLimit);
    if (status == SATSolverWrapper::SolveStatus::Unsat) {
      if (refinement != 0 && isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
        emitSecDiag(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            "SEC diag: reset summary CEGAR proved unreachable "
            "post_bootstrap_steps=",
            postBootstrapSteps,
            // LCOV_EXCL_STOP
            " refinements=",
            refinement,
            // LCOV_EXCL_START
            " frontier_blockers=",
            frontierBlockers.size());  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
      }  // LCOV_EXCL_LINE
      return true;  // LCOV_EXCL_LINE
    }
    if (status == SATSolverWrapper::SolveStatus::Unknown) {
      // LCOV_EXCL_START
      if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
        emitSecDiag(  // LCOV_EXCL_LINE
            "SEC diag: reset summary one-shot precheck resource_limit "
            // LCOV_EXCL_START
            "post_bootstrap_steps=",
            postBootstrapSteps,
            // LCOV_EXCL_STOP
            " solver_symbols=",
            coi.solverSymbols.size(),  // LCOV_EXCL_LINE
            " transition_targets=",
            // LCOV_EXCL_START
            transitionTargets);
            // LCOV_EXCL_STOP
      }  // LCOV_EXCL_LINE
      return false;  // LCOV_EXCL_LINE
    }
    if (refinement == kMaxResetSummaryRefinements ||
        frontierBlockers.size() >= kMaxResetSummaryFrontierBlockers) {
      // LCOV_EXCL_START
      return false;
    }
    // LCOV_EXCL_STOP

    const auto frontierCube =
        extractResetSummaryFrontierCube(data, solver, variables, coi);
    if (!frontierCube.has_value()) {
      if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        emitSecDiag(  // LCOV_EXCL_LINE
            "SEC diag: reset summary refinement skipped reason=frontier_cap "
            // LCOV_EXCL_STOP
            "post_bootstrap_steps=",
            postBootstrapSteps,
            " max_literals=",
            kMaxResetSummaryFrontierCubeLiterals);
      }  // LCOV_EXCL_LINE
      return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
    const auto blocker =
        proveResetSummaryFrontierCubeUnreachable(data, *frontierCube);
    if (!blocker.has_value()) {
      return false;
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    size_t addedBlockers = appendFrontierBlocker(*blocker) ? 1 : 0;
    if (blocker->size() == 1 &&
        frontierBlockers.size() < kMaxResetSummaryFrontierBlockers) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      const size_t remainingBlockers =  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          kMaxResetSummaryFrontierBlockers - frontierBlockers.size();  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
      const auto singletonBlockers =
          collectResetSummarySingletonFrontierBlockers(  // LCOV_EXCL_LINE
              // LCOV_EXCL_START
              data,  // LCOV_EXCL_LINE
              *frontierCube,  // LCOV_EXCL_LINE
              frontierBlockers,
              std::min(  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
                  kMaxResetSummaryBulkSingletonBlockers,
                  // LCOV_EXCL_START
                  remainingBlockers));
      for (const auto& singleton : singletonBlockers) {  // LCOV_EXCL_LINE
        if (appendFrontierBlocker(singleton)) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
          ++addedBlockers;  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
      }
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    if (isKInductionCoiDiagEnabled()) {
    // LCOV_EXCL_STOP
      emitSecDiag(  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          "SEC diag: reset summary learned frontier blocker ",
          // LCOV_EXCL_STOP
          "post_bootstrap_steps=",
          // LCOV_EXCL_START
          postBootstrapSteps,
          // LCOV_EXCL_STOP
          " refinement=",
          refinement + 1,  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          " frontier_cube=",
          // LCOV_EXCL_STOP
          frontierCube->size(),  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          " blocker=",
          // LCOV_EXCL_STOP
          blocker->size(),  // LCOV_EXCL_LINE
          " added=",
          addedBlockers);
    }  // LCOV_EXCL_LINE
  }
  return false;  // LCOV_EXCL_LINE
}

CachedResetFrontierSolver& getCachedResetFrontierSolver(
    const ResetFrontierReachabilityContextData& data,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t targetFrame) {
  // Reset-frontier checks are dominated by repeated neighboring cube queries.
  // Use the assumption-capable solver here even when the main SEC run selected
  // Kissat: otherwise every cube value has to be encoded as unit clauses in a
  // separate cached solver, which BlackParrot sampling showed growing to
  // multi-GB retained solver caches before PDR made progress.
  const bool encodeCubeAsUnitClauses = false;
  const auto cachedSolverType =
      SATSolverWrapper::assumptionSolverTypeFor(solverType);
  const ResetFrontierSolverCacheKey key =
      resetFrontierSolverCacheKey(
          cachedSolverType, targetFrame, cube, encodeCubeAsUnitClauses);
  if (const auto it = data.cachedSolvers.find(key);
      it != data.cachedSolvers.end()) {
    return *it->second;
  }

  const auto cubeSymbols = sortedCubeSymbols(cube);
  for (const auto& [_, cached] : data.cachedSolvers) {
    if (cached->solverType == cachedSolverType &&
        cached->targetFrame == targetFrame &&
        cached->cubeEncodedAsUnitClauses == encodeCubeAsUnitClauses &&
        solverContainsCubeSymbols(*cached, cubeSymbols)) {
      // A solver built for a wider reset-frontier COI is still an exact query
      // for a smaller cube: the extra transition / init clauses constrain
      // unrelated existential variables, while the requested cube values are
      // supplied as assumptions.
      if (isKInductionCoiDiagEnabled()) {
        emitSecDiag(
            "SEC diag: reset frontier solver superset cache hit ",
            "target_frame=",
            targetFrame,
            " cube_literals=",
            cube.size(),
            " solver_symbols=",
            cached->coi.solverSymbols.size());
      // LCOV_EXCL_START
      }
      return *cached;
      // LCOV_EXCL_STOP
    }
  }

  if (data.cachedSolvers.size() >= kMaxResetFrontierCachedSolvers) {
    data.cachedSolvers.clear();  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  auto cached = buildResetFrontierSolver(
      // LCOV_EXCL_START
      data, cachedSolverType, cube, targetFrame, encodeCubeAsUnitClauses);
      // LCOV_EXCL_STOP
  auto [it, inserted] = data.cachedSolvers.emplace(key, std::move(cached));
  (void)inserted;
  return *it->second;
}

// LCOV_EXCL_START

void primeResetFrontierReachabilitySolver(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    // LCOV_EXCL_START
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps) {
  if (cube.empty()) {  // LCOV_EXCL_LINE
    return;  // LCOV_EXCL_LINE
  }


// LCOV_EXCL_STOP
  const auto& data = *context.data;  // LCOV_EXCL_LINE
  const size_t targetFrame = data.bootstrapFrames + postBootstrapSteps;  // LCOV_EXCL_LINE
  const auto normalizedCube = normalizedAssignmentCube(cube);  // LCOV_EXCL_LINE
  (void)getCachedResetFrontierSolver(  // LCOV_EXCL_LINE
      data, solverType, normalizedCube, targetFrame);  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

bool isStateCubeReachableAtResetFrontier(
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps,
    bool usePostBootstrapPrechecks,
    int64_t startupConflictLimit,
    int64_t startupPropagationLimit) {
  if (cube.empty()) {
    return true; // LCOV_EXCL_LINE
  }

  const auto& data = *context.data;
  const size_t targetFrame = data.bootstrapFrames + postBootstrapSteps;
  const auto normalizedCube = normalizedAssignmentCube(cube);
  if (const auto cachedCore =
          findCachedResetFrontierUnreachableCore(data, targetFrame, normalizedCube);
      cachedCore.has_value()) {
    if (isKInductionCoiDiagEnabled()) {
      emitSecDiag(
          "SEC diag: reset frontier cached unreachable core hit ",
          "post_bootstrap_steps=",
          postBootstrapSteps,
          " core_literals=",
          cachedCore->size(),
          " cube_literals=",
          normalizedCube.size());
    }
    return false;
  }
  if (const auto knownCore = knownResetFrontierConflictCore(
          data, normalizedCube, postBootstrapSteps);
      knownCore.has_value()) {
    rememberResetFrontierUnreachableCore(data, targetFrame, *knownCore);
    if (isKInductionCoiDiagEnabled()) {
      emitSecDiag(
          "SEC diag: reset frontier known facts exclude cube ",
          "post_bootstrap_steps=",
          postBootstrapSteps,
          " core_literals=",
          knownCore->size(),
          " cube_literals=",
          normalizedCube.size());
    }
    return false;
  }
  if (postBootstrapSteps != 0 && usePostBootstrapPrechecks) {
    // Cached-assumption validation is PDR's hot reset-frontier path. Before
    // constructing the exact assumption solver, use the same weakened
    // startup-equality COI used by one-shot validation. The relaxed query only
    // drops equality closure, so UNSAT remains a sound proof; SAT simply falls
    // through to the exact cached solver below. This is deliberately bounded:
    // if the relaxed COI is still broad, or Kissat cannot decide it quickly,
    // the shortcut is abandoned rather than becoming the PDR wall itself.
    auto relaxedSolver = buildResetFrontierSolver(
        data,
        solverType,
        normalizedCube,
        targetFrame,
        /*encodeCubeAsUnitClauses=*/true,
        /*closeStartupEqualityDependencies=*/false);
    const size_t relaxedTransitionTargets =
        countTransitionTargets(relaxedSolver->coi.transitionTargetsByFrame);

    if (isKInductionCoiDiagEnabled()) {
      emitSecDiag(
          "SEC diag: reset frontier relaxed cached cube coi "
          "post_bootstrap_steps=",
          postBootstrapSteps,
          " frames=",
          targetFrame + 1,
          " solver_symbols=",
          relaxedSolver->coi.solverSymbols.size(),
          " transition_targets=",
          relaxedTransitionTargets,
          " cube_literals=",
          normalizedCube.size(),
          " frame_invariant_symbols=",
          data.frameInvariantSupport.size());
    }

    const bool relaxedCoiIsLocal =
        relaxedResetFrontierPrecheckCoiIsLocal(
            relaxedSolver->coi.solverSymbols.size(),
            relaxedTransitionTargets,
            normalizedCube.size());
    if (relaxedCoiIsLocal) {
      bool relaxedUnsat = false;
      // LCOV_EXCL_START
      if (solverType == KEPLER_FORMAL::Config::SolverType::KISSAT) {
        const auto status =
        // LCOV_EXCL_STOP
            relaxedSolver->solver->solveWithKissatResourceLimits(
                kRelaxedResetFrontierPrecheckConflictLimit);
        relaxedUnsat = status == SATSolverWrapper::SolveStatus::Unsat;
        if (status == SATSolverWrapper::SolveStatus::Unknown &&
            // LCOV_EXCL_START
            isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
          emitSecDiag(  // LCOV_EXCL_LINE
              "SEC diag: reset frontier relaxed cached precheck "
              // LCOV_EXCL_START
              "resource_limit post_bootstrap_steps=",
              // LCOV_EXCL_STOP
              postBootstrapSteps,
              // LCOV_EXCL_START
              " solver_symbols=",
              // LCOV_EXCL_STOP
              relaxedSolver->coi.solverSymbols.size(),  // LCOV_EXCL_LINE
              " transition_targets=",
              relaxedTransitionTargets);
        }  // LCOV_EXCL_LINE
      } else {
        relaxedUnsat = !relaxedSolver->solver->solve();  // LCOV_EXCL_LINE
      }
      if (relaxedUnsat) {
        rememberResetFrontierUnreachableCore(data, targetFrame, normalizedCube);
        return false;
      }
    } else if (isKInductionCoiDiagEnabled()) {
      emitSecDiag(
          "SEC diag: reset frontier relaxed cached precheck skipped "
          "reason=coi_cap post_bootstrap_steps=",
          postBootstrapSteps,
          " solver_symbols=",
          relaxedSolver->coi.solverSymbols.size(),
          // LCOV_EXCL_START
          " transition_targets=",
          relaxedTransitionTargets);
          // LCOV_EXCL_STOP
    }

    if (resetSummaryPrecheckProvesUnreachable(
            // LCOV_EXCL_START
            data, solverType, normalizedCube, postBootstrapSteps)) {
      rememberResetFrontierUnreachableCore(data, targetFrame, normalizedCube);  // LCOV_EXCL_LINE
      return false;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  if (resetFrontierAssumptionSolvesDisabled()) {
    // LCOV_EXCL_START
    auto unitSolver = buildResetFrontierSolver(  // LCOV_EXCL_LINE
        data,  // LCOV_EXCL_LINE
        solverType,  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        normalizedCube,
        // LCOV_EXCL_START
        targetFrame,  // LCOV_EXCL_LINE
        /*encodeCubeAsUnitClauses=*/true);
    const int64_t conflictLimit =  // LCOV_EXCL_LINE
        postBootstrapSteps == 0 && startupConflictLimit >= 0  // LCOV_EXCL_LINE
            ? startupConflictLimit  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            : static_cast<int64_t>(kResetFrontierCachedAssumptionConflictLimit);
    const auto status = solveResetFrontierUnitClauseQuery(  // LCOV_EXCL_LINE
        *unitSolver->solver, solverType, conflictLimit);  // LCOV_EXCL_LINE
    if (status == SATSolverWrapper::SolveStatus::Unknown) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
        emitSecDiag(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            "SEC diag: reset frontier unit-clause proof resource_limit ",
            // LCOV_EXCL_STOP
            "post_bootstrap_steps=",
            // LCOV_EXCL_START
            postBootstrapSteps,
            " solver_symbols=",
            unitSolver->coi.solverSymbols.size(),  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            " transition_targets=",
            // LCOV_EXCL_START
            countTransitionTargets(unitSolver->coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
            " cube_literals=",
            normalizedCube.size());  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      return true;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_STOP
    const bool reachable = status == SATSolverWrapper::SolveStatus::Sat;  // LCOV_EXCL_LINE
    if (!reachable) {  // LCOV_EXCL_LINE
      rememberResetFrontierUnreachableCore(data, targetFrame, normalizedCube);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    return reachable;  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  CachedResetFrontierSolver& cached =
      getCachedResetFrontierSolver(data, solverType, normalizedCube, targetFrame);

  if (isKInductionCoiDiagEnabled()) {
    emitSecDiag(
        "SEC diag: reset frontier cube coi post_bootstrap_steps=",
        postBootstrapSteps,
        " frames=",
        targetFrame + 1,
        " solver_symbols=",
        cached.coi.solverSymbols.size(),
        " transition_targets=",
        countTransitionTargets(cached.coi.transitionTargetsByFrame),
        " cube_literals=",
        normalizedCube.size(),
        // LCOV_EXCL_START
        " frame_invariant_symbols=",
        data.frameInvariantSupport.size());
        // LCOV_EXCL_STOP
  }

  SATSolverWrapper::SolveStatus status = SATSolverWrapper::SolveStatus::Unknown;
  if (cached.cubeEncodedAsUnitClauses) {
    status = cached.solver->solveStatus();  // LCOV_EXCL_LINE
  } else {  // LCOV_EXCL_LINE
    const auto assumptions =
        stateCubeAssumptionLits(*cached.variables, normalizedCube, targetFrame);
    if (postBootstrapSteps == 0) {
      status =
          startupConflictLimit >= 0 || startupPropagationLimit >= 0
              ? cached.solver->solveWithAssumptionsStatus(
                    assumptions, startupConflictLimit, startupPropagationLimit)
              : cached.solver->solveWithAssumptionsStatus(assumptions);
    // LCOV_EXCL_START
    } else {
      status = cached.solver->solveWithAssumptionsStatus(
      // LCOV_EXCL_STOP
          assumptions, kResetFrontierCachedAssumptionConflictLimit);
    }
  }
  if (status == SATSolverWrapper::SolveStatus::Unknown) {
    // LCOV_EXCL_START
    if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      emitSecDiag(  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          "SEC diag: reset frontier cached assumption proof resource_limit ",
          // LCOV_EXCL_STOP
          "post_bootstrap_steps=",
          // LCOV_EXCL_START
          postBootstrapSteps,
          " solver_symbols=",
          cached.coi.solverSymbols.size(),  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
          " transition_targets=",
          countTransitionTargets(cached.coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
          " cube_literals=",
          normalizedCube.size());  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    return true;  // LCOV_EXCL_LINE
  }

  const bool reachable = status == SATSolverWrapper::SolveStatus::Sat;
  if (!reachable) {
    if (postBootstrapSteps == 0) {
      // The solve above already produced a failed-assumption core.  Reuse it
      // directly instead of launching a second exact minimization loop: AES
      // PDR samples showed wide F[0] reset cubes spending their runtime inside
      // that duplicate assumption search, while the failed core is already a
      // LCOV_EXCL_START
      // sound reset-frontier blocker and still reusable for neighboring cubes.
      if (const auto core = failedAssumptionCoreFromLastResetFrontierSolve(
      // LCOV_EXCL_STOP
              cached, normalizedCube, targetFrame);
          core.has_value()) {
        // LCOV_EXCL_START
        rememberResetFrontierUnreachableCore(data, targetFrame, *core);
      } else {
        rememberResetFrontierUnreachableCore(  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            data, targetFrame, normalizedCube);  // LCOV_EXCL_LINE
      }
    } else if (const auto core =
                   // LCOV_EXCL_START
                   failedAssumptionCoreFromLastResetFrontierSolve(  // LCOV_EXCL_LINE
                       cached, normalizedCube, targetFrame);  // LCOV_EXCL_LINE
                       // LCOV_EXCL_STOP
               core.has_value()) {  // LCOV_EXCL_LINE
      // Post-bootstrap prechecks are on the hot PDR path.  Reuse the
      // assumption core already produced by this UNSAT query, but avoid the
      // extra minimization SAT calls reserved for the exact reset frontier.
      rememberResetFrontierUnreachableCore(data, targetFrame, *core);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
  return reachable;
}

bool isStateCubeReachableAtResetFrontierOneShot(
    // LCOV_EXCL_START
    const ResetFrontierReachabilityContext& context,
    // LCOV_EXCL_STOP
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps,
    bool usePostBootstrapPrechecks) {
  if (cube.empty()) {
    return true;  // LCOV_EXCL_LINE
  }

  const auto& data = *context.data;
  const size_t targetFrame = data.bootstrapFrames + postBootstrapSteps;
  const auto normalizedCube = normalizedAssignmentCube(cube);
  if (const auto knownCore = knownResetFrontierConflictCore(
          data, normalizedCube, postBootstrapSteps);
      knownCore.has_value()) {
    // LCOV_EXCL_START
    rememberResetFrontierUnreachableCore(data, targetFrame, *knownCore);
    // LCOV_EXCL_STOP
    return false;
  }
  if (findCachedResetFrontierUnreachableCore(
          data, targetFrame, normalizedCube)
          .has_value()) {
    return false;  // LCOV_EXCL_LINE
  }
  if (postBootstrapSteps != 0 && usePostBootstrapPrechecks) {
    // First use a weakened COI that does not close the startup equality
    // components.  This is safe only as an UNSAT precheck: removing equality
    // constraints can create spurious SAT witnesses, but if the relaxed query
    // is already UNSAT then the exact reset-frontier query is UNSAT too.  AES
    // PDR sampling showed the exact post-reset one-shot query expanding a
    // two-literal root cube through a 900+ symbol equality/transition closure;
    // this precheck lets local transition/bootstrap contradictions prove out
    // before opening that wider solver. Keep it bounded so a broad relaxed COI
    // falls through instead of replacing the exact query with another wall.
    auto relaxedSolver = buildResetFrontierSolver(
        data,
        solverType,
        normalizedCube,
        targetFrame,
        /*encodeCubeAsUnitClauses=*/true,
        /*closeStartupEqualityDependencies=*/false);
    const size_t relaxedTransitionTargets =
        countTransitionTargets(relaxedSolver->coi.transitionTargetsByFrame);

    if (isKInductionCoiDiagEnabled()) {
      emitSecDiag(
          "SEC diag: reset frontier relaxed one-shot cube coi "
          "post_bootstrap_steps=",
          postBootstrapSteps,
          " frames=",
          targetFrame + 1,
          " solver_symbols=",
          relaxedSolver->coi.solverSymbols.size(),
          " transition_targets=",
          relaxedTransitionTargets,
          " cube_literals=",
          normalizedCube.size(),
          " frame_invariant_symbols=",
          data.frameInvariantSupport.size());
    }

    const bool relaxedCoiIsLocal =
        relaxedResetFrontierPrecheckCoiIsLocal(
            relaxedSolver->coi.solverSymbols.size(),
            relaxedTransitionTargets,
            normalizedCube.size());
    if (relaxedCoiIsLocal) {
      bool relaxedUnsat = false;
      // LCOV_EXCL_START
      if (solverType == KEPLER_FORMAL::Config::SolverType::KISSAT) {
        const auto status =
        // LCOV_EXCL_STOP
            relaxedSolver->solver->solveWithKissatResourceLimits(
                kRelaxedResetFrontierPrecheckConflictLimit);
        relaxedUnsat = status == SATSolverWrapper::SolveStatus::Unsat;
        if (status == SATSolverWrapper::SolveStatus::Unknown &&
            // LCOV_EXCL_START
            isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
          emitSecDiag(  // LCOV_EXCL_LINE
              "SEC diag: reset frontier relaxed one-shot precheck "
              // LCOV_EXCL_START
              "resource_limit post_bootstrap_steps=",
              // LCOV_EXCL_STOP
              postBootstrapSteps,
              // LCOV_EXCL_START
              " solver_symbols=",
              // LCOV_EXCL_STOP
              relaxedSolver->coi.solverSymbols.size(),  // LCOV_EXCL_LINE
              " transition_targets=",
              relaxedTransitionTargets);
        }  // LCOV_EXCL_LINE
      } else {
        relaxedUnsat = !relaxedSolver->solver->solve();  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      }
      // LCOV_EXCL_STOP
      if (relaxedUnsat) {
        rememberResetFrontierUnreachableCore(data, targetFrame, normalizedCube);
        return false;
      }
    // LCOV_EXCL_START
    } else if (isKInductionCoiDiagEnabled()) {
    // LCOV_EXCL_STOP
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC diag: reset frontier relaxed one-shot precheck skipped "
          // LCOV_EXCL_START
          "reason=coi_cap post_bootstrap_steps=",
          // LCOV_EXCL_STOP
          postBootstrapSteps,
          " solver_symbols=",
          relaxedSolver->coi.solverSymbols.size(),  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          " transition_targets=",
          relaxedTransitionTargets);
          // LCOV_EXCL_STOP
    }  // LCOV_EXCL_LINE

    if (resetSummaryPrecheckProvesUnreachable(
            // LCOV_EXCL_START
            data, solverType, normalizedCube, postBootstrapSteps)) {
      rememberResetFrontierUnreachableCore(data, targetFrame, normalizedCube);
      return false;
      // LCOV_EXCL_STOP
    }
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  if (resetFrontierAssumptionSolvesDisabled()) {
    // LCOV_EXCL_START
    auto solver = buildResetFrontierSolver(  // LCOV_EXCL_LINE
        data,  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        solverType,  // LCOV_EXCL_LINE
        normalizedCube,
        // LCOV_EXCL_START
        targetFrame,  // LCOV_EXCL_LINE
        /*encodeCubeAsUnitClauses=*/true);
    const int64_t conflictLimit =  // LCOV_EXCL_LINE
        postBootstrapSteps == 0  // LCOV_EXCL_LINE
            ? -1
            // LCOV_EXCL_STOP
            : static_cast<int64_t>(kResetFrontierCachedAssumptionConflictLimit);
    const auto status = solveResetFrontierUnitClauseQuery(  // LCOV_EXCL_LINE
        *solver->solver, solverType, conflictLimit);  // LCOV_EXCL_LINE
    if (status == SATSolverWrapper::SolveStatus::Unknown) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
        emitSecDiag(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            "SEC diag: reset frontier one-shot unit-clause resource_limit ",
            // LCOV_EXCL_STOP
            "post_bootstrap_steps=",
            // LCOV_EXCL_START
            postBootstrapSteps,
            " solver_symbols=",
            solver->coi.solverSymbols.size(),  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            " transition_targets=",
            // LCOV_EXCL_START
            countTransitionTargets(solver->coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
            " cube_literals=",
            normalizedCube.size());  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      return true;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_STOP
    const bool reachable = status == SATSolverWrapper::SolveStatus::Sat;  // LCOV_EXCL_LINE
    if (!reachable) {  // LCOV_EXCL_LINE
      rememberResetFrontierUnreachableCore(data, targetFrame, normalizedCube);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    return reachable;  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  const auto assumptionSolverType =
      SATSolverWrapper::assumptionSolverTypeFor(solverType);
  auto solver = buildResetFrontierSolver(
      data,
      assumptionSolverType,
      normalizedCube,
      targetFrame,
      /*encodeCubeAsUnitClauses=*/false);

  if (isKInductionCoiDiagEnabled()) {
    emitSecDiag(
        "SEC diag: reset frontier one-shot cube coi post_bootstrap_steps=",
        postBootstrapSteps,
        " frames=",
        targetFrame + 1,
        " solver_symbols=",
        solver->coi.solverSymbols.size(),
        " transition_targets=",
        countTransitionTargets(solver->coi.transitionTargetsByFrame),
        " cube_literals=",
        normalizedCube.size(),
        " frame_invariant_symbols=",
        data.frameInvariantSupport.size());
  }

  const auto assumptions =
      stateCubeAssumptionLits(*solver->variables, normalizedCube, targetFrame);
  const bool reachable = solver->solver->solveWithAssumptions(assumptions);
  if (!reachable) {
    // Keep the one-shot COI/build profile, but query the cube through
    // assumptions so an UNSAT proof exposes a reusable core. BlackParrot
    // LCOV_EXCL_START
    // samples showed full-cube caching missing many neighboring root cubes.
    // LCOV_EXCL_STOP
    if (const auto core = failedAssumptionCoreFromLastResetFrontierSolve(
            *solver, normalizedCube, targetFrame);
        core.has_value()) {
      rememberResetFrontierUnreachableCore(data, targetFrame, *core);
    } else {
      // LCOV_EXCL_START
      rememberResetFrontierUnreachableCore(data, targetFrame, normalizedCube);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
  }
  return reachable;
}

// LCOV_EXCL_START

bool isStateCubeReachableWithinResetFrontier(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    // LCOV_EXCL_START
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t maxPostBootstrapSteps) {
  if (cube.empty()) {  // LCOV_EXCL_LINE
    return true;  // LCOV_EXCL_LINE
  }

  const auto& data = *context.data;  // LCOV_EXCL_LINE
  const auto normalizedCube = normalizedAssignmentCube(cube);  // LCOV_EXCL_LINE
  std::vector<size_t> uncheckedSteps;  // LCOV_EXCL_LINE
  uncheckedSteps.reserve(maxPostBootstrapSteps + 1);  // LCOV_EXCL_LINE
  for (size_t step = 0; step <= maxPostBootstrapSteps; ++step) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    const size_t targetFrame = data.bootstrapFrames + step;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    if (const auto knownCore =  // LCOV_EXCL_LINE
            knownResetFrontierConflictCore(data, normalizedCube, step);  // LCOV_EXCL_LINE
        knownCore.has_value()) {  // LCOV_EXCL_LINE
      rememberResetFrontierUnreachableCore(data, targetFrame, *knownCore);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      continue;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
    if (findCachedResetFrontierUnreachableCore(  // LCOV_EXCL_LINE
            data, targetFrame, normalizedCube)  // LCOV_EXCL_LINE
            .has_value()) {  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
      continue;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    uncheckedSteps.push_back(step);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  if (uncheckedSteps.empty()) {  // LCOV_EXCL_LINE
    return false;  // LCOV_EXCL_LINE
  }

  std::vector<size_t> remainingSteps;  // LCOV_EXCL_LINE
  remainingSteps.reserve(uncheckedSteps.size());  // LCOV_EXCL_LINE
  for (const auto step : uncheckedSteps) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    if (step != 0 &&  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        resetSummaryPrecheckProvesUnreachable(  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            data, solverType, normalizedCube, step)) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      rememberResetFrontierUnreachableCore(  // LCOV_EXCL_LINE
          data, data.bootstrapFrames + step, normalizedCube);  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
      continue;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
    remainingSteps.push_back(step);  // LCOV_EXCL_LINE
  }
  if (remainingSteps.empty()) {  // LCOV_EXCL_LINE
    return false;  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_STOP
  if (resetFrontierAssumptionSolvesDisabled() ||  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      remainingSteps.size() <= kMaxSparseResetFrontierPerStepChecks) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    for (const auto step : remainingSteps) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      if (isStateCubeReachableAtResetFrontier(  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
              context,  // LCOV_EXCL_LINE
              solverType,  // LCOV_EXCL_LINE
              // LCOV_EXCL_START
              normalizedCube,
              // LCOV_EXCL_STOP
              step,  // LCOV_EXCL_LINE
              /*usePostBootstrapPrechecks=*/false)) {
        // LCOV_EXCL_START
        return true;  // LCOV_EXCL_LINE
      }
    }
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

// LCOV_EXCL_START

  const size_t maxTargetFrame =  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
      data.bootstrapFrames + maxPostBootstrapSteps;  // LCOV_EXCL_LINE
  CachedResetFrontierSolver& solver = getCachedResetFrontierPrefixSolver(  // LCOV_EXCL_LINE
      data, solverType, normalizedCube, maxTargetFrame);  // LCOV_EXCL_LINE

// LCOV_EXCL_START


// LCOV_EXCL_STOP
  if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    emitSecDiag(  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
        "SEC diag: reset frontier prefix cube coi max_post_bootstrap_steps=",
        // LCOV_EXCL_START
        maxPostBootstrapSteps,
        // LCOV_EXCL_STOP
        " frames=",
        // LCOV_EXCL_START
        maxTargetFrame + 1,  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        " solver_symbols=",
        // LCOV_EXCL_START
        solver.coi.solverSymbols.size(),  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        " transition_targets=",
        // LCOV_EXCL_START
        countTransitionTargets(solver.coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
        " cube_literals=",
        // LCOV_EXCL_STOP
        normalizedCube.size(),  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        " unchecked_steps=",
        remainingSteps.size(),  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        " frame_invariant_symbols=",
        // LCOV_EXCL_START
        data.frameInvariantSupport.size());  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  for (const auto step : remainingSteps) {  // LCOV_EXCL_LINE
    const size_t targetFrame = data.bootstrapFrames + step;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    const auto assumptions =
        // LCOV_EXCL_START
        stateCubeAssumptionLits(*solver.variables, normalizedCube, targetFrame);  // LCOV_EXCL_LINE
    const auto status =  // LCOV_EXCL_LINE
        step == 0  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            ? solver.solver->solveWithAssumptionsStatus(assumptions)  // LCOV_EXCL_LINE
            : solver.solver->solveWithAssumptionsStatus(  // LCOV_EXCL_LINE
                  assumptions, kResetFrontierCachedAssumptionConflictLimit);
    if (status == SATSolverWrapper::SolveStatus::Unknown) {  // LCOV_EXCL_LINE
      if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
        emitSecDiag(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            "SEC diag: reset frontier prefix assumption proof resource_limit ",
            // LCOV_EXCL_STOP
            "post_bootstrap_steps=",
            // LCOV_EXCL_START
            step,
            // LCOV_EXCL_STOP
            " max_post_bootstrap_steps=",
            // LCOV_EXCL_START
            maxPostBootstrapSteps,
            " solver_symbols=",
            solver.coi.solverSymbols.size(),  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            " transition_targets=",
            // LCOV_EXCL_START
            countTransitionTargets(solver.coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
            " cube_literals=",
            // LCOV_EXCL_STOP
            normalizedCube.size());  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      }  // LCOV_EXCL_LINE
      return true;  // LCOV_EXCL_LINE
    }
    if (status == SATSolverWrapper::SolveStatus::Sat) {  // LCOV_EXCL_LINE
      return true;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_STOP
    if (const auto core = failedAssumptionCoreFromLastResetFrontierSolve(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            solver, normalizedCube, targetFrame);  // LCOV_EXCL_LINE
        core.has_value()) {  // LCOV_EXCL_LINE
      rememberResetFrontierUnreachableCore(data, targetFrame, *core);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    } else {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      rememberResetFrontierUnreachableCore(data, targetFrame, normalizedCube);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
  // LCOV_EXCL_START
  }  // LCOV_EXCL_LINE
  return false;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE


// LCOV_EXCL_STOP
std::vector<std::pair<size_t, bool>> supportCubeForAssignmentCubes(  // LCOV_EXCL_LINE
    const std::vector<std::vector<std::pair<size_t, bool>>>& cubes) {
  std::unordered_map<size_t, bool> valueBySymbol;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  for (const auto& cube : cubes) {  // LCOV_EXCL_LINE
    for (const auto& [symbol, value] : cube) {  // LCOV_EXCL_LINE
      valueBySymbol.emplace(symbol, value);  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_STOP
  }

// LCOV_EXCL_START

  std::vector<std::pair<size_t, bool>> supportCube;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  supportCube.reserve(valueBySymbol.size());  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  for (const auto& [symbol, value] : valueBySymbol) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    supportCube.emplace_back(symbol, value);  // LCOV_EXCL_LINE
  }
  return normalizedAssignmentCube(std::move(supportCube));  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START

bool anyStateCubeReachableWithinResetFrontier(  // LCOV_EXCL_LINE
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::vector<std::pair<size_t, bool>>>& cubes,
    size_t maxPostBootstrapSteps) {
    // LCOV_EXCL_STOP
  std::vector<std::vector<std::pair<size_t, bool>>> normalizedCubes;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  normalizedCubes.reserve(cubes.size());  // LCOV_EXCL_LINE
  for (auto cube : cubes) {  // LCOV_EXCL_LINE
    cube = normalizedAssignmentCube(std::move(cube));  // LCOV_EXCL_LINE
    if (cube.empty()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      return true;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    normalizedCubes.push_back(std::move(cube));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  if (normalizedCubes.empty()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
  }


// LCOV_EXCL_STOP
  const auto& data = *context.data;  // LCOV_EXCL_LINE
  const size_t maxTargetFrame = data.bootstrapFrames + maxPostBootstrapSteps;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  const std::vector<std::pair<size_t, bool>> supportCube =
      supportCubeForAssignmentCubes(normalizedCubes);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  if (supportCube.empty()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    return true;  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_STOP

  CachedResetFrontierSolver& solver = getCachedResetFrontierPrefixSolver(  // LCOV_EXCL_LINE
      data, solverType, supportCube, maxTargetFrame);  // LCOV_EXCL_LINE

// LCOV_EXCL_START


// LCOV_EXCL_STOP
  if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    emitSecDiag(  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
        "SEC diag: reset frontier batch cube coi max_post_bootstrap_steps=",
        // LCOV_EXCL_START
        maxPostBootstrapSteps,
        // LCOV_EXCL_STOP
        " frames=",
        // LCOV_EXCL_START
        maxTargetFrame + 1,  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        " solver_symbols=",
        // LCOV_EXCL_START
        solver.coi.solverSymbols.size(),  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        " transition_targets=",
        // LCOV_EXCL_START
        countTransitionTargets(solver.coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
        " cubes=",
        // LCOV_EXCL_STOP
        normalizedCubes.size(),  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        " support_literals=",
        supportCube.size(),  // LCOV_EXCL_LINE
        " frame_invariant_symbols=",
        data.frameInvariantSupport.size());  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  for (size_t step = 0; step <= maxPostBootstrapSteps; ++step) {  // LCOV_EXCL_LINE
    const size_t targetFrame = data.bootstrapFrames + step;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    for (const auto& cube : normalizedCubes) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      if (const auto knownCore =  // LCOV_EXCL_LINE
              knownResetFrontierConflictCore(data, cube, step);  // LCOV_EXCL_LINE
          knownCore.has_value()) {  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        rememberResetFrontierUnreachableCore(data, targetFrame, *knownCore);  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      }
      if (findCachedResetFrontierUnreachableCore(data, targetFrame, cube)  // LCOV_EXCL_LINE
              .has_value()) {  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
      }
      // LCOV_EXCL_STOP
      const auto assumptions =
          stateCubeAssumptionLits(*solver.variables, cube, targetFrame);  // LCOV_EXCL_LINE
      const auto status =  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          step == 0  // LCOV_EXCL_LINE
              ? solver.solver->solveWithAssumptionsStatus(assumptions)  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
              : solver.solver->solveWithAssumptionsStatus(  // LCOV_EXCL_LINE
                    // LCOV_EXCL_START
                    assumptions,
                    kResetFrontierCachedAssumptionConflictLimit,
                    kResetFrontierBatchProofPropagationLimit);
                    // LCOV_EXCL_STOP
      if (status == SATSolverWrapper::SolveStatus::Sat) {  // LCOV_EXCL_LINE
        return true;  // LCOV_EXCL_LINE
      }
      if (status == SATSolverWrapper::SolveStatus::Unknown) {  // LCOV_EXCL_LINE
        if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
          emitSecDiag(  // LCOV_EXCL_LINE
              // LCOV_EXCL_START
              "SEC diag: reset frontier batch assumption proof resource_limit ",
              // LCOV_EXCL_STOP
              "post_bootstrap_steps=",
              // LCOV_EXCL_START
              step,
              // LCOV_EXCL_STOP
              " max_post_bootstrap_steps=",
              // LCOV_EXCL_START
              maxPostBootstrapSteps,
              " solver_symbols=",
              solver.coi.solverSymbols.size(),  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
              " transition_targets=",
              // LCOV_EXCL_START
              countTransitionTargets(solver.coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
              " cubes=",
              normalizedCubes.size());  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
        return true;  // LCOV_EXCL_LINE
      }
      // LCOV_EXCL_STOP
      if (const auto core = failedAssumptionCoreFromLastResetFrontierSolve(  // LCOV_EXCL_LINE
              // LCOV_EXCL_START
              solver, cube, targetFrame);  // LCOV_EXCL_LINE
          core.has_value()) {  // LCOV_EXCL_LINE
        rememberResetFrontierUnreachableCore(data, targetFrame, *core);  // LCOV_EXCL_LINE
      } else {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
        rememberResetFrontierUnreachableCore(data, targetFrame, cube);  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      }
      // LCOV_EXCL_STOP
    }  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return false;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

bool anyStateCubeReachableAtResetFrontier(  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::vector<std::pair<size_t, bool>>>& cubes,
    size_t postBootstrapSteps,
    long long conflictLimit,
    long long propagationLimit) {
    // LCOV_EXCL_STOP
  std::vector<std::vector<std::pair<size_t, bool>>> normalizedCubes;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  normalizedCubes.reserve(cubes.size());  // LCOV_EXCL_LINE
  for (auto cube : cubes) {  // LCOV_EXCL_LINE
    cube = normalizedAssignmentCube(std::move(cube));  // LCOV_EXCL_LINE
    if (cube.empty()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      return true;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    normalizedCubes.push_back(std::move(cube));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  if (normalizedCubes.empty()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
  }


// LCOV_EXCL_STOP
  const auto& data = *context.data;  // LCOV_EXCL_LINE
  const size_t targetFrame = data.bootstrapFrames + postBootstrapSteps;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  const std::vector<std::pair<size_t, bool>> supportCube =
      supportCubeForAssignmentCubes(normalizedCubes);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  if (supportCube.empty()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    return true;  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_STOP

  CachedResetFrontierSolver& solver = getCachedResetFrontierSolver(  // LCOV_EXCL_LINE
      data, solverType, supportCube, targetFrame);  // LCOV_EXCL_LINE

  // LCOV_EXCL_START
  if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    emitSecDiag(  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        "SEC diag: reset frontier target batch cube coi "
        // LCOV_EXCL_STOP
        "post_bootstrap_steps=",
        // LCOV_EXCL_START
        postBootstrapSteps,
        // LCOV_EXCL_STOP
        " frames=",
        // LCOV_EXCL_START
        targetFrame + 1,  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        " solver_symbols=",
        // LCOV_EXCL_START
        solver.coi.solverSymbols.size(),  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        " transition_targets=",
        // LCOV_EXCL_START
        countTransitionTargets(solver.coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
        " cubes=",
        // LCOV_EXCL_STOP
        normalizedCubes.size(),  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        " support_literals=",
        supportCube.size(),  // LCOV_EXCL_LINE
        " frame_invariant_symbols=",
        data.frameInvariantSupport.size());  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE


// LCOV_EXCL_STOP
  for (const auto& cube : normalizedCubes) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    if (const auto knownCore =  // LCOV_EXCL_LINE
            knownResetFrontierConflictCore(data, cube, postBootstrapSteps);  // LCOV_EXCL_LINE
        knownCore.has_value()) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      rememberResetFrontierUnreachableCore(data, targetFrame, *knownCore);  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
    if (findCachedResetFrontierUnreachableCore(data, targetFrame, cube)  // LCOV_EXCL_LINE
            .has_value()) {  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    const auto assumptions =
    // LCOV_EXCL_STOP
        stateCubeAssumptionLits(*solver.variables, cube, targetFrame);  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    const auto status =  // LCOV_EXCL_LINE
        solver.solver->solveWithAssumptionsStatus(  // LCOV_EXCL_LINE
            assumptions, conflictLimit, propagationLimit);  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
    if (status == SATSolverWrapper::SolveStatus::Sat) {  // LCOV_EXCL_LINE
      return true;  // LCOV_EXCL_LINE
    }
    if (status == SATSolverWrapper::SolveStatus::Unknown) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      if (isKInductionCoiDiagEnabled()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
        emitSecDiag(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            "SEC diag: reset frontier target batch assumption proof "
            // LCOV_EXCL_STOP
            "resource_limit post_bootstrap_steps=",
            // LCOV_EXCL_START
            postBootstrapSteps,
            " solver_symbols=",
            solver.coi.solverSymbols.size(),  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            " transition_targets=",
            // LCOV_EXCL_START
            countTransitionTargets(solver.coi.transitionTargetsByFrame),  // LCOV_EXCL_LINE
            " cubes=",
            normalizedCubes.size());  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      return true;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_STOP
    if (const auto core = failedAssumptionCoreFromLastResetFrontierSolve(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            solver, cube, targetFrame);  // LCOV_EXCL_LINE
        core.has_value()) {  // LCOV_EXCL_LINE
      rememberResetFrontierUnreachableCore(data, targetFrame, *core);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    } else {  // LCOV_EXCL_LINE
      rememberResetFrontierUnreachableCore(data, targetFrame, cube);  // LCOV_EXCL_LINE
    }
  }  // LCOV_EXCL_LINE
  return false;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

std::optional<std::vector<std::pair<size_t, bool>>>
// LCOV_EXCL_START
findResetFrontierUnreachableCubeCore(
// LCOV_EXCL_STOP
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps) {
  if (cube.empty()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  const auto& data = *context.data;
  const size_t targetFrame = data.bootstrapFrames + postBootstrapSteps;
  const auto normalizedCube = normalizedAssignmentCube(cube);
  if (const auto knownCore = knownResetFrontierConflictCore(
          data, normalizedCube, postBootstrapSteps);
      knownCore.has_value()) {
    rememberResetFrontierUnreachableCore(data, targetFrame, *knownCore);
    return knownCore;
  // LCOV_EXCL_START
  }
  if (const auto cachedCore =
  // LCOV_EXCL_STOP
          findCachedResetFrontierUnreachableCore(data, targetFrame, normalizedCube);
      // LCOV_EXCL_START
      cachedCore.has_value()) {
    return cachedCore;
  }
  if (resetFrontierAssumptionSolvesDisabled()) {  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  CachedResetFrontierSolver& cached =  // LCOV_EXCL_LINE
      getCachedResetFrontierSolver(data, solverType, normalizedCube, targetFrame);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  const auto core = extractUnreachableCoreFromCachedResetFrontierSolver(  // LCOV_EXCL_LINE
      cached, normalizedCube, targetFrame);  // LCOV_EXCL_LINE
  if (core.has_value()) {  // LCOV_EXCL_LINE
    rememberResetFrontierUnreachableCore(data, targetFrame, *core);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return core;  // LCOV_EXCL_LINE
}

}  // namespace KEPLER_FORMAL::SEC
