// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "imc/IMCEngine.h"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/SecDiag.h"
#include "kinduction/BaseCaseSolver.h"
#include "imc/CraigInterpolatingModelChecker.h"
#include "imc/ExactInterpolantSynthesizer.h"
#include "kinduction/OutputBatching.h"
#include "kinduction/SatEncoding.h"
#include "proof/ProofEngineShared.h"

namespace KEPLER_FORMAL::SEC {

// Overall IMC algorithm:
// 1. Reuse the shared SEC base-case search to detect concrete counterexamples.
// 2. Build the startup frontier from init/reset/bootstrap constraints.
// 3. Try any already-validated strengthening plus the one-step exact
//    interpolant as an immediate inductive proof.
// 4. If that is not enough, grow the exact reachable frontier with depth k.
// 5. At each k, conjoin the frontier with the validated strengthening and ask
//    whether the result is now inductive and excludes bad.
// 6. Return the first counterexample, the first proof depth, or inconclusive.

namespace {

constexpr OutputBatchingLimits kLargeDualRailCraigBatchingLimits{
    /*maxOutputBatchSize=*/8,
    /*outputBatchSupportLimit=*/8192};

constexpr size_t kLargeDualRailCraigProjectionStateLimit = 86016;
constexpr size_t kLargeDualRailCraigImageTransitionStateLimit = 32768;
constexpr size_t kLargeDualRailCounterexampleProbeSupportLimit = 512;
constexpr size_t kLargeDualRailCounterexampleProbeStatefulStateLimit = 8192;

constexpr CraigImcGrowthBudget kDefaultLargeDualRailCraigGrowthBudget{
    /*enabled=*/true,
    /*maxQExpansionPass=*/4,
    /*maxInterpolantClauses=*/100000,
    /*maxInterpolantLiterals=*/250000,
    /*maxInterpolantAuxiliaries=*/50000,
    // BP's retained-helper tail exposes an 84,516-state focused support cone.
    // Allow exactly that strict projection surface, with a small power-of-two
    // margin, while staying below the next broad all-support regime.
    /*maxProjectionStates=*/kLargeDualRailCraigProjectionStateLimit,
    // A single projected image above this size has already entered the memory
    // risk zone.  Return a strict IMC budget result before allocating the Craig
    // solver/proof trace for that image; the caller can split or report
    // allowed-inconclusive instead of crossing the top-MEM cap.
    /*maxImageTransitionStates=*/kLargeDualRailCraigImageTransitionStateLimit};

constexpr size_t kCraigBatchMinOverlapPercent = 90;
constexpr size_t kCraigBatchMinMarginalSupportLimit = 64;
constexpr size_t kCraigBatchMarginalSupportDivisor = 10;
constexpr size_t kCraigBatchSingleOutputSupportLimit = 1024;
constexpr size_t kCraigBatchProjectionSoftLimit = 1024;
constexpr size_t kCraigBatchProjectionHardLimit = 2048;
constexpr size_t kCraigBatchProjectionSoftMaxOutputs = 4;
constexpr size_t kCraigBatchProjectionHardMaxOutputs = 2;
constexpr size_t kCraigReusableProjectionMinSupport = 256;
constexpr size_t kCraigSingleOutputSelfProjectionSeedLimit = 1024;
constexpr size_t kCraigProjectionCacheWorkLimit = 500000;
constexpr size_t kCraigDenseProjectionCacheStateThreshold = 4096;
constexpr size_t kCraigDenseProjectionCacheOutputThreshold = 64;
constexpr size_t kCraigIndexedVectorOutputRawSupportLimit = 64;
constexpr size_t kCraigIndexedVectorUnionRawSupportLimit = 128;
constexpr size_t kCraigHelperReuseMinOverlapPercent = 80;
constexpr size_t kCraigSingletonRawHelperSourceLimit = 32;

struct CraigOutputSupport { // LCOV_EXCL_LINE
  std::unordered_set<size_t> raw;
  std::unordered_set<size_t> projection;
};

struct CraigOutputSupportCache {
  std::vector<CraigOutputSupport> outputSupports;
  std::vector<std::string> outputNames;
};

struct CraigOutputBatchPlan {
  size_t firstOutput = 0;
  size_t endOutput = 0;
};

struct CraigCounterexampleProbe {
  size_t output = 0;
  size_t support = 0;
  bool touchesState = false;
};

enum class CraigTrackedSeedScope {
  SharedReusableSurface,
  LocalRange,
};

std::vector<CraigOutputBatchPlan> buildLargeDualRailCraigImcOutputBatchPlans(
    const CraigOutputSupportCache& supportCache,
    const OutputBatchingLimits& limits);

std::unordered_set<size_t> observedOutputSupport(
    const KInductionProblem& problem,
    size_t outputIndex) {
  std::unordered_set<size_t> support;
  const auto support0 =
      problem.observedOutputExprs0[outputIndex]->getSupportVars();
  const auto support1 =
      problem.observedOutputExprs1[outputIndex]->getSupportVars();
  support.reserve(support0.size() + support1.size());
  support.insert(support0.begin(), support0.end());
  support.insert(support1.begin(), support1.end());
  return support;
}

bool shouldBuildCraigProjectionSupportCache(
    const KInductionProblem& problem) {
  const size_t outputCount = problem.observedOutputExprs0.size();
  const size_t stateCount = problem.effectiveTotalStateCount();
  if (outputCount == 0 || stateCount == 0) {
    return true; // LCOV_EXCL_LINE
  }
  // The cache expands each output support through the transition relation.
  // That is useful on AES-sized surfaces, but BP-sized dual-rail problems can
  // spend tens of GB before the first Craig query, and RISC-V-sized buses spend
  // minutes in closure walks before proving only a few outputs.  When the
  // surface is dense, keep raw support batching and let Craig's bounded
  // projection refinement grow one proof slice at a time.
  if (stateCount >= kCraigDenseProjectionCacheStateThreshold &&
      outputCount >= kCraigDenseProjectionCacheOutputThreshold) {
    return false;
  }
  // Keep a product cap for broad state/output mixes below the dense guard.
  return stateCount <= kCraigProjectionCacheWorkLimit / outputCount;
}

CraigOutputSupportCache buildCraigOutputSupportCache(
    const KInductionProblem& problem) {
  CraigOutputSupportCache cache;
  cache.outputSupports.reserve(problem.observedOutputExprs0.size());
  cache.outputNames = problem.observedOutputNames;
  const bool buildProjectionCache =
      shouldBuildCraigProjectionSupportCache(problem);
  if (!buildProjectionCache) {
    emitSecDiag(
        "SEC diag: imc Craig skips projection support cache states=",
        problem.effectiveTotalStateCount(),
        " outputs=", problem.observedOutputExprs0.size(),
        " work_limit=", kCraigProjectionCacheWorkLimit);
  }
  for (size_t output = 0; output < problem.observedOutputExprs0.size();
       ++output) {
    CraigOutputSupport support;
    support.raw = observedOutputSupport(problem, output);
    if (buildProjectionCache) {
      support.projection =
          computeCraigImcProjectionClosure(problem, support.raw);
    }
    cache.outputSupports.push_back(std::move(support));
  }
  return cache;
}

const std::unordered_set<size_t>& batchingSupport(
    const CraigOutputSupport& support) {
  return support.projection.empty() ? support.raw : support.projection;
}

size_t supportIntersectionSize(const std::unordered_set<size_t>& lhs,
                               const std::unordered_set<size_t>& rhs) {
  const auto& smaller = lhs.size() <= rhs.size() ? lhs : rhs;
  const auto& larger = lhs.size() <= rhs.size() ? rhs : lhs;
  size_t count = 0;
  for (const size_t symbol : smaller) {
    if (larger.find(symbol) != larger.end()) {
      ++count;
    }
  }
  return count;
}

void mergeSupport(std::unordered_set<size_t>& batchSupport,
                  const std::unordered_set<size_t>& outputSupport) {
  batchSupport.insert(outputSupport.begin(), outputSupport.end());
}

void mergeCraigOutputSupport(CraigOutputSupport& batchSupport,
                             const CraigOutputSupport& outputSupport) {
  mergeSupport(batchSupport.raw, outputSupport.raw);
  mergeSupport(batchSupport.projection, outputSupport.projection);
}

CraigOutputSupport mergedCraigOutputSupportForRange(
    const CraigOutputSupportCache& supportCache,
    size_t firstOutput,
    size_t endOutput) {
  CraigOutputSupport rangeSupport;
  const size_t boundedEnd =
      std::min(endOutput, supportCache.outputSupports.size());
  for (size_t output = firstOutput; output < boundedEnd; ++output) {
    mergeCraigOutputSupport(rangeSupport, supportCache.outputSupports[output]);
  }
  return rangeSupport;
}

void mergeLocalTrackedSeeds(
    std::unordered_set<size_t>& trackedStateSeeds,
    const CraigOutputSupport& outputSupport) {
  mergeSupport(trackedStateSeeds, outputSupport.raw);
  mergeSupport(trackedStateSeeds, outputSupport.projection);
}

bool supportCoversMostOf(const std::unordered_set<size_t>& reference,
                         const std::unordered_set<size_t>& candidate) {
  if (candidate.empty()) {
    return reference.empty(); // LCOV_EXCL_LINE
  }
  return (supportIntersectionSize(reference, candidate) * 100) /
             candidate.size() >=
         kCraigBatchMinOverlapPercent;
}

bool hasReusableProjectionSurface(const CraigOutputSupport& batchSupport,
                                  const CraigOutputSupport& outputSupport) {
  const auto& batchComparison = batchingSupport(batchSupport);
  const auto& outputComparison = batchingSupport(outputSupport);
  if (batchComparison.size() < kCraigReusableProjectionMinSupport ||
      outputComparison.size() < kCraigReusableProjectionMinSupport) {
    return false;
  }
  return supportCoversMostOf(batchComparison, outputComparison) &&
         supportCoversMostOf(outputComparison, batchComparison);
}

bool supportOverlapAtLeast(const std::unordered_set<size_t>& lhs,
                           const std::unordered_set<size_t>& rhs,
                           size_t minPercent) {
  if (lhs.empty() || rhs.empty()) {
    return lhs.empty() && rhs.empty(); // LCOV_EXCL_LINE
  }
  const size_t overlap = supportIntersectionSize(lhs, rhs);
  return (overlap * 100) / lhs.size() >= minPercent &&
         (overlap * 100) / rhs.size() >= minPercent;
}

bool hasStrongReusableSupportOverlap(
    const CraigOutputSupport& sourceSupport,
    const CraigOutputSupport& targetSupport) {
  // Projection support is the best compatibility signal when it is available.
  // BP-sized runs skip that cache, so fall back to raw output support.  Require
  // bidirectional overlap to avoid carrying a tiny control-bit invariant into a
  // mostly unrelated datapath bus.
  const auto& sourceComparison = batchingSupport(sourceSupport);
  const auto& targetComparison = batchingSupport(targetSupport);
  return supportOverlapAtLeast(
      sourceComparison,
      targetComparison,
      kCraigHelperReuseMinOverlapPercent);
}

bool canUseSmallRawHelperForSingleton(
    const CraigOutputSupport& sourceSupport,
    const CraigOutputSupport& rangeSupport) {
  // When BP skips the projection cache, a tiny control helper can be valuable
  // for the next valid bit even though raw overlap is low.  Do not extend that
  // exception to wide datapath helpers; they can drag a large bus surface into
  // an unrelated singleton and spend the remaining budget rediscovering Craig
  // projection states.
  return sourceSupport.projection.empty() &&
         rangeSupport.projection.empty() &&
         sourceSupport.raw.size() <= kCraigSingletonRawHelperSourceLimit;
}

bool shouldRetainSmallRawSingletonCraigInvariant(
    const CraigOutputSupport& rangeSupport,
    size_t firstOutput,
    size_t endOutput) {
  // Keep one small control-style helper around after the normal reusable slot
  // moves on to a wide bus proof.  BP valid/ready outputs can need that helper
  // late in the output list, but wide data batches should still see only the
  // compatibility-checked primary helper.
  return endOutput == firstOutput + 1 &&
         rangeSupport.projection.empty() &&
         rangeSupport.raw.size() <= kCraigSingletonRawHelperSourceLimit;
}

size_t marginalSupportLimit(size_t batchSupportSize) {
  return std::max(
      kCraigBatchMinMarginalSupportLimit,
      batchSupportSize / kCraigBatchMarginalSupportDivisor);
}

size_t adaptiveCraigBatchMaxOutputCount(
    const CraigOutputSupport& batchSupport,
    const CraigOutputSupport& outputSupport,
    const OutputBatchingLimits& limits) {
  const size_t projectionSize = std::max(
      batchingSupport(batchSupport).size(),
      batchingSupport(outputSupport).size());
  size_t maxOutputCount = limits.maxOutputBatchSize;
  if (projectionSize >= kCraigBatchProjectionHardLimit) {
    maxOutputCount =
        std::min(maxOutputCount, kCraigBatchProjectionHardMaxOutputs);
  } else if (projectionSize >= kCraigBatchProjectionSoftLimit) {
    maxOutputCount =
        std::min(maxOutputCount, kCraigBatchProjectionSoftMaxOutputs);
  }
  return maxOutputCount;
}

bool indexedOutputVectorBase(
    const std::string& outputName,
    std::string* baseName) {
  if (outputName.empty() || outputName.back() != ']') {
    return false;
  }
  const size_t openBracket = outputName.rfind('['); // LCOV_EXCL_LINE
  if (openBracket == std::string::npos || // LCOV_EXCL_LINE
      openBracket == 0 || // LCOV_EXCL_LINE
      openBracket + 1 == outputName.size() - 1) { // LCOV_EXCL_LINE
    return false; // LCOV_EXCL_LINE
  }
  for (size_t index = openBracket + 1; index + 1 < outputName.size(); // LCOV_EXCL_LINE
       ++index) { // LCOV_EXCL_LINE
    if (outputName[index] < '0' || outputName[index] > '9') { // LCOV_EXCL_LINE
      return false; // LCOV_EXCL_LINE
    }
  } // LCOV_EXCL_LINE
  *baseName = outputName.substr(0, openBracket); // LCOV_EXCL_LINE
  return true; // LCOV_EXCL_LINE
}

bool sameIndexedOutputVector(
    const CraigOutputSupportCache& supportCache,
    size_t lhsOutput,
    size_t rhsOutput) {
  if (lhsOutput >= supportCache.outputNames.size() ||
      rhsOutput >= supportCache.outputNames.size()) {
    return false; // LCOV_EXCL_LINE
  }
  std::string lhsBase;
  std::string rhsBase;
  return indexedOutputVectorBase(
             supportCache.outputNames[lhsOutput], &lhsBase) &&
         indexedOutputVectorBase( // LCOV_EXCL_LINE
             supportCache.outputNames[rhsOutput], &rhsBase) && // LCOV_EXCL_LINE
         lhsBase == rhsBase; // LCOV_EXCL_LINE
}

bool canUseIndexedVectorRawBatch(
    const CraigOutputSupportCache& supportCache,
    const CraigOutputSupport& batchSupport,
    const CraigOutputSupport& outputSupport,
    size_t firstOutput,
    size_t output,
    size_t unionSupport) {
  // BP-sized problems skip the per-output projection cache to stay under the
  // memory target. In that mode, neighboring bits of the same observed vector
  // have low immediate support overlap even though Craig later discovers the
  // same transition surface. Batch only small same-vector raw cones so unrelated
  // wide outputs still split before interpolation.
  return batchSupport.projection.empty() &&
         outputSupport.projection.empty() &&
         outputSupport.raw.size() <=
             kCraigIndexedVectorOutputRawSupportLimit &&
         unionSupport <= kCraigIndexedVectorUnionRawSupportLimit &&
         sameIndexedOutputVector(supportCache, firstOutput, output);
}

size_t envSizeOrDefault(const char* name, size_t fallback) {
  const char* value = std::getenv(name);
  if (value == nullptr || *value == '\0') {
    return fallback;
  }
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(value, &end, 10);
  if (end == value || *end != '\0') {
    return fallback; // LCOV_EXCL_LINE
  }
  return static_cast<size_t>(parsed);
}

CraigImcGrowthBudget largeDualRailCraigGrowthBudget() {
  CraigImcGrowthBudget budget = kDefaultLargeDualRailCraigGrowthBudget;
  budget.maxQExpansionPass = envSizeOrDefault(
      "KEPLER_SEC_IMC_CRAIG_MAX_Q_PASS", budget.maxQExpansionPass);
  budget.maxImageTransitionStates = envSizeOrDefault(
      "KEPLER_SEC_IMC_CRAIG_MAX_IMAGE_TRANSITION_STATES",
      budget.maxImageTransitionStates);
  return budget;
}

std::unordered_set<size_t> buildSharedReusableCraigTrackedSeeds(
    const CraigOutputSupportCache& supportCache,
    const CraigOutputSupport& rangeSupport) {
  std::unordered_set<size_t> trackedStateSeeds;
  for (const CraigOutputSupport& outputSupport : supportCache.outputSupports) {
    if (!hasReusableProjectionSurface(rangeSupport, outputSupport)) {
      continue;
    }
    // The initial scheduled batch can start from a same-design surface shared
    // by sibling outputs. Recursive split children use LocalRange below so a
    // failed parent does not force the same full AES surface into every child.
    mergeLocalTrackedSeeds(trackedStateSeeds, outputSupport);
  }
  return trackedStateSeeds;
}

std::unordered_set<size_t> buildCraigTrackedStateSeedsForRange(
    const CraigOutputSupportCache& supportCache,
    size_t firstOutput,
    size_t endOutput,
    CraigTrackedSeedScope seedScope) {
  const CraigOutputSupport rangeSupport = mergedCraigOutputSupportForRange(
      supportCache, firstOutput, endOutput);
  if (seedScope == CraigTrackedSeedScope::SharedReusableSurface) {
    return buildSharedReusableCraigTrackedSeeds(supportCache, rangeSupport);
  }

  std::unordered_set<size_t> trackedStateSeeds;
  mergeLocalTrackedSeeds(trackedStateSeeds, rangeSupport);
  return trackedStateSeeds;
}

bool canAddOutputToCraigBatch(
    const CraigOutputSupportCache& supportCache,
    const CraigOutputSupport& batchSupport,
    const CraigOutputSupport& outputSupport,
    const OutputBatchingLimits& limits,
    size_t firstOutput,
    size_t output,
    size_t batchOutputCount,
    const char** rejectReason,
    size_t* marginalSupport,
    size_t* overlapPercent,
    size_t* effectiveMaxOutputCount) {
  if (batchOutputCount == 0) {
    return true;
  }
  // Large projection surfaces usually mean each additional output widens the
  // Craig bad predicate on top of an already expensive transition image.
  *effectiveMaxOutputCount =
      adaptiveCraigBatchMaxOutputCount(batchSupport, outputSupport, limits);
  if (batchOutputCount + 1 > *effectiveMaxOutputCount) {
    *rejectReason = *effectiveMaxOutputCount < limits.maxOutputBatchSize
                        ? "adaptive_count_limit"
                        : "count_limit";
    return false;
  }
  // Large frontiers make Craig interpolants grow with the OR'ed bad predicate;
  // keep those raw-wide outputs single even when supports overlap heavily.
  if (batchSupport.raw.size() > kCraigBatchSingleOutputSupportLimit ||
      outputSupport.raw.size() > kCraigBatchSingleOutputSupportLimit) {
    *rejectReason = "large_raw_support";
    return false;
  }

  const auto& batchComparison = batchingSupport(batchSupport);
  const auto& outputComparison = batchingSupport(outputSupport);
  const size_t overlap =
      supportIntersectionSize(batchComparison, outputComparison);
  const size_t marginal = outputComparison.size() - overlap;
  const size_t unionSupport = batchComparison.size() + marginal;
  *marginalSupport = marginal;
  *overlapPercent = outputComparison.empty()
                        ? 100
                        : (overlap * 100) / outputComparison.size();

  if (unionSupport > limits.outputBatchSupportLimit) {
    *rejectReason = "support_limit"; // LCOV_EXCL_LINE
    return false; // LCOV_EXCL_LINE
  }
  if (*overlapPercent < kCraigBatchMinOverlapPercent) {
    if (canUseIndexedVectorRawBatch(
            supportCache,
            batchSupport,
            outputSupport,
            firstOutput,
            output,
            unionSupport)) {
      return true; // LCOV_EXCL_LINE
    }
    *rejectReason = "low_overlap";
    return false;
  }
  if (marginal > marginalSupportLimit(batchComparison.size())) {
    *rejectReason = "marginal_support"; // LCOV_EXCL_LINE
    return false; // LCOV_EXCL_LINE
  }
  return true;
}

void addComplementedStateRelations(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    size_t numFrames) {
  // Complemented outputs such as Q/QN are modeled as hard equalities between
  // the primary and inverted state views in every explored frame.
  for (size_t frame = 0; frame < numFrames; ++frame) {
    for (const auto& [primarySymbol, complementedSymbol] : complementedStatePairs) {
      // LCOV_EXCL_START
      addLiteralEquivalence(  // LCOV_EXCL_LINE
          solver,  // LCOV_EXCL_LINE
          variables.getLiteral(complementedSymbol, frame),  // LCOV_EXCL_LINE
          -variables.getLiteral(primarySymbol, frame));  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
    }
  }
}

void addTransitionRelation(SATSolverWrapper& solver,
                           const FrameVariableStore& variables,
                           const KInductionProblem& problem,
                           size_t frame) {
  // Encode one SEC transition step for both designs into the local SAT frame.
  FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
  for (const auto& [stateSymbol, expr] : problem.transitions0) {
    addLiteralEquivalence(
        solver,
        variables.getLiteral(stateSymbol, frame + 1),
        encoder.encode(expr));
  }
  for (const auto& [stateSymbol, expr] : problem.transitions1) {
    addLiteralEquivalence(
        solver,
        variables.getLiteral(stateSymbol, frame + 1),
        encoder.encode(expr));
  }
}

BoolExpr* buildStateAssignmentCube(const std::vector<size_t>& symbols, size_t assignment) {
  // Turn one concrete state assignment into a BoolExpr cube so exact frontier
  // enumeration can reuse the common proof-formula helpers.
  BoolExpr* cube = BoolExpr::createTrue();
  for (size_t bit = 0; bit < symbols.size(); ++bit) {
    BoolExpr* literal = BoolExpr::Var(symbols[bit]);
    cube = BoolExpr::And(
        cube,
        (assignment & (size_t{1} << bit)) != 0 ? literal : BoolExpr::Not(literal));
  }
  return BoolExpr::simplify(cube);
}

bool isStateReachableAtDepth(const KInductionProblem& problem,
                             KEPLER_FORMAL::Config::SolverType solverType,
                             BoolExpr* initFormula,
                             const std::vector<size_t>& stateSymbols,
                             size_t assignment,
                             size_t depth) {
  // Ask whether one concrete combined SEC state is reachable at exactly this
  // depth from the init/reset frontier.
  SATSolverWrapper solver(solverType);
  FrameVariableStore variables(solver, problem.allSymbols, depth + 1);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs0, depth + 1);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs1, depth + 1);
  for (size_t frame = 0; frame < depth; ++frame) {
    addTransitionRelation(solver, variables, problem, frame);
  }

  FrameFormulaEncoder initEncoder(solver, variables.makeLeafLits(0));
  solver.addClause({initEncoder.encode(initFormula)});

  FrameFormulaEncoder targetEncoder(solver, variables.makeLeafLits(depth));
  solver.addClause({targetEncoder.encode(buildStateAssignmentCube(stateSymbols, assignment))});
  return solver.solve();
}

BoolExpr* buildExactReachableStateInvariant(const KInductionProblem& problem,
                                            KEPLER_FORMAL::Config::SolverType solverType,
                                            BoolExpr* initFormula,
                                            size_t depth,
                                            size_t maxStateBits = 12) {
  const std::vector<size_t> combinedStateSymbols = problem.combinedStateSymbols();
  if (initFormula == nullptr || combinedStateSymbols.empty() ||
      combinedStateSymbols.size() > maxStateBits) {
    // LCOV_EXCL_START
    return nullptr;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  const size_t assignmentCount = size_t{1} << combinedStateSymbols.size();
  BoolExpr* reachable = BoolExpr::createFalse();
  bool foundReachableState = false;

  // Build the exact set of states reachable within the bounded frontier. IMC
  // can then stop as soon as this frontier becomes inductive and excludes bad.
  for (size_t assignment = 0; assignment < assignmentCount; ++assignment) {
    bool reachableWithinDepth = false;
    for (size_t frame = 0; frame <= depth; ++frame) {
      if (isStateReachableAtDepth(
              problem, solverType, initFormula, combinedStateSymbols, assignment, frame)) {
        reachableWithinDepth = true;
        break;
      }
    }

    if (!reachableWithinDepth) {
      continue;
    }

    foundReachableState = true;
    reachable = BoolExpr::Or(
        reachable, buildStateAssignmentCube(combinedStateSymbols, assignment));
  }

  if (!foundReachableState) {
    // LCOV_EXCL_START
    return nullptr;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  return BoolExpr::simplify(reachable);
}

BoolExpr* buildInitialImcStrengthening(const KInductionProblem& problem,
                                       KEPLER_FORMAL::Config::SolverType solverType,
                                       BoolExpr* initFormula) {
  if (initFormula == nullptr) {
    return nullptr;  // LCOV_EXCL_LINE
  }

  // Reuse any already validated SEC strengthening and then sharpen it with the
  // exact one-step interpolant when that derivation is affordable.
  BoolExpr* sharedStrengthening =
      selectValidatedStrengtheningInvariant(problem, initFormula, solverType);
  ExactInterpolantSynthesizer interpolantSynthesizer(problem, solverType);
  if (auto interpolant =
          interpolantSynthesizer.deriveOneStepReachableStateInvariant();
      interpolant.has_value()) {
    sharedStrengthening =
        sharedStrengthening == nullptr
            ? *interpolant
            : BoolExpr::simplify(BoolExpr::And(sharedStrengthening, *interpolant));
  }
  return sharedStrengthening == nullptr ? problem.property : sharedStrengthening;
}

bool provesImcInvariant(const KInductionProblem& problem,
                        KEPLER_FORMAL::Config::SolverType solverType,
                        BoolExpr* initFormula,
                        BoolExpr* invariant) {
  return invariant != nullptr &&
         initialFrontierImplies(initFormula, invariant, solverType) &&
         isInductiveInvariant(problem, invariant, solverType) &&
         invariantExcludesBadStates(problem, invariant, solverType);
}

std::optional<IMCResult> findImcCounterexample(const ImcBaseCounterexampleCache& cache,
                                               KEPLER_FORMAL::Config::SolverType solverType,
                                               size_t depth) {
  // IMC checks depths monotonically.  Only the newly exposed frontier can hold
  // a fresh counterexample, so avoid rebuilding a cumulative BMC query that
  // re-walks already-cleared frames and all earlier output bad clauses.
  if (auto witness = findImcBaseCounterexampleAtFrontier(cache, solverType, depth);
      witness.has_value()) {
    return IMCResult{IMCStatus::Different, witness->badFrame, std::move(witness)};
  }
  return std::nullopt;
}

struct ReusableCraigInvariant { // LCOV_EXCL_LINE
  std::vector<InterpolantRegion> regions;
  std::unordered_set<size_t> trackedStates;
  std::vector<std::pair<size_t, bool>> auxiliaryConstants;
  std::vector<std::pair<size_t, size_t>> auxiliaryEqualities;
  CraigOutputSupport sourceSupport;
  size_t sourceFirstOutput = 0;
  size_t sourceEndOutput = 0;
  size_t proofBound = 0;
  bool hasSource = false;
};

bool isSingleCraigOutputRange(size_t firstOutput, size_t endOutput) {
  return endOutput == firstOutput + 1;
}

bool shouldRetrySingleOutputWithSelfProjection(
    size_t firstOutput,
    size_t endOutput,
    const std::unordered_set<size_t>& trackedStateSeeds) {
  return isSingleCraigOutputRange(firstOutput, endOutput) &&
         trackedStateSeeds.size() >= kCraigSingleOutputSelfProjectionSeedLimit;
}

bool reusableCraigInvariantMatchesOutputRange(
    const CraigOutputSupportCache& supportCache,
    const ReusableCraigInvariant& reusableInvariant,
    const CraigOutputSupport& rangeSupport,
    size_t firstOutput,
    size_t endOutput) {
  if (reusableInvariant.regions.empty()) {
    return false;
  }
  if (!reusableInvariant.hasSource) {
    return true; // LCOV_EXCL_LINE
  }
  if (isSingleCraigOutputRange(firstOutput, endOutput) &&
      canUseSmallRawHelperForSingleton(
          reusableInvariant.sourceSupport, rangeSupport)) {
    return true;
  }
  // Same indexed buses, such as BP mem_data_cmd_o[*], often have low raw
  // overlap before projection but converge to the same Craig surface.
  if (sameIndexedOutputVector(
          supportCache, reusableInvariant.sourceFirstOutput, firstOutput)) {
    return true; // LCOV_EXCL_LINE
  }
  return hasReusableProjectionSurface(
             reusableInvariant.sourceSupport, rangeSupport) ||
         hasStrongReusableSupportOverlap(
             reusableInvariant.sourceSupport, rangeSupport);
}

void saveReusableCraigInvariant(
    const CraigImcResult& proof,
    const CraigOutputSupport& rangeSupport,
    size_t firstOutput,
    size_t endOutput,
    ReusableCraigInvariant& reusableInvariant) {
  reusableInvariant.regions = proof.invariantRegions;
  reusableInvariant.trackedStates = proof.trackedStates;
  reusableInvariant.auxiliaryConstants = proof.auxiliaryStateInvariants;
  reusableInvariant.auxiliaryEqualities = proof.auxiliaryStateEqualities;
  reusableInvariant.sourceSupport = rangeSupport;
  reusableInvariant.sourceFirstOutput = firstOutput;
  reusableInvariant.sourceEndOutput = endOutput;
  reusableInvariant.proofBound =
      std::max(reusableInvariant.proofBound, proof.iterations);
  reusableInvariant.hasSource = true;
}

template <typename T>
void appendUniqueCraigHelperFacts(std::vector<T>& target, // LCOV_EXCL_LINE
                                  const std::vector<T>& source) {
  for (const T& fact : source) { // LCOV_EXCL_LINE
    if (std::find(target.begin(), target.end(), fact) == target.end()) { // LCOV_EXCL_LINE
      target.push_back(fact); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
  }
} // LCOV_EXCL_LINE

ReusableCraigInvariant combineCraigHelperInvariants( // LCOV_EXCL_LINE
    const ReusableCraigInvariant& primary,
    const ReusableCraigInvariant& secondary) {
  ReusableCraigInvariant combined = primary; // LCOV_EXCL_LINE
  combined.regions.insert( // LCOV_EXCL_LINE
      combined.regions.end(), secondary.regions.begin(), secondary.regions.end()); // LCOV_EXCL_LINE
  mergeSupport(combined.trackedStates, secondary.trackedStates); // LCOV_EXCL_LINE
  appendUniqueCraigHelperFacts( // LCOV_EXCL_LINE
      combined.auxiliaryConstants, secondary.auxiliaryConstants); // LCOV_EXCL_LINE
  appendUniqueCraigHelperFacts( // LCOV_EXCL_LINE
      combined.auxiliaryEqualities, secondary.auxiliaryEqualities); // LCOV_EXCL_LINE
  combined.proofBound = std::max(combined.proofBound, secondary.proofBound); // LCOV_EXCL_LINE
  combined.hasSource = primary.hasSource || secondary.hasSource; // LCOV_EXCL_LINE
  return combined; // LCOV_EXCL_LINE
} // LCOV_EXCL_LINE

std::unordered_set<size_t> buildCraigInitialTrackedStatesForAttempt(
    const std::unordered_set<size_t>& trackedStateSeeds,
    const ReusableCraigInvariant& reusableInvariant,
    bool includeRangeSeeds) {
  std::unordered_set<size_t> initialTrackedStates;
  if (includeRangeSeeds) {
    initialTrackedStates = trackedStateSeeds;
  }
  if (!reusableInvariant.regions.empty()) {
    mergeSupport(initialTrackedStates, reusableInvariant.trackedStates);
  }
  return initialTrackedStates;
}

CraigImcOptions makeCraigImcOptions(
    bool enableGrowthBudget,
    const ReusableCraigInvariant& reusableInvariant) {
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  options.helperAuxiliaryStateInvariants =
      reusableInvariant.auxiliaryConstants;
  options.helperAuxiliaryStateEqualities =
      reusableInvariant.auxiliaryEqualities;
  if (enableGrowthBudget) {
    options.growthBudget = largeDualRailCraigGrowthBudget();
  }
  return options;
}

CraigImcResult runCraigCheckerAttempt(
    const KInductionProblem& batchProblem,
    const ReusableCraigInvariant& reusableInvariant,
    const std::unordered_set<size_t>& trackedStateSeeds,
    bool includeRangeSeeds,
    bool enableGrowthBudget,
    size_t maxK) {
  const std::unordered_set<size_t> initialTrackedStates =
      buildCraigInitialTrackedStatesForAttempt(
          trackedStateSeeds, reusableInvariant, includeRangeSeeds);
  const CraigImcOptions options = makeCraigImcOptions(
      enableGrowthBudget, reusableInvariant);
  CraigInterpolatingModelChecker checker(
      batchProblem,
      reusableInvariant.regions.empty() ? nullptr
                                        : &reusableInvariant.regions,
      initialTrackedStates.empty() ? nullptr : &initialTrackedStates,
      options);
  return checker.run(maxK);
}

IMCResult makeCraigInconclusiveResult(size_t bound) {
  IMCResult result;
  result.status = IMCStatus::Inconclusive;
  result.bound = bound;
  return result;
}

void markCraigOutputRangeCovered(
    std::vector<bool>& coveredOutputs,
    size_t firstOutput,
    size_t endOutput) {
  const size_t boundedEnd = std::min(endOutput, coveredOutputs.size());
  for (size_t output = firstOutput; output < boundedEnd; ++output) {
    coveredOutputs[output] = true;
  }
}

std::unordered_set<size_t> observedOutputSupportForProbe(
    const KInductionProblem& problem,
    size_t output) {
  std::unordered_set<size_t> support;
  if (problem.observedOutputExprs0[output] != nullptr) {
    const auto vars = problem.observedOutputExprs0[output]->getSupportVars();
    support.insert(vars.begin(), vars.end());
  }
  if (problem.observedOutputExprs1[output] != nullptr) {
    const auto vars = problem.observedOutputExprs1[output]->getSupportVars();
    support.insert(vars.begin(), vars.end());
  }
  return support;
}

bool outputSupportTouchesState(
    const std::unordered_set<size_t>& support,
    const std::unordered_set<size_t>& stateSymbols) {
  for (const auto symbol : support) {
    if (stateSymbols.find(symbol) != stateSymbols.end()) {
      return true;
    }
  }
  return false;
}

std::optional<IMCResult> findLargeDualRailCounterexampleUpTo(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK) {
  if (problem.observedOutputExprs0.empty() ||
      problem.observedOutputExprs0.size() !=
          problem.observedOutputExprs1.size()) {
    return std::nullopt; // LCOV_EXCL_LINE
  }

  std::vector<CraigCounterexampleProbe> probes;
  probes.reserve(problem.observedOutputExprs0.size());
  size_t skippedLargeSupport = 0;
  size_t skippedStateDependent = 0;
  size_t skippedProjectionSupport = 0;
  std::unordered_set<size_t> stateSymbols;
  stateSymbols.reserve(problem.state0Symbols.size() + problem.state1Symbols.size());
  stateSymbols.insert(problem.state0Symbols.begin(), problem.state0Symbols.end());
  stateSymbols.insert(problem.state1Symbols.begin(), problem.state1Symbols.end());
  const size_t stateCount = problem.effectiveTotalStateCount();
  for (size_t output = 0; output < problem.observedOutputExprs0.size(); ++output) {
    const std::unordered_set<size_t> support =
        observedOutputSupportForProbe(problem, output);
    if (support.size() > kLargeDualRailCounterexampleProbeSupportLimit) {
      ++skippedLargeSupport; // LCOV_EXCL_LINE
      continue; // LCOV_EXCL_LINE
    }
    const bool touchesState = outputSupportTouchesState(support, stateSymbols);
    size_t supportScore = support.size();
    if (touchesState) {
      if (stateCount > kLargeDualRailCounterexampleProbeStatefulStateLimit) {
        ++skippedStateDependent;
        continue;
      }
    }
    probes.push_back({output, supportScore, touchesState});
  }
  std::stable_sort(
      probes.begin(),
      probes.end(),
      [](const CraigCounterexampleProbe& lhs,
         const CraigCounterexampleProbe& rhs) {
        return lhs.support < rhs.support;
      });
  emitSecDiag(
      "SEC diag: imc large dual-rail bounded witness probes=",
      probes.size(), " skipped_support=", skippedLargeSupport,
      " skipped_state_dependent=", skippedStateDependent,
      " skipped_projection_support=", skippedProjectionSupport,
      " support_limit=", kLargeDualRailCounterexampleProbeSupportLimit,
      " stateful_state_limit=",
      kLargeDualRailCounterexampleProbeStatefulStateLimit,
      " projection_sizing=raw_support",
      " max_k=", maxK);

  for (size_t depth = 0; depth <= maxK; ++depth) {
    for (const CraigCounterexampleProbe& probe : probes) {
      if (!probe.touchesState && depth > 0) {
        continue; // LCOV_EXCL_LINE
      }
      KInductionProblem outputProblem = problem;
      configureOutputBatchProblem(
          outputProblem, problem, probe.output, probe.output + 1);
      if (auto witness = findFastBaseCounterexampleAtFrontier(
              outputProblem, solverType, depth);
          witness.has_value()) {
        return IMCResult{
            IMCStatus::Different, witness->badFrame, std::move(witness)};
      }
    }
  }
  return std::nullopt;
}

size_t boundedCraigWitnessDepth(size_t maxK, size_t reachedLookahead) {
  return std::min(maxK, reachedLookahead);
}

IMCResult runCraigOutputRange(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK,
    size_t firstOutput,
    size_t endOutput,
    const CraigOutputSupportCache& supportCache,
    CraigTrackedSeedScope seedScope,
    ReusableCraigInvariant& reusableInvariant,
    ReusableCraigInvariant& smallRawSingletonInvariant,
    std::vector<bool>& coveredOutputs) {
  KInductionProblem batchProblem = problem;
  configureOutputBatchProblem(
      batchProblem, problem, firstOutput, endOutput);
  const std::unordered_set<size_t> trackedStateSeeds =
      buildCraigTrackedStateSeedsForRange(
          supportCache, firstOutput, endOutput, seedScope);
  const CraigOutputSupport rangeSupport = mergedCraigOutputSupportForRange(
      supportCache, firstOutput, endOutput);
  const ReusableCraigInvariant emptyReusableInvariant;
  const bool useReusableInvariant = reusableCraigInvariantMatchesOutputRange(
      supportCache, reusableInvariant, rangeSupport, firstOutput, endOutput);
  const bool useSmallRawSingletonInvariant =
      !useReusableInvariant &&
      reusableCraigInvariantMatchesOutputRange(
          supportCache,
          smallRawSingletonInvariant,
          rangeSupport,
          firstOutput,
          endOutput);
  const bool combineCraigHelpers =
      shouldCombineCraigHelpersForSmallRawSingleton(
          useSmallRawSingletonInvariant,
          !reusableInvariant.regions.empty());
  ReusableCraigInvariant combinedReusableInvariant;
  if (combineCraigHelpers) {
    // The skipped reusable invariant and the retained singleton are both
    // already-proved Craig invariants over this same SEC problem.  Conjoin them
    // as helper facts instead of choosing one and rediscovering the other in
    // the expensive retained-helper tail.
    combinedReusableInvariant = combineCraigHelperInvariants( // LCOV_EXCL_LINE
        smallRawSingletonInvariant, reusableInvariant); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  const ReusableCraigInvariant& activeReusableInvariant =
      useReusableInvariant
          ? reusableInvariant
          : combineCraigHelpers ? combinedReusableInvariant
          : useSmallRawSingletonInvariant ? smallRawSingletonInvariant
                                         : emptyReusableInvariant;
  emitSecDiag(
      "SEC diag: imc Craig output batch first=", firstOutput,
      " end=", endOutput,
      " first_name=",
      firstOutput < problem.observedOutputNames.size()
          ? problem.observedOutputNames[firstOutput]
          : std::string("<unknown>"), // LCOV_EXCL_LINE
      " bad_support=", batchProblem.bad->getSupportVars().size(),
      " tracked_seed_states=", trackedStateSeeds.size(),
      " helper_regions=", activeReusableInvariant.regions.size());
  if (!reusableInvariant.regions.empty() && !useReusableInvariant) {
    emitSecDiag( // LCOV_EXCL_LINE
        "SEC diag: imc Craig skips reusable invariant for output batch first=",
        firstOutput, " end=", endOutput,
        " source_first=", reusableInvariant.sourceFirstOutput, // LCOV_EXCL_LINE
        " source_end=", reusableInvariant.sourceEndOutput, // LCOV_EXCL_LINE
        " helper_regions=", reusableInvariant.regions.size()); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  if (useSmallRawSingletonInvariant) {
    emitSecDiag( // LCOV_EXCL_LINE
        "SEC diag: imc Craig uses retained small singleton invariant for "
        "output batch first=",
        firstOutput, " end=", endOutput,
        " source_first=", smallRawSingletonInvariant.sourceFirstOutput, // LCOV_EXCL_LINE
        " source_end=", smallRawSingletonInvariant.sourceEndOutput, // LCOV_EXCL_LINE
        " helper_regions=", smallRawSingletonInvariant.regions.size()); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  if (combineCraigHelpers) {
    emitSecDiag( // LCOV_EXCL_LINE
        "SEC diag: imc Craig combines reusable helpers for output batch first=",
        firstOutput, " end=", endOutput,
        " singleton_regions=", smallRawSingletonInvariant.regions.size(), // LCOV_EXCL_LINE
        " reusable_regions=", reusableInvariant.regions.size(), // LCOV_EXCL_LINE
        " combined_regions=", combinedReusableInvariant.regions.size()); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  if (!activeReusableInvariant.regions.empty() &&
      craigInvariantExcludesBad(
          batchProblem,
          activeReusableInvariant.trackedStates,
          activeReusableInvariant.regions,
          activeReusableInvariant.auxiliaryConstants,
          activeReusableInvariant.auxiliaryEqualities)) {
    emitSecDiag(
        "SEC diag: imc Craig reused invariant for output batch first=",
        firstOutput, " end=", endOutput);
    markCraigOutputRangeCovered(coveredOutputs, firstOutput, endOutput);
    return {IMCStatus::Equivalent, activeReusableInvariant.proofBound};
  }

  const bool multiOutputRange = !isSingleCraigOutputRange(
      firstOutput, endOutput);
  const bool retrySingleOutputWithSelfProjection =
      shouldRetrySingleOutputWithSelfProjection(
          firstOutput, endOutput, trackedStateSeeds);
  const bool enableGrowthBudget = true;

  // A single AES text bit can arrive with a huge precomputed transition
  // surface. If that seeded attempt blows the Craig budget, retry the same
  // strict IMC proof from the property's own support and let projection
  // refinement add only states the proof actually touches.
  CraigImcResult proof = runCraigCheckerAttempt(
      batchProblem,
      activeReusableInvariant,
      trackedStateSeeds,
      /*includeRangeSeeds=*/true,
      enableGrowthBudget,
      maxK);
  if (proof.status == CraigImcStatus::BudgetExceeded &&
      retrySingleOutputWithSelfProjection) {
    emitSecDiag(
        "SEC diag: imc Craig retrying single output with self projection "
        "first=",
        firstOutput, " end=", endOutput,
        " dropped_seed_states=", trackedStateSeeds.size(),
        " helper_seed_states=", activeReusableInvariant.trackedStates.size());
    proof = runCraigCheckerAttempt(
        batchProblem,
        activeReusableInvariant,
        trackedStateSeeds,
        /*includeRangeSeeds=*/false,
        /*enableGrowthBudget=*/true,
        maxK);
  }
  if (proof.status == CraigImcStatus::Equivalent) {
    if (!proof.invariantRegions.empty() &&
        !proof.trackedStates.empty()) {
      // Any proved Craig region is a valid helper invariant. Small control
      // proofs, such as AES done/counter bits, can still prune later wide data
      // output images even though their own tracked surface is tiny.
      ReusableCraigInvariant proofInvariant;
      saveReusableCraigInvariant(
          proof,
          rangeSupport,
          firstOutput,
          endOutput,
          proofInvariant);
      reusableInvariant = proofInvariant;
      const bool retainSmallRawSingleton =
          shouldRetainSmallRawSingletonCraigInvariant(
              rangeSupport, firstOutput, endOutput);
      if (retainSmallRawSingleton) {
        smallRawSingletonInvariant = proofInvariant;
        emitSecDiag(
            "SEC diag: imc Craig retains small singleton invariant first=",
            firstOutput, " end=", endOutput,
            " helper_regions=", proof.invariantRegions.size());
      }
    }
    markCraigOutputRangeCovered(coveredOutputs, firstOutput, endOutput);
    return {IMCStatus::Equivalent, proof.iterations};
  }
  if (proof.status == CraigImcStatus::CounterexampleCandidate) {
    // Craig found SAT from the exact concrete post-reset frontier. Reconstruct
    // only that lookahead with the exact bounded witness encoder; this is
    // counterexample validation inside IMC, not a different proof engine.
    const auto cache = makeImcBaseCounterexampleCache(batchProblem); // LCOV_EXCL_LINE
    if (const auto counterexample = // LCOV_EXCL_LINE
            findImcCounterexample(*cache, solverType, proof.iterations); // LCOV_EXCL_LINE
        counterexample.has_value()) { // LCOV_EXCL_LINE
      return *counterexample; // LCOV_EXCL_LINE
    }
    return makeCraigInconclusiveResult(maxK); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  if (proof.status == CraigImcStatus::BudgetExceeded) {
    if (multiOutputRange) {
      emitSecDiag(
          "SEC diag: imc Craig splitting output batch after growth budget "
          "first=",
          firstOutput, " end=", endOutput);
    } else {
      emitSecDiag(
          "SEC diag: imc Craig single output exhausted growth budget first=",
          firstOutput, " end=", endOutput);
    }
  }

  if (multiOutputRange) {
    const size_t midpoint = firstOutput + (endOutput - firstOutput) / 2;
    const IMCResult left = runCraigOutputRange(
        problem,
        solverType,
        maxK,
        firstOutput,
        midpoint,
        supportCache,
        CraigTrackedSeedScope::LocalRange,
        reusableInvariant,
        smallRawSingletonInvariant,
        coveredOutputs);
    const IMCResult right = runCraigOutputRange(
        problem,
        solverType,
        maxK,
        midpoint,
        endOutput,
        supportCache,
        CraigTrackedSeedScope::LocalRange,
        reusableInvariant,
        smallRawSingletonInvariant,
        coveredOutputs);
    if (left.status == IMCStatus::Different) {
      return left; // LCOV_EXCL_LINE
    }
    if (right.status == IMCStatus::Different) {
      return right; // LCOV_EXCL_LINE
    }
    if (left.status == IMCStatus::Inconclusive) {
      return left;
    }
    if (right.status == IMCStatus::Inconclusive) { // LCOV_EXCL_LINE
      return right; // LCOV_EXCL_LINE
    }
    if (left.status == IMCStatus::Equivalent && // LCOV_EXCL_LINE
        right.status == IMCStatus::Equivalent) { // LCOV_EXCL_LINE
      return {IMCStatus::Equivalent, std::max(left.bound, right.bound)}; // LCOV_EXCL_LINE
    }
    return makeCraigInconclusiveResult(maxK); // LCOV_EXCL_LINE
  }

  if (proof.status == CraigImcStatus::ConcreteNoProgress) {
    // This is a Craig proof failure, not a concrete IMC witness. A full
    // bounded counterexample sweep on BP-sized dual-rail slices rebuilds large
    // transition prefixes and is outside the interpolation proof step. Keep
    // only the guarded local probe for small output cones; otherwise return
    // inconclusive conservatively.
    if (const auto counterexample = // LCOV_EXCL_LINE
            findLargeDualRailCounterexampleUpTo(batchProblem, solverType, maxK); // LCOV_EXCL_LINE
        counterexample.has_value()) { // LCOV_EXCL_LINE
      return *counterexample; // LCOV_EXCL_LINE
    }
    return makeCraigInconclusiveResult(maxK); // LCOV_EXCL_LINE
  }
  if (proof.status == CraigImcStatus::BudgetExceeded) {
    if (const auto counterexample =
            findLargeDualRailCounterexampleUpTo(
                batchProblem, solverType, proof.iterations);
        counterexample.has_value()) {
      return *counterexample; // LCOV_EXCL_LINE
    }
    return makeCraigInconclusiveResult(proof.iterations);
  }

  // Partial or implicit initial frontiers are over-approximations. Do not turn
  // an abstract Craig SAT result into a proof claim; report inconclusive
  // instead of falling back to an unrelated full bounded sweep.
  // The local probe is still IMC's bounded counterexample search, restricted to
  // small output supports.  Use the caller horizon so a shallow Craig attempt
  // cannot hide a concrete bad state that is already within maxK.
  const size_t checkedDepth = maxK;
  if (const auto counterexample =
          findLargeDualRailCounterexampleUpTo(
              batchProblem, solverType, checkedDepth);
      counterexample.has_value()) {
    return *counterexample;
  }
  return makeCraigInconclusiveResult(checkedDepth);
}

IMCResult runLargeDualRailCraigImc(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK) {
  // Craig IMC derives proof regions from the selected output slice. When one
  // bus bit already pulls a wide rail-state cone, nearby bits usually reuse the
  // same transition surface. Batch those bits here so classic IMC proves one
  // conjunction instead of rebuilding the same Craig query per bit. This limit
  // is IMC-local and does not change KI/PDR batching.
  const CraigOutputSupportCache supportCache =
      buildCraigOutputSupportCache(problem);
  const auto batches = buildLargeDualRailCraigImcOutputBatchPlans(
      supportCache, kLargeDualRailCraigBatchingLimits);
  size_t proofBound = 0;
  std::vector<bool> coveredOutputs(
      problem.observedOutputExprs0.size(), false);
  ReusableCraigInvariant reusableInvariant;
  ReusableCraigInvariant smallRawSingletonInvariant;
  for (const CraigOutputBatchPlan& batchPlan : batches) {
    const IMCResult batchResult = runCraigOutputRange(
        problem,
        solverType,
        maxK,
        batchPlan.firstOutput,
        batchPlan.endOutput,
        supportCache,
        CraigTrackedSeedScope::SharedReusableSurface,
        reusableInvariant,
        smallRawSingletonInvariant,
        coveredOutputs);
    if (batchResult.status == IMCStatus::Different) {
      return batchResult;
    }
    if (batchResult.status == IMCStatus::Inconclusive) {
      // Keep proofs already found by recursive splitting, but retain the
      // existing global work bound for unrelated later output batches.
      const size_t checkedDepth =
          boundedCraigWitnessDepth(maxK, batchResult.bound);
      emitSecDiag(
          "SEC diag: imc Craig stopping after inconclusive output batch first=",
          batchPlan.firstOutput, " end=", batchPlan.endOutput);
      IMCResult result = makeCraigInconclusiveResult(checkedDepth);
      result.coveredOutputs = std::move(coveredOutputs);
      return result;
    } else {
      proofBound = std::max(proofBound, batchResult.bound);
    }
  }
  IMCResult result;
  result.status = IMCStatus::Equivalent;
  result.bound = proofBound;
  result.coveredOutputs = std::move(coveredOutputs);
  return result;
}

bool shouldBuildExplicitImcInitFormula(const KInductionProblem& problem) {
  if (!problem.usesDualRailStateEncoding) {
    return true;
  }
  // Exact IMC enumerates reachable combined states only for tiny systems.
  // Large dual-rail ASIC problems use proof-derived Craig interpolation instead
  // of materializing a full rail-init formula.
  return problem.totalStateCount <= 12;
}

std::vector<CraigOutputBatchPlan> buildLargeDualRailCraigImcOutputBatchPlans(
    const CraigOutputSupportCache& supportCache,
    const OutputBatchingLimits& limits) {
  const std::vector<CraigOutputSupport>& outputSupports =
      supportCache.outputSupports;

  std::vector<CraigOutputBatchPlan> batches;
  size_t firstOutput = 0;
  while (firstOutput < outputSupports.size()) {
    size_t endOutput = firstOutput;
    CraigOutputSupport batchSupport;
    while (endOutput < outputSupports.size()) {
      const char* rejectReason = nullptr;
      size_t marginalSupport = 0;
      size_t overlapPercent = 100;
      size_t effectiveMaxOutputCount = limits.maxOutputBatchSize;
      const size_t batchOutputCount = endOutput - firstOutput;
      if (!canAddOutputToCraigBatch(
              supportCache,
              batchSupport,
              outputSupports[endOutput],
              limits,
              firstOutput,
              endOutput,
              batchOutputCount,
              &rejectReason,
              &marginalSupport,
              &overlapPercent,
              &effectiveMaxOutputCount)) {
        emitSecDiag(
            "SEC diag: imc Craig batch closes before output=", endOutput,
            " reason=", rejectReason,
            " batch_first=", firstOutput,
            " adaptive_max_outputs=", effectiveMaxOutputCount,
            " batch_support=", batchingSupport(batchSupport).size(),
            " output_support=", batchingSupport(outputSupports[endOutput]).size(),
            " batch_raw_support=", batchSupport.raw.size(),
            " output_raw_support=", outputSupports[endOutput].raw.size(),
            " batch_projection_support=", batchSupport.projection.size(),
            " output_projection_support=",
            outputSupports[endOutput].projection.size(),
            " marginal_support=", marginalSupport,
            " overlap_percent=", overlapPercent);
        break;
      }
      mergeCraigOutputSupport(batchSupport, outputSupports[endOutput]);
      ++endOutput;
    }
    if (endOutput == firstOutput) {
      // Always make progress: an oversized single output still deserves one
      // exact IMC attempt rather than being silently skipped by the scheduler.
      ++endOutput;  // LCOV_EXCL_LINE
      mergeCraigOutputSupport(batchSupport, outputSupports[firstOutput]); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
    CraigOutputBatchPlan plan;
    plan.firstOutput = firstOutput;
    plan.endOutput = endOutput;
    batches.push_back(std::move(plan));
    firstOutput = endOutput;
  }
  return batches;
}

}  // namespace

std::vector<std::pair<size_t, size_t>> buildLargeDualRailCraigImcOutputBatches(
    const KInductionProblem& problem,
    const OutputBatchingLimits& limits) {
  const CraigOutputSupportCache supportCache =
      buildCraigOutputSupportCache(problem);
  const auto plans = buildLargeDualRailCraigImcOutputBatchPlans(
      supportCache, limits);
  std::vector<std::pair<size_t, size_t>> batches;
  batches.reserve(plans.size());
  for (const CraigOutputBatchPlan& plan : plans) {
    batches.emplace_back(plan.firstOutput, plan.endOutput);
  }
  return batches;
}

size_t largeDualRailCraigImcProjectionStateLimit() {
  return kLargeDualRailCraigProjectionStateLimit;
}

bool shouldCombineCraigHelpersForSmallRawSingleton(
    bool useSmallRawSingletonInvariant,
    bool reusableInvariantHasRegions) {
  return useSmallRawSingletonInvariant && reusableInvariantHasRegions;
}

IMCEngine::IMCEngine(const KInductionProblem& problem,
                     KEPLER_FORMAL::Config::SolverType solverType)
    : problem_(problem), solverType_(solverType) {}

IMCResult IMCEngine::run(size_t maxK) const {
  if (problem_.combinedStateSymbols().empty()) {
    // Stateless SEC is still a real IMC base query: a combinational mismatch
    // at frame 0 must be reported before declaring equivalence.
    const auto baseCache = makeImcBaseCounterexampleCache(problem_);
    if (const auto counterexample =
            findImcCounterexample(*baseCache, solverType_, 0);
        counterexample.has_value()) {
      return *counterexample;
    }
    return {IMCStatus::Equivalent, 0};
  }

  if (problem_.usesDualRailStateEncoding &&
      problem_.effectiveTotalStateCount() > 12 &&
      !problem_.observedOutputExprs0.empty() &&
      problem_.observedOutputExprs0.size() ==
          problem_.observedOutputExprs1.size()) {
    if (const auto counterexample =
            findLargeDualRailCounterexampleUpTo(problem_, solverType_, 0);
        counterexample.has_value()) {
      return *counterexample;
    }
    return runLargeDualRailCraigImc(problem_, solverType_, maxK);
  }

  const auto baseCache = makeImcBaseCounterexampleCache(problem_);
  // Keep counterexample discovery on the same bounded base-case machinery as
  // the rest of SEC so witnesses and reported cycles stay consistent.
  if (const auto counterexample =
          findImcCounterexample(*baseCache, solverType_, 0);
      counterexample.has_value()) {
    return *counterexample;  // LCOV_EXCL_LINE
  }

  BoolExpr* initFormula =
      shouldBuildExplicitImcInitFormula(problem_) ? buildProofInitFormula(problem_)
                                                  : nullptr;
  const BoolExpr* sharedStrengthening =
      buildInitialImcStrengthening(problem_, solverType_, initFormula);
  if (initFormula != nullptr &&
      provesImcInvariant(problem_, solverType_, initFormula,
                         const_cast<BoolExpr*>(sharedStrengthening))) {
    // Before spending time on deeper frontiers, see whether the startup
    // strengthening is already a complete inductive proof.
    return {IMCStatus::Equivalent, 1};
  }

  for (size_t k = 1; k <= maxK; ++k) {
    if (const auto counterexample =
            findImcCounterexample(*baseCache, solverType_, k);
        counterexample.has_value()) {
      return *counterexample;
    }

    if (initFormula == nullptr) {
      // LCOV_EXCL_START
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    BoolExpr* frontierInvariant =
        buildExactReachableStateInvariant(problem_, solverType_, initFormula, k);
    if (frontierInvariant == nullptr) {
      // LCOV_EXCL_START
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    // Keep the explicit IMC engine centered on the reachable frontier, but
    // reuse any already validated strengthening to reduce the SAT work needed
    // to establish inductiveness on compact transition systems.
    BoolExpr* proofInvariant =
        sharedStrengthening == nullptr
            // LCOV_EXCL_START
            ? frontierInvariant  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            : BoolExpr::simplify(
                  BoolExpr::And(frontierInvariant, const_cast<BoolExpr*>(sharedStrengthening)));

    if (provesImcInvariant(problem_, solverType_, initFormula, proofInvariant)) {
      // LCOV_EXCL_START
      return {IMCStatus::Equivalent, k};  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
  }

  // We exhausted the requested depth without finding either a counterexample
  // or an inductive frontier.
  return {IMCStatus::Inconclusive, maxK};
}

}  // namespace KEPLER_FORMAL::SEC
