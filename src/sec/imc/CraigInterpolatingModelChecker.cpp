// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "imc/CraigInterpolatingModelChecker.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory_resource>
#include <optional>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/SecDiag.h"
#include "common/BoolExprUtils.h"
#include "kinduction/SatEncoding.h"
#include "proof/TransitionExprResolver.h"

namespace KEPLER_FORMAL::SEC {

namespace {

using VariablePartition = SATSolverWrapper::CraigVariablePartition;
using ClausePartition = SATSolverWrapper::CraigClausePartition;
using SteadyClock = std::chrono::steady_clock;

constexpr size_t kAuxiliaryInvariantSupportLimit = 256;
constexpr size_t kAuxiliaryInvariantCandidateLimit = 10000;
constexpr size_t kLocalAuxiliaryInvariantCandidateLimit = 4096;
constexpr size_t kLocalAuxiliaryInvariantLargeCandidateThreshold = 8192;
constexpr size_t kLocalAuxiliaryInvariantLargeCandidateLimit = 1024;
constexpr size_t kLocalAuxiliaryInvariantSupportMiningLimit = 65536;
constexpr size_t kLocalAuxiliaryPreservationScoreBonus = 1ULL << 32;
constexpr size_t kLocalAuxiliaryPreservationScoreSupportLimit = 512;
constexpr size_t kLocalAuxiliaryPreservedCandidateExtraLimit = 512;
constexpr size_t kAuxiliaryEqualityCandidateLimit = 1024;
constexpr size_t kCraigSemanticSimplifyClauseLimit = 256;
constexpr size_t kCraigSemanticSimplifyVariableLimit = 128;
constexpr size_t kCraigSubsumptionClauseLimit = 4096;
constexpr size_t kCraigRegionCompactionStart = 8;
constexpr size_t kCraigRegionCompactionCandidateLimit = 4;
constexpr size_t kCraigModelGuidedProjectionRefinementLimit = 64;
constexpr size_t kCraigLargeProjectionRefinementCandidateThreshold = 8192;
constexpr size_t kCraigLargeProjectionRefinementLimit = 512;
constexpr size_t kCraigTinyModelGuidedBackfillSupportThreshold =
    kCraigLargeProjectionRefinementCandidateThreshold;
constexpr size_t kCraigLowScoreBackfillSupportThreshold = 100000;
constexpr size_t kCraigLowScoreBackfillScoreLimit = 64;
constexpr size_t kCraigLowScoreProjectionRefinementLimit = 128;
constexpr size_t kCraigTightLowScoreBackfillSupportThreshold = 104000;
constexpr size_t kCraigTightLowScoreProjectionRefinementLimit = 64;
constexpr size_t kCraigHighSupportRefinementThreshold = 65536;
constexpr size_t kCraigHighSupportProjectionRefinementLimit = 4096;
constexpr size_t kCraigVeryHighSupportRefinementThreshold = 90000;
constexpr size_t kCraigVeryHighSupportHighScoreProjectionRefinementLimit = 1024;
constexpr size_t kCraigVeryHighSupportProjectionRefinementLimit = 128;
constexpr size_t kCraigDirectCdclProjectionQueryStateLimit = 4096;
constexpr size_t kCraigLocalAuxiliaryRetryLimit = 2;
constexpr size_t kCraigProjectedLocalSemanticsSupportLimit = 100000;
constexpr size_t kCraigProjectedLocalSemanticsTrackedStateLimit = 4096;
constexpr size_t kCraigFocusedImageTransitionTrackedStateLimit = 1024;
constexpr size_t kCraigFocusedImageTransitionRequestLimit = 12000;
constexpr size_t kCraigRetainedHelperFocusedImageTransitionStateThreshold =
    48000;
constexpr size_t kCraigRetainedHelperFocusedImageTransitionRequestLimit = 10000;
constexpr size_t
    kCraigBroadRetainedHelperFocusedImageTransitionStateThreshold = 48000;
constexpr size_t
    kCraigBroadRetainedHelperFocusedImageTransitionRegionThreshold = 6;
constexpr size_t
    kCraigBroadRetainedHelperFocusedImageTransitionRequestLimit = 8192;
constexpr size_t kCraigFocusedProjectionRefinementSupportLimit = 65536;
constexpr size_t kCraigFocusedProjectionRefinementLimit = 4096;
constexpr size_t kCraigProjectedTransitionBuildSupportLimit = 131072;
constexpr size_t kCraigProjectedTransitionBuildTargetLimit = 12000;
constexpr size_t kCraigProjectedTransitionBuildNodeLimit = 1000000;
constexpr size_t kCraigMaxSolverTseitinReserveHint = 65536;
constexpr size_t kCraigFocusedProjectionBulkSupportLimit = 49152;
constexpr size_t kCraigFocusedProjectionBulkCandidateLimit = 32768;
// Saturated focused projections are a last chance after projection stops
// moving.  Keep the extra Q expansion bounded; BP's q7 proof alone can cross
// the 10GB target, so advance lookahead after six retained regions.
constexpr size_t kCraigFocusedSaturatedQExpansionPassLimit = 6;
constexpr size_t kCraigRetainedHelperFocusedQExpansionPassLimit = 3;
constexpr size_t kCraigRetainedHelperFocusedProjectionQExpansionPassLimit = 3;
// Advancing lookahead discards the current interpolant instead of adding it to
// Q.  Permit a larger one-time proof than the normal retained-region budget,
// but still stop before runaway traces threaten the memory target.
constexpr size_t kCraigFocusedLookaheadAdvanceMaxClauses = 750000;
constexpr size_t kCraigFocusedLookaheadAdvanceMaxLiterals = 1750000;
constexpr size_t kCraigFocusedLookaheadAdvanceMaxAuxiliaries = 250000;
// The BP tail showed that waiting until the saturated q-pass limit can build a
// second huge proof and cross the memory target.  Once a focused proof reaches
// this size, try the next strict-IMC lookahead before growing Q again.
constexpr size_t kCraigFocusedLookaheadAdvanceMinClauses = 200000;
constexpr size_t kCraigFocusedLookaheadAdvanceMinLiterals = 500000;
constexpr size_t kCraigFocusedLookaheadAdvanceMinAuxiliaries = 70000;
constexpr size_t kCraigHelperAuxiliaryCarryLimit = 4096;
constexpr size_t kCraigRetainedHelperLocalAuxiliarySkipRegions = 6;
constexpr size_t kCraigRetainedHelperLocalAuxiliarySkipStateThreshold = 32768;
constexpr size_t kCraigNearSaturatedProjectionRemainderLimit = 256;

struct AuxiliaryStateInvariants {
  std::vector<std::pair<size_t, bool>> constants;
  std::vector<std::pair<size_t, size_t>> equalities;

  bool empty() const {
    return constants.empty() && equalities.empty();
  }
};

struct IndexedStatePairSemantics {
  size_t rhs = 0;
  bool complemented = false;
};

struct IndexedDualRailSemantics {
  size_t mayBeZero = 0;
};

struct StateSemanticsIndex {
  std::unordered_map<size_t, std::vector<IndexedStatePairSemantics>>
      pairSemanticsByLhs;
  std::unordered_map<size_t, std::vector<IndexedDualRailSemantics>>
      dualRailSemanticsByMayOne;
};

struct CraigProblemStaticIndex {
  std::unordered_set<size_t> states;
  StateSemanticsIndex stateSemantics;
};

struct ProjectedTransitionFrame {
  std::unordered_set<size_t> requests;
  std::vector<size_t> targets;
  size_t expressionNodes = 0;
};

struct StateConstantPreservationKey {
  BoolExpr* expr = nullptr;
  bool value = false;

  bool operator==(const StateConstantPreservationKey& other) const {
    return expr == other.expr && value == other.value;
  }
};

struct StateConstantPreservationKeyHash {
  size_t operator()(const StateConstantPreservationKey& key) const {
    size_t hash = std::hash<BoolExpr*>{}(key.expr);
    if (key.value) {
      hash ^= 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
    return hash;
  }
};

using StateConstantPreservationCache =
    std::unordered_map<
        StateConstantPreservationKey,
        bool,
        StateConstantPreservationKeyHash>;

struct ProjectedTransitionBuildEstimate {
  size_t stateSupport = 0;
  size_t encodedTargets = 0;
  size_t expressionNodes = 0;
};

struct ProjectedTransitionPlan {
  std::vector<ProjectedTransitionFrame> frames;
  ProjectedTransitionBuildEstimate estimate;
  size_t largestTransitionRequestCount = 0;
};

using ProjectedTransitionPlanCache =
    std::unordered_map<size_t, ProjectedTransitionPlan>;

struct RegionVariableKey {
  bool isState = false;
  size_t index = 0;

  bool operator==(const RegionVariableKey& other) const {
    return isState == other.isState && index == other.index;
  }
};

struct RegionVariableKeyHash {
  size_t operator()(const RegionVariableKey& key) const {
    return std::hash<size_t>{}((key.index << 1) ^ (key.isState ? 1u : 0u));
  }
};

struct ExprBootstrapValueKey {
  BoolExpr* expr = nullptr;
  bool value = false;

  bool operator==(const ExprBootstrapValueKey& other) const {
    return expr == other.expr && value == other.value;
  }
};

struct ExprBootstrapValueKeyHash {
  size_t operator()(const ExprBootstrapValueKey& key) const {
    return std::hash<BoolExpr*>{}(key.expr) ^
           (key.value ? 0x9e3779b97f4a7c15ULL : 0ULL);
  }
};

struct ProjectionRefinementCandidate {
  size_t symbol = 0;
  size_t score = 0;
};

struct FrontierResult;

std::unordered_map<size_t, size_t> primaryByComplement(
    const KInductionProblem& problem);

size_t transitionTargetFor(
    size_t symbol,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary);

bool projectionRefinementCandidateBetter(
    const ProjectionRefinementCandidate& lhs,
    const ProjectionRefinementCandidate& rhs) {
  if (lhs.score != rhs.score) {
    return lhs.score > rhs.score;
  }
  return lhs.symbol < rhs.symbol;
}

size_t localAuxiliaryCandidateLimit(size_t candidateCount) {
  // Large dual-rail cones can expose tens of thousands of local bootstrap
  // constants.  Keep the full slice for moderate supports, but switch to the
  // highest-scored quarter slice once the proof cost is dominated by thousands
  // of per-candidate preservation SAT calls.
  return candidateCount > kLocalAuxiliaryInvariantLargeCandidateThreshold
             ? kLocalAuxiliaryInvariantLargeCandidateLimit
             : kLocalAuxiliaryInvariantCandidateLimit;
}

size_t largeSupportProjectionRefinementLimit(size_t supportSize) {
  // Once the support pool is BP-sized, many 64-symbol projection rounds spend
  // more time rebuilding Craig queries than proving.  Take a wider scored
  // slice while still avoiding the old all-support import that caused blowups.
  return supportSize > kCraigLargeProjectionRefinementCandidateThreshold
             ? kCraigLargeProjectionRefinementLimit
             : kCraigModelGuidedProjectionRefinementLimit;
}

size_t boundedProjectionRefinementLimit(size_t candidateCount) {
  return largeSupportProjectionRefinementLimit(candidateCount);
}

size_t modelGuidedProjectionRefinementLimit(size_t transitionSupportSize) {
  return largeSupportProjectionRefinementLimit(transitionSupportSize);
}

bool shouldMineLocalAuxiliaryInvariants(size_t transitionSupportSize) {
  return transitionSupportSize <= kLocalAuxiliaryInvariantSupportMiningLimit;
}

bool shouldUseDirectCdclCraigProjectionQuery(size_t trackedStateCount) {
  return trackedStateCount > kCraigDirectCdclProjectionQueryStateLimit;
}

bool shouldEncodeProjectedLocalStateSemantics(
    size_t leafCount,
    size_t trackedStateCount) {
  return leafCount <= kCraigProjectedLocalSemanticsSupportLimit &&
         trackedStateCount <= kCraigProjectedLocalSemanticsTrackedStateLimit;
}

bool shouldFocusImageTransitionRequests(
    size_t trackedStateCount,
    size_t lookahead) {
  return lookahead > 0 &&
         trackedStateCount > kCraigFocusedImageTransitionTrackedStateLimit;
}

bool shouldUseRetainedHelperFocusedImageRequestLimit(
    size_t trackedStateCount,
    size_t helperInvariantRegionCount) {
  return trackedStateCount >=
             kCraigRetainedHelperFocusedImageTransitionStateThreshold &&
         helperInvariantRegionCount > 0;
}

bool shouldUseBroadRetainedHelperFocusedImageRequestLimit(
    size_t trackedStateCount,
    size_t helperInvariantRegionCount) {
  return trackedStateCount >=
             kCraigBroadRetainedHelperFocusedImageTransitionStateThreshold &&
         helperInvariantRegionCount >=
             kCraigBroadRetainedHelperFocusedImageTransitionRegionThreshold;
}

size_t focusedImageTransitionRequestLimit(
    size_t trackedStateCount,
    size_t helperInvariantRegionCount) {
  // Retained helper tails already carry proof-derived pruning facts.  Keep
  // later broad retained-helper tails on a smaller strict over-approximation
  // instead of rebuilding the same saturated suffix cone.
  if (shouldUseBroadRetainedHelperFocusedImageRequestLimit(
          trackedStateCount, helperInvariantRegionCount)) {
    return kCraigBroadRetainedHelperFocusedImageTransitionRequestLimit;
  }
  if (shouldUseRetainedHelperFocusedImageRequestLimit(
          trackedStateCount, helperInvariantRegionCount)) {
    return kCraigRetainedHelperFocusedImageTransitionRequestLimit;
  }
  return kCraigFocusedImageTransitionRequestLimit;
}

bool shouldCapFocusedImageTransitionRequests(
    size_t expandedRequestCount,
    size_t requestLimit) {
  return expandedRequestCount > requestLimit;
}

size_t cappedFocusedImageTransitionRequestCount(
    size_t currentRequestCount,
    size_t expandedRequestCount,
    size_t requestLimit) {
  if (!shouldCapFocusedImageTransitionRequests(
          expandedRequestCount, requestLimit)) {
    return expandedRequestCount; // LCOV_EXCL_LINE
  }
  return std::max(currentRequestCount, requestLimit);
}

void configureCraigProjectionSolver(
    SATSolverWrapper& solver,
    size_t trackedStateCount,
    const char* phase) {
  if (!shouldUseDirectCdclCraigProjectionQuery(trackedStateCount)) {
    return;
  }
  // Reuse the shared direct-CaDiCaL SEC profile locally.  Craig tracing keeps
  // every preprocessing resolvent, so large BP projections must avoid variable
  // elimination instead of building an exponential proof trace before CDCL.
  solver.configureForSecPdrQuery(trackedStateCount);
  emitSecDiag(
      "SEC diag: imc Craig uses direct CaDiCaL profile phase=", phase,
      " tracked_states=", trackedStateCount,
      " state_limit=", kCraigDirectCdclProjectionQueryStateLimit);
}

void configureCraigProjectionSolverForTransitionRequests(
    SATSolverWrapper& solver,
    size_t transitionRequestCount,
    const char* phase) {
  if (!shouldUseDirectCdclCraigProjectionQuery(transitionRequestCount)) {
    return;
  }
  solver.configureForSecPdrQuery(transitionRequestCount);
  emitSecDiag(
      "SEC diag: imc Craig uses direct CaDiCaL profile phase=", phase,
      " transition_requests=", transitionRequestCount,
      " state_limit=", kCraigDirectCdclProjectionQueryStateLimit);
}

bool modelGuidedProjectionNeedsBoundedBackfill(
    const FrontierResult& frontier);

bool imcAuxiliaryInvariantsEnabled() {
  const char* enabled = std::getenv("KEPLER_SEC_IMC_AUX_INVARIANTS");
  return enabled != nullptr && std::strcmp(enabled, "1") == 0;
}

bool imcDirectCubeSourceEnabled(bool forceEnabled) {
  if (forceEnabled) {
    return true;
  }
  const char* enabled = std::getenv("KEPLER_SEC_IMC_DIRECT_CUBE_SOURCE");
  return enabled != nullptr && std::strcmp(enabled, "1") == 0;
}

int64_t elapsedMilliseconds(SteadyClock::time_point start) {
  return std::chrono::duration_cast<std::chrono::milliseconds>(
             SteadyClock::now() - start)
      .count();
}

struct TransitionEncodingResult {
  std::unordered_map<size_t, int> currentLits;
};

struct TransitionEncodingCache {
  std::unordered_map<size_t, bool> constantAssignments;
};

size_t regionClauseCount(const InterpolantRegion& region) {
  return region.definitionClauseEnds.size();
}

size_t regionLiteralCount(const InterpolantRegion& region) {
  return region.definitionLiterals.size();
}

bool sameRegionVariable(const RegionLiteral& lhs, const RegionLiteral& rhs) {
  return lhs.isState == rhs.isState && lhs.index == rhs.index;
}

bool sameRegionLiteral(const RegionLiteral& lhs, const RegionLiteral& rhs) {
  return sameRegionVariable(lhs, rhs) && lhs.positive == rhs.positive;
}

bool regionLiteralLess(const RegionLiteral& lhs, const RegionLiteral& rhs) {
  if (lhs.isState != rhs.isState) {
    return lhs.isState < rhs.isState;
  }
  if (lhs.index != rhs.index) {
    return lhs.index < rhs.index;
  }
  return lhs.positive < rhs.positive;
}

bool regionClauseSubsumes(const std::vector<RegionLiteral>& lhs,
                          const std::vector<RegionLiteral>& rhs) {
  return std::includes(
      rhs.begin(), rhs.end(), lhs.begin(), lhs.end(), regionLiteralLess);
}

bool regionClauseLess(const std::vector<RegionLiteral>& lhs,
                      const std::vector<RegionLiteral>& rhs) {
  return std::lexicographical_compare(
      lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), regionLiteralLess);
}

bool sameRegionClause(const std::vector<RegionLiteral>& lhs,
                      const std::vector<RegionLiteral>& rhs) {
  return lhs.size() == rhs.size() &&
         std::equal(lhs.begin(), lhs.end(), rhs.begin(), sameRegionLiteral);
}

int localRegionLiteral(
    SATSolverWrapper& solver,
    const RegionLiteral& literal,
    std::unordered_map<RegionVariableKey, int, RegionVariableKeyHash>& literals) {
  const RegionVariableKey key{literal.isState != 0, literal.index};
  auto [it, inserted] = literals.emplace(key, 0);
  if (inserted) {
    it->second = solver.newVar() + 2;
  }
  return literal.positive ? it->second : -it->second;
}

std::vector<std::vector<RegionLiteral>> normalizedRegionClauses(
    const InterpolantRegion& region) {
  std::vector<std::vector<RegionLiteral>> clauses;
  clauses.reserve(region.definitionClauseEnds.size());

  size_t clauseBegin = 0;
  for (const size_t clauseEnd : region.definitionClauseEnds) {
    std::vector<RegionLiteral> clause(
        region.definitionLiterals.begin() + clauseBegin,
        region.definitionLiterals.begin() + clauseEnd);
    clauseBegin = clauseEnd;
    std::sort(clause.begin(), clause.end(), regionLiteralLess);
    clause.erase(
        std::unique(clause.begin(), clause.end(), sameRegionLiteral),
        clause.end());

    bool tautology = false;
    for (size_t i = 1; i < clause.size(); ++i) {
      if (sameRegionVariable(clause[i - 1], clause[i]) &&
          clause[i - 1].positive != clause[i].positive) {
        tautology = true;
        break;
      }
    }
    if (!tautology) {
      clauses.push_back(std::move(clause));
    }
  }

  std::sort(clauses.begin(), clauses.end(), regionClauseLess);
  clauses.erase(
      std::unique(clauses.begin(), clauses.end(), sameRegionClause),
      clauses.end());
  return clauses;
}

void removeSubsumedRegionClauses(
    std::vector<std::vector<RegionLiteral>>& clauses) {
  if (clauses.size() > kCraigSubsumptionClauseLimit) {
    return; // LCOV_EXCL_LINE
  }

  std::vector<bool> removed(clauses.size(), false);
  for (size_t i = 0; i < clauses.size(); ++i) {
    if (removed[i]) {
      continue;
    }
    for (size_t j = 0; j < clauses.size(); ++j) {
      if (i == j || removed[j]) {
        continue;
      }
      if (regionClauseSubsumes(clauses[i], clauses[j])) {
        removed[j] = true;
      }
    }
  }

  size_t write = 0;
  for (size_t read = 0; read < clauses.size(); ++read) {
    if (!removed[read]) {
      if (write != read) {
        clauses[write] = std::move(clauses[read]); // LCOV_EXCL_LINE
      } // LCOV_EXCL_LINE
      ++write;
    }
  }
  clauses.resize(write);
}

size_t regionVariableCount(
    const std::vector<std::vector<RegionLiteral>>& clauses) {
  std::unordered_set<RegionVariableKey, RegionVariableKeyHash> variables;
  for (const auto& clause : clauses) {
    for (const RegionLiteral& literal : clause) {
      variables.insert({literal.isState != 0, literal.index});
    }
  }
  return variables.size();
}

void removeSemanticallyRedundantRegionClauses(
    std::vector<std::vector<RegionLiteral>>& clauses) {
  if (clauses.size() > kCraigSemanticSimplifyClauseLimit ||
      regionVariableCount(clauses) > kCraigSemanticSimplifyVariableLimit) {
    return; // LCOV_EXCL_LINE
  }

  std::vector<bool> removed(clauses.size(), false);
  for (size_t candidate = 0; candidate < clauses.size(); ++candidate) {
    if (removed[candidate]) {
      continue; // LCOV_EXCL_LINE
    }

    SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::KISSAT);
    std::unordered_map<RegionVariableKey, int, RegionVariableKeyHash> literals;
    for (size_t clauseIndex = 0; clauseIndex < clauses.size(); ++clauseIndex) {
      if (clauseIndex == candidate || removed[clauseIndex]) {
        continue;
      }
      std::vector<int> clause;
      clause.reserve(clauses[clauseIndex].size());
      for (const RegionLiteral& literal : clauses[clauseIndex]) {
        clause.push_back(localRegionLiteral(solver, literal, literals));
      }
      solver.addClause(clause);
    }
    for (const RegionLiteral& literal : clauses[candidate]) {
      solver.addClause({-localRegionLiteral(solver, literal, literals)});
    }
    if (solver.solveStatus() == SATSolverWrapper::SolveStatus::Unsat) {
      removed[candidate] = true;
    }
  }

  size_t write = 0;
  for (size_t read = 0; read < clauses.size(); ++read) {
    if (!removed[read]) {
      if (write != read) {
        clauses[write] = std::move(clauses[read]);
      }
      ++write;
    }
  }
  clauses.resize(write);
}

void rebuildRegionClauses(
    InterpolantRegion& region,
    const std::vector<std::vector<RegionLiteral>>& clauses) {
  region.definitionLiterals.clear();
  region.definitionClauseEnds.clear();
  size_t literalCount = 0;
  for (const auto& clause : clauses) {
    literalCount += clause.size();
  }
  region.definitionLiterals.reserve(literalCount);
  region.definitionClauseEnds.reserve(clauses.size());
  for (const auto& clause : clauses) {
    region.definitionLiterals.insert(
        region.definitionLiterals.end(), clause.begin(), clause.end());
    region.definitionClauseEnds.push_back(region.definitionLiterals.size());
  }
}

InterpolantRegion simplifyCraigInterpolantRegionImpl(InterpolantRegion region) {
  if (region.type != InterpolantRegion::Type::Normal) {
    return region; // LCOV_EXCL_LINE
  }

  const size_t oldClauses = regionClauseCount(region);
  const size_t oldLiterals = regionLiteralCount(region);
  auto clauses = normalizedRegionClauses(region);
  if (std::any_of(clauses.begin(), clauses.end(),
                  [](const auto& clause) { return clause.empty(); })) {
    return {InterpolantRegion::Type::False}; // LCOV_EXCL_LINE
  }

  // McMillan's original IMC implementation reduced redundant interpolant logic
  // with small BDDs. Keep the same role local to Craig IMC: first remove cheap
  // syntactic redundancy, then run a bounded exact SAT implication cleanup for
  // small interpolants where the extra solver calls are predictable.
  removeSubsumedRegionClauses(clauses);
  removeSemanticallyRedundantRegionClauses(clauses);
  if (clauses.empty() && !region.root.isState) {
    return {InterpolantRegion::Type::True}; // LCOV_EXCL_LINE
  }
  rebuildRegionClauses(region, clauses);

  emitSecDiag(
      "SEC diag: imc Craig interpolant simplified clauses=",
      oldClauses,
      "->",
      regionClauseCount(region),
      " literals=",
      oldLiterals,
      "->",
      regionLiteralCount(region));
  return region;
}

int instantiateRegionLiteral(
    const RegionLiteral& literal,
    const std::unordered_map<size_t, int>& stateLits,
    const std::vector<int>& auxiliaryLits) {
  const int positive = literal.isState
                           ? stateLits.at(literal.index)
                           : auxiliaryLits.at(literal.index);
  return literal.positive ? positive : -positive;
}

std::unordered_set<size_t> stateSymbolSet(
    const KInductionProblem& problem) {
  std::unordered_set<size_t> states;
  states.reserve(problem.state0Symbols.size() + problem.state1Symbols.size());
  states.insert(problem.state0Symbols.begin(), problem.state0Symbols.end());
  states.insert(problem.state1Symbols.begin(), problem.state1Symbols.end());
  return states;
}

void closePairDependencies(
    const std::vector<std::pair<size_t, size_t>>& pairs,
    std::unordered_set<size_t>& symbols) {
  bool changed = true;
  while (changed) {
    changed = false;
    for (const auto& [lhs, rhs] : pairs) {
      if (!symbols.contains(lhs) && !symbols.contains(rhs)) {
        continue;
      }
      changed |= symbols.insert(lhs).second;
      changed |= symbols.insert(rhs).second;
    }
  }
}

void closeSameDesignStateSemantics(
    const KInductionProblem& problem,
    std::unordered_set<size_t>& symbols) {
  closePairDependencies(problem.complementedStatePairs0, symbols);
  closePairDependencies(problem.complementedStatePairs1, symbols);
  closePairDependencies(problem.sameFrameStateEqualityPairs0, symbols);
  closePairDependencies(problem.sameFrameStateEqualityPairs1, symbols);
  for (const auto& rails : problem.dualRailStatePairs) {
    if (symbols.contains(rails.mayBeOne) ||
        symbols.contains(rails.mayBeZero)) {
      symbols.insert(rails.mayBeOne);
      symbols.insert(rails.mayBeZero);
    }
  }
}

std::vector<size_t> sortedSymbols(
    const std::unordered_set<size_t>& symbols) {
  std::vector<size_t> sorted(symbols.begin(), symbols.end());
  std::sort(sorted.begin(), sorted.end());
  return sorted;
}

void closeProjectionRefinementTransitionTargets(
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    std::unordered_set<size_t>& symbols) {
  // Refinement does not need to track every same-design semantic partner.  If a
  // selected complemented rail maps to a primary transition target, add that
  // target so the next-state relation is still encoded; leaving other partners
  // untracked only over-approximates the Craig query, which is sound for proofs.
  for (const size_t requested : sortedSymbols(symbols)) {
    const size_t target =
        transitionTargetFor(requested, resolver, complementPrimary);
    if (resolver.contains(target)) {
      symbols.insert(target);
    }
  }
}

void insertFocusedTransitionRequestWithTarget( // LCOV_EXCL_LINE
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& trackedStates,
    size_t symbol,
    std::unordered_set<size_t>& requests) {
  if (!trackedStates.contains(symbol)) { // LCOV_EXCL_LINE
    return; // LCOV_EXCL_LINE
  }
  const size_t target = // LCOV_EXCL_LINE
      transitionTargetFor(symbol, resolver, complementPrimary); // LCOV_EXCL_LINE
  if (resolver.contains(target) && trackedStates.contains(target)) { // LCOV_EXCL_LINE
    requests.insert(target); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  requests.insert(symbol); // LCOV_EXCL_LINE
} // LCOV_EXCL_LINE

std::unordered_set<size_t> selectCappedFocusedImageTransitionRequests( // LCOV_EXCL_LINE
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& trackedStates,
    const std::unordered_set<size_t>& currentRequests,
    const std::unordered_set<size_t>& expandedRequests,
    size_t requestLimit) {
  std::unordered_set<size_t> requests = currentRequests; // LCOV_EXCL_LINE
  requests.reserve(std::min(expandedRequests.size(), requestLimit)); // LCOV_EXCL_LINE
  if (requests.size() >= requestLimit) { // LCOV_EXCL_LINE
    return requests; // LCOV_EXCL_LINE
  }
  for (const size_t symbol : sortedSymbols(expandedRequests)) { // LCOV_EXCL_LINE
    if (requests.size() >= requestLimit) { // LCOV_EXCL_LINE
      break; // LCOV_EXCL_LINE
    }
    if (requests.contains(symbol)) { // LCOV_EXCL_LINE
      continue; // LCOV_EXCL_LINE
    }
    insertFocusedTransitionRequestWithTarget( // LCOV_EXCL_LINE
        resolver, complementPrimary, trackedStates, symbol, requests); // LCOV_EXCL_LINE
  }
  return requests; // LCOV_EXCL_LINE
} // LCOV_EXCL_LINE

std::unordered_set<size_t> initialTrackedStates(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& states) {
  std::unordered_set<size_t> tracked;
  for (const size_t symbol : problem.bad->getSupportVars()) {
    if (states.contains(symbol)) {
      tracked.insert(symbol);
    }
  }
  closeSameDesignStateSemantics(problem, tracked);
  return tracked;
}

std::unordered_set<size_t> focusedImageTransitionRequests(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& states,
    const std::unordered_set<size_t>& trackedStates,
    size_t helperInvariantRegionCount,
    size_t lookahead) {
  if (!shouldFocusImageTransitionRequests(trackedStates.size(), lookahead)) {
    return trackedStates;
  }
  const size_t requestLimit = focusedImageTransitionRequestLimit(
      trackedStates.size(), helperInvariantRegionCount);

  std::unordered_set<size_t> requests;
  for (const size_t symbol : problem.bad->getSupportVars()) {
    if (states.contains(symbol) && trackedStates.contains(symbol)) {
      requests.insert(symbol);
    }
  }
  closeSameDesignStateSemantics(problem, requests);
  closeProjectionRefinementTransitionTargets(
      resolver, complementPrimary, requests);
  for (size_t depth = 1; depth < lookahead; ++depth) {
    // A multi-step image needs the preimage support of the suffix-observed bad
    // slice.  Grow the focused request set one transition layer per remaining
    // suffix step; unrelated projected globals still stay unconstrained.
    std::unordered_set<size_t> expandedRequests = requests; // LCOV_EXCL_LINE
    for (const size_t requested : sortedSymbols(requests)) { // LCOV_EXCL_LINE
      const size_t target = // LCOV_EXCL_LINE
          transitionTargetFor(requested, resolver, complementPrimary); // LCOV_EXCL_LINE
      if (!resolver.contains(target)) { // LCOV_EXCL_LINE
        continue; // LCOV_EXCL_LINE
      }
      for (const size_t symbol : resolver.support(target)) { // LCOV_EXCL_LINE
        if (states.contains(symbol) && trackedStates.contains(symbol)) { // LCOV_EXCL_LINE
          expandedRequests.insert(symbol); // LCOV_EXCL_LINE
        } // LCOV_EXCL_LINE
      }
    }
    closeSameDesignStateSemantics(problem, expandedRequests); // LCOV_EXCL_LINE
    closeProjectionRefinementTransitionTargets( // LCOV_EXCL_LINE
        resolver, complementPrimary, expandedRequests); // LCOV_EXCL_LINE
    const size_t cappedRequestCount = cappedFocusedImageTransitionRequestCount( // LCOV_EXCL_LINE
        requests.size(), expandedRequests.size(), requestLimit); // LCOV_EXCL_LINE
    if (cappedRequestCount < expandedRequests.size()) { // LCOV_EXCL_LINE
      // Keeping only a deterministic prefix of the expanded request layer
      // leaves the remaining projected next-state functions unconstrained.
      // That weakens the image query, preserving strict Craig proof soundness
      // while avoiding BP's 130K-leaf lookahead-4 proof spike.
      requests = selectCappedFocusedImageTransitionRequests( // LCOV_EXCL_LINE
          resolver, // LCOV_EXCL_LINE
          complementPrimary, // LCOV_EXCL_LINE
          trackedStates, // LCOV_EXCL_LINE
          requests,
          expandedRequests,
          cappedRequestCount); // LCOV_EXCL_LINE
      emitSecDiag( // LCOV_EXCL_LINE
          "SEC diag: imc Craig caps focused image transition requests "
          "tracked_states=",
          trackedStates.size(), // LCOV_EXCL_LINE
          " lookahead=", lookahead,
          " depth=", depth + 1, // LCOV_EXCL_LINE
          " requests=", requests.size(), // LCOV_EXCL_LINE
          " expanded=", expandedRequests.size(), // LCOV_EXCL_LINE
          " limited=", cappedRequestCount,
          " request_limit=", requestLimit);
      break; // LCOV_EXCL_LINE
    }
    requests = std::move(expandedRequests); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  for (auto it = requests.begin(); it != requests.end();) {
    if (!trackedStates.contains(*it)) {
      it = requests.erase(it); // LCOV_EXCL_LINE
    } else { // LCOV_EXCL_LINE
      ++it;
    }
  }
  if (requests.empty()) {
    return trackedStates;
  }

  // A safety-image query only observes the bad predicate on the suffix frames.
  // Encoding transitions for unrelated tracked globals makes a stronger image
  // but is not required for a strict Craig proof; the smaller request set is a
  // sound over-approximation and keeps BP-sized proof traces from materializing
  // every projected next-state cone.
  emitSecDiag(
      "SEC diag: imc Craig focuses image transition requests tracked_states=",
      trackedStates.size(),
      " requests=", requests.size(),
      " lookahead=", lookahead,
      " tracked_limit=", kCraigFocusedImageTransitionTrackedStateLimit);
  return requests;
}

std::vector<size_t> projectedTransitionTargets(
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& trackedStates,
    const std::unordered_set<size_t>& transitionRequests) {
  std::vector<size_t> targets;
  targets.reserve(transitionRequests.size());
  std::unordered_set<size_t> encodedTargets;
  encodedTargets.reserve(transitionRequests.size());
  for (const size_t requested : sortedSymbols(transitionRequests)) {
    if (!trackedStates.contains(requested)) {
      continue; // LCOV_EXCL_LINE
    }
    const size_t target =
        transitionTargetFor(requested, resolver, complementPrimary);
    if (!trackedStates.contains(target) || !resolver.contains(target) ||
        !encodedTargets.insert(target).second) {
      continue;
    }
    targets.push_back(target);
  }
  return targets;
}

std::unordered_map<size_t, int> allocateLeafLits(
    SATSolverWrapper& solver,
    const std::vector<size_t>& symbols,
    VariablePartition partition) {
  solver.setCraigVariablePartition(partition);
  std::unordered_map<size_t, int> leaves;
  leaves.reserve(symbols.size());
  for (const size_t symbol : symbols) {
    leaves.emplace(symbol, solver.newVar() + 2);
  }
  return leaves;
}

class ConstantAwareFrameFormulaEncoder {
 public:
  ConstantAwareFrameFormulaEncoder(
      SATSolverWrapper& solver,
      std::unordered_map<size_t, int> leafLits,
      const std::unordered_map<size_t, bool>& constantAssignments,
      VariablePartition variablePartition,
      ClausePartition clausePartition,
      size_t expectedNodeHint = 0)
      : solver_(solver),
        leafLits_(std::move(leafLits)),
        constantAssignments_(constantAssignments),
        variablePartition_(variablePartition),
        clausePartition_(clausePartition),
        nodeArena_(nodeArenaBuffer_.data(), nodeArenaBuffer_.size()),
        nodeToLit_(&nodeArena_) {
    size_t expectedNodes = std::max(
        static_cast<size_t>(256),
        leafLits_.size() * static_cast<size_t>(8));
    if (expectedNodeHint != 0) {
      const size_t hintedNodes =
          expectedNodeHint +
          std::max(expectedNodeHint / 8, static_cast<size_t>(256));
      expectedNodes = std::max(expectedNodes, hintedNodes);
      leafLits_.reserve(
          leafLits_.size() +
          std::min(
              std::max(expectedNodeHint / 16, static_cast<size_t>(256)),
              static_cast<size_t>(65536)));
    }
    nodeToLit_.reserve(expectedNodes);
    constantValueByNode_.reserve(expectedNodes);
    solver_.reserveAdditionalVars(
        std::min(expectedNodes, kCraigMaxSolverTseitinReserveHint));
  }

  int encode(BoolExpr* expr) {
    if (expr == nullptr) {
      throw std::invalid_argument( // LCOV_EXCL_LINE
          "ConstantAwareFrameFormulaEncoder::encode: null expr");
    }
    return encodeImpl(expr);
  } // LCOV_EXCL_LINE

  const std::unordered_map<size_t, int>& leafLits() const {
    return leafLits_;
  }

 private:
  int newLiteral() {
    solver_.setCraigVariablePartition(variablePartition_);
    return solver_.newVar() + 2;
  }

  void addClause(const std::vector<int>& clause) {
    solver_.setCraigClausePartition(clausePartition_);
    solver_.addClause(clause);
  }

  int constLit(bool value) {
    if (trueLit_ == 0) {
      trueLit_ = newLiteral();
      addClause({trueLit_});
    }
    return value ? trueLit_ : -trueLit_;
  } // LCOV_EXCL_LINE

  bool isConstLit(int lit, bool value) const {
    return trueLit_ != 0 && lit == (value ? trueLit_ : -trueLit_);
  }

  struct BoolExprStackFrame {
    BoolExpr* expr = nullptr;
    bool visited = false;
  };

  std::optional<bool> cachedConstantValue(BoolExpr* expr) const {
    const auto cached = constantValueByNode_.find(expr);
    if (cached == constantValueByNode_.end()) {
      return std::nullopt;
    }
    return cached->second;
  }

  std::optional<bool> assignedConstantValue(BoolExpr* expr) {
    if (const auto cached = constantValueByNode_.find(expr);
        cached != constantValueByNode_.end()) {
      return cached->second;
    }

    std::vector<BoolExprStackFrame> stack;
    stack.reserve(512);
    stack.push_back({expr, false});

    while (!stack.empty()) {
      const BoolExprStackFrame current = stack.back();
      stack.pop_back();
      BoolExpr* node = current.expr;
      if (constantValueByNode_.find(node) != constantValueByNode_.end()) {
        continue;
      }

      if (node->getOp() == Op::VAR) {
        std::optional<bool> value;
        const size_t symbol = node->getId();
        if (symbol == 0) {
          value = false;
        } else if (symbol == 1) {
          value = true; // LCOV_EXCL_LINE
        } else if (const auto constant = constantAssignments_.find(symbol);
                   constant != constantAssignments_.end()) {
          value = constant->second;
        }
        constantValueByNode_.emplace(node, value);
        continue;
      }

      if (!current.visited) {
        stack.push_back({node, true});
        if (node->getRight() != nullptr &&
            constantValueByNode_.find(node->getRight()) ==
                constantValueByNode_.end()) {
          stack.push_back({node->getRight(), false});
        }
        if (node->getLeft() != nullptr &&
            constantValueByNode_.find(node->getLeft()) ==
                constantValueByNode_.end()) {
          stack.push_back({node->getLeft(), false});
        }
        continue;
      }

      const auto left =
          node->getLeft() == nullptr ? std::optional<bool>{}
                                     : cachedConstantValue(node->getLeft());
      const auto right =
          node->getRight() == nullptr ? std::optional<bool>{}
                                      : cachedConstantValue(node->getRight());
      std::optional<bool> value;
      switch (node->getOp()) {
        case Op::NOT:
          if (left.has_value()) {
            value = !*left;
          }
          break;
        case Op::AND:
          if ((left.has_value() && !*left) ||
              (right.has_value() && !*right)) {
            value = false;
          } else if (left.has_value() && right.has_value()) {
            value = *left && *right; // LCOV_EXCL_LINE
          } else if (left.has_value() && *left) {
            value = right;
          } else if (right.has_value() && *right) {
            value = left; // LCOV_EXCL_LINE
          } // LCOV_EXCL_LINE
          break;
        case Op::OR:
          if ((left.has_value() && *left) ||
              (right.has_value() && *right)) {
            value = true; // LCOV_EXCL_LINE
          } else if (left.has_value() && right.has_value()) {
            value = *left || *right;
          } else if (left.has_value() && !*left) {
            value = right;
          } else if (right.has_value() && !*right) {
            value = left;
          }
          break;
        case Op::XOR:
          if (left.has_value() && right.has_value()) { // LCOV_EXCL_LINE
            value = *left != *right; // LCOV_EXCL_LINE
          } // LCOV_EXCL_LINE
          break; // LCOV_EXCL_LINE
        case Op::VAR: // LCOV_EXCL_LINE
        case Op::NONE: // LCOV_EXCL_LINE
        default:
          throw std::runtime_error( // LCOV_EXCL_LINE
              "Unsupported BoolExpr operator in constant-aware evaluation");
      }

      constantValueByNode_.emplace(node, value);
    }

    return cachedConstantValue(expr);
  }

  int encodeImpl(BoolExpr* expr) {
    std::vector<BoolExprStackFrame> stack;
    stack.reserve(512);
    stack.push_back({expr, false});

    while (!stack.empty()) {
      const BoolExprStackFrame current = stack.back();
      stack.pop_back();
      BoolExpr* node = current.expr;
      if (nodeToLit_.find(node) != nodeToLit_.end()) {
        continue;
      }

      if (const auto value = assignedConstantValue(node); value.has_value()) {
        nodeToLit_.emplace(node, constLit(*value));
        continue;
      }

      if (node->getOp() == Op::VAR) {
        const size_t symbol = node->getId();
        if (symbol == 0) {
          nodeToLit_.emplace(node, constLit(false)); // LCOV_EXCL_LINE
          continue; // LCOV_EXCL_LINE
        }
        if (symbol == 1) {
          nodeToLit_.emplace(node, constLit(true)); // LCOV_EXCL_LINE
          continue; // LCOV_EXCL_LINE
        }
        if (const auto constant = constantAssignments_.find(symbol);
            constant != constantAssignments_.end()) {
          nodeToLit_.emplace(node, constLit(constant->second)); // LCOV_EXCL_LINE
          continue; // LCOV_EXCL_LINE
        }
        auto leaf = leafLits_.find(symbol);
        if (leaf == leafLits_.end()) {
          leaf = leafLits_.emplace(symbol, newLiteral()).first;
        }
        nodeToLit_.emplace(node, leaf->second);
        continue;
      }

      if (!current.visited) {
        stack.push_back({node, true});
        if (node->getRight() != nullptr &&
            nodeToLit_.find(node->getRight()) == nodeToLit_.end()) {
          stack.push_back({node->getRight(), false});
        }
        if (node->getLeft() != nullptr &&
            nodeToLit_.find(node->getLeft()) == nodeToLit_.end()) {
          stack.push_back({node->getLeft(), false});
        }
        continue;
      }

      const auto leftIt =
          node->getLeft() == nullptr ? nodeToLit_.end()
                                     : nodeToLit_.find(node->getLeft());
      const auto rightIt =
          node->getRight() == nullptr ? nodeToLit_.end()
                                      : nodeToLit_.find(node->getRight());
      const int leftLit =
          leftIt == nodeToLit_.end() ? 0 : leftIt->second;
      const int rightLit =
          rightIt == nodeToLit_.end() ? 0 : rightIt->second;
      int lit = 0;
      switch (node->getOp()) {
        case Op::NOT:
          lit = -leftLit;
          break;
        case Op::AND:
          if (isConstLit(leftLit, false)) {
            lit = constLit(false); // LCOV_EXCL_LINE
          } else if (isConstLit(leftLit, true)) {
            lit = rightLit;
          } else if (leftLit == rightLit || isConstLit(rightLit, true)) {
            lit = leftLit; // LCOV_EXCL_LINE
          } else if (leftLit == -rightLit || isConstLit(rightLit, false)) {
            lit = constLit(false); // LCOV_EXCL_LINE
          } else { // LCOV_EXCL_LINE
            lit = newLiteral();
            addClause({-lit, leftLit});
            addClause({-lit, rightLit});
            addClause({lit, -leftLit, -rightLit});
          }
          break;
        case Op::OR:
          if (isConstLit(leftLit, true)) {
            lit = constLit(true); // LCOV_EXCL_LINE
          } else if (isConstLit(leftLit, false)) {
            lit = rightLit;
          } else if (leftLit == rightLit || isConstLit(rightLit, false)) {
            lit = leftLit;
          } else if (leftLit == -rightLit || isConstLit(rightLit, true)) {
            lit = constLit(true); // LCOV_EXCL_LINE
          } else { // LCOV_EXCL_LINE
            lit = newLiteral(); // LCOV_EXCL_LINE
            addClause({-leftLit, lit}); // LCOV_EXCL_LINE
            addClause({-rightLit, lit}); // LCOV_EXCL_LINE
            addClause({-lit, leftLit, rightLit}); // LCOV_EXCL_LINE
          }
          break;
        case Op::XOR:
          if (isConstLit(leftLit, false)) { // LCOV_EXCL_LINE
            lit = rightLit; // LCOV_EXCL_LINE
          } else if (isConstLit(leftLit, true)) { // LCOV_EXCL_LINE
            lit = -rightLit; // LCOV_EXCL_LINE
          } else if (leftLit == rightLit) { // LCOV_EXCL_LINE
            lit = constLit(false); // LCOV_EXCL_LINE
          } else if (leftLit == -rightLit) { // LCOV_EXCL_LINE
            lit = constLit(true); // LCOV_EXCL_LINE
          } else if (isConstLit(rightLit, false)) { // LCOV_EXCL_LINE
            lit = leftLit; // LCOV_EXCL_LINE
          } else if (isConstLit(rightLit, true)) { // LCOV_EXCL_LINE
            lit = -leftLit; // LCOV_EXCL_LINE
          } else { // LCOV_EXCL_LINE
            lit = newLiteral(); // LCOV_EXCL_LINE
            addClause({-lit, -leftLit, -rightLit}); // LCOV_EXCL_LINE
            addClause({-lit, leftLit, rightLit}); // LCOV_EXCL_LINE
            addClause({lit, -leftLit, rightLit}); // LCOV_EXCL_LINE
            addClause({lit, leftLit, -rightLit}); // LCOV_EXCL_LINE
          }
          break; // LCOV_EXCL_LINE
        case Op::VAR: // LCOV_EXCL_LINE
        case Op::NONE: // LCOV_EXCL_LINE
        default:
          throw std::runtime_error( // LCOV_EXCL_LINE
              "Unsupported BoolExpr operator in constant-aware encoding");
      }

      nodeToLit_.emplace(node, lit);
    }

    return nodeToLit_.at(expr);
  }

  SATSolverWrapper& solver_;
  std::unordered_map<size_t, int> leafLits_;
  const std::unordered_map<size_t, bool>& constantAssignments_;
  VariablePartition variablePartition_;
  ClausePartition clausePartition_;
  std::array<std::byte, 16 * 1024> nodeArenaBuffer_{};
  std::pmr::monotonic_buffer_resource nodeArena_;
  std::pmr::unordered_map<BoolExpr*, int> nodeToLit_;
  std::unordered_map<BoolExpr*, std::optional<bool>> constantValueByNode_;
  int trueLit_ = 0;
};

void addLiteralEquivalenceForPartition(
    SATSolverWrapper& solver,
    int lhs,
    int rhs,
    ClausePartition partition) {
  solver.setCraigClausePartition(partition);
  addLiteralEquivalence(solver, lhs, rhs);
}

void addPairEqualities(
    SATSolverWrapper& solver,
    const std::unordered_map<size_t, int>& leaves,
    const std::vector<std::pair<size_t, size_t>>& pairs,
    bool complemented,
    ClausePartition partition) {
  for (const auto& [lhsSymbol, rhsSymbol] : pairs) {
    const auto lhs = leaves.find(lhsSymbol);
    const auto rhs = leaves.find(rhsSymbol);
    if (lhs == leaves.end() || rhs == leaves.end()) {
      continue;
    }
    addLiteralEquivalenceForPartition(
        solver,
        rhs->second,
        complemented ? -lhs->second : lhs->second,
        partition);
  }
}

void indexPairEqualities(
    StateSemanticsIndex& index,
    const std::vector<std::pair<size_t, size_t>>& pairs,
    bool complemented) {
  for (const auto& [lhsSymbol, rhsSymbol] : pairs) {
    index.pairSemanticsByLhs[lhsSymbol].push_back(
        {rhsSymbol, complemented});
  }
}

StateSemanticsIndex buildStateSemanticsIndex(
    const KInductionProblem& problem) {
  StateSemanticsIndex index;
  indexPairEqualities(
      index, problem.complementedStatePairs0, /*complemented=*/true);
  indexPairEqualities(
      index, problem.complementedStatePairs1, /*complemented=*/true);
  indexPairEqualities(
      index, problem.sameFrameStateEqualityPairs0, /*complemented=*/false);
  indexPairEqualities(
      index, problem.sameFrameStateEqualityPairs1, /*complemented=*/false);
  for (const auto& rails : problem.dualRailStatePairs) {
    index.dualRailSemanticsByMayOne[rails.mayBeOne].push_back(
        {rails.mayBeZero});
  }
  return index;
}

void addIndexedStateSemantics(
    SATSolverWrapper& solver,
    const StateSemanticsIndex& index,
    const std::unordered_map<size_t, int>& leaves,
    ClausePartition partition) {
  for (const auto& [lhsSymbol, lhsLit] : leaves) {
    const auto relations = index.pairSemanticsByLhs.find(lhsSymbol);
    if (relations == index.pairSemanticsByLhs.end()) {
      continue;
    }
    for (const IndexedStatePairSemantics& relation : relations->second) {
      const auto rhs = leaves.find(relation.rhs);
      if (rhs == leaves.end()) {
        continue;
      }
      addLiteralEquivalenceForPartition(
          solver,
          rhs->second,
          relation.complemented ? -lhsLit : lhsLit,
          partition);
    }
  }

  solver.setCraigClausePartition(partition);
  for (const auto& [mayOneSymbol, mayOneLit] : leaves) {
    const auto relations = index.dualRailSemanticsByMayOne.find(mayOneSymbol);
    if (relations == index.dualRailSemanticsByMayOne.end()) {
      continue;
    }
    for (const IndexedDualRailSemantics& relation : relations->second) {
      const auto mayZero = leaves.find(relation.mayBeZero);
      if (mayZero != leaves.end()) {
        solver.addClause({mayOneLit, mayZero->second});
      }
    }
  }
}

void addStateSemantics(
    SATSolverWrapper& solver,
    const KInductionProblem& problem,
    const std::unordered_map<size_t, int>& leaves,
    ClausePartition partition) {
  addPairEqualities(
      solver, leaves, problem.complementedStatePairs0, true, partition);
  addPairEqualities(
      solver, leaves, problem.complementedStatePairs1, true, partition);
  addPairEqualities(
      solver, leaves, problem.sameFrameStateEqualityPairs0, false, partition);
  addPairEqualities(
      solver, leaves, problem.sameFrameStateEqualityPairs1, false, partition);
  solver.setCraigClausePartition(partition);
  for (const auto& rails : problem.dualRailStatePairs) {
    const auto mayOne = leaves.find(rails.mayBeOne);
    const auto mayZero = leaves.find(rails.mayBeZero);
    if (mayOne != leaves.end() && mayZero != leaves.end()) {
      solver.addClause({mayOne->second, mayZero->second});
    }
  }
}

void addAuxiliaryStateInvariants(
    SATSolverWrapper& solver,
    const std::unordered_map<size_t, int>& leaves,
    const AuxiliaryStateInvariants& invariants,
    ClausePartition partition) {
  if (invariants.empty()) {
    return;
  }
  solver.setCraigClausePartition(partition);
  for (const auto& [symbol, value] : invariants.constants) {
    if (const auto leaf = leaves.find(symbol); leaf != leaves.end()) {
      solver.addClause({value ? leaf->second : -leaf->second});
    }
  }
  for (const auto& [lhsSymbol, rhsSymbol] : invariants.equalities) {
    const auto lhs = leaves.find(lhsSymbol);
    const auto rhs = leaves.find(rhsSymbol);
    if (lhs == leaves.end() || rhs == leaves.end()) {
      continue;
    }
    addLiteralEquivalenceForPartition(
        solver, lhs->second, rhs->second, partition);
  }
}

std::unordered_map<size_t, bool> auxiliaryConstantAssignments(
    const AuxiliaryStateInvariants& invariants) {
  std::unordered_map<size_t, bool> assignments;
  assignments.reserve(invariants.constants.size());
  for (const auto& [symbol, value] : invariants.constants) {
    assignments.emplace(symbol, value);
  }
  return assignments;
}

std::unordered_map<size_t, bool> transitionConstantAssignmentsForEncoding(
    const KInductionProblem& problem,
    const AuxiliaryStateInvariants& invariants) {
  std::unordered_map<size_t, bool> assignments =
      auxiliaryConstantAssignments(invariants);
  assignments.reserve(assignments.size() + problem.resetBootstrapInputs.size());
  for (const auto& [symbol, assertedValue] : problem.resetBootstrapInputs) {
    // Projected transition queries already constrain reset inputs to their
    // inactive value after the bootstrap frontier. Substitute the same value
    // before Tseitin encoding so reset-gated BP cones never enter the Craig
    // proof trace.
    assignments[symbol] = !assertedValue;
  }
  return assignments;
}

void resetTransitionEncodingCache(
    TransitionEncodingCache& cache,
    const KInductionProblem& problem,
    const AuxiliaryStateInvariants& invariants) {
  cache.constantAssignments =
      transitionConstantAssignmentsForEncoding(problem, invariants);
}

std::unordered_map<size_t, bool> bootstrapStateConstants(
    const KInductionProblem& problem) {
  std::unordered_map<size_t, bool> constants;
  constants.reserve(problem.bootstrapStateAssignments.size());
  const auto states = stateSymbolSet(problem);
  for (const auto& [symbol, value] : problem.bootstrapStateAssignments) {
    if (states.contains(symbol)) {
      constants[symbol] = value;
    }
  }
  return constants;
}

void closeComplementConstantsForPairs(
    const std::vector<std::pair<size_t, size_t>>& pairs,
    std::unordered_map<size_t, bool>& constants) {
  for (const auto& [primary, complement] : pairs) {
    const auto primaryValue = constants.find(primary);
    const auto complementValue = constants.find(complement);
    if (primaryValue != constants.end() &&
        complementValue == constants.end()) {
      constants.emplace(complement, !primaryValue->second); // LCOV_EXCL_LINE
    } else if (primaryValue == constants.end() &&
               complementValue != constants.end()) {
      constants.emplace(primary, !complementValue->second); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
  }
}

void closeComplementStateConstants(
    const KInductionProblem& problem,
    std::unordered_map<size_t, bool>& constants) {
  // Q/QN rail pairs are local to one design and are already encoded as state
  // semantics in every Craig query.  Closing auxiliary constants over those
  // pairs lets a proof of Q=0 also prune support that only mentions QN=1.
  closeComplementConstantsForPairs(problem.complementedStatePairs0, constants);
  closeComplementConstantsForPairs(problem.complementedStatePairs1, constants);
}

void closeAuxiliaryInvariantConstantComplements(
    const KInductionProblem& problem,
    AuxiliaryStateInvariants& invariants) {
  std::unordered_map<size_t, bool> constants;
  constants.reserve(invariants.constants.size());
  for (const auto& [symbol, value] : invariants.constants) {
    constants.emplace(symbol, value);
  }
  closeComplementStateConstants(problem, constants);
  invariants.constants.assign(constants.begin(), constants.end());
  std::sort(invariants.constants.begin(), invariants.constants.end());
}

std::unordered_map<size_t, bool> bootstrapAssignmentsForSymbols(
    const std::vector<std::pair<size_t, bool>>& assignments,
    const std::unordered_set<size_t>& symbols) {
  std::unordered_map<size_t, bool> constants;
  if (symbols.empty()) {
    return constants; // LCOV_EXCL_LINE
  }
  constants.reserve(std::min(assignments.size(), symbols.size()));
  for (const auto& [symbol, value] : assignments) {
    if (symbols.contains(symbol)) {
      constants[symbol] = value;
    }
  }
  return constants;
}

std::optional<bool> constantBoolExprValue(BoolExpr* expr) {
  if (expr == nullptr || expr->getOp() != Op::VAR || expr->getId() > 1) {
    return std::nullopt;
  }
  return expr->getId() == 1;
}

std::optional<bool> evaluateBoolExprWithKnownSupport(
    BoolExpr* expr,
    const std::set<size_t>& support,
    const std::unordered_map<size_t, bool>& constants) {
  if (expr == nullptr) {
    return std::nullopt; // LCOV_EXCL_LINE
  }
  for (const size_t supportSymbol : support) {
    if (supportSymbol >= 2 && !constants.contains(supportSymbol)) {
      return std::nullopt;
    }
  }
  try {
    return expr->evaluate(constants);
  } catch (const std::exception&) { // LCOV_EXCL_LINE
    return std::nullopt; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
}

std::optional<std::pair<size_t, bool>> transitionConstantTargetFor(
    size_t symbol,
    bool value,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary) {
  if (resolver.contains(symbol)) {
    return std::make_pair(symbol, value);
  }
  const auto primary = complementPrimary.find(symbol);
  if (primary != complementPrimary.end() && resolver.contains(primary->second)) {
    // Complement rails may not have an explicit transition expression.  Because
    // same-design state semantics enforces complement = !primary, proving the
    // primary next-state has the opposite value proves the complement constant.
    return std::make_pair(primary->second, !value);
  }
  return std::nullopt;
}

size_t localAuxiliaryPreservationScore(
    size_t symbol,
    bool value,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_map<size_t, bool>& constants) {
  const auto target =
      transitionConstantTargetFor(symbol, value, resolver, complementPrimary);
  if (!target.has_value()) {
    return 0;
  }
  BoolExpr* expr = resolver.at(target->first);
  if (expr == nullptr) {
    return 0; // LCOV_EXCL_LINE
  }
  if (const auto constantValue = constantBoolExprValue(expr);
      constantValue.has_value()) {
    return *constantValue == target->second
               ? kLocalAuxiliaryPreservationScoreBonus
               : 0;
  }
  const auto& support = resolver.support(target->first); // LCOV_EXCL_LINE
  if (support.size() > kLocalAuxiliaryPreservationScoreSupportLimit) { // LCOV_EXCL_LINE
    return 0; // LCOV_EXCL_LINE
  }
  if (const auto knownValue = // LCOV_EXCL_LINE
          evaluateBoolExprWithKnownSupport(expr, support, constants); // LCOV_EXCL_LINE
      knownValue.has_value() && *knownValue == target->second) { // LCOV_EXCL_LINE
    return kLocalAuxiliaryPreservationScoreBonus; // LCOV_EXCL_LINE
  }
  return 0; // LCOV_EXCL_LINE
}

bool transitionPreservesStateConstant(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const StateSemanticsIndex& stateSemanticsIndex,
    size_t symbol,
    bool value,
    const std::unordered_map<size_t, bool>& constants,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    StateConstantPreservationCache* preservationCache) {
  const auto target =
      transitionConstantTargetFor(symbol, value, resolver, complementPrimary);
  if (!target.has_value()) {
    return false;
  }
  const size_t targetSymbol = target->first;
  const bool targetValue = target->second;
  if (resolver.at(targetSymbol) == nullptr) {
    return false; // LCOV_EXCL_LINE
  }
  BoolExpr* targetExpr = resolver.at(targetSymbol);
  const StateConstantPreservationKey cacheKey{targetExpr, targetValue};
  if (preservationCache != nullptr) {
    const auto cached = preservationCache->find(cacheKey);
    if (cached != preservationCache->end()) {
      return cached->second;
    }
  }
  const auto recordResult = [&](bool result) {
    if (preservationCache != nullptr) {
      preservationCache->emplace(cacheKey, result);
    }
    return result;
  };

  if (const auto constantValue = constantBoolExprValue(targetExpr);
      constantValue.has_value()) {
    // Constant next-state functions do not need a SAT proof. This case is
    // common in reset-pruned dual-rail cones, and avoiding thousands of tiny
    // KISSAT instances keeps local auxiliary mining from becoming the bottleneck.
    return recordResult(*constantValue == targetValue);
  }
  const auto& support = resolver.support(targetSymbol);
  if (const auto knownValue =
          evaluateBoolExprWithKnownSupport(
              targetExpr, support, constants);
      knownValue.has_value()) {
    // If the current candidate invariant already determines the transition
    // expression, the fixed-point local auxiliary loop can prove the constant
    // without opening another SAT solver.  Candidates that only hold because of
    // an invalid assumption are removed on the next fixed-point pass.
    return recordResult(*knownValue == targetValue);
  }
  if (support.size() > kAuxiliaryInvariantSupportLimit) {
    return recordResult(false);
  }

  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.configureForSecLocalBooleanCheck(support.size());
  std::unordered_map<size_t, int> leaves;
  leaves.reserve(support.size());
  for (const size_t supportSymbol : support) {
    if (supportSymbol >= 2) {
      leaves.emplace(supportSymbol, solver.newVar() + 2);
    }
  }

  AuxiliaryStateInvariants assumptions;
  assumptions.constants.assign(constants.begin(), constants.end());
  addIndexedStateSemantics(
      solver, stateSemanticsIndex, leaves, ClausePartition::A);
  addAuxiliaryStateInvariants(
      solver, leaves, assumptions, ClausePartition::A);

  FrameFormulaEncoder encoder(
      solver, std::move(leaves), /*createMissingLeaves=*/true);
  const int nextValue = encoder.encode(targetExpr);
  solver.addClause({targetValue ? -nextValue : nextValue});
  return recordResult(
      solver.solveStatus() == SATSolverWrapper::SolveStatus::Unsat);
}

bool addTransitionSupportWithinLimit(
    const TransitionExprResolver& resolver,
    size_t symbol,
    std::unordered_set<size_t>& support) {
  if (!resolver.contains(symbol)) {
    return false; // LCOV_EXCL_LINE
  }
  for (const size_t supportSymbol : resolver.support(symbol)) {
    if (supportSymbol < 2) {
      continue; // LCOV_EXCL_LINE
    }
    support.insert(supportSymbol);
    if (support.size() > kAuxiliaryInvariantSupportLimit) {
      return false; // LCOV_EXCL_LINE
    }
  }
  return true;
}

std::pair<size_t, size_t> canonicalStateEqualityPair(size_t lhs, size_t rhs) {
  return lhs < rhs ? std::make_pair(lhs, rhs) : std::make_pair(rhs, lhs);
}

bool bothStatesHaveSameConstant(
    const std::unordered_map<size_t, bool>& constants,
    size_t lhs,
    size_t rhs) {
  const auto lhsIt = constants.find(lhs);
  const auto rhsIt = constants.find(rhs);
  return lhsIt != constants.end() && rhsIt != constants.end() &&
         lhsIt->second == rhsIt->second;
}

size_t appendAuxiliaryEqualityCandidatesForDesign(
    const TransitionExprResolver& resolver,
    const std::vector<size_t>& designStates,
    const std::unordered_set<size_t>* eligibleSymbols,
    const std::unordered_map<size_t, bool>& bootstrapValues,
    const std::unordered_map<size_t, bool>& constants,
    std::vector<std::pair<size_t, size_t>>& candidates) {
  std::unordered_map<
      ExprBootstrapValueKey,
      std::vector<size_t>,
      ExprBootstrapValueKeyHash>
      buckets;
  for (const size_t symbol : designStates) {
    if (eligibleSymbols != nullptr && !eligibleSymbols->contains(symbol)) {
      continue;
    }
    const auto value = bootstrapValues.find(symbol);
    if (value == bootstrapValues.end() || !resolver.contains(symbol)) {
      continue;
    }
    if (resolver.support(symbol).size() > kAuxiliaryInvariantSupportLimit) {
      continue;
    }
    BoolExpr* expr = resolver.at(symbol);
    if (expr == nullptr) {
      continue; // LCOV_EXCL_LINE
    }
    buckets[{expr, value->second}].push_back(symbol);
  }

  size_t candidateCount = 0;
  for (auto& [key, bucket] : buckets) {
    (void)key;
    if (bucket.size() < 2) {
      continue;
    }
    std::sort(bucket.begin(), bucket.end());
    const size_t anchor = bucket.front();
    for (size_t index = 1; index < bucket.size(); ++index) {
      if (bothStatesHaveSameConstant(constants, anchor, bucket[index])) {
        continue;
      }
      ++candidateCount;
      if (candidates.size() < kAuxiliaryEqualityCandidateLimit) {
        candidates.push_back(
            canonicalStateEqualityPair(anchor, bucket[index]));
      }
    }
  }
  return candidateCount;
}

size_t appendAuxiliaryEqualityCandidates(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const std::unordered_set<size_t>* eligibleSymbols,
    const std::unordered_map<size_t, bool>& bootstrapValues,
    const std::unordered_map<size_t, bool>& constants,
    std::vector<std::pair<size_t, size_t>>& candidates) {
  size_t candidateCount = appendAuxiliaryEqualityCandidatesForDesign(
      resolver, problem.state0Symbols, eligibleSymbols, bootstrapValues,
      constants, candidates);
  candidateCount += appendAuxiliaryEqualityCandidatesForDesign(
      resolver, problem.state1Symbols, eligibleSymbols, bootstrapValues,
      constants, candidates);
  std::sort(candidates.begin(), candidates.end());
  candidates.erase(
      std::unique(candidates.begin(), candidates.end()), candidates.end());
  return candidateCount;
}

bool transitionPreservesStateEquality(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const StateSemanticsIndex& stateSemanticsIndex,
    size_t lhsSymbol,
    size_t rhsSymbol,
    const std::unordered_map<size_t, bool>& constants,
    const std::vector<std::pair<size_t, size_t>>& equalities) {
  if (lhsSymbol == rhsSymbol) {
    return true; // LCOV_EXCL_LINE
  }
  if (!resolver.contains(lhsSymbol) || !resolver.contains(rhsSymbol)) {
    return false; // LCOV_EXCL_LINE
  }
  BoolExpr* lhsExpr = resolver.at(lhsSymbol);
  BoolExpr* rhsExpr = resolver.at(rhsSymbol);
  if (lhsExpr == nullptr || rhsExpr == nullptr) {
    return false; // LCOV_EXCL_LINE
  }
  std::unordered_set<size_t> support;
  if (!addTransitionSupportWithinLimit(resolver, lhsSymbol, support) ||
      !addTransitionSupportWithinLimit(resolver, rhsSymbol, support)) {
    return false; // LCOV_EXCL_LINE
  }

  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.configureForSecLocalBooleanCheck(support.size());
  auto leaves = allocateLeafLits(
      solver, sortedSymbols(support), VariablePartition::ALocal);

  AuxiliaryStateInvariants assumptions;
  assumptions.constants.assign(constants.begin(), constants.end());
  assumptions.equalities = equalities;
  // Auxiliary equalities are only trusted after this local induction check.
  // Same-design rail semantics are unconditional, so the proof may rely on
  // them just like every Craig image/fixed-point query does.
  addIndexedStateSemantics(
      solver, stateSemanticsIndex, leaves, ClausePartition::A);
  addAuxiliaryStateInvariants(
      solver, leaves, assumptions, ClausePartition::A);

  FrameFormulaEncoder encoder(
      solver, std::move(leaves), /*createMissingLeaves=*/true);
  const int lhsNext = encoder.encode(lhsExpr);
  const int rhsNext = encoder.encode(rhsExpr);
  solver.addClause({lhsNext, rhsNext});
  solver.addClause({-lhsNext, -rhsNext});
  return solver.solveStatus() == SATSolverWrapper::SolveStatus::Unsat;
}

AuxiliaryStateInvariants deriveAuxiliaryStateInvariants(
    const KInductionProblem& problem,
    bool forceEnabled) {
  if ((!forceEnabled && !imcAuxiliaryInvariantsEnabled()) ||
      problem.bootstrapStateAssignments.empty()) {
    return {};
  }

  const size_t candidateEstimate = problem.bootstrapStateAssignments.size();
  if (candidateEstimate > kAuxiliaryInvariantCandidateLimit) {
    // Auxiliary invariants are a pruning optimization.  Probing millions of
    // bootstrap entries can exceed CI memory before Craig IMC starts, so huge
    // dual-rail problems skip this optional pass and let the strict IMC
    // projection loop refine only the current proof slice.  Use the raw
    // assignment count as an upper bound so the skip itself stays allocation
    // free on BP-sized designs.
    emitSecDiag(
        "SEC diag: imc Craig skips auxiliary invariants candidates=",
        candidateEstimate,
        " candidate_limit=", kAuxiliaryInvariantCandidateLimit);
    return {};
  }

  auto constants = bootstrapStateConstants(problem);
  const size_t candidateCount = constants.size();
  const TransitionExprResolver resolver(problem);
  const auto complementPrimary = primaryByComplement(problem);
  const StateSemanticsIndex stateSemanticsIndex =
      buildStateSemanticsIndex(problem);
  closeComplementStateConstants(problem, constants);
  bool changed = true;
  while (changed) {
    changed = false;
    StateConstantPreservationCache preservationCache;
    preservationCache.reserve(constants.size());
    for (auto it = constants.begin(); it != constants.end();) {
      if (transitionPreservesStateConstant(
              problem,
              resolver,
              stateSemanticsIndex,
              it->first,
              it->second,
              constants,
              complementPrimary,
              &preservationCache)) {
        ++it;
        continue;
      }
      it = constants.erase(it);
      preservationCache.clear();
      changed = true;
    }
  }

  std::vector<std::pair<size_t, size_t>> equalityCandidates;
  const size_t equalityCandidateCount = appendAuxiliaryEqualityCandidates(
      problem, resolver, /*eligibleSymbols=*/nullptr,
      bootstrapStateConstants(problem), constants,
      equalityCandidates);
  std::vector<std::pair<size_t, size_t>> equalities = equalityCandidates;
  changed = true;
  while (changed) {
    changed = false;
    for (auto it = equalities.begin(); it != equalities.end();) {
      if (transitionPreservesStateEquality(
              problem,
              resolver,
              stateSemanticsIndex,
              it->first,
              it->second,
              constants,
              equalities)) {
        ++it;
        continue;
      }
      it = equalities.erase(it); // LCOV_EXCL_LINE
      changed = true; // LCOV_EXCL_LINE
    }
  }

  AuxiliaryStateInvariants invariants;
  invariants.constants.assign(constants.begin(), constants.end());
  invariants.equalities = std::move(equalities);
  closeAuxiliaryInvariantConstantComplements(problem, invariants);
  std::sort(invariants.equalities.begin(), invariants.equalities.end());
  emitSecDiag(
      "SEC diag: imc Craig auxiliary constants=", invariants.constants.size(),
      " candidates=", candidateCount,
      " equalities=", invariants.equalities.size(),
      " equality_candidates=", equalityCandidateCount,
      " equality_candidate_limit=", kAuxiliaryEqualityCandidateLimit,
      " support_limit=", kAuxiliaryInvariantSupportLimit);
  return invariants;
}

std::unordered_map<size_t, bool> localAuxiliaryCandidateConstants(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& candidateSymbols,
    const AuxiliaryStateInvariants& existingInvariants,
    size_t& candidateCount,
    size_t& selectedCandidateCount,
    size_t& selectedCandidateLimit,
    size_t& topCandidateScore) {
  candidateCount = 0;
  selectedCandidateCount = 0;
  selectedCandidateLimit = 0;
  topCandidateScore = 0;
  if (candidateSymbols.empty() ||
      problem.bootstrapStateAssignments.empty()) {
    return {};
  }

  std::unordered_set<size_t> existingConstants;
  existingConstants.reserve(existingInvariants.constants.size());
  for (const auto& [symbol, value] : existingInvariants.constants) {
    (void)value;
    existingConstants.insert(symbol);
  }

  std::unordered_map<size_t, bool> constants;
  constants.reserve(std::min(
      candidateSymbols.size(), problem.bootstrapStateAssignments.size()));
  for (const auto& [symbol, value] : problem.bootstrapStateAssignments) {
    if (candidateSymbols.contains(symbol) &&
        !existingConstants.contains(symbol)) {
      constants.emplace(symbol, value);
    }
  }
  candidateCount = constants.size();
  selectedCandidateLimit = localAuxiliaryCandidateLimit(candidateCount);
  if (constants.size() <= selectedCandidateLimit) {
    selectedCandidateCount = constants.size();
    closeComplementStateConstants(problem, constants);
    return constants;
  }

  std::unordered_map<size_t, bool> scoringConstants = constants;
  closeComplementStateConstants(problem, scoringConstants);
  const bool usePreservedExtraSlice =
      candidateCount > kLocalAuxiliaryInvariantLargeCandidateThreshold;
  std::unordered_map<size_t, size_t> faninScore;
  size_t scoredRequestCount = 0;
  for (const size_t requested : sortedSymbols(candidateSymbols)) {
    const size_t target =
        transitionTargetFor(requested, resolver, complementPrimary);
    if (!resolver.contains(target)) {
      continue;
    }
    ++scoredRequestCount;
    for (const size_t supportSymbol : resolver.support(target)) {
      if (constants.contains(supportSymbol)) {
        ++faninScore[supportSymbol];
      }
    }
  }

  std::vector<ProjectionRefinementCandidate> candidates;
  std::vector<ProjectionRefinementCandidate> preservedCandidates;
  candidates.reserve(constants.size());
  for (const auto& [symbol, value] : constants) {
    const auto fanin = faninScore.find(symbol);
    const size_t score = fanin == faninScore.end() ? 0 : fanin->second;
    candidates.push_back({symbol, score});
    if (usePreservedExtraSlice && score != 0) {
      const size_t preservationScore = localAuxiliaryPreservationScore(
          symbol, value, resolver, complementPrimary, scoringConstants);
      if (preservationScore != 0) {
        preservedCandidates.push_back({symbol, preservationScore + score});
      }
    }
  }
  std::sort(
      candidates.begin(),
      candidates.end(),
      projectionRefinementCandidateBetter);
  topCandidateScore = candidates.empty() ? 0 : candidates.front().score;

  std::unordered_map<size_t, bool> selected;
  selected.reserve(
      selectedCandidateLimit +
      (usePreservedExtraSlice ? kLocalAuxiliaryPreservedCandidateExtraLimit
                              : 0));
  for (const ProjectionRefinementCandidate& candidate : candidates) {
    if (selected.size() >= selectedCandidateLimit) {
      break;
    }
    selected.emplace(candidate.symbol, constants.at(candidate.symbol));
  }
  if (usePreservedExtraSlice) {
    std::sort(
        preservedCandidates.begin(),
        preservedCandidates.end(),
        projectionRefinementCandidateBetter);
    const size_t selectedWithPreservedLimit =
        selectedCandidateLimit + kLocalAuxiliaryPreservedCandidateExtraLimit;
    for (const ProjectionRefinementCandidate& candidate : preservedCandidates) {
      if (selected.size() >= selectedWithPreservedLimit) {
        break; // LCOV_EXCL_LINE
      }
      if (selected.contains(candidate.symbol)) {
        continue;
      }
      // Preserve the original high-fan-in slice, then add candidates that the
      // cheap constant/known-support screen expects to prove without opening
      // many new SAT instances. The later induction filter remains the proof
      // gate.
      selected.emplace(candidate.symbol, constants.at(candidate.symbol));
    }
  }
  selectedCandidateCount = selected.size();
  closeComplementStateConstants(problem, selected);
  return selected;
}

AuxiliaryStateInvariants deriveLocalAuxiliaryStateInvariants(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const StateSemanticsIndex& stateSemanticsIndex,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& candidateSymbols,
    const AuxiliaryStateInvariants& existingInvariants) {
  if (!shouldMineLocalAuxiliaryInvariants(candidateSymbols.size())) {
    // At BP scale the local auxiliary pass can spend gigabytes validating a
    // small handful of constants.  Once the Craig transition support is this
    // large, direct projection refinement is the cheaper strict-IMC move.
    emitSecDiag(
        "SEC diag: imc Craig skips local auxiliary invariants support=",
        candidateSymbols.size(),
        " support_limit=", kLocalAuxiliaryInvariantSupportMiningLimit);
    return {};
  }

  size_t candidateCount = 0;
  size_t selectedCandidateCount = 0;
  size_t selectedCandidateLimit = 0;
  size_t topCandidateScore = 0;
  std::unordered_map<size_t, bool> localConstants =
      localAuxiliaryCandidateConstants(
          problem,
          resolver,
          complementPrimary,
          candidateSymbols,
          existingInvariants,
          candidateCount,
          selectedCandidateCount,
          selectedCandidateLimit,
          topCandidateScore);
  if (localConstants.empty()) {
    return {};
  }
  std::unordered_set<size_t> selectedSymbols;
  selectedSymbols.reserve(localConstants.size());
  for (const auto& [symbol, value] : localConstants) {
    (void)value;
    selectedSymbols.insert(symbol);
  }

  std::unordered_map<size_t, bool> assumptions;
  assumptions.reserve(
      existingInvariants.constants.size() + localConstants.size());
  for (const auto& [symbol, value] : existingInvariants.constants) {
    assumptions[symbol] = value;
  }
  for (const auto& [symbol, value] : localConstants) {
    assumptions[symbol] = value;
  }
  closeComplementStateConstants(problem, assumptions);

  bool changed = true;
  while (changed) {
    changed = false;
    StateConstantPreservationCache preservationCache;
    preservationCache.reserve(localConstants.size());
    for (auto it = localConstants.begin(); it != localConstants.end();) {
      if (transitionPreservesStateConstant(
              problem,
              resolver,
              stateSemanticsIndex,
              it->first,
              it->second,
              assumptions,
              complementPrimary,
              &preservationCache)) {
        ++it;
        continue;
      }
      assumptions.erase(it->first);
      it = localConstants.erase(it);
      preservationCache.clear();
      changed = true;
    }
  }

  std::vector<std::pair<size_t, size_t>> equalityCandidates;
  const size_t equalityCandidateCount = appendAuxiliaryEqualityCandidates(
      problem,
      resolver,
      &selectedSymbols,
      bootstrapAssignmentsForSymbols(
          problem.bootstrapStateAssignments, selectedSymbols),
      assumptions,
      equalityCandidates);
  std::set<std::pair<size_t, size_t>> existingEqualities;
  for (const auto& equality : existingInvariants.equalities) {
    existingEqualities.insert(canonicalStateEqualityPair(
        equality.first, equality.second));
  }
  std::vector<std::pair<size_t, size_t>> localEqualities;
  localEqualities.reserve(equalityCandidates.size());
  for (const auto& equality : equalityCandidates) {
    const auto canonical =
        canonicalStateEqualityPair(equality.first, equality.second);
    if (!existingEqualities.contains(canonical)) {
      localEqualities.push_back(canonical);
    }
  }
  std::sort(localEqualities.begin(), localEqualities.end());
  localEqualities.erase(
      std::unique(localEqualities.begin(), localEqualities.end()),
      localEqualities.end());

  std::vector<std::pair<size_t, size_t>> equalityAssumptions =
      existingInvariants.equalities;
  equalityAssumptions.insert(
      equalityAssumptions.end(), localEqualities.begin(),
      localEqualities.end());
  std::sort(equalityAssumptions.begin(), equalityAssumptions.end());
  equalityAssumptions.erase(
      std::unique(equalityAssumptions.begin(), equalityAssumptions.end()),
      equalityAssumptions.end());
  changed = true;
  while (changed) {
    changed = false;
    for (auto it = localEqualities.begin(); it != localEqualities.end();) {
      if (transitionPreservesStateEquality(
              problem,
              resolver,
              stateSemanticsIndex,
              it->first,
              it->second,
              assumptions,
              equalityAssumptions)) {
        ++it;
        continue;
      }
      const auto removed = *it; // LCOV_EXCL_LINE
      it = localEqualities.erase(it); // LCOV_EXCL_LINE
      const auto assumption = std::lower_bound( // LCOV_EXCL_LINE
          equalityAssumptions.begin(), equalityAssumptions.end(), removed); // LCOV_EXCL_LINE
      if (assumption != equalityAssumptions.end() && // LCOV_EXCL_LINE
          *assumption == removed) { // LCOV_EXCL_LINE
        equalityAssumptions.erase(assumption); // LCOV_EXCL_LINE
      } // LCOV_EXCL_LINE
      changed = true; // LCOV_EXCL_LINE
    }
  }

  AuxiliaryStateInvariants invariants;
  invariants.constants.assign(localConstants.begin(), localConstants.end());
  invariants.equalities = std::move(localEqualities);
  closeAuxiliaryInvariantConstantComplements(problem, invariants);
  std::sort(invariants.equalities.begin(), invariants.equalities.end());
  if (!invariants.empty()) {
    // BP-sized designs cannot afford global auxiliary mining over millions of
    // bootstrap assignments.  These local constants/equalities are still
    // transition-proven, but only for the current Craig projection support.
    emitSecDiag(
        "SEC diag: imc Craig local auxiliary constants=",
        invariants.constants.size(),
        " candidates=", candidateCount,
        " selected=", selectedCandidateCount,
        " candidate_limit=", selectedCandidateLimit,
        " top_score=", topCandidateScore,
        " support_limit=", kAuxiliaryInvariantSupportLimit,
        " equalities=", invariants.equalities.size(),
        " equality_candidates=", equalityCandidateCount,
        " equality_candidate_limit=", kAuxiliaryEqualityCandidateLimit);
  }
  return invariants;
}

size_t mergeAuxiliaryStateInvariants(
    AuxiliaryStateInvariants& target,
    const AuxiliaryStateInvariants& additions) {
  size_t added = 0;
  std::unordered_set<size_t> existingConstants;
  existingConstants.reserve(target.constants.size());
  for (const auto& [symbol, value] : target.constants) {
    (void)value;
    existingConstants.insert(symbol);
  }
  for (const auto& [symbol, value] : additions.constants) {
    if (existingConstants.insert(symbol).second) {
      target.constants.push_back({symbol, value});
      ++added;
    }
  }
  std::sort(target.constants.begin(), target.constants.end());
  std::set<std::pair<size_t, size_t>> existingEqualities;
  for (const auto& equality : target.equalities) {
    existingEqualities.insert(canonicalStateEqualityPair(
        equality.first, equality.second));
  }
  for (const auto& equality : additions.equalities) {
    const auto canonical =
        canonicalStateEqualityPair(equality.first, equality.second);
    if (existingEqualities.insert(canonical).second) {
      target.equalities.push_back(canonical);
      ++added;
    }
  }
  std::sort(target.equalities.begin(), target.equalities.end());
  return added;
}

bool promoteLocalAuxiliaryInvariants(
    AuxiliaryStateInvariants& activeInvariants,
    AuxiliaryStateInvariants& projectionInvariants,
    const AuxiliaryStateInvariants& localInvariants) {
  if (mergeAuxiliaryStateInvariants(activeInvariants, localInvariants) == 0) {
    return false;
  }
  const size_t promotedCount =
      mergeAuxiliaryStateInvariants(projectionInvariants, localInvariants);
  if (promotedCount != 0) {
    // Local auxiliaries have already passed the same transition-preservation
    // filter as global auxiliaries. Keeping them for later projection rounds
    // avoids re-mining identical BP slices.
    emitSecDiag(
        "SEC diag: imc Craig promotes local auxiliary invariants added=",
        promotedCount,
        " total_constants=", projectionInvariants.constants.size(),
        " total_equalities=", projectionInvariants.equalities.size());
  }
  return true;
}

AuxiliaryStateInvariants helperAuxiliaryStateInvariantsFromOptions(
    const KInductionProblem& problem,
    const CraigImcOptions& options) {
  AuxiliaryStateInvariants invariants;
  invariants.constants = options.helperAuxiliaryStateInvariants;
  invariants.equalities = options.helperAuxiliaryStateEqualities;
  const size_t rawCount = invariants.constants.size() +
                          invariants.equalities.size();
  if (rawCount > kCraigHelperAuxiliaryCarryLimit) {
    // Helper auxiliaries are useful when they are still a compact proof
    // certificate.  BP's retained-helper tail needs its moderate mined packet,
    // but very broad packets make every Craig image carry a second proof.
    emitSecDiag(
        "SEC diag: imc Craig skips broad helper auxiliary invariants "
        "constants=",
        invariants.constants.size(),
        " equalities=", invariants.equalities.size(),
        " limit=", kCraigHelperAuxiliaryCarryLimit);
    return {};
  }
  closeAuxiliaryInvariantConstantComplements(problem, invariants);
  std::sort(invariants.constants.begin(), invariants.constants.end());
  std::sort(invariants.equalities.begin(), invariants.equalities.end());
  invariants.constants.erase(
      std::unique(invariants.constants.begin(), invariants.constants.end()),
      invariants.constants.end());
  invariants.equalities.erase(
      std::unique(invariants.equalities.begin(), invariants.equalities.end()),
      invariants.equalities.end());
  const size_t normalizedCount = invariants.constants.size() +
                                 invariants.equalities.size();
  if (normalizedCount > kCraigHelperAuxiliaryCarryLimit) {
    emitSecDiag( // LCOV_EXCL_LINE
        "SEC diag: imc Craig skips broad helper auxiliary invariants "
        "constants=",
        invariants.constants.size(), // LCOV_EXCL_LINE
        " equalities=", invariants.equalities.size(), // LCOV_EXCL_LINE
        " limit=", kCraigHelperAuxiliaryCarryLimit);
    return {}; // LCOV_EXCL_LINE
  }
  return invariants;
}

CraigImcResult makeCraigEquivalentResult(
    size_t iterations,
    std::vector<InterpolantRegion> invariantRegions,
    std::unordered_set<size_t> trackedStates,
    const AuxiliaryStateInvariants& auxiliaryStateInvariants) {
  CraigImcResult result;
  result.status = CraigImcStatus::Equivalent;
  result.iterations = iterations;
  result.invariantRegions = std::move(invariantRegions);
  result.trackedStates = std::move(trackedStates);
  result.auxiliaryStateInvariants = auxiliaryStateInvariants.constants;
  result.auxiliaryStateEqualities = auxiliaryStateInvariants.equalities;
  return result;
}

void addResetInputValue(
    SATSolverWrapper& solver,
    const KInductionProblem& problem,
    std::unordered_map<size_t, int>& leaves,
    bool asserted,
    VariablePartition variablePartition,
    ClausePartition clausePartition) {
  solver.setCraigVariablePartition(variablePartition);
  solver.setCraigClausePartition(clausePartition);
  for (const auto& [symbol, assertedValue] : problem.resetBootstrapInputs) {
    auto [it, inserted] = leaves.emplace(symbol, 0);
    if (inserted) {
      it->second = solver.newVar() + 2;
    }
    const bool value = asserted ? assertedValue : !assertedValue;
    solver.addClause({value ? it->second : -it->second});
  }
}

std::unordered_map<size_t, size_t> primaryByComplement(
    const KInductionProblem& problem) {
  std::unordered_map<size_t, size_t> result;
  for (const auto& [primary, complement] : problem.complementedStatePairs0) {
    result.emplace(complement, primary);
  }
  for (const auto& [primary, complement] : problem.complementedStatePairs1) {
    result.emplace(complement, primary); // LCOV_EXCL_LINE
  }
  return result;
}

size_t transitionTargetFor(
    size_t symbol,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary) {
  if (resolver.contains(symbol)) {
    return symbol;
  }
  const auto primary = complementPrimary.find(symbol);
  if (primary != complementPrimary.end() && resolver.contains(primary->second)) {
    return primary->second;
  }
  return symbol;
}

std::unordered_set<size_t> transitionProjectionSupport(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& trackedStates) {
  const auto states = stateSymbolSet(problem);
  std::unordered_set<size_t> support = trackedStates;
  for (const size_t requested : sortedSymbols(trackedStates)) {
    const size_t target =
        transitionTargetFor(requested, resolver, complementPrimary);
    // addProjectedTransition only encodes targets that also have a next-state
    // leaf in the current projection. Match that behavior here so batching
    // sees the same closure Craig IMC will later refine toward.
    if (!trackedStates.contains(target) || !resolver.contains(target)) {
      continue;
    }
    for (const size_t symbol : resolver.support(target)) {
      if (states.contains(symbol)) {
        support.insert(symbol);
      }
    }
  }
  closeSameDesignStateSemantics(problem, support);
  return support;
}

size_t countTransitionExpressionNodes(
    BoolExpr* root,
    std::unordered_set<BoolExpr*>& visited) {
  if (root == nullptr) {
    return 0; // LCOV_EXCL_LINE
  }
  const size_t before = visited.size();
  std::vector<BoolExpr*> stack{root};
  while (!stack.empty()) {
    BoolExpr* node = stack.back();
    stack.pop_back();
    if (node == nullptr || !visited.insert(node).second) {
      continue;
    }
    if (node->getLeft() != nullptr) {
      stack.push_back(node->getLeft());
    }
    if (node->getRight() != nullptr) {
      stack.push_back(node->getRight());
    }
  }
  return visited.size() - before;
}

ProjectedTransitionBuildEstimate estimateProjectedTransitionBuild(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const std::unordered_set<size_t>& states,
    const std::unordered_set<size_t>& trackedStates,
    std::vector<ProjectedTransitionFrame>& transitionFrames) {
  std::unordered_set<size_t> support = trackedStates;
  ProjectedTransitionBuildEstimate estimate;
  for (ProjectedTransitionFrame& frame : transitionFrames) {
    frame.expressionNodes = 0;
    std::unordered_set<BoolExpr*> frameExpressionNodes;
    for (const size_t target : frame.targets) {
      ++estimate.encodedTargets;
      const size_t expressionNodes = countTransitionExpressionNodes(
          resolver.at(target), frameExpressionNodes);
      frame.expressionNodes += expressionNodes;
      estimate.expressionNodes += expressionNodes;
      for (const size_t symbol : resolver.support(target)) {
        if (states.contains(symbol)) {
          support.insert(symbol);
        }
      }
    }
  }
  closeSameDesignStateSemantics(problem, support);
  estimate.stateSupport = support.size();
  return estimate;
}

const ProjectedTransitionPlan& projectedTransitionPlanForLookahead(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& states,
    const std::unordered_set<size_t>& trackedStates,
    size_t helperInvariantRegionCount,
    size_t lookahead,
    ProjectedTransitionPlanCache& cache) {
  auto [planIt, inserted] = cache.try_emplace(lookahead);
  ProjectedTransitionPlan& plan = planIt->second;
  if (!inserted) {
    return plan;
  }

  plan.frames.reserve(lookahead);
  for (size_t frame = 0; frame < lookahead; ++frame) {
    ProjectedTransitionFrame transitionFrame;
    transitionFrame.requests = focusedImageTransitionRequests(
        problem,
        resolver,
        complementPrimary,
        states,
        trackedStates,
        helperInvariantRegionCount,
        lookahead - frame);
    transitionFrame.targets = projectedTransitionTargets(
        resolver,
        complementPrimary,
        trackedStates,
        transitionFrame.requests);
    plan.largestTransitionRequestCount = std::max(
        plan.largestTransitionRequestCount,
        transitionFrame.requests.size());
    plan.frames.push_back(std::move(transitionFrame));
  }
  plan.estimate = estimateProjectedTransitionBuild(
      problem,
      resolver,
      states,
      trackedStates,
      plan.frames);
  return plan;
}

TransitionEncodingResult addProjectedTransition(
    SATSolverWrapper& solver,
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const StateSemanticsIndex& stateSemanticsIndex,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<size_t>& transitionTargets,
    const AuxiliaryStateInvariants& auxiliaryInvariants,
    TransitionEncodingCache& transitionEncodingCache,
    const std::unordered_map<size_t, int>& currentLits,
    const std::unordered_map<size_t, int>& nextStateLits,
    VariablePartition localVariablePartition,
    ClausePartition clausePartition,
    size_t transitionExpressionNodeHint = 0) {
  solver.setCraigVariablePartition(localVariablePartition);
  // FrameFormulaEncoder emits Tseitin clauses while encoding each transition.
  // Select the transition's Craig side before encoding, not only before adding
  // the final next-state equivalence.
  solver.setCraigClausePartition(clausePartition);
  std::unordered_map<size_t, int> leaves;
  if (transitionEncodingCache.constantAssignments.empty()) {
    FrameFormulaEncoder encoder(
        solver,
        std::unordered_map<size_t, int>(currentLits),
        /*createMissingLeaves=*/true,
        transitionExpressionNodeHint);
    for (const size_t target : transitionTargets) {
      const auto next = nextStateLits.find(target);
      if (next == nextStateLits.end()) {
        continue; // LCOV_EXCL_LINE
      }
      solver.setCraigVariablePartition(localVariablePartition);
      solver.setCraigClausePartition(clausePartition);
      const int transitionLit = encoder.encode(resolver.at(target));
      addLiteralEquivalenceForPartition(
          solver, next->second, transitionLit, clausePartition);
    }
    leaves = encoder.leafLits();
  } else {
    ConstantAwareFrameFormulaEncoder encoder(
        solver,
        std::unordered_map<size_t, int>(currentLits),
        transitionEncodingCache.constantAssignments,
        localVariablePartition,
        clausePartition,
        transitionExpressionNodeHint);
    for (const size_t target : transitionTargets) {
      const auto next = nextStateLits.find(target);
      if (next == nextStateLits.end()) {
        continue; // LCOV_EXCL_LINE
      }
      solver.setCraigVariablePartition(localVariablePartition);
      solver.setCraigClausePartition(clausePartition);
      const int transitionLit = encoder.encode(resolver.at(target));
      addLiteralEquivalenceForPartition(
          solver, next->second, transitionLit, clausePartition);
    }
    leaves = encoder.leafLits();
  }

  addResetInputValue(
      solver,
      problem,
      leaves,
      /*asserted=*/false,
      localVariablePartition,
      clausePartition);
  if (shouldEncodeProjectedLocalStateSemantics(
          leaves.size(), trackedStates.size())) {
    addIndexedStateSemantics(
        solver, stateSemanticsIndex, leaves, clausePartition);
  } else {
    // The transition support leaves are local to one Craig side.  Dropping
    // local same-frame rail semantics weakens the projected query, so every
    // UNSAT interpolant remains a strict IMC proof while avoiding a large
    // proof-trace clause burst on BP-sized cones.
    emitSecDiag(
        "SEC diag: imc Craig weakens projected local state semantics leaves=",
        leaves.size(),
        " leaf_limit=", kCraigProjectedLocalSemanticsSupportLimit,
        " tracked_states=", trackedStates.size(),
        " tracked_limit=", kCraigProjectedLocalSemanticsTrackedStateLimit);
  }
  addAuxiliaryStateInvariants(
      solver, leaves, auxiliaryInvariants, clausePartition);
  return {std::move(leaves)};
}

RegionLiteral convertInterpolantLiteral(
    int literal,
    const std::unordered_map<int, size_t>& stateByVariable,
    int firstAuxiliaryVariable,
    std::unordered_map<int, size_t>& auxiliaryByVariable) {
  const int variable = std::abs(literal);
  if (const auto state = stateByVariable.find(variable);
      state != stateByVariable.end()) {
    return {true, state->second, literal > 0};
  }
  if (variable < firstAuxiliaryVariable) {
    throw std::runtime_error( // LCOV_EXCL_LINE
        "Craig interpolant contains a non-global original variable");
  }
  auto [auxiliary, inserted] =
      auxiliaryByVariable.emplace(variable, auxiliaryByVariable.size());
  (void)inserted;
  return {false, auxiliary->second, literal > 0};
}

InterpolantRegion convertInterpolant(
    const SATSolverWrapper::CraigInterpolantCnf& cnf,
    const std::unordered_map<int, size_t>& stateByVariable) {
  if (cnf.type ==
      SATSolverWrapper::CraigInterpolantCnf::Type::ConstantFalse) {
    return {InterpolantRegion::Type::False}; // LCOV_EXCL_LINE
  }
  if (cnf.type ==
      SATSolverWrapper::CraigInterpolantCnf::Type::ConstantTrue) {
    return {InterpolantRegion::Type::True};
  }
  if (cnf.type != SATSolverWrapper::CraigInterpolantCnf::Type::Normal ||
      cnf.clauses.empty() || cnf.clauses.back().size() != 1) {
    throw std::runtime_error("CaDiCaL returned an invalid Craig interpolant CNF"); // LCOV_EXCL_LINE
  }

  InterpolantRegion region;
  region.type = InterpolantRegion::Type::Normal;
  std::unordered_map<int, size_t> auxiliaryByVariable;
  region.definitionClauseEnds.reserve(cnf.clauses.size() - 1);
  size_t inputLiteralCount = 0;
  for (size_t clauseIndex = 0; clauseIndex + 1 < cnf.clauses.size();
       ++clauseIndex) {
    inputLiteralCount += cnf.clauses[clauseIndex].size();
    for (const int literal : cnf.clauses[clauseIndex]) {
      region.definitionLiterals.push_back(convertInterpolantLiteral(
          literal,
          stateByVariable,
          cnf.firstAuxiliaryVariable,
          auxiliaryByVariable));
    }
    region.definitionClauseEnds.push_back(region.definitionLiterals.size());
  }
  region.root = convertInterpolantLiteral(
      cnf.clauses.back().front(),
      stateByVariable,
      cnf.firstAuxiliaryVariable,
      auxiliaryByVariable);
  region.auxiliaryCount = auxiliaryByVariable.size();
  emitSecDiag(
      "SEC diag: imc Craig interpolant converted clauses=",
      regionClauseCount(region),
      " literals=", regionLiteralCount(region),
      " input_literals=", inputLiteralCount,
      " auxiliaries=", region.auxiliaryCount);
  return simplifyCraigInterpolantRegion(std::move(region));
}

int instantiateRegion(
    SATSolverWrapper& solver,
    const InterpolantRegion& region,
    const std::unordered_map<size_t, int>& stateLits,
    VariablePartition variablePartition,
    ClausePartition clausePartition) {
  if (region.type == InterpolantRegion::Type::True) {
    solver.setCraigVariablePartition(variablePartition);
    const int literal = solver.newVar() + 2;
    solver.setCraigClausePartition(clausePartition);
    solver.addClause({literal});
    return literal;
  }
  if (region.type == InterpolantRegion::Type::False) {
    solver.setCraigVariablePartition(variablePartition); // LCOV_EXCL_LINE
    const int literal = solver.newVar() + 2; // LCOV_EXCL_LINE
    solver.setCraigClausePartition(clausePartition); // LCOV_EXCL_LINE
    solver.addClause({-literal}); // LCOV_EXCL_LINE
    return literal; // LCOV_EXCL_LINE
  }

  solver.setCraigVariablePartition(variablePartition);
  std::vector<int> auxiliaryLits;
  auxiliaryLits.reserve(region.auxiliaryCount);
  for (size_t index = 0; index < region.auxiliaryCount; ++index) {
    auxiliaryLits.push_back(solver.newVar() + 2);
  }

  solver.setCraigClausePartition(clausePartition);
  size_t clauseBegin = 0;
  for (const size_t clauseEnd : region.definitionClauseEnds) {
    std::vector<int> clause;
    clause.reserve(clauseEnd - clauseBegin);
    for (size_t index = clauseBegin; index < clauseEnd; ++index) {
      clause.push_back(
          instantiateRegionLiteral(
              region.definitionLiterals[index], stateLits, auxiliaryLits));
    }
    solver.addClause(clause);
    clauseBegin = clauseEnd;
  }
  return instantiateRegionLiteral(region.root, stateLits, auxiliaryLits);
}

void addRegionUnionConstraint(
    SATSolverWrapper& solver,
    const std::vector<InterpolantRegion>& regions,
    const std::unordered_map<size_t, int>& stateLits,
    VariablePartition auxiliaryVariablePartition,
    ClausePartition clausePartition) {
  if (regions.empty()) {
    return;
  }
  std::vector<int> roots;
  roots.reserve(regions.size());
  for (const InterpolantRegion& region : regions) {
    roots.push_back(instantiateRegion(
        solver, region, stateLits, auxiliaryVariablePartition, clausePartition));
  }
  solver.setCraigClausePartition(clausePartition);
  solver.addClause(roots);
}

std::optional<InterpolantRegion> buildConcreteAssignmentRegion(
    const std::vector<std::pair<size_t, bool>>& assignments,
    const std::unordered_set<size_t>& trackedStates) {
  const std::unordered_map<size_t, bool> assignmentBySymbol =
      bootstrapAssignmentsForSymbols(assignments, trackedStates);
  if (assignmentBySymbol.size() < trackedStates.size()) {
    return std::nullopt;
  }

  std::vector<RegionLiteral> cubeLiterals;
  cubeLiterals.reserve(trackedStates.size());
  for (const size_t symbol : sortedSymbols(trackedStates)) {
    const auto assignment = assignmentBySymbol.find(symbol);
    if (assignment == assignmentBySymbol.end()) {
      // A partial assignment is an over-approximation, not the exact concrete
      // post-reset cube required by the bounded counterexample fast path.
      return std::nullopt; // LCOV_EXCL_LINE
    }
    cubeLiterals.push_back({true, symbol, assignment->second});
  }
  if (cubeLiterals.empty()) {
    return std::nullopt; // LCOV_EXCL_LINE
  }

  InterpolantRegion region;
  region.type = InterpolantRegion::Type::Normal;
  region.auxiliaryCount = 1;
  region.root = {false, 0, true};
  region.definitionLiterals.reserve(cubeLiterals.size() * 2 + 1);
  region.definitionClauseEnds.reserve(cubeLiterals.size() + 1);

  // Encode root <-> conjunction(assignments). Both directions are required:
  // the positive root selects the exact initial cube, while a negated root in
  // an inductiveness query denotes the exact complement of that cube.
  for (const RegionLiteral& literal : cubeLiterals) {
    region.definitionLiterals.push_back({false, 0, false});
    region.definitionLiterals.push_back(literal);
    region.definitionClauseEnds.push_back(region.definitionLiterals.size());
  }
  region.definitionLiterals.push_back(region.root);
  for (RegionLiteral literal : cubeLiterals) {
    literal.positive = !literal.positive;
    region.definitionLiterals.push_back(literal);
  }
  region.definitionClauseEnds.push_back(region.definitionLiterals.size());
  return region;
}

int encodeBadFormulaRoot(
    SATSolverWrapper& solver,
    const KInductionProblem& problem,
    std::unordered_map<size_t, int> nextLeaves) {
  solver.setCraigVariablePartition(VariablePartition::BLocal);
  for (const size_t symbol : problem.inputSymbols) {
    if (!nextLeaves.contains(symbol)) {
      nextLeaves.emplace(symbol, solver.newVar() + 2);
    }
  }
  addResetInputValue(
      solver,
      problem,
      nextLeaves,
      /*asserted=*/false,
      VariablePartition::BLocal,
      ClausePartition::B);
  FrameFormulaEncoder encoder(
      solver, std::move(nextLeaves), /*createMissingLeaves=*/true);
  solver.setCraigVariablePartition(VariablePartition::BLocal);
  return encoder.encode(problem.bad);
}

void addBadFormula(
    SATSolverWrapper& solver,
    const KInductionProblem& problem,
    std::unordered_map<size_t, int> nextLeaves) {
  const int bad = encodeBadFormulaRoot(solver, problem, std::move(nextLeaves));
  solver.setCraigClausePartition(ClausePartition::B);
  solver.addClause({bad});
}

void addSuffixBadFormula(
    SATSolverWrapper& solver,
    const KInductionProblem& problem,
    const std::vector<std::unordered_map<size_t, int>>& frameLits,
    size_t firstFrame,
    size_t lastFrame) {
  std::vector<int> badRoots;
  badRoots.reserve(lastFrame - firstFrame + 1);
  for (size_t frame = firstFrame; frame <= lastFrame; ++frame) {
    badRoots.push_back(
        encodeBadFormulaRoot(solver, problem, frameLits[frame]));
  }
  solver.setCraigClausePartition(ClausePartition::B);
  solver.addClause(badRoots);
}

void addInitialFrontierConstraint(
    SATSolverWrapper& solver,
    const KInductionProblem& problem,
    std::unordered_map<size_t, int> leaves) {
  BoolExpr* initial = BoolExpr::createTrue();
  bool hasInitialConstraint = false;
  for (const auto& [symbol, value] : problem.initialStateAssignments) {
    initial = BoolExpr::And(
        initial,
        value ? BoolExpr::Var(symbol)
              : BoolExpr::Not(BoolExpr::Var(symbol)));
    hasInitialConstraint = true;
  }
  if (!problem.hasExplicitInitialState() && !hasInitialConstraint) {
    // Resetless SEC observes both designs from an already-matching top-level
    // frontier. This is an interface property, not an internal state relation.
    initial = problem.property; // LCOV_EXCL_LINE
    hasInitialConstraint = initial != nullptr; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
  if (!hasInitialConstraint) {
    return;
  }

  solver.setCraigVariablePartition(VariablePartition::ALocal);
  FrameFormulaEncoder encoder(
      solver, std::move(leaves), /*createMissingLeaves=*/true);
  const int initialLiteral = encoder.encode(initial);
  solver.setCraigClausePartition(ClausePartition::A);
  solver.addClause({initialLiteral});
}

struct FrontierResult {
  std::optional<InterpolantRegion> region;
  std::unordered_set<size_t> transitionStateSupport;
  std::unordered_set<size_t> modelGuidedTransitionStateSupport;
  size_t modelGuidedTransitionStateCandidates = 0;
  SATSolverWrapper::SolveStatus solveStatus =
      SATSolverWrapper::SolveStatus::Unknown;
  std::int64_t solveElapsedMilliseconds = 0;
  std::int64_t interpolationElapsedMilliseconds = 0;
  bool focusedTransitionProjection = false;
  size_t largestFocusedTransitionRequestCount = 0;
  bool buildBudgetExceeded = false;
  ProjectedTransitionBuildEstimate estimatedTransitionBuild;

  void recordTransitionRequests(
      size_t trackedStateCount,
      const std::unordered_set<size_t>& requests) {
    if (requests.size() >= trackedStateCount) {
      return;
    }
    focusedTransitionProjection = true;
    largestFocusedTransitionRequestCount =
        std::max(largestFocusedTransitionRequestCount, requests.size());
  }

  bool usesFocusedTransitionProjection() const {
    return focusedTransitionProjection;
  }
};

bool shouldTrySelectedLocalAuxiliaryInvariants(
    const FrontierResult& frontier,
    const std::unordered_set<size_t>& selectedSupport,
    size_t retryCount) {
  // Full local mining is disabled on BP-sized supports.  The refinement slice
  // is already capped and proof-relevant, so it is a cheap last chance to prove
  // selected bootstrap constants before permanently growing tracked state.
  return retryCount < kCraigLocalAuxiliaryRetryLimit &&
         !selectedSupport.empty() &&
         frontier.transitionStateSupport.size() >
             kLocalAuxiliaryInvariantSupportMiningLimit;
}

std::unordered_set<size_t> modelGuidedBootstrapProjectionSupport(
    const KInductionProblem& problem,
    const SATSolverWrapper& solver,
    const std::unordered_map<size_t, int>& currentLits,
    const std::unordered_set<size_t>& trackedStates,
    const std::unordered_set<size_t>& transitionStateSupport,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    size_t& candidateCount) {
  candidateCount = 0;
  if (problem.bootstrapStateAssignments.empty()) {
    return {}; // LCOV_EXCL_LINE
  }

  const std::unordered_map<size_t, bool> bootstrapAssignments =
      bootstrapAssignmentsForSymbols(
          problem.bootstrapStateAssignments, transitionStateSupport);
  const size_t refinementLimit =
      modelGuidedProjectionRefinementLimit(transitionStateSupport.size());
  std::unordered_set<size_t> selectedSupport;
  for (const size_t symbol : sortedSymbols(transitionStateSupport)) {
    if (trackedStates.contains(symbol)) {
      continue;
    }
    const auto expected = bootstrapAssignments.find(symbol);
    const auto literal = currentLits.find(symbol);
    if (expected == bootstrapAssignments.end() ||
        literal == currentLits.end()) {
      continue;
    }
    if (solver.getLiteralValue(literal->second) == expected->second) {
      continue;
    }
    ++candidateCount;
    if (selectedSupport.size() < refinementLimit) {
      selectedSupport.insert(symbol);
    }
  }
  closeProjectionRefinementTransitionTargets(
      resolver, complementPrimary, selectedSupport);
  return selectedSupport;
}

const std::unordered_set<size_t>& projectionRefinementSupport(
    const FrontierResult& frontier) {
  return frontier.modelGuidedTransitionStateSupport.empty()
             ? frontier.transitionStateSupport
             : frontier.modelGuidedTransitionStateSupport;
}

std::unordered_set<size_t> nearSaturatedProjectionRemainderSupport(
    const std::unordered_set<size_t>& transitionStateSupport,
    const std::unordered_set<size_t>& trackedStates) {
  std::unordered_set<size_t> remainder;
  for (const size_t symbol : transitionStateSupport) {
    if (trackedStates.contains(symbol)) {
      continue;
    }
    remainder.insert(symbol);
    if (remainder.size() > kCraigNearSaturatedProjectionRemainderLimit) {
      return {};
    }
  }
  return remainder;
}

bool modelGuidedProjectionNeedsBoundedBackfill(
    const FrontierResult& frontier) {
  // A SAT model may point at only a couple of bootstrap mismatches even though
  // the untracked transition cone is BP-sized.  Backfill those tiny slices with
  // the scored projection picker so one or two model bits do not force a full
  // Craig rebuild round.
  return !frontier.modelGuidedTransitionStateSupport.empty() &&
         frontier.transitionStateSupport.size() >
             kCraigTinyModelGuidedBackfillSupportThreshold &&
         frontier.modelGuidedTransitionStateSupport.size() <
             modelGuidedProjectionRefinementLimit(
                 frontier.transitionStateSupport.size());
}

bool shouldUseSmallLowScoreRefinementStride(
    size_t transitionSupportSize,
    size_t topCandidateScore) {
  // Once the remaining BP support surface has only tiny fan-in scores, smaller
  // strict Craig projection slices are more valuable than fewer rebuild rounds:
  // each imported state stays live in the later proof queries.
  return transitionSupportSize > kCraigLowScoreBackfillSupportThreshold &&
         topCandidateScore <= kCraigLowScoreBackfillScoreLimit;
}

bool shouldKeepOnlyModelGuidedRefinement(
    const FrontierResult& frontier,
    size_t boundedTopCandidateScore) {
  return !frontier.modelGuidedTransitionStateSupport.empty() &&
         shouldUseSmallLowScoreRefinementStride(
             frontier.transitionStateSupport.size(), boundedTopCandidateScore);
}

class ProjectionRefinementScorer {
 public:
  ProjectionRefinementScorer(
      const TransitionExprResolver& resolver,
      const std::unordered_map<size_t, size_t>& complementPrimary)
      : resolver_(resolver), complementPrimary_(complementPrimary) {}

  std::unordered_set<size_t> selectBoundedSupport(
      const std::unordered_set<size_t>& trackedStates,
      const std::unordered_set<size_t>& transitionStateSupport,
      bool focusedTransitionProjection,
      size_t& candidateCount,
      size_t& topCandidateScore,
      size_t& refinementLimit,
      bool& frozeScoreUpdates) {
    candidateCount = 0;
    topCandidateScore = 0;
    refinementLimit = 0;
    frozeScoreUpdates = false;

    std::vector<size_t> candidateSymbols;
    for (const size_t symbol : sortedSymbols(transitionStateSupport)) {
      if (!trackedStates.contains(symbol)) {
        candidateSymbols.push_back(symbol);
      }
    }
    candidateCount = candidateSymbols.size();
    const size_t defaultRefinementLimit = defaultBoundedRefinementLimit(
        candidateCount, transitionStateSupport.size(),
        focusedTransitionProjection);
    refinementLimit = defaultRefinementLimit;
    if (candidateCount <= defaultRefinementLimit) {
      return {};
    }

    if (!scoreUpdatesFrozen_) {
      scoreNewTrackedTargets(trackedStates);
    }

    std::vector<ProjectionRefinementCandidate> candidates;
    candidates.reserve(candidateSymbols.size());
    for (const size_t symbol : candidateSymbols) {
      candidates.push_back({symbol, scoreFor(symbol)});
    }
    std::sort(
        candidates.begin(),
        candidates.end(),
        projectionRefinementCandidateBetter);
    topCandidateScore = candidates.empty() ? 0 : candidates.front().score;
    if (shouldFreezeScoreUpdates(
            transitionStateSupport.size(), topCandidateScore)) {
      scoreUpdatesFrozen_ = true;
      frozeScoreUpdates = true;
    }
    refinementLimit = supportAwareRefinementLimit(
        transitionStateSupport.size(), topCandidateScore,
        defaultRefinementLimit);

    std::unordered_set<size_t> selectedSupport;
    for (const ProjectionRefinementCandidate& candidate : candidates) {
      if (selectedSupport.size() >= refinementLimit) {
        break;
      }
      selectedSupport.insert(candidate.symbol);
    }
    closeProjectionRefinementTransitionTargets(
        resolver_, complementPrimary_, selectedSupport);
    return selectedSupport;
  }

  std::unordered_set<size_t> selectModelGuidedBackfill(
      const std::unordered_set<size_t>& trackedStates,
      const std::unordered_set<size_t>& transitionStateSupport,
      const std::unordered_set<size_t>& modelGuidedSupport,
      size_t refinementLimit,
      size_t& candidateCount,
      size_t& topCandidateScore) {
    candidateCount = 0;
    topCandidateScore = 0;
    if (modelGuidedSupport.empty() || refinementLimit == 0) {
      return {}; // LCOV_EXCL_LINE
    }

    std::unordered_map<size_t, size_t> supportScore;
    for (const size_t requested : sortedSymbols(modelGuidedSupport)) {
      const size_t target =
          transitionTargetFor(requested, resolver_, complementPrimary_);
      if (!resolver_.contains(target)) {
        continue;
      }
      for (const size_t symbol : resolver_.support(target)) { // LCOV_EXCL_LINE
        if (symbol < 2 || trackedStates.contains(symbol) || // LCOV_EXCL_LINE
            !transitionStateSupport.contains(symbol)) { // LCOV_EXCL_LINE
          continue; // LCOV_EXCL_LINE
        }
        ++supportScore[symbol]; // LCOV_EXCL_LINE
      }
    }

    std::vector<ProjectionRefinementCandidate> candidates;
    candidates.reserve(supportScore.size());
    for (const auto& [symbol, score] : supportScore) {
      candidates.push_back({symbol, score}); // LCOV_EXCL_LINE
    }
    std::sort(
        candidates.begin(),
        candidates.end(),
        projectionRefinementCandidateBetter);
    candidateCount = candidates.size();
    topCandidateScore = candidates.empty() ? 0 : candidates.front().score;

    std::unordered_set<size_t> selectedSupport;
    for (const ProjectionRefinementCandidate& candidate : candidates) {
      if (selectedSupport.size() >= refinementLimit) { // LCOV_EXCL_LINE
        break; // LCOV_EXCL_LINE
      }
      selectedSupport.insert(candidate.symbol); // LCOV_EXCL_LINE
    }
    // The SAT model identifies currently missing states; adding their local
    // transition fanin keeps the strict projection step focused on that witness.
    closeProjectionRefinementTransitionTargets(
        resolver_, complementPrimary_, selectedSupport);
    return selectedSupport;
  }

 private:
  bool shouldFreezeScoreUpdates(
      size_t transitionSupportSize,
      size_t topCandidateScore) const {
    // BP-sized cones eventually reach a flat fan-in surface where every
    // remaining candidate has a tiny score.  Continuing to rescore newly
    // imported targets only materializes huge resolver support sets; freezing
    // the score map keeps the strict Craig projection heuristic bounded.
    return !scoreUpdatesFrozen_ &&
           transitionSupportSize > kCraigLowScoreBackfillSupportThreshold &&
           topCandidateScore <= kCraigLowScoreBackfillScoreLimit;
  }

  size_t supportAwareRefinementLimit(
      size_t transitionSupportSize,
      size_t topCandidateScore,
      size_t defaultRefinementLimit) const {
    if (scoreUpdatesFrozen_ &&
        shouldUseSmallLowScoreRefinementStride(
            transitionSupportSize, topCandidateScore)) {
      if (transitionSupportSize > kCraigTightLowScoreBackfillSupportThreshold) {
        // BP's first hard output enters a long flat-score plateau above 104K
        // support.  At that point each extra tracked state lives across many
        // Craig rebuilds, so use a tighter strict-refinement slice.
        return std::min( // LCOV_EXCL_LINE
            defaultRefinementLimit,
            kCraigTightLowScoreProjectionRefinementLimit);
      }
      return std::min(
          defaultRefinementLimit, kCraigLowScoreProjectionRefinementLimit);
    }
    if (transitionSupportSize > kCraigHighSupportRefinementThreshold) {
      // Once local auxiliary mining is disabled, BP-sized cones keep the scored
      // order but take bounded slices so every Craig rebuild stays below the
      // physical memory target.  BP's 84K partial-cap tail was still dominated
      // by repeated rebuilds at 1024 states, so use the focused projection
      // slice there while keeping the >90K low-score plateau tight below.
      if (topCandidateScore > kCraigLowScoreBackfillScoreLimit) {
        // The partial focused image cap exposes BP's hard tail as an 84K-state
        // cone, below the old very-high threshold but still with a strong
        // fan-in signal.  Use the wider strict slice while that signal exists
        // so the proof does not rebuild the same capped query per 256 states.
        if (transitionSupportSize <= kCraigVeryHighSupportRefinementThreshold) {
          return std::max(
              defaultRefinementLimit,
              kCraigHighSupportProjectionRefinementLimit);
        }
        return std::max(
            defaultRefinementLimit,
            kCraigVeryHighSupportHighScoreProjectionRefinementLimit);
      }
      if (transitionSupportSize > kCraigVeryHighSupportRefinementThreshold) {
        return std::min( // LCOV_EXCL_LINE
            defaultRefinementLimit,
            kCraigVeryHighSupportProjectionRefinementLimit);
      }
      return std::max(
          defaultRefinementLimit,
          kCraigHighSupportProjectionRefinementLimit);
    }
    return defaultRefinementLimit;
  }

  size_t defaultBoundedRefinementLimit(
      size_t candidateCount,
      size_t transitionSupportSize,
      bool focusedTransitionProjection) const {
    return craigBoundedProjectionRefinementLimit(
        candidateCount, transitionSupportSize, focusedTransitionProjection);
  }

  void scoreNewTrackedTargets(
      const std::unordered_set<size_t>& trackedStates) {
    // Projection only grows.  Remember already-scored requested states so each
    // bounded refinement round pays for support fan-in of newly imported
    // targets, not every target imported by earlier rounds.
    for (const size_t requested : sortedSymbols(trackedStates)) {
      const size_t target =
          transitionTargetFor(requested, resolver_, complementPrimary_);
      if (!trackedStates.contains(target) || !resolver_.contains(target)) {
        continue;
      }
      if (!scoredRequestedStates_.insert(requested).second) {
        continue;
      }
      for (const size_t symbol : resolver_.support(target)) {
        if (symbol >= 2) {
          ++faninScore_[symbol];
        }
      }
    }
  }

  size_t scoreFor(size_t symbol) const {
    const auto score = faninScore_.find(symbol);
    return score == faninScore_.end() ? 0 : score->second;
  }

  const TransitionExprResolver& resolver_;
  const std::unordered_map<size_t, size_t>& complementPrimary_;
  std::unordered_set<size_t> scoredRequestedStates_;
  std::unordered_map<size_t, size_t> faninScore_;
  bool scoreUpdatesFrozen_ = false;
};

std::unordered_set<size_t> boundedProjectionRefinementSupport(
    ProjectionRefinementScorer& scorer,
    const std::unordered_set<size_t>& trackedStates,
    const std::unordered_set<size_t>& transitionStateSupport,
    bool focusedTransitionProjection,
    size_t& candidateCount,
    size_t& topCandidateScore,
    size_t& refinementLimit,
    bool& frozeScoreUpdates) {
  return scorer.selectBoundedSupport(
      trackedStates,
      transitionStateSupport,
      focusedTransitionProjection,
      candidateCount,
      topCandidateScore,
      refinementLimit,
      frozeScoreUpdates);
}

size_t frontierRegionClauseCount(const FrontierResult& frontier) {
  return frontier.region.has_value() ? regionClauseCount(*frontier.region) : 0;
}

size_t frontierRegionLiteralCount(const FrontierResult& frontier) {
  return frontier.region.has_value() ? regionLiteralCount(*frontier.region) : 0;
}

size_t frontierRegionAuxiliaryCount(const FrontierResult& frontier) {
  return frontier.region.has_value() ? frontier.region->auxiliaryCount : 0;
}

bool isLargeFocusedLookaheadAdvanceProof(
    size_t interpolantClauses,
    size_t interpolantLiterals,
    size_t interpolantAuxiliaries) {
  return interpolantClauses >= kCraigFocusedLookaheadAdvanceMinClauses ||
         interpolantLiterals >= kCraigFocusedLookaheadAdvanceMinLiterals ||
         interpolantAuxiliaries >=
             kCraigFocusedLookaheadAdvanceMinAuxiliaries;
}

bool craigGrowthBudgetExceeded(
    const CraigImcGrowthBudget& budget,
    size_t qExpansionPass,
    const FrontierResult& frontier,
    const char** reason,
    size_t qExpansionPassLimitOverride = 0) {
  if (!budget.enabled) {
    return false;
  }
  const size_t qExpansionPassLimit =
      qExpansionPassLimitOverride > 0 ? qExpansionPassLimitOverride
                                      : budget.maxQExpansionPass;
  if (qExpansionPassLimit > 0 && qExpansionPass >= qExpansionPassLimit) {
    *reason = "q_pass";
    return true;
  }
  if (budget.maxInterpolantClauses > 0 &&
      frontierRegionClauseCount(frontier) > budget.maxInterpolantClauses) {
    *reason = "clauses"; // LCOV_EXCL_LINE
    return true; // LCOV_EXCL_LINE
  }
  if (budget.maxInterpolantLiterals > 0 &&
      frontierRegionLiteralCount(frontier) > budget.maxInterpolantLiterals) {
    *reason = "literals"; // LCOV_EXCL_LINE
    return true; // LCOV_EXCL_LINE
  }
  if (budget.maxInterpolantAuxiliaries > 0 &&
      frontierRegionAuxiliaryCount(frontier) >
          budget.maxInterpolantAuxiliaries) {
    *reason = "auxiliaries"; // LCOV_EXCL_LINE
    return true; // LCOV_EXCL_LINE
  }
  return false;
}

void emitCraigGrowthBudgetExceeded(
    size_t lookahead,
    size_t qExpansionPass,
    const FrontierResult& frontier,
    const char* reason) {
  // A budget trip is not a proof result. It only tells the large dual-rail
  // caller to bisect the current output batch and retry strict Craig IMC on
  // narrower bad predicates.
  emitSecDiag(
      "SEC diag: imc Craig growth budget exceeded reason=", reason,
      " lookahead=", lookahead,
      " q_pass=", qExpansionPass,
      " solve_ms=", frontier.solveElapsedMilliseconds,
      " interpolant_ms=", frontier.interpolationElapsedMilliseconds,
      " clauses=", frontierRegionClauseCount(frontier),
      " literals=", frontierRegionLiteralCount(frontier),
      " auxiliaries=", frontierRegionAuxiliaryCount(frontier));
}

bool craigProjectionBudgetExceeded(
    const CraigImcGrowthBudget& budget,
    size_t projectionStates) {
  return budget.enabled && budget.maxProjectionStates > 0 &&
         projectionStates > budget.maxProjectionStates;
}

bool craigTransitionBuildBudgetExceeded(
    const CraigImcGrowthBudget& budget,
    const ProjectedTransitionBuildEstimate& estimate) {
  if (!budget.enabled) {
    return false;
  }
  const size_t supportLimit = kCraigProjectedTransitionBuildSupportLimit;
  return (budget.maxImageTransitionStates > 0 &&
          estimate.stateSupport > budget.maxImageTransitionStates) ||
         estimate.stateSupport > supportLimit ||
         estimate.encodedTargets > kCraigProjectedTransitionBuildTargetLimit ||
         estimate.expressionNodes > kCraigProjectedTransitionBuildNodeLimit;
}

bool shouldRefineFocusedProjectionAfterGrowthBudget(
    const FrontierResult& frontier,
    const std::unordered_set<size_t>& trackedStates) {
  if (!frontier.usesFocusedTransitionProjection()) {
    return false;
  }
  for (const size_t symbol : frontier.transitionStateSupport) {
    if (!trackedStates.contains(symbol)) {
      return true; // LCOV_EXCL_LINE
    }
  }
  return false;
}

bool shouldContinueSaturatedFocusedQExpansion(
    const FrontierResult& frontier,
    const std::unordered_set<size_t>& trackedStates,
    const char* budgetReason,
    size_t qExpansionPass,
    size_t qExpansionPassLimit) {
  if (std::strcmp(budgetReason, "q_pass") != 0 ||
      qExpansionPass >= qExpansionPassLimit) {
    return false; // LCOV_EXCL_LINE
  }
  if (isLargeFocusedLookaheadAdvanceProof(
          frontierRegionClauseCount(frontier),
          frontierRegionLiteralCount(frontier),
          frontierRegionAuxiliaryCount(frontier))) {
    return false; // LCOV_EXCL_LINE
  }
  // Once the focused transition support is fully tracked, projection
  // refinement cannot make the Craig image stronger.  Spend a bounded number of
  // extra McMillan Q-expansion passes before reporting the batch as budgeted.
  return frontier.usesFocusedTransitionProjection() &&
         !shouldRefineFocusedProjectionAfterGrowthBudget( // LCOV_EXCL_LINE
             frontier, trackedStates); // LCOV_EXCL_LINE
}

bool shouldAdvanceLookaheadAfterSaturatedFocusedQBudget(
    const FrontierResult& frontier,
    const std::unordered_set<size_t>& trackedStates,
    const char* budgetReason,
    size_t qExpansionPass,
    size_t lookahead,
    size_t maxLookahead,
    const CraigImcGrowthBudget& growthBudget,
    size_t qExpansionPassLimit) {
  return shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      frontier.usesFocusedTransitionProjection(),
      shouldRefineFocusedProjectionAfterGrowthBudget(frontier, trackedStates),
      budgetReason,
      qExpansionPass,
      lookahead,
      maxLookahead,
      growthBudget,
      frontierRegionClauseCount(frontier),
      frontierRegionLiteralCount(frontier),
      frontierRegionAuxiliaryCount(frontier),
      qExpansionPassLimit);
}

bool shouldAdvanceLookaheadAfterBudgetedFocusedSat(
    const FrontierResult& frontier,
    const std::unordered_set<size_t>& trackedStates,
    const char* budgetReason,
    size_t lookahead,
    size_t maxLookahead,
    const CraigImcGrowthBudget& growthBudget) {
  return shouldAdvanceCraigLookaheadAfterBudgetedFocusedSat(
      frontier.usesFocusedTransitionProjection(),
      shouldRefineFocusedProjectionAfterGrowthBudget(frontier, trackedStates),
      budgetReason,
      lookahead,
      maxLookahead,
      growthBudget);
}

void emitCraigProjectionBudgetExceeded(
    size_t projectionRound,
    size_t projectionStates,
    const CraigImcGrowthBudget& budget) {
  // Projection growth is another Craig proof-size budget. For a single hard
  // output there is nothing left to split, so the caller reports strict IMC as
  // inconclusive instead of letting the process run into the memory cap.
  emitSecDiag(
      "SEC diag: imc Craig growth budget exceeded reason=projection_states",
      " projection_round=", projectionRound,
      " states=", projectionStates,
      " state_limit=", budget.maxProjectionStates);
}

void emitCraigTransitionBuildBudgetExceeded(
    size_t lookahead,
    size_t qExpansionPass,
    size_t trackedStates,
    const ProjectedTransitionBuildEstimate& estimate,
    size_t largestTransitionRequestCount,
    const CraigImcGrowthBudget& budget) {
  const size_t supportLimit = kCraigProjectedTransitionBuildSupportLimit;
  const size_t activeSupportLimit =
      budget.maxImageTransitionStates > 0
          ? std::min(supportLimit, budget.maxImageTransitionStates)
          : supportLimit;
  emitSecDiag(
      "SEC diag: imc Craig growth budget exceeded "
      "reason=transition_build",
      " lookahead=", lookahead,
      " q_pass=", qExpansionPass,
      " tracked_states=", trackedStates,
      " transition_support=", estimate.stateSupport,
      " transition_targets=", estimate.encodedTargets,
      " transition_nodes=", estimate.expressionNodes,
      " largest_transition_requests=", largestTransitionRequestCount,
      " support_limit=", activeSupportLimit,
      " target_limit=", kCraigProjectedTransitionBuildTargetLimit,
      " node_limit=", kCraigProjectedTransitionBuildNodeLimit);
}

void addStateAssignments(
    SATSolverWrapper& solver,
    const std::unordered_map<size_t, int>& leaves,
    const std::vector<std::pair<size_t, bool>>& assignments) {
  solver.setCraigClausePartition(ClausePartition::A);
  for (const auto& [symbol, value] : assignments) {
    if (const auto literal = leaves.find(symbol); literal != leaves.end()) {
      solver.addClause({value ? literal->second : -literal->second});
    }
  }
}

void addConcreteCubeOrRegionRoots(
    SATSolverWrapper& solver,
    const std::unordered_map<size_t, int>& leaves,
    const std::vector<std::pair<size_t, bool>>& cubeAssignments,
    const std::vector<int>& regionRoots) {
  solver.setCraigClausePartition(ClausePartition::A);
  for (const auto& [symbol, value] : cubeAssignments) {
    const auto leaf = leaves.find(symbol);
    if (leaf == leaves.end()) {
      continue;
    }
    std::vector<int> clause = regionRoots;
    clause.push_back(value ? leaf->second : -leaf->second);
    solver.addClause(clause);
  }
}

FrontierResult deriveBoundedFrontierRegion(
    const KInductionProblem& problem,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& trackedStates,
    const AuxiliaryStateInvariants& auxiliaryInvariants,
    size_t proofDepth) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.enableCraigInterpolation();
  configureCraigProjectionSolver(solver, trackedStates.size(), "bounded");
  const bool startsAtConcreteBootstrapFrontier =
      problem.resetBootstrapCycles != 0 &&
      !problem.bootstrapStateAssignments.empty();
  // Concrete bootstrap assignments are independently derived for each design
  // by reset simulation. Starting from that post-reset frontier avoids
  // rebuilding the reset transition cone and introduces no cross-design state
  // relation.
  const size_t bootstrapDepth =
      startsAtConcreteBootstrapFrontier ? 0 : problem.resetBootstrapCycles;
  const size_t depth = bootstrapDepth + proofDepth;
  std::vector<std::unordered_map<size_t, int>> frameLits(depth + 1);
  frameLits[depth] = allocateLeafLits(
      solver, sortedSymbols(trackedStates), VariablePartition::Global);
  std::unordered_map<int, size_t> stateByVariable;
  for (const auto& [symbol, literal] : frameLits[depth]) {
    stateByVariable.emplace(std::abs(literal), symbol);
  }

  FrontierResult result;
  const auto buildStart = SteadyClock::now();
  emitSecDiag(
      "SEC diag: imc Craig bounded build begin depth=", proofDepth,
      " tracked_states=", trackedStates.size(),
      " total_depth=", depth);
  if (depth == 0 && !startsAtConcreteBootstrapFrontier) {
    addInitialFrontierConstraint(solver, problem, frameLits[0]);
  }

  std::unordered_set<size_t> required = trackedStates;
  for (size_t reverseFrame = depth; reverseFrame > 0; --reverseFrame) {
    const size_t frame = reverseFrame - 1;
    frameLits[frame] = allocateLeafLits(
        solver, problem.inputSymbols, VariablePartition::ALocal);
    FrameFormulaEncoder encoder(
        solver, std::move(frameLits[frame]), /*createMissingLeaves=*/true);
    std::unordered_set<size_t> encodedTargets;
    for (const size_t requested : required) {
      const size_t target =
          transitionTargetFor(requested, resolver, complementPrimary);
      if (!resolver.contains(target) || !encodedTargets.insert(target).second) {
        continue; // LCOV_EXCL_LINE
      }
      const auto next = frameLits[frame + 1].find(target);
      if (next == frameLits[frame + 1].end()) {
        continue; // LCOV_EXCL_LINE
      }
      solver.setCraigVariablePartition(VariablePartition::ALocal);
      const int transition = encoder.encode(resolver.at(target));
      addLiteralEquivalenceForPartition(
          solver, next->second, transition, ClausePartition::A);
    }
    frameLits[frame] = encoder.leafLits();
    addResetInputValue(
        solver,
        problem,
        frameLits[frame],
        /*asserted=*/frame < bootstrapDepth,
        VariablePartition::ALocal,
        ClausePartition::A);
    addStateSemantics(
        solver, problem, frameLits[frame], ClausePartition::A);
    addAuxiliaryStateInvariants(
        solver, frameLits[frame], auxiliaryInvariants, ClausePartition::A);

    required.clear();
    const auto states = stateSymbolSet(problem);
    for (const auto& [symbol, literal] : frameLits[frame]) {
      (void)literal;
      if (states.contains(symbol)) {
        required.insert(symbol); // LCOV_EXCL_LINE
      } // LCOV_EXCL_LINE
    }
    closeSameDesignStateSemantics(problem, required);
    if (frame > 0) {
      for (const size_t symbol : required) { // LCOV_EXCL_LINE
        if (!frameLits[frame].contains(symbol)) { // LCOV_EXCL_LINE
          solver.setCraigVariablePartition(VariablePartition::ALocal); // LCOV_EXCL_LINE
          frameLits[frame].emplace(symbol, solver.newVar() + 2); // LCOV_EXCL_LINE
        } // LCOV_EXCL_LINE
      }
    } // LCOV_EXCL_LINE
  }
  emitSecDiag(
      "SEC diag: imc Craig bounded build after_unroll depth=", proofDepth,
      " elapsed_ms=", elapsedMilliseconds(buildStart));

  if (!startsAtConcreteBootstrapFrontier) {
    addStateAssignments(
        solver, frameLits[0], problem.initialStateAssignments);
  }
  if (bootstrapDepth <= depth) {
    addStateAssignments(
        solver,
        frameLits[bootstrapDepth],
        problem.bootstrapStateAssignments);
  }
  if (startsAtConcreteBootstrapFrontier && proofDepth == 0) {
    emitSecDiag(
        "SEC diag: imc Craig starts at concrete post-reset frontier "
        "assignments=",
        problem.bootstrapStateAssignments.size(),
        " reset_cycles=", problem.resetBootstrapCycles);
  }
  addStateSemantics(
      solver, problem, frameLits[depth], ClausePartition::A);
  addAuxiliaryStateInvariants(
      solver, frameLits[depth], auxiliaryInvariants, ClausePartition::A);
  const auto states = stateSymbolSet(problem);
  for (const auto& leaves : frameLits) {
    for (const auto& [symbol, literal] : leaves) {
      (void)literal;
      if (states.contains(symbol)) {
        result.transitionStateSupport.insert(symbol);
      }
    }
  }
  closeSameDesignStateSemantics(problem, result.transitionStateSupport);
  addBadFormula(solver, problem, frameLits[depth]);
  emitSecDiag(
      "SEC diag: imc Craig bounded build end depth=", proofDepth,
      " transition_states=", result.transitionStateSupport.size(),
      " elapsed_ms=", elapsedMilliseconds(buildStart));

  // Keep phase timing diagnostic-only. It identifies whether a difficult IMC
  // batch is dominated by SAT or by Craig-interpolant materialization.
  emitSecDiag(
      "SEC diag: imc Craig bounded solve begin depth=", proofDepth,
      " tracked_states=", trackedStates.size(),
      " transition_states=", result.transitionStateSupport.size());
  const auto solveStart = SteadyClock::now();
  result.solveStatus = solver.solveStatus();
  result.solveElapsedMilliseconds = elapsedMilliseconds(solveStart);
  emitSecDiag(
      "SEC diag: imc Craig bounded solve end depth=", proofDepth,
      " status=", static_cast<int>(result.solveStatus),
      " elapsed_ms=", result.solveElapsedMilliseconds);
  if (result.solveStatus != SATSolverWrapper::SolveStatus::Unsat) {
    return result;
  }
  const auto interpolationStart = SteadyClock::now();
  result.region =
      convertInterpolant(solver.createCraigInterpolant(), stateByVariable);
  result.interpolationElapsedMilliseconds =
      elapsedMilliseconds(interpolationStart);
  emitSecDiag(
      "SEC diag: imc Craig interpolant built depth=", proofDepth,
      " type=", static_cast<int>(result.region->type),
      " clauses=", regionClauseCount(*result.region),
      " literals=", regionLiteralCount(*result.region),
      " auxiliaries=", result.region->auxiliaryCount,
      " elapsed_ms=", result.interpolationElapsedMilliseconds);
  return result;
}

FrontierResult deriveLookaheadFrontierRegion(
    const KInductionProblem& problem,
    const CraigProblemStaticIndex& staticIndex,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& reachableRegions,
    const std::vector<InterpolantRegion>& helperInvariantRegions,
    const AuxiliaryStateInvariants& auxiliaryInvariants,
    TransitionEncodingCache& transitionEncodingCache,
    ProjectedTransitionPlanCache& transitionPlanCache,
    const CraigImcGrowthBudget& growthBudget,
    bool sourceIncludesConcreteBootstrapCube,
    bool directConcreteCubeSourceEnabled,
    size_t lookahead,
    size_t qExpansionPass) {
  const ProjectedTransitionPlan& transitionPlan =
      projectedTransitionPlanForLookahead(
          problem,
          resolver,
          complementPrimary,
          staticIndex.states,
          trackedStates,
          helperInvariantRegions.size(),
          lookahead,
          transitionPlanCache);
  const std::vector<ProjectedTransitionFrame>& transitionFrames =
      transitionPlan.frames;
  const ProjectedTransitionBuildEstimate& estimatedTransitionBuild =
      transitionPlan.estimate;
  const size_t largestTransitionRequestCount =
      transitionPlan.largestTransitionRequestCount;
  if (craigTransitionBuildBudgetExceeded(
          growthBudget, estimatedTransitionBuild)) {
    FrontierResult result;
    result.buildBudgetExceeded = true;
    result.estimatedTransitionBuild = estimatedTransitionBuild;
    result.recordTransitionRequests(
        trackedStates.size(), transitionFrames.front().requests);
    emitCraigTransitionBuildBudgetExceeded(
        lookahead,
        qExpansionPass,
        trackedStates.size(),
        estimatedTransitionBuild,
        largestTransitionRequestCount,
        growthBudget);
    return result;
  }

  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.enableCraigInterpolation();
  configureCraigProjectionSolver(solver, trackedStates.size(), "image");
  configureCraigProjectionSolverForTransitionRequests(
      solver, largestTransitionRequestCount, "image_transition");

  const std::vector<size_t> tracked = sortedSymbols(trackedStates);
  std::vector<std::unordered_map<size_t, int>> frameLits(lookahead + 1);
  frameLits[0] = allocateLeafLits(
      solver, tracked, VariablePartition::ALocal);
  frameLits[1] = allocateLeafLits(
      solver, tracked, VariablePartition::Global);
  for (size_t frame = 2; frame <= lookahead; ++frame) {
    frameLits[frame] = allocateLeafLits(
        solver, tracked, VariablePartition::BLocal);
  }
  addRegionUnionConstraint(
      solver,
      helperInvariantRegions,
      frameLits[0],
      VariablePartition::ALocal,
      ClausePartition::A);
  addRegionUnionConstraint(
      solver,
      helperInvariantRegions,
      frameLits[1],
      VariablePartition::ALocal,
      ClausePartition::A);
  for (size_t frame = 2; frame <= lookahead; ++frame) {
    addRegionUnionConstraint(
        solver,
        helperInvariantRegions,
        frameLits[frame],
        VariablePartition::BLocal,
        ClausePartition::B);
  }
  addAuxiliaryStateInvariants(
      solver, frameLits[0], auxiliaryInvariants, ClausePartition::A);
  // Frame 1 receives the same A-side auxiliary clauses below together with
  // state semantics. Avoid adding the duplicate copy before transition build.
  for (size_t frame = 2; frame <= lookahead; ++frame) {
    addAuxiliaryStateInvariants(
        solver, frameLits[frame], auxiliaryInvariants, ClausePartition::B);
  }
  std::unordered_map<int, size_t> stateByVariable;
  for (const auto& [symbol, literal] : frameLits[1]) {
    stateByVariable.emplace(std::abs(literal), symbol);
  }

  const bool useDirectConcreteCube =
      sourceIncludesConcreteBootstrapCube &&
      imcDirectCubeSourceEnabled(directConcreteCubeSourceEnabled);
  if (useDirectConcreteCube) {
    std::vector<int> regionRoots;
    regionRoots.reserve(
        reachableRegions.empty() ? 0 : reachableRegions.size() - 1);
    for (size_t index = 1; index < reachableRegions.size(); ++index) {
      regionRoots.push_back(instantiateRegion(
          solver,
          reachableRegions[index],
          frameLits[0],
          VariablePartition::ALocal,
          ClausePartition::A));
    }
    // The concrete post-reset cube is always the first source region. Encode
    // cube OR other_regions directly instead of rebuilding an auxiliary cube
    // root and then OR'ing that root with the interpolant roots.
    addConcreteCubeOrRegionRoots(
        solver,
        frameLits[0],
        problem.bootstrapStateAssignments,
        regionRoots);
  } else {
    std::vector<int> currentRegionRoots;
    currentRegionRoots.reserve(reachableRegions.size());
    for (const auto& region : reachableRegions) {
      currentRegionRoots.push_back(instantiateRegion(
          solver,
          region,
          frameLits[0],
          VariablePartition::ALocal,
          ClausePartition::A));
    }
    solver.setCraigClausePartition(ClausePartition::A);
    solver.addClause(currentRegionRoots);
  }

  const auto buildStart = SteadyClock::now();
  emitSecDiag(
      "SEC diag: imc Craig image build begin lookahead=", lookahead,
      " q_pass=", qExpansionPass,
      " regions=", reachableRegions.size(),
      " tracked_states=", trackedStates.size(),
      " direct_cube_source=", useDirectConcreteCube ? 1 : 0);
  const std::unordered_set<size_t>& imageTransitionRequests =
      transitionFrames.front().requests;
  const TransitionEncodingResult transition = addProjectedTransition(
      solver,
      problem,
      resolver,
      staticIndex.stateSemantics,
      trackedStates,
      transitionFrames.front().targets,
      auxiliaryInvariants,
      transitionEncodingCache,
      frameLits[0],
      frameLits[1],
      VariablePartition::ALocal,
      ClausePartition::A,
      transitionFrames.front().expressionNodes);
  emitSecDiag(
      "SEC diag: imc Craig image build after_a_transition lookahead=",
      lookahead,
      " q_pass=", qExpansionPass,
      " elapsed_ms=", elapsedMilliseconds(buildStart));

  FrontierResult result;
  result.recordTransitionRequests(trackedStates.size(), imageTransitionRequests);
  for (const auto& [symbol, literal] : transition.currentLits) {
    (void)literal;
    if (staticIndex.states.contains(symbol)) {
      result.transitionStateSupport.insert(symbol);
    }
  }
  for (size_t frame = 1; frame < lookahead; ++frame) {
    const std::unordered_set<size_t>& suffixTransitionRequests =
        transitionFrames[frame].requests;
    result.recordTransitionRequests(
        trackedStates.size(), suffixTransitionRequests);
    const TransitionEncodingResult suffixTransition = addProjectedTransition(
        solver,
        problem,
        resolver,
        staticIndex.stateSemantics,
        trackedStates,
        transitionFrames[frame].targets,
        auxiliaryInvariants,
        transitionEncodingCache,
        frameLits[frame],
        frameLits[frame + 1],
        VariablePartition::BLocal,
        ClausePartition::B,
        transitionFrames[frame].expressionNodes);
    for (const auto& [symbol, literal] : suffixTransition.currentLits) {
      (void)literal;
      if (staticIndex.states.contains(symbol)) {
        result.transitionStateSupport.insert(symbol);
      }
    }
  }
  emitSecDiag(
      "SEC diag: imc Craig image build after_b_suffix lookahead=",
      lookahead,
      " q_pass=", qExpansionPass,
      " suffix_frames=", lookahead > 0 ? lookahead - 1 : 0,
      " largest_transition_requests=",
      largestTransitionRequestCount,
      " elapsed_ms=", elapsedMilliseconds(buildStart));
  closeSameDesignStateSemantics(problem, result.transitionStateSupport);
  addIndexedStateSemantics(
      solver, staticIndex.stateSemantics, frameLits[1], ClausePartition::A);
  addAuxiliaryStateInvariants(
      solver, frameLits[1], auxiliaryInvariants, ClausePartition::A);
  addIndexedStateSemantics(
      solver,
      staticIndex.stateSemantics,
      frameLits[lookahead],
      ClausePartition::B);
  addAuxiliaryStateInvariants(
      solver, frameLits[lookahead], auxiliaryInvariants, ClausePartition::B);
  // McMillan's suffix checks whether the bad set appears at any suffix frame,
  // not only at the last unrolled frame. This keeps each interpolant itself
  // outside Bad, which is required for the R' => R fixed-point test below.
  addSuffixBadFormula(
      solver, problem, frameLits, /*firstFrame=*/1, /*lastFrame=*/lookahead);
  emitSecDiag(
      "SEC diag: imc Craig image build end lookahead=", lookahead,
      " q_pass=", qExpansionPass,
      " transition_states=", result.transitionStateSupport.size(),
      " elapsed_ms=", elapsedMilliseconds(buildStart));

  emitSecDiag(
      "SEC diag: imc Craig image solve begin lookahead=", lookahead,
      " q_pass=", qExpansionPass,
      " regions=", reachableRegions.size(),
      " tracked_states=", trackedStates.size());
  const auto solveStart = SteadyClock::now();
  result.solveStatus = solver.solveStatus();
  result.solveElapsedMilliseconds = elapsedMilliseconds(solveStart);
  emitSecDiag(
      "SEC diag: imc Craig image solve end lookahead=", lookahead,
      " q_pass=", qExpansionPass,
      " status=", static_cast<int>(result.solveStatus),
      " elapsed_ms=", result.solveElapsedMilliseconds);
  if (result.solveStatus != SATSolverWrapper::SolveStatus::Unsat) {
    if (result.solveStatus == SATSolverWrapper::SolveStatus::Sat &&
        sourceIncludesConcreteBootstrapCube) {
      // A SAT result under a projected transition is often spurious because
      // untracked source states are unconstrained. Prefer the reset-known
      // state bits contradicted by this model before pulling the whole cone.
      result.modelGuidedTransitionStateSupport =
          modelGuidedBootstrapProjectionSupport(
              problem,
              solver,
              transition.currentLits,
              trackedStates,
              result.transitionStateSupport,
              resolver,
              complementPrimary,
              result.modelGuidedTransitionStateCandidates);
    }
    return result;
  }

  const auto interpolationStart = SteadyClock::now();
  result.region =
      convertInterpolant(solver.createCraigInterpolant(), stateByVariable);
  result.interpolationElapsedMilliseconds =
      elapsedMilliseconds(interpolationStart);
  emitSecDiag(
      "SEC diag: imc Craig image interpolant built lookahead=", lookahead,
      " q_pass=", qExpansionPass,
      " type=", static_cast<int>(result.region->type),
      " clauses=", regionClauseCount(*result.region),
      " literals=", regionLiteralCount(*result.region),
      " auxiliaries=", result.region->auxiliaryCount,
      " elapsed_ms=", result.interpolationElapsedMilliseconds);
  return result;
}

bool regionContainedInReachableUnionSkipping(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& reachableRegions,
    const std::vector<InterpolantRegion>& helperInvariantRegions,
    const AuxiliaryStateInvariants& auxiliaryInvariants,
    const InterpolantRegion& candidateRegion,
    std::optional<size_t> skippedReachableRegion) {
  if (candidateRegion.type == InterpolantRegion::Type::False) {
    return true; // LCOV_EXCL_LINE
  }
  if (reachableRegions.empty()) {
    return false; // LCOV_EXCL_LINE
  }

  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  auto stateLits = allocateLeafLits(
      solver, sortedSymbols(trackedStates), VariablePartition::ALocal);
  addRegionUnionConstraint(
      solver,
      helperInvariantRegions,
      stateLits,
      VariablePartition::ALocal,
      ClausePartition::A);
  addStateSemantics(solver, problem, stateLits, ClausePartition::A);
  addAuxiliaryStateInvariants(
      solver, stateLits, auxiliaryInvariants, ClausePartition::A);

  const int candidateRoot = instantiateRegion(
      solver,
      candidateRegion,
      stateLits,
      VariablePartition::ALocal,
      ClausePartition::A);
  solver.setCraigClausePartition(ClausePartition::A);
  solver.addClause({candidateRoot});
  size_t unionRegionCount = 0;
  for (size_t regionIndex = 0; regionIndex < reachableRegions.size();
       ++regionIndex) {
    if (skippedReachableRegion.has_value() &&
        *skippedReachableRegion == regionIndex) {
      continue;
    }
    const InterpolantRegion& region = reachableRegions[regionIndex];
    const int root = instantiateRegion(
        solver,
        region,
        stateLits,
        VariablePartition::ALocal,
        ClausePartition::A);
    solver.addClause({-root});
    ++unionRegionCount;
  }
  if (unionRegionCount == 0) {
    return false; // LCOV_EXCL_LINE
  }

  const auto solveStart = SteadyClock::now();
  emitSecDiag(
      "SEC diag: imc Craig fixedpoint containment begin regions=",
      unionRegionCount,
      " tracked_states=", trackedStates.size());
  const auto status = solver.solveStatus();
  emitSecDiag(
      "SEC diag: imc Craig fixedpoint containment end status=",
      static_cast<int>(status),
      " elapsed_ms=", elapsedMilliseconds(solveStart));
  return status == SATSolverWrapper::SolveStatus::Unsat;
}

bool regionContainedInReachableUnion(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& reachableRegions,
    const std::vector<InterpolantRegion>& helperInvariantRegions,
    const AuxiliaryStateInvariants& auxiliaryInvariants,
    const InterpolantRegion& candidateRegion) {
  return regionContainedInReachableUnionSkipping(
      problem,
      trackedStates,
      reachableRegions,
      helperInvariantRegions,
      auxiliaryInvariants,
      candidateRegion,
      std::nullopt);
}

size_t compactReachableRegionsImpl(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& helperInvariantRegions,
    const AuxiliaryStateInvariants& auxiliaryInvariants,
    std::vector<InterpolantRegion>& reachableRegions,
    size_t compactionStart,
    size_t candidateLimit) {
  if (reachableRegions.size() < compactionStart || candidateLimit == 0) {
    return 0;
  }

  size_t checked = 0;
  std::optional<size_t> removedRegionIndex;
  for (size_t regionIndex = 0;
       regionIndex < reachableRegions.size() &&
           checked < candidateLimit;
       ++regionIndex) {
    ++checked;
    if (regionContainedInReachableUnionSkipping(
            problem,
            trackedStates,
            reachableRegions,
            helperInvariantRegions,
            auxiliaryInvariants,
            reachableRegions[regionIndex],
            regionIndex)) {
      removedRegionIndex = regionIndex;
      break;
    }
  }
  if (!removedRegionIndex.has_value()) {
    return 0;
  }

  const size_t oldRegionCount = reachableRegions.size();
  std::vector<InterpolantRegion> compactedRegions;
  compactedRegions.reserve(oldRegionCount - 1);
  for (size_t regionIndex = 0; regionIndex < reachableRegions.size();
       ++regionIndex) {
    if (regionIndex != *removedRegionIndex) {
      compactedRegions.push_back(std::move(reachableRegions[regionIndex]));
    }
  }
  reachableRegions = std::move(compactedRegions);
  // Craig IMC's Q is a union of regions.  Removing a region already contained
  // in the remaining union preserves Q while keeping later image queries from
  // carrying stale roots indefinitely.
  emitSecDiag(
      "SEC diag: imc Craig compacted reachable regions ",
      oldRegionCount, "->", reachableRegions.size(),
      " checked=", checked);
  return 1;
}

bool craigInvariantExcludesBadInternal(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& invariantRegions,
    const AuxiliaryStateInvariants& auxiliaryStateInvariants) {
  if (trackedStates.empty() || invariantRegions.empty()) {
    return false; // LCOV_EXCL_LINE
  }

  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  auto stateLits = allocateLeafLits(
      solver, sortedSymbols(trackedStates), VariablePartition::ALocal);
  addStateSemantics(solver, problem, stateLits, ClausePartition::A);
  addAuxiliaryStateInvariants(
      solver, stateLits, auxiliaryStateInvariants, ClausePartition::A);

  std::vector<int> regionRoots;
  regionRoots.reserve(invariantRegions.size());
  for (const auto& region : invariantRegions) {
    regionRoots.push_back(instantiateRegion(
        solver,
        region,
        stateLits,
        VariablePartition::ALocal,
        ClausePartition::A));
  }
  solver.setCraigClausePartition(ClausePartition::A);
  solver.addClause(regionRoots);
  addBadFormula(solver, problem, stateLits);

  // Reusing a Craig IMC invariant is sound because it is already proven
  // reachable-frontier inductive for the concrete transition system.  This
  // query only asks whether the new top-level bad predicate intersects that
  // invariant; it does not add cross-design state assumptions.
  return solver.solveStatus() == SATSolverWrapper::SolveStatus::Unsat;
}

CraigImcResult runWithProjection(
    const KInductionProblem& problem,
    const CraigProblemStaticIndex& staticIndex,
    std::unordered_set<size_t>& trackedStates,
    const TransitionExprResolver& resolver,
    const std::unordered_map<size_t, size_t>& complementPrimary,
    ProjectionRefinementScorer& refinementScorer,
    const std::vector<InterpolantRegion>& helperInvariantRegions,
    AuxiliaryStateInvariants& projectionAuxiliaryInvariants,
    const CraigImcGrowthBudget& growthBudget,
    bool directConcreteCubeSourceEnabled,
    size_t maxLookahead) {
  FrontierResult frontier;
  AuxiliaryStateInvariants activeAuxiliaryInvariants =
      projectionAuxiliaryInvariants;
  std::optional<InterpolantRegion> initialRegion =
      buildConcreteAssignmentRegion(
          problem.bootstrapStateAssignments, trackedStates);
  const bool hasConcreteInitialCube = initialRegion.has_value();
  if (!initialRegion.has_value()) {
    frontier = deriveBoundedFrontierRegion(
        problem,
        resolver,
        complementPrimary,
        trackedStates,
        activeAuxiliaryInvariants,
        0);
    initialRegion = std::move(frontier.region);
    if (!initialRegion.has_value()) {
      const size_t oldSize = trackedStates.size();
      trackedStates.insert(
          frontier.transitionStateSupport.begin(),
          frontier.transitionStateSupport.end());
      closeSameDesignStateSemantics(problem, trackedStates);
      return {
          CraigImcStatus::NoProgress,
          trackedStates.size() == oldSize ? 1u : 0u};
    }
  } else {
    emitSecDiag(
        "SEC diag: imc Craig uses exact post-reset initial cube states=",
        trackedStates.size(),
        " reset_cycles=", problem.resetBootstrapCycles);
  }
  const InterpolantRegion concreteInitialRegion = std::move(*initialRegion);
  if (!craigInvariantExcludesBadInternal(
          problem, trackedStates, {concreteInitialRegion},
          activeAuxiliaryInvariants)) {
    return { // LCOV_EXCL_LINE
        hasConcreteInitialCube ? CraigImcStatus::CounterexampleCandidate // LCOV_EXCL_LINE
                               : CraigImcStatus::NoProgress,
        hasConcreteInitialCube ? 0u : 1u}; // LCOV_EXCL_LINE
  }
  size_t lookahead = 1;
  size_t qExpansionPass = 1;
  size_t localAuxiliaryRetryCount = 0;
  std::vector<InterpolantRegion> reachableRegions{concreteInitialRegion};
  TransitionEncodingCache transitionEncodingCache{
      transitionConstantAssignmentsForEncoding(
          problem, activeAuxiliaryInvariants)};
  ProjectedTransitionPlanCache transitionPlanCache;
  while (lookahead <= maxLookahead) {
    frontier = deriveLookaheadFrontierRegion(
        problem,
        staticIndex,
        resolver,
        complementPrimary,
        trackedStates,
        reachableRegions,
        helperInvariantRegions,
        activeAuxiliaryInvariants,
        transitionEncodingCache,
        transitionPlanCache,
        growthBudget,
        hasConcreteInitialCube,
        directConcreteCubeSourceEnabled,
        lookahead,
        qExpansionPass);
    if (!frontier.region.has_value()) {
      if (frontier.buildBudgetExceeded) {
        return {CraigImcStatus::BudgetExceeded, lookahead};
      }
      const bool skipLocalAuxiliaryMining =
          shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
              frontier.usesFocusedTransitionProjection(),
              trackedStates.size(),
              frontier.transitionStateSupport.size(),
              helperInvariantRegions.size());
      if (frontier.solveStatus == SATSolverWrapper::SolveStatus::Sat) {
        if (skipLocalAuxiliaryMining) {
          emitSecDiag( // LCOV_EXCL_LINE
              "SEC diag: imc Craig skips local auxiliary mining for retained "
              "helper tail lookahead=",
              lookahead,
              " q_pass=", qExpansionPass,
              " tracked_states=", trackedStates.size(), // LCOV_EXCL_LINE
              " transition_states=", frontier.transitionStateSupport.size(), // LCOV_EXCL_LINE
              " helper_regions=", helperInvariantRegions.size()); // LCOV_EXCL_LINE
        } else if (localAuxiliaryRetryCount < kCraigLocalAuxiliaryRetryLimit) {
          const AuxiliaryStateInvariants localInvariants =
              deriveLocalAuxiliaryStateInvariants(
                  problem,
                  resolver,
                  staticIndex.stateSemantics,
                  complementPrimary,
                  frontier.transitionStateSupport,
                  activeAuxiliaryInvariants);
          if (promoteLocalAuxiliaryInvariants(
                  activeAuxiliaryInvariants,
                  projectionAuxiliaryInvariants,
                  localInvariants)) {
            ++localAuxiliaryRetryCount;
            resetTransitionEncodingCache(
                transitionEncodingCache, problem, activeAuxiliaryInvariants);
            // A projected SAT model may only exist because support state bits
            // are unconstrained.  Retry a tiny bounded number of local slices
            // before projection refinement so BP/AES cones can harvest the next
            // useful constants without turning auxiliary mining into a proof
            // search of its own.
            continue;
          }
        } else {
          emitSecDiag( // LCOV_EXCL_LINE
              "SEC diag: imc Craig local auxiliary retry limit reached "
              "lookahead=",
              lookahead,
              " q_pass=", qExpansionPass,
              " limit=", kCraigLocalAuxiliaryRetryLimit);
        }
      }
      const size_t oldSize = trackedStates.size();
      size_t boundedCandidateCount = 0;
      size_t boundedTopCandidateScore = 0;
      size_t boundedRefinementLimit = 0;
      bool frozeBoundedScoreUpdates = false;
      const bool hasModelGuidedRefinement =
          !frontier.modelGuidedTransitionStateSupport.empty();
      const bool needsBoundedBackfill =
          modelGuidedProjectionNeedsBoundedBackfill(frontier);
      std::unordered_set<size_t> boundedRefinementSupport;
      if (!hasModelGuidedRefinement || needsBoundedBackfill) {
        boundedRefinementSupport = boundedProjectionRefinementSupport(
            refinementScorer,
            trackedStates,
            frontier.transitionStateSupport,
            frontier.usesFocusedTransitionProjection(),
            boundedCandidateCount,
            boundedTopCandidateScore,
            boundedRefinementLimit,
            frozeBoundedScoreUpdates);
        if (frozeBoundedScoreUpdates) {
          emitSecDiag(
              "SEC diag: imc Craig freezes low-score fanin scoring "
              "candidates=",
              boundedCandidateCount,
              " top_score=", boundedTopCandidateScore,
              " support=", frontier.transitionStateSupport.size(),
              " score_limit=", kCraigLowScoreBackfillScoreLimit);
        }
        if (!boundedRefinementSupport.empty() &&
            boundedRefinementLimit <
                boundedProjectionRefinementLimit(boundedCandidateCount)) {
          const bool cappedByLowScore =
              boundedRefinementLimit <= kCraigLowScoreProjectionRefinementLimit &&
              frontier.transitionStateSupport.size() >
                  kCraigLowScoreBackfillSupportThreshold &&
              boundedTopCandidateScore <= kCraigLowScoreBackfillScoreLimit;
          if (cappedByLowScore) {
            emitSecDiag(
                "SEC diag: imc Craig caps low-score bounded refinement "
                "candidates=",
                boundedCandidateCount,
                " selected_limit=", boundedRefinementLimit,
                " top_score=", boundedTopCandidateScore,
                " support=", frontier.transitionStateSupport.size(),
                " score_limit=", kCraigLowScoreBackfillScoreLimit);
          } else {
            const size_t supportLimit = // LCOV_EXCL_LINE
                boundedRefinementLimit <= // LCOV_EXCL_LINE
                        kCraigVeryHighSupportProjectionRefinementLimit
                    ? kCraigVeryHighSupportRefinementThreshold
                    : kCraigHighSupportRefinementThreshold;
            emitSecDiag( // LCOV_EXCL_LINE
                "SEC diag: imc Craig caps high-support bounded refinement "
                "candidates=",
                boundedCandidateCount,
                " selected_limit=", boundedRefinementLimit,
                " top_score=", boundedTopCandidateScore,
                " support=", frontier.transitionStateSupport.size(), // LCOV_EXCL_LINE
                " support_limit=", supportLimit);
          }
        }
        if (needsBoundedBackfill &&
            shouldKeepOnlyModelGuidedRefinement(
                frontier, boundedTopCandidateScore)) {
          // In huge BP cones the bounded fan-in score eventually collapses.
          // At that point the SAT model is the precise missing slice; importing
          // another low-score 512-state backfill only grows Craig memory.
          emitSecDiag( // LCOV_EXCL_LINE
              "SEC diag: imc Craig skips low-score bounded backfill "
              "model_selected=",
              frontier.modelGuidedTransitionStateSupport.size(), // LCOV_EXCL_LINE
              " candidates=", boundedCandidateCount,
              " top_score=", boundedTopCandidateScore,
              " support=", frontier.transitionStateSupport.size(), // LCOV_EXCL_LINE
              " score_limit=", kCraigLowScoreBackfillScoreLimit);
          boundedRefinementSupport.clear(); // LCOV_EXCL_LINE
        } else if (needsBoundedBackfill) {
          if (boundedTopCandidateScore <= kCraigLowScoreBackfillScoreLimit) {
            size_t guidedBackfillCandidateCount = 0;
            size_t guidedBackfillTopScore = 0;
            std::unordered_set<size_t> guidedBackfillSupport =
                refinementScorer.selectModelGuidedBackfill(
                    trackedStates,
                    frontier.transitionStateSupport,
                    frontier.modelGuidedTransitionStateSupport,
                    boundedRefinementLimit,
                    guidedBackfillCandidateCount,
                    guidedBackfillTopScore);
            if (!guidedBackfillSupport.empty()) {
              emitSecDiag( // LCOV_EXCL_LINE
                  "SEC diag: imc Craig model-guided bounded backfill "
                  "candidates=",
                  guidedBackfillCandidateCount,
                  " selected=", guidedBackfillSupport.size(), // LCOV_EXCL_LINE
                  " top_score=", guidedBackfillTopScore,
                  " support=", frontier.transitionStateSupport.size()); // LCOV_EXCL_LINE
              boundedRefinementSupport = std::move(guidedBackfillSupport); // LCOV_EXCL_LINE
            } // LCOV_EXCL_LINE
          }
          boundedRefinementSupport.insert(
              frontier.modelGuidedTransitionStateSupport.begin(),
              frontier.modelGuidedTransitionStateSupport.end());
        }
      }
      std::unordered_set<size_t> nearSaturatedRefinementSupport =
          nearSaturatedProjectionRemainderSupport(
              frontier.transitionStateSupport, trackedStates);
      const std::unordered_set<size_t>* refinementSupport =
          !boundedRefinementSupport.empty()
              ? &boundedRefinementSupport
              : &projectionRefinementSupport(frontier);
      if (!nearSaturatedRefinementSupport.empty()) {
        emitSecDiag(
            "SEC diag: imc Craig imports near-saturated projection remainder "
            "selected=",
            nearSaturatedRefinementSupport.size(),
            " tracked_states=", trackedStates.size(),
            " full=", frontier.transitionStateSupport.size(),
            " limit=", kCraigNearSaturatedProjectionRemainderLimit);
        refinementSupport = &nearSaturatedRefinementSupport;
      }
      if (!skipLocalAuxiliaryMining &&
          shouldTrySelectedLocalAuxiliaryInvariants(
              frontier, *refinementSupport, localAuxiliaryRetryCount)) {
        const AuxiliaryStateInvariants localInvariants =
            deriveLocalAuxiliaryStateInvariants(
                problem,
                resolver,
                staticIndex.stateSemantics,
                complementPrimary,
                *refinementSupport,
                activeAuxiliaryInvariants);
        if (promoteLocalAuxiliaryInvariants(
                activeAuxiliaryInvariants,
                projectionAuxiliaryInvariants,
                localInvariants)) {
          ++localAuxiliaryRetryCount; // LCOV_EXCL_LINE
          resetTransitionEncodingCache( // LCOV_EXCL_LINE
              transitionEncodingCache, problem, activeAuxiliaryInvariants); // LCOV_EXCL_LINE
          emitSecDiag( // LCOV_EXCL_LINE
              "SEC diag: imc Craig retries selected local auxiliary "
              "invariants selected=",
              refinementSupport->size(), // LCOV_EXCL_LINE
              " full_support=", frontier.transitionStateSupport.size(), // LCOV_EXCL_LINE
              " retry=", localAuxiliaryRetryCount,
              " retry_limit=", kCraigLocalAuxiliaryRetryLimit);
          continue; // LCOV_EXCL_LINE
        }
      }
      // Do not close refinement picks over every rail/equality partner here.
      // Missing partners weaken the projected Craig query but keep any proof
      // sound, while avoiding thousands of extra global interpolant variables.
      trackedStates.insert(
          refinementSupport->begin(),
          refinementSupport->end());
      if (trackedStates.size() != oldSize) {
        if (hasModelGuidedRefinement) {
          emitSecDiag(
              "SEC diag: imc Craig model-guided projection refinement "
              "candidates=",
              frontier.modelGuidedTransitionStateCandidates,
              " selected=",
              frontier.modelGuidedTransitionStateSupport.size(),
              " full=", frontier.transitionStateSupport.size());
        }
        if (!boundedRefinementSupport.empty()) {
          emitSecDiag(
              "SEC diag: imc Craig bounded projection refinement candidates=",
              boundedCandidateCount,
              " selected=", boundedRefinementSupport.size(),
              " full=", frontier.transitionStateSupport.size(),
              " top_score=", boundedTopCandidateScore);
        }
        emitSecDiag(
            "SEC diag: imc Craig refines transition projection states=",
            oldSize, "->", trackedStates.size());
        return {CraigImcStatus::NoProgress, 0};
      }
      if (hasConcreteInitialCube &&
          qExpansionPass == 1 &&
          frontier.solveStatus == SATSolverWrapper::SolveStatus::Sat) { // LCOV_EXCL_LINE
        // With Q == S0, SAT is a concrete bounded candidate. After Q grows,
        // SAT only means the over-approximation was too coarse.
        return { // LCOV_EXCL_LINE
            CraigImcStatus::CounterexampleCandidate,
            lookahead}; // LCOV_EXCL_LINE
      }
      const char* budgetReason = nullptr;
      if (craigGrowthBudgetExceeded(
              growthBudget, qExpansionPass, frontier, &budgetReason)) {
        if (shouldAdvanceLookaheadAfterBudgetedFocusedSat(
                frontier,
                trackedStates,
                budgetReason,
                lookahead,
                maxLookahead,
                growthBudget)) {
          // A focused SAT frontier with every transition-support bit already
          // tracked has no projection work left.  Let strict IMC perform its
          // ordinary SAT action, k := k + 1, instead of treating the q-pass
          // guard as a proof failure.
          emitSecDiag( // LCOV_EXCL_LINE
              "SEC diag: imc Craig advances focused lookahead after budgeted "
              "sat frontier lookahead=",
              lookahead,
              " next=", lookahead + 1, // LCOV_EXCL_LINE
              " q_pass=", qExpansionPass,
              " tracked_states=", trackedStates.size(), // LCOV_EXCL_LINE
              " solve_ms=", frontier.solveElapsedMilliseconds); // LCOV_EXCL_LINE
          ++lookahead; // LCOV_EXCL_LINE
          qExpansionPass = 1; // LCOV_EXCL_LINE
          localAuxiliaryRetryCount = 0; // LCOV_EXCL_LINE
          reachableRegions = {concreteInitialRegion}; // LCOV_EXCL_LINE
          continue; // LCOV_EXCL_LINE
        } else {
          emitCraigGrowthBudgetExceeded(
              lookahead, qExpansionPass, frontier, budgetReason);
          return {CraigImcStatus::BudgetExceeded, lookahead};
        }
      }
      // McMillan SAT branch: increase k and restart from Q := S0.
      ++lookahead;
      qExpansionPass = 1;
      localAuxiliaryRetryCount = 0;
      reachableRegions = {concreteInitialRegion};
      continue;
    }
    const char* budgetReason = nullptr;
    const size_t projectionQExpansionPassLimit =
        craigFocusedProjectionRefinementQExpansionPassLimit(
            frontier.usesFocusedTransitionProjection(),
            trackedStates.size(),
            frontier.transitionStateSupport.size(),
            helperInvariantRegions.size());
    const bool budgetExceeded = craigGrowthBudgetExceeded(
        growthBudget,
        qExpansionPass,
        frontier,
        &budgetReason,
        projectionQExpansionPassLimit);
    const InterpolantRegion& nextRegion = *frontier.region;
    if (regionContainedInReachableUnion(
            problem,
            trackedStates,
            reachableRegions,
            helperInvariantRegions,
            activeAuxiliaryInvariants,
            nextRegion)) {
      return makeCraigEquivalentResult(
          lookahead,
          std::move(reachableRegions),
          trackedStates,
          activeAuxiliaryInvariants);
    }
    if (budgetExceeded) {
      // The Craig proof has already paid to build this interpolant.  Let a large
      // final interpolant close the strict IMC fixed point, but stop before
      // adding it to Q when it only grows the next proof obligation.  For a
      // focused BP/AES image, use the proof's transition support as the next
      // strict projection slice instead of spending the budget on deeper suffix
      // frames over the same projected state surface.
      if (shouldRefineFocusedProjectionAfterGrowthBudget(
              frontier, trackedStates)) {
        const size_t oldSize = trackedStates.size(); // LCOV_EXCL_LINE
        size_t boundedCandidateCount = 0; // LCOV_EXCL_LINE
        size_t boundedTopCandidateScore = 0; // LCOV_EXCL_LINE
        size_t boundedRefinementLimit = 0; // LCOV_EXCL_LINE
        bool frozeBoundedScoreUpdates = false; // LCOV_EXCL_LINE
        std::unordered_set<size_t> boundedRefinementSupport =
            boundedProjectionRefinementSupport( // LCOV_EXCL_LINE
                refinementScorer, // LCOV_EXCL_LINE
                trackedStates, // LCOV_EXCL_LINE
                frontier.transitionStateSupport, // LCOV_EXCL_LINE
                /*focusedTransitionProjection=*/true,
                boundedCandidateCount,
                boundedTopCandidateScore,
                boundedRefinementLimit,
                frozeBoundedScoreUpdates);
        const std::unordered_set<size_t>& refinementSupport = // LCOV_EXCL_LINE
            boundedRefinementSupport.empty() ? frontier.transitionStateSupport // LCOV_EXCL_LINE
                                             : boundedRefinementSupport;
        trackedStates.insert( // LCOV_EXCL_LINE
            refinementSupport.begin(), refinementSupport.end()); // LCOV_EXCL_LINE
        if (trackedStates.size() != oldSize) { // LCOV_EXCL_LINE
          if (frozeBoundedScoreUpdates) { // LCOV_EXCL_LINE
            emitSecDiag( // LCOV_EXCL_LINE
                "SEC diag: imc Craig freezes low-score fanin scoring "
                "candidates=",
                boundedCandidateCount,
                " top_score=", boundedTopCandidateScore,
                " support=", frontier.transitionStateSupport.size(), // LCOV_EXCL_LINE
                " score_limit=", kCraigLowScoreBackfillScoreLimit);
          } // LCOV_EXCL_LINE
          if (!boundedRefinementSupport.empty()) { // LCOV_EXCL_LINE
            emitSecDiag( // LCOV_EXCL_LINE
                "SEC diag: imc Craig refines focused projection after growth "
                "budget reason=",
                budgetReason,
                " candidates=", boundedCandidateCount,
                " selected=", boundedRefinementSupport.size(), // LCOV_EXCL_LINE
                " full=", frontier.transitionStateSupport.size(), // LCOV_EXCL_LINE
                " top_score=", boundedTopCandidateScore);
          } // LCOV_EXCL_LINE
          emitSecDiag( // LCOV_EXCL_LINE
              "SEC diag: imc Craig refines transition projection states=",
              oldSize, "->", trackedStates.size()); // LCOV_EXCL_LINE
          return {CraigImcStatus::NoProgress, 0}; // LCOV_EXCL_LINE
        }
      } // LCOV_EXCL_LINE
      const size_t qExpansionPassLimit =
          craigFocusedSaturatedQExpansionPassLimit(
              frontier.usesFocusedTransitionProjection(),
              trackedStates.size(),
              frontier.transitionStateSupport.size(),
              helperInvariantRegions.size());
      if (shouldContinueSaturatedFocusedQExpansion(
              frontier,
              trackedStates,
              budgetReason,
              qExpansionPass,
              qExpansionPassLimit)) {
        emitSecDiag( // LCOV_EXCL_LINE
            "SEC diag: imc Craig continues saturated focused q expansion "
            "lookahead=",
            lookahead,
            " q_pass=", qExpansionPass,
            " limit=", qExpansionPassLimit,
            " tracked_states=", trackedStates.size()); // LCOV_EXCL_LINE
      } else if (shouldAdvanceLookaheadAfterSaturatedFocusedQBudget(
                     frontier,
                     trackedStates,
                     budgetReason,
                     qExpansionPass,
                     lookahead,
                     maxLookahead,
                     growthBudget,
                     qExpansionPassLimit)) {
        // The q-pass limit is an optimization guard, not a semantic bound.
        // When projection is saturated, avoid q13's proof explosion by trying
        // the next strict IMC unroll depth from the concrete frontier.
        emitSecDiag( // LCOV_EXCL_LINE
            "SEC diag: imc Craig advances focused lookahead after q budget "
            "lookahead=",
            lookahead,
            " next=", lookahead + 1, // LCOV_EXCL_LINE
            " q_pass=", qExpansionPass,
            " tracked_states=", trackedStates.size(), // LCOV_EXCL_LINE
            " clauses=", frontierRegionClauseCount(frontier), // LCOV_EXCL_LINE
            " literals=", frontierRegionLiteralCount(frontier), // LCOV_EXCL_LINE
            " auxiliaries=", frontierRegionAuxiliaryCount(frontier)); // LCOV_EXCL_LINE
        ++lookahead; // LCOV_EXCL_LINE
        qExpansionPass = 1; // LCOV_EXCL_LINE
        localAuxiliaryRetryCount = 0; // LCOV_EXCL_LINE
        reachableRegions = {concreteInitialRegion}; // LCOV_EXCL_LINE
        continue; // LCOV_EXCL_LINE
      } else {
        emitCraigGrowthBudgetExceeded(
            lookahead, qExpansionPass, frontier, budgetReason);
        return {CraigImcStatus::BudgetExceeded, lookahead};
      }
    } // LCOV_EXCL_LINE

    // McMillan UNSAT branch: Q := Q OR I, then repeat the same loop without
    // increasing k.
    reachableRegions.push_back(std::move(*frontier.region));
    compactReachableRegionsImpl(
        problem,
        trackedStates,
        helperInvariantRegions,
        activeAuxiliaryInvariants,
        reachableRegions,
        kCraigRegionCompactionStart,
        kCraigRegionCompactionCandidateLimit);
    ++qExpansionPass;
    localAuxiliaryRetryCount = 0;
  }
  return {
      hasConcreteInitialCube ? CraigImcStatus::ConcreteNoProgress
                             : CraigImcStatus::NoProgress,
      maxLookahead};
}

}  // namespace

InterpolantRegion simplifyCraigInterpolantRegion(InterpolantRegion region) {
  return simplifyCraigInterpolantRegionImpl(std::move(region));
} // LCOV_EXCL_LINE

size_t compactCraigReachableRegions(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& helperInvariantRegions,
    std::vector<InterpolantRegion>& reachableRegions,
    size_t compactionStart,
    size_t candidateLimit) {
  return compactReachableRegionsImpl(
      problem,
      trackedStates,
      helperInvariantRegions,
      /*auxiliaryInvariants=*/{},
      reachableRegions,
      compactionStart,
      candidateLimit);
} // LCOV_EXCL_LINE

bool shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
    bool focusedTransitionProjection,
    bool hasUntrackedTransitionSupport,
    const char* budgetReason,
    size_t qExpansionPass,
    size_t lookahead,
    size_t maxLookahead,
    const CraigImcGrowthBudget& budget,
    size_t interpolantClauses,
    size_t interpolantLiterals,
    size_t interpolantAuxiliaries,
    size_t qExpansionPassLimit) {
  (void)budget;
  const bool saturatedQPass =
      qExpansionPass >= qExpansionPassLimit;
  const bool largeProof = isLargeFocusedLookaheadAdvanceProof(
      interpolantClauses, interpolantLiterals, interpolantAuxiliaries);
  return budgetReason != nullptr &&
         std::strcmp(budgetReason, "q_pass") == 0 &&
         (saturatedQPass || largeProof) &&
         lookahead < maxLookahead &&
         focusedTransitionProjection &&
         !hasUntrackedTransitionSupport &&
         interpolantClauses <= kCraigFocusedLookaheadAdvanceMaxClauses &&
         interpolantLiterals <= kCraigFocusedLookaheadAdvanceMaxLiterals &&
         interpolantAuxiliaries <=
             kCraigFocusedLookaheadAdvanceMaxAuxiliaries;
}

bool shouldAdvanceCraigLookaheadAfterBudgetedFocusedSat(
    bool focusedTransitionProjection,
    bool hasUntrackedTransitionSupport,
    const char* budgetReason,
    size_t lookahead,
    size_t maxLookahead,
    const CraigImcGrowthBudget& budget) {
  (void)budget;
  return budgetReason != nullptr &&
         std::strcmp(budgetReason, "q_pass") == 0 &&
         lookahead < maxLookahead &&
         focusedTransitionProjection &&
         !hasUntrackedTransitionSupport;
}

size_t craigBoundedProjectionRefinementLimit(
    size_t candidateCount,
    size_t transitionSupportSize,
    bool focusedTransitionProjection) {
  const size_t defaultLimit = boundedProjectionRefinementLimit(candidateCount);
  if (!focusedTransitionProjection ||
      transitionSupportSize > kCraigFocusedProjectionRefinementSupportLimit) {
    return defaultLimit;
  }
  if (transitionSupportSize <= kCraigFocusedProjectionBulkSupportLimit &&
      candidateCount <= kCraigFocusedProjectionBulkCandidateLimit) {
    // Focused image queries already cap the encoded transition request.  When
    // the remaining support is BP-modest, import it in one strict refinement
    // step instead of rebuilding the same saturated Craig proof per 4K states.
    return candidateCount;
  }
  return std::max(defaultLimit, kCraigFocusedProjectionRefinementLimit);
}

bool shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
    bool focusedTransitionProjection,
    size_t trackedStateCount,
    size_t transitionSupportSize,
    size_t helperInvariantRegionCount) {
  // Helper-backed BP/AES singleton outputs need their small and mid-size local
  // auxiliary packets; those transition-proven facts can close a strict Craig
  // fixed point before the projection grows to the full suffix surface.  Skip
  // only the broad retained-helper surface where validation probes dominate.
  return focusedTransitionProjection &&
         helperInvariantRegionCount >=
             kCraigRetainedHelperLocalAuxiliarySkipRegions &&
         trackedStateCount >=
             kCraigRetainedHelperLocalAuxiliarySkipStateThreshold &&
         transitionSupportSize >=
             kCraigFocusedProjectionRefinementSupportLimit;
}

bool shouldShortenRetainedHelperQReplay(
    bool focusedTransitionProjection,
    size_t trackedStateCount,
    size_t transitionSupportSize,
    size_t helperInvariantRegionCount) {
  // The first helper-backed singleton still needs q6 to close its Craig fixed
  // point.  Shorten only later retained-helper tails, where q4+ mostly replays
  // the same saturated suffix image before each strict projection import.
  return focusedTransitionProjection &&
         helperInvariantRegionCount >=
             kCraigRetainedHelperLocalAuxiliarySkipRegions &&
         trackedStateCount >=
             kCraigRetainedHelperLocalAuxiliarySkipStateThreshold &&
         transitionSupportSize >=
             kCraigRetainedHelperLocalAuxiliarySkipStateThreshold;
}

size_t craigFocusedSaturatedQExpansionPassLimit(
    bool focusedTransitionProjection,
    size_t trackedStateCount,
    size_t transitionSupportSize,
    size_t helperInvariantRegionCount) {
  // The first singleton still needs the default q6 guard.  Large helper-backed
  // focused tails replay the same saturated proof before each strict projection
  // import, so stop that replay at q3 once the support has already reached the
  // BP/AES-sized surface.
  if (shouldShortenRetainedHelperQReplay(
          focusedTransitionProjection,
          trackedStateCount,
          transitionSupportSize,
          helperInvariantRegionCount)) {
    return kCraigRetainedHelperFocusedQExpansionPassLimit;
  }
  return kCraigFocusedSaturatedQExpansionPassLimit;
}

size_t craigFocusedProjectionRefinementQExpansionPassLimit(
    bool focusedTransitionProjection,
    size_t trackedStateCount,
    size_t transitionSupportSize,
    size_t helperInvariantRegionCount) {
  // Projection refinement only needs enough Q expansion to rank the next
  // strict state slice.  In BP's unpruned retained-helper tail q4 mostly
  // repeats the same capped lookahead query, so use q3 there; zero keeps the
  // caller's normal growth budget everywhere else.
  if (shouldShortenRetainedHelperQReplay(
          focusedTransitionProjection,
          trackedStateCount,
          transitionSupportSize,
          helperInvariantRegionCount)) {
    return kCraigRetainedHelperFocusedProjectionQExpansionPassLimit;
  }
  return 0;
}

bool shouldCapCraigFocusedImageTransitionRequests(size_t expandedRequestCount) {
  return shouldCapFocusedImageTransitionRequests(
      expandedRequestCount, kCraigFocusedImageTransitionRequestLimit);
}

size_t craigFocusedImageTransitionRequestLimit(
    size_t trackedStateCount,
    size_t helperInvariantRegionCount) {
  return focusedImageTransitionRequestLimit(
      trackedStateCount, helperInvariantRegionCount);
}

size_t cappedCraigFocusedImageTransitionRequestCount(
    size_t currentRequestCount,
    size_t expandedRequestCount) {
  return cappedFocusedImageTransitionRequestCount(
      currentRequestCount,
      expandedRequestCount,
      kCraigFocusedImageTransitionRequestLimit);
}

CraigInterpolatingModelChecker::CraigInterpolatingModelChecker(
    const KInductionProblem& problem,
    const std::vector<InterpolantRegion>* helperInvariantRegions,
    const std::unordered_set<size_t>* initialTrackedStates,
    CraigImcOptions options)
    : problem_(problem),
      helperInvariantRegions_(helperInvariantRegions),
      initialTrackedStates_(initialTrackedStates),
      options_(options) {}

CraigImcResult CraigInterpolatingModelChecker::run(
    size_t maxLookahead) const {
  const CraigProblemStaticIndex staticIndex{
      stateSymbolSet(problem_), buildStateSemanticsIndex(problem_)};
  std::unordered_set<size_t> trackedStates =
      initialTrackedStates(problem_, staticIndex.states);
  if (initialTrackedStates_ != nullptr) {
    for (const size_t symbol : *initialTrackedStates_) {
      if (staticIndex.states.contains(symbol)) {
        trackedStates.insert(symbol);
      }
    }
    closeSameDesignStateSemantics(problem_, trackedStates);
  }
  if (trackedStates.empty()) {
    return {CraigImcStatus::Equivalent, 0};
  }
  const std::vector<InterpolantRegion> emptyHelperRegions;
  const std::vector<InterpolantRegion>& helperRegions =
      helperInvariantRegions_ == nullptr ? emptyHelperRegions
                                         : *helperInvariantRegions_;
  AuxiliaryStateInvariants projectionAuxiliaryInvariants =
      deriveAuxiliaryStateInvariants(
          problem_, options_.enableAuxiliaryInvariants);
  const AuxiliaryStateInvariants helperAuxiliaryInvariants =
      helperAuxiliaryStateInvariantsFromOptions(problem_, options_);
  if (!helperAuxiliaryInvariants.empty()) {
    const size_t added = mergeAuxiliaryStateInvariants(
        projectionAuxiliaryInvariants, helperAuxiliaryInvariants);
    emitSecDiag(
        "SEC diag: imc Craig seeds helper auxiliary invariants constants=",
        helperAuxiliaryInvariants.constants.size(),
        " equalities=", helperAuxiliaryInvariants.equalities.size(),
        " added=", added);
  }
  const TransitionExprResolver resolver(problem_);
  const auto complementPrimary = primaryByComplement(problem_);
  ProjectionRefinementScorer refinementScorer(resolver, complementPrimary);

  for (size_t projectionRound = 0;
       projectionRound <= problem_.state0Symbols.size() +
                              problem_.state1Symbols.size();
       ++projectionRound) {
    const size_t projectionSize = trackedStates.size();
    if (craigProjectionBudgetExceeded(
            options_.growthBudget, projectionSize)) {
      emitCraigProjectionBudgetExceeded( // LCOV_EXCL_LINE
          projectionRound, projectionSize, options_.growthBudget); // LCOV_EXCL_LINE
      return {CraigImcStatus::BudgetExceeded, 0}; // LCOV_EXCL_LINE
    }
    emitSecDiag(
        "SEC diag: imc Craig projection round=", projectionRound,
        " states=", projectionSize);
    const CraigImcResult result =
        runWithProjection(
            problem_,
            staticIndex,
            trackedStates,
            resolver,
            complementPrimary,
            refinementScorer,
            helperRegions,
            projectionAuxiliaryInvariants,
            options_.growthBudget,
            options_.enableDirectConcreteCubeSource,
            maxLookahead);
    if (result.status == CraigImcStatus::Equivalent ||
        result.iterations != 0 ||
        trackedStates.size() == projectionSize) {
      return result;
    }
    if (craigProjectionBudgetExceeded(
            options_.growthBudget, trackedStates.size())) {
      emitCraigProjectionBudgetExceeded(
          projectionRound, trackedStates.size(), options_.growthBudget);
      return {CraigImcStatus::BudgetExceeded, 0};
    }
  }
  return {}; // LCOV_EXCL_LINE
}

bool craigInvariantExcludesBad(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& invariantRegions,
    const std::vector<std::pair<size_t, bool>>& auxiliaryStateInvariants,
    const std::vector<std::pair<size_t, size_t>>& auxiliaryStateEqualities) {
  AuxiliaryStateInvariants auxiliaryInvariants;
  auxiliaryInvariants.constants = auxiliaryStateInvariants;
  auxiliaryInvariants.equalities = auxiliaryStateEqualities;
  return craigInvariantExcludesBadInternal(
      problem, trackedStates, invariantRegions, auxiliaryInvariants);
}

std::unordered_set<size_t> computeCraigImcProjectionClosure(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& seedSupport) {
  const auto states = stateSymbolSet(problem);
  std::unordered_set<size_t> trackedStates;
  for (const size_t symbol : seedSupport) {
    if (states.contains(symbol)) {
      trackedStates.insert(symbol);
    }
  }
  closeSameDesignStateSemantics(problem, trackedStates);
  if (trackedStates.empty()) {
    return trackedStates;
  }

  const TransitionExprResolver resolver(problem);
  const auto complementPrimary = primaryByComplement(problem);
  for (size_t round = 0; round <= states.size(); ++round) {
    const size_t oldSize = trackedStates.size();
    const auto support = transitionProjectionSupport(
        problem, resolver, complementPrimary, trackedStates);
    trackedStates.insert(support.begin(), support.end());
    closeSameDesignStateSemantics(problem, trackedStates);
    if (trackedStates.size() == oldSize) {
      break;
    }
  }
  return trackedStates;
}

}  // namespace KEPLER_FORMAL::SEC
