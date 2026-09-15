// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <vector>

#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

struct RegionLiteral { // LCOV_EXCL_LINE
  // Large Craig proofs contain millions of literals. Bit fields keep each
  // literal in one word instead of padding two booleans around a size_t.
  size_t isState : 1 = false; // LCOV_EXCL_LINE
  size_t index : (sizeof(size_t) * 8 - 2) = 0; // LCOV_EXCL_LINE
  size_t positive : 1 = true; // LCOV_EXCL_LINE
};
static_assert(sizeof(RegionLiteral) == sizeof(size_t));

struct InterpolantRegion { // LCOV_EXCL_LINE
  enum class Type {
    False,
    True,
    Normal,
  };

  Type type = Type::False; // LCOV_EXCL_LINE
  size_t auxiliaryCount = 0; // LCOV_EXCL_LINE
  // Flat clause storage avoids one heap allocation and one vector object per
  // Tseitin clause while preserving the exact interpolant CNF.
  std::vector<RegionLiteral> definitionLiterals;
  std::vector<size_t> definitionClauseEnds;
  RegionLiteral root;
};

enum class CraigImcStatus {
  Equivalent,
  CounterexampleCandidate,
  ConcreteNoProgress,
  BudgetExceeded,
  NoProgress,
};

struct CraigImcGrowthBudget { // LCOV_EXCL_LINE
  bool enabled = false; // LCOV_EXCL_LINE
  size_t maxQExpansionPass = 0; // LCOV_EXCL_LINE
  size_t maxInterpolantClauses = 0; // LCOV_EXCL_LINE
  size_t maxInterpolantLiterals = 0; // LCOV_EXCL_LINE
  size_t maxInterpolantAuxiliaries = 0; // LCOV_EXCL_LINE
  size_t maxProjectionStates = 0; // LCOV_EXCL_LINE
  size_t maxImageTransitionStates = 0; // LCOV_EXCL_LINE
};

struct CraigImcOptions { // LCOV_EXCL_LINE
  // Large dual-rail IMC may derive transition-proven constants from the
  // bootstrap cube to prune Craig image queries. The environment switch still
  // enables the same path for direct checker experiments.
  bool enableAuxiliaryInvariants = false; // LCOV_EXCL_LINE
  // Avoid reifying the concrete post-reset cube as an ordinary region when the
  // caller already has an exact bootstrap assignment.
  bool enableDirectConcreteCubeSource = false; // LCOV_EXCL_LINE
  // Helper Craig invariants can depend on transition-proven auxiliary facts.
  // Seed them into the next strict IMC batch instead of re-mining the same
  // constants/equalities from scratch.
  std::vector<std::pair<size_t, bool>> helperAuxiliaryStateInvariants;
  std::vector<std::pair<size_t, size_t>> helperAuxiliaryStateEqualities;
  CraigImcGrowthBudget growthBudget;
};

struct CraigImcResult { // LCOV_EXCL_LINE
  CraigImcStatus status = CraigImcStatus::NoProgress;
  size_t iterations = 0;
  std::vector<InterpolantRegion> invariantRegions;
  std::unordered_set<size_t> trackedStates;
  std::vector<std::pair<size_t, bool>> auxiliaryStateInvariants;
  std::vector<std::pair<size_t, size_t>> auxiliaryStateEqualities;
};

// IMC stores proof interpolants in a compact CNF-like region form. This cleanup
// is intentionally IMC-local: it removes redundant interpolant clauses before
// later Craig reachability iterations re-instantiate them.
InterpolantRegion simplifyCraigInterpolantRegion(InterpolantRegion region);

// Drops at most one reachable Craig region already contained in the remaining
// union.  This preserves Q exactly while bounding later image queries.
size_t compactCraigReachableRegions(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& helperInvariantRegions,
    std::vector<InterpolantRegion>& reachableRegions,
    size_t compactionStart,
    size_t candidateLimit);

// Policy hook for the focused Craig q-pass guard.  A true result means the
// checker should restart from S0 at the next IMC lookahead instead of trying a
// larger saturated q-pass on the same projection.
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
    size_t qExpansionPassLimit = 6);

// SAT on a fully tracked focused frontier is the normal strict-IMC signal to
// increase k.  The q-pass budget is only a growth guard, so it must not turn
// that SAT step into a terminal inconclusive result.
bool shouldAdvanceCraigLookaheadAfterBudgetedFocusedSat(
    bool focusedTransitionProjection,
    bool hasUntrackedTransitionSupport,
    const char* budgetReason,
    size_t lookahead,
    size_t maxLookahead,
    const CraigImcGrowthBudget& budget);

// Projection refinement policy shared by the checker and tests. Focused Craig
// image queries encode a narrow transition request; modest focused support
// cones can be imported at once, while large BP-like cones keep a bounded
// stride to protect memory.
size_t craigBoundedProjectionRefinementLimit(
    size_t candidateCount,
    size_t transitionSupportSize,
    bool focusedTransitionProjection);

// Local auxiliary mining is optional.  Skip it when reusable helpers already
// cover a small singleton, or when a broad retained-helper tail would turn the
// validation pass into the bottleneck.
bool shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
    bool focusedTransitionProjection,
    size_t trackedStateCount,
    size_t transitionSupportSize,
    size_t helperInvariantRegionCount);

// Helper-backed singleton tails replay the same saturated q proof before the
// capped lookahead import.  This hook shortens that replay once the strict
// focused support is already BP/AES-sized.
size_t craigFocusedSaturatedQExpansionPassLimit(
    bool focusedTransitionProjection,
    size_t trackedStateCount,
    size_t transitionSupportSize,
    size_t helperInvariantRegionCount);

// Return zero when the normal growth-budget q-pass limit should be used.
size_t craigFocusedProjectionRefinementQExpansionPassLimit(
    bool focusedTransitionProjection,
    size_t trackedStateCount,
    size_t transitionSupportSize,
    size_t helperInvariantRegionCount);

// Focused multi-step images grow the requested transition slice by suffix
// preimage layers.  A true result means the next layer is large enough that the
// checker should keep the previous strict over-approximation for this query.
bool shouldCapCraigFocusedImageTransitionRequests(size_t expandedRequestCount);

size_t craigFocusedImageTransitionRequestLimit(
    size_t trackedStateCount,
    size_t helperInvariantRegionCount);

// Test hook for the focused request cap. A capped query keeps the old request
// layer and admits a deterministic prefix of the expanded layer instead of
// materializing the full BP/AES suffix preimage.
size_t cappedCraigFocusedImageTransitionRequestCount(
    size_t currentRequestCount,
    size_t expandedRequestCount);

// Large-state IMC based on proof-derived Craig interpolants. Internal state in
// the two designs remains independent; the only relational clauses retained by
// this checker are clauses emitted by CaDiCaL's UNSAT proof. A SAT result is
// reported with its concrete lookahead so the caller can reconstruct exactly
// that counterexample instead of repeating the complete bounded sweep.
class CraigInterpolatingModelChecker { // LCOV_EXCL_LINE
 public:
  explicit CraigInterpolatingModelChecker(
      const KInductionProblem& problem,
      const std::vector<InterpolantRegion>* helperInvariantRegions = nullptr,
      const std::unordered_set<size_t>* initialTrackedStates = nullptr,
      CraigImcOptions options = {});

  CraigImcResult run(size_t maxLookahead) const;

 private:
  const KInductionProblem& problem_;
  const std::vector<InterpolantRegion>* helperInvariantRegions_;
  const std::unordered_set<size_t>* initialTrackedStates_;
  CraigImcOptions options_;
};

bool craigInvariantExcludesBad(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& trackedStates,
    const std::vector<InterpolantRegion>& invariantRegions,
    const std::vector<std::pair<size_t, bool>>& auxiliaryStateInvariants = {},
    const std::vector<std::pair<size_t, size_t>>& auxiliaryStateEqualities = {});

// Computes the same state projection closure that Craig IMC will discover
// before SAT solving. This is an IMC batching aid only: it follows each
// design's own transition support and never relates internal state across
// designs.
std::unordered_set<size_t> computeCraigImcProjectionClosure(
    const KInductionProblem& problem,
    const std::unordered_set<size_t>& seedSupport);

}  // namespace KEPLER_FORMAL::SEC
