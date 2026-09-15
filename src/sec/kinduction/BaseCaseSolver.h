// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "../../config/Config.h"
#include "kinduction/KInductionEngine.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

class TransitionExprResolver;
struct ResetFrontierReachabilityContext;

// Chooses the SAT backend used by KI/IMC/PDR one-shot base-case validation
// queries. These BMC checks do not use incremental assumptions, so they keep
// the requested proof solver; assumption-based reset/frontier helpers use
// SATSolverWrapper::assumptionSolverTypeFor at their call sites.
KEPLER_FORMAL::Config::SolverType baseCaseValidationSolverType(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType);

// Returns true when the selected base-case validation backend should use the
// short-lived local-query profile instead of the proof-oriented cone profile.
// This keeps KI base BMC probes on a direct CDCL path and gives unit tests a
// stable policy hook instead of a timing-based regression check.
bool baseCaseValidationUsesLocalQueryProfile(
    KEPLER_FORMAL::Config::SolverType solverType);

enum class BaseCounterexampleCheckStatus {
  NoCounterexample,
  Counterexample,
  Unknown,
};

struct BaseCounterexampleCheckResult {
  BaseCounterexampleCheckStatus status =
      BaseCounterexampleCheckStatus::Unknown;
  std::optional<KInductionResult::CounterexampleWitness> witness;
};

// Solves the bounded SEC base case for a single horizon k and reconstructs the
// first concrete counterexample if the property can fail within that prefix.
std::optional<KInductionResult::CounterexampleWitness> findBaseCounterexample(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k);

// Resource-bounded base proof for localized recovery paths.  A true UNSAT
// answer is required before an output may be covered; timeout stays Unknown so
// callers can conservatively split or skip the hard residual.
BaseCounterexampleCheckResult checkBaseCounterexampleWithFastValidation(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k);

// Same bounded transition prefix as findBaseCounterexample(k), but the bad
// predicate is asserted only on the newest frontier frame. K-induction calls
// this after previous smaller frontiers were already proved safe, avoiding a
// repeated OR over old bad frames on every depth.
std::optional<KInductionResult::CounterexampleWitness>
findBaseCounterexampleAtFrontier(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k);

struct ImcBaseCounterexampleCache;
using KInductionBaseCounterexampleCache = ImcBaseCounterexampleCache;

// KI/IMC-only newest-frontier query cache.  PDR keeps using the ordinary exact
// validation entry points above; this cache stores immutable COI indexes only
// for engines that sweep the same output slice across increasing depths.
std::shared_ptr<ImcBaseCounterexampleCache> makeImcBaseCounterexampleCache(
    const KInductionProblem& problem);

std::optional<KInductionResult::CounterexampleWitness>
findImcBaseCounterexampleAtFrontier(
    const ImcBaseCounterexampleCache& cache,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k);

std::shared_ptr<KInductionBaseCounterexampleCache>
makeKInductionBaseCounterexampleCache(const KInductionProblem& problem);

std::optional<KInductionResult::CounterexampleWitness>
findKInductionBaseCounterexampleAtFrontier(
    const KInductionBaseCounterexampleCache& cache,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k);

// Counterexample-search variant of the newest-frontier base query.  It keeps
// the same exact bounded-prefix semantics and witness reconstruction, but uses
// the fast SAT-validation solver profile instead of the UNSAT/proof-oriented
// profile.  Regressions enable this for cases that are expected to be different
// so KI does not spend minutes trying to prove a frontier that should produce a
// concrete witness.
std::optional<KInductionResult::CounterexampleWitness>
findFastBaseCounterexampleAtFrontier(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k);

// Exact one-sided PDR refinement helper: returns true only when an
// over-approximate, startup-pruned base query proves the frontier bad predicate
// unreachable. A false result is inconclusive, not a counterexample, because
// the pruned query may admit abstract startup states that the full SEC base
// query would reject.
bool provesNoBaseCounterexampleAtFrontier(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k);

// Checks whether a concrete reset/bootstrap prefix can reach a state cube at
// the requested post-bootstrap depth. PDR uses this as an exact CEGAR check for
// abstract startup-frontier obligations; it shares the base-case COI machinery
// so ASIC-sized memories do not force a full transition unroll.
bool isStateCubeReachableAtResetFrontier(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps);

// PDR already owns a transition resolver for its blocking loop. Reusing it here
// keeps the reset-frontier CEGAR query exact without rebuilding the large
// state-to-transition index for every blocked cube.
bool isStateCubeReachableAtResetFrontier(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    const TransitionExprResolver& transitionByState,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps);

// Builds the immutable indexing needed by repeated reset-frontier cube checks.
// PDR uses this for level-zero CEGAR refinements so thousands of tiny cubes do
// not rescan the same initial/bootstrap equality tables.
std::shared_ptr<ResetFrontierReachabilityContext>
makeResetFrontierReachabilityContext(
    const KInductionProblem& problem,
    const TransitionExprResolver& transitionByState,
    BoolExpr* frameInvariant = nullptr);

// Records an externally-proved unreachable reset-frontier cube in the same
// cache used by exact reset-frontier SAT queries. PDR uses this when a cheaper
// reset-specialized proof has already established unreachability, allowing
// later post-bootstrap queries to reuse that fact as a safe-prefix blocker.
void rememberResetFrontierUnreachableCube(
    const ResetFrontierReachabilityContext& context,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps);

bool isStateCubeReachableAtResetFrontier(
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps,
    bool usePostBootstrapPrechecks = true,
    int64_t startupConflictLimit = -1,
    int64_t startupPropagationLimit = -1);

// Prebuilds the shared assumption-capable reset-frontier solver for a cube
// support. Later exact cube queries whose symbols are covered by this support
// reuse the solver instead of rebuilding the same bounded reset prefix.
void primeResetFrontierReachabilitySolver(
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps);

// Same exact bounded-prefix check as the cached API above, but encodes the
// queried cube as unit clauses in a fresh solver. This is intentionally used
// for isolated final PDR candidate validation, where BlackParrot sampling
// showed a single incremental assumption query dominating runtime.
bool isStateCubeReachableAtResetFrontierOneShot(
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps,
    bool usePostBootstrapPrechecks = true);

// Checks the same cube against every concrete reset/bootstrap frontier from
// post-bootstrap step 0 through maxPostBootstrapSteps using one shared unroll.
// This is for PDR final-candidate validation where rebuilding the same reset
// prefix once per frame dominates runtime.
bool isStateCubeReachableWithinResetFrontier(
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t maxPostBootstrapSteps);

// Checks whether any one of several state cubes can appear at any concrete
// reset/bootstrap frontier up to maxPostBootstrapSteps. This is the batched
// form of the cube API above: it encodes one shared prefix and one disjunction
// over the requested cubes, so PDR can validate small state-only bad CNFs
// without solving one reset-frontier query per assignment.
bool anyStateCubeReachableWithinResetFrontier(
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::vector<std::pair<size_t, bool>>>& cubes,
    size_t maxPostBootstrapSteps);

// Checks whether any one of several state cubes can appear at exactly the
// requested post-bootstrap frontier. Unknown/resource-limited SAT queries are
// treated as reachable, so callers can use this as a bounded refinement proof:
// false is a sound unreachable result, true is reachable or inconclusive.
bool anyStateCubeReachableAtResetFrontier(
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::vector<std::pair<size_t, bool>>>& cubes,
    size_t postBootstrapSteps,
    long long conflictLimit,
    long long propagationLimit);

// Returns a smaller cube that is still unreachable at the requested
// reset/bootstrap frontier, when the assumption-capable reset solver can
// extract such a core. std::nullopt means the original cube is reachable.
std::optional<std::vector<std::pair<size_t, bool>>>
findResetFrontierUnreachableCubeCore(
    const ResetFrontierReachabilityContext& context,
    KEPLER_FORMAL::Config::SolverType solverType,
    const std::vector<std::pair<size_t, bool>>& cube,
    size_t postBootstrapSteps);

}  // namespace KEPLER_FORMAL::SEC
