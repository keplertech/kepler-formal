// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "imc/ExactInterpolantSynthesizer.h"

#include <vector>

#include "proof/ProofEngineShared.h"

namespace KEPLER_FORMAL::SEC {

// Overall exact-interpolation algorithm:
// 1. Build the one-step reachable frontier from Init and the SEC transition.
// 2. Build the one-step bad region over the next-state variables.
// 3. Enumerate all shared-state assignments exactly to synthesize an
//    interpolant between reachable and bad.
// 4. Remap that interpolant back onto the current-state symbols.
// 5. Accept it only if it is inductive for the full SEC transition system.

namespace {

BoolExpr* buildAssignmentCube(const std::vector<size_t>& symbols, size_t assignment) {
  // Enumerate one complete valuation of the shared state support.
  BoolExpr* cube = BoolExpr::createTrue();
  for (size_t bit = 0; bit < symbols.size(); ++bit) {
    BoolExpr* literal = BoolExpr::Var(symbols[bit]);
    cube = BoolExpr::And(
        cube,
        (assignment & (size_t{1} << bit)) != 0 ? literal : BoolExpr::Not(literal));
  }
  return BoolExpr::simplify(cube);
}

BoolExpr* computeExactInterpolant(
    BoolExpr* lhs,
    BoolExpr* rhs,
    const std::vector<size_t>& sharedSymbols,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (lhs == nullptr || rhs == nullptr || sharedSymbols.empty()) {
    return nullptr;
  }
  if (isProofFormulaSatisfiable(BoolExpr::And(lhs, rhs), solverType)) {
    return nullptr;
  }

  BoolExpr* interpolant = BoolExpr::createFalse();
  const size_t assignmentCount = size_t{1} << sharedSymbols.size();
  // Exact interpolation here is brute-force over the shared support: keep the
  // assignments that remain compatible with the reachable side, then verify
  // that the resulting formula really separates lhs from rhs.
  for (size_t assignment = 0; assignment < assignmentCount; ++assignment) {
    BoolExpr* cube = buildAssignmentCube(sharedSymbols, assignment);
    if (isProofFormulaSatisfiable(BoolExpr::And(lhs, cube), solverType)) {
      interpolant = BoolExpr::Or(interpolant, cube);
    }
  }

  interpolant = BoolExpr::simplify(interpolant);
  if (isProofFormulaSatisfiable(
          BoolExpr::And(lhs, BoolExpr::Not(interpolant)), solverType)) {
    return nullptr;
  }
  if (isProofFormulaSatisfiable(BoolExpr::And(interpolant, rhs), solverType)) {
    return nullptr;
  }
  return interpolant;
}

}  // namespace

ExactInterpolantSynthesizer::ExactInterpolantSynthesizer(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType)
    : problem_(problem), solverType_(solverType) {}

std::optional<BoolExpr*> ExactInterpolantSynthesizer::deriveOneStepReachableStateInvariant(
    size_t maxSharedStateBits) const {
  const std::vector<size_t> combinedStateSymbols = problem_.combinedStateSymbols();
  if (combinedStateSymbols.empty() || combinedStateSymbols.size() > maxSharedStateBits) {
    return std::nullopt;
  }

  BoolExpr* init = buildProofInitFormula(problem_);
  if (init == nullptr) {
    return std::nullopt;
  }

  // Build the one-step frontier over fresh "next-state" symbols so the shared
  // support used for interpolation is the post-state space, not the pre-state.
  size_t nextSymbol = nextFreshProofSymbol(problem_);
  const auto nextStateSymbols =
      allocateFreshProofSymbols(combinedStateSymbols, nextSymbol);
  const auto badInputSymbols =
      allocateFreshProofSymbols(problem_.inputSymbols, nextSymbol);

  BoolExpr* lhs = BoolExpr::And(
      init, buildOneStepTransitionFormula(problem_, nextStateSymbols));

  std::unordered_map<size_t, size_t> badRemap = nextStateSymbols;
  badRemap.insert(badInputSymbols.begin(), badInputSymbols.end());
  BoolExpr* rhs = remapProofFormula(problem_.bad, badRemap);

  std::vector<size_t> sharedSymbols;
  sharedSymbols.reserve(combinedStateSymbols.size());
  for (const auto symbol : combinedStateSymbols) {
    sharedSymbols.push_back(nextStateSymbols.at(symbol));
  }

  BoolExpr* interpolant =
      computeExactInterpolant(lhs, rhs, sharedSymbols, solverType_);
  if (interpolant == nullptr) {
    return std::nullopt;
  }

  // Translate the interpolant back to the current-state symbol space so it can
  // be reused as a normal SEC invariant candidate.
  std::unordered_map<size_t, size_t> restoreMap;
  restoreMap.reserve(nextStateSymbols.size());
  for (const auto& [originalSymbol, freshSymbol] : nextStateSymbols) {
    restoreMap.emplace(freshSymbol, originalSymbol);
  }

  BoolExpr* restored =
      BoolExpr::simplify(remapProofFormula(interpolant, restoreMap));
  if (!isInductiveInvariant(problem_, restored, solverType_)) {
    // Separation alone is not enough for SEC; the candidate must also be an
    // inductive invariant for the full transition relation.
    return std::nullopt;
  }
  return restored;
}

}  // namespace KEPLER_FORMAL::SEC
