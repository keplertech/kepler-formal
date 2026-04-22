// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "imc/IMCEngine.h"

#include <vector>

#include "kinduction/BaseCaseSolver.h"
#include "imc/ExactInterpolantSynthesizer.h"
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

void addComplementedStateRelations(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    size_t numFrames) {
  // Complemented outputs such as Q/QN are modeled as hard equalities between
  // the primary and inverted state views in every explored frame.
  for (size_t frame = 0; frame < numFrames; ++frame) {
    for (const auto& [primarySymbol, complementedSymbol] : complementedStatePairs) {
      addLiteralEquivalence(
          solver,
          variables.getLiteral(complementedSymbol, frame),
          -variables.getLiteral(primarySymbol, frame));
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
    return nullptr;
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
    return nullptr;
  }
  return BoolExpr::simplify(reachable);
}

}  // namespace

IMCEngine::IMCEngine(const KInductionProblem& problem,
                     KEPLER_FORMAL::Config::SolverType solverType)
    : problem_(problem), solverType_(solverType) {}

IMCResult IMCEngine::run(size_t maxK) const {
  // Keep counterexample discovery on the same bounded base-case machinery as
  // the rest of SEC so witnesses and reported cycles stay consistent.
  if (auto witness = findBaseCounterexample(problem_, solverType_, 0); witness.has_value()) {
    return {IMCStatus::Different, witness->badFrame, std::move(witness)};
  }

  if (problem_.combinedStateSymbols().empty()) {
    return {IMCStatus::Equivalent, 0};
  }

  BoolExpr* initFormula = buildProofInitFormula(problem_);
  BoolExpr* sharedStrengthening = nullptr;
  if (initFormula != nullptr) {
    // Start from any validated shared strengthening the strategy already built
    // from reset/bootstrap and aligned state structure, then let interpolation
    // sharpen it with a reachable one-step frontier when that exact derivation
    // is affordable.
    sharedStrengthening =
        selectValidatedStrengtheningInvariant(problem_, initFormula, solverType_);
    ExactInterpolantSynthesizer interpolantSynthesizer(problem_, solverType_);
    if (auto interpolant =
            interpolantSynthesizer.deriveOneStepReachableStateInvariant();
        interpolant.has_value()) {
      sharedStrengthening =
          sharedStrengthening == nullptr
              ? *interpolant
              : BoolExpr::simplify(BoolExpr::And(sharedStrengthening, *interpolant));
    }
    if (sharedStrengthening == nullptr) {
      sharedStrengthening = problem_.property;
    }

    // Before spending time on deeper frontiers, see whether the startup
    // strengthening is already a complete inductive proof.
    if (initialFrontierImplies(initFormula, sharedStrengthening, solverType_) &&
        isInductiveInvariant(problem_, sharedStrengthening, solverType_) &&
        invariantExcludesBadStates(problem_, sharedStrengthening, solverType_)) {
      return {IMCStatus::Equivalent, 1};
    }
  }

  for (size_t k = 1; k <= maxK; ++k) {
    // IMC keeps counterexample discovery and proof growth in lockstep by depth:
    // first rule out a real bug at k, then try to turn the reachable frontier
    // up to k into an inductive invariant.
    if (auto witness = findBaseCounterexample(problem_, solverType_, k);
        witness.has_value()) {
      return {IMCStatus::Different, witness->badFrame, std::move(witness)};
    }

    if (initFormula == nullptr) {
      continue;
    }

    BoolExpr* frontierInvariant =
        buildExactReachableStateInvariant(problem_, solverType_, initFormula, k);
    if (frontierInvariant == nullptr) {
      continue;
    }

    // Keep the explicit IMC engine centered on the reachable frontier, but
    // reuse any already validated strengthening to reduce the SAT work needed
    // to establish inductiveness on compact transition systems.
    BoolExpr* proofInvariant =
        sharedStrengthening == nullptr
            ? frontierInvariant
            : BoolExpr::simplify(BoolExpr::And(frontierInvariant, sharedStrengthening));

    if (initialFrontierImplies(initFormula, proofInvariant, solverType_) &&
        isInductiveInvariant(problem_, proofInvariant, solverType_) &&
        invariantExcludesBadStates(problem_, proofInvariant, solverType_)) {
      return {IMCStatus::Equivalent, k};
    }
  }

  // We exhausted the requested depth without finding either a counterexample
  // or an inductive frontier.
  return {IMCStatus::Inconclusive, maxK};
}

}  // namespace KEPLER_FORMAL::SEC
