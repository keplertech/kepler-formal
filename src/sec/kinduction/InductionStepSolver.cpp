// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "kinduction/InductionStepSolver.h"

#include "kinduction/SatEncoding.h"

namespace KEPLER_FORMAL::SEC {

namespace {

void addComplementedStateRelations(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    size_t numFrames) {
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

void addInductiveStateEqualities(SATSolverWrapper& solver,
                                 const FrameVariableStore& variables,
                                 const KInductionProblem& problem,
                                 size_t firstFrame,
                                 size_t lastFrame) {
  if (problem.inductiveStateEqualityPairs.empty() || firstFrame > lastFrame) {
    return;
  }

  for (size_t frame = firstFrame; frame <= lastFrame; ++frame) {
    for (const auto& [lhsSymbol, rhsSymbol] : problem.inductiveStateEqualityPairs) {
      addLiteralEquivalence(
          solver,
          variables.getLiteral(lhsSymbol, frame),
          variables.getLiteral(rhsSymbol, frame));
    }
  }
}

}  // namespace

bool provesByInduction(const KInductionProblem& problem,
                       KEPLER_FORMAL::Config::SolverType solverType,
                       size_t k) {
  const bool hasExplicitInductionInvariant = problem.inductionProperty != nullptr;
  BoolExpr* inductionProperty =
      hasExplicitInductionInvariant ? problem.inductionProperty : problem.property;
  BoolExpr* inductionBad =
      problem.inductionBad != nullptr ? problem.inductionBad : problem.bad;

  SATSolverWrapper solver(solverType);
  FrameVariableStore variables(solver, problem.allSymbols, k + 1);
  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs0, k + 1);
  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs1, k + 1);

  for (size_t frame = 0; frame < k; ++frame) {
    addTransitionRelation(solver, variables, problem, frame);
  }

  for (size_t frame = 0; frame < k; ++frame) {
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
    solver.addClause({encoder.encode(inductionProperty)});
  }
  if (!hasExplicitInductionInvariant) {
    addInductiveStateEqualities(solver, variables, problem, 0, k - 1);
  }

  addSimplePathConstraint(
      solver, variables, problem.combinedStateSymbols(), k + 1);

  FrameFormulaEncoder lastFrameEncoder(solver, variables.makeLeafLits(k));
  solver.addClause({lastFrameEncoder.encode(inductionBad)});
  return !solver.solve();
}

}  // namespace KEPLER_FORMAL::SEC
