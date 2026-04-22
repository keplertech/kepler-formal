// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "kinduction/BaseCaseSolver.h"

#include <stdexcept>
#include <unordered_map>

#include "kinduction/SatEncoding.h"

namespace KEPLER_FORMAL::SEC {

namespace {

enum class InitialConstraintMode {
  None,
  ObservationOnly,
  PartialInit,
  CompleteInit,
};

size_t resetBootstrapFrames(const KInductionProblem& problem) {
  return (!problem.hasCompleteInitialState() && problem.hasResetBootstrap())
             ? problem.resetBootstrapCycles
             : 0u;
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

void addBootstrapStateEqualities(SATSolverWrapper& solver,
                                 const FrameVariableStore& variables,
                                 const KInductionProblem& problem,
                                 size_t frame) {
  for (const auto& [lhsSymbol, rhsSymbol] : problem.bootstrapStateEqualityPairs) {
    addLiteralEquivalence(
        solver,
        variables.getLiteral(lhsSymbol, frame),
        variables.getLiteral(rhsSymbol, frame));
  }
}

void addInitialStateEqualities(SATSolverWrapper& solver,
                               const FrameVariableStore& variables,
                               const KInductionProblem& problem) {
  for (const auto& [lhsSymbol, rhsSymbol] : problem.initialStateEqualityPairs) {
    addLiteralEquivalence(
        solver,
        variables.getLiteral(lhsSymbol, 0),
        variables.getLiteral(rhsSymbol, 0));
  }
}

void addBootstrapStateAssignments(SATSolverWrapper& solver,
                                  const FrameVariableStore& variables,
                                  const KInductionProblem& problem,
                                  size_t frame) {
  for (const auto& [symbol, value] : problem.bootstrapStateAssignments) {
    solver.addClause({value ? variables.getLiteral(symbol, frame)
                            : -variables.getLiteral(symbol, frame)});
  }
}

InitialConstraintMode addInitialConstraints(SATSolverWrapper& solver,
                                            const FrameVariableStore& variables,
                                            const KInductionProblem& problem) {
  if (!problem.hasSequentialState()) {
    return InitialConstraintMode::None;
  }

  const bool hasInitialStateRelation = !problem.initialStateEqualityPairs.empty();
  FrameFormulaEncoder encoder(solver, variables.makeLeafLits(0));
  if (problem.hasCompleteInitialState()) {
    solver.addClause({encoder.encode(problem.initialCondition)});
    return InitialConstraintMode::CompleteInit;
  }

  if (problem.hasExplicitInitialState()) {
    solver.addClause({encoder.encode(problem.initialCondition)});
    if (hasInitialStateRelation) {
      return InitialConstraintMode::CompleteInit;
    }
    solver.addClause({encoder.encode(problem.property)});
    return InitialConstraintMode::PartialInit;
  }

  if (hasInitialStateRelation) {
    return InitialConstraintMode::CompleteInit;
  }

  solver.addClause({encoder.encode(problem.property)});
  return InitialConstraintMode::ObservationOnly;
}

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

std::unordered_map<size_t, bool> buildFrameEnvironment(
    const SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<size_t>& symbols,
    size_t frame) {
  std::unordered_map<size_t, bool> environment;
  environment.reserve(symbols.size());
  for (const auto symbol : symbols) {
    environment.emplace(
        symbol, solver.getLiteralValue(variables.getLiteral(symbol, frame)));
  }
  return environment;
}

size_t findFirstBadFrame(const SATSolverWrapper& solver,
                         const FrameVariableStore& variables,
                         const KInductionProblem& problem,
                         size_t firstBadFrame,
                         size_t maxFrame) {
  for (size_t frame = firstBadFrame; frame <= maxFrame; ++frame) {
    if (problem.bad->evaluate(
            buildFrameEnvironment(solver, variables, problem.allSymbols, frame))) {
      return frame;
    }
  }
  throw std::runtime_error("SAT model does not satisfy any bad frame");
}

std::vector<KInductionResult::FrameInputAssignments> buildInputTrace(
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
      frameAssignments.assignments.push_back(
          {problem.environmentInputNames[i],
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
  const auto environment =
      buildFrameEnvironment(solver, variables, problem.allSymbols, frame);
  std::vector<KInductionResult::SignalMismatch> mismatches;
  for (size_t i = 0; i < problem.observedOutputExprs0.size(); ++i) {
    const bool value0 = problem.observedOutputExprs0[i]->evaluate(environment);
    const bool value1 = problem.observedOutputExprs1[i]->evaluate(environment);
    if (value0 != value1) {
      mismatches.push_back(
          {problem.observedOutputNames[i], value0, value1});
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

}  // namespace

std::optional<KInductionResult::CounterexampleWitness> findBaseCounterexample(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k) {
  const size_t bootstrapFrames = resetBootstrapFrames(problem);
  const size_t internalK = k + bootstrapFrames;
  SATSolverWrapper solver(solverType);
  FrameVariableStore variables(solver, problem.allSymbols, internalK + 1);
  addResetBootstrapConstraints(solver, variables, problem, internalK + 1);
  const InitialConstraintMode initialMode =
      bootstrapFrames == 0 ? addInitialConstraints(solver, variables, problem)
                           : InitialConstraintMode::None;

  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs0, internalK + 1);
  addComplementedStateRelations(
      solver, variables, problem.complementedStatePairs1, internalK + 1);
  addInitialStateEqualities(solver, variables, problem);

  for (size_t frame = 0; frame < internalK; ++frame) {
    addTransitionRelation(solver, variables, problem, frame);
  }
  if (bootstrapFrames != 0) {
    addBootstrapStateAssignments(solver, variables, problem, bootstrapFrames);
    addBootstrapStateEqualities(solver, variables, problem, bootstrapFrames);
  }

  size_t firstBadFrame = 0;
  if (bootstrapFrames != 0) {
    firstBadFrame = bootstrapFrames;
  } else if (initialMode == InitialConstraintMode::ObservationOnly ||
             initialMode == InitialConstraintMode::PartialInit) {
    firstBadFrame = 1;
  }
  if (firstBadFrame > internalK) {
    return std::nullopt;
  }

  std::vector<int> badClause;
  badClause.reserve(internalK - firstBadFrame + 1);
  for (size_t frame = firstBadFrame; frame <= internalK; ++frame) {
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
    badClause.push_back(encoder.encode(problem.bad));
  }
  solver.addClause(badClause);

  if (!solver.solve()) {
    return std::nullopt;
  }
  return buildCounterexampleWitness(
      solver, variables, problem, firstBadFrame, internalK, bootstrapFrames);
}

}  // namespace KEPLER_FORMAL::SEC
