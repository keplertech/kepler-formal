// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "pdr/PDREngine.h"

#include <algorithm>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

#include "proof/ProofEngineShared.h"
#include "kinduction/SatEncoding.h"

namespace KEPLER_FORMAL::SEC {

// Overall PDR algorithm:
// 1. Build Init from the SEC startup constraints and reuse any already
//    validated strengthening invariant when it is sound to do so.
// 2. Maintain frames F[0], F[1], ... where each frame stores clauses known to
//    hold for all states reachable within that many steps.
// 3. At each level, ask whether a bad state still survives the current frame.
// 4. If so, recursively search for predecessors until either Init is reached
//    (real counterexample) or the bad cube is blocked by a learned clause.
// 5. Generalize learned blocking clauses, add them to all earlier frames, and
//    then propagate them forward when the transition relation preserves them.
// 6. Stop once two adjacent frames converge, when a real bug is found, or when
//    the requested frame budget is exhausted.

namespace {

// Cubes represent a concrete bad/predecessor state, while clauses are the
// blocked generalization of such a state stored in a PDR frame.
struct CubeLiteral {
  size_t symbol = 0;
  bool value = false;

  bool operator==(const CubeLiteral& other) const {
    return symbol == other.symbol && value == other.value;
  }
};

struct ClauseLiteral {
  size_t symbol = 0;
  bool positive = false;

  bool operator==(const ClauseLiteral& other) const {
    return symbol == other.symbol && positive == other.positive;
  }
};

using StateCube = std::vector<CubeLiteral>;
using StateClause = std::vector<ClauseLiteral>;

struct FrameClauses {
  // F[i] stores clauses known to hold for all states reachable within i steps.
  std::vector<StateClause> clauses;
};

struct ProofObligation {
  // "cube is bad at level" requests either a predecessor or a blocking clause.
  StateCube cube;
  size_t level = 0;
  size_t badFrame = 0;
};

void addComplementedStateRelations(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    size_t numFrames);

void addTransitionRelation(SATSolverWrapper& solver,
                           const FrameVariableStore& variables,
                           const KInductionProblem& problem,
                           size_t frame);

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

void normalizeCube(StateCube& cube) {
  // Canonical ordering lets us compare cubes structurally and avoid learning
  // the same obligation more than once with a different literal order.
  std::sort(cube.begin(), cube.end(), [](const CubeLiteral& lhs, const CubeLiteral& rhs) {
    if (lhs.symbol != rhs.symbol) {
      return lhs.symbol < rhs.symbol;
    }
    return lhs.value < rhs.value;
  });
  cube.erase(std::unique(cube.begin(), cube.end()), cube.end());
}

void normalizeClause(StateClause& clause) {
  // Clauses are canonicalized for the same reason: later subsumption and
  // convergence checks depend on stable ordering and deduplication.
  std::sort(
      clause.begin(), clause.end(), [](const ClauseLiteral& lhs, const ClauseLiteral& rhs) {
        if (lhs.symbol != rhs.symbol) {
          return lhs.symbol < rhs.symbol;
        }
        return lhs.positive < rhs.positive;
      });
  clause.erase(std::unique(clause.begin(), clause.end()), clause.end());
}

StateClause clauseFromCube(const StateCube& cube) {
  StateClause clause;
  clause.reserve(cube.size());
  for (const auto& literal : cube) {
    clause.push_back({literal.symbol, !literal.value});
  }
  normalizeClause(clause);
  return clause;
}

StateCube cubeFromClauseNegation(const StateClause& clause) {
  StateCube cube;
  cube.reserve(clause.size());
  for (const auto& literal : clause) {
    cube.push_back({literal.symbol, !literal.positive});
  }
  normalizeCube(cube);
  return cube;
}

bool clauseSubsumes(const StateClause& lhs, const StateClause& rhs) {
  return std::includes(rhs.begin(), rhs.end(), lhs.begin(), lhs.end(),
                       [](const ClauseLiteral& a, const ClauseLiteral& b) {
                         if (a.symbol != b.symbol) {
                           return a.symbol < b.symbol;
                         }
                         return a.positive < b.positive;
                       });
}

bool frameHasSubsumingClause(const FrameClauses& frame, const StateClause& clause) {
  for (const auto& existingClause : frame.clauses) {
    if (clauseSubsumes(existingClause, clause)) {
      return true;
    }
  }
  return false;
}

void addClauseToFrame(FrameClauses& frame, StateClause clause) {
  normalizeClause(clause);
  if (frameHasSubsumingClause(frame, clause)) {
    return;
  }

  // Keep each frame minimal so later SAT queries do not carry redundant facts.
  frame.clauses.erase(
      std::remove_if(
          frame.clauses.begin(),
          frame.clauses.end(),
          [&](const StateClause& existingClause) {
            return clauseSubsumes(clause, existingClause);
          }),
      frame.clauses.end());
  frame.clauses.push_back(std::move(clause));
}

void addClauseToFrames(std::vector<FrameClauses>& frames,
                       const StateClause& clause,
                       size_t maxLevel) {
  for (size_t level = 1; level <= maxLevel; ++level) {
    addClauseToFrame(frames[level], clause);
  }
}

void addStateClause(SATSolverWrapper& solver,
                    const FrameVariableStore& variables,
                    const StateClause& clause,
                    size_t frame) {
  std::vector<int> satClause;
  satClause.reserve(clause.size());
  for (const auto& literal : clause) {
    const int satLiteral = variables.getLiteral(literal.symbol, frame);
    satClause.push_back(literal.positive ? satLiteral : -satLiteral);
  }
  solver.addClause(satClause);
}

void addCubeAssumptions(SATSolverWrapper& solver,
                        const FrameVariableStore& variables,
                        const StateCube& cube,
                        size_t frame) {
  for (const auto& literal : cube) {
    solver.addClause(
        {literal.value ? variables.getLiteral(literal.symbol, frame)
                       : -variables.getLiteral(literal.symbol, frame)});
  }
}

void addNegatedCubeClause(SATSolverWrapper& solver,
                          const FrameVariableStore& variables,
                          const StateCube& cube,
                          size_t frame) {
  std::vector<int> satClause;
  satClause.reserve(cube.size());
  for (const auto& literal : cube) {
    const int satLiteral = variables.getLiteral(literal.symbol, frame);
    satClause.push_back(literal.value ? -satLiteral : satLiteral);
  }
  solver.addClause(satClause);
}

void addFrameConstraints(SATSolverWrapper& solver,
                         const FrameVariableStore& variables,
                         BoolExpr* initFormula,
                         BoolExpr* frameInvariant,
                         const std::vector<FrameClauses>& frames,
                         size_t level,
                         size_t frame) {
  if (level == 0) {
    // F[0] is Init, so the SAT query is anchored directly in the startup
    // frontier rather than in learned blocking clauses.
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
    solver.addClause({encoder.encode(initFormula)});
    return;
  }

  // For higher frames, materialize the currently learned blocking clauses and
  // any validated strengthening invariant the strategy handed to PDR.
  for (const auto& clause : frames[level].clauses) {
    addStateClause(solver, variables, clause, frame);
  }
  if (frameInvariant != nullptr) {
    // The optional strengthening is treated exactly like a frame fact, but it
    // is validated before we allow the engine to rely on it.
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
    solver.addClause({encoder.encode(frameInvariant)});
  }
}

StateCube extractStateCube(const SATSolverWrapper& solver,
                           const FrameVariableStore& variables,
                           const std::vector<size_t>& stateSymbols,
                           size_t frame) {
  StateCube cube;
  cube.reserve(stateSymbols.size());
  for (const auto symbol : stateSymbols) {
    cube.push_back({symbol, solver.getLiteralValue(variables.getLiteral(symbol, frame))});
  }
  normalizeCube(cube);
  return cube;
}

std::optional<StateCube> findBadCube(const KInductionProblem& problem,
                                     KEPLER_FORMAL::Config::SolverType solverType,
                                     BoolExpr* initFormula,
                                     BoolExpr* frameInvariant,
                                     const std::vector<FrameClauses>& frames,
                                     size_t level) {
  // Search the current frame for a concrete state that still satisfies bad
  // after all learned blocking clauses and optional strengthening are applied.
  SATSolverWrapper solver(solverType);
  FrameVariableStore variables(solver, problem.allSymbols, 1);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs0, 1);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs1, 1);
  addFrameConstraints(
      solver, variables, initFormula, frameInvariant, frames, level, 0);
  FrameFormulaEncoder encoder(solver, variables.makeLeafLits(0));
  solver.addClause({encoder.encode(problem.bad)});
  if (!solver.solve()) {
    return std::nullopt;
  }
  return extractStateCube(solver, variables, problem.combinedStateSymbols(), 0);
}

std::optional<StateCube> findPredecessorCube(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    BoolExpr* initFormula,
    BoolExpr* frameInvariant,
    const std::vector<FrameClauses>& frames,
    size_t level,
    const StateCube& targetCube,
    bool excludeTargetOnCurrentFrame) {
  // This is the one-step predecessor query at the heart of PDR: does some
  // state in F[level] transition into the target cube on the next frame?
  SATSolverWrapper solver(solverType);
  FrameVariableStore variables(solver, problem.allSymbols, 2);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs0, 2);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs1, 2);
  addFrameConstraints(
      solver, variables, initFormula, frameInvariant, frames, level, 0);
  addTransitionRelation(solver, variables, problem, 0);
  addCubeAssumptions(solver, variables, targetCube, 1);
  if (excludeTargetOnCurrentFrame) {
    addNegatedCubeClause(solver, variables, targetCube, 0);
  }
  if (!solver.solve()) {
    return std::nullopt;
  }
  return extractStateCube(solver, variables, problem.combinedStateSymbols(), 0);
}

bool cubeIntersectsInit(const KInductionProblem& problem,
                        KEPLER_FORMAL::Config::SolverType solverType,
                        BoolExpr* initFormula,
                        const StateCube& cube) {
  // A clause is only safe to learn if its negated cube stays outside Init.
  SATSolverWrapper solver(solverType);
  FrameVariableStore variables(solver, problem.allSymbols, 1);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs0, 1);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs1, 1);
  FrameFormulaEncoder encoder(solver, variables.makeLeafLits(0));
  solver.addClause({encoder.encode(initFormula)});
  addCubeAssumptions(solver, variables, cube, 0);
  return solver.solve();
}

StateCube generalizeBlockedCube(const KInductionProblem& problem,
                                KEPLER_FORMAL::Config::SolverType solverType,
                                BoolExpr* initFormula,
                                BoolExpr* frameInvariant,
                                const std::vector<FrameClauses>& frames,
                                size_t level,
                                const StateCube& cube) {
  // Greedy clause generalization: drop literals that are unnecessary for
  // blocking and not needed to keep the cube outside Init.
  StateCube candidate = cube;
  size_t index = 0;
  while (index < candidate.size()) {
    StateCube reduced = candidate;
    reduced.erase(reduced.begin() + static_cast<std::ptrdiff_t>(index));
    if (cubeIntersectsInit(problem, solverType, initFormula, reduced)) {
      ++index;
      continue;
    }
    if (!findPredecessorCube(
             problem,
             solverType,
             initFormula,
             frameInvariant,
             frames,
             level - 1,
             reduced,
             true)
             .has_value()) {
      candidate = std::move(reduced);
      continue;
    }
    ++index;
  }
  return candidate;
}

bool framesConverged(const FrameClauses& lhs, const FrameClauses& rhs) {
  if (lhs.clauses.size() != rhs.clauses.size()) {
    return false;
  }
  for (const auto& clause : lhs.clauses) {
    if (!frameHasSubsumingClause(rhs, clause)) {
      return false;
    }
  }
  for (const auto& clause : rhs.clauses) {
    if (!frameHasSubsumingClause(lhs, clause)) {
      return false;
    }
  }
  return true;
}

bool obligationAlreadyBlocked(const std::vector<FrameClauses>& frames,
                              const ProofObligation& obligation) {
  return frameHasSubsumingClause(frames[obligation.level], clauseFromCube(obligation.cube));
}

size_t popNextObligationIndex(const std::vector<ProofObligation>& queue) {
  size_t bestIndex = 0;
  for (size_t i = 1; i < queue.size(); ++i) {
    if (queue[i].level < queue[bestIndex].level ||
        (queue[i].level == queue[bestIndex].level &&
         queue[i].badFrame < queue[bestIndex].badFrame)) {
      bestIndex = i;
    }
  }
  return bestIndex;
}

bool blockProofObligations(const KInductionProblem& problem,
                           KEPLER_FORMAL::Config::SolverType solverType,
                           BoolExpr* initFormula,
                           BoolExpr* frameInvariant,
                           std::vector<FrameClauses>& frames,
                           const StateCube& rootCube,
                           size_t rootLevel,
                           size_t& badFrame) {
  // This is the paper's recursive blocking idea expressed as an explicit queue
  // so we do not depend on deep recursion for large obligation stacks.
  std::vector<ProofObligation> queue = {{{rootCube, rootLevel, rootLevel}}};

  while (!queue.empty()) {
    const size_t obligationIndex = popNextObligationIndex(queue);
    const ProofObligation obligation = queue[obligationIndex];
    queue.erase(queue.begin() + static_cast<std::ptrdiff_t>(obligationIndex));

    if (obligationAlreadyBlocked(frames, obligation)) {
      continue;
    }

    if (obligation.level == 0) {
      badFrame = obligation.badFrame;
      return false;
    }

    if (const auto predecessor = findPredecessorCube(
            problem,
            solverType,
            initFormula,
            frameInvariant,
            frames,
            obligation.level - 1,
            obligation.cube,
            false);
        predecessor.has_value()) {
      queue.push_back(obligation);
      queue.push_back({*predecessor, obligation.level - 1, obligation.badFrame});
      continue;
    }

    // No predecessor survives F[level-1], so the cube can be blocked at every
    // frame up to "level" after a small literal-dropping generalization pass.
    const StateCube generalizedCube = generalizeBlockedCube(
        problem,
        solverType,
        initFormula,
        frameInvariant,
        frames,
        obligation.level,
        obligation.cube);
    addClauseToFrames(frames, clauseFromCube(generalizedCube), obligation.level);
    if (obligation.level < obligation.badFrame) {
      queue.push_back({generalizedCube, obligation.level + 1, obligation.badFrame});
    }
  }

  return true;
}

std::vector<StateClause> buildSeedClauses(const KInductionProblem& problem,
                                          KEPLER_FORMAL::Config::SolverType solverType,
                                          BoolExpr* initFormula) {
  std::vector<StateClause> seedClauses;
  // Seed the first learned frame with state equalities that are already
  // guaranteed by Init/bootstrap, so PDR starts from facts that are known
  // reachable-state invariants instead of rediscovering them from scratch.
  for (const auto& [lhsSymbol, rhsSymbol] : problem.inductiveStateEqualityPairs) {
    StateClause clause0 = {{lhsSymbol, false}, {rhsSymbol, true}};
    StateClause clause1 = {{lhsSymbol, true}, {rhsSymbol, false}};
    normalizeClause(clause0);
    normalizeClause(clause1);

    // Promote already-anchored state equalities into initial frame facts when
    // they are guaranteed by Init/bootstrap instead of guessed from structure.
    if (!cubeIntersectsInit(problem, solverType, initFormula, cubeFromClauseNegation(clause0))) {
      seedClauses.push_back(clause0);
    }
    if (!cubeIntersectsInit(problem, solverType, initFormula, cubeFromClauseNegation(clause1))) {
      seedClauses.push_back(clause1);
    }
  }
  return seedClauses;
}

void propagateClauses(const KInductionProblem& problem,
                      KEPLER_FORMAL::Config::SolverType solverType,
                      BoolExpr* initFormula,
                      BoolExpr* frameInvariant,
                      std::vector<FrameClauses>& frames,
                      size_t maxLevel) {
  // Standard PDR propagation: if F[i] /\ T implies a clause on the next frame,
  // move that clause forward into F[i+1].
  for (size_t level = 1; level <= maxLevel; ++level) {
    const auto snapshot = frames[level].clauses;
    for (const auto& clause : snapshot) {
      if (frameHasSubsumingClause(frames[level + 1], clause)) {
        continue;
      }
      const StateCube violatingCube = cubeFromClauseNegation(clause);
      if (!findPredecessorCube(
               problem,
               solverType,
               initFormula,
               frameInvariant,
               frames,
               level,
               violatingCube,
               false)
               .has_value()) {
        addClauseToFrame(frames[level + 1], clause);
      }
    }
  }
}

}  // namespace

PDREngine::PDREngine(const KInductionProblem& problem,
                     KEPLER_FORMAL::Config::SolverType solverType)
    : problem_(problem), solverType_(solverType) {}

PDRResult PDREngine::run(size_t maxFrames) const {
  // Build the SEC startup frontier once so every frame query shares the same
  // interpretation of reset/bootstrap and frame-0 equality constraints.
  BoolExpr* initFormula = buildProofInitFormula(problem_);
  if (initFormula == nullptr) {
    return {PDRStatus::Inconclusive, 0};
  }

  // Before entering the clause loop, try any reusable proof candidates the SEC
  // strategy already built. A candidate is only safe to inject into every PDR
  // frame once it is known to be inductive; otherwise it would over-constrain
  // the search and could hide a real counterexample.
  BoolExpr* frameInvariant = nullptr;
  const auto tryProofCandidate = [&](BoolExpr* candidate) -> std::optional<PDRResult> {
    if (candidate == nullptr ||
        !initialFrontierImplies(initFormula, candidate, solverType_) ||
        !isInductiveInvariant(problem_, candidate, solverType_)) {
      return std::nullopt;
    }

    if (invariantExcludesBadStates(problem_, candidate, solverType_)) {
      // If the strategy already synthesized a sound inductive strengthening, the
      // new engine can accept that proof immediately instead of rediscovering it.
      return PDRResult{PDRStatus::Equivalent, 1};
    }

    // Otherwise the candidate is a safe frame fact, even if it is not by
    // itself enough to prove the full SEC property.
    if (frameInvariant == nullptr) {
      frameInvariant = candidate;
    }
    return std::nullopt;
  };

  if (const auto proof = tryProofCandidate(
          selectValidatedStrengtheningInvariant(problem_, initFormula, solverType_));
      proof.has_value()) {
    return *proof;
  }

  // Keep PDR aligned with IMC: if the plain SEC property is already inductive
  // from the startup frontier, there is no reason to spend time blocking cubes.
  if (const auto proof = tryProofCandidate(problem_.property); proof.has_value()) {
    return *proof;
  }

  std::vector<FrameClauses> frames(1);

  // Before growing any frame sequence, check whether Init itself already
  // contains a bad state.
  if (auto badCube = findBadCube(
          problem_, solverType_, initFormula, frameInvariant, frames, 0);
      badCube.has_value()) {
    return {PDRStatus::Different, 0};
  }

  if (maxFrames == 0) {
    return {PDRStatus::Inconclusive, 0};
  }

  const auto seedClauses = buildSeedClauses(problem_, solverType_, initFormula);
  frames.emplace_back(FrameClauses{seedClauses});
  for (size_t level = 1; level <= maxFrames; ++level) {
    // Phase 1: exhaust the proof obligations created by bad states that still
    // survive in the current frontier.
    while (true) {
      const auto badCube =
          findBadCube(
              problem_, solverType_, initFormula, frameInvariant, frames, level);
      if (!badCube.has_value()) {
        break;
      }
      size_t badFrame = level;
      if (!blockProofObligations(
              problem_,
              solverType_,
              initFormula,
              frameInvariant,
              frames,
              *badCube,
              level,
              badFrame)) {
        return {PDRStatus::Different, badFrame};
      }
    }

    // Phase 2: create the next frame, seed it with already-known startup
    // facts, and then try to push learned clauses forward.
    frames.emplace_back(FrameClauses{seedClauses});
    propagateClauses(
        problem_, solverType_, initFormula, frameInvariant, frames, level);

    // Phase 3: convergence means F[i] == F[i+1], so the frame has become an
    // inductive invariant and the SEC property is proved.
    for (size_t convergenceLevel = 1; convergenceLevel <= level; ++convergenceLevel) {
      if (framesConverged(frames[convergenceLevel], frames[convergenceLevel + 1])) {
        return {PDRStatus::Equivalent, convergenceLevel};
      }
    }
  }

  return {PDRStatus::Inconclusive, maxFrames};
}

}  // namespace KEPLER_FORMAL::SEC
