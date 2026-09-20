// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "proof/InternalRelations.h"

#include <array>
#include <cstdio>
#include <set>
#include <unordered_map>

#include "model/SequentialDesignModel.h"
#include "kinduction/SatEncoding.h"
#include "proof/TransitionExprResolver.h"

namespace KEPLER_FORMAL::SEC {

std::vector<std::pair<size_t, size_t>> proveInternalRelations(
    const KInductionProblem& problem,
    const std::vector<InternalRelationCandidate>& candidates,
    const InternalRelationOptions& options,
    Config::SolverType solverType) {
  if (!options.learnInternalRelations || candidates.empty()) {
    return {};
  }
  const std::unordered_map<size_t, bool> initial(
      problem.initialStateAssignments.begin(), problem.initialStateAssignments.end());
  TransitionExprResolver transitions(problem);
  std::vector<const InternalRelationCandidate*> active;
  std::set<size_t> targets;
  size_t nodeBudget = 0;
  for (const auto& candidate : candidates) {
    if (active.size() >= 4096 || nodeBudget >= 250000) {
      break;
    }
    bool baseHolds = !candidate.equalities.empty();
    for (const auto& [lhs, rhs] : candidate.equalities) {
      baseHolds &= initial.contains(lhs) && initial.contains(rhs) &&
                   initial.at(lhs) == initial.at(rhs) &&
                   transitions.contains(lhs) && transitions.contains(rhs);
    }
    if (!options.allowXEqualityInInternalRelations) {
      for (const auto& value : candidate.values) {
        baseHolds &= initial.contains(value.mayBeOne) &&
                     initial.contains(value.mayBeZero) &&
                     initial.at(value.mayBeOne) != initial.at(value.mayBeZero);
      }
    }
    if (!baseHolds) {
      continue;
    }
    active.push_back(&candidate);
    for (const auto& [lhs, rhs] : candidate.equalities) {
      for (size_t symbol : {lhs, rhs}) {
        if (targets.insert(symbol).second) {
          nodeBudget += transitions.nodeCount(symbol);
        }
      }
    }
  }
  if (active.empty() || nodeBudget > 250000) {
    return {};
  }

  // Incremental assumptions let refinement remove hypotheses. Keeping a
  // rejected hypothesis as a clause would make later proofs unsound.
  SATSolverWrapper solver(SATSolverWrapper::assumptionSolverTypeFor(solverType));
  FrameVariableStore variables(solver, problem.allSymbols, 2);
  FrameFormulaEncoder current(solver, variables.makeLeafLits(0));
  FrameFormulaEncoder next(solver, variables.makeLeafLits(1));
  for (const auto& value : problem.dualRailStatePairs) {
    for (size_t frame : {0u, 1u}) {
      solver.addClause({variables.getLiteral(value.mayBeOne, frame),
                        variables.getLiteral(value.mayBeZero, frame)});
    }
  }
  for (size_t symbol : targets) {
    addLiteralEquivalence(solver, variables.getLiteral(symbol, 1),
                          current.encode(transitions.at(symbol)));
  }
  std::vector<int> hypotheses;
  std::vector<int> conclusions;
  for (const auto* candidate : active) {
    BoolExpr* relation = BoolExpr::createTrue();
    for (const auto& [lhs, rhs] : candidate->equalities) {
      relation = BoolExpr::And(relation, BoolExpr::Not(
          BoolExpr::Xor(BoolExpr::Var(lhs), BoolExpr::Var(rhs))));
    }
    if (!options.allowXEqualityInInternalRelations) {
      for (const auto& value : candidate->values) {
        relation = BoolExpr::And(relation, BoolExpr::Xor(
            BoolExpr::Var(value.mayBeOne), BoolExpr::Var(value.mayBeZero)));
      }
    }
    hypotheses.push_back(current.encode(relation));
    conclusions.push_back(next.encode(relation));
  }
  SATSolverWrapper::CadicalWorkBudget budget(100000, 1000000, 10000000);
  SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(budget);
  for (size_t round = 0; round < 64 && !active.empty(); ++round) {
    const int selector = solver.newVar() + 2;
    std::vector<int> badClause{-selector};
    for (int conclusion : conclusions) {
      badClause.push_back(-conclusion);
    }
    solver.addClause(badClause);
    auto assumptions = hypotheses;
    assumptions.push_back(selector);
    const auto status = solver.solveWithAssumptionsStatus(
        assumptions, 10000, 100000, 1000000);
    if (status == SATSolverWrapper::SolveStatus::Unsat) {
      std::set<std::pair<size_t, size_t>> proved;
      for (const auto* candidate : active) {
        proved.insert(candidate->equalities.begin(), candidate->equalities.end());
      }
      return {proved.begin(), proved.end()};
    }
    if (status == SATSolverWrapper::SolveStatus::Unknown) {
      return {};
    }
    size_t kept = 0;
    for (size_t i = 0; i < active.size(); ++i) {
      if (solver.getLiteralValue(conclusions[i])) {
        active[kept] = active[i];
        hypotheses[kept] = hypotheses[i];
        conclusions[kept++] = conclusions[i];
      }
    }
    active.resize(kept);
    hypotheses.resize(kept);
    conclusions.resize(kept);
    solver.addClause({-selector});
  }
  // An unfinished fixed point is not a certificate for any surviving guess.
  return {};
}

void learnInternalStateRelations(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    KInductionProblem& problem,
    const InternalRelationOptions& options,
    Config::SolverType solverType,
    bool diagnostics) {
  if (!options.learnInternalRelations) {
    return;
  }
  const size_t width = problem.usesDualRailStateEncoding ? 2 : 1;
  const std::array<const SequentialDesignModel*, 2> models{&model0, &model1};
  const std::array<const std::vector<size_t>*, 2> symbols{
      &problem.state0Symbols, &problem.state1Symbols};
  using Location = std::pair<size_t, size_t>;
  std::vector<InternalRelationCandidate> candidates;
  auto append = [&](Location lhs, Location rhs) {
    InternalRelationCandidate candidate;
    for (size_t rail = 0; rail < width; ++rail) {
      candidate.equalities.emplace_back(
          symbols[lhs.first]->at(width * lhs.second + rail),
          symbols[rhs.first]->at(width * rhs.second + rail));
    }
    if (width == 2) {
      for (const auto& location : {lhs, rhs}) {
        candidate.values.push_back({
            symbols[location.first]->at(2 * location.second),
            symbols[location.first]->at(2 * location.second + 1)});
      }
    }
    candidates.push_back(std::move(candidate));
  };
  // Names and shared drivers are guesses only. Neither participates in the
  // certificate: the SAT learner checks the original transitions and boot.
  std::unordered_map<std::string, Location> namedStates;
  for (size_t side = 0; side < models.size(); ++side) {
    const auto& model = *models[side];
    std::unordered_map<BoolExpr*, Location> sharedDrivers;
    for (size_t i = 0; i < model.stateBits.size(); ++i) {
      const auto& key = model.stateBits[i];
      const Location location{side, i};
      if (auto name = model.displayNameByKey.find(key);
          name != model.displayNameByKey.end() && !name->second.empty()) {
        auto [existing, inserted] = namedStates.emplace(name->second, location);
        if (!inserted) {
          append(existing->second, location);
        }
      }
      auto [driver, inserted] = sharedDrivers.emplace(
          model.nextStateExprByStateKey.at(key), location);
      if (!inserted) {
        append(driver->second, location);
      }
    }
  }
  const auto proved = proveInternalRelations(problem, candidates, options, solverType);
  problem.sameFrameStateEqualityPairs0.insert(
      problem.sameFrameStateEqualityPairs0.end(), proved.begin(), proved.end());
  for (const auto& [lhs, rhs] : proved) {
    BoolExpr* equality = BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(lhs), BoolExpr::Var(rhs)));
    problem.learnedInternalRelationInvariant = problem.learnedInternalRelationInvariant == nullptr
        ? equality : BoolExpr::And(problem.learnedInternalRelationInvariant, equality);
  }
  if (diagnostics) {
    printf("SEC summary: internal_relation_candidates=%zu proved_rail_equalities=%zu "
           "allow_x_equality_in_internal_relations=%d\n",
           candidates.size(), proved.size(), options.allowXEqualityInInternalRelations);
    fflush(stdout);
  }
}

}  // namespace KEPLER_FORMAL::SEC
