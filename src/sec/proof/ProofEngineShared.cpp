// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "proof/ProofEngineShared.h"

#include <algorithm>
#include <unordered_set>

#include "../../config/Config.h"
#include "common/BoolExprUtils.h"
#include "kinduction/SatEncoding.h"

namespace KEPLER_FORMAL::SEC {

namespace {

BoolExpr* appendStructuredAssignmentFacts(
    BoolExpr* init,
    const std::vector<std::pair<size_t, bool>>& assignments,
    bool& hasConstraint) {
  for (const auto& [symbol, value] : assignments) {
    init = BoolExpr::And(
        init, value ? BoolExpr::Var(symbol) : BoolExpr::Not(BoolExpr::Var(symbol)));
    hasConstraint = true;
  }
  return init;
}

std::vector<size_t> sortUniqueSymbols(std::vector<size_t> symbols) {
  std::sort(symbols.begin(), symbols.end());
  symbols.erase(std::unique(symbols.begin(), symbols.end()), symbols.end());
  return symbols;
}

std::vector<size_t> sortUniqueSymbols(std::unordered_set<size_t> symbols) {
  return sortUniqueSymbols(
      std::vector<size_t>(symbols.begin(), symbols.end()));
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

std::vector<size_t> buildFormulaSupportVector(BoolExpr* formula) {
  std::vector<size_t> support;
  if (formula == nullptr) {
    // LCOV_EXCL_START
    return support;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  const auto supportSet = formula->getSupportVars();
  support.reserve(supportSet.size());
  for (const auto symbol : supportSet) {
    support.push_back(symbol);
  }
  return support;
}

const std::vector<size_t>& cachedFormulaSupport(
    BoolExpr* formula,
    FormulaSupportCache& supportCache) {
  auto [it, inserted] = supportCache.emplace(formula, std::vector<size_t>{});
  if (inserted) {
    // Invariant validation reuses the same large transition and strengthening
    // DAGs across several PDR candidates. Cache their support locally to avoid
    // repeatedly walking and allocating for identical BoolExpr subgraphs.
    it->second = buildFormulaSupportVector(formula);
  }
  return it->second;
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

void addSupportSymbols(const std::vector<size_t>& support,
                       std::unordered_set<size_t>& symbols) {
  for (const auto symbol : support) {
    if (symbol >= 2) {
      symbols.insert(symbol);
    }
  }
}

std::unordered_set<size_t> buildCombinedStateSymbolSet(
    const KInductionProblem& problem);

std::vector<size_t> collectStateSupportSymbols(
    const KInductionProblem& problem,
    const std::vector<size_t>& formulaSupport) {
  std::vector<size_t> support;
  const auto stateSymbolSet = buildCombinedStateSymbolSet(problem);
  for (const auto symbol : formulaSupport) {
    if (stateSymbolSet.find(symbol) != stateSymbolSet.end()) {
      support.push_back(symbol);
    }
  }
  std::sort(support.begin(), support.end());
  support.erase(std::unique(support.begin(), support.end()), support.end());
  return support;
}

void addComplementedStateRelations(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    size_t numFrames) {
  for (size_t frame = 0; frame < numFrames; ++frame) {
    for (const auto& [primarySymbol, complementedSymbol] : complementedStatePairs) {
      if (!variables.hasSymbol(primarySymbol) ||
          !variables.hasSymbol(complementedSymbol)) {
        continue; // LCOV_EXCL_LINE
      }
      addLiteralEquivalence(
          solver,
          variables.getLiteral(complementedSymbol, frame),
          -variables.getLiteral(primarySymbol, frame));
    }
  }
}

void addDualRailStateValidity(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const std::vector<DualRailSymbolPair>& railPairs,
    size_t numFrames) {
  for (size_t frame = 0; frame < numFrames; ++frame) {
    for (const auto& rails : railPairs) {
      if (!variables.hasSymbol(rails.mayBeOne) ||
          !variables.hasSymbol(rails.mayBeZero)) {
        // LCOV_EXCL_START
        continue;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // Valid dual-rail values are 0, 1, and X.  The empty set is not a
      // reachable ternary value and must not be considered by invariant checks.
      solver.addClause({
          variables.getLiteral(rails.mayBeOne, frame),
          variables.getLiteral(rails.mayBeZero, frame)});
    }
  }
}

std::unordered_map<size_t, BoolExpr*> buildTransitionExprByStateSymbol(
    const KInductionProblem& problem) {
  std::unordered_map<size_t, BoolExpr*> transitionExprByStateSymbol;
  transitionExprByStateSymbol.reserve(
      problem.transitions0.size() + problem.transitions1.size() +
      problem.auxiliaryTransitions.size());
  for (const auto& [stateSymbol, expr] : problem.transitions0) {
    transitionExprByStateSymbol.emplace(stateSymbol, expr);
  }
  for (const auto& [stateSymbol, expr] : problem.transitions1) {
    transitionExprByStateSymbol.emplace(stateSymbol, expr);
  }
  for (const auto& [stateSymbol, expr] : problem.auxiliaryTransitions) {
    transitionExprByStateSymbol.emplace(stateSymbol, expr);
  }
  return transitionExprByStateSymbol;
}

std::unordered_map<size_t, size_t> buildComplementPrimaryByStateSymbol(
    const KInductionProblem& problem) {
  std::unordered_map<size_t, size_t> primaryByComplement;
  primaryByComplement.reserve(
      problem.complementedStatePairs0.size() +
      problem.complementedStatePairs1.size());
  for (const auto& [primarySymbol, complementedSymbol] :
       problem.complementedStatePairs0) {
    primaryByComplement.emplace(complementedSymbol, primarySymbol);
  }
  for (const auto& [primarySymbol, complementedSymbol] :
       problem.complementedStatePairs1) {
    primaryByComplement.emplace(complementedSymbol, primarySymbol);
  }
  return primaryByComplement;
}

std::unordered_set<size_t> buildCombinedStateSymbolSet(
    const KInductionProblem& problem) {
  std::unordered_set<size_t> stateSymbols;
  stateSymbols.reserve(
      problem.state0Symbols.size() + problem.state1Symbols.size() +
      problem.auxiliaryStateSymbols.size());
  stateSymbols.insert(problem.state0Symbols.begin(), problem.state0Symbols.end());
  stateSymbols.insert(problem.state1Symbols.begin(), problem.state1Symbols.end());
  stateSymbols.insert(
      problem.auxiliaryStateSymbols.begin(),
      problem.auxiliaryStateSymbols.end());
  return stateSymbols;
}

void addRelevantComplementedStatePartners(
    const std::vector<std::pair<size_t, size_t>>& complementedStatePairs,
    std::unordered_set<size_t>& symbols) {
  for (const auto& [primarySymbol, complementedSymbol] : complementedStatePairs) {
    if (symbols.find(primarySymbol) != symbols.end() ||
        symbols.find(complementedSymbol) != symbols.end()) { // LCOV_EXCL_LINE
      symbols.insert(primarySymbol);
      symbols.insert(complementedSymbol);
    }
  }
}

void addRelevantDualRailPartners(
    const std::vector<DualRailSymbolPair>& railPairs,
    std::unordered_set<size_t>& symbols) {
  for (const auto& rails : railPairs) {
    if (symbols.find(rails.mayBeOne) != symbols.end() ||
        // LCOV_EXCL_START
        symbols.find(rails.mayBeZero) != symbols.end()) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      symbols.insert(rails.mayBeOne);
      symbols.insert(rails.mayBeZero);
    }
  }
}

std::vector<size_t> expandTransitionTargets(
    const KInductionProblem& problem,
    const std::vector<size_t>& requestedTargets,
    const std::unordered_map<size_t, BoolExpr*>& transitionExprByStateSymbol) {
  const auto primaryByComplement = buildComplementPrimaryByStateSymbol(problem);
  std::vector<size_t> targets;
  targets.reserve(requestedTargets.size());

  for (const auto symbol : requestedTargets) {
    if (transitionExprByStateSymbol.find(symbol) !=
        transitionExprByStateSymbol.end()) {
      targets.push_back(symbol);
      continue;
    }
    // LCOV_EXCL_START
    if (const auto primaryIt = primaryByComplement.find(symbol);  // LCOV_EXCL_LINE
        primaryIt != primaryByComplement.end() &&  // LCOV_EXCL_LINE
        transitionExprByStateSymbol.find(primaryIt->second) !=  // LCOV_EXCL_LINE
            transitionExprByStateSymbol.end()) {  // LCOV_EXCL_LINE
      targets.push_back(primaryIt->second);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  return sortUniqueSymbols(std::move(targets));
}

std::vector<size_t> inductiveInvariantQuerySymbols(
    const KInductionProblem& problem,
    BoolExpr* invariant,
    const std::unordered_map<size_t, BoolExpr*>& transitionExprByStateSymbol,
    FormulaSupportCache& supportCache) {
  std::unordered_set<size_t> symbols;
  const auto& invariantSupport = cachedFormulaSupport(invariant, supportCache);
  addSupportSymbols(invariantSupport, symbols);

  const auto targets = expandTransitionTargets(
      problem,
      collectStateSupportSymbols(problem, invariantSupport),
      transitionExprByStateSymbol);
  for (const auto stateSymbol : targets) {
    symbols.insert(stateSymbol);
    addSupportSymbols(
        cachedFormulaSupport(
            transitionExprByStateSymbol.at(stateSymbol), supportCache),
        symbols);
  }

  // Reset-bootstrap PDR frames are post-reset frames.  If the invariant or the
  // transition cone mentions reset controls, the inductiveness query must see
  // the same deasserted reset environment used by PDR's blocking queries.
  if (problem.resetBootstrapCycles != 0) {
    for (const auto& [symbol, _] : problem.resetBootstrapInputs) {
      symbols.insert(symbol);  // LCOV_EXCL_LINE
    }
  }

  addRelevantComplementedStatePartners(problem.complementedStatePairs0, symbols);
  addRelevantComplementedStatePartners(problem.complementedStatePairs1, symbols);
  addRelevantDualRailPartners(problem.dualRailStatePairs, symbols);
  return sortUniqueSymbols(std::move(symbols));
}

void addTransitionRelationForTargets(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    const std::unordered_map<size_t, BoolExpr*>& transitionExprByStateSymbol,
    size_t frame,
    const std::vector<size_t>& requestedTargets) {
  const auto encodedTargets = expandTransitionTargets(
      problem, requestedTargets, transitionExprByStateSymbol);

  FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
  for (const auto stateSymbol : encodedTargets) {
    addLiteralEquivalence(
        solver,
        variables.getLiteral(stateSymbol, frame + 1),
        encoder.encode(transitionExprByStateSymbol.at(stateSymbol)));
  }
}

void addPostBootstrapResetInputConstraints(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    size_t frame) {
  if (problem.resetBootstrapCycles == 0) {
    return;
  }

  for (const auto& [symbol, assertedValue] : problem.resetBootstrapInputs) {
    if (!variables.hasSymbol(symbol)) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    solver.addClause(  // LCOV_EXCL_LINE
        {assertedValue ? -variables.getLiteral(symbol, frame)  // LCOV_EXCL_LINE
                       // LCOV_EXCL_START
                       : variables.getLiteral(symbol, frame)});  // LCOV_EXCL_LINE
                       // LCOV_EXCL_STOP
  }
}

}  // namespace

BoolExpr* buildProofInitFormula(const KInductionProblem& problem) {
  BoolExpr* init = BoolExpr::createTrue();
  bool hasConstraint = false;

  if (problem.resetBootstrapCycles != 0) {
    for (const auto& [symbol, value] : problem.bootstrapStateAssignments) {
      init = BoolExpr::And(
          init, value ? BoolExpr::Var(symbol) : BoolExpr::Not(BoolExpr::Var(symbol)));
      hasConstraint = true;
    }
  } else {
    if (problem.initialCondition == BoolExpr::createTrue() &&
        !problem.initialStateAssignments.empty()) {
      // Dual-rail SEC keeps the boot rails as structured unit facts so PDR and
      // k-induction can encode only the local COI.  Formula-based callers such
      // as exact IMC still need a real BoolExpr frontier, so materialize it
      // here rather than during SEC problem construction.
      init = appendStructuredAssignmentFacts(
          init, problem.initialStateAssignments, hasConstraint);
    } else if (problem.initialCondition != nullptr) {
      init = BoolExpr::And(init, problem.initialCondition);
      hasConstraint = true;
    }
    const bool needsObservationFrontier =
        problem.hasSequentialState() && problem.property != nullptr &&
        ((!problem.hasExplicitInitialState()) ||
         (problem.hasExplicitInitialState() &&
          !problem.hasCompleteInitialState()));
    if (needsObservationFrontier) {
      init = BoolExpr::And(init, problem.property);
      hasConstraint = true;
    }
  }

  if (!hasConstraint) {
    return nullptr;
  }
  return BoolExpr::simplify(init);
}

size_t nextFreshProofSymbol(const KInductionProblem& problem) {
  size_t nextSymbol = 2;
  for (const auto symbol : problem.allSymbols) {
    nextSymbol = std::max(nextSymbol, symbol + 1);
  }
  return nextSymbol;
}

std::unordered_map<size_t, size_t> allocateFreshProofSymbols(
    const std::vector<size_t>& originalSymbols,
    size_t& nextSymbol) {
  std::unordered_map<size_t, size_t> symbolMap;
  symbolMap.reserve(originalSymbols.size());
  for (const auto symbol : originalSymbols) {
    symbolMap.emplace(symbol, nextSymbol++);
  }
  return symbolMap;
}

BoolExpr* buildOneStepTransitionFormula(
    const KInductionProblem& problem,
    const std::unordered_map<size_t, size_t>& nextStateSymbols) {
  BoolExpr* transition = BoolExpr::createTrue();
  for (const auto& [stateSymbol, expr] : problem.transitions0) {
    transition = BoolExpr::And(
        transition,
        makeEqualityExpr(BoolExpr::Var(nextStateSymbols.at(stateSymbol)), expr));
  }
  for (const auto& [stateSymbol, expr] : problem.transitions1) {
    transition = BoolExpr::And(
        transition,
        makeEqualityExpr(BoolExpr::Var(nextStateSymbols.at(stateSymbol)), expr));
  }
  for (const auto& [stateSymbol, expr] : problem.auxiliaryTransitions) {
    transition = BoolExpr::And(
        transition,
        makeEqualityExpr(BoolExpr::Var(nextStateSymbols.at(stateSymbol)), expr));
  }
  for (const auto& [primarySymbol, complementedSymbol] : problem.complementedStatePairs0) {
    transition = BoolExpr::And(
        transition,
        makeEqualityExpr(
            BoolExpr::Var(nextStateSymbols.at(complementedSymbol)),
            BoolExpr::Not(BoolExpr::Var(nextStateSymbols.at(primarySymbol)))));
  }
  for (const auto& [primarySymbol, complementedSymbol] : problem.complementedStatePairs1) {
    transition = BoolExpr::And(
        transition,
        makeEqualityExpr(
            BoolExpr::Var(nextStateSymbols.at(complementedSymbol)),
            BoolExpr::Not(BoolExpr::Var(nextStateSymbols.at(primarySymbol)))));
  }
  return BoolExpr::simplify(transition);
}

BoolExpr* buildCurrentStateLegalityFormula(const KInductionProblem& problem) {
  BoolExpr* legality = BoolExpr::createTrue();
  for (const auto& [primarySymbol, complementedSymbol] : problem.complementedStatePairs0) {
    // LCOV_EXCL_START
    legality = BoolExpr::And(  // LCOV_EXCL_LINE
        legality,  // LCOV_EXCL_LINE
        makeEqualityExpr(  // LCOV_EXCL_LINE
            BoolExpr::Var(complementedSymbol), BoolExpr::Not(BoolExpr::Var(primarySymbol))));  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
  }
  for (const auto& [primarySymbol, complementedSymbol] : problem.complementedStatePairs1) {
    legality = BoolExpr::And(
        legality,
        makeEqualityExpr(
            BoolExpr::Var(complementedSymbol), BoolExpr::Not(BoolExpr::Var(primarySymbol))));
  }
  return BoolExpr::simplify(legality);
}

BoolExpr* remapProofFormula(
    BoolExpr* formula,
    const std::unordered_map<size_t, size_t>& symbolMap) {
  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  return remapBoolExprVariables(formula, symbolMap, memo);
}

bool isProofFormulaSatisfiable(
    BoolExpr* formula,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (formula == nullptr) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  SATSolverWrapper solver(solverType);
  const auto support = formula->getSupportVars();
  std::unordered_map<size_t, int> leafLits;
  leafLits.reserve(support.size());
  for (const auto symbol : support) {
    if (symbol < 2) {
      continue;
    }
    leafLits.emplace(symbol, solver.newVar() + 2);
  }

  FrameFormulaEncoder encoder(solver, std::move(leafLits));
  solver.addClause({encoder.encode(formula)});
  return solver.solve();
}

bool initialFrontierImplies(
    BoolExpr* initFormula,
    BoolExpr* invariant,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (initFormula == nullptr || invariant == nullptr) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  return !isProofFormulaSatisfiable(
      BoolExpr::And(initFormula, BoolExpr::Not(invariant)), solverType);
}

BoolExpr* selectValidatedStrengtheningInvariant(
    const KInductionProblem& problem,
    BoolExpr* initFormula,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (problem.inductionProperty == nullptr || initFormula == nullptr ||
      problem.inductionProperty == problem.property) {
    return nullptr;
  }

  if (!initialFrontierImplies(initFormula, problem.inductionProperty, solverType)) {
    return nullptr; // LCOV_EXCL_LINE
  }
  return problem.inductionProperty;
}

bool invariantExcludesBadStates(
    const KInductionProblem& problem,
    BoolExpr* invariant,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (invariant == nullptr || problem.bad == nullptr) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  return !isProofFormulaSatisfiable(
      BoolExpr::And(invariant, problem.bad), solverType);
}

bool isInductiveInvariant(
    const KInductionProblem& problem,
    BoolExpr* invariant,
    KEPLER_FORMAL::Config::SolverType solverType) {
  FormulaSupportCache supportCache;
  return isInductiveInvariant(problem, invariant, solverType, supportCache);
}

bool isInductiveInvariant(
    const KInductionProblem& problem,
    BoolExpr* invariant,
    KEPLER_FORMAL::Config::SolverType solverType,
    FormulaSupportCache& supportCache) {
  if (invariant == nullptr) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  const auto transitionExprByStateSymbol =
      buildTransitionExprByStateSymbol(problem);
  const auto querySymbols =
      inductiveInvariantQuerySymbols(
          problem, invariant, transitionExprByStateSymbol, supportCache);
  const auto& invariantSupport = cachedFormulaSupport(invariant, supportCache);
  const auto invariantStateSupport =
      collectStateSupportSymbols(problem, invariantSupport);

  SATSolverWrapper solver(solverType);
  FrameVariableStore variables(solver, querySymbols, 2);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs0, 2);
  addComplementedStateRelations(solver, variables, problem.complementedStatePairs1, 2);
  addDualRailStateValidity(solver, variables, problem.dualRailStatePairs, 2);
  addPostBootstrapResetInputConstraints(solver, variables, problem, 0);
  // Only next-state symbols read by the candidate invariant need transition
  // equations. Encoding every flop transition here made PDR's immediate
  // invariant validation scale like a full-design induction proof even when
  // the candidate touched a small cone.
  addTransitionRelationForTargets(
      solver,
      variables,
      problem,
      transitionExprByStateSymbol,
      0,
      invariantStateSupport);

  FrameFormulaEncoder currentEncoder(
      solver, variables.makeLeafLits(0, invariantSupport));
  FrameFormulaEncoder nextEncoder(
      solver, variables.makeLeafLits(1, invariantSupport));
  solver.addClause({currentEncoder.encode(invariant)});
  solver.addClause({nextEncoder.encode(BoolExpr::Not(invariant))});
  return !solver.solve();
}

}  // namespace KEPLER_FORMAL::SEC
