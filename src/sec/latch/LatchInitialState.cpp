// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/LatchInitialState.h"

#include <stdexcept>
#include <algorithm>
#include "common/BoolExprUtils.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
void append(KInductionProblem& problem, BoolExpr* expression) {
  problem.initialCondition = BoolExpr::And(
      problem.initialCondition ? problem.initialCondition : BoolExpr::createTrue(),
      expression);
}

void addRelation(const SequentialDesignModel& model,
                 const std::unordered_map<size_t, size_t>& symbols,
                 KInductionProblem& problem) {
  if (!model.initialCondition) return;
  for (auto symbol : model.initialCondition->getSupportVars())
    if (symbol >= 2 && !symbols.count(symbol))
      throw std::runtime_error("Event initial relation references an undeclared symbol");
  append(problem, remapBoolExprVariables(model.initialCondition, symbols));
}
}  // namespace

void reuseInitialInputHistory(SequentialDesignModel& model,
    const std::vector<SignalKey>& inputs, const std::vector<BoolExpr*>& history) {
  std::unordered_map<size_t, size_t> replacements;
  for (size_t i = 0; i < inputs.size(); ++i) {
    const auto origin = model.initialInputStateKeyByInputKey.find(inputs[i]);
    if (origin == model.initialInputStateKeyByInputKey.end()) continue;
    if (!history.at(i) || history[i]->getOp() != Op::VAR || history[i]->getId() < 2)
      throw std::runtime_error("Symbolic initial input lacks a remembered state bit");
    const auto originalKey = origin->second;
    const auto id = model.inputVarByKey.at(originalKey);
    const auto current = std::find_if(model.stateBits.begin(), model.stateBits.end(),
        [&](const auto& key) { return model.inputVarByKey.at(key) == history[i]->getId(); });
    if (current == model.stateBits.end() || history[i]->getId() == id)
      throw std::runtime_error("Symbolic input history is not a distinct state bit");
    origin->second = *current;
    replacements.emplace(id, history[i]->getId());
    std::erase(model.stateBits, originalKey);
    model.nextStateExprByStateKey.erase(originalKey);
    model.inputVarByKey.erase(originalKey);
    model.displayNameByKey.erase(originalKey);
  }
  if (!replacements.empty()) {
    for (auto symbol : model.initialCondition->getSupportVars())
      if (symbol >= 2) replacements.try_emplace(symbol, symbol);
    model.initialCondition = remapBoolExprVariables(model.initialCondition, replacements);
  }
  // A BOOT origin erased by transparency/reset and absent from the entire
  // initial relation is existentially irrelevant. Macro transitions/outputs
  // contain only boundary state and event inputs, never BOOT parameters.
  const auto initialSupport = model.initialCondition->getSupportVars();
  std::erase_if(model.stateBits, [&](const auto& key) {
    if (key.first.size() != 3 || key.first[0] != (uint64_t(1) << 61) || key.first[1] != 8 ||
        initialSupport.count(model.inputVarByKey.at(key))) return false;
    model.nextStateExprByStateKey.erase(key);
    model.inputVarByKey.erase(key);
    model.displayNameByKey.erase(key);
    return true;
  });
}

void integrateEventInitialState(
    const SequentialDesignModel& first, const SequentialDesignModel& second,
    const AlignedSignals& inputs,
    const std::unordered_map<size_t, size_t>& symbols0,
    const std::unordered_map<size_t, size_t>& symbols1,
    KInductionProblem& problem) {
  if (!first.initialCondition && !second.initialCondition) return;
  problem.hasExactRelationalInitialState = true;
  addRelation(first, symbols0, problem);
  addRelation(second, symbols1, problem);
  for (const auto& [symbol, value] : problem.initialStateAssignments)
    append(problem, value ? BoolExpr::Var(symbol) : BoolExpr::Not(BoolExpr::Var(symbol)));
  for (size_t i = 0; i < inputs.keys0.size(); ++i) {
    const auto lhs = first.initialInputStateKeyByInputKey.find(inputs.keys0[i]);
    const auto rhs = second.initialInputStateKeyByInputKey.find(inputs.keys1[i]);
    if (lhs == first.initialInputStateKeyByInputKey.end() ||
        rhs == second.initialInputStateKeyByInputKey.end()) continue;
    append(problem, makeEqualityExpr(
        BoolExpr::Var(symbols0.at(first.inputVarByKey.at(lhs->second))),
        BoolExpr::Var(symbols1.at(second.inputVarByKey.at(rhs->second)))));
  }
}

void constrainEventInitialRails(
    const SequentialDesignModel& model,
    const std::unordered_map<size_t, DualRailSymbolPair>& rails,
    KInductionProblem& problem, bool secondDesign) {
  if (!model.initialCondition) return;
  auto& complements = secondDesign ? problem.complementedStatePairs1 : problem.complementedStatePairs0;
  for (const auto& [symbol, pair] : rails) {
    append(problem, BoolExpr::Xor(BoolExpr::Var(pair.mayBeOne),
                                  BoolExpr::Var(pair.mayBeZero)));
    // BOOT is Boolean and every event next-state expression is Boolean over
    // Boolean state/inputs. Dual-rail lifting therefore preserves complements
    // inductively. This excludes unreachable X states, not genuine starts.
    complements.emplace_back(pair.mayBeOne, pair.mayBeZero);
  }
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
