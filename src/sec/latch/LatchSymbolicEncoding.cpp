// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/LatchSymbolicEncoding.h"
#include <limits>
#include <set>
#include <unordered_map>

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
BoolExpr* selected(const SymbolicBits& selector, size_t value) {
  auto* result = BoolExpr::createTrue();
  for (size_t i = 0; i < selector.size(); ++i)
    result = BoolExpr::And(result, (value >> i) & 1 ? selector[i] : BoolExpr::Not(selector[i]));
  return result;
}
// Iterative, shared-DAG traversal: long unfolded latch chains must not consume
// the native stack, and substituting a variable is simultaneous, not recursive
// rewriting of its replacement (the destination may reuse a numeric symbol).
SymbolicBits rewrite(const SymbolicBits& roots, const std::unordered_map<size_t, BoolExpr*>& replacements) {
  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  SymbolicBits result;
  for (auto* root : roots) {
    std::vector<std::pair<BoolExpr*, bool>> pending{{root, false}};
    while (!pending.empty()) {
      const auto [node, ready] = pending.back();
      pending.pop_back();
      if (!node || !node->isValid()) throw std::invalid_argument("invalid symbolic macro expression");
      if (memo.count(node)) continue;
      if (node->getOp() == Op::VAR) {
        if (node->getId() < 2) memo.emplace(node, node);
        else {
          const auto found = replacements.find(node->getId());
          if (found == replacements.end()) throw std::invalid_argument("uncertified auxiliary symbol in macro");
          memo.emplace(node, found->second);
        }
        continue;
      }
      if (!ready) {
        pending.emplace_back(node, true);
        if (node->getRight()) pending.emplace_back(node->getRight(), false);
        pending.emplace_back(node->getLeft(), false);
        continue;
      }
      auto* left = memo.at(node->getLeft());
      auto* right = node->getRight() ? memo.at(node->getRight()) : nullptr;
      switch (node->getOp()) {
        case Op::NOT: memo.emplace(node, BoolExpr::Not(left)); break;
        case Op::AND: memo.emplace(node, BoolExpr::And(left, right)); break;
        case Op::OR: memo.emplace(node, BoolExpr::Or(left, right)); break;
        case Op::XOR: memo.emplace(node, BoolExpr::Xor(left, right)); break;
        default: throw std::invalid_argument("unsupported symbolic macro expression");
      }
    }
    result.push_back(memo.at(root));
  }
  return result;
}
}  // namespace

SymbolicBits encodeSymbolicInitialState(const SymbolicMacro& macro,
                                      const SymbolicBits& parameters) {
  if (parameters.size() != macro.initialParameterSymbols.size() ||
      (!macro.initialState.empty() && macro.initialState.size() != macro.stateSymbols.size()) ||
      (!macro.initialStateExpressions.empty() && macro.initialStateExpressions.size() != macro.stateSymbols.size()))
    throw std::invalid_argument("symbolic initial interface mismatch");
  if (macro.initialStateExpressions.empty()) {
    if (macro.initialState.size() != macro.stateSymbols.size() || !parameters.empty())
      throw std::invalid_argument("symbolic macro lacks its initialization relation");
    SymbolicBits result;
    for (auto bit : macro.initialState) {
      if (bit > 1) throw std::invalid_argument("non-Boolean initial state");
      result.push_back(BoolExpr::Var(bit));
    }
    return result;
  }
  std::set<size_t> ordinary(macro.stateSymbols.begin(), macro.stateSymbols.end());
  ordinary.insert(macro.inputSymbols.begin(), macro.inputSymbols.end());
  std::unordered_map<size_t, BoolExpr*> replacements;
  for (size_t i = 0; i < parameters.size(); ++i) {
    const size_t id = macro.initialParameterSymbols[i];
    if (id < 2 || ordinary.count(id) || !parameters[i] || !parameters[i]->isValid() ||
        !replacements.emplace(id, parameters[i]).second)
      throw std::invalid_argument("invalid or duplicate symbolic initial parameter");
  }
  return rewrite(macro.initialStateExpressions, replacements);
}

BoolExpr* encodeSymbolicInitialRelation(const SymbolicMacro& macro,
    const SymbolicBits& state, const SymbolicBits& parameters) {
  if (state.size() != macro.stateSymbols.size())
    throw std::invalid_argument("symbolic initial state width mismatch");
  const auto initial = encodeSymbolicInitialState(macro, parameters);
  auto* relation = BoolExpr::createTrue();
  for (size_t i = 0; i < state.size(); ++i) {
    if (!state[i] || !state[i]->isValid()) throw std::invalid_argument("invalid symbolic initial state");
    relation = BoolExpr::And(relation, BoolExpr::Not(BoolExpr::Xor(state[i], initial[i])));
  }
  return relation;
}

BoundaryEncoding encodeSymbolicMacro(const SymbolicMacro& macro, const Network& network,
    const SymbolicBits& state, const SymbolicBits& inputs, bool singleInputChange,
    const SymbolicBits& selector, BoolExpr* eventValue,
    const std::vector<size_t>& globalInputIndices) {
  if (state.size() != macro.stateSymbols.size() || state.size() < network.netCount ||
      macro.nextState.size() != state.size() ||
      (macro.initialState.size() != state.size() && macro.initialStateExpressions.size() != state.size()) ||
      inputs.size() != macro.inputSymbols.size() || inputs.size() != network.externalInputs.size() ||
      macro.observedNets.size() != network.netCount ||
      macro.singleExternalInputChange != singleInputChange ||
      macro.externalInputNets != network.externalInputs ||
      selector.size() >= std::numeric_limits<size_t>::digits ||
      (singleInputChange && (!eventValue || globalInputIndices.size() != inputs.size())))
    throw std::invalid_argument("symbolic macro interface mismatch");
  for (auto* bit : selector)
    if (!bit || !bit->isValid()) throw std::invalid_argument("invalid symbolic event selector");
  if (singleInputChange && !eventValue->isValid())
    throw std::invalid_argument("invalid symbolic event value");
  std::unordered_map<size_t, BoolExpr*> replacements;
  for (size_t i = 0; i < state.size(); ++i)
    if (macro.stateSymbols[i] < 2 || !state[i] || !state[i]->isValid() ||
        !replacements.emplace(macro.stateSymbols[i], state[i]).second)
      throw std::invalid_argument("invalid or duplicate symbolic state");
  std::set<size_t> indices;
  for (size_t i = 0; i < inputs.size(); ++i) {
    auto* input = inputs[i];
    if (singleInputChange) {
      const size_t global = globalInputIndices[i];
      if (global >= (size_t(1) << selector.size()) || !indices.insert(global).second)
        throw std::invalid_argument("invalid symbolic event selector index");
      input = symbolicMux(selected(selector, global), eventValue, state.at(network.externalInputs[i]));
    }
    if (macro.inputSymbols[i] < 2 || !input || !input->isValid() ||
        !replacements.emplace(macro.inputSymbols[i], input).second)
      throw std::invalid_argument("invalid or duplicate symbolic input");
  }
  auto roots = macro.nextState;
  roots.insert(roots.end(), macro.observedNets.begin(), macro.observedNets.end());
  const auto encoded = rewrite(roots, replacements);
  return {{encoded.begin(), encoded.begin() + state.size()},
          {encoded.begin() + state.size(), encoded.end()}};
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
