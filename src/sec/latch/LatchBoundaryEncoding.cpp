// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/LatchBoundaryEncoding.h"
#include <stdexcept>
#include <limits>
#include <set>

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
BoolExpr* equals(const std::vector<BoolExpr*>& bits, size_t value) {
  BoolExpr* result = BoolExpr::createTrue();
  for (size_t i = 0; i < bits.size(); ++i)
    result = BoolExpr::And(result, (value >> i) & 1 ? bits[i] : BoolExpr::Not(bits[i]));
  return result;
}
BoolExpr* literal(BoolExpr* expression, bool value) {
  return value ? expression : BoolExpr::Not(expression);
}
}  // namespace

size_t boundaryEncodingBits(size_t states) {
  if (!states) throw std::invalid_argument("empty boundary table");
  size_t bits = 1;
  for (size_t remaining = states - 1; remaining >>= 1;) ++bits;
  return bits;
}

BoundaryEncoding encodeBoundaryTable(
    const TransitionTable& table, const Network& network,
    const std::vector<BoolExpr*>& state, const std::vector<BoolExpr*>& inputs,
    const std::vector<BoolExpr*>& selector, BoolExpr* eventValue,
    const std::vector<size_t>& globalInputIndices) {
  if (state.size() != boundaryEncodingBits(table.boundaries.size()) ||
      inputs.size() != network.externalInputs.size() ||
      (table.singleExternalInputChange &&
       (!eventValue || globalInputIndices.size() != inputs.size())))
    throw std::invalid_argument("boundary encoding interface mismatch");
  if (selector.size() >= std::numeric_limits<size_t>::digits)
    throw std::invalid_argument("event selector exceeds index width");
  std::set<size_t> uniqueIndices;
  for (auto index : globalInputIndices)
    if (index >= (size_t(1) << selector.size()) || !uniqueIndices.insert(index).second)
      throw std::invalid_argument("aliased or out-of-range event selector index");
  BoundaryEncoding encoded;
  encoded.nextState.assign(state.size(), BoolExpr::createFalse());
  encoded.observedNets.assign(network.netCount, BoolExpr::createFalse());
  std::vector<BoolExpr*> fromConditions;
  for (size_t i = 0; i < table.boundaries.size(); ++i)
    fromConditions.push_back(equals(state, i));
  std::vector<BoolExpr*> selected;
  BoolExpr* selectsThisComponent = BoolExpr::createFalse();
  if (table.singleExternalInputChange) {
    for (auto global : globalInputIndices) {
      auto* condition = equals(selector, global);
      selected.push_back(condition);
      selectsThisComponent = BoolExpr::Or(selectsThisComponent, condition);
    }
  }
  for (const auto& row : table.rows) {
    if (row.from >= table.boundaries.size() || row.to >= table.boundaries.size() ||
        row.input.size() != inputs.size())
      throw std::invalid_argument("invalid certified boundary row");
    BoolExpr* inputCondition = BoolExpr::createTrue();
    if (!table.singleExternalInputChange) {
      for (size_t i = 0; i < inputs.size(); ++i)
        inputCondition = BoolExpr::And(inputCondition, literal(inputs[i], row.input[i]));
    } else {
      inputCondition = BoolExpr::Not(selectsThisComponent);
      size_t differences = 0;
      BoolExpr* changed = nullptr;
      for (size_t i = 0; i < inputs.size(); ++i) {
        const bool old = table.boundaries[row.from].current.at(network.externalInputs[i]);
        auto* choose = BoolExpr::And(selected[i], literal(eventValue, row.input[i]));
        if (row.input[i] != old) { ++differences; changed = choose; }
        else inputCondition = BoolExpr::Or(inputCondition, choose);
      }
      if (differences > 1) throw std::invalid_argument("non-single-input certified row");
      if (differences == 1) inputCondition = changed;
    }
    auto* condition = BoolExpr::And(fromConditions[row.from], inputCondition);
    for (size_t i = 0; i < state.size(); ++i)
      if ((row.to >> i) & 1)
        encoded.nextState[i] = BoolExpr::Or(encoded.nextState[i], condition);
    const auto& next = table.boundaries[row.to];
    for (size_t i = 0; i < network.netCount; ++i)
      if (next.current.at(i))
        encoded.observedNets[i] = BoolExpr::Or(encoded.observedNets[i], condition);
  }
  // Out-of-range IDs are unreachable from the explicitly initialized state and
  // closed certified table. Their total encoding is zero, never an assumption.
  return encoded;
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
