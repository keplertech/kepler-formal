// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <algorithm>
#include <unordered_map>

#include "BoolExprCache.h"
#include "latch/LatchBoundaryEncoding.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

class LatchBoundaryEncodingTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }

  TransitionTable compile(const Network& network, bool singleChange) {
    CompileOptions options;
    options.initialInputs.assign(network.externalInputs.size(), 0);
    options.singleExternalInputChange = singleChange;
    for (const auto& primitive : network.primitives)
      options.initialStorage.emplace_back(primitive.storageBits, uint8_t{0});
    auto result = compileTransitionTable(EventModel(network), options);
    EXPECT_TRUE(result.certified()) << result.detail;
    return result.table.value();
  }

  std::vector<BoolExpr*> variables(size_t count, size_t base) {
    std::vector<BoolExpr*> result;
    for (size_t i = 0; i < count; ++i) result.push_back(BoolExpr::Var(base + i));
    return result;
  }

  void assign(std::unordered_map<size_t, bool>& environment,
              const std::vector<BoolExpr*>& variables, size_t value) {
    for (size_t i = 0; i < variables.size(); ++i)
      environment[variables[i]->getId()] = (value >> i) & 1;
  }

  size_t evaluateId(const BoundaryEncoding& encoding,
                    const std::unordered_map<size_t, bool>& environment) {
    size_t result = 0;
    for (size_t i = 0; i < encoding.nextState.size(); ++i)
      result |= size_t{encoding.nextState[i]->evaluate(environment)} << i;
    return result;
  }

  void expectRow(const BoundaryEncoding& encoding, const TransitionTable& table,
                 size_t target, const std::unordered_map<size_t, bool>& environment) {
    EXPECT_EQ(evaluateId(encoding, environment), target);
    ASSERT_EQ(encoding.observedNets.size(), table.boundaries[target].current.size());
    for (size_t net = 0; net < encoding.observedNets.size(); ++net)
      EXPECT_EQ(encoding.observedNets[net]->evaluate(environment),
                table.boundaries[target].current[net] != 0) << "net=" << net;
  }
};

TEST_F(LatchBoundaryEncodingTests, StateWidthCoversEveryIdentifierIncludingSingleton) {
  EXPECT_THROW(boundaryEncodingBits(0), std::invalid_argument);
  EXPECT_EQ(boundaryEncodingBits(1), 1u);
  EXPECT_EQ(boundaryEncodingBits(2), 1u);
  EXPECT_EQ(boundaryEncodingBits(3), 2u);
  EXPECT_EQ(boundaryEncodingBits(4), 2u);
  EXPECT_EQ(boundaryEncodingBits(5), 3u);
  EXPECT_EQ(boundaryEncodingBits(256), 8u);
  EXPECT_EQ(boundaryEncodingBits(257), 9u);
}

TEST_F(LatchBoundaryEncodingTests, EncodesEveryAllChangeTruthTableRowExactly) {
  Network network{3, {0, 1}, {combinational("xor", {0, 1}, {2}, [](const Bits& pins) {
    return Bits{static_cast<uint8_t>(pins[0] ^ pins[1])};
  })}};
  const auto table = compile(network, false);
  const auto state = variables(boundaryEncodingBits(table.boundaries.size()), 2);
  const auto input = variables(2, 20);
  const auto encoding = encodeBoundaryTable(table, network, state, input);
  for (const auto& row : table.rows) {
    std::unordered_map<size_t, bool> environment;
    assign(environment, state, row.from);
    for (size_t i = 0; i < row.input.size(); ++i) environment[input[i]->getId()] = row.input[i];
    expectRow(encoding, table, row.to, environment);
  }
}

TEST_F(LatchBoundaryEncodingTests, SharedSelectorChangesOnlyItsNamedGlobalInput) {
  Network network{3, {0, 1}, {latch("l", 0, 1, 2)}};
  const auto table = compile(network, true);
  const auto state = variables(boundaryEncodingBits(table.boundaries.size()), 2);
  const auto inputs = variables(2, 20);
  const auto selector = variables(3, 30);
  auto* value = BoolExpr::Var(40);
  const std::vector<size_t> global{2, 5};
  const auto encoding = encodeBoundaryTable(table, network, state, inputs, selector, value, global);
  for (size_t from = 0; from < table.boundaries.size(); ++from) {
    for (size_t selected = 0; selected < 8; ++selected) {
      for (size_t event = 0; event < 2; ++event) {
        Bits expectedInput;
        for (size_t i = 0; i < global.size(); ++i)
          expectedInput.push_back(global[i] == selected ? event :
              table.boundaries[from].current[network.externalInputs[i]]);
        const auto row = std::find_if(table.rows.begin(), table.rows.end(), [&](const auto& row) {
          return row.from == from && row.input == expectedInput;
        });
        ASSERT_NE(row, table.rows.end());
        std::unordered_map<size_t, bool> environment;
        assign(environment, state, from);
        assign(environment, selector, selected);
        environment[value->getId()] = event;
        // Original full-valued input variables must not constrain selector mode.
        environment[inputs[0]->getId()] = !expectedInput[0];
        environment[inputs[1]->getId()] = !expectedInput[1];
        expectRow(encoding, table, row->to, environment);
      }
    }
  }
}

TEST_F(LatchBoundaryEncodingTests, SharedSelectorMaintainsIndependentComponents) {
  const Network network{2, {0}, {combinational("buffer", {0}, {1}, [](const Bits& pins) { return pins; })}};
  const auto table = compile(network, true);
  const auto selector = variables(2, 30);
  auto* value = BoolExpr::Var(40);
  const auto firstState = variables(boundaryEncodingBits(table.boundaries.size()), 2);
  const auto secondState = variables(boundaryEncodingBits(table.boundaries.size()), 10);
  const auto first = encodeBoundaryTable(table, network, firstState, {BoolExpr::createFalse()}, selector, value, {0});
  const auto second = encodeBoundaryTable(table, network, secondState, {BoolExpr::createFalse()}, selector, value, {1});
  std::unordered_map<size_t, bool> environment;
  assign(environment, firstState, table.initials.front().boundary);
  assign(environment, secondState, table.initials.front().boundary);
  assign(environment, selector, 0);
  environment[40] = true;
  EXPECT_TRUE(first.observedNets[1]->evaluate(environment));
  EXPECT_FALSE(second.observedNets[1]->evaluate(environment));
  assign(environment, selector, 1);
  EXPECT_FALSE(first.observedNets[1]->evaluate(environment));
  EXPECT_TRUE(second.observedNets[1]->evaluate(environment));
}

TEST_F(LatchBoundaryEncodingTests, StateIdsPreserveHiddenInputHistoryNotOnlyOutput) {
  const Network network{3, {0, 1}, {latch("l", 0, 1, 2)}};
  const auto table = compile(network, true);
  size_t closedZero = table.boundaries.size(), closedOne = table.boundaries.size();
  for (size_t i = 0; i < table.boundaries.size(); ++i) {
    if (table.boundaries[i].current == Bits{0, 0, 0}) closedZero = i;
    if (table.boundaries[i].current == Bits{1, 0, 0}) closedOne = i;
  }
  ASSERT_LT(closedZero, table.boundaries.size());
  ASSERT_LT(closedOne, table.boundaries.size());
  ASSERT_NE(closedZero, closedOne);
  const auto state = variables(boundaryEncodingBits(table.boundaries.size()), 2);
  const auto selector = variables(2, 20);
  const auto encoding = encodeBoundaryTable(table, network, state,
      {BoolExpr::createFalse(), BoolExpr::createFalse()}, selector, BoolExpr::Var(30), {0, 1});
  std::unordered_map<size_t, bool> environment{{30, true}};
  assign(environment, selector, 1);  // Open without changing remembered data input.
  assign(environment, state, closedZero);
  EXPECT_FALSE(encoding.observedNets[2]->evaluate(environment));
  assign(environment, state, closedOne);
  EXPECT_TRUE(encoding.observedNets[2]->evaluate(environment));
}

TEST_F(LatchBoundaryEncodingTests, UnreachableStateIdsHaveTotalZeroEncoding) {
  Network network{3, {0, 1}, {latch("l", 0, 1, 2)}};
  const auto table = compile(network, true);
  const auto state = variables(boundaryEncodingBits(table.boundaries.size()), 2);
  const auto selector = variables(2, 20);
  const auto encoding = encodeBoundaryTable(table, network, state,
      {BoolExpr::createFalse(), BoolExpr::createFalse()}, selector, BoolExpr::Var(30), {0, 1});
  ASSERT_LT(table.boundaries.size(), size_t{1} << state.size());
  for (size_t id = table.boundaries.size(); id < (size_t{1} << state.size()); ++id) {
    std::unordered_map<size_t, bool> environment{{30, true}};
    assign(environment, state, id);
    assign(environment, selector, 1);
    EXPECT_EQ(evaluateId(encoding, environment), 0u);
    for (auto* expression : encoding.observedNets) EXPECT_FALSE(expression->evaluate(environment));
  }
}

TEST_F(LatchBoundaryEncodingTests, RejectsMismatchedInterfaces) {
  const Network network{3, {0, 1}, {latch("l", 0, 1, 2)}};
  const auto table = compile(network, true);
  const auto state = variables(boundaryEncodingBits(table.boundaries.size()), 2);
  const auto inputs = variables(2, 20);
  const auto selector = variables(2, 30);
  auto* value = BoolExpr::Var(40);
  EXPECT_THROW(encodeBoundaryTable(table, network, {}, inputs, selector, value, {0, 1}), std::invalid_argument);
  EXPECT_THROW(encodeBoundaryTable(table, network, state, {}, selector, value, {}), std::invalid_argument);
  EXPECT_THROW(encodeBoundaryTable(table, network, state, inputs, selector, nullptr, {0, 1}), std::invalid_argument);
  EXPECT_THROW(encodeBoundaryTable(table, network, state, inputs, selector, value, {0}), std::invalid_argument);
}

TEST_F(LatchBoundaryEncodingTests, RejectsAliasingGlobalSelectorIndices) {
  const Network network{3, {0, 1}, {latch("l", 0, 1, 2)}};
  const auto table = compile(network, true);
  const auto state = variables(boundaryEncodingBits(table.boundaries.size()), 2);
  const auto inputs = variables(2, 20);
  EXPECT_THROW(encodeBoundaryTable(table, network, state, inputs, variables(2, 30),
      BoolExpr::Var(40), {1, 1}), std::invalid_argument);
  EXPECT_THROW(encodeBoundaryTable(table, network, state, inputs, variables(1, 30),
      BoolExpr::Var(40), {0, 2}), std::invalid_argument);
}

TEST_F(LatchBoundaryEncodingTests, RejectsMalformedCertifiedRows) {
  const Network network{3, {0, 1}, {latch("l", 0, 1, 2)}};
  const auto table = compile(network, true);
  const auto state = variables(boundaryEncodingBits(table.boundaries.size()), 2);
  const auto inputs = variables(2, 20);
  const auto selector = variables(2, 30);
  auto bad = table;
  bad.rows.front().to = table.boundaries.size();
  EXPECT_THROW(encodeBoundaryTable(bad, network, state, inputs, selector, BoolExpr::Var(40), {0, 1}), std::invalid_argument);
  bad = table;
  bad.rows.front().input.clear();
  EXPECT_THROW(encodeBoundaryTable(bad, network, state, inputs, selector, BoolExpr::Var(40), {0, 1}), std::invalid_argument);
  bad = table;
  const auto& source = bad.boundaries[bad.rows.front().from];
  bad.rows.front().input = {static_cast<uint8_t>(!source.current[0]), static_cast<uint8_t>(!source.current[1])};
  EXPECT_THROW(encodeBoundaryTable(bad, network, state, inputs, selector, BoolExpr::Var(40), {0, 1}), std::invalid_argument);
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
