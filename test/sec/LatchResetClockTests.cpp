// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include <gtest/gtest.h>
#include <algorithm>
#include "BoolExprCache.h"
#include "latch/LatchResetClock.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
class LatchResetClockTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }
  size_t external() {
    const size_t net = network.netCount++;
    network.externalInputs.push_back(net);
    return net;
  }
  size_t constantNet(bool value) {
    const size_t net = network.netCount++;
    network.constantByNet.resize(network.netCount);
    network.constantByNet[net] = value;
    return net;
  }
  size_t add(std::vector<size_t> inputs, ResetClockPrimitive meta) {
    const size_t output = network.netCount++;
    Primitive primitive;
    primitive.name = "cell" + std::to_string(network.primitives.size());
    primitive.inputs = std::move(inputs);
    primitive.outputs = {output};
    primitive.storageBits = meta.kind == ResetClockPrimitive::Kind::Combinational ? 0 : 1;
    network.primitives.push_back(std::move(primitive));
    metadata.push_back(std::move(meta));
    if (!network.constantByNet.empty()) network.constantByNet.resize(network.netCount);
    return output;
  }
  size_t ff(size_t clock, bool inverse = false) {
    auto* expression = BoolExpr::Var(2);
    if (inverse) expression = BoolExpr::Not(expression);
    return add({clock}, {ResetClockPrimitive::Kind::FlipFlop, expression, {}});
  }
  size_t gate(std::vector<size_t> inputs, BoolExpr* output) {
    return add(std::move(inputs), {ResetClockPrimitive::Kind::Combinational, nullptr, {output}});
  }
  ResetClockDiscovery discover() { return discoverResetClock(network, metadata); }
  void expectUnsupported(const std::string& fragment) {
    const auto result = discover();
    EXPECT_EQ(result.status, ResetClockDiscovery::Status::Unsupported) << result.detail;
    EXPECT_FALSE(result.rootNet);
    EXPECT_FALSE(result.externalInputIndex);
    EXPECT_NE(result.detail.find(fragment), std::string::npos) << result.detail;
  }
  Network network;
  std::vector<ResetClockPrimitive> metadata;
};

TEST_F(LatchResetClockTests, DirectPositiveClockKeepsExternalInputIdentity) {
  external();
  const auto clock = external();
  ff(clock);
  const auto result = discover();
  ASSERT_TRUE(result.resolved()) << result.detail;
  EXPECT_EQ(result.rootNet, clock);
  EXPECT_EQ(result.externalInputIndex, 1u);
}

TEST_F(LatchResetClockTests, NegativeEdgeAndMixedEdgesUseOneRoot) {
  const auto clock = external();
  ff(clock, true);
  ff(clock, false);
  const auto result = discover();
  ASSERT_TRUE(result.resolved()) << result.detail;
  EXPECT_EQ(result.rootNet, clock);
}

TEST_F(LatchResetClockTests, BuffersAndInvertersPreserveClockRoot) {
  const auto clock = external();
  const auto inverted = gate({clock}, BoolExpr::Not(BoolExpr::Var(2)));
  const auto buffered = gate({inverted}, BoolExpr::Var(2));
  ff(buffered);
  ff(clock);
  const auto result = discover();
  ASSERT_TRUE(result.resolved()) << result.detail;
  EXPECT_EQ(result.rootNet, clock);
}

TEST_F(LatchResetClockTests, ConstantGateInputsAreFoldedExactly) {
  const auto clock = external();
  const auto one = constantNet(true);
  const auto zero = constantNet(false);
  const auto enabled = gate({clock, one}, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));
  const auto routed = gate({enabled, zero}, BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(3)));
  const auto inverted = gate({routed, one}, BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3)));
  ff(inverted);
  const auto result = discover();
  ASSERT_TRUE(result.resolved()) << result.detail;
  EXPECT_EQ(result.rootNet, clock);
}

TEST_F(LatchResetClockTests, AliasedPinsCollapseToTheSameCarrier) {
  const auto clock = external();
  const auto copy = gate({clock}, BoolExpr::Var(2));
  const auto repeated = gate({copy, clock}, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));
  ff(repeated);
  EXPECT_TRUE(discover().resolved());
}

TEST_F(LatchResetClockTests, LatchEnableDoesNotBecomeAnExtraClockDomain) {
  const auto clock = external();
  const auto enable = external();
  add({enable}, {ResetClockPrimitive::Kind::Latch, BoolExpr::Var(2), {}});
  ff(clock);
  const auto result = discover();
  ASSERT_TRUE(result.resolved()) << result.detail;
  EXPECT_EQ(result.rootNet, clock);
}

TEST_F(LatchResetClockTests, LatchOnlyCircuitHasNoInferredClock) {
  const auto enable = external();
  add({enable}, {ResetClockPrimitive::Kind::Latch, BoolExpr::Var(2), {}});
  network.primitives.back().name = "CLOCK_GATE";
  const auto result = discover();
  EXPECT_EQ(result.status, ResetClockDiscovery::Status::NoEdgeClock);
  EXPECT_FALSE(result.rootNet);
  EXPECT_NE(result.detail.find("latch enables"), std::string::npos);
}

TEST_F(LatchResetClockTests, EmptyCircuitCannotDefineClockCycles) {
  EXPECT_EQ(discover().status, ResetClockDiscovery::Status::NoEdgeClock);
}

TEST_F(LatchResetClockTests, IndependentClockRootsRequireProtocol) {
  ff(external());
  ff(external());
  expectUnsupported("multiple independent");
}

TEST_F(LatchResetClockTests, ArbitraryClockGateDoesNotGuessCarrierFromNames) {
  const auto clock = external();
  const auto gateInput = external();
  const auto gated = gate({clock, gateInput}, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));
  network.primitives.back().name = "CLOCK_GATE";
  ff(gated);
  expectUnsupported("no unique external carrier");
}

TEST_F(LatchResetClockTests, ValidDomainDoesNotHideAnotherUnresolvedClock) {
  const auto clock = external();
  const auto enable = external();
  ff(clock);
  const auto gated = gate({clock, enable}, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));
  ff(gated);
  expectUnsupported("no unique external carrier");
}

TEST_F(LatchResetClockTests, MultiPinClockExpressionCannotGuessRoot) {
  const auto a = external(), b = external();
  add({a, b}, {ResetClockPrimitive::Kind::FlipFlop,
      BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(3)), {}});
  expectUnsupported("no unique external carrier");
}

TEST_F(LatchResetClockTests, StateGeneratedClockIsNotExternalCarrier) {
  const auto clock = external();
  const auto divided = ff(clock);
  ff(divided);
  expectUnsupported("state-generated");
}

TEST_F(LatchResetClockTests, ConstantClockCannotCountResetEdges) {
  ff(constantNet(false));
  expectUnsupported("clock is constant");
}

TEST_F(LatchResetClockTests, DisabledClockGateIsRejectedEvenWithKnownRoot) {
  const auto clock = external();
  const auto zero = constantNet(false);
  ff(gate({clock, zero}, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3))));
  expectUnsupported("clock is constant");
}

TEST_F(LatchResetClockTests, UnrelatedCombinationalCycleDoesNotChangeKnownClock) {
  const auto clock = external();
  ff(clock);
  const size_t first = network.netCount, second = first + 1;
  gate({second}, BoolExpr::Var(2));
  gate({first}, BoolExpr::Var(2));
  EXPECT_TRUE(discover().resolved());
}

TEST_F(LatchResetClockTests, CyclicClockRoutingTerminatesWithoutGuessing) {
  const size_t first = network.netCount, second = first + 1;
  gate({second}, BoolExpr::Var(2));
  gate({first}, BoolExpr::Var(2));
  ff(first);
  expectUnsupported("cyclic");
}

TEST_F(LatchResetClockTests, ReverseOrderedDeepRoutingDoesNotRecurse) {
  const auto clock = external();
  size_t routed = clock;
  for (size_t i = 0; i < 4096; ++i) routed = gate({routed}, BoolExpr::Not(BoolExpr::Var(2)));
  std::reverse(network.primitives.begin(), network.primitives.end());
  std::reverse(metadata.begin(), metadata.end());
  ff(routed);
  const auto result = discover();
  ASSERT_TRUE(result.resolved()) << result.detail;
  EXPECT_EQ(result.rootNet, clock);
}

TEST_F(LatchResetClockTests, MissingMetadataDoesNotSilentlyDropAClockDomain) {
  ff(external());
  metadata.clear();
  expectUnsupported("metadata for every primitive");
}

TEST_F(LatchResetClockTests, UnknownPrimitiveCannotBeAssumedCombinational) {
  ff(external());
  add({}, {});
  expectUnsupported("unclassified primitive");
}

TEST_F(LatchResetClockTests, MissingClockExpressionIsUnsupported) {
  ff(external());
  metadata.back().clock = nullptr;
  expectUnsupported("invalid clock routing expression");
}

TEST_F(LatchResetClockTests, ClockExpressionMayNotReferenceUndeclaredPins) {
  ff(external());
  metadata.back().clock = BoolExpr::Var(3);
  expectUnsupported("non-input symbol");
}

TEST_F(LatchResetClockTests, IncompleteGateMetadataIsRejected) {
  const auto clock = external();
  const auto buffered = gate({clock}, BoolExpr::Var(2));
  metadata.back().outputs.clear();
  ff(buffered);
  expectUnsupported("incomplete combinational");
}

TEST_F(LatchResetClockTests, DuplicateOrDrivenExternalInputsAreRejected) {
  const auto clock = external();
  ff(clock);
  network.externalInputs.push_back(clock);
  expectUnsupported("duplicate external");
  network.externalInputs.pop_back();
  network.primitives[0].outputs[0] = clock;
  expectUnsupported("multiple drivers");
}

TEST_F(LatchResetClockTests, InvalidConstantsAndPinRangesAreRejected) {
  const auto clock = external();
  ff(clock);
  network.constantByNet = {false};
  expectUnsupported("invalid constant net table");
  network.constantByNet.resize(network.netCount);
  expectUnsupported("also a constant");
  network.constantByNet.clear();
  network.primitives[0].inputs[0] = network.netCount;
  expectUnsupported("out-of-range input");
}
}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
