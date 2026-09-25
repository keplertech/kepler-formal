// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include <gtest/gtest.h>
#include <set>
#include <unordered_map>

#include "BoolExprCache.h"
#include "latch/LatchSymbolicModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using Environment = std::unordered_map<size_t, bool>;

class LatchSymbolicModelTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }
  SymbolicNetwork lifted(Network network) {
    SymbolicNetwork result;
    result.reference = std::move(network);
    for (const auto& primitive : result.reference.primitives) result.primitives.push_back(liftPrimitive(primitive));
    return result;
  }
  SymbolicBits constants(const Bits& bits) {
    SymbolicBits result;
    for (auto bit : bits) result.push_back(bit ? BoolExpr::createTrue() : BoolExpr::createFalse());
    return result;
  }
  SymbolicState constants(const State& state) {
    SymbolicState result;
    result.current = constants(state.current); result.previous = constants(state.previous);
    result.active = constants(state.active);
    for (const auto& storage : state.storage) result.storage.push_back(constants(storage));
    result.bootstrap = BoolExpr::Var(state.bootstrap); result.error = BoolExpr::Var(state.error);
    return result;
  }
  SymbolicBits variables(size_t count, size_t& next) {
    SymbolicBits result;
    for (size_t i = 0; i < count; ++i) result.push_back(BoolExpr::Var(next++));
    return result;
  }
  Bits evaluate(const SymbolicBits& bits, const Environment& env) {
    Bits result;
    for (auto* bit : bits) result.push_back(bit->evaluate(env));
    return result;
  }
  State evaluate(const SymbolicState& state, const Environment& env) {
    State result;
    result.current = evaluate(state.current, env); result.previous = evaluate(state.previous, env);
    result.active = evaluate(state.active, env);
    for (const auto& storage : state.storage) result.storage.push_back(evaluate(storage, env));
    result.bootstrap = state.bootstrap->evaluate(env); result.error = state.error->evaluate(env);
    return result;
  }
  State withoutReason(State state) { state.errorReason.clear(); return state; }
  Environment assignment(size_t first, size_t count, size_t value, Environment base = {}) {
    for (size_t bit = 0; bit < count; ++bit) base[first + bit] = (value >> bit) & 1;
    return base;
  }
  void expectSuccessors(const EventModel& concrete, const SymbolicEventModel& symbolic, const State& state) {
    size_t next = 1000;
    const auto formula = symbolic.wave(constants(state), [&] { return BoolExpr::Var(next++); });
    ASSERT_LT(next - 1000, 12u);
    std::set<State> actual, expected;
    for (size_t value = 0; value < (size_t{1} << (next - 1000)); ++value)
      actual.insert(evaluate(formula, assignment(1000, next - 1000, value)));
    for (const auto& successor : concrete.successors(state)) expected.insert(withoutReason(successor));
    EXPECT_EQ(actual, expected);
  }
  State settle(const EventModel& model, State state) {
    for (size_t i = 0; i < 32 && !model.stable(state); ++i) {
      const auto next = model.successors(state);
      EXPECT_EQ(next.size(), 1u);
      state = next.front();
    }
    EXPECT_TRUE(model.stable(state));
    return state;
  }
};

TEST_F(LatchSymbolicModelTests, ExhaustiveLatchWaveMatchesEveryCompleteConcreteState) {
  const auto network = lifted({3, {0, 1}, {latch("l", 0, 1, 2)}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  size_t next = 2;
  SymbolicState state;
  state.current = variables(3, next); state.previous = variables(3, next);
  state.storage = {variables(1, next)}; state.active = variables(1, next);
  state.bootstrap = variables(1, next).front(); state.error = variables(1, next).front();
  size_t choice = 1000;
  const auto result = symbolic.wave(state, [&] { return BoolExpr::Var(choice++); });
  for (size_t value = 0; value < (size_t{1} << (next - 2)); ++value) {
    const auto env = assignment(2, next - 2, value);
    const auto concreteState = evaluate(state, env);
    std::set<State> actual, expected;
    for (size_t order = 0; order < (size_t{1} << (choice - 1000)); ++order)
      actual.insert(evaluate(result, assignment(1000, choice - 1000, order, env)));
    for (const auto& successor : concrete.successors(concreteState)) expected.insert(withoutReason(successor));
    EXPECT_EQ(actual, expected) << "assignment=" << value;
  }
}

TEST_F(LatchSymbolicModelTests, ExhaustiveFlipFlopWaveConsumesClockOnlyOnClockVisit) {
  const auto network = lifted({3, {0, 1}, {flipFlop("f", 0, 1, 2)}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  for (size_t value = 0; value < 128; ++value) {
    State state;
    state.current = {uint8_t(value & 1), uint8_t((value >> 1) & 1), uint8_t((value >> 2) & 1)};
    state.previous = {uint8_t((value >> 3) & 1), uint8_t((value >> 4) & 1), uint8_t((value >> 5) & 1)};
    state.storage = {{uint8_t((value >> 6) & 1)}};
    state.active = {1};
    expectSuccessors(concrete, symbolic, state);
  }
}

TEST_F(LatchSymbolicModelTests, BootstrapMatchesAllInitialInputsStorageAndSeeds) {
  const auto network = lifted({3, {0, 1}, {latch("l", 0, 1, 2)}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  size_t next = 2;
  const auto inputs = variables(2, next), storage = variables(1, next), seeds = variables(3, next);
  const auto boot = symbolic.bootstrap(inputs, {storage}, seeds);
  for (size_t value = 0; value < 64; ++value) {
    const auto env = assignment(2, 6, value);
    const auto expected = concrete.bootstrap(evaluate(inputs, env), {evaluate(storage, env)}, evaluate(seeds, env));
    EXPECT_EQ(evaluate(boot, env), withoutReason(expected));
    expectSuccessors(concrete, symbolic, expected);
  }
}

TEST_F(LatchSymbolicModelTests, InputDependentBootstrapUsesFrozenSeedSnapshot) {
  auto first = latch("a", 0, 1, 2), second = latch("b", 2, 1, 3);
  first.initialOutputValues = second.initialOutputValues = [](const Bits& storage, const Bits& pins) {
    return Bits{static_cast<uint8_t>(storage[0] & pins[0])};
  };
  const auto network = lifted({4, {0, 1}, {first, second}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  for (uint8_t seed : {0, 1}) {
    const auto expected = concrete.bootstrap({1, 0}, {{1}, {1}}, {0, 0, seed, 0});
    const auto result = symbolic.bootstrap(constants({1, 0}), {constants({1}), constants({1})}, constants({0, 0, seed, 0}));
    EXPECT_EQ(evaluate(result, {}), withoutReason(expected));
    EXPECT_EQ(expected.current[3], seed);
  }
}

TEST_F(LatchSymbolicModelTests, AdmissionMatchesAllTransactionsIncludingNonstableErrors) {
  const auto network = lifted({3, {0, 1}, {latch("l", 0, 1, 2)}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  for (size_t value = 0; value < 8; ++value) {
    auto state = settle(concrete, concrete.bootstrap({uint8_t(value & 1), uint8_t((value >> 1) & 1)},
        {{uint8_t((value >> 2) & 1)}}, Bits(3)));
    for (size_t tx = 0; tx < 4; ++tx) {
      const Bits inputs{uint8_t(tx & 1), uint8_t((tx >> 1) & 1)};
      EXPECT_EQ(evaluate(symbolic.admit(constants(state), constants(inputs)), {}), withoutReason(concrete.admit(state, inputs)));
    }
    state.active = {1};
    EXPECT_EQ(evaluate(symbolic.admit(constants(state), constants({0, 0})), {}), withoutReason(concrete.admit(state, {0, 0})));
  }
}

TEST_F(LatchSymbolicModelTests, FourOpenLatchesAndSelfFeedbackPreserveHistory) {
  const auto network = lifted({6, {0, 1}, {latch("a", 0, 1, 2), latch("b", 2, 1, 3),
      latch("c", 3, 1, 4), latch("d", 4, 1, 5)}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network, {}, 4);
  auto state = settle(concrete, concrete.bootstrap({0, 1}, {{0}, {0}, {0}, {0}}, Bits(6)));
  state = concrete.admit(state, {1, 1});
  for (size_t i = 0; i < 5; ++i) {
    expectSuccessors(concrete, symbolic, state);
    state = concrete.successors(state).front();
  }
  EXPECT_TRUE(concrete.stable(state));
  EXPECT_EQ(state.current[5], 1);
  const auto self = lifted({2, {0}, {latch("self", 1, 0, 1)}});
  EventModel concreteSelf(self.reference);
  SymbolicEventModel symbolicSelf(self);
  for (uint8_t bit : {0, 1}) expectSuccessors(concreteSelf, symbolicSelf, concreteSelf.bootstrap({1}, {{bit}}, Bits(2)));
}

TEST_F(LatchSymbolicModelTests, IntermediateEnablePulseSurvivesEveryWave) {
  Network reference{6, {0}, {
    combinational("a", {0}, {1}, [](const Bits& bits) { return bits; }),
    combinational("b", {1}, {2}, [](const Bits& bits) { return bits; }),
    combinational("xor", {0, 2}, {3}, [](const Bits& bits) { return Bits{uint8_t(bits[0] ^ bits[1])}; }),
    latch("l", 4, 3, 5)}};
  reference.constantByNet.resize(6); reference.constantByNet[4] = true;
  const auto network = lifted(reference);
  EventModel concrete(reference);
  SymbolicEventModel symbolic(network);
  auto state = settle(concrete, concrete.bootstrap({0}, {{}, {}, {}, {0}}, Bits(6)));
  state = concrete.admit(state, {1});
  bool pulse = false;
  for (size_t i = 0; i < 8 && !concrete.stable(state); ++i) {
    expectSuccessors(concrete, symbolic, state);
    state = concrete.successors(state).front(); pulse |= state.current[3];
  }
  EXPECT_TRUE(pulse); EXPECT_EQ(state.current[3], 0); EXPECT_EQ(state.current[5], 1);
}

TEST_F(LatchSymbolicModelTests, SharedFanoutChoiceAndWorkerOrderingAreIdentical) {
  const auto network = lifted({5, {0, 1}, {flipFlop("f", 0, 1, 2),
    combinational("a", {2}, {3}, [](const Bits& bits) { return bits; }),
    combinational("b", {2}, {4}, [](const Bits& bits) { return bits; })}});
  EventModel concrete(network.reference);
  SymbolicEventModel serial(network, {}, 1), parallel(network, {}, 4);
  auto initial = settle(concrete, concrete.bootstrap({0, 0}, {{0}, {}, {}}, Bits(5)));
  auto state = constants(concrete.admit(initial, {1, 1}));
  size_t first = 100, second = 100;
  auto a = serial.wave(state, [&] { return BoolExpr::Var(first++); });
  auto b = parallel.wave(state, [&] { return BoolExpr::Var(second++); });
  EXPECT_EQ(first, second);
  EXPECT_EQ(flattenSymbolicBoundary(a), flattenSymbolicBoundary(b));
  a = serial.wave(a, [&] { return BoolExpr::Var(first++); });
  b = parallel.wave(b, [&] { return BoolExpr::Var(second++); });
  EXPECT_EQ(flattenSymbolicBoundary(a), flattenSymbolicBoundary(b));
  for (size_t choice = 0; choice < (size_t{1} << (first - 100)); ++choice) {
    const auto result = evaluate(a, assignment(100, first - 100, choice));
    EXPECT_EQ(result.current[2], result.current[3]); EXPECT_EQ(result.current[3], result.current[4]);
  }
}

TEST_F(LatchSymbolicModelTests, UnusedChoiceEncodingsSelectLegalOrders) {
  auto primitive = flipFlop("f", 0, 1, 3);
  primitive.inputs.push_back(2);  // Six total permutations use three choice bits.
  const auto network = lifted({4, {0, 1, 2}, {primitive}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  auto state = settle(concrete, concrete.bootstrap({0, 0, 0}, {{0}}, Bits(4)));
  state = concrete.admit(state, {1, 1, 1});
  expectSuccessors(concrete, symbolic, state);
}

TEST_F(LatchSymbolicModelTests, StickyErrorsRemainVisibleWithSuccessfulAlternative) {
  auto primitive = latch("l", 0, 1, 2);
  primitive.react = [](const Bits& storage, const Bits&, const Bits& pins, std::optional<size_t>, bool) -> Reaction {
    if (pins[0] && pins[1]) return {{}, {}, true, "conflict"};
    const Bits next{pins[0] ? uint8_t{0} : pins[1] ? uint8_t{1} : storage[0]};
    return {next, next};
  };
  const auto network = lifted({3, {0, 1}, {primitive}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  auto state = settle(concrete, concrete.bootstrap({1, 0}, {{0}}, Bits(3)));
  state = concrete.admit(state, {0, 1});
  expectSuccessors(concrete, symbolic, state);
  for (const auto& outcome : concrete.successors(state)) {
    expectSuccessors(concrete, symbolic, outcome);
    EXPECT_EQ(symbolic.stable(constants(outcome))->evaluate({}), concrete.stable(outcome));
  }
}

TEST_F(LatchSymbolicModelTests, BoundaryInvariantChecksConstantsLevelsAndNoFakeFfCapture) {
  const auto network = lifted({3, {0, 1}, {latch("l", 0, 1, 2)}});
  SymbolicEventModel latchModel(network);
  EXPECT_TRUE(latchModel.boundaryInvariant(latchModel.boundary(constants({0, 0, 1}), {constants({1})}))->evaluate({}));
  EXPECT_FALSE(latchModel.boundaryInvariant(latchModel.boundary(constants({0, 1, 1}), {constants({1})}))->evaluate({}));
  EXPECT_FALSE(latchModel.boundaryInvariant(latchModel.boundary(constants({0, 0, 0}), {constants({1})}))->evaluate({}));
  SymbolicEventModel flopModel(lifted({3, {0, 1}, {flipFlop("f", 0, 1, 2)}}));
  EXPECT_TRUE(flopModel.boundaryInvariant(flopModel.boundary(constants({1, 1, 0}), {constants({0})}))->evaluate({}));
  Network reference{1, {}, {}, {true}};
  SymbolicEventModel constantModel(lifted(reference));
  EXPECT_FALSE(constantModel.boundaryInvariant(constantModel.boundary(constants({0}), {}))->evaluate({}));
  EXPECT_TRUE(constantModel.boundaryInvariant(constantModel.boundary(constants({1}), {}))->evaluate({}));
}

TEST_F(LatchSymbolicModelTests, ResourceAndMalformedCallbackFailuresNeverPruneOutcomes) {
  const auto network = lifted({3, {0, 1}, {flipFlop("f", 0, 1, 2)}});
  SymbolicEventModel limited(network, {1, 100});
  EventModel concrete(network.reference);
  const auto state = concrete.admit(settle(concrete, concrete.bootstrap({0, 0}, {{0}}, Bits(3))), {1, 1});
  size_t next = 100;
  EXPECT_THROW(limited.wave(constants(state), [&] { return BoolExpr::Var(next++); }), Limit);
  EXPECT_THROW(liftPrimitive(network.reference.primitives[0], 1), Limit);
  auto malformed = network;
  malformed.primitives[0].react = [](const SymbolicBits&, const SymbolicBits&, const SymbolicBits&,
      std::optional<size_t>, bool) { return SymbolicReaction{}; };
  SymbolicEventModel bad(malformed);
  const auto result = bad.wave(constants(state), [&] { return BoolExpr::Var(next++); });
  EXPECT_TRUE(result.error->evaluate({}));
  EXPECT_FALSE(bad.stable(result)->evaluate({}));
  EXPECT_EQ(evaluate(bad.wave(result, {}), {}), evaluate(result, {}));
}

TEST_F(LatchSymbolicModelTests, FlattenedBoundaryIncludesAllNetAndStorageBitsInFixedOrder) {
  SymbolicEventModel model(lifted({3, {0, 1}, {latch("l", 0, 1, 2)}}));
  size_t next = 2;
  const auto nets = variables(3, next), storage = variables(1, next);
  const auto boundary = model.boundary(nets, {storage});
  EXPECT_EQ(flattenSymbolicBoundary(boundary), (SymbolicBits{nets[0], nets[1], nets[2], storage[0]}));
  EXPECT_EQ(boundary.previous, nets);
  EXPECT_TRUE(model.stable(boundary)->evaluate({}));
}

TEST_F(LatchSymbolicModelTests, MissingExplicitInitialProjectionCannotFallBackToWrongIdentity) {
  auto primitive = latch("complement", 0, 1, 2);
  primitive.initialOutputs = [](const Bits& storage) { return Bits{uint8_t(!storage[0])}; };
  auto network = lifted({3, {0, 1}, {primitive}});
  network.primitives.front().initialOutputs = {};
  SymbolicEventModel model(network);
  const auto state = model.bootstrap(constants({0, 0}), {constants({0})}, constants({0, 0, 0}));
  EXPECT_TRUE(state.error->evaluate({}));
  EXPECT_EQ(evaluate(model.wave(state, {}), {}), evaluate(state, {}));
}

TEST_F(LatchSymbolicModelTests, SymbolicBootstrapDoesNotFabricateFlopStartupEdge) {
  const auto network = lifted({3, {0, 1}, {flipFlop("f", 0, 1, 2)}});
  SymbolicEventModel model(network);
  auto state = model.bootstrap(constants({1, 1}), {constants({0})}, constants({0, 0, 1}));
  state = model.wave(state, {});
  EXPECT_TRUE(model.stable(state)->evaluate({}));
  EXPECT_FALSE(state.current[2]->evaluate({}));
}

TEST_F(LatchSymbolicModelTests, ZeroChangedActivePrimitiveUsesNullVisitFallback) {
  auto primitive = latch("null_visit", 0, 1, 2);
  primitive.react = [](const Bits& storage, const Bits&, const Bits&,
      std::optional<size_t> changed, bool bootstrap) {
    const Bits next{uint8_t(!bootstrap && !changed ? !storage[0] : storage[0])};
    return Reaction{next, next};
  };
  const auto network = lifted({3, {0, 1}, {primitive}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  State state{{0, 0, 0}, {0, 0, 0}, {{0}}, {1}, false, false, {}};
  expectSuccessors(concrete, symbolic, state);
}

TEST_F(LatchSymbolicModelTests, PureCombinationalSnapshotsNeverInventXorPulse) {
  const auto network = lifted({4, {0}, {
    combinational("a", {0}, {1}, [](const Bits& bits) { return bits; }),
    combinational("b", {0}, {2}, [](const Bits& bits) { return bits; }),
    combinational("xor", {1, 2}, {3}, [](const Bits& bits) { return Bits{uint8_t(bits[0] ^ bits[1])}; })}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network, {}, 4);
  auto state = settle(concrete, concrete.bootstrap({0}, {{}, {}, {}}, Bits(4)));
  state = concrete.admit(state, {1});
  for (size_t i = 0; i < 3; ++i) {
    expectSuccessors(concrete, symbolic, state);
    state = concrete.successors(state).front();
    EXPECT_EQ(state.current[3], 0);
  }
}

TEST_F(LatchSymbolicModelTests, IndependentPrimitiveChoicesAreAllocatedSeriallyAndRemainIndependent) {
  const auto network = lifted({4, {0, 1}, {flipFlop("a", 0, 1, 2), flipFlop("b", 0, 1, 3)}});
  EventModel concrete(network.reference);
  SymbolicEventModel serial(network, {}, 1), parallel(network, {}, 4);
  const auto initial = settle(concrete, concrete.bootstrap({0, 0}, {{0}, {0}}, Bits(4)));
  const auto state = concrete.admit(initial, {1, 1});
  size_t first = 100, second = 100;
  const auto a = serial.wave(constants(state), [&] { return BoolExpr::Var(first++); });
  const auto b = parallel.wave(constants(state), [&] { return BoolExpr::Var(second++); });
  ASSERT_EQ(first, 102u); ASSERT_EQ(second, first);
  EXPECT_EQ(flattenSymbolicBoundary(a), flattenSymbolicBoundary(b));
  std::set<State> outcomes;
  for (size_t bits = 0; bits < 4; ++bits) outcomes.insert(evaluate(a, assignment(100, 2, bits)));
  EXPECT_EQ(outcomes.size(), 4u);
  expectSuccessors(concrete, parallel, state);
}

TEST_F(LatchSymbolicModelTests, FixedChoiceModeSelectsLegalConcreteOrderWithoutFreeSymbols) {
  const auto network = lifted({3, {0, 1}, {flipFlop("f", 0, 1, 2)}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  const auto state = concrete.admit(settle(concrete, concrete.bootstrap({0, 0}, {{0}}, Bits(3))), {1, 1});
  std::set<State> permitted;
  for (auto successor : concrete.successors(state)) permitted.insert(withoutReason(successor));
  for (bool value : {false, true}) {
    const auto fixed = symbolic.wave(constants(state), [value] { return BoolExpr::Var(value); });
    EXPECT_TRUE(permitted.count(evaluate(fixed, {})));
  }
}

TEST_F(LatchSymbolicModelTests, SymbolicStableValuationsPadEveryStateFieldByIdentity) {
  SymbolicEventModel model(lifted({4, {0, 1},
      {latch("l", 0, 1, 2), flipFlop("f", 0, 1, 3)}}));
  size_t next = 2;
  SymbolicState state;
  state.current = variables(4, next);
  state.previous = variables(4, next);
  state.storage = {variables(1, next), variables(1, next)};
  state.active = variables(2, next);
  state.bootstrap = variables(1, next).front();
  state.error = variables(1, next).front();
  ASSERT_NE(model.stable(state), BoolExpr::createTrue());
  size_t choice = 1000;
  const auto result = model.wave(state, [&] { return BoolExpr::Var(choice++); });
  for (size_t value = 0; value < 64; ++value) {
    Environment env;
    for (size_t bit = 0; bit < 4; ++bit) {
      env[state.current[bit]->getId()] = (value >> bit) & 1;
      env[state.previous[bit]->getId()] = (value >> bit) & 1;
    }
    env[state.storage[0][0]->getId()] = (value >> 4) & 1;
    env[state.storage[1][0]->getId()] = (value >> 5) & 1;
    for (auto* active : state.active) env[active->getId()] = false;
    env[state.bootstrap->getId()] = false;
    env[state.error->getId()] = false;
    ASSERT_TRUE(model.stable(state)->evaluate(env));
    for (size_t order = 0; order < (size_t{1} << (choice - 1000)); ++order)
      EXPECT_EQ(evaluate(result, assignment(1000, choice - 1000, order, env)), evaluate(state, env));
  }
}

TEST_F(LatchSymbolicModelTests, ConsumedAndSinkHistoryMatchesConcreteAdmissionIncludingErrors) {
  const auto network = lifted({5, {0, 1, 2},
      {latch("first", 0, 1, 3), latch("second", 3, 1, 4)}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  size_t next = 2;
  auto state = symbolic.boundary(variables(5, next), {variables(1, next), variables(1, next)});
  const auto input = variables(3, next);
  state.error = variables(1, next).front();
  const auto admitted = symbolic.admit(state, input);
  for (size_t value = 0; value < (size_t{1} << (next - 2)); ++value) {
    const auto env = assignment(2, next - 2, value);
    EXPECT_EQ(evaluate(admitted, env), withoutReason(concrete.admit(evaluate(state, env), evaluate(input, env))));
  }
}

TEST_F(LatchSymbolicModelTests, ConsumedAndSinkHistoryMatchesConcreteWavesIncludingErrors) {
  const auto network = lifted({5, {0, 1, 2},
      {latch("first", 0, 1, 3), latch("second", 3, 1, 4)}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  for (size_t value = 0; value < 512; ++value) {
    State state;
    for (size_t bit = 0; bit < 5; ++bit) state.current.push_back((value >> bit) & 1);
    state.previous = state.current;
    // Include stale sink history, changed consumed feedback and a changed
    // unconsumed external input. Supplied activation need not be consistent:
    // the optimization applies only to the newly computed activation set.
    for (auto net : {2u, 3u, 4u}) state.previous[net] ^= 1;
    state.storage = {{uint8_t((value >> 5) & 1)}, {uint8_t((value >> 6) & 1)}};
    state.active = {uint8_t((value >> 7) & 1), uint8_t((value >> 8) & 1)};
    expectSuccessors(concrete, symbolic, state);
    state.error = true;
    expectSuccessors(concrete, symbolic, state);
  }
}

TEST_F(LatchSymbolicModelTests, EntirelyUnconsumedNetworkStillNormalizesHistory) {
  const auto network = lifted({1, {0}, {}});
  EventModel concrete(network.reference);
  SymbolicEventModel symbolic(network);
  auto state = settle(concrete, concrete.bootstrap({0}, {}, {0}));
  const auto result = symbolic.admit(constants(state), constants({1}));
  EXPECT_EQ(evaluate(result, {}), withoutReason(concrete.admit(state, {1})));
  EXPECT_EQ(evaluate(result.previous, {}), Bits{1});
  state.previous = {1};
  expectSuccessors(concrete, symbolic, state);
}
}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
