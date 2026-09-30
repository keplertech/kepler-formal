// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <algorithm>
#include <set>

#include "latch/LatchEventModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

Primitive buffer(size_t input, size_t output) {
  return combinational("buffer", {input}, {output}, [](const Bits& bits) { return bits; });
}

Primitive inverter(size_t input, size_t output) {
  return combinational("inverter", {input}, {output}, [](const Bits& bits) {
    return Bits{static_cast<uint8_t>(!bits[0])};
  });
}

std::vector<State> settle(const EventModel& model, State state, size_t bound = 32) {
  std::vector<State> states{std::move(state)};
  for (size_t wave = 0; wave <= bound; ++wave) {
    if (std::all_of(states.begin(), states.end(), [&](const auto& item) { return model.stable(item); })) {
      return states;
    }
    std::set<State> next;
    for (const auto& item : states) {
      const auto successors = model.successors(item);
      next.insert(successors.begin(), successors.end());
    }
    states.assign(next.begin(), next.end());
  }
  ADD_FAILURE() << "Expected bounded settling";
  return states;
}

State simpleBoundary(const EventModel& model, Bits inputs, std::vector<Bits> storage) {
  const auto states = settle(model, model.bootstrap(inputs, storage, Bits(model.network().netCount)));
  EXPECT_EQ(states.size(), 1);
  return states.front();
}

TEST(LatchEventModelTests, EmptyNetworkHasOneForcedBootstrapWave) {
  EventModel model(Network{});
  const auto initial = model.bootstrap({}, {}, {});
  EXPECT_FALSE(model.stable(initial));
  const auto next = model.successors(initial);
  ASSERT_EQ(next.size(), 1);
  EXPECT_TRUE(model.stable(next.front()));
  EXPECT_EQ(model.successors(next.front()), next);
}

TEST(LatchEventModelTests, BootstrapCopiesInitialStorageBeforeForcedEvaluation) {
  EventModel model({3, {0, 1}, {latch("l", 0, 1, 2)}});
  auto initial = model.bootstrap({0, 0}, {{1}}, {1, 1, 0});
  EXPECT_EQ(initial.current, (Bits{0, 0, 1}));
  EXPECT_EQ(initial.previous, initial.current);
  EXPECT_EQ(settle(model, initial).front().storage[0], Bits{1});
}

TEST(LatchEventModelTests, LatchFollowsDataWhileOpenAndHoldsWhenClosed) {
  EventModel model({3, {0, 1}, {latch("l", 0, 1, 2)}});
  auto state = simpleBoundary(model, {0, 1}, {{0}});
  state = settle(model, model.admit(state, {1, 1})).front();
  EXPECT_EQ(state.current[2], 1);
  state = settle(model, model.admit(state, {1, 0})).front();
  state = settle(model, model.admit(state, {0, 0})).front();
  EXPECT_EQ(state.current[2], 1);
}

TEST(LatchEventModelTests, ActiveLowLatchUsesExplicitPolarity) {
  EventModel model({3, {0, 1}, {latch("l", 0, 1, 2, false)}});
  auto state = simpleBoundary(model, {1, 0}, {{0}});
  EXPECT_EQ(state.storage[0], Bits{1});
  state = settle(model, model.admit(state, {1, 1})).front();
  state = settle(model, model.admit(state, {0, 1})).front();
  EXPECT_EQ(state.storage[0], Bits{1});
}

TEST(LatchEventModelTests, FourOpenLatchesPropagateInOneExternalEpisode) {
  EventModel model({6, {0, 1}, {latch("l0", 0, 1, 2), latch("l1", 2, 1, 3),
                               latch("l2", 3, 1, 4), latch("l3", 4, 1, 5)}});
  auto state = simpleBoundary(model, {0, 1}, {{0}, {0}, {0}, {0}});
  state = model.admit(state, {1, 1});
  for (size_t wave = 0; wave < 4; ++wave) {
    ASSERT_FALSE(model.stable(state));
    auto next = model.successors(state);
    ASSERT_EQ(next.size(), 1);
    state = next.front();
    for (size_t output = 0; output < 4; ++output) {
      EXPECT_EQ(state.current[output + 2], output <= wave);
    }
  }
  EXPECT_TRUE(model.stable(state));
}

TEST(LatchEventModelTests, SimultaneousCloseAndDataPreservesBothPinOrders) {
  EventModel model({3, {0, 1}, {latch("l", 0, 1, 2)}});
  const auto state = simpleBoundary(model, {0, 1}, {{0}});
  const auto results = settle(model, model.admit(state, {1, 0}));
  ASSERT_EQ(results.size(), 2);
  EXPECT_EQ(results[0].current[2], 0);
  EXPECT_EQ(results[1].current[2], 1);
}

TEST(LatchEventModelTests, OpeningWithDataChangeDeduplicatesIdenticalFinalOutcomes) {
  EventModel model({3, {0, 1}, {latch("l", 0, 1, 2)}});
  const auto state = simpleBoundary(model, {0, 0}, {{0}});
  const auto results = settle(model, model.admit(state, {1, 1}));
  ASSERT_EQ(results.size(), 1);
  EXPECT_EQ(results.front().current[2], 1);
}

TEST(LatchEventModelTests, BootstrapDoesNotInventClockEdge) {
  EventModel model({3, {0, 1}, {flipFlop("f", 0, 1, 2)}});
  const auto state = simpleBoundary(model, {1, 1}, {{0}});
  EXPECT_EQ(state.storage[0], Bits{0});
}

TEST(LatchEventModelTests, FlipFlopCapturesOnlySpecifiedEdge) {
  EventModel model({3, {0, 1}, {flipFlop("f", 0, 1, 2)}});
  auto state = simpleBoundary(model, {1, 0}, {{0}});
  state = settle(model, model.admit(state, {1, 1})).front();
  EXPECT_EQ(state.current[2], 1);
  state = settle(model, model.admit(state, {0, 1})).front();
  state = settle(model, model.admit(state, {0, 0})).front();
  EXPECT_EQ(state.current[2], 1);
}

TEST(LatchEventModelTests, FallingEdgeFlipFlopUsesExplicitPolarity) {
  EventModel model({3, {0, 1}, {flipFlop("f", 0, 1, 2, false)}});
  auto state = simpleBoundary(model, {1, 1}, {{0}});
  state = settle(model, model.admit(state, {1, 0})).front();
  EXPECT_EQ(state.current[2], 1);
}

TEST(LatchEventModelTests, ClockDataRaceDoesNotConsumeClockEdgeTwice) {
  EventModel model({3, {0, 1}, {flipFlop("f", 0, 1, 2)}});
  const auto state = simpleBoundary(model, {0, 0}, {{0}});
  const auto results = settle(model, model.admit(state, {1, 1}));
  // Clock first keeps old data; a subsequent data-pin visit MUST NOT re-use it.
  ASSERT_EQ(results.size(), 2);
  EXPECT_EQ(results[0].current[2], 0);
  EXPECT_EQ(results[1].current[2], 1);
}

TEST(LatchEventModelTests, GeneratedClockAfterBootstrapCommitIsRealEvent) {
  EventModel model({4, {0, 1}, {buffer(1, 2), flipFlop("f", 0, 2, 3)}});
  auto states = settle(model, model.bootstrap({1, 1}, {{}, {0}}, {0, 0, 0, 0}));
  ASSERT_EQ(states.size(), 1);
  EXPECT_EQ(states[0].current[3], 1);
}

TEST(LatchEventModelTests, GeneratedClockBootstrapSeedsCanChangeRememberedOutcome) {
  EventModel model({4, {0, 1}, {buffer(1, 2), flipFlop("f", 0, 2, 3)}});
  const auto lowSeed = settle(model, model.bootstrap({1, 1}, {{}, {0}}, {0, 0, 0, 0}));
  const auto highSeed = settle(model, model.bootstrap({1, 1}, {{}, {0}}, {0, 0, 1, 0}));
  EXPECT_EQ(lowSeed.front().current[3], 1);
  EXPECT_EQ(highSeed.front().current[3], 0);
}

TEST(LatchEventModelTests, SelfFeedbackRetainsHistoryInsteadOfChoosingFixedPoint) {
  EventModel model({2, {0}, {latch("self", 1, 0, 1)}});
  for (uint8_t value : {0, 1}) {
    const auto state = simpleBoundary(model, {1}, {{value}});
    EXPECT_EQ(state.current[1], value);
    EXPECT_EQ(model.successors(state), std::vector<State>{state});
  }
}

TEST(LatchEventModelTests, InvertingTransparentFeedbackKeepsExecuting) {
  EventModel model({3, {0}, {inverter(2, 1), latch("loop", 1, 0, 2)}});
  auto state = model.bootstrap({1}, {{}, {0}}, {0, 1, 0});
  std::set<State> seen;
  bool repeated = false;
  for (size_t i = 0; i < 20; ++i) {
    EXPECT_FALSE(model.stable(state));
    if (!seen.insert(state).second) repeated = true;
    auto next = model.successors(state);
    ASSERT_EQ(next.size(), 1);
    state = next.front();
  }
  EXPECT_TRUE(repeated);
}

TEST(LatchEventModelTests, ClosedLatchBreaksInvertingFeedback) {
  EventModel model({3, {0}, {inverter(2, 1), latch("loop", 1, 0, 2)}});
  const auto state = simpleBoundary(model, {0}, {{}, {0}});
  EXPECT_EQ(state.current, (Bits{0, 1, 0}));
}

TEST(LatchEventModelTests, IntermediateEnablePulseIsNotLostAtRegionBoundary) {
  Network network{6, {0}, {buffer(0, 1), buffer(1, 2),
      combinational("xor", {0, 2}, {3}, [](const Bits& bits) {
        return Bits{static_cast<uint8_t>(bits[0] ^ bits[1])};
      }), latch("capture", 4, 3, 5)}};
  network.constantByNet.resize(6);
  network.constantByNet[4] = true;
  EventModel model(network);
  auto state = simpleBoundary(model, {0}, {{}, {}, {}, {0}});
  state = model.admit(state, {1});
  bool sawPulse = false;
  for (size_t i = 0; !model.stable(state) && i < 16; ++i) {
    state = model.successors(state).front();
    sawPulse |= state.current[3] != 0;
  }
  EXPECT_TRUE(sawPulse);
  EXPECT_TRUE(model.stable(state));
  EXPECT_EQ(state.current[3], 0);
  EXPECT_EQ(state.current[5], 1);
}

TEST(LatchEventModelTests, EpochBarrierDoesNotInventReconvergentXorPulse) {
  Network network{6, {0}, {buffer(0, 1), buffer(0, 2),
      combinational("xor", {1, 2}, {3}, [](const Bits& bits) {
        return Bits{static_cast<uint8_t>(bits[0] ^ bits[1])};
      }), latch("capture", 4, 3, 5)}};
  network.constantByNet.resize(6);
  network.constantByNet[4] = true;
  for (size_t workers : {1, 2, 4}) {
    EventModel model(network, {}, workers);
    auto state = simpleBoundary(model, {0}, {{}, {}, {}, {0}});
    state = model.admit(state, {1});
    while (!model.stable(state)) {
      state = model.successors(state).front();
      EXPECT_EQ(state.current[3], 0);
      EXPECT_EQ(state.current[5], 0);
    }
  }
}

TEST(LatchEventModelTests, NondeterministicProducerChoiceIsSharedByAllFanout) {
  EventModel model({5, {0, 1}, {flipFlop("f", 0, 1, 2), buffer(2, 3), buffer(2, 4)}});
  auto state = simpleBoundary(model, {0, 0}, {{0}, {}, {}});
  const auto results = settle(model, model.admit(state, {1, 1}));
  ASSERT_EQ(results.size(), 2);
  for (const auto& result : results) {
    EXPECT_EQ(result.current[2], result.current[3]);
    EXPECT_EQ(result.current[3], result.current[4]);
  }
}

TEST(LatchEventModelTests, IndependentPrimitiveRacesKeepCartesianProduct) {
  Network network{4, {0, 1}, {flipFlop("f0", 0, 1, 2), flipFlop("f1", 0, 1, 3)}};
  std::vector<State> expected;
  for (size_t workers : {1, 2, 4}) {
    EventModel model(network, {}, workers);
    const auto state = simpleBoundary(model, {0, 0}, {{0}, {0}});
    const auto results = settle(model, model.admit(state, {1, 1}));
    ASSERT_EQ(results.size(), 4);
    if (expected.empty()) expected = results;
    EXPECT_EQ(results, expected);
  }
}

TEST(LatchEventModelTests, DuplicateNetPinsAreDistinctVisitedPositions) {
  auto primitive = latch("duplicated", 0, 0, 1);
  EventModel model({2, {0}, {primitive}});
  auto state = simpleBoundary(model, {1}, {{1}});
  const auto results = settle(model, model.admit(state, {0}));
  // Both positions change: enable first holds 1, data first captures 0.
  ASSERT_EQ(results.size(), 2);
  EXPECT_EQ(results[0].current[1], 0);
  EXPECT_EQ(results[1].current[1], 1);
}

TEST(LatchEventModelTests, EnumeratesAllSixOrdersOfThreeChangedPinPositions) {
  Primitive primitive;
  primitive.name = "order_recorder";
  primitive.inputs = {0, 1, 2};
  primitive.outputs = {3, 4, 5, 6, 7, 8};
  primitive.storageBits = 6;
  primitive.react = [](const Bits& storage, const Bits&, const Bits&,
                        std::optional<size_t> changedPin, bool boot) {
    if (boot || !changedPin) return Reaction{storage, storage};
    size_t code = 0;
    for (size_t bit = 0; bit < storage.size(); ++bit) code |= size_t{storage[bit]} << bit;
    code = 4 * code + *changedPin + 1;
    Bits next(6);
    for (size_t bit = 0; bit < next.size(); ++bit) next[bit] = (code >> bit) & 1;
    return Reaction{next, next};
  };
  EventModel model({9, {0, 1, 2}, {primitive}});
  const auto state = simpleBoundary(model, {0, 0, 0}, {Bits(6)});
  const auto results = settle(model, model.admit(state, {1, 1, 1}));
  EXPECT_EQ(results.size(), 6);
  std::set<size_t> codes;
  for (const auto& result : results) {
    size_t code = 0;
    for (size_t bit = 0; bit < 6; ++bit) code |= size_t{result.storage[0][bit]} << bit;
    codes.insert(code);
  }
  EXPECT_EQ(codes, (std::set<size_t>{27, 30, 39, 45, 54, 57}));
}

TEST(LatchEventModelTests, AsyncClearWorksWithoutClockEdgeIncludingBootstrap) {
  auto primitive = flipFlop("clear_ff", 0, 1, 3);
  const auto baseReaction = primitive.react;
  primitive.inputs.push_back(2);
  primitive.react = [baseReaction](const Bits& storage, const Bits& before, const Bits& pins,
                                  std::optional<size_t> changedPin, bool boot) {
    if (pins[2]) return Reaction{{0}, {0}};
    return baseReaction(storage, before, pins, changedPin, boot);
  };
  EventModel model({4, {0, 1, 2}, {primitive}});
  auto state = simpleBoundary(model, {1, 0, 0}, {{1}});
  state = settle(model, model.admit(state, {1, 0, 1})).front();
  EXPECT_EQ(state.current[3], 0);
  EXPECT_EQ(simpleBoundary(model, {1, 1, 1}, {{1}}).current[3], 0);
}

TEST(LatchEventModelTests, ConflictingControlsProduceStickyNonstableError) {
  auto primitive = latch("controls", 0, 1, 4);
  const auto baseReaction = primitive.react;
  primitive.inputs.insert(primitive.inputs.end(), {2, 3});
  primitive.react = [baseReaction](const Bits& storage, const Bits& before, const Bits& pins,
                                  std::optional<size_t> changedPin, bool boot) {
    if (pins[2] && pins[3]) return Reaction{{}, {}, true, "clear and preset conflict"};
    if (pins[2]) return Reaction{{0}, {0}};
    if (pins[3]) return Reaction{{1}, {1}};
    return baseReaction(storage, before, pins, changedPin, boot);
  };
  EventModel model({5, {0, 1, 2, 3}, {primitive}});
  auto state = model.bootstrap({0, 0, 1, 1}, {{0}}, Bits(5));
  state = model.successors(state).front();
  EXPECT_TRUE(state.error);
  EXPECT_FALSE(model.stable(state));
  EXPECT_NE(state.errorReason.find("conflict"), std::string::npos);
  EXPECT_EQ(model.successors(state), std::vector<State>{state});
}

TEST(LatchEventModelTests, OneFailingPinOrderCannotVanishAmongSuccessfulOrders) {
  Primitive primitive;
  primitive.name = "clear_preset";
  primitive.inputs = {0, 1};
  primitive.outputs = {2};
  primitive.storageBits = 1;
  primitive.react = [](const Bits& storage, const Bits&, const Bits& pins,
                        std::optional<size_t>, bool) -> Reaction {
    if (pins[0] && pins[1]) return {{}, {}, true, "conflicting controls"};
    const Bits value{pins[0] ? uint8_t{0} : pins[1] ? uint8_t{1} : storage[0]};
    return {value, value};
  };
  EventModel model({3, {0, 1}, {primitive}});
  const auto state = simpleBoundary(model, {1, 0}, {{0}});
  const auto outcomes = model.successors(model.admit(state, {0, 1}));
  ASSERT_EQ(outcomes.size(), 2);
  size_t failures = 0;
  size_t successes = 0;
  for (const auto& outcome : outcomes) {
    if (outcome.error) {
      ++failures;
      EXPECT_FALSE(model.stable(outcome));
      EXPECT_EQ(model.successors(outcome), std::vector<State>{outcome});
    } else {
      ++successes;
      EXPECT_EQ(outcome.current[2], 1);
    }
  }
  EXPECT_EQ(failures, 1);
  EXPECT_EQ(successes, 1);
}

TEST(LatchEventModelTests, MultiOutputTupleAndComplementStorageMappingAreAtomic) {
  auto primitive = latch("pair", 0, 1, 2);
  primitive.outputs.push_back(3);
  primitive.initialOutputs = [](const Bits& storage) {
    return Bits{storage[0], static_cast<uint8_t>(!storage[0])};
  };
  const auto baseReaction = primitive.react;
  primitive.react = [baseReaction](const Bits& storage, const Bits& before, const Bits& pins,
                                  std::optional<size_t> changedPin, bool boot) {
    auto result = baseReaction(storage, before, pins, changedPin, boot);
    result.outputs.push_back(static_cast<uint8_t>(!result.outputs[0]));
    return result;
  };
  EventModel model({4, {0, 1}, {primitive}});
  const auto initial = model.bootstrap({0, 0}, {{1}}, Bits(4));
  EXPECT_EQ(initial.current, (Bits{0, 0, 1, 0}));
  auto state = settle(model, initial).front();
  state = settle(model, model.admit(state, {0, 1})).front();
  EXPECT_EQ(state.current, (Bits{0, 1, 0, 1}));
}

TEST(LatchEventModelTests, InputDependentInitialProjectionsReadOneFrozenSnapshot) {
  auto first = latch("first", 1, 0, 2);
  auto second = latch("second", 2, 0, 3);
  first.initialOutputValues = second.initialOutputValues = [](const Bits& storage, const Bits& pins) {
    return Bits{static_cast<uint8_t>(storage[0] & pins[0])};
  };
  // First's projected output becomes 1, but second must read the original seed 0.
  // Swapping primitive order cannot change this simultaneous projection.
  EventModel forward({4, {0, 1}, {first, second}});
  EventModel reversed({4, {0, 1}, {second, first}});
  const auto a = forward.bootstrap({0, 1}, {{1}, {1}}, {0, 0, 0, 0});
  const auto b = reversed.bootstrap({0, 1}, {{1}, {1}}, {0, 0, 0, 0});
  EXPECT_EQ(a.current, (Bits{0, 1, 1, 0}));
  EXPECT_EQ(a.current, b.current);
  EXPECT_EQ(a.previous, a.current);
  EXPECT_EQ(b.previous, b.current);
}

TEST(LatchEventModelTests, AdmissionActivatesExactlyChangedConsumersAndReplacesOldWork) {
  EventModel model({4, {0, 1}, {buffer(0, 2), buffer(1, 3)}});
  auto state = simpleBoundary(model, {0, 0}, {{}, {}});
  state = model.admit(state, {1, 0});
  EXPECT_EQ(state.active, (Bits{1, 0}));
  state = model.successors(state).front();
  EXPECT_EQ(state.active, (Bits{0, 0}));
  EXPECT_EQ(state.previous, state.current);
  EXPECT_TRUE(model.stable(state));
  EXPECT_EQ(model.admit(state, {1, 0}), state);
}

TEST(LatchEventModelTests, UnconsumedExternalChangeStillNormalizesHistory) {
  EventModel model({1, {0}, {}});
  auto state = simpleBoundary(model, {0}, {});
  state = model.admit(state, {1});
  EXPECT_TRUE(model.stable(state));
  EXPECT_EQ(state.previous, state.current);
}

TEST(LatchEventModelTests, InvalidAdmissionDoesNotRemoveTransaction) {
  EventModel model({3, {0, 1}, {latch("l", 0, 1, 2)}});
  const auto initial = model.bootstrap({0, 0}, {{0}}, Bits(3));
  EXPECT_TRUE(model.admit(initial, {1, 0}).error);
  const auto state = settle(model, initial).front();
  EXPECT_TRUE(model.admit(state, {}).error);
  EXPECT_TRUE(model.admit(state, {2, 0}).error);
  const auto bad = model.admit(state, {});
  EXPECT_EQ(model.admit(bad, {0, 0}), bad);
  EXPECT_FALSE(model.stable(bad));
}

TEST(LatchEventModelTests, MissingReactionIsExplicitAbsorbingError) {
  auto primitive = latch("missing", 0, 1, 2);
  primitive.react = {};
  EventModel model({3, {0, 1}, {primitive}});
  const auto next = model.successors(model.bootstrap({0, 0}, {{0}}, Bits(3)));
  ASSERT_EQ(next.size(), 1);
  EXPECT_TRUE(next[0].error);
  EXPECT_FALSE(model.stable(next[0]));
  EXPECT_EQ(model.successors(next[0]), next);
}

TEST(LatchEventModelTests, MalformedAndThrowingReactionRemainVisibleErrors) {
  for (int variant = 0; variant < 4; ++variant) {
    auto primitive = latch("bad", 0, 1, 2);
    primitive.react = [variant](const Bits&, const Bits&, const Bits&, std::optional<size_t>, bool) -> Reaction {
      if (variant == 0) return {{}, {0}};
      if (variant == 1) return {{0}, {}};
      if (variant == 2) return {{0}, {2}};
      throw std::runtime_error("undefined case");
    };
    EventModel model({3, {0, 1}, {primitive}});
    const auto next = model.successors(model.bootstrap({0, 0}, {{0}}, Bits(3)));
    ASSERT_EQ(next.size(), 1);
    EXPECT_TRUE(next[0].error);
    EXPECT_FALSE(model.stable(next[0]));
  }
}

TEST(LatchEventModelTests, MissingOrMalformedInitialOutputMappingIsError) {
  auto primitive = latch("map", 0, 1, 2);
  primitive.outputs.push_back(3);
  EventModel missing({4, {0, 1}, {primitive}});
  EXPECT_TRUE(missing.bootstrap({0, 0}, {{0}}, Bits(4)).error);
  primitive.initialOutputs = [](const Bits&) { return Bits{0, 2}; };
  EventModel malformed({4, {0, 1}, {primitive}});
  EXPECT_TRUE(malformed.bootstrap({0, 0}, {{0}}, Bits(4)).error);
}

TEST(LatchEventModelTests, ConstantsAreFixedAcrossBootstrapAndTransactions) {
  Network network{3, {0}, {latch("constant_enable", 0, 1, 2)}};
  network.constantByNet.resize(3);
  network.constantByNet[1] = true;
  EventModel model(network);
  auto state = simpleBoundary(model, {0}, {{0}});
  state = settle(model, model.admit(state, {1})).front();
  EXPECT_EQ(state.current, (Bits{1, 1, 1}));
  state.current[1] = 0;
  EXPECT_THROW(model.successors(state), std::invalid_argument);
}

TEST(LatchEventModelTests, RejectsMultipleDriversAndUndrivenOrOutOfRangeNets) {
  EXPECT_THROW(EventModel((Network{2, {0}, {buffer(0, 0)}})), std::invalid_argument);
  EXPECT_THROW(EventModel((Network{2, {0}, {buffer(2, 1)}})), std::invalid_argument);
  EXPECT_THROW(EventModel((Network{2, {0}, {buffer(0, 2)}})), std::invalid_argument);
  EXPECT_THROW(EventModel((Network{2, {0}, {}})), std::invalid_argument);
  EXPECT_THROW(EventModel((Network{1, {0, 0}, {}})), std::invalid_argument);
  Network constantConflict{1, {0}, {}, {true}};
  EXPECT_THROW(EventModel{constantConflict}, std::invalid_argument);
  Network badConstants{1, {0}, {}, {true, false}};
  EXPECT_THROW(EventModel{badConstants}, std::invalid_argument);
}

TEST(LatchEventModelTests, RejectsMalformedBootstrapAndStateShapes) {
  EventModel model({3, {0, 1}, {latch("l", 0, 1, 2)}});
  EXPECT_THROW(model.bootstrap({}, {{0}}, Bits(3)), std::invalid_argument);
  EXPECT_THROW(model.bootstrap({2, 0}, {{0}}, Bits(3)), std::invalid_argument);
  EXPECT_THROW(model.bootstrap({0, 0}, {}, Bits(3)), std::invalid_argument);
  EXPECT_THROW(model.bootstrap({0, 0}, {{2}}, Bits(3)), std::invalid_argument);
  EXPECT_THROW(model.bootstrap({0, 0}, {{0}}, Bits(2)), std::invalid_argument);
  auto state = model.bootstrap({0, 0}, {{0}}, Bits(3));
  state.active.push_back(0);
  EXPECT_THROW(model.successors(state), std::invalid_argument);
}

TEST(LatchEventModelTests, PinOrderingResourceLimitThrowsWithoutTruncatingOutcomes) {
  EventModel model({3, {0, 1}, {flipFlop("f", 0, 1, 2)}}, {1, 100});
  auto state = simpleBoundary(model, {0, 0}, {{0}});
  EXPECT_THROW(model.successors(model.admit(state, {1, 1})), Limit);
}

TEST(LatchEventModelTests, SuccessorResourceLimitThrowsWithoutTruncatingOutcomes) {
  EventModel model({4, {0, 1}, {flipFlop("f0", 0, 1, 2), flipFlop("f1", 0, 1, 3)}}, {100, 3});
  auto state = simpleBoundary(model, {0, 0}, {{0}, {0}});
  EXPECT_THROW(model.successors(model.admit(state, {1, 1})), Limit);
}

TEST(LatchEventModelTests, RejectsZeroResourceLimit) {
  EXPECT_THROW(EventModel(Network{}, Limits{0, 1}), std::invalid_argument);
  EXPECT_THROW(EventModel(Network{}, Limits{1, 0}), std::invalid_argument);
}

TEST(LatchEventModelTests, StateKeysIncludeHistoryWorkStorageBootstrapAndErrors) {
  State state{{0}, {0}, {{0}}, {0}, false, false, {}};
  std::set<State> states{state};
  std::set<std::string> keys{state.key()};
  const auto insert = [&](State different) {
    EXPECT_NE(different, state);
    states.insert(different);
    keys.insert(different.key());
  };
  auto different = state; different.current[0] = 1; insert(different);
  different = state; different.previous[0] = 1; insert(different);
  different = state; different.storage[0][0] = 1; insert(different);
  different = state; different.active[0] = 1; insert(different);
  different = state; different.bootstrap = true; insert(different);
  different = state; different.error = true; insert(different);
  different = state; different.errorReason = std::string("a\0b", 3); insert(different);
  EXPECT_EQ(states.size(), 8);
  EXPECT_EQ(keys.size(), 8);
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
