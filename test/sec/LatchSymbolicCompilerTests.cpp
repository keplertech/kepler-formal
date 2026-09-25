// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "BoolExprCache.h"
#include "latch/LatchSymbolicCompiler.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

class LatchSymbolicCompilerTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }

  SymbolicNetwork lifted(Network reference) {
    SymbolicNetwork result;
    for (const auto& primitive : reference.primitives) {
      result.primitives.push_back(liftPrimitive(primitive));
    }
    result.reference = std::move(reference);
    return result;
  }

  SymbolicCompileOptions options(const Network& network, bool single = true) {
    SymbolicCompileOptions result;
    result.initialInputs.assign(network.externalInputs.size(), 0);
    for (const auto& primitive : network.primitives) {
      result.initialStorage.emplace_back(primitive.storageBits, 0);
    }
    result.singleExternalInputChange = single;
    result.maxWaves = 32;
    return result;
  }

  Network latchNetwork() {
    return Network{3, {0, 1}, {latch("latch", 0, 1, 2)}};
  }

  Bits flatten(const State& state) {
    Bits result = state.current;
    for (const auto& storage : state.storage) {
      result.insert(result.end(), storage.begin(), storage.end());
    }
    return result;
  }

  Bits evaluate(const SymbolicBits& expressions, const SymbolicMacro& macro,
                const Bits& state, const Bits& inputs) {
    std::unordered_map<size_t, bool> values;
    for (size_t i = 0; i < state.size(); ++i) values.emplace(macro.stateSymbols.at(i), state[i]);
    for (size_t i = 0; i < inputs.size(); ++i) values.emplace(macro.inputSymbols.at(i), inputs[i]);
    Bits result;
    for (auto* expression : expressions) result.push_back(expression->evaluate(values));
    return result;
  }

  Bits initialValues(const SymbolicMacro& macro, const Bits& parameters) {
    std::unordered_map<size_t, bool> values;
    for (size_t i = 0; i < parameters.size(); ++i)
      values.emplace(macro.initialParameterSymbols.at(i), parameters[i]);
    Bits result;
    for (auto* expression : macro.initialStateExpressions)
      result.push_back(expression->evaluate(values));
    return result;
  }

  void compareWithFinite(const Network& network, bool single = true) {
    const auto settings = options(network, single);
    const auto symbolic = compileSymbolicNetwork(lifted(network), settings);
    ASSERT_TRUE(symbolic.certified()) << symbolic.detail;
    CompileOptions exactSettings;
    exactSettings.initialInputs = settings.initialInputs;
    exactSettings.singleExternalInputChange = single;
    for (const auto& storage : settings.initialStorage) {
      exactSettings.initialStorage.emplace_back(storage.begin(), storage.end());
    }
    const auto exact = compileTransitionTable(EventModel(network), exactSettings);
    ASSERT_TRUE(exact.certified()) << exact.detail;
    const auto& macro = *symbolic.model;
    EXPECT_EQ(macro.singleExternalInputChange, single);
    EXPECT_EQ(macro.externalInputNets, network.externalInputs);
    EXPECT_EQ(macro.initialState, flatten(exact.table->boundaries.at(exact.table->initials[0].boundary)));
    for (const auto& row : exact.table->rows) {
      const auto old = flatten(exact.table->boundaries[row.from]);
      const auto expected = flatten(exact.table->boundaries[row.to]);
      EXPECT_EQ(evaluate(macro.nextState, macro, old, row.input), expected);
      EXPECT_EQ(evaluate(macro.observedNets, macro, old, row.input),
                exact.table->boundaries[row.to].current);
    }
  }
};

TEST_F(LatchSymbolicCompilerTests, SingleInputLatchMatchesEveryExactTableRow) {
  compareWithFinite(latchNetwork());
}

TEST_F(LatchSymbolicCompilerTests, IndependentFlopAndLatchMatchExactHistory) {
  auto network = latchNetwork();
  network.netCount = 5;
  network.externalInputs.push_back(3);
  network.primitives.push_back(flipFlop("ff", 2, 3, 4));
  compareWithFinite(network);
}

TEST_F(LatchSymbolicCompilerTests, AllChangeCombinationalModelMatchesExactTable) {
  Network network{3, {0, 1}, {combinational("xor", {0, 1}, {2}, [](const Bits& pins) {
    return Bits{uint8_t(pins[0] ^ pins[1])};
  })}};
  compareWithFinite(network, false);
}

TEST_F(LatchSymbolicCompilerTests, StoredSelfFeedbackIsNotClassifiedAsOscillation) {
  Network network{2, {}, {latch("self", 0, 1, 0)}};
  network.constantByNet = {std::nullopt, true};
  compareWithFinite(network, false);
}

TEST_F(LatchSymbolicCompilerTests, CapturesNegativeClockEdgesExactly) {
  Network network{3, {0, 1}, {flipFlop("falling", 0, 1, 2, false)}};
  compareWithFinite(network);
}

TEST_F(LatchSymbolicCompilerTests, WideStateAndInputSpacesNeedNoValuationEnumeration) {
  Network network;
  constexpr size_t width = 16;
  network.netCount = width * 3;
  for (size_t input = 0; input < width * 2; ++input) network.externalInputs.push_back(input);
  for (size_t bit = 0; bit < width; ++bit) {
    network.primitives.push_back(latch("l" + std::to_string(bit),
                                      bit * 2, bit * 2 + 1, width * 2 + bit));
  }
  const auto result = compileSymbolicNetwork(lifted(network), options(network));
  ASSERT_TRUE(result.certified()) << result.detail;
  EXPECT_EQ(result.model->inputSymbols.size(), 32);
  EXPECT_EQ(result.model->stateSymbols.size(), 64);
  Bits state = result.model->initialState;
  Bits input(network.externalInputs.size(), 0);
  input[0] = 1;
  state = evaluate(result.model->nextState, *result.model, state, input);
  EXPECT_EQ(state[32], 0);
  input[1] = 1;
  state = evaluate(result.model->nextState, *result.model, state, input);
  EXPECT_EQ(state[32], 1);
  EXPECT_EQ(state[48], 1);
  EXPECT_EQ(state[33], 0);
}

TEST_F(LatchSymbolicCompilerTests, GrowingBoundHandlesAChainBeyondTwelveStateBits) {
  Network network;
  constexpr size_t length = 13;
  network.netCount = length + 2;
  network.externalInputs = {0, 1};
  for (size_t stage = 0; stage < length; ++stage) {
    network.primitives.push_back(latch("stage" + std::to_string(stage),
        stage ? stage + 1 : 0, 1, stage + 2));
  }
  auto settings = options(network);
  settings.maxWaves = 16;
  const auto result = compileSymbolicNetwork(lifted(network), settings);
  ASSERT_TRUE(result.certified()) << result.detail;
  EXPECT_GE(result.model->transitionWaves, length);
  EXPECT_LE(result.model->transitionWaves, settings.maxWaves);
  auto state = result.model->initialState;
  state = evaluate(result.model->nextState, *result.model, state, {1, 0});
  state = evaluate(result.model->nextState, *result.model, state, {1, 1});
  EXPECT_EQ(state[length + 1], 1);
}

TEST_F(LatchSymbolicCompilerTests, SimultaneousClosingDataRaceIsNotAProvenTransition) {
  const auto network = latchNetwork();
  const auto result = compileSymbolicNetwork(lifted(network), options(network, false));
  EXPECT_EQ(result.status, CertificationStatus::UnprovedBound) << result.detail;
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, NonsettlingBootstrapIsUnprovedNotSilentlyExcluded) {
  Network network{3, {}, {latch("self", 2, 1, 0),
      combinational("invert", {0}, {2}, [](const Bits& pins) { return Bits{uint8_t(!pins[0])}; })}};
  network.constantByNet = {std::nullopt, true, std::nullopt};
  auto settings = options(network);
  settings.maxWaves = 8;
  const auto result = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(result.status, CertificationStatus::UnprovedBound) << result.detail;
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, ATooSmallCandidateDoesNotClaimConvergence) {
  Network network{5, {0, 1}, {latch("first", 0, 1, 2),
      latch("second", 2, 1, 3), latch("third", 3, 1, 4)}};
  auto settings = options(network);
  settings.maxWaves = 1;
  const auto result = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(result.status, CertificationStatus::UnprovedBound) << result.detail;
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, BootstrapErrorsAreNotPrunedOrAbsorbedAsSuccess) {
  auto network = latchNetwork();
  auto symbolic = lifted(network);
  // Reference callback error must remain represented through every wave.
  symbolic.primitives[0].react = [](const SymbolicBits& storage,
      const SymbolicBits&, const SymbolicBits&, std::optional<size_t>, bool) {
    return SymbolicReaction{storage, storage, BoolExpr::createTrue()};
  };
  auto settings = options(network);
  settings.maxWaves = 2;
  const auto result = compileSymbolicNetwork(symbolic, settings);
  EXPECT_EQ(result.status, CertificationStatus::UnprovedBound) << result.detail;
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, CandidateInvariantFailureIsNotAReachableBugReport) {
  Primitive cell;
  cell.name = "unreachable_error";
  cell.inputs = {0};
  cell.outputs = {1};
  cell.storageBits = 1;
  cell.react = [](const Bits& storage, const Bits&, const Bits& pins,
                  std::optional<size_t>, bool) {
    return Reaction{storage, storage, bool(storage[0] && pins[0])};
  };
  Network network{2, {0}, {cell}};
  const auto settings = options(network);
  const auto symbolic = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(symbolic.status, CertificationStatus::UnprovedBound) << symbolic.detail;
  EXPECT_FALSE(symbolic.model);
  // The error is not reachable from intended storage=0. The generic invariant
  // nevertheless includes storage=1/input=0, so rejecting the symbolic proof
  // cannot be reported as a reachable circuit bug. The exact fallback succeeds.
  CompileOptions exactSettings;
  exactSettings.initialInputs = {0};
  exactSettings.initialStorage = {{false}};
  exactSettings.singleExternalInputChange = true;
  EXPECT_TRUE(compileTransitionTable(EventModel(network), exactSettings).certified());
}

TEST_F(LatchSymbolicCompilerTests, UniquenessIncludesHiddenStorageNotOnlyOutputs) {
  auto cell = latch("hidden", 0, 1, 2);
  const auto reaction = cell.react;
  cell.react = [reaction](const Bits& storage, const Bits& before, const Bits& after,
                          std::optional<size_t> changedPin, bool bootstrap) {
    auto result = reaction(storage, before, after, changedPin, bootstrap);
    result.outputs = {0};
    return result;
  };
  cell.initialOutputs = [](const Bits&) { return Bits{0}; };
  Network network{3, {0, 1}, {cell}};
  const auto result = compileSymbolicNetwork(lifted(network), options(network, false));
  EXPECT_EQ(result.status, CertificationStatus::UnprovedBound) << result.detail;
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, NonPowerOfTwoLimitIsCheckedWithoutOvershooting) {
  Network network{5, {0, 1}, {latch("first", 0, 1, 2),
      latch("second", 2, 1, 3), latch("third", 3, 1, 4)}};
  auto settings = options(network);
  settings.maxWaves = 3;
  const auto result = compileSymbolicNetwork(lifted(network), settings);
  ASSERT_TRUE(result.certified()) << result.detail;
  EXPECT_EQ(result.model->transitionWaves, 3);
}

TEST_F(LatchSymbolicCompilerTests, ZeroWaveLimitDoesNotSkipBootstrap) {
  const auto network = latchNetwork();
  auto settings = options(network);
  settings.maxWaves = 0;
  const auto result = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(result.status, CertificationStatus::UnprovedBound) << result.detail;
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, AuxiliaryBootstrapSeedsRemainIndependent) {
  Network network{4, {0, 1}, {combinational("clock", {1}, {2},
      [](const Bits& pins) { return pins; }), flipFlop("capture", 0, 2, 3)}};
  auto settings = options(network);
  settings.initialInputs = {1, 1};
  const auto result = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(result.status, CertificationStatus::OrderDependent) << result.detail;
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, InsufficientNodeAndSatBudgetsNeverReturnPartialModels) {
  const auto network = latchNetwork();
  auto settings = options(network);
  settings.maxNodes = 2;
  auto result = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
  EXPECT_FALSE(result.model);
  settings = options(network);
  settings.maxSatConflicts = 0;
  result = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
  EXPECT_FALSE(result.model);
  settings = options(network);
  settings.maxSatDecisions = 0;
  result = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, InvalidInitializationIsRejected) {
  const auto network = latchNetwork();
  auto settings = options(network);
  settings.initialInputs = {0};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
  settings = options(network);
  settings.initialStorage = {{2}};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
  settings = options(network);
  settings.initialStorage = {{}};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
}

TEST_F(LatchSymbolicCompilerTests, UnspecifiedInputsAndStorageRetainEveryBooleanBootstrapOrigin) {
  const auto network = latchNetwork();
  SymbolicCompileOptions settings;
  settings.singleExternalInputChange = true;
  const auto compiled = compileSymbolicNetwork(lifted(network), settings);
  ASSERT_TRUE(compiled.certified()) << compiled.detail;
  const auto& macro = *compiled.model;
  EXPECT_TRUE(macro.initialState.empty());
  ASSERT_EQ(macro.initialParameterSymbols.size(), 3u);
  ASSERT_EQ(macro.initialInputSymbols.size(), 2u);
  ASSERT_EQ(macro.initialStorageSymbols.size(), 1u);
  EXPECT_EQ(macro.initialInputSymbols[0], macro.initialParameterSymbols[0]);
  EXPECT_EQ(macro.initialInputSymbols[1], macro.initialParameterSymbols[1]);
  EXPECT_EQ(macro.initialStorageSymbols[0][0], macro.initialParameterSymbols[2]);
  for (size_t code = 0; code < 8; ++code) {
    const Bits input{uint8_t(code & 1), uint8_t((code >> 1) & 1)};
    const Bits stored{uint8_t((code >> 2) & 1)};
    const auto expected = certifyBootstrap(EventModel(network), input, {stored});
    ASSERT_TRUE(expected.certified()) << expected.detail;
    const auto state = initialValues(macro, {input[0], input[1], stored[0]});
    EXPECT_EQ(state, flatten(expected.stableStates.front()));
    EXPECT_EQ(evaluate(macro.nextState, macro, state, input), state);
    EXPECT_EQ(state[2], input[1] ? input[0] : stored[0]);
  }
}

TEST_F(LatchSymbolicCompilerTests, MixedPerBitRestrictionsDoNotFixUnspecifiedOrigins) {
  SymbolicCompileOptions settings;
  settings.singleExternalInputChange = true;
  settings.initialInputValues = {std::nullopt, uint8_t(0)};
  settings.initialStorageValues = {{uint8_t(1)}};
  const auto compiled = compileSymbolicNetwork(lifted(latchNetwork()), settings);
  ASSERT_TRUE(compiled.certified()) << compiled.detail;
  const auto& macro = *compiled.model;
  ASSERT_EQ(macro.initialParameterSymbols.size(), 1u);
  EXPECT_TRUE(macro.initialInputSymbols[0]);
  EXPECT_FALSE(macro.initialInputSymbols[1]);
  EXPECT_FALSE(macro.initialStorageSymbols[0][0]);
  EXPECT_EQ(initialValues(macro, {0}), (Bits{0, 0, 1, 1}));
  EXPECT_EQ(initialValues(macro, {1}), (Bits{1, 0, 1, 1}));
}

TEST_F(LatchSymbolicCompilerTests, UnknownSelfHeldStorageIsHistoryNotSchedulerNondeterminism) {
  Network network{2, {}, {latch("self", 0, 1, 0)}};
  network.constantByNet = {std::nullopt, true};
  const auto compiled = compileSymbolicNetwork(lifted(network), {});
  ASSERT_TRUE(compiled.certified()) << compiled.detail;
  const auto& macro = *compiled.model;
  ASSERT_EQ(macro.initialParameterSymbols.size(), 1u);
  EXPECT_TRUE(macro.initialInputSymbols.empty());
  for (bool value : {false, true}) {
    const auto state = initialValues(macro, {uint8_t(value)});
    EXPECT_EQ(state, (Bits{uint8_t(value), 1, uint8_t(value)}));
    EXPECT_EQ(evaluate(macro.nextState, macro, state, {}), state);
  }
}

TEST_F(LatchSymbolicCompilerTests, UnknownInitialClockDoesNotFabricateAnEdge) {
  Network network{3, {0, 1}, {flipFlop("ff", 0, 1, 2)}};
  SymbolicCompileOptions settings;
  settings.singleExternalInputChange = true;
  const auto compiled = compileSymbolicNetwork(lifted(network), settings);
  ASSERT_TRUE(compiled.certified()) << compiled.detail;
  const auto& macro = *compiled.model;
  for (size_t code = 0; code < 8; ++code) {
    const uint8_t data = code & 1, clock = (code >> 1) & 1, stored = (code >> 2) & 1;
    const auto state = initialValues(macro, {data, clock, stored});
    EXPECT_EQ(state, (Bits{data, clock, stored, stored}));
    EXPECT_EQ(evaluate(macro.nextState, macro, state, {data, clock}), state);
    const auto changed = evaluate(macro.nextState, macro, state, {data, uint8_t(!clock)});
    EXPECT_EQ(changed[2], clock ? stored : data);
  }
}

TEST_F(LatchSymbolicCompilerTests, BootstrapResetCanEliminateUnknownStorageWithoutChoosingItsValue) {
  Primitive cell;
  cell.name = "async_clear";
  cell.inputs = {0}; cell.outputs = {1}; cell.storageBits = 1;
  cell.react = [](const Bits& storage, const Bits&, const Bits& pins,
                  std::optional<size_t>, bool) {
    const Bits next{uint8_t(pins[0] ? 0 : storage[0])};
    return Reaction{next, next};
  };
  Network network{2, {0}, {cell}};
  SymbolicCompileOptions settings;
  settings.initialInputValues = {uint8_t(1)};
  const auto compiled = compileSymbolicNetwork(lifted(network), settings);
  ASSERT_TRUE(compiled.certified()) << compiled.detail;
  const auto& macro = *compiled.model;
  ASSERT_EQ(macro.initialParameterSymbols.size(), 1u);
  EXPECT_EQ(macro.initialState, (Bits{1, 0, 0}));
  EXPECT_EQ(initialValues(macro, {0}), macro.initialState);
  EXPECT_EQ(initialValues(macro, {1}), macro.initialState);
}

TEST_F(LatchSymbolicCompilerTests, AnInvalidUnknownBootstrapOriginCannotDisappear) {
  auto cell = latch("reject_one", 0, 1, 2);
  const auto react = cell.react;
  cell.react = [react](const Bits& storage, const Bits& before, const Bits& pins,
                      std::optional<size_t> pin, bool boot) {
    auto result = react(storage, before, pins, pin, boot);
    result.error = boot && storage[0];
    return result;
  };
  Network network{3, {0, 1}, {cell}};
  auto settings = options(network);
  settings.initialStorage.clear();
  settings.maxWaves = 4;
  const auto compiled = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(compiled.status, CertificationStatus::UnprovedBound) << compiled.detail;
  EXPECT_FALSE(compiled.model);
}

TEST_F(LatchSymbolicCompilerTests, OriginRestrictionShapesAndConflictingFormsAreRejected) {
  const auto network = latchNetwork();
  auto settings = options(network);
  settings.initialInputValues = {std::nullopt, std::nullopt};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
  settings = options(network);
  settings.initialStorageValues = {{std::nullopt}};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
  settings = {};
  settings.initialInputValues = {std::nullopt};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
  settings.initialInputValues = {uint8_t(2), std::nullopt};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
  settings = {};
  settings.initialStorageValues = {{}};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
  settings.initialStorageValues = {{uint8_t(2)}};
  EXPECT_EQ(compileSymbolicNetwork(lifted(network), settings).status, CertificationStatus::Invalid);
}

TEST_F(LatchSymbolicCompilerTests, NonzeroSatBudgetExhaustionDoesNotBecomeAProof) {
  Network network{5, {0, 1}, {latch("first", 0, 1, 2),
      latch("second", 2, 1, 3), latch("third", 3, 1, 4)}};
  auto settings = options(network);
  settings.maxSatConflicts = 1;
  settings.maxSatDecisions = 1;
  const auto result = compileSymbolicNetwork(lifted(network), settings);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit) << result.detail;
  EXPECT_FALSE(result.model);
}

TEST_F(LatchSymbolicCompilerTests, ParallelWorkersPreserveTheCompiledFunction) {
  Network network{4, {0, 1}, {latch("a", 0, 1, 2), latch("b", 0, 1, 3)}};
  auto settings = options(network);
  settings.workers = 1;
  const auto serial = compileSymbolicNetwork(lifted(network), settings);
  settings.workers = 4;
  const auto parallel = compileSymbolicNetwork(lifted(network), settings);
  ASSERT_TRUE(serial.certified()) << serial.detail;
  ASSERT_TRUE(parallel.certified()) << parallel.detail;
  EXPECT_EQ(serial.model->initialState, parallel.model->initialState);
  EXPECT_EQ(serial.model->stateSymbols, parallel.model->stateSymbols);
  EXPECT_EQ(serial.model->inputSymbols, parallel.model->inputSymbols);
  auto state = serial.model->initialState;
  for (const auto& input : std::vector<Bits>{{1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}}) {
    const auto a = evaluate(serial.model->nextState, *serial.model, state, input);
    const auto b = evaluate(parallel.model->nextState, *parallel.model, state, input);
    EXPECT_EQ(a, b);
    state = a;
  }
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
