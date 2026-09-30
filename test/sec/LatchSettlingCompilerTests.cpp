// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "latch/LatchSettlingCompiler.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

State graphState(size_t id, bool stable = false) {
  State state;
  for (size_t bit = 0; bit < 4; ++bit) {
    state.current.push_back(static_cast<uint8_t>((id >> bit) & 1));
  }
  state.previous = state.current;
  state.storage = {{0}};
  state.active = {static_cast<uint8_t>(!stable)};
  return state;
}

EpisodeRelation graphRelation(std::vector<State> states,
                              std::vector<std::vector<size_t>> edges) {
  return {
      [](const State& state) { return state.active == Bits{0} && !state.error; },
      [states = std::move(states), edges = std::move(edges)](const State& state) {
        const auto found = std::find(states.begin(), states.end(), state);
        if (found == states.end()) {
          throw std::invalid_argument("Unknown synthetic reference state");
        }
        std::vector<State> result;
        for (const auto next : edges[static_cast<size_t>(found - states.begin())]) {
          result.push_back(states[next]);
        }
        return result;
      }};
}

Network simpleLatchNetwork() {
  Network network;
  network.netCount = 3;
  network.externalInputs = {0, 1};
  network.primitives = {latch("latch", 0, 1, 2)};
  return network;
}

CompileOptions simpleLatchOptions() {
  CompileOptions options;
  options.initialInputs = {0, 0};
  options.initialStorage = {{uint8_t{0}}};
  options.singleExternalInputChange = true;
  return options;
}

Network selfLatchNetwork(bool inverted) {
  Network network;
  network.netCount = inverted ? 3 : 2;
  network.constantByNet.resize(network.netCount);
  network.constantByNet[1] = true;
  network.primitives.push_back(latch("self", inverted ? 2 : 0, 1, 0));
  if (inverted) {
    network.primitives.push_back(combinational(
        "invert", {0}, {2}, [](const Bits& inputs) { return Bits{uint8_t(!inputs[0])}; }));
  }
  return network;
}

TEST(LatchSettlingCompilerTests, EmptyEntryRelationIsNotAVacuousCertificate) {
  const auto result = certifyEpisode({}, {});
  EXPECT_EQ(result.status, CertificationStatus::Invalid);
  EXPECT_NE(result.detail.find("empty"), std::string::npos);
}

TEST(LatchSettlingCompilerTests, MissingRelationCallbacksAreInvalid) {
  EXPECT_EQ(certifyEpisode({graphState(0)}, {}).status, CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, StableEntryNeedsZeroWavesButIdentitySuccessor) {
  const auto stable = graphState(0, true);
  CompilerLimits limits;
  limits.maxWaves = 0;
  const auto result = certifyEpisode({stable}, graphRelation({stable}, {{0}}), limits);
  ASSERT_TRUE(result.certified()) << result.detail;
  EXPECT_EQ(result.maxWaves, 0);
  EXPECT_EQ(result.stableStates, (std::vector<State>{stable}));
}

TEST(LatchSettlingCompilerTests, MissingNonstableSuccessorCannotDisappear) {
  const auto state = graphState(0);
  EXPECT_EQ(certifyEpisode({state}, graphRelation({state}, {{}})).status,
            CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, MissingStableSuccessorCannotDisappear) {
  const auto state = graphState(0, true);
  EXPECT_EQ(certifyEpisode({state}, graphRelation({state}, {{}})).status,
            CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, StableSuccessorMustPreserveTheEntireState) {
  const auto first = graphState(0, true);
  auto changed = first;
  changed.storage[0][0] = 1;
  EXPECT_EQ(certifyEpisode({first}, graphRelation({first, changed}, {{1}, {1}})).status,
            CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, ExplicitErrorRejectsEvenIfAnotherBranchSettles) {
  const auto entry = graphState(0);
  const auto stable = graphState(1, true);
  auto error = graphState(2);
  error.error = true;
  error.errorReason = "invalid clear/preset combination";
  const auto result = certifyEpisode(
      {entry}, graphRelation({entry, stable, error}, {{1, 2}, {1}, {2}}));
  EXPECT_EQ(result.status, CertificationStatus::Invalid);
  EXPECT_NE(result.detail.find("clear/preset"), std::string::npos);
}

TEST(LatchSettlingCompilerTests, ErrorCannotBeRelabeledStable) {
  auto error = graphState(0, true);
  error.error = true;
  const EpisodeRelation relation{
      [](const State&) { return true; },
      [](const State& state) { return std::vector<State>{state}; }};
  EXPECT_EQ(certifyEpisode({error}, relation).status, CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, LongestPathNotShortestDistanceDefinesTheBound) {
  const auto entry = graphState(0);
  const auto middle = graphState(1);
  const auto stable = graphState(2, true);
  const auto result = certifyEpisode(
      {entry}, graphRelation({entry, middle, stable}, {{2, 1}, {2}, {2}}));
  ASSERT_TRUE(result.certified()) << result.detail;
  EXPECT_EQ(result.maxWaves, 2);
  EXPECT_EQ(result.exploredStates, 3);
  EXPECT_EQ(result.exploredTransitions, 4);
}

TEST(LatchSettlingCompilerTests, AllEntryStatesContributeToTheUniformBound) {
  const auto first = graphState(0);
  const auto second = graphState(1);
  const auto stable = graphState(2, true);
  const auto result = certifyEpisode(
      {stable, second, first}, graphRelation({first, second, stable}, {{1}, {2}, {2}}));
  ASSERT_TRUE(result.certified()) << result.detail;
  EXPECT_EQ(result.maxWaves, 2);
}

TEST(LatchSettlingCompilerTests, NonstableCycleIsNotASufficientlyLargeBound) {
  const auto first = graphState(0);
  const auto second = graphState(1);
  const auto result = certifyEpisode(
      {first}, graphRelation({first, second}, {{1}, {0}}));
  EXPECT_EQ(result.status, CertificationStatus::NonSettling);
}

TEST(LatchSettlingCompilerTests, CycleWithAnExitStillFailsUniversalSettling) {
  const auto entry = graphState(0);
  const auto stable = graphState(1, true);
  const auto result = certifyEpisode(
      {entry}, graphRelation({entry, stable}, {{0, 1}, {1}}));
  EXPECT_EQ(result.status, CertificationStatus::NonSettling);
}

TEST(LatchSettlingCompilerTests, TwoStableEntriesAreNotCherryPicked) {
  const auto first = graphState(0, true);
  const auto second = graphState(1, true);
  const auto result = certifyEpisode(
      {first, second}, graphRelation({first, second}, {{0}, {1}}));
  EXPECT_EQ(result.status, CertificationStatus::OrderDependent);
  EXPECT_EQ(result.stableStates.size(), 2);
}

TEST(LatchSettlingCompilerTests, EqualOutputsDoNotHideDifferentRetainedStorage) {
  const auto entry = graphState(0);
  auto first = graphState(1, true);
  auto second = first;
  second.storage[0][0] = 1;
  const auto result = certifyEpisode(
      {entry}, graphRelation({entry, first, second}, {{1, 2}, {1}, {2}}));
  EXPECT_EQ(result.status, CertificationStatus::OrderDependent);
  ASSERT_EQ(result.stableStates.size(), 2);
  EXPECT_EQ(result.stableStates[0].current, result.stableStates[1].current);
}

TEST(LatchSettlingCompilerTests, HistoryIsPartOfBoundaryUniqueness) {
  const auto entry = graphState(0);
  auto first = graphState(1, true);
  auto second = first;
  second.previous[0] = 0;
  EXPECT_EQ(certifyEpisode({entry}, graphRelation(
                    {entry, first, second}, {{1, 2}, {1}, {2}})).status,
            CertificationStatus::OrderDependent);
}

TEST(LatchSettlingCompilerTests, EquivalentDuplicateSuccessorsDoNotMakeACycle) {
  const auto entry = graphState(0);
  const auto stable = graphState(1, true);
  const auto result = certifyEpisode(
      {entry, entry}, graphRelation({entry, stable}, {{1, 1}, {1, 1}}));
  ASSERT_TRUE(result.certified()) << result.detail;
  EXPECT_EQ(result.maxWaves, 1);
}

TEST(LatchSettlingCompilerTests, ExhaustiveSmallGraphsMatchBoundedUniversalSemantics) {
  const std::vector<State> states = {
      graphState(0), graphState(1), graphState(2, true), graphState(3, true)};
  // Every nonempty successor subset for each of two nonstable states. The two
  // stable states pad identically. Four waves suffice for any acyclic graph on
  // these four states; surviving nonstable paths therefore witness a cycle.
  for (size_t firstMask = 1; firstMask < 16; ++firstMask) {
    for (size_t secondMask = 1; secondMask < 16; ++secondMask) {
      SCOPED_TRACE("masks " + std::to_string(firstMask) + "," +
                   std::to_string(secondMask));
      std::vector<std::vector<size_t>> edges(4);
      for (size_t state = 0; state < 4; ++state) {
        if (firstMask & (size_t{1} << state)) edges[0].push_back(state);
        if (secondMask & (size_t{1} << state)) edges[1].push_back(state);
      }
      edges[2] = {2};
      edges[3] = {3};
      size_t frontier = 1;
      size_t expectedDepth = 0;
      for (size_t depth = 1; depth <= 4; ++depth) {
        size_t next = 0;
        for (size_t state = 0; state < 4; ++state) {
          if (frontier & (size_t{1} << state)) {
            for (const auto successor : edges[state]) next |= size_t{1} << successor;
          }
        }
        frontier = next;
        if ((frontier & 3) == 0 && expectedDepth == 0) expectedDepth = depth;
      }
      const auto expected = (frontier & 3)
                                ? CertificationStatus::NonSettling
                                : frontier == 12 ? CertificationStatus::OrderDependent
                                                 : CertificationStatus::Certified;
      const auto result = certifyEpisode({states[0]}, graphRelation(states, edges));
      EXPECT_EQ(result.status, expected) << result.detail;
      if (expected != CertificationStatus::NonSettling) {
        EXPECT_EQ(result.maxWaves, expectedDepth);
      }
    }
  }
}

TEST(LatchSettlingCompilerTests, StateLayoutCannotChangeDuringPropagation) {
  const auto entry = graphState(0);
  auto malformed = graphState(1, true);
  malformed.current.push_back(0);
  EXPECT_EQ(certifyEpisode({entry}, graphRelation(
                    {entry, malformed}, {{1}, {1}})).status,
            CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, NonBooleanReferenceValuesAreUnsupported) {
  auto malformed = graphState(0, true);
  malformed.storage[0][0] = 2;
  EXPECT_EQ(certifyEpisode({malformed}, graphRelation({malformed}, {{0}})).status,
            CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, InsufficientWaveBudgetIsNotNonSettling) {
  const auto entry = graphState(0);
  const auto stable = graphState(1, true);
  CompilerLimits limits;
  limits.maxWaves = 0;
  const auto result = certifyEpisode(
      {entry}, graphRelation({entry, stable}, {{1}, {1}}), limits);
  EXPECT_EQ(result.status, CertificationStatus::UnprovedBound);
  EXPECT_EQ(result.maxWaves, 1);
}

TEST(LatchSettlingCompilerTests, GraphLimitsNeverTruncateToASuccessfulProof) {
  const auto entry = graphState(0);
  const auto stable = graphState(1, true);
  const auto relation = graphRelation({entry, stable}, {{1}, {1}});
  CompilerLimits limits;
  limits.maxEpisodeStates = 1;
  EXPECT_EQ(certifyEpisode({entry}, relation, limits).status,
            CertificationStatus::ResourceLimit);
  limits = {};
  limits.maxEpisodeTransitions = 1;
  EXPECT_EQ(certifyEpisode({entry}, relation, limits).status,
            CertificationStatus::ResourceLimit);
  limits = {};
  limits.maxStateBits = 1;
  EXPECT_EQ(certifyEpisode({entry}, relation, limits).status,
            CertificationStatus::ResourceLimit);
}

TEST(LatchSettlingCompilerTests, PrimitiveEnumerationLimitIsNotASampledProof) {
  const EpisodeRelation relation{
      [](const State&) { return false; },
      [](const State&) -> std::vector<State> { throw Limit("pin order limit"); }};
  EXPECT_EQ(certifyEpisode({graphState(0)}, relation).status,
            CertificationStatus::ResourceLimit);
}

TEST(LatchSettlingCompilerTests, ReferenceExceptionCannotBecomeAnEmptyRelation) {
  const EpisodeRelation relation{
      [](const State&) { return false; },
      [](const State&) -> std::vector<State> {
        throw std::invalid_argument("missing primitive rule");
      }};
  EXPECT_EQ(certifyEpisode({graphState(0)}, relation).status,
            CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, OpenSelfLatchPreservesBothInitialHistories) {
  EventModel model(selfLatchNetwork(false));
  CompileOptions options;
  const auto result = compileTransitionTable(model, options);
  ASSERT_TRUE(result.certified()) << result.detail;
  ASSERT_TRUE(result.table);
  EXPECT_EQ(result.table->initials.size(), 2);
  EXPECT_EQ(result.table->boundaries.size(), 2);
  EXPECT_EQ(result.table->rows.size(), 2);
  for (const auto& origin : result.table->initials) {
    EXPECT_EQ(origin.storage[0][0], result.table->boundaries[origin.boundary].current[0]);
  }
}

TEST(LatchSettlingCompilerTests, InvertingOpenLatchLoopCannotBeCompiled) {
  EventModel model(selfLatchNetwork(true));
  CompileOptions options;
  options.initialStorage = {{uint8_t{0}}, {}};
  const auto result = compileTransitionTable(model, options);
  EXPECT_EQ(result.status, CertificationStatus::NonSettling) << result.detail;
  EXPECT_FALSE(result.table);
}

TEST(LatchSettlingCompilerTests, AllInitialOriginsRemainWhenTheyMerge) {
  Network network;
  network.netCount = 3;
  network.constantByNet = {true, true, std::nullopt};
  network.primitives = {latch("capture_one", 0, 1, 2)};
  EventModel model(network);
  const auto result = compileTransitionTable(model, {});
  ASSERT_TRUE(result.certified()) << result.detail;
  ASSERT_TRUE(result.table);
  ASSERT_EQ(result.table->initials.size(), 2);
  EXPECT_EQ(result.table->boundaries.size(), 1);
  EXPECT_EQ(result.table->initials[0].boundary, result.table->initials[1].boundary);
  EXPECT_NE(result.table->initials[0].storage, result.table->initials[1].storage);
}

TEST(LatchSettlingCompilerTests, DefaultContractRetainsConcurrentClosingDataRace) {
  EventModel model(simpleLatchNetwork());
  auto options = simpleLatchOptions();
  options.singleExternalInputChange = false;
  const auto result = compileTransitionTable(model, options);
  EXPECT_EQ(result.status, CertificationStatus::OrderDependent) << result.detail;
  EXPECT_FALSE(result.table);
}

TEST(LatchSettlingCompilerTests, ExplicitSingleInputContractHasCompleteClosure) {
  EventModel model(simpleLatchNetwork());
  const auto result = compileTransitionTable(model, simpleLatchOptions());
  ASSERT_TRUE(result.certified()) << result.detail;
  ASSERT_TRUE(result.table);
  EXPECT_EQ(result.table->boundaries.size(), 6);
  EXPECT_EQ(result.table->rows.size(), 18);
  for (size_t from = 0; from < result.table->boundaries.size(); ++from) {
    size_t rowCount = 0;
    size_t stutters = 0;
    const auto& boundary = result.table->boundaries[from];
    for (const auto& row : result.table->rows) {
      if (row.from != from) {
        continue;
      }
      ++rowCount;
      const size_t changes = (row.input[0] != boundary.current[0]) +
                             (row.input[1] != boundary.current[1]);
      EXPECT_LE(changes, 1);
      if (changes == 0) {
        ++stutters;
        EXPECT_EQ(row.to, from);
      }
      const auto replay = certifyEpisode(model, {model.admit(boundary, row.input)});
      ASSERT_TRUE(replay.certified()) << replay.detail;
      EXPECT_EQ(replay.stableStates.front(), result.table->boundaries[row.to]);
    }
    EXPECT_EQ(rowCount, 3);
    EXPECT_EQ(stutters, 1);
  }
}

TEST(LatchSettlingCompilerTests, WorkerScheduleDoesNotChangeCompiledTransitions) {
  auto network = simpleLatchNetwork();
  network.netCount = 4;
  network.primitives.push_back(latch("second_latch", 0, 1, 3));
  auto options = simpleLatchOptions();
  options.initialStorage.push_back({uint8_t{1}});
  const auto serial = compileTransitionTable(EventModel(network, {}, 1), options);
  const auto parallel = compileTransitionTable(EventModel(network, {}, 4), options);
  ASSERT_TRUE(serial.certified()) << serial.detail;
  ASSERT_TRUE(parallel.certified()) << parallel.detail;
  EXPECT_EQ(serial.table->boundaries, parallel.table->boundaries);
  ASSERT_EQ(serial.table->rows.size(), parallel.table->rows.size());
  for (size_t index = 0; index < serial.table->rows.size(); ++index) {
    EXPECT_EQ(serial.table->rows[index].from, parallel.table->rows[index].from);
    EXPECT_EQ(serial.table->rows[index].input, parallel.table->rows[index].input);
    EXPECT_EQ(serial.table->rows[index].to, parallel.table->rows[index].to);
  }
  EXPECT_EQ(serial.table->maxWaves, parallel.table->maxWaves);
}

TEST(LatchSettlingCompilerTests, BootstrapUniversallyChecksAuxiliaryGateSeeds) {
  Network network;
  network.netCount = 4;
  network.externalInputs = {0, 1};
  network.primitives.push_back(combinational(
      "buffer", {0}, {2}, [](const Bits& inputs) { return inputs; }));
  network.primitives.push_back(latch("held", 2, 1, 3));
  const EventModel model(network);
  const auto result = certifyBootstrap(model, {1, 0}, {{}, {0}});
  ASSERT_TRUE(result.certified()) << result.detail;
  ASSERT_EQ(result.stableStates.size(), 1);
  EXPECT_EQ(result.stableStates.front().current, (Bits{1, 0, 1, 0}));
}

TEST(LatchSettlingCompilerTests, GeneratedStartupClockCannotChooseConvenientSeed) {
  Network network;
  network.netCount = 4;
  network.externalInputs = {0, 1};
  network.primitives.push_back(combinational(
      "clock_buffer", {1}, {2}, [](const Bits& inputs) { return inputs; }));
  network.primitives.push_back(flipFlop("capture", 0, 2, 3));
  const EventModel model(network);
  const auto result = certifyBootstrap(model, {1, 1}, {{}, {0}});
  EXPECT_EQ(result.status, CertificationStatus::OrderDependent) << result.detail;
  ASSERT_EQ(result.stableStates.size(), 2);
  EXPECT_NE(result.stableStates[0].storage[1], result.stableStates[1].storage[1]);
}

TEST(LatchSettlingCompilerTests, ProjectionReadsOverwrittenStorageSeedsUniversally) {
  Network network;
  network.netCount = 3;
  network.externalInputs = {2};
  Primitive projected;
  projected.name = "projected_clock";
  projected.inputs = {1};
  projected.outputs = {0};
  projected.storageBits = 1;
  projected.initialOutputValues = [](const Bits&, const Bits& pins) { return pins; };
  projected.react = [](const Bits& storage, const Bits&, const Bits&,
                       std::optional<size_t>, bool) {
    return Reaction{storage, {0}, false, {}};
  };
  network.primitives = {projected, flipFlop("fall_capture", 2, 0, 1, false)};
  // The second output is initialized from stored zero, but its PRE-projection
  // seed is read by the first output's projection. Seed one creates a falling
  // clock event after BOOT; seed zero does not. Enumerating only gate outputs
  // would miss this distinction because this network has no zero-storage gates.
  const auto result = certifyBootstrap(EventModel(network), {1}, {{0}, {0}});
  EXPECT_EQ(result.status, CertificationStatus::OrderDependent) << result.detail;
}

TEST(LatchSettlingCompilerTests, SingleExternalChangeDoesNotExcludeInternalRaces) {
  Network network;
  network.netCount = 4;
  network.externalInputs = {0};
  network.primitives.push_back(combinational(
      "data", {0}, {1}, [](const Bits& inputs) { return inputs; }));
  network.primitives.push_back(combinational(
      "enable", {0}, {2}, [](const Bits& inputs) { return Bits{uint8_t(!inputs[0])}; }));
  network.primitives.push_back(latch("capture", 1, 2, 3));
  // BOOT at input zero ends with an open latch and known data zero, independently
  // of auxiliary gate seeds. A single external rise then makes data rise and enable fall in
  // the same internal wave. The external restriction must not discard either
  // ordering at the latch pins.
  CompileOptions options;
  options.initialInputs = {0};
  options.initialStorage = {{}, {}, {uint8_t{0}}};
  options.singleExternalInputChange = true;
  const auto result = compileTransitionTable(EventModel(network), options);
  EXPECT_EQ(result.status, CertificationStatus::OrderDependent) << result.detail;
  EXPECT_FALSE(result.table);
}

TEST(LatchSettlingCompilerTests, BootstrapGateSeedLimitNeverMeansAllZeroSeeds) {
  Network network;
  network.netCount = 2;
  network.externalInputs = {0};
  network.primitives.push_back(combinational(
      "buffer", {0}, {1}, [](const Bits& inputs) { return inputs; }));
  CompilerLimits limits;
  limits.maxBootstrapSeedBits = 0;
  const auto result = certifyBootstrap(EventModel(network), {0}, {{}}, limits);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
}

TEST(LatchSettlingCompilerTests, ConstantStorageOutputDoesNotAddAuxiliarySeeds) {
  CompilerLimits limits;
  limits.maxBootstrapSeedBits = 0;
  const auto result = certifyBootstrap(EventModel(selfLatchNetwork(false)), {}, {{1}}, limits);
  ASSERT_TRUE(result.certified()) << result.detail;
  EXPECT_EQ(result.stableStates.front().storage, (std::vector<Bits>{{1}}));
}

TEST(LatchSettlingCompilerTests, InvalidInitializationIsNotAnEmptyInitialRelation) {
  const EventModel model(simpleLatchNetwork());
  auto options = simpleLatchOptions();
  options.initialInputs = {0};
  EXPECT_EQ(compileTransitionTable(model, options).status, CertificationStatus::Invalid);
  options = simpleLatchOptions();
  options.initialStorage = {{uint8_t{2}}};
  EXPECT_EQ(compileTransitionTable(model, options).status, CertificationStatus::Invalid);
  options.initialStorage = {{}};
  EXPECT_EQ(compileTransitionTable(model, options).status, CertificationStatus::Invalid);
  options.initialStorage = {{uint8_t{0}}, {}};
  EXPECT_EQ(compileTransitionTable(model, options).status, CertificationStatus::Invalid);
}

TEST(LatchSettlingCompilerTests, InitializationBudgetsDoNotSelectASubset) {
  const EventModel model(selfLatchNetwork(false));
  CompileOptions options;
  options.limits.maxInitialStorageBits = 0;
  EXPECT_EQ(compileTransitionTable(model, options).status, CertificationStatus::ResourceLimit);
  options = {};
  options.limits.maxInitialConfigurations = 1;
  const auto result = compileTransitionTable(model, options);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
  EXPECT_FALSE(result.table);
}

TEST(LatchSettlingCompilerTests, InitializationLimitCoversStorageTimesAuxiliarySeeds) {
  Network network;
  network.netCount = 4;
  network.externalInputs = {0, 1};
  network.primitives.push_back(combinational(
      "buffer", {0}, {2}, [](const Bits& inputs) { return inputs; }));
  network.primitives.push_back(latch("held", 2, 1, 3));
  CompileOptions options;
  options.initialInputs = {0, 0};
  options.singleExternalInputChange = true;
  options.limits.maxInitialConfigurations = 2;
  // Two initial storage choices TIMES two independent auxiliary seed choices.
  const auto result = compileTransitionTable(EventModel(network), options);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
  EXPECT_FALSE(result.table);
}

TEST(LatchSettlingCompilerTests, BoundaryTransactionAndInputLimitsRejectPartialTables) {
  const EventModel model(simpleLatchNetwork());
  auto options = simpleLatchOptions();
  options.limits.maxBoundaryStates = 1;
  auto result = compileTransitionTable(model, options);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
  EXPECT_FALSE(result.table);
  options = simpleLatchOptions();
  options.limits.maxTransactions = 1;
  result = compileTransitionTable(model, options);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
  EXPECT_FALSE(result.table);
  options = simpleLatchOptions();
  options.limits.maxExternalBits = 1;
  result = compileTransitionTable(model, options);
  EXPECT_EQ(result.status, CertificationStatus::ResourceLimit);
  EXPECT_FALSE(result.table);
}

TEST(LatchSettlingCompilerTests, MissingBootstrapRuleCannotBeIgnored) {
  auto network = simpleLatchNetwork();
  network.primitives[0].react = [](const Bits&, const Bits&, const Bits&,
                                 std::optional<size_t>, bool) {
    return Reaction{{}, {}, true, "missing reaction"};
  };
  const auto result = compileTransitionTable(EventModel(network), simpleLatchOptions());
  EXPECT_EQ(result.status, CertificationStatus::Invalid);
  EXPECT_FALSE(result.table);
}

TEST(LatchSettlingCompilerTests, StatusNamesDistinguishUnknownFromNonsettling) {
  EXPECT_STREQ(certificationStatusName(CertificationStatus::Certified), "certified");
  EXPECT_STREQ(certificationStatusName(CertificationStatus::NonSettling), "non-settling");
  EXPECT_STREQ(certificationStatusName(CertificationStatus::UnprovedBound), "unproved-bound");
  EXPECT_STREQ(certificationStatusName(CertificationStatus::ResourceLimit), "resource-limit");
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
