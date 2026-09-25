// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <limits>
#include <memory>
#include <unordered_map>

#include "BoolExprCache.h"
#include "latch/LatchResetAdapter.h"
#include "strategy/SequentialEquivalenceStrategy.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using Values = std::unordered_map<size_t, bool>;

BoolExpr* choose(BoolExpr* condition, BoolExpr* yes, BoolExpr* no) {
  return BoolExpr::Or(BoolExpr::And(condition, yes), BoolExpr::And(BoolExpr::Not(condition), no));
}
SignalKey key(size_t group, size_t index) { return {{group, index}, {0}}; }

// Four remembered PIs plus resettable FF, non-resettable FF and open latch.
// This direct event relation makes the reset wrapper independently testable,
// without relying on Naja extraction or the settling compiler to build its oracle.
SequentialDesignModel synthetic(bool single, bool activeHigh = true, bool falling = false) {
  SequentialDesignModel model;
  model.eventContract = single ? "test-events;single" : "test-events;any";
  auto interface = std::make_shared<EventResetInterface>();
  interface->singleInputChange = single;
  interface->clockInputIndex = 0;
  interface->inputNames = {"clock[0]", "data[0]", "enable[0]", "reset[0]"};
  for (size_t i = 0; i < 4; ++i) {
    const auto inputKey = key(0, i);
    model.environmentInputs.push_back(inputKey);
    model.inputVarByKey.emplace(inputKey, i + 2);
    model.displayNameByKey.emplace(inputKey, interface->inputNames[i]);
    interface->inputKeys.push_back(inputKey);
    interface->currentInputs.push_back(BoolExpr::Var(20 + i));
  }
  if (single) {
    for (size_t i = 0; i < 4; ++i) {
      const auto inputKey = key(1, i);
      model.environmentInputs.push_back(inputKey);
      model.inputVarByKey.emplace(inputKey, 6 + i);
      model.displayNameByKey.emplace(inputKey, "selector" + std::to_string(i));
      if (i < 3) interface->selectorSymbols.push_back(6 + i);
      else interface->valueSymbol = 9;
    }
  }
  std::vector<BoolExpr*> inputs;
  for (size_t i = 0; i < 4; ++i) {
    auto* value = BoolExpr::Var(2 + i);
    if (single) {
      auto* selected = BoolExpr::createTrue();
      for (size_t bit = 0; bit < 3; ++bit) {
        auto* variable = BoolExpr::Var(6 + bit);
        selected = BoolExpr::And(selected, ((i >> bit) & 1) ? variable : BoolExpr::Not(variable));
      }
      value = choose(selected, BoolExpr::Var(9), BoolExpr::Var(20 + i));
    }
    inputs.push_back(value);
  }
  auto* edge = falling ? BoolExpr::And(BoolExpr::Var(20), BoolExpr::Not(inputs[0]))
                      : BoolExpr::And(BoolExpr::Not(BoolExpr::Var(20)), inputs[0]);
  auto* reset = activeHigh ? inputs[3] : BoolExpr::Not(inputs[3]);
  inputs.push_back(choose(edge, choose(reset, BoolExpr::createFalse(), inputs[1]), BoolExpr::Var(24)));
  inputs.push_back(choose(edge, inputs[1], BoolExpr::Var(25)));
  inputs.push_back(choose(inputs[2], inputs[1], BoolExpr::Var(26)));
  for (size_t i = 0; i < inputs.size(); ++i) {
    const auto stateKey = key(2, i);
    model.stateBits.push_back(stateKey);
    model.inputVarByKey.emplace(stateKey, 20 + i);
    model.nextStateExprByStateKey.emplace(stateKey, inputs[i]);
    model.initialStateValueByKey.emplace(stateKey, i >= 4);
    model.displayNameByKey.emplace(stateKey, "state" + std::to_string(i));
  }
  for (size_t i = 0; i < 3; ++i) {
    const auto outputKey = key(3, i);
    model.observedOutputs.push_back(outputKey);
    model.allObservedOutputs.push_back(outputKey);
    model.observedOutputExprByKey.emplace(outputKey, inputs[4 + i]);
    model.displayNameByKey.emplace(outputKey, "out" + std::to_string(i));
  }
  model.eventResetInterface = std::move(interface);
  return model;
}

Values initial(const SequentialDesignModel& model) {
  Values result;
  for (const auto& state : model.stateBits)
    result.emplace(model.inputVarByKey.at(state), model.initialStateValueByKey.at(state));
  return result;
}

Values transaction(bool single, size_t selected, bool value, bool clock = false,
                   bool data = false, bool enable = false, bool reset = false) {
  Values result{{2, clock}, {3, data}, {4, enable}, {5, reset}};
  if (single) {
    for (size_t bit = 0; bit < 3; ++bit) result.emplace(6 + bit, (selected >> bit) & 1);
    result.emplace(9, value);
  }
  return result;
}

std::vector<bool> step(const SequentialDesignModel& model, Values& state, const Values& input) {
  auto environment = state;
  environment.insert(input.begin(), input.end());
  // New reset-order inputs default to code zero unless a test selects an order.
  for (const auto& key : model.environmentInputs)
    environment.try_emplace(model.inputVarByKey.at(key), false);
  std::vector<bool> outputs;
  for (const auto& key : model.observedOutputs)
    outputs.push_back(model.observedOutputExprByKey.at(key)->evaluate(environment));
  for (const auto& key : model.stateBits)
    state[model.inputVarByKey.at(key)] = model.nextStateExprByStateKey.at(key)->evaluate(environment);
  return outputs;
}

size_t remaining(const SequentialDesignModel& model, const Values& state) {
  size_t result = 0;
  for (const auto& key : model.stateBits) {
    const auto& name = model.displayNameByKey.at(key);
    const std::string prefix = "$event.reset.remaining[";
    if (name.find(prefix) == 0 && state.at(model.inputVarByKey.at(key)))
      result |= size_t(1) << std::stoul(name.substr(prefix.size()));
  }
  return result;
}

class LatchResetAdapterTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }
};

TEST_F(LatchResetAdapterTests, PrefixRunsExactCyclesSamplesSymbolicDataAndReleasesReset) {
  for (const bool single : {false, true}) {
    for (const bool activeHigh : {false, true}) {
      for (const bool falling : {false, true}) {
        SCOPED_TRACE(::testing::Message() << "single=" << single << " active=" << activeHigh
            << " falling=" << falling);
        const auto original = synthetic(single, activeHigh, falling);
        const auto adapted = adaptResetCycles(original, {3, {{"reset", activeHigh}}});
        ASSERT_TRUE(adapted.model) << adapted.error;
        const auto& model = *adapted.model;
        EXPECT_EQ(model.eventResetCycles, 3u);
        EXPECT_NE(model.eventContract.find("reset_protocol=low-sample-high-low-release"), std::string::npos);
        EXPECT_EQ(original.stateBits.size(), 7u);
        EXPECT_EQ(original.eventResetCycles, 0u);
        EXPECT_TRUE(original.eventResetInterface);
        auto state = initial(model);
        for (size_t cycle = 0; cycle < 3; ++cycle) {
          const bool data = cycle != 1;
          EXPECT_EQ(remaining(model, state), 3u - cycle);
          EXPECT_EQ(step(model, state, transaction(single, 1, data, true, data, false, !activeHigh)),
                    (std::vector<bool>{false, false, false}));
          EXPECT_FALSE(state.at(20));  // full cycle ends at low carrier
          EXPECT_EQ(state.at(21), data);  // arbitrary sample was retained
          EXPECT_EQ(state.at(23), cycle == 2 ? !activeHigh : activeHigh);
          EXPECT_FALSE(state.at(24));  // synchronous reset observed a real edge
          EXPECT_EQ(state.at(25), data);  // non-reset FF captures the same free sample
        }
        EXPECT_EQ(remaining(model, state), 0u);
        const auto outputs = step(model, state, transaction(single, 4, false, false, true, false, activeHigh));
        EXPECT_EQ(outputs, (std::vector<bool>{false, true, true}));
        EXPECT_EQ(remaining(model, state), 0u);
        EXPECT_EQ(state.at(23), !activeHigh);
      }
    }
  }
}

TEST_F(LatchResetAdapterTests, ResetAndClockSelectorsCannotOverridePrefixAndResetStaysInactiveAfterward) {
  const auto result = adaptResetCycles(synthetic(true), {2, {{"reset", true}}});
  ASSERT_TRUE(result.model) << result.error;
  auto state = initial(*result.model);
  step(*result.model, state, transaction(true, 0, true));  // blocked during sample
  EXPECT_FALSE(state.at(20));
  EXPECT_TRUE(state.at(23));
  step(*result.model, state, transaction(true, 3, false));
  EXPECT_FALSE(state.at(20));
  EXPECT_FALSE(state.at(23));
  step(*result.model, state, transaction(true, 3, true));  // reset permanently inactive
  EXPECT_FALSE(state.at(23));
  step(*result.model, state, transaction(true, 1, true));
  const auto output = step(*result.model, state, transaction(true, 0, true));
  EXPECT_TRUE(output[0]);  // post-prefix real clock events are no longer blocked
  EXPECT_TRUE(output[1]);
  EXPECT_TRUE(state.at(20));
}

TEST_F(LatchResetAdapterTests, LatchTransparencyContinuesDuringPrefixAndHoldAfterRelease) {
  for (const bool single : {false, true}) {
    const auto result = adaptResetCycles(synthetic(single), {3, {{"reset[0]", true}}});
    ASSERT_TRUE(result.model) << result.error;
    auto state = initial(*result.model);
    step(*result.model, state, transaction(single, 2, true, false, false, true));
    EXPECT_FALSE(state.at(26));  // open while D=0
    step(*result.model, state, transaction(single, 1, true, false, true, true));
    EXPECT_TRUE(state.at(26));
    step(*result.model, state, transaction(single, 2, false, false, true, false));
    EXPECT_TRUE(state.at(26));
    const auto output = step(*result.model, state, transaction(single, 1, false, false, false, false));
    EXPECT_TRUE(output[2]);
    EXPECT_TRUE(state.at(26));
  }
}

TEST_F(LatchResetAdapterTests, FirstNormalObservationIsNotMaskedAfterOneCycle) {
  const auto result = adaptResetCycles(synthetic(false), {1, {{"reset", true}}});
  ASSERT_TRUE(result.model) << result.error;
  auto state = initial(*result.model);
  EXPECT_EQ(step(*result.model, state, transaction(false, 0, false, false, false, false)),
            (std::vector<bool>{false, false, false}));
  EXPECT_EQ(step(*result.model, state, transaction(false, 0, false, true, true, true, true)),
            (std::vector<bool>{true, true, true}));
}

TEST_F(LatchResetAdapterTests, ComposedRelationMatchesExplicitCertifiedEventSequenceForEverySample) {
  for (const bool single : {false, true}) {
    for (const bool activeHigh : {false, true}) {
      for (const bool falling : {false, true}) {
        const auto source = synthetic(single, activeHigh, falling);
        const auto result = adaptResetCycles(source, {2, {{"reset", activeHigh}}});
        ASSERT_TRUE(result.model) << result.error;
        // Enumerate all data/enable levels and all reset-order input encodings.
        // Then compare ALL retained source state, not merely observed outputs.
        for (size_t sample = 0; sample < (single ? 16u : 4u); ++sample) {
          auto actual = initial(*result.model);
          auto reference = initial(source);
          for (size_t cycle = 0; cycle < 2; ++cycle) {
            const auto force = [&](size_t pin, bool value) {
              auto input = transaction(single, pin, value, reference.at(20),
                  reference.at(21), reference.at(22), reference.at(23));
              if (!single) input[2 + pin] = value;
              step(source, reference, input);
            };
            force(3, activeHigh);
            force(0, false);
            auto input = transaction(single, 0, false, true, sample & 1, sample & 2, !activeHigh);
            auto sampled = input;
            if (single) {
              for (const auto& key : result.model->environmentInputs) {
                const auto& name = result.model->displayNameByKey.at(key);
                if (name == "$event.reset.order[0][0]") input[result.model->inputVarByKey.at(key)] = sample & 4;
                if (name == "$event.reset.order[1][0]") input[result.model->inputVarByKey.at(key)] = sample & 8;
              }
              const size_t first = sample & 4 ? 2 : 1;
              const size_t second = first == 1 ? 2 : 1;
              force(first, input.at(2 + first));
              force(second, input.at(2 + second));
            } else {
              sampled[2] = false;
              sampled[5] = activeHigh;
              step(source, reference, sampled);
            }
            force(0, true);
            force(0, false);
            if (cycle == 1) force(3, !activeHigh);
            EXPECT_EQ(step(*result.model, actual, input), (std::vector<bool>{false, false, false}));
            for (const auto& key : source.stateBits) {
              const auto id = source.inputVarByKey.at(key);
              EXPECT_EQ(actual.at(id), reference.at(id))
                  << "single=" << single << " polarity=" << activeHigh << " falling=" << falling
                  << " sample=" << sample << " cycle=" << cycle << " variable=" << id;
            }
          }
        }
      }
    }
  }
}

TEST_F(LatchResetAdapterTests, AllNoncontrolLevelsAreFreeBeforeFirstResetClock) {
  auto source = synthetic(true);
  auto* nextClock = source.nextStateExprByStateKey.at(key(2, 0));
  auto* edge = BoolExpr::And(BoolExpr::Not(BoolExpr::Var(20)), nextClock);
  auto* bothInputs = BoolExpr::And(source.nextStateExprByStateKey.at(key(2, 1)),
      source.nextStateExprByStateKey.at(key(2, 2)));
  auto* captured = choose(edge, bothInputs, BoolExpr::Var(25));
  source.nextStateExprByStateKey[key(2, 5)] = captured;
  source.observedOutputExprByKey[key(3, 1)] = captured;
  source.initialStateValueByKey[key(2, 5)] = false;
  const auto result = adaptResetCycles(source, {1, {{"reset", true}}});
  ASSERT_TRUE(result.model) << result.error;
  auto state = initial(*result.model);
  step(*result.model, state, transaction(true, 4, false, false, true, true));
  EXPECT_TRUE(state.at(21));
  EXPECT_TRUE(state.at(22));
  EXPECT_TRUE(state.at(25));  // Both initially-low free pins became high before capture.
}

TEST_F(LatchResetAdapterTests, BothArrivalOrdersRemainSharedEnvironmentChoices) {
  auto source = synthetic(true);
  source.initialStateValueByKey[key(2, 1)] = true;
  source.initialStateValueByKey[key(2, 2)] = true;
  const auto result = adaptResetCycles(source, {1, {{"reset", true}}});
  ASSERT_TRUE(result.model) << result.error;
  for (const bool enableFirst : {false, true}) {
    auto state = initial(*result.model);
    auto input = transaction(true, 4, false, false, false, false);
    for (const auto& key : result.model->environmentInputs)
      if (result.model->displayNameByKey.at(key) == "$event.reset.order[0][0]")
        input[result.model->inputVarByKey.at(key)] = enableFirst;
    step(*result.model, state, input);
    EXPECT_FALSE(state.at(21));
    EXPECT_FALSE(state.at(22));
    EXPECT_EQ(state.at(26), enableFirst);  // E closes first: hold 1; D falls first: capture 0.
  }
}

TEST_F(LatchResetAdapterTests, CounterWidthAndSaturationDoNotAddOrDropCycles) {
  for (const size_t cycles : {size_t(1), size_t(2), size_t(3), size_t(4), size_t(7), size_t(8), size_t(16)}) {
    const auto result = adaptResetCycles(synthetic(true), {cycles, {{"reset", true}}});
    ASSERT_TRUE(result.model) << result.error;
    auto state = initial(*result.model);
    for (size_t i = 0; i < cycles + 2; ++i) {
      EXPECT_EQ(remaining(*result.model, state), i < cycles ? cycles - i : 0u);
      step(*result.model, state, transaction(true, 4, false));
    }
  }
}

TEST_F(LatchResetAdapterTests, RejectsMissingOrAmbiguousContractRatherThanGuessingClock) {
  const SecResetSpec reset{2, {{"reset", true}}};
  auto model = synthetic(true);
  model.eventResetInterface.reset();
  EXPECT_NE(adaptResetCycles(model, reset).error.find("metadata"), std::string::npos);
  model = synthetic(true);
  auto interface = std::make_shared<EventResetInterface>(*model.eventResetInterface);
  model.eventResetInterface = interface;
  interface->clockInputIndex.reset();
  EXPECT_NE(adaptResetCycles(model, reset).error.find("unambiguous"), std::string::npos);
  interface->clockError = "multiple clock carriers";
  EXPECT_EQ(adaptResetCycles(model, reset).error, "multiple clock carriers");
  EXPECT_FALSE(adaptResetCycles(synthetic(true), {0, {{"reset", true}}}).model);
  EXPECT_FALSE(adaptResetCycles(synthetic(true), {std::numeric_limits<size_t>::max(), {{"reset", true}}}).model);
  EXPECT_FALSE(adaptResetCycles(synthetic(true), {1, {}}).model);
  EXPECT_FALSE(adaptResetCycles(synthetic(true), {1, {{"reset", true}, {"enable", true}}}).model);
  EXPECT_FALSE(adaptResetCycles(synthetic(true), {1, {{"missing", true}}}).model);
  EXPECT_FALSE(adaptResetCycles(synthetic(true), {1, {{"clock", true}}}).model);
}

TEST_F(LatchResetAdapterTests, RejectsMalformedRememberedInputsAndSelectorMetadata) {
  const auto check = [](auto mutate) {
    auto model = synthetic(true);
    auto interface = std::make_shared<EventResetInterface>(*model.eventResetInterface);
    model.eventResetInterface = interface;
    mutate(model, *interface);
    const auto result = adaptResetCycles(model, {1, {{"reset", true}}});
    EXPECT_FALSE(result.model);
    EXPECT_FALSE(result.error.empty());
  };
  check([](auto&, auto& i) { i.currentInputs.pop_back(); });
  check([](auto&, auto& i) { i.currentInputs[0] = nullptr; });
  check([](auto&, auto& i) { i.currentInputs[0] = BoolExpr::Var(6); });
  check([](auto&, auto& i) { i.selectorSymbols = {6}; });
  check([](auto&, auto& i) { i.selectorSymbols = {6, 6, 8}; });
  check([](auto&, auto& i) { i.valueSymbol.reset(); });
  check([](auto&, auto& i) { i.inputKeys[1] = i.inputKeys[2]; });
  check([](auto& m, auto&) { m.initialStateValueByKey.clear(); });
  check([](auto& m, auto&) { m.nextStateExprByStateKey.clear(); });
}

TEST_F(LatchResetAdapterTests, BareBusResetNameCannotSilentlySelectBitZero) {
  auto source = synthetic(true);
  auto interface = std::make_shared<EventResetInterface>(*source.eventResetInterface);
  interface->inputNames[2] = "reset[1]";
  source.eventResetInterface = interface;
  const auto ambiguous = adaptResetCycles(source, {1, {{"reset", true}}});
  EXPECT_FALSE(ambiguous.model);
  EXPECT_NE(ambiguous.error.find("bus bit explicitly"), std::string::npos);
  EXPECT_TRUE(adaptResetCycles(source, {1, {{"reset[0]", true}}}).model);
}

TEST_F(LatchResetAdapterTests, AdaptationIsIdempotenceGuardedAndOriginalRemainsUnmodified) {
  const auto source = synthetic(true);
  const auto result = adaptResetCycles(source, {2, {{"reset", true}}});
  ASSERT_TRUE(result.model) << result.error;
  EXPECT_FALSE(result.model->eventResetInterface);
  EXPECT_FALSE(adaptResetCycles(*result.model, {2, {{"reset", true}}}).model);
  EXPECT_EQ(source.eventContract, "test-events;single");
  EXPECT_EQ(source.stateBits.size(), 7u);
  EXPECT_TRUE(source.eventResetInterface);
  EXPECT_TRUE(adaptResetCycles(source, {2, {{"reset", true}}}).model);
}

TEST_F(LatchResetAdapterTests, CompositionBudgetFailsWithoutPublishingPartialModel) {
  for (const size_t limit : {size_t(0), size_t(1), size_t(40)}) {
    auto source = synthetic(true);
    auto interface = std::make_shared<EventResetInterface>(*source.eventResetInterface);
    interface->maxCompositionNodes = limit;
    source.eventResetInterface = interface;
    const auto result = adaptResetCycles(source, {1, {{"reset", true}}});
    EXPECT_FALSE(result.model);
    EXPECT_NE(result.error.find("budget"), std::string::npos);
    EXPECT_EQ(source.stateBits.size(), 7u);
    EXPECT_EQ(source.environmentInputs.size(), 8u);
  }
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
