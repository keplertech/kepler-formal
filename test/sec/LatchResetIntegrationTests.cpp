// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <tuple>
#include <unordered_map>
#include "BoolExprCache.h"
#include "DNL.h"
#include "NLDB.h"
#include "NLDB0.h"
#include "NLLibrary.h"
#include "NLName.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "latch/LatchResetAdapter.h"
#include "latch/LatchSupportOptions.h"
#include "model/SequentialDesignModel.h"
#include "strategy/SequentialEquivalenceStrategy.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using namespace naja::NL;
using Modeling = SNLDesignModeling;
using Expression = Modeling::BooleanExpression;
using Operator = Expression::Operator;

Expression termExpression(SNLBitTerm* term) {
  Expression expression;
  expression.root = expression.addTerm(term);
  return expression;
}

enum class ClockRoute { Direct, Buffered, Inverted, Constant, Gated, Multiple };
struct CircuitOptions {
  ClockRoute clock = ClockRoute::Direct;
  bool latchChain = true;
  bool mixedEdges = false;
  bool asyncThroughLatch = false;
  bool invertData = false;
  bool invertOnlyDuringReset = false;
  bool ignoreReset = false;
  bool andDataEnable = false;
  bool constantData = false;
};

class LatchResetIntegrationTests : public ::testing::Test {
 protected:
  void SetUp() override {
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    designs_ = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
    primitives_ = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("primitives"));
  }
  void TearDown() override {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) universe->destroy();
    BoolExprCache::destroy();
  }
  SupportOptions settings(bool initialStorage = false) {
    SupportOptions options;
    options.enabled = true;
    options.singleInputChange = true;
    options.initialInputs = false;
    options.initialStorage = initialStorage;
    options.workers = 2;
    return options;
  }
  SNLScalarNet* port(SNLDesign* top, const std::string& name, SNLTerm::Direction direction) {
    auto* term = SNLScalarTerm::create(top, direction, NLName(name));
    auto* net = SNLScalarNet::create(top, NLName(name));
    term->setNet(net);
    return net;
  }
  void instance(SNLDesign* top, SNLDesign* primitive, const std::string& name,
                const std::vector<std::pair<std::string, SNLScalarNet*>>& pins) {
    auto* cell = SNLInstance::create(top, primitive, NLName(name));
    for (const auto& [pin, net] : pins)
      cell->getInstTerm(primitive->getScalarTerm(NLName(pin)))->setNet(net);
  }
  SNLDesign* flop(const std::string& name, const CircuitOptions& options, bool falling = false) {
    auto* cell = SNLDesign::create(primitives_, SNLDesign::Type::Primitive, NLName(name));
    auto* data = options.constantData ? nullptr
        : SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("D"));
    auto* clock = SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("C"));
    auto* reset = SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("R"));
    auto* output = SNLScalarTerm::create(cell, SNLTerm::Direction::Output, NLName("Q"));
    Modeling::SequentialModel model;
    model.kind = Modeling::SequentialModel::Kind::FlipFlop;
    model.clockedOn = termExpression(clock);
    if (falling) model.clockedOn.root = model.clockedOn.addOperation(Operator::Not, {model.clockedOn.root});
    Modeling::SequentialState state;
    if (options.constantData) state.nextState.root = state.nextState.addConstant(false);
    else state.nextState = termExpression(data);
    if (options.invertData)
      state.nextState.root = state.nextState.addOperation(Operator::Not, {state.nextState.root});
    if (options.asyncThroughLatch) state.clear = termExpression(reset);
    else if (!options.ignoreReset) {
      auto inactiveReset = state.nextState.addTerm(reset);
      inactiveReset = state.nextState.addOperation(Operator::Not, {inactiveReset});
      state.nextState.root = state.nextState.addOperation(Operator::And, {state.nextState.root, inactiveReset});
    }
    model.states.push_back(state);
    Expression physical;
    physical.root = physical.addState(0);
    if (options.invertOnlyDuringReset) {
      const auto asserted = physical.addTerm(reset);
      physical.root = physical.addOperation(Operator::Xor, {physical.root, asserted});
    }
    model.outputs.push_back({output, physical});
    Modeling::setSequentialModel(cell, model);
    return cell;
  }
  SNLDesign* routeCell(const std::string& name, ClockRoute route) {
    auto* cell = SNLDesign::create(primitives_, SNLDesign::Type::Primitive, NLName(name));
    auto* a = SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("A"));
    Modeling::BitTerms inputs{a};
    if (route == ClockRoute::Gated)
      inputs.push_back(SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("B")));
    auto* y = SNLScalarTerm::create(cell, SNLTerm::Direction::Output, NLName("Y"));
    Modeling::addCombinatorialArcs(inputs, {y});
    const auto arity = inputs.size();
    Modeling::setTruthTable(cell, SNLTruthTable(arity,
        route == ClockRoute::Gated ? 8 : route == ClockRoute::Inverted ? 1 : 2,
        SNLTruthTable::fullDependencies(arity)));
    return cell;
  }
  SNLDesign* circuit(const std::string& name, CircuitOptions options = {}) {
    auto* top = SNLDesign::create(designs_, SNLDesign::Type::Standard, NLName(name));
    auto* clock = port(top, "clock", SNLTerm::Direction::Input);
    auto* data = port(top, "data", SNLTerm::Direction::Input);
    auto* enable = port(top, "enable", SNLTerm::Direction::Input);
    auto* reset = port(top, "reset", SNLTerm::Direction::Input);
    auto* output = port(top, "out", SNLTerm::Direction::Output);
    auto* localClock = clock;
    if (options.clock == ClockRoute::Constant) {
      localClock = SNLScalarNet::create(top, NLName("constant_clock"));
      localClock->setType(SNLNet::Type::Assign0);
    } else if (options.clock != ClockRoute::Direct && options.clock != ClockRoute::Multiple) {
      localClock = SNLScalarNet::create(top, NLName("routed_clock"));
      auto* primitive = routeCell(name + "_route", options.clock);
      std::vector<std::pair<std::string, SNLScalarNet*>> pins{{"A", clock}, {"Y", localClock}};
      if (options.clock == ClockRoute::Gated) pins.emplace_back("B", enable);
      instance(top, primitive, "clock_route", pins);
    }
    auto* localReset = reset;
    if (options.asyncThroughLatch) {
      localReset = SNLScalarNet::create(top, NLName("latched_reset"));
      instance(top, NLDB0::getDLatch(), "reset_latch", {{"D", reset}, {"E", enable}, {"Q", localReset}});
    }
    auto* held = options.latchChain || options.mixedEdges
        ? SNLScalarNet::create(top, NLName("held")) : output;
    auto* localData = data;
    if (options.andDataEnable) {
      localData = SNLScalarNet::create(top, NLName("qualified_data"));
      instance(top, routeCell(name + "_data_gate", ClockRoute::Gated), "data_gate",
          {{"A", data}, {"B", enable}, {"Y", localData}});
    }
    std::vector<std::pair<std::string, SNLScalarNet*>> firstPins{
        {"C", localClock}, {"R", localReset}, {"Q", held}};
    if (!options.constantData) firstPins.emplace_back("D", localData);
    instance(top, flop(name + "_ff", options), "ff", firstPins);
    if (options.mixedEdges) {
      auto* falling = options.latchChain ? SNLScalarNet::create(top, NLName("falling")) : output;
      std::vector<std::pair<std::string, SNLScalarNet*>> fallingPins{
          {"C", localClock}, {"R", localReset}, {"Q", falling}};
      if (!options.constantData) fallingPins.emplace_back("D", held);
      instance(top, flop(name + "_falling", options, true), "falling_ff", fallingPins);
      held = falling;
    }
    if (options.latchChain) {
      auto* middle = SNLScalarNet::create(top, NLName("middle"));
      instance(top, NLDB0::getDLatch(), "first_latch", {{"D", held}, {"E", enable}, {"Q", middle}});
      instance(top, NLDB0::getDLatch(), "second_latch", {{"D", middle}, {"E", enable}, {"Q", output}});
    }
    if (options.clock == ClockRoute::Multiple) {
      auto* otherClock = port(top, "other_clock", SNLTerm::Direction::Input);
      auto* otherOutput = port(top, "other_out", SNLTerm::Direction::Output);
      instance(top, flop(name + "_other", options), "other_ff",
          {{"D", data}, {"C", otherClock}, {"R", reset}, {"Q", otherOutput}});
    }
    return top;
  }
  SNLDesign* resetObservationCircuit(const std::string& name, bool exposeReset) {
    auto* top = SNLDesign::create(designs_, SNLDesign::Type::Standard, NLName(name));
    auto* clock = port(top, "clock", SNLTerm::Direction::Input);
    auto* reset = port(top, "reset", SNLTerm::Direction::Input);
    auto* output = port(top, "out", SNLTerm::Direction::Output);
    if (exposeReset) top->getScalarTerm(NLName("out"))->setNet(reset);
    else output->setType(SNLNet::Type::Assign0);
    auto* hidden = SNLScalarNet::create(top, NLName("hidden"));
    CircuitOptions options;
    options.constantData = true;
    // A real sequential cell supplies clock discovery, but its data state is
    // irrelevant to the output-mask and permanent-reset-clamp property.
    instance(top, flop(name + "_ff", options), "ff",
        {{"C", clock}, {"R", reset}, {"Q", hidden}});
    return top;
  }
  SequentialEquivalenceResult compare(SNLDesign* first, SNLDesign* second,
      SecEngine engine = SecEngine::KInduction, SecEncoding encoding = SecEncoding::Binary,
      SecResetSpec reset = {2, {{"reset", true}}}) {
    ScopedSupportOptions scope(settings());
    const auto left = SequentialDesignModel::extract(first);
    const auto right = SequentialDesignModel::extract(second);
    EXPECT_FALSE(left.hasUnsupportedFeatures());
    EXPECT_FALSE(right.hasUnsupportedFeatures());
    for (const auto* extracted : {&left, &right}) {
      for (const auto& skipped : extracted->skippedObservedOutputs)
        ADD_FAILURE() << (extracted == &left ? "first: " : "second: ")
                      << extracted->connectivitySkipInfoByKey.at(skipped).detail;
    }
    SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
                                           engine, encoding, reset);
    return strategy.runExtractedModels(left, right, 48);
  }
  NLLibrary* designs_ = nullptr;
  NLLibrary* primitives_ = nullptr;
};

class LatchResetEngineTests : public LatchResetIntegrationTests,
                             public ::testing::WithParamInterface<std::tuple<SecEngine, SecEncoding>> {};

TEST_P(LatchResetEngineTests, SynchronousResetThenIndependentEnableLatchChain) {
  auto* first = circuit("first");
  auto* second = circuit("second");
  const auto [engine, encoding] = GetParam();
  const auto result = compare(first, second, engine, encoding);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_P(LatchResetEngineTests, RealPostResetMismatchIsNotMaskedForever) {
  auto* first = circuit("first");
  CircuitOptions options;
  options.invertData = true;
  auto* second = circuit("second", options);
  const auto [engine, encoding] = GetParam();
  const auto result = compare(first, second, engine, encoding);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_P(LatchResetEngineTests, ResetOnlyOutputDifferenceIsMaskedAndResetStaysInactiveAfterRelease) {
  auto* first = resetObservationCircuit("first", false);
  auto* second = resetObservationCircuit("second", true);
  const auto [engine, encoding] = GetParam();
  // Two cycles expose active reset at the first prefix boundary: missing
  // masking would fail there. After release, every selector/value sequence is
  // checked, including attempts to reassert reset. The irrelevant FF data
  // state is intentionally outside this focused property's observable cone.
  const auto result = compare(first, second, engine, encoding);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchResetIntegrationTests, ResetSensitivePhysicalSequentialOutputIsMasked) {
  CircuitOptions options;
  options.latchChain = false;
  options.constantData = true;
  auto* first = circuit("first", options);
  options.invertOnlyDuringReset = true;
  auto* second = circuit("second", options);
  // This richer output depends on both hidden storage and current reset. The
  // existing dual-rail IMC Craig interpolant-growth budget cannot prove this
  // fixture; the smaller reset-observation property above is required to pass
  // all six engine/encoding combinations. Do not turn a resource-limited
  // result into an expected equivalence or silently relax that assertion.
  for (auto engine : {SecEngine::KInduction, SecEngine::Imc, SecEngine::Pdr}) {
    for (auto encoding : {SecEncoding::Binary, SecEncoding::DualRailSteady}) {
      if (engine == SecEngine::Imc && encoding == SecEncoding::DualRailSteady) continue;
      SCOPED_TRACE(::testing::Message() << "engine=" << int(engine) << " encoding=" << int(encoding));
      const auto result = compare(first, second, engine, encoding, {1, {{"reset", true}}});
      EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
      EXPECT_EQ(result.coveredOutputs, 1u);
    }
  }
}

INSTANTIATE_TEST_SUITE_P(AllEnginesAndEncodings, LatchResetEngineTests,
    ::testing::Combine(::testing::Values(SecEngine::KInduction, SecEngine::Imc, SecEngine::Pdr),
                       ::testing::Values(SecEncoding::Binary, SecEncoding::DualRailSteady)));

TEST_F(LatchResetIntegrationTests, BufferedAndInvertedClockRoutesSupportFullCycles) {
  size_t index = 0;
  for (auto route : {ClockRoute::Buffered, ClockRoute::Inverted}) {
    CircuitOptions options;
    options.clock = route;
    auto* first = circuit("first" + std::to_string(index), options);
    auto* second = circuit("second" + std::to_string(index++), options);
    const auto result = compare(first, second);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
    EXPECT_EQ(result.coveredOutputs, 1u);
  }
}

TEST_F(LatchResetIntegrationTests, MixedEdgeFlopsUseTheSameExternalCarrier) {
  CircuitOptions options;
  options.mixedEdges = true;
  const auto result = compare(circuit("first", options), circuit("second", options));
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchResetIntegrationTests, OneResetCycleActuallyVisitsBothEdgesAndClearsStoredOnes) {
  CircuitOptions options;
  options.latchChain = false;
  options.mixedEdges = true;
  auto* top = circuit("top", options);
  ScopedSupportOptions scope(settings(true));
  const auto original = SequentialDesignModel::extract(top);
  ASSERT_EQ(original.observedOutputs.size(), 1u);
  ASSERT_TRUE(original.skippedObservedOutputs.empty());
  const auto wrapped = adaptResetCycles(original, {1, {{"reset", true}}});
  ASSERT_TRUE(wrapped.model.has_value()) << wrapped.error;
  const auto& model = *wrapped.model;
  std::unordered_map<size_t, bool> state;
  for (const auto& key : model.stateBits)
    state[model.inputVarByKey.at(key)] = model.initialStateValueByKey.at(key);
  const auto step = [&] {
    auto environment = state;
    // Selector zero and value zero request a clock-low stutter. The prefix
    // itself must supply both edges; the environment provides no clock tick.
    for (const auto& key : model.environmentInputs)
      environment[model.inputVarByKey.at(key)] = false;
    const bool output = model.observedOutputExprByKey.at(model.observedOutputs.front())->evaluate(environment);
    for (const auto& key : model.stateBits)
      state[model.inputVarByKey.at(key)] = model.nextStateExprByStateKey.at(key)->evaluate(environment);
    return output;
  };
  EXPECT_FALSE(step());  // Reset cycle is masked.
  EXPECT_FALSE(step());  // Normal observation: falling-edge state really cleared.
  EXPECT_FALSE(step());
}

TEST_F(LatchResetIntegrationTests, AsynchronousResetCanPropagateThroughALatch) {
  CircuitOptions options;
  options.asyncThroughLatch = true;
  const auto result = compare(circuit("first", options), circuit("second", options));
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 1u);
}

TEST_F(LatchResetIntegrationTests, ResetCycleSamplesAllFreeInputsRatherThanOnlyOneInput) {
  CircuitOptions options;
  options.latchChain = false;
  options.ignoreReset = true;
  options.andDataEnable = true;
  auto* top = circuit("top", options);
  ScopedSupportOptions scope(settings());
  const auto original = SequentialDesignModel::extract(top);
  ASSERT_EQ(original.observedOutputs.size(), 1u);
  ASSERT_TRUE(original.eventResetInterface);
  const auto wrapped = adaptResetCycles(original, {1, {{"reset", true}}});
  ASSERT_TRUE(wrapped.model) << wrapped.error;
  const auto& model = *wrapped.model;
  std::unordered_map<size_t, bool> state;
  for (const auto& key : model.stateBits)
    state[model.inputVarByKey.at(key)] = model.initialStateValueByKey.at(key);
  const auto step = [&] {
    auto environment = state;
    for (const auto& key : model.environmentInputs)
      environment[model.inputVarByKey.at(key)] = false;
    const auto& interface = *original.eventResetInterface;
    for (size_t i = 0; i < interface.inputNames.size(); ++i) {
      if (interface.inputNames[i] == "data[0]" || interface.inputNames[i] == "enable[0]")
        environment[model.inputVarByKey.at(interface.inputKeys[i])] = true;
    }
    const bool output = model.observedOutputExprByKey.at(model.observedOutputs.front())->evaluate(environment);
    for (const auto& key : model.stateBits)
      state[model.inputVarByKey.at(key)] = model.nextStateExprByStateKey.at(key)->evaluate(environment);
    return output;
  };
  EXPECT_FALSE(step());  // Masked reset cycle, with both free inputs driven high.
  EXPECT_TRUE(step());   // AND(data,enable) was captured on the forced rising edge.
}

TEST_F(LatchResetIntegrationTests, MultipleIndependentClockRootsHaveASpecificDiagnostic) {
  CircuitOptions options;
  options.clock = ClockRoute::Multiple;
  const auto result = compare(circuit("first", options), circuit("second", options));
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported) << result.reason;
  EXPECT_NE(result.reason.find("multiple independent"), std::string::npos) << result.reason;
}

TEST_F(LatchResetIntegrationTests, ConstantAndGatedClocksAreNotInventedAsPeriodicCarriers) {
  size_t index = 0;
  for (auto route : {ClockRoute::Constant, ClockRoute::Gated}) {
    CircuitOptions options;
    options.clock = route;
    const auto suffix = std::to_string(index++);
    const auto result = compare(circuit("first" + suffix, options),
                                circuit("second" + suffix, options));
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported) << result.reason;
    EXPECT_NE(result.reason.find(route == ClockRoute::Constant ? "constant" : "carrier"),
              std::string::npos) << result.reason;
  }
}

TEST_F(LatchResetIntegrationTests, MultipleResetPortsAreExplicitlyUnsupported) {
  const auto result = compare(circuit("first"), circuit("second"), SecEngine::KInduction,
      SecEncoding::Binary, {2, {{"reset", true}, {"enable", false}}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported) << result.reason;
  EXPECT_NE(result.reason.find("reset"), std::string::npos) << result.reason;
  EXPECT_EQ(result.reason.find("cannot use clock-cycle"), std::string::npos) << result.reason;
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
