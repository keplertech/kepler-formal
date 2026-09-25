// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
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
#include "latch/LatchNetlistAdapter.h"
#include "latch/LatchSupportOptions.h"
#include "latch/NajaEventPrimitive.h"
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

Expression stateExpression(bool invert = false) {
  Expression expression;
  auto root = expression.addState(0);
  if (invert) root = expression.addOperation(Operator::Not, {root});
  expression.root = root;
  return expression;
}

class LatchNetlistAdapterTests : public ::testing::Test {
 protected:
  void SetUp() override {
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    designs = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
    primitives = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("primitives"));
  }
  void TearDown() override {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) universe->destroy();
    BoolExprCache::destroy();
  }

  SupportOptions options(bool single = true) {
    SupportOptions result;
    result.enabled = true;
    result.singleInputChange = single;
    result.initialInputs = false;
    result.initialStorage = false;
    result.workers = 2;
    return result;
  }

  SNLScalarNet* port(SNLDesign* top, const char* name, SNLTerm::Direction direction) {
    auto* term = SNLScalarTerm::create(top, direction, NLName(name));
    auto* net = SNLScalarNet::create(top, NLName(name));
    term->setNet(net);
    return net;
  }

  SNLDesign* directTop(const char* name, SNLDesign* primitive, const char* control = "E") {
    auto* top = SNLDesign::create(designs, SNLDesign::Type::Standard, NLName(name));
    auto* data = port(top, "data", SNLTerm::Direction::Input);
    auto* enable = port(top, "enable", SNLTerm::Direction::Input);
    auto* cell = SNLInstance::create(top, primitive, NLName("cell"));
    cell->getInstTerm(primitive->getScalarTerm(NLName("D")))->setNet(data);
    cell->getInstTerm(primitive->getScalarTerm(NLName(control)))->setNet(enable);
    size_t outputIndex = 0;
    for (auto* term : primitive->getBitTerms()) {
      if (term->getDirection() == SNLTerm::Direction::Output) {
        auto* output = port(top, outputIndex++ ? "out_n" : "out", SNLTerm::Direction::Output);
        cell->getInstTerm(term)->setNet(output);
      }
    }
    return top;
  }

  SNLDesign* explicitLatch(const char* name, bool inverted = false, bool pair = false,
                           bool gatedOutput = false) {
    auto* cell = SNLDesign::create(primitives, SNLDesign::Type::Primitive, NLName(name));
    auto* data = SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("D"));
    auto* enable = SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("E"));
    auto* output = SNLScalarTerm::create(cell, SNLTerm::Direction::Output, NLName("Q"));
    Modeling::SequentialModel model;
    model.kind = Modeling::SequentialModel::Kind::Latch;
    model.clockedOn = termExpression(enable);
    if (gatedOutput) {
      model.clockedOn.root = model.clockedOn.addOperation(Operator::Not, {model.clockedOn.root});
    }
    Modeling::SequentialState state;
    state.nextState = termExpression(data);
    model.states.push_back(state);
    auto outputExpression = stateExpression(inverted);
    if (gatedOutput) {
      const auto phase = outputExpression.addTerm(enable);
      outputExpression.root = outputExpression.addOperation(Operator::And, {outputExpression.root, phase});
    }
    model.outputs.push_back({output, outputExpression});
    if (pair) {
      auto* complement = SNLScalarTerm::create(cell, SNLTerm::Direction::Output, NLName("QN"));
      model.outputs.push_back({complement, stateExpression(!inverted)});
    }
    Modeling::setSequentialModel(cell, model);
    return cell;
  }

  SNLDesign* inverterModel() {
    auto* cell = SNLDesign::create(primitives, SNLDesign::Type::Primitive, NLName("inverter"));
    auto* input = SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("A"));
    auto* output = SNLScalarTerm::create(cell, SNLTerm::Direction::Output, NLName("Y"));
    Modeling::addCombinatorialArcs({input}, {output});
    Modeling::setTruthTable(cell, SNLTruthTable(1, 1, SNLTruthTable::fullDependencies(1)));
    return cell;
  }

  std::unordered_map<size_t, bool> initialState(const SequentialDesignModel& model) {
    std::unordered_map<size_t, bool> result;
    for (const auto& key : model.stateBits) {
      EXPECT_EQ(model.initialStateValueByKey.count(key), 1u);
      result[model.inputVarByKey.at(key)] = model.initialStateValueByKey.at(key);
    }
    return result;
  }

  std::map<std::string, bool> step(const SequentialDesignModel& model,
      std::unordered_map<size_t, bool>& state, size_t selector, bool value) {
    auto environment = state;
    for (const auto& key : model.environmentInputs) {
      const auto& name = model.displayNameByKey.at(key);
      bool bit = false;
      if (name == "$event.value") bit = value;
      else if (name.find("$event.select[") == 0) {
        const auto index = std::stoul(name.substr(std::string("$event.select[").size()));
        bit = (selector >> index) & 1;
      }
      environment[model.inputVarByKey.at(key)] = bit;
    }
    std::map<std::string, bool> outputs;
    for (const auto& key : model.observedOutputs)
      outputs[model.displayNameByKey.at(key)] = model.observedOutputExprByKey.at(key)->evaluate(environment);
    for (const auto& key : model.stateBits)
      state[model.inputVarByKey.at(key)] = model.nextStateExprByStateKey.at(key)->evaluate(environment);
    return outputs;
  }

  void expectPublishedSupport(const SequentialDesignModel& model) {
    std::set<size_t> variables;
    for (const auto& [key, variable] : model.inputVarByKey) variables.insert(variable);
    const auto check = [&](BoolExpr* expression) {
      for (auto variable : expression->getSupportVars()) {
        if (variable > 1) EXPECT_TRUE(variables.count(variable)) << variable;
      }
    };
    for (const auto& [key, expression] : model.nextStateExprByStateKey) check(expression);
    for (const auto& [key, expression] : model.observedOutputExprByKey) check(expression);
  }

  NLLibrary* designs = nullptr;
  NLLibrary* primitives = nullptr;
};

TEST_F(LatchNetlistAdapterTests, DefaultOptionsLeaveExistingExtractionUntouched) {
  auto* top = directTop("top", NLDB0::getDLatch());
  EXPECT_FALSE(extractEventDesign(top, {}, 0).has_value());
  const auto model = SequentialDesignModel::extract(top);
  EXPECT_TRUE(model.observedOutputs.empty());
  EXPECT_EQ(model.skippedObservedOutputs.size(), 1u);
}

TEST_F(LatchNetlistAdapterTests, ExplicitContractRejectsUnspecifiedInitialization) {
  auto setting = options();
  setting.initialInputs.reset();
  ScopedSupportOptions scope(setting);
  const auto model = SequentialDesignModel::extract(directTop("top", NLDB0::getDLatch()));
  EXPECT_TRUE(model.hasUnsupportedFeatures());
  EXPECT_TRUE(model.observedOutputs.empty());
}

TEST_F(LatchNetlistAdapterTests, ScopedOptionsRestorePreviousSemanticContract) {
  EXPECT_FALSE(supportOptions().enabled);
  {
    ScopedSupportOptions outer(options());
    EXPECT_TRUE(supportOptions().enabled);
    {
      ScopedSupportOptions inner(SupportOptions{});
      EXPECT_FALSE(supportOptions().enabled);
    }
    EXPECT_TRUE(supportOptions().enabled);
    EXPECT_TRUE(supportOptions().singleInputChange);
  }
  EXPECT_FALSE(supportOptions().enabled);
}

TEST_F(LatchNetlistAdapterTests, Db0LatchFollowsOpenDataAndRetainsClosingValue) {
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", NLDB0::getDLatch()));
  ASSERT_FALSE(model.hasUnsupportedFeatures());
  ASSERT_EQ(model.observedOutputs.size(), 1u);
  EXPECT_TRUE(model.skippedObservedOutputs.empty());
  EXPECT_FALSE(model.stateBits.empty());
  expectPublishedSupport(model);
  auto state = initialState(model);
  EXPECT_FALSE(step(model, state, 0, true).at("out[0]"));  // Closed.
  EXPECT_TRUE(step(model, state, 1, true).at("out[0]"));   // Open.
  EXPECT_FALSE(step(model, state, 0, false).at("out[0]"));
  EXPECT_TRUE(step(model, state, 0, true).at("out[0]"));
  EXPECT_TRUE(step(model, state, 1, false).at("out[0]"));  // Close.
  EXPECT_TRUE(step(model, state, 0, false).at("out[0]"));  // Hold.
}

TEST_F(LatchNetlistAdapterTests, AnyChangeLatchRaceRemainsOpaqueRatherThanSelectingOrder) {
  ScopedSupportOptions scope(options(false));
  const auto model = SequentialDesignModel::extract(directTop("top", NLDB0::getDLatch()));
  ASSERT_FALSE(model.hasUnsupportedFeatures());
  EXPECT_TRUE(model.observedOutputs.empty());
  ASSERT_EQ(model.skippedObservedOutputs.size(), 1u);
  const auto& detail = model.connectivitySkipInfoByKey.at(model.skippedObservedOutputs.front()).detail;
  EXPECT_NE(detail.find("order"), std::string::npos);
}

TEST_F(LatchNetlistAdapterTests, Db0FlipFlopConsumesClockEventOnlyOnce) {
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", NLDB0::getDFF(), "C"));
  ASSERT_EQ(model.observedOutputs.size(), 1u);
  auto state = initialState(model);
  EXPECT_FALSE(step(model, state, 0, true).at("out[0]"));
  EXPECT_TRUE(step(model, state, 1, true).at("out[0]"));
  EXPECT_TRUE(step(model, state, 0, false).at("out[0]"));
  EXPECT_TRUE(step(model, state, 1, false).at("out[0]"));
  EXPECT_FALSE(step(model, state, 1, true).at("out[0]"));
}

TEST_F(LatchNetlistAdapterTests, InvertedPhysicalOutputDoesNotInvertRememberedStorage) {
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", explicitLatch("inverse", true)));
  ASSERT_EQ(model.observedOutputs.size(), 1u);
  auto state = initialState(model);
  EXPECT_TRUE(step(model, state, 3, false).at("out[0]"));
  EXPECT_TRUE(step(model, state, 0, true).at("out[0]"));
  EXPECT_FALSE(step(model, state, 1, true).at("out[0]"));
  EXPECT_FALSE(step(model, state, 1, false).at("out[0]"));
  EXPECT_FALSE(step(model, state, 0, false).at("out[0]"));
}

TEST_F(LatchNetlistAdapterTests, ComplementaryOutputsCommitConsistentlyAndHavePublishedSymbols) {
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", explicitLatch("pair", false, true)));
  ASSERT_EQ(model.observedOutputs.size(), 2u);
  expectPublishedSupport(model);
  auto state = initialState(model);
  for (const auto event : {std::pair<size_t, bool>{3, false}, {0, true}, {1, true}, {1, false}, {0, false}}) {
    const auto outputs = step(model, state, event.first, event.second);
    EXPECT_NE(outputs.at("out[0]"), outputs.at("out_n[0]"));
  }
}

TEST_F(LatchNetlistAdapterTests, InputDependentPhysicalOutputTracksClockGatePhase) {
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", explicitLatch("icg", false, false, true)));
  ASSERT_EQ(model.observedOutputs.size(), 1u);
  auto state = initialState(model);
  EXPECT_FALSE(step(model, state, 0, true).at("out[0]"));  // Enable remembered during low carrier.
  EXPECT_TRUE(step(model, state, 1, true).at("out[0]"));
  EXPECT_TRUE(step(model, state, 0, false).at("out[0]"));  // Held during high carrier.
  EXPECT_FALSE(step(model, state, 1, false).at("out[0]"));
  EXPECT_FALSE(step(model, state, 1, true).at("out[0]"));
}

TEST_F(LatchNetlistAdapterTests, LatchGeneratedClockDrivesFlipFlopThroughInternalEvents) {
  auto* top = SNLDesign::create(designs, SNLDesign::Type::Standard, NLName("clock_gated"));
  auto* clock = port(top, "clock", SNLTerm::Direction::Input);
  auto* data = port(top, "data", SNLTerm::Direction::Input);
  auto* gate = port(top, "gate", SNLTerm::Direction::Input);
  auto* output = port(top, "out", SNLTerm::Direction::Output);
  auto* generatedClock = SNLScalarNet::create(top, NLName("generated_clock"));
  auto* gateModel = explicitLatch("clock_gate", false, false, true);
  auto* gateCell = SNLInstance::create(top, gateModel, NLName("gate_cell"));
  gateCell->getInstTerm(gateModel->getScalarTerm(NLName("D")))->setNet(gate);
  gateCell->getInstTerm(gateModel->getScalarTerm(NLName("E")))->setNet(clock);
  gateCell->getInstTerm(gateModel->getScalarTerm(NLName("Q")))->setNet(generatedClock);
  auto* flop = SNLInstance::create(top, NLDB0::getDFF(), NLName("flop"));
  flop->getInstTerm(NLDB0::getDFFData())->setNet(data);
  flop->getInstTerm(NLDB0::getDFFClock())->setNet(generatedClock);
  flop->getInstTerm(NLDB0::getDFFOutput())->setNet(output);

  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(top);
  ASSERT_FALSE(model.hasUnsupportedFeatures());
  ASSERT_EQ(model.observedOutputs.size(), 1u);
  EXPECT_TRUE(model.skippedObservedOutputs.empty());
  expectPublishedSupport(model);
  auto state = initialState(model);
  // Selector order follows the common alphabetical interface: clock, data, gate.
  EXPECT_FALSE(step(model, state, 1, true).at("out[0]"));
  EXPECT_FALSE(step(model, state, 2, true).at("out[0]"));
  EXPECT_TRUE(step(model, state, 0, true).at("out[0]"));  // Generated rising edge captures 1.
  EXPECT_TRUE(step(model, state, 1, false).at("out[0]")); // Data change is not another edge.
  EXPECT_TRUE(step(model, state, 2, false).at("out[0]"));
  EXPECT_TRUE(step(model, state, 0, false).at("out[0]"));
  EXPECT_TRUE(step(model, state, 0, true).at("out[0]"));  // Gate disabled: no generated edge.
  EXPECT_TRUE(step(model, state, 2, true).at("out[0]"));  // High-phase gate change stays held.
  EXPECT_TRUE(step(model, state, 0, false).at("out[0]"));
  EXPECT_FALSE(step(model, state, 0, true).at("out[0]")); // Next generated edge captures 0.
}

TEST_F(LatchNetlistAdapterTests, StateDependentLatchDataIsOpaqueNotFalseSettling) {
  auto* primitive = explicitLatch("hidden_feedback");
  auto sequential = Modeling::getSequentialModel(primitive);
  sequential.clockedOn = Expression{};
  sequential.clockedOn.root = sequential.clockedOn.addConstant(true);
  sequential.states.front().nextState = stateExpression(true);
  Modeling::setSequentialModel(primitive, sequential);
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", primitive));
  EXPECT_TRUE(model.observedOutputs.empty());
  EXPECT_EQ(model.skippedObservedOutputs.size(), 1u);
}

TEST_F(LatchNetlistAdapterTests, StateDependentAsyncControlIsOpaqueForLatchAndFlop) {
  for (bool flop : {false, true}) {
    auto* primitive = explicitLatch(flop ? "hidden_ff_clear" : "hidden_latch_clear");
    auto sequential = Modeling::getSequentialModel(primitive);
    if (flop) sequential.kind = Modeling::SequentialModel::Kind::FlipFlop;
    sequential.states.front().clear = stateExpression(true);
    Modeling::setSequentialModel(primitive, sequential);
    ScopedSupportOptions scope(options());
    const auto model = SequentialDesignModel::extract(directTop(flop ? "ff_top" : "latch_top", primitive));
    EXPECT_TRUE(model.observedOutputs.empty());
    EXPECT_EQ(model.skippedObservedOutputs.size(), 1u);
  }
}

TEST_F(LatchNetlistAdapterTests, UnsupportedLoopDoesNotDiscardIndependentSupportedOutput) {
  auto* top = directTop("top", NLDB0::getDLatch());
  auto* enable = top->getScalarNet(NLName("enable"));
  auto* loop = SNLScalarNet::create(top, NLName("loop"));
  auto* inverted = SNLScalarNet::create(top, NLName("inverted"));
  auto* bad = port(top, "bad", SNLTerm::Direction::Output);
  auto* invModel = inverterModel();
  auto* inv = SNLInstance::create(top, invModel, NLName("inv"));
  inv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(loop);
  inv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(inverted);
  auto* latch = SNLInstance::create(top, NLDB0::getDLatch(), NLName("bad_latch"));
  latch->getInstTerm(NLDB0::getDLatchEnable())->setNet(enable);
  latch->getInstTerm(NLDB0::getDLatchData())->setNet(inverted);
  latch->getInstTerm(NLDB0::getDLatchOutput())->setNet(loop);
  top->getScalarTerm(NLName("bad"))->setNet(loop);
  (void)bad;
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(top);
  ASSERT_FALSE(model.hasUnsupportedFeatures());
  ASSERT_EQ(model.observedOutputs.size(), 1u);
  ASSERT_EQ(model.skippedObservedOutputs.size(), 1u);
  EXPECT_EQ(model.displayNameByKey.at(model.observedOutputs.front()), "out[0]");
  EXPECT_EQ(model.displayNameByKey.at(model.skippedObservedOutputs.front()), "bad[0]");
}

TEST_F(LatchNetlistAdapterTests, UnsupportedControlSourceMakesDependentComponentOpaque) {
  auto* unknown = SNLDesign::create(primitives, SNLDesign::Type::Primitive, NLName("opaque_source"));
  auto* sourceData = SNLScalarTerm::create(unknown, SNLTerm::Direction::Input, NLName("D"));
  auto* sourceOut = SNLScalarTerm::create(unknown, SNLTerm::Direction::Output, NLName("Q"));
  auto* top = directTop("top", NLDB0::getDLatch());
  auto* net = SNLScalarNet::create(top, NLName("control"));
  auto* source = SNLInstance::create(top, unknown, NLName("source"));
  source->getInstTerm(sourceData)->setNet(top->getScalarNet(NLName("enable")));
  source->getInstTerm(sourceOut)->setNet(net);
  top->getInstance(NLName("cell"))->getInstTerm(NLDB0::getDLatchEnable())->setNet(net);
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(top);
  EXPECT_TRUE(model.observedOutputs.empty());
  ASSERT_EQ(model.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(model.connectivitySkipInfoByKey.at(model.skippedObservedOutputs[0]).detail.find("source"), std::string::npos);
}

TEST_F(LatchNetlistAdapterTests, ExtractionRestoresBorrowedTopAndFlattenedGraph) {
  auto* saved = directTop("saved", NLDB0::getDLatch());
  auto* target = directTop("target", NLDB0::getDLatch());
  NLUniverse::get()->setTopDesign(saved);
  auto* dnl = naja::DNL::get();
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(target);
  EXPECT_EQ(model.observedOutputs.size(), 1u);
  EXPECT_EQ(NLUniverse::get()->getTopDesign(), saved);
  EXPECT_EQ(naja::DNL::get(), dnl);
}

TEST_F(LatchNetlistAdapterTests, ExtractionRestoresNondefaultSourceOrderingIds) {
  auto* primitive = explicitLatch("ordering");
  auto* top = directTop("top", primitive);
  std::vector<std::pair<SNLBitTerm*, NLID::DesignObjectID>> orders;
  NLID::DesignObjectID order = 51;
  for (auto* design : {top, primitive}) {
    for (auto* term : design->getBitTerms()) {
      term->setOrderID(order++);
      orders.emplace_back(term, term->getOrderID());
    }
  }
  auto* instance = top->getInstance(NLName("cell"));
  instance->setOrderID(77);
  const auto* previousDb = NLUniverse::get()->getTopDB();
  const auto* previousTop = top->getDB()->getTopDesign();
  ScopedSupportOptions scope(options());
  const auto model = extractEventDesign(top, {}, 0);
  ASSERT_TRUE(model.has_value());
  EXPECT_EQ(model->observedOutputs.size(), 1u);
  for (const auto& [term, savedOrder] : orders) EXPECT_EQ(term->getOrderID(), savedOrder);
  EXPECT_EQ(instance->getOrderID(), 77u);
  EXPECT_EQ(NLUniverse::get()->getTopDB(), previousDb);
  EXPECT_EQ(top->getDB()->getTopDesign(), previousTop);
}

TEST_F(LatchNetlistAdapterTests, LatchFreeSideUsesSameEventEnvironmentForComparison) {
  auto* top = SNLDesign::create(designs, SNLDesign::Type::Standard, NLName("wire"));
  auto* data = port(top, "data", SNLTerm::Direction::Input);
  port(top, "enable", SNLTerm::Direction::Input);
  auto* output = SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  output->setNet(data);
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(top);
  ASSERT_EQ(model.observedOutputs.size(), 1u);
  EXPECT_FALSE(model.stateBits.empty());
  auto state = initialState(model);
  EXPECT_TRUE(step(model, state, 0, true).at("out[0]"));
  EXPECT_TRUE(step(model, state, 1, true).at("out[0]"));
}

TEST_F(LatchNetlistAdapterTests, ProofEnginesAcceptSelfEquivalentLatchInBothEncodings) {
  auto* first = directTop("first", NLDB0::getDLatch());
  auto* second = directTop("second", NLDB0::getDLatch());
  ScopedSupportOptions scope(options());
  for (auto engine : {SecEngine::KInduction, SecEngine::Imc, SecEngine::Pdr}) {
    for (auto encoding : {SecEncoding::Binary, SecEncoding::DualRailSteady}) {
      SequentialEquivalenceStrategy strategy(first, second, Config::SolverType::KISSAT, engine, encoding);
      const auto result = strategy.run(16);
      EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
      EXPECT_EQ(result.coveredOutputs, 1u);
    }
  }
}

TEST_F(LatchNetlistAdapterTests, ProofEnginesFindDifferentPhysicalLatchOutputs) {
  auto* first = directTop("first", NLDB0::getDLatch());
  auto* second = directTop("second", explicitLatch("inverted", true));
  ScopedSupportOptions scope(options());
  for (auto engine : {SecEngine::KInduction, SecEngine::Pdr}) {
    SequentialEquivalenceStrategy strategy(first, second, Config::SolverType::KISSAT, engine, SecEncoding::Binary);
    const auto result = strategy.run(8);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
    EXPECT_NE(result.reason.find("event transaction"), std::string::npos);
    EXPECT_NE(result.reason.find("Event contract: boolean-epochs-v1"), std::string::npos);
    EXPECT_EQ(result.reason.find("at cycle"), std::string::npos);
  }
}

TEST_F(LatchNetlistAdapterTests, SameSizeDifferentNamedInterfacesCannotShareSelectorAccidentally) {
  auto* first = directTop("first", NLDB0::getDLatch());
  auto* second = directTop("second", NLDB0::getDLatch());
  second->getScalarTerm(NLName("data"))->setName(NLName("different_data"));
  ScopedSupportOptions scope(options());
  SequentialEquivalenceStrategy strategy(first, second, Config::SolverType::KISSAT,
      SecEngine::KInduction, SecEncoding::Binary);
  EXPECT_THROW(strategy.run(8), std::runtime_error);
}

TEST_F(LatchNetlistAdapterTests, CopiedPrimitiveReactionsSurviveReleaseOfNajaNetlist) {
  auto* top = directTop("top", explicitLatch("copied", false, true));
  auto* instance = top->getInstance(NLName("cell"));
  std::map<const SNLBitTerm*, size_t> netByTerm;
  Network network;
  for (auto* term : instance->getModel()->getBitTerms()) {
    netByTerm.emplace(term, network.netCount);
    if (term->getDirection() == SNLTerm::Direction::Input)
      network.externalInputs.push_back(network.netCount);
    ++network.netCount;
  }
  network.primitives.push_back(makeNajaEventPrimitive(instance, "cell", netByTerm));
  NLUniverse::get()->destroy();
  EventModel model(network, {}, 2);
  auto state = model.bootstrap({1, 1}, {{0}}, Bits(network.netCount));
  for (size_t wave = 0; !model.stable(state) && wave < 16; ++wave)
    state = model.successors(state).front();
  ASSERT_TRUE(model.stable(state));
  const auto& outputs = network.primitives.front().outputs;
  EXPECT_EQ(state.current[outputs[0]], 1);
  EXPECT_EQ(state.current[outputs[1]], 0);
}

TEST_F(LatchNetlistAdapterTests, ConflictingConstantAnnotationsStayOpaqueWithExactlyOneDriver) {
  auto* top = directTop("top", NLDB0::getDLatch());
  auto* data = top->getScalarNet(NLName("data"));
  data->setType(SNLNet::Type::Assign0);
  auto* child = SNLDesign::create(designs, SNLDesign::Type::Standard, NLName("annotation"));
  auto* input = SNLScalarTerm::create(child, SNLTerm::Direction::Input, NLName("i"));
  auto* childNet = SNLScalarNet::create(child, NLName("conflicting_one"));
  childNet->setType(SNLNet::Type::Assign1);
  input->setNet(childNet);
  auto* instance = SNLInstance::create(top, child, NLName("annotation"));
  instance->getInstTerm(input)->setNet(data);

  NLUniverse::get()->setTopDesign(top);
  const auto* dnl = naja::DNL::get();
  const auto& terminal = dnl->getTop().getTerminalFromBitTerm(top->getScalarTerm(NLName("data")));
  const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(terminal.getIsoID());
  ASSERT_EQ(iso.getType(), naja::DNL::DNLIso::AMBIGUOUS);
  ASSERT_EQ(iso.getDrivers().size(), 1u);  // Driver-count validation alone misses this case.

  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(top);
  EXPECT_TRUE(model.observedOutputs.empty());
  ASSERT_EQ(model.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(model.connectivitySkipInfoByKey.at(model.skippedObservedOutputs.front()).detail.find(
      "conflicting constant"), std::string::npos);
}

TEST_F(LatchNetlistAdapterTests, DuplicateSequentialOutputAssociationsAreNotChosenArbitrarily) {
  auto* primitive = explicitLatch("duplicate_outputs");
  auto sequential = Modeling::getSequentialModel(primitive);
  sequential.outputs.push_back({sequential.outputs.front().term, stateExpression(true)});
  Modeling::setSequentialModel(primitive, sequential);
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", primitive));
  EXPECT_TRUE(model.observedOutputs.empty());
  ASSERT_EQ(model.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(model.connectivitySkipInfoByKey.at(model.skippedObservedOutputs.front()).detail.find(
      "duplicate sequential output"), std::string::npos);
}

TEST_F(LatchNetlistAdapterTests, ReachableClearPresetToggleConflictIsOpaqueNotFalseSettling) {
  for (bool flop : {false, true}) {
    auto* primitive = explicitLatch(flop ? "toggle_flop" : "toggle_latch");
    auto sequential = Modeling::getSequentialModel(primitive);
    if (flop) sequential.kind = Modeling::SequentialModel::Kind::FlipFlop;
    sequential.states.front().clear = termExpression(primitive->getScalarTerm(NLName("D")));
    sequential.states.front().preset = termExpression(primitive->getScalarTerm(NLName("E")));
    sequential.states.front().clearPresetValue = Modeling::SequentialState::ClearPresetValue::Toggle;
    Modeling::setSequentialModel(primitive, sequential);
    ScopedSupportOptions scope(options());
    const auto model = SequentialDesignModel::extract(directTop(flop ? "ff_top" : "latch_top", primitive));
    EXPECT_TRUE(model.observedOutputs.empty());
    ASSERT_EQ(model.skippedObservedOutputs.size(), 1u);
    EXPECT_NE(model.connectivitySkipInfoByKey.at(model.skippedObservedOutputs.front()).detail.find(
        "toggle"), std::string::npos);
  }
}

TEST_F(LatchNetlistAdapterTests, UnreachableToggleConflictDoesNotRejectUsableCell) {
  auto* primitive = explicitLatch("never_conflicting");
  auto sequential = Modeling::getSequentialModel(primitive);
  sequential.states.front().clear = termExpression(primitive->getScalarTerm(NLName("D")));
  auto preset = *sequential.states.front().clear;
  preset.root = preset.addOperation(Operator::Not, {preset.root});
  sequential.states.front().preset = preset;
  sequential.states.front().clearPresetValue = Modeling::SequentialState::ClearPresetValue::Toggle;
  Modeling::setSequentialModel(primitive, sequential);
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", primitive));
  ASSERT_EQ(model.observedOutputs.size(), 1u);
  EXPECT_TRUE(model.skippedObservedOutputs.empty());
  auto state = initialState(model);
  EXPECT_TRUE(step(model, state, 3, false).at("out[0]"));
  EXPECT_FALSE(step(model, state, 0, true).at("out[0]"));
}

TEST_F(LatchNetlistAdapterTests, ExtractedModelsCannotMixEventAndCycleContracts) {
  ScopedSupportOptions scope(options());
  const auto eventModel = SequentialDesignModel::extract(directTop("top", NLDB0::getDLatch()));
  ASSERT_FALSE(eventModel.eventContract.empty());
  auto legacyModel = eventModel;
  legacyModel.eventContract.clear();
  SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
      SecEngine::KInduction, SecEncoding::Binary);
  const auto result = strategy.runExtractedModels(eventModel, legacyModel, 8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("different clock/event"), std::string::npos);
}

TEST_F(LatchNetlistAdapterTests, ExtractedModelsCannotMixBooleanInitializationContracts) {
  auto* top = directTop("top", NLDB0::getDLatch());
  SequentialDesignModel zero, one;
  {
    ScopedSupportOptions scope(options());
    zero = SequentialDesignModel::extract(top);
  }
  {
    auto setting = options();
    setting.initialStorage = true;
    ScopedSupportOptions scope(setting);
    one = SequentialDesignModel::extract(top);
  }
  SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
      SecEngine::KInduction, SecEncoding::Binary);
  const auto result = strategy.runExtractedModels(zero, one, 8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("initialization contracts"), std::string::npos);
}

TEST_F(LatchNetlistAdapterTests, ResetCyclesRequireDiscoverableClockNotAnArbitraryLatchEnable) {
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", NLDB0::getDLatch()));
  SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
      SecEngine::KInduction, SecEncoding::Binary, SecResetSpec{1, {{"data", false}}});
  const auto result = strategy.runExtractedModels(model, model, 8);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("clock"), std::string::npos);
}

TEST_F(LatchNetlistAdapterTests, Btor2ExportRetainsEventStepAndInitializationMetadata) {
  struct TemporaryDirectory {
    std::filesystem::path path;
    TemporaryDirectory() {
      for (size_t attempt = 0; attempt < 32; ++attempt) {
        path = std::filesystem::temp_directory_path() /
            ("kf-event-metadata-" + std::to_string(std::random_device{}()) + "-" + std::to_string(attempt));
        if (std::filesystem::create_directory(path)) return;
      }
      throw std::runtime_error("Cannot create event export test directory");
    }
    ~TemporaryDirectory() {
      std::error_code ignored;
      std::filesystem::remove_all(path, ignored);
    }
  } directory;
  ScopedSupportOptions scope(options());
  const auto model = SequentialDesignModel::extract(directTop("top", NLDB0::getDLatch()));
  const auto output = directory.path / "events.btor2";
  SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
      SecEngine::KInduction, SecEncoding::Binary, {}, Btor2ExportOptions{output.string(), true});
  const auto result = strategy.runExtractedModels(model, model, 8);
  ASSERT_EQ(result.status, SequentialEquivalenceStatus::Exported) << result.reason;
  std::ifstream file(output);
  ASSERT_TRUE(file.good());
  std::ostringstream contents;
  contents << file.rdbuf();
  EXPECT_NE(contents.str().find("step_semantics=" + model.eventContract), std::string::npos);
  EXPECT_NE(contents.str().find("initial_inputs=0;initial_storage=0"), std::string::npos);
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
