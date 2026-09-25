// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include <gtest/gtest.h>
#include <iomanip>
#include <set>
#include <sstream>
#include <unordered_map>
#include "BoolExprCache.h"
#include "DNL.h"
#include "NLDB.h"
#include "NLDB0.h"
#include "NLLibrary.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLBusTerm.h"
#include "SNLBusTermBit.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "latch/LatchSupportOptions.h"
#include "latch/LatchSymbolicEncoding.h"
#include "latch/NajaEventPrimitive.h"
#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using namespace naja::NL;
class LatchSymbolicIntegrationTests : public ::testing::Test {
 protected:
  void SetUp() override {
    auto* db = NLDB::create(NLUniverse::create());
    designs = NLLibrary::create(db, NLName("designs"));
    primitives = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("primitives"));
  }
  void TearDown() override {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) universe->destroy();
    BoolExprCache::destroy();
  }
  SNLScalarNet* port(SNLDesign* top, const std::string& name, SNLTerm::Direction dir) {
    auto* net = SNLScalarNet::create(top, NLName(name));
    SNLScalarTerm::create(top, dir, NLName(name))->setNet(net);
    return net;
  }
  SNLDesign* wideLatchCone(size_t width) {
    auto* gate = NLDB0::getOrCreateNInputGate(NLDB0::GateType::And, width);
    std::vector<SNLBitTerm*> pins;
    for (auto* bit : NLDB0::getGateNTerms(gate)->getBusBits()) pins.push_back(bit);
    auto* result = NLDB0::getGateSingleTerm(gate);
    auto* top = SNLDesign::create(designs, NLName("top"));
    auto* enable = port(top, "enable", SNLTerm::Direction::Input);
    auto* reduction = SNLInstance::create(top, gate, NLName("reduction"));
    reduction->getInstTerm(result)->setNet(port(top, "out", SNLTerm::Direction::Output));
    for (size_t i = 0; i < width; ++i) {
      std::ostringstream name;
      name << "data" << std::setw(4) << std::setfill('0') << i;
      auto* latch = SNLInstance::create(top, NLDB0::getDLatch(), NLName("latch" + std::to_string(i)));
      latch->getInstTerm(NLDB0::getDLatchData())->setNet(port(top, name.str(), SNLTerm::Direction::Input));
      latch->getInstTerm(NLDB0::getDLatchEnable())->setNet(enable);
      auto* q = SNLScalarNet::create(top, NLName("q" + std::to_string(i)));
      latch->getInstTerm(NLDB0::getDLatchOutput())->setNet(q);
      reduction->getInstTerm(pins[i])->setNet(q);
    }
    return top;
  }
  SupportOptions options() {
    SupportOptions result;
    result.enabled = true; result.singleInputChange = true;
    result.initialInputs = false; result.initialStorage = false; result.workers = 2;
    return result;
  }
  std::unordered_map<size_t, bool> initial(const SequentialDesignModel& model) {
    std::unordered_map<size_t, bool> result;
    for (const auto& key : model.stateBits)
      result.emplace(model.inputVarByKey.at(key), model.initialStateValueByKey.at(key));
    return result;
  }
  bool step(const SequentialDesignModel& model, std::unordered_map<size_t, bool>& state,
            size_t selector, bool value) {
    auto environment = state;
    for (const auto& key : model.environmentInputs) {
      const auto& name = model.displayNameByKey.at(key);
      bool bit = name == "$event.value" && value;
      if (name.starts_with("$event.select[")) bit = (selector >> std::stoul(name.substr(14))) & 1;
      environment[model.inputVarByKey.at(key)] = bit;
    }
    const auto output = model.observedOutputExprByKey.at(model.observedOutputs.at(0))->evaluate(environment);
    for (const auto& key : model.stateBits)
      state[model.inputVarByKey.at(key)] = model.nextStateExprByStateKey.at(key)->evaluate(environment);
    return output;
  }
  NLLibrary* designs = nullptr;
  NLLibrary* primitives = nullptr;
};

TEST_F(LatchSymbolicIntegrationTests, ModelsWideConnectedLatchConeBeyondFiniteEnumerationLimits) {
  constexpr size_t width = 32;
  auto* top = wideLatchCone(width);
  auto config = options();
  // Neither a state/input table nor bootstrap seed enumeration may prove this.
  config.limits.maxBoundaryStates = 1;
  config.limits.maxExternalBits = 1;
  config.limits.maxBootstrapSeedBits = 1;
  ScopedSupportOptions scope(config);
  const auto model = SequentialDesignModel::extract(top);
  ASSERT_TRUE(model.unsupportedReasons.empty());
  ASSERT_EQ(model.observedOutputs.size(), 1u) << (model.connectivitySkipInfoByKey.empty()
      ? "no skip diagnostic" : model.connectivitySkipInfoByKey.begin()->second.detail);
  EXPECT_TRUE(model.skippedObservedOutputs.empty());
  EXPECT_GT(model.stateBits.size(), 64u);
  auto state = initial(model);
  EXPECT_FALSE(step(model, state, width, true));
  for (size_t i = 0; i < width; ++i) EXPECT_EQ(step(model, state, i, true), i + 1 == width);
  EXPECT_TRUE(step(model, state, width, false));
  for (size_t i = 0; i < width; ++i) EXPECT_TRUE(step(model, state, i, false));
  EXPECT_FALSE(step(model, state, width, true));
}

TEST_F(LatchSymbolicIntegrationTests, ExhaustingBothCompilersLeavesConeOpaque) {
  auto* top = wideLatchCone(16);
  auto config = options();
  config.maxSymbolicNodes = 1;
  config.limits.maxExternalBits = 1;
  ScopedSupportOptions scope(config);
  const auto model = SequentialDesignModel::extract(top);
  EXPECT_TRUE(model.observedOutputs.empty());
  ASSERT_EQ(model.skippedObservedOutputs.size(), 1u);
  const auto& reason = model.connectivitySkipInfoByKey.at(model.skippedObservedOutputs[0]).detail;
  EXPECT_NE(reason.find("symbolic"), std::string::npos) << reason;
}

TEST_F(LatchSymbolicIntegrationTests, NativeNajaSymbolicCallbacksMatchConcreteStorageAndOutputs) {
  for (auto* cell : {NLDB0::getDLatch(), NLDB0::getDFF()}) {
    auto* top = SNLDesign::create(designs);
    auto* instance = SNLInstance::create(top, cell);
    std::map<const SNLBitTerm*, size_t> nets;
    for (auto* term : cell->getBitTerms()) nets.emplace(term, nets.size());
    SymbolicPrimitive symbolic;
    const auto concrete = makeNajaEventPrimitive(instance, "cell", nets, &symbolic);
    ASSERT_EQ(concrete.inputs.size(), 2u);
    SymbolicBits old{BoolExpr::Var(2)}, before{BoolExpr::Var(3), BoolExpr::Var(4)}, now{BoolExpr::Var(5), BoolExpr::Var(6)};
    for (bool boot : {false, true}) for (int pin = -1; pin < 2; ++pin) {
      const auto changed = pin < 0 ? std::optional<size_t>{} : size_t(pin);
      const auto reaction = symbolic.react(old, before, now, changed, boot);
      for (size_t code = 0; code < 32; ++code) {
        std::unordered_map<size_t, bool> environment;
        for (size_t i = 0; i < 5; ++i) environment.emplace(i + 2, (code >> i) & 1);
        const auto expected = concrete.react({uint8_t(code & 1)},
            {uint8_t((code >> 1) & 1), uint8_t((code >> 2) & 1)},
            {uint8_t((code >> 3) & 1), uint8_t((code >> 4) & 1)}, changed, boot);
        EXPECT_EQ(reaction.error->evaluate(environment), expected.error);
        EXPECT_EQ(reaction.storage[0]->evaluate(environment), expected.storage[0]);
        EXPECT_EQ(reaction.outputs[0]->evaluate(environment), expected.outputs[0]);
      }
    }
  }
}

TEST_F(LatchSymbolicIntegrationTests, MacroEncodingSubstitutesSimultaneouslyAndRejectsHiddenChoices) {
  Network network{1, {0}, {}, {}};
  SymbolicMacro macro;
  macro.externalInputNets = network.externalInputs;
  macro.stateSymbols = {2}; macro.inputSymbols = {3}; macro.initialState = {0};
  macro.nextState = {BoolExpr::Var(3)}; macro.observedNets = macro.nextState;
  const auto encoded = encodeSymbolicMacro(macro, network, {BoolExpr::Var(3)}, {BoolExpr::Var(2)}, false);
  EXPECT_EQ(encoded.nextState[0], BoolExpr::Var(2));
  macro.nextState[0] = BoolExpr::Var(99);
  EXPECT_THROW(encodeSymbolicMacro(macro, network, {BoolExpr::Var(3)}, {BoolExpr::Var(2)}, false), std::invalid_argument);
}

TEST_F(LatchSymbolicIntegrationTests, EncoderCannotBroadenCertificateOrRelabelInputs) {
  Network network{2, {0}, {}, {}};
  SymbolicMacro macro;
  macro.stateSymbols = {2, 3}; macro.inputSymbols = {4}; macro.initialState = {0, 0};
  macro.nextState = {BoolExpr::Var(4), BoolExpr::Var(3)};
  macro.observedNets = macro.nextState;
  macro.externalInputNets = {0}; macro.singleExternalInputChange = true;
  const SymbolicBits state{BoolExpr::Var(10), BoolExpr::Var(11)}, input{BoolExpr::Var(12)};
  EXPECT_THROW(encodeSymbolicMacro(macro, network, state, input, false), std::invalid_argument);
  const SymbolicBits selector{BoolExpr::Var(13)};
  EXPECT_NO_THROW(encodeSymbolicMacro(macro, network, state, input, true, selector, BoolExpr::Var(14), {0}));
  network.externalInputs = {1};
  EXPECT_THROW(encodeSymbolicMacro(macro, network, state, input, true, selector, BoolExpr::Var(14), {0}), std::invalid_argument);
  network.externalInputs = {0};
  EXPECT_THROW(encodeSymbolicMacro(macro, network, state, input, true, {nullptr}, BoolExpr::Var(14), {0}), std::invalid_argument);
}

TEST_F(LatchSymbolicIntegrationTests, InitialRelationPreservesAllOriginsAndSharedInputCorrelation) {
  SymbolicMacro macro;
  macro.stateSymbols = {2, 3};
  macro.inputSymbols = {4};
  macro.initialParameterSymbols = {5, 6};
  macro.initialStateExpressions = {BoolExpr::Var(5), BoolExpr::Xor(BoolExpr::Var(5), BoolExpr::Var(6))};
  const SymbolicBits state{BoolExpr::Var(10), BoolExpr::Var(11)};
  const SymbolicBits origins{BoolExpr::Var(12), BoolExpr::Var(13)};
  auto* relation = encodeSymbolicInitialRelation(macro, state, origins);
  for (size_t code = 0; code < 16; ++code) {
    const bool first = code & 1, second = (code >> 1) & 1;
    const bool input = (code >> 2) & 1, stored = (code >> 3) & 1;
    EXPECT_EQ(relation->evaluate({{10, first}, {11, second}, {12, input}, {13, stored}}),
        first == input && second == (input != stored));
  }
  // Two components can share input origin 12 while their storage origins 13/14
  // stay distinct; the compiler never imposes an unintended storage equality.
  const auto other = encodeSymbolicInitialState(macro, {origins[0], BoolExpr::Var(14)});
  EXPECT_EQ(other[0], origins[0]);
  EXPECT_EQ(other[1]->getSupportVars(), (std::set<size_t>{12, 14}));
  // Replacement IDs may collide with local IDs; substitution is simultaneous.
  const auto swapped = encodeSymbolicInitialState(macro, {BoolExpr::Var(6), BoolExpr::Var(5)});
  EXPECT_EQ(swapped[0], BoolExpr::Var(6));
  EXPECT_EQ(swapped[1], BoolExpr::Xor(BoolExpr::Var(6), BoolExpr::Var(5)));
}

TEST_F(LatchSymbolicIntegrationTests, InitialEncodingRejectsHiddenChoicesAndOrdinaryFrameInputs) {
  SymbolicMacro macro;
  macro.stateSymbols = {2}; macro.inputSymbols = {3};
  macro.initialParameterSymbols = {4};
  macro.initialStateExpressions = {BoolExpr::Var(4)};
  EXPECT_THROW(encodeSymbolicInitialState(macro, {}), std::invalid_argument);
  EXPECT_THROW(encodeSymbolicInitialState(macro, {nullptr}), std::invalid_argument);
  for (size_t illegal : {size_t(2), size_t(3), size_t(99)}) {
    macro.initialStateExpressions[0] = BoolExpr::Var(illegal);
    EXPECT_THROW(encodeSymbolicInitialState(macro, {BoolExpr::Var(20)}), std::invalid_argument);
  }
  macro.initialStateExpressions = {BoolExpr::Var(4)};
  macro.initialParameterSymbols = {2};
  EXPECT_THROW(encodeSymbolicInitialState(macro, {BoolExpr::Var(20)}), std::invalid_argument);
  macro.initialParameterSymbols = {4, 4};
  EXPECT_THROW(encodeSymbolicInitialState(macro, {BoolExpr::Var(20), BoolExpr::Var(21)}), std::invalid_argument);
  macro.initialParameterSymbols.clear(); macro.initialStateExpressions.clear();
  macro.initialState = {1};
  EXPECT_EQ(encodeSymbolicInitialState(macro, {}), (SymbolicBits{BoolExpr::createTrue()}));
}
}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
