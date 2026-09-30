// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <string>
#include <unordered_map>
#include <vector>

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
#include "latch/LatchConstantNet.h"
#include "latch/LatchSupportOptions.h"
#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using namespace naja::NL;

class LatchConstantNetTests : public ::testing::Test {
 protected:
  void SetUp() override {
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    designs_ = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  }
  void TearDown() override {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) universe->destroy();
    BoolExprCache::destroy();
  }
  SNLDesign* design(const std::string& name) {
    return SNLDesign::create(designs_, SNLDesign::Type::Standard, NLName(name));
  }
  SNLScalarNet* output(SNLDesign* top, const std::string& name, SNLNet::Type type) {
    auto* term = SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName(name));
    auto* net = SNLScalarNet::create(top, NLName(name));
    net->setType(type);
    term->setNet(net);
    return net;
  }
  SequentialDesignModel extract(SNLDesign* top) {
    SupportOptions options;
    options.enabled = true;
    options.singleInputChange = true;
    options.initialInputs = false;
    options.initialStorage = false;
    options.workers = 2;
    ScopedSupportOptions scope(options);
    return SequentialDesignModel::extract(top);
  }
  void expectConstant(const SequentialDesignModel& model, bool expected) {
    EXPECT_FALSE(model.hasUnsupportedFeatures());
    ASSERT_EQ(model.allObservedOutputs.size(), 1u);
    ASSERT_EQ(model.observedOutputs.size(), 1u);
    EXPECT_TRUE(model.skippedObservedOutputs.empty());
    std::unordered_map<size_t, bool> environment;
    for (const auto& key : model.stateBits)
      environment[model.inputVarByKey.at(key)] = model.initialStateValueByKey.at(key);
    for (const auto& key : model.environmentInputs)
      environment[model.inputVarByKey.at(key)] = false;
    EXPECT_EQ(model.observedOutputExprByKey.at(model.observedOutputs.front())->evaluate(environment), expected);
  }
  void expectOpaque(const SequentialDesignModel& model) {
    EXPECT_TRUE(model.observedOutputs.empty());
    ASSERT_EQ(model.allObservedOutputs.size(), 1u);
    ASSERT_EQ(model.skippedObservedOutputs.size(), 1u);
    const auto& reason = model.connectivitySkipInfoByKey.at(model.skippedObservedOutputs.front()).detail;
    EXPECT_FALSE(reason.empty());
  }
  SNLInstance* constantLatch(SNLDesign* top, SNLScalarNet* source, SNLScalarNet* target) {
    auto* primitive = NLDB0::getDLatch();
    auto* cell = SNLInstance::create(top, primitive, NLName("latch"));
    cell->getInstTerm(primitive->getScalarTerm(NLName("D")))->setNet(source);
    cell->getInstTerm(primitive->getScalarTerm(NLName("E")))->setNet(source);
    cell->getInstTerm(primitive->getScalarTerm(NLName("Q")))->setNet(target);
    return cell;
  }
  SNLDesign* constantChild(const std::string& name, SNLNet::Type type) {
    auto* child = design(name);
    auto* constant = output(child, "out", type);
    auto* held = SNLScalarNet::create(child, NLName("held"));
    // A real nested cell makes this a hierarchical design, while no primitive
    // drives the annotated output. The annotation also feeds its D/E pins.
    constantLatch(child, constant, held);
    return child;
  }
  SNLDesign* wrapper(const std::string& name, SNLDesign* child, SNLNet::Type type) {
    auto* top = design(name);
    auto* wire = output(top, "out", type);
    auto* cell = SNLInstance::create(top, child, NLName("nested"));
    cell->getInstTerm(child->getScalarTerm(NLName("out")))->setNet(wire);
    return top;
  }
  NLLibrary* designs_ = nullptr;
};

TEST_F(LatchConstantNetTests, DriverlessTopZeroAndOneAreObservedConstants) {
  for (bool value : {false, true}) {
    auto* top = design(value ? "one" : "zero");
    output(top, "out", value ? SNLNet::Type::Assign1 : SNLNet::Type::Assign0);
    expectConstant(extract(top), value);
  }
}

TEST_F(LatchConstantNetTests, DisconnectedOutputIsNotInventedAsZero) {
  auto* top = design("disconnected");
  SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  expectOpaque(extract(top));
}

TEST_F(LatchConstantNetTests, ConnectedButFloatingWireIsNotInventedAsZero) {
  auto* top = design("floating");
  output(top, "out", SNLNet::Type::Standard);
  expectOpaque(extract(top));
}

TEST_F(LatchConstantNetTests, DriverlessUnknownAndHighImpedanceRemainOpaque) {
  for (const auto type : {SNLNet::Type::AssignX, SNLNet::Type::AssignZ}) {
    auto* top = design(type == SNLNet::Type::AssignX ? "unknown" : "high_impedance");
    output(top, "out", type);
    expectOpaque(extract(top));
  }
}

TEST_F(LatchConstantNetTests, MatchingHierarchicalConstantsRemainKnown) {
  for (bool value : {false, true}) {
    const auto type = value ? SNLNet::Type::Assign1 : SNLNet::Type::Assign0;
    const std::string suffix = value ? "one" : "zero";
    auto* child = constantChild("child_" + suffix, type);
    auto* middle = wrapper("middle_" + suffix, child, type);
    auto* top = wrapper("top_" + suffix, middle, type);
    expectConstant(extract(top), value);
  }
}

TEST_F(LatchConstantNetTests, HierarchicalConstantFlowsThroughUnannotatedParentNets) {
  auto* child = constantChild("child", SNLNet::Type::Assign1);
  auto* middle = wrapper("middle", child, SNLNet::Type::Standard);
  auto* top = wrapper("top", middle, SNLNet::Type::Standard);
  expectConstant(extract(top), true);
}

TEST_F(LatchConstantNetTests, ConflictingDriverlessHierarchicalAnnotationsRemainOpaque) {
  for (bool parentValue : {false, true}) {
    const std::string suffix = parentValue ? "one" : "zero";
    auto* child = constantChild("child_" + suffix,
        parentValue ? SNLNet::Type::Assign0 : SNLNet::Type::Assign1);
    auto* middle = wrapper("middle_" + suffix, child, SNLNet::Type::Standard);
    auto* top = wrapper("top_" + suffix, middle,
        parentValue ? SNLNet::Type::Assign1 : SNLNet::Type::Assign0);
    expectOpaque(extract(top));
  }
}

TEST_F(LatchConstantNetTests, HierarchicalUnknownAnnotationCannotBeOverriddenByBooleanParent) {
  for (const auto type : {SNLNet::Type::AssignX, SNLNet::Type::AssignZ}) {
    const std::string suffix = type == SNLNet::Type::AssignX ? "x" : "z";
    auto* child = constantChild("child_" + suffix, type);
    auto* top = wrapper("top_" + suffix, child, SNLNet::Type::Assign0);
    expectOpaque(extract(top));
  }
}

TEST_F(LatchConstantNetTests, DriverlessConstantOnLatchPinsIsResolvedWithoutDriverCell) {
  auto* top = design("top");
  auto* one = SNLScalarNet::create(top, NLName("one"));
  one->setType(SNLNet::Type::Assign1);
  auto* out = output(top, "out", SNLNet::Type::Standard);
  constantLatch(top, one, out);
  // Explicit storage initialization is zero; transparent D=E=1 must settle to
  // one. Neither a dangling-input default nor the stored value is correct.
  expectConstant(extract(top), true);
}

TEST_F(LatchConstantNetTests, SharedHierarchyAndPrimitiveRemainUnchangedAcrossExtraction) {
  auto* child = constantChild("shared_child", SNLNet::Type::Assign1);
  auto* first = wrapper("first", child, SNLNet::Type::Standard);
  auto* second = wrapper("second", child, SNLNet::Type::Assign1);
  auto* childNet = child->getScalarNet(NLName("out"));
  auto* cell = child->getInstance(NLName("latch"));
  auto* primitive = cell->getModel();
  auto* data = primitive->getScalarTerm(NLName("D"));
  auto* enable = primitive->getScalarTerm(NLName("E"));
  auto* outputTerm = child->getScalarTerm(NLName("out"));
  NLUniverse::get()->setTopDesign(first);
  auto* originalDnl = naja::DNL::get();
  const auto dataOrder = data->getOrderID();
  const auto cellOrder = cell->getOrderID();
  const auto outputOrder = outputTerm->getOrderID();
  const auto childNetCount = child->getNets().size();
  const auto primitiveTermCount = primitive->getBitTerms().size();

  expectConstant(extract(second), true);
  expectConstant(extract(first), true);
  expectConstant(extract(second), true);

  EXPECT_EQ(naja::DNL::get(), originalDnl);
  EXPECT_EQ(NLUniverse::get()->getTopDesign(), first);
  EXPECT_EQ(cell->getModel(), primitive);
  EXPECT_EQ(child->getInstances().size(), 1u);
  EXPECT_EQ(child->getNets().size(), childNetCount);
  EXPECT_EQ(primitive->getBitTerms().size(), primitiveTermCount);
  EXPECT_EQ(childNet->getType(), SNLNet::Type::Assign1);
  EXPECT_EQ(outputTerm->getNet(), childNet);
  EXPECT_EQ(cell->getInstTerm(data)->getNet(), childNet);
  EXPECT_EQ(cell->getInstTerm(enable)->getNet(), childNet);
  EXPECT_EQ(data->getOrderID(), dataOrder);
  EXPECT_EQ(cell->getOrderID(), cellOrder);
  EXPECT_EQ(outputTerm->getOrderID(), outputOrder);
  EXPECT_EQ(SNLDesignModeling::getSequentialModel(primitive).kind,
      SNLDesignModeling::SequentialModel::Kind::Latch);
}

TEST_F(LatchConstantNetTests, ResolverTraversesAllHierarchicalAnnotationsWithoutChangingIsoIds) {
  auto* top = design("top");
  const std::vector<std::pair<SNLNet::Type, SNLNet::Type>> annotations{
      {SNLNet::Type::Assign0, SNLNet::Type::Assign0},
      {SNLNet::Type::Assign1, SNLNet::Type::Assign1},
      {SNLNet::Type::Assign0, SNLNet::Type::Assign1},
      {SNLNet::Type::Assign1, SNLNet::Type::Assign0},
      {SNLNet::Type::Assign0, SNLNet::Type::AssignX},
      {SNLNet::Type::Assign0, SNLNet::Type::AssignZ}};
  std::vector<SNLScalarTerm*> terms;
  for (size_t i = 0; i < annotations.size(); ++i) {
    const auto suffix = std::to_string(i);
    auto* child = constantChild("child" + suffix, annotations[i].second);
    auto* wire = output(top, "out" + suffix, annotations[i].first);
    auto* cell = SNLInstance::create(top, child, NLName("child" + suffix));
    cell->getInstTerm(child->getScalarTerm(NLName("out")))->setNet(wire);
    terms.push_back(top->getScalarTerm(NLName("out" + suffix)));
  }
  NLUniverse::get()->setTopDesign(top);
  const auto* dnl = naja::DNL::get();
  std::vector<naja::DNL::DNLID> originalIds;
  for (const auto& term : dnl->getDNLTerms()) originalIds.push_back(term.getIsoID());
  DriverlessConstantResolver resolver(*dnl);
  // Call the resolver directly even where DNL already supplies a valid iso:
  // these assertions exercise its traversal, not the adapter's normal path.
  for (size_t repetition = 0; repetition < 3; ++repetition) {
    for (size_t i = 0; i < terms.size(); ++i) {
      const auto& terminal = dnl->getTop().getTerminalFromBitTerm(terms[i]);
      if (i < 2) EXPECT_NE(terminal.getIsoID(), naja::DNL::DNLID_MAX);
      const auto result = resolver.resolve(terminal);
      if (i < 2) {
        ASSERT_TRUE(result.value.has_value());
        EXPECT_EQ(*result.value, i == 1);
        EXPECT_FALSE(result.conflictingOrUnknown);
      } else {
        EXPECT_FALSE(result.value.has_value());
        EXPECT_TRUE(result.conflictingOrUnknown);
      }
    }
  }
  ASSERT_EQ(dnl->getDNLTerms().size(), originalIds.size());
  for (size_t i = 0; i < originalIds.size(); ++i)
    EXPECT_EQ(dnl->getDNLTerms()[i].getIsoID(), originalIds[i]);
  EXPECT_EQ(naja::DNL::get(), dnl);
  EXPECT_EQ(NLUniverse::get()->getTopDesign(), top);
}

TEST_F(LatchConstantNetTests, ResolverRejectsConstantAnnotatedNetWithActualDriver) {
  auto* top = design("top");
  auto* wire = output(top, "out", SNLNet::Type::Assign0);
  auto* driver = SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("driver"));
  driver->setNet(wire);
  NLUniverse::get()->setTopDesign(top);
  const auto* dnl = naja::DNL::get();
  const auto& terminal = dnl->getTop().getTerminalFromBitTerm(top->getScalarTerm(NLName("out")));
  const auto originalId = terminal.getIsoID();
  ASSERT_NE(originalId, naja::DNL::DNLID_MAX);
  ASSERT_FALSE(dnl->getDNLIsoDB().getIsoFromIsoIDconst(originalId).getDrivers().empty());
  DriverlessConstantResolver resolver(*dnl);
  for (size_t repetition = 0; repetition < 3; ++repetition) {
    const auto result = resolver.resolve(terminal);
    EXPECT_FALSE(result.value.has_value());
    EXPECT_FALSE(result.conflictingOrUnknown);
    EXPECT_EQ(terminal.getIsoID(), originalId);
  }
  EXPECT_EQ(wire->getType(), SNLNet::Type::Assign0);
  EXPECT_EQ(driver->getNet(), wire);
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
