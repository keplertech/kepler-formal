// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <tuple>

#include "BoolExprCache.h"
#include "DesignBoundary.h"
#include "DNL.h"
#include "NLDB0.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "Tree2BoolExpr.h"
#include "model/SequentialDesignModel.h"
#include "strategy/SequentialEquivalenceStrategy.h"

namespace KEPLER_FORMAL::SEC {
namespace {

using namespace naja::NL;

struct DestroyBoundaryDnl {
  ~DestroyBoundaryDnl() { naja::DNL::destroy(); }
};

class BoundaryEquivalenceTests
    : public ::testing::TestWithParam<std::tuple<SecEngine, SecEncoding>> {
 protected:
  void SetUp() override {
    auto* universe = NLUniverse::create();
    auto* db = NLDB::create(universe);
    designs_ = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
    auto* primitives =
        NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("primitives"));
    opaque_ = SNLDesign::create(primitives, SNLDesign::Type::Primitive,
                                NLName("OPAQUE"));
    SNLScalarTerm::create(opaque_, SNLTerm::Direction::Input, NLName("A"));
    SNLScalarTerm::create(opaque_, SNLTerm::Direction::Output, NLName("Y"));
    inverter_ = SNLDesign::create(primitives, SNLDesign::Type::Primitive,
                                  NLName("INV"));
    auto* a = SNLScalarTerm::create(inverter_, SNLTerm::Direction::Input,
                                   NLName("A"));
    auto* y = SNLScalarTerm::create(inverter_, SNLTerm::Direction::Output,
                                   NLName("Y"));
    SNLDesignModeling::addCombinatorialArcs({a}, {y});
    SNLDesignModeling::setTruthTable(inverter_, SNLTruthTable::Inv());
  }

  void TearDown() override {
    naja::DNL::destroy();
    NLUniverse::get()->destroy();
    Tree2BoolExpr::iso2boolExpr_.clear();
    BoolExprCache::destroy();
  }

  SNLScalarNet* port(SNLDesign* top, const char* name,
                     SNLTerm::Direction direction) {
    auto* net = SNLScalarNet::create(top, NLName(name));
    SNLScalarTerm::create(top, direction, NLName(name))->setNet(net);
    return net;
  }

  SNLScalarNet* invert(SNLDesign* top, SNLScalarNet* input,
                       const char* name) {
    auto* output = SNLScalarNet::create(top, NLName(name));
    auto* instance = SNLInstance::create(top, inverter_, NLName(name));
    instance->getInstTerm(inverter_->getScalarTerm(NLName("A")))->setNet(input);
    instance->getInstTerm(inverter_->getScalarTerm(NLName("Y")))->setNet(output);
    return output;
  }

  SNLDesign* wrapper(const char* name, const char* instanceName,
                     bool invertInput = false, bool invertOutput = false,
                     bool registered = false, SNLDesign* blockModel = nullptr) {
    auto* top = SNLDesign::create(designs_, NLName(name));
    auto* a = port(top, "a", SNLTerm::Direction::Input);
    auto* clk = port(top, "clk", SNLTerm::Direction::Input);
    auto* rst = port(top, "rst", SNLTerm::Direction::Input);
    auto* y = SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"));
    auto* model = blockModel ? blockModel : opaque_;
    auto* block = SNLInstance::create(top, model, NLName(instanceName));
    auto* raw = SNLScalarNet::create(top, NLName("raw"));
    block->getInstTerm(model->getScalarTerm(NLName("A")))
        ->setNet(invertInput ? invert(top, a, "input_inv") : a);
    block->getInstTerm(model->getScalarTerm(NLName("Y")))->setNet(raw);
    auto* value = invertOutput ? invert(top, raw, "output_inv") : raw;
    if (registered) {
      auto* ff = SNLInstance::create(top, NLDB0::getDFFR(), NLName("ff"));
      auto* q = SNLScalarNet::create(top, NLName("q"));
      ff->getInstTerm(NLDB0::getDFFRClock())->setNet(clk);
      ff->getInstTerm(NLDB0::getDFFRReset())->setNet(rst);
      ff->getInstTerm(NLDB0::getDFFRData())->setNet(value);
      ff->getInstTerm(NLDB0::getDFFROutput())->setNet(q);
      value = q;
    }
    y->setNet(value);
    return top;
  }

  SequentialEquivalenceResult prove(SNLDesign* left, SNLDesign* right,
                                    BoundaryPairs pairs,
                                    bool registered = false) {
    BoundaryDesign first(left, pairs, 0);
    BoundaryDesign second(right, pairs, 1);
    DestroyBoundaryDnl cleanup;
    validateBoundaryInterfaces(first.getPorts(), second.getPorts());
    SecResetSpec reset;
    if (registered) {
      reset.cycles = 1;
      reset.ports.push_back({"rst", true});
    }
    return SequentialEquivalenceStrategy(
               first.getTop(), second.getTop(), Config::SolverType::KISSAT,
               std::get<0>(GetParam()), std::get<1>(GetParam()), reset)
        .run(4);
  }

  SNLDesign* sequentialBlock(const char* name, bool invertClock = false) {
    auto* top = SNLDesign::create(designs_, NLName(name));
    auto* a = port(top, "a", SNLTerm::Direction::Input);
    auto* clk = port(top, "clk", SNLTerm::Direction::Input);
    auto* rst = port(top, "rst", SNLTerm::Direction::Input);
    auto* y = port(top, "y", SNLTerm::Direction::Output);
    auto* ff = SNLInstance::create(top, NLDB0::getDFFR(), NLName("block"));
    ff->getInstTerm(NLDB0::getDFFRClock())
        ->setNet(invertClock ? invert(top, clk, "clock_inv") : clk);
    ff->getInstTerm(NLDB0::getDFFRReset())->setNet(rst);
    ff->getInstTerm(NLDB0::getDFFRData())->setNet(a);
    ff->getInstTerm(NLDB0::getDFFROutput())->setNet(y);
    return top;
  }

  SNLDesign* ordinaryRegisteredInterface(const char* name) {
    auto* top = wrapper(name, "macro", false, false, true);
    auto* block = top->getInstance(NLName("macro"));
    auto* data = block->getInstTerm(opaque_->getScalarTerm(NLName("Y")))->getNet();
    block->destroy();
    SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("external_data"))
        ->setNet(data);
    SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("checked_a"))
        ->setNet(top->getScalarTerm(NLName("a"))->getNet());
    return top;
  }

  SNLDesign* twoBlocks(const char* name, bool swapOutputs) {
    auto* top = SNLDesign::create(designs_, NLName(name));
    auto* a = port(top, "a", SNLTerm::Direction::Input);
    auto* y0 = port(top, "y0", SNLTerm::Direction::Output);
    auto* y1 = port(top, "y1", SNLTerm::Direction::Output);
    for (size_t index = 0; index < 2; ++index) {
      auto* block = SNLInstance::create(
          top, opaque_, NLName(index == 0 ? "block0" : "block1"));
      block->getInstTerm(opaque_->getScalarTerm(NLName("A")))->setNet(a);
      block->getInstTerm(opaque_->getScalarTerm(NLName("Y")))
          ->setNet((index == 0) != swapOutputs ? y0 : y1);
    }
    return top;
  }

  NLLibrary* designs_ = nullptr;
  SNLDesign* opaque_ = nullptr;
  SNLDesign* inverter_ = nullptr;
};

TEST_P(BoundaryEquivalenceTests, OpaquePairUsesNormalEquivalentStatus) {
  auto* left = wrapper("left", "macro");
  auto* right = wrapper("right", "renamed_macro");
  const auto result = prove(left, right, {{"macro", "renamed_macro"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.totalOutputs, 2u);
  EXPECT_EQ(result.coveredOutputs, 2u);
  EXPECT_TRUE(result.opaqueCellSkippedOutputs.empty());
  EXPECT_NE(left->getInstance(NLName("macro")), nullptr);
  EXPECT_EQ(left->getTerms().size(), 4u);
}

TEST_P(BoundaryEquivalenceTests, DifferentBoundaryInputsFailEvenWithEqualTopOutput) {
  const auto result = prove(wrapper("left", "macro"),
                            wrapper("right", "macro", true),
                            {{"macro", "macro"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
}

TEST_P(BoundaryEquivalenceTests, DifferentLogicAfterBoundaryFails) {
  const auto result = prove(wrapper("left", "macro"),
                            wrapper("right", "macro", false, true),
                            {{"macro", "macro"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
}

TEST_P(BoundaryEquivalenceTests, ConstantBoundaryInputRemainsAnObservedOutput) {
  auto* left = wrapper("left", "macro");
  auto* right = wrapper("right", "macro");
  auto* constant = SNLScalarNet::create(right, NLName("constant"));
  constant->setType(SNLNet::Type::Assign0);
  right->getInstance(NLName("macro"))
      ->getInstTerm(opaque_->getScalarTerm(NLName("A")))->setNet(constant);
  const auto result = prove(left, right, {{"macro", "macro"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
  EXPECT_EQ(result.totalOutputs, 2u);
}

TEST_P(BoundaryEquivalenceTests, IdenticalConstantBoundaryInputsAreProved) {
  auto* left = wrapper("left", "macro");
  auto* right = wrapper("right", "macro");
  for (auto* top : {left, right}) {
    auto* constant = SNLScalarNet::create(top, NLName("constant"));
    constant->setType(SNLNet::Type::Assign1);
    top->getInstance(NLName("macro"))
        ->getInstTerm(opaque_->getScalarTerm(NLName("A")))->setNet(constant);
  }
  const auto result = prove(left, right, {{"macro", "macro"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.totalOutputs, 2u);
  EXPECT_EQ(result.coveredOutputs, 2u);
}

TEST_P(BoundaryEquivalenceTests, DifferentConstantBoundaryInputsFail) {
  auto* left = wrapper("left", "macro");
  auto* right = wrapper("right", "macro");
  for (auto* top : {left, right}) {
    auto* constant = SNLScalarNet::create(top, NLName("constant"));
    constant->setType(top == left ? SNLNet::Type::Assign0 : SNLNet::Type::Assign1);
    top->getInstance(NLName("macro"))
        ->getInstTerm(opaque_->getScalarTerm(NLName("A")))->setNet(constant);
  }
  const auto result = prove(left, right, {{"macro", "macro"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
  EXPECT_EQ(result.totalOutputs, 2u);
}

TEST_P(BoundaryEquivalenceTests, ModelledAndOpaqueInternalsAreBothRemoved) {
  const auto result = prove(wrapper("left", "rtl", false, false, false, inverter_),
                            wrapper("right", "gate"), {{"rtl", "gate"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
}

TEST_P(BoundaryEquivalenceTests, SharedBoundaryOutputFeedsSequentialLogic) {
  const auto result = prove(wrapper("left", "macro", false, false, true),
                            wrapper("right", "macro", false, false, true),
                            {{"macro", "macro"}}, true);
  // Compare with an ordinary, manually exposed top interface. In particular,
  // current small-state IMC can be inconclusive with reset bootstrap; adding
  // boundaries must preserve that engine's normal semantics and coverage.
  auto* reference = ordinaryRegisteredInterface("external_left");
  auto* implementation = ordinaryRegisteredInterface("external_right");
  const auto baseline = SequentialEquivalenceStrategy(
                            reference, implementation, Config::SolverType::KISSAT,
                            std::get<0>(GetParam()), std::get<1>(GetParam()),
                            SecResetSpec{1, {{"rst", true}}})
                            .run(4);
  EXPECT_EQ(result.status, baseline.status) << result.reason;
  EXPECT_EQ(result.coveredOutputs, baseline.coveredOutputs);
  EXPECT_EQ(result.totalOutputs, 2u);
  if (std::get<0>(GetParam()) != SecEngine::Imc) {
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
    EXPECT_EQ(result.coveredOutputs, 2u);
  }
}

TEST_P(BoundaryEquivalenceTests, SequentialMutationAfterBoundaryFails) {
  const auto result = prove(wrapper("left", "macro", false, false, true),
                            wrapper("right", "macro", false, true, true),
                            {{"macro", "macro"}}, true);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
}

TEST_P(BoundaryEquivalenceTests, ExtractedModelsSurviveBoundaryDestruction) {
  auto* left = wrapper("left", "macro");
  auto* right = wrapper("right", "other");
  const BoundaryPairs pairs{{"macro", "other"}};
  SequentialDesignModel first;
  SequentialDesignModel second;
  std::vector<BoundaryPort> firstPorts;
  {
    BoundaryDesign boundary(left, pairs, 0);
    DestroyBoundaryDnl cleanup;
    firstPorts = boundary.getPorts();
    first = SequentialDesignModel::extract(boundary.getTop());
  }
  {
    BoundaryDesign boundary(right, pairs, 1);
    DestroyBoundaryDnl cleanup;
    validateBoundaryInterfaces(firstPorts, boundary.getPorts());
    second = SequentialDesignModel::extract(boundary.getTop());
  }
  const auto result = SequentialEquivalenceStrategy(
                          nullptr, nullptr, Config::SolverType::KISSAT,
                          std::get<0>(GetParam()), std::get<1>(GetParam()))
                          .runExtractedModels(first, second, 4);
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 2u);
}

TEST_P(BoundaryEquivalenceTests, StatefulBoundaryComparesClockResetAndData) {
  const auto result = prove(sequentialBlock("left"), sequentialBlock("right"),
                            {{"block", "block"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 4u);
}

TEST_P(BoundaryEquivalenceTests, StatefulBoundaryDetectsChangedClock) {
  const auto result = prove(sequentialBlock("left"),
                            sequentialBlock("right", true),
                            {{"block", "block"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
}

TEST_P(BoundaryEquivalenceTests, MultiplePairsExposeAllInputChecks) {
  const auto result = prove(twoBlocks("left", false), twoBlocks("right", false),
                            {{"block0", "block0"}, {"block1", "block1"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.coveredOutputs, 4u);
}

TEST_P(BoundaryEquivalenceTests, DifferentPairsHaveIndependentSharedInputs) {
  const auto result = prove(twoBlocks("left", false), twoBlocks("right", true),
                            {{"block0", "block0"}, {"block1", "block1"}});
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different) << result.reason;
}

INSTANTIATE_TEST_SUITE_P(
    AllEncodingsAndEngines, BoundaryEquivalenceTests,
    ::testing::Combine(
        ::testing::Values(SecEngine::KInduction, SecEngine::Imc, SecEngine::Pdr),
        ::testing::Values(SecEncoding::Binary, SecEncoding::DualRailSteady)));

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
