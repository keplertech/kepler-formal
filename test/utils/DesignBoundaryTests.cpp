// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "DesignBoundary.h"
#include "NLDB.h"
#include "NLDB0.h"
#include "NLLibrary.h"
#include "NLName.h"
#include "NLUniverse.h"
#include "SNLBundleTerm.h"
#include "SNLBusTerm.h"
#include "SNLBusTermBit.h"
#include "SNLDesign.h"
#include "SNLInstance.h"
#include "SNLInstTerm.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"

namespace KEPLER_FORMAL {
namespace {

using namespace naja::NL;

class DesignBoundaryTests : public ::testing::Test {
 protected:
  void SetUp() override {
    universe_ = NLUniverse::create();
    db_ = NLDB::create(universe_);
    designs_ = NLLibrary::create(
        db_, NLLibrary::Type::Standard, NLName("designs"));
    primitives_ = NLLibrary::create(
        db_, NLLibrary::Type::Primitives, NLName("primitives"));
  }

  void TearDown() override {
    if (NLUniverse::get() != nullptr) {
      NLUniverse::get()->destroy();
    }
  }

  SNLScalarNet* addPort(SNLDesign* design,
                        const std::string& name,
                        SNLTerm::Direction direction) {
    auto* net = SNLScalarNet::create(design, NLName(name + "_net"));
    SNLScalarTerm::create(design, direction, NLName(name))->setNet(net);
    return net;
  }

  SNLDesign* createScalarBlock(const std::string& name,
                               const std::vector<std::string>& inputs,
                               const std::vector<std::string>& outputs) {
    auto* model = SNLDesign::create(
        primitives_, SNLDesign::Type::Primitive, NLName(name));
    for (const auto& input : inputs) {
      SNLScalarTerm::create(
          model, SNLTerm::Direction::Input, NLName(input));
    }
    for (const auto& output : outputs) {
      SNLScalarTerm::create(
          model, SNLTerm::Direction::Output, NLName(output));
    }
    return model;
  }

  static const BoundaryPort* findPort(const BoundaryDesign& boundary,
                                      size_t pairIndex,
                                      const std::string& pinName,
                                      int32_t bit = 0) {
    const auto& ports = boundary.getPorts();
    const auto found = std::find_if(
        ports.begin(), ports.end(), [&](const BoundaryPort& port) {
          return port.pairIndex == pairIndex && port.pinName == pinName &&
                 port.bit == bit;
        });
    return found == ports.end() ? nullptr : &*found;
  }

  NLUniverse* universe_ = nullptr;
  NLDB* db_ = nullptr;
  NLLibrary* designs_ = nullptr;
  NLLibrary* primitives_ = nullptr;
};

TEST_F(DesignBoundaryTests, PromotesEveryPinBitAndLeavesSourceUntouched) {
  auto* block = SNLDesign::create(
      primitives_, SNLDesign::Type::Primitive, NLName("WIDE_BLOCK"));
  auto* input = SNLBusTerm::create(
      block, SNLTerm::Direction::Input, 3, 2, NLName("A"));
  auto* output = SNLBusTerm::create(
      block, SNLTerm::Direction::Output, 1, 0, NLName("Y"));

  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  for (int bit : {3, 2}) {
    auto* net = addPort(
        top, "a" + std::to_string(bit), SNLTerm::Direction::Input);
    instance->getInstTerm(input->getBit(bit))->setNet(net);
  }
  for (int bit : {1, 0}) {
    auto* net = addPort(
        top, "y" + std::to_string(bit), SNLTerm::Direction::Output);
    instance->getInstTerm(output->getBit(bit))->setNet(net);
  }
  db_->setTopDesign(top);
  const size_t sourceTermCount = top->getTerms().size();
  const size_t libraryCount = db_->getGlobalLibraries().size();

  {
    BoundaryDesign boundary(top, {{"u", "u"}}, 0);
    ASSERT_NE(boundary.getTop(), top);
    EXPECT_EQ(top, db_->getTopDesign());
    EXPECT_NE(nullptr, top->getInstance(NLName("u")));
    EXPECT_EQ(sourceTermCount, top->getTerms().size());
    EXPECT_EQ(nullptr, boundary.getTop()->getInstance(NLName("u")));
    ASSERT_EQ(4u, boundary.getPorts().size());

    for (int bit : {3, 2}) {
      const auto* port = findPort(boundary, 0, "A", bit);
      ASSERT_NE(nullptr, port);
      EXPECT_TRUE(port->isInput);
      EXPECT_EQ(2u, port->width);
      EXPECT_EQ(3, port->msb);
      EXPECT_EQ(2, port->lsb);
      auto* promoted =
          boundary.getTop()->getScalarTerm(NLName(port->topTermName));
      ASSERT_NE(nullptr, promoted);
      EXPECT_EQ(SNLTerm::Direction::Output, promoted->getDirection());
      EXPECT_NE(nullptr, promoted->getNet());
    }
    for (int bit : {1, 0}) {
      const auto* port = findPort(boundary, 0, "Y", bit);
      ASSERT_NE(nullptr, port);
      EXPECT_FALSE(port->isInput);
      EXPECT_EQ(2u, port->width);
      EXPECT_EQ(1, port->msb);
      EXPECT_EQ(0, port->lsb);
      auto* promoted =
          boundary.getTop()->getScalarTerm(NLName(port->topTermName));
      ASSERT_NE(nullptr, promoted);
      EXPECT_EQ(SNLTerm::Direction::Input, promoted->getDirection());
      EXPECT_NE(nullptr, promoted->getNet());
    }

    // If a client temporarily makes the scratch design active, RAII restores
    // the caller's top before releasing the scratch library.
    universe_->setTopDesign(boundary.getTop());
  }
  EXPECT_EQ(top, db_->getTopDesign());
  EXPECT_EQ(libraryCount, db_->getGlobalLibraries().size());
}

TEST_F(DesignBoundaryTests, EmptyBoundaryListClonesAnEmptyDesign) {
  auto* top = SNLDesign::create(designs_, NLName("empty_top"));

  BoundaryDesign boundary(top, {}, 0);
  ASSERT_NE(nullptr, boundary.getTop());
  EXPECT_NE(top, boundary.getTop());
  EXPECT_TRUE(boundary.getPorts().empty());
  EXPECT_EQ(0u, boundary.getTop()->getInstances().size());
  EXPECT_EQ(0u, boundary.getTop()->getTerms().size());
}

TEST_F(DesignBoundaryTests, HierarchicalCutUniquifiesOnlySelectedOccurrence) {
  auto* block = createScalarBlock("BLOCK", {"A"}, {"Y"});
  auto* wrapper = SNLDesign::create(designs_, NLName("wrapper"));
  auto* wrapperInput = SNLScalarTerm::create(
      wrapper, SNLTerm::Direction::Input, NLName("I"));
  auto* wrapperOutput = SNLScalarTerm::create(
      wrapper, SNLTerm::Direction::Output, NLName("O"));
  auto* inputNet = SNLScalarNet::create(wrapper, NLName("input"));
  auto* outputNet = SNLScalarNet::create(wrapper, NLName("output"));
  wrapperInput->setNet(inputNet);
  wrapperOutput->setNet(outputNet);
  auto* leaf = SNLInstance::create(wrapper, block, NLName("leaf"));
  leaf->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(inputNet);
  leaf->getInstTerm(block->getScalarTerm(NLName("Y")))->setNet(outputNet);

  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* topInput = addPort(top, "a", SNLTerm::Direction::Input);
  auto* topOutput = addPort(top, "y", SNLTerm::Direction::Output);
  auto* selected = SNLInstance::create(top, wrapper, NLName("selected"));
  auto* untouched = SNLInstance::create(top, wrapper, NLName("untouched"));
  selected->getInstTerm(wrapperInput)->setNet(topInput);
  selected->getInstTerm(wrapperOutput)->setNet(topOutput);

  BoundaryDesign boundary(top, {{"selected/leaf", "selected/leaf"}}, 0);
  auto* clonedSelected = boundary.getTop()->getInstance(NLName("selected"));
  auto* clonedUntouched = boundary.getTop()->getInstance(NLName("untouched"));
  ASSERT_NE(nullptr, clonedSelected);
  ASSERT_NE(nullptr, clonedUntouched);
  EXPECT_NE(wrapper, clonedSelected->getModel());
  EXPECT_EQ(wrapper, clonedUntouched->getModel());
  EXPECT_EQ(nullptr, clonedSelected->getModel()->getInstance(NLName("leaf")));
  EXPECT_NE(nullptr, wrapper->getInstance(NLName("leaf")));

  const auto* inputPort = findPort(boundary, 0, "A");
  const auto* outputPort = findPort(boundary, 0, "Y");
  ASSERT_NE(nullptr, inputPort);
  ASSERT_NE(nullptr, outputPort);
  auto* innerInput = clonedSelected->getModel()->getScalarTerm(
      NLName(inputPort->topTermName));
  auto* outerInput =
      boundary.getTop()->getScalarTerm(NLName(inputPort->topTermName));
  ASSERT_NE(nullptr, innerInput);
  ASSERT_NE(nullptr, outerInput);
  EXPECT_EQ(
      outerInput->getNet(), clonedSelected->getInstTerm(innerInput)->getNet());

  // The wrapper's pre-existing output alias remains connected while the new
  // boundary PI replaces the removed leaf's driver on the same inner net.
  auto* clonedWrapperOutput =
      clonedSelected->getModel()->getScalarTerm(NLName("O"));
  auto* innerOutput = clonedSelected->getModel()->getScalarTerm(
      NLName(outputPort->topTermName));
  ASSERT_NE(nullptr, clonedWrapperOutput);
  ASSERT_NE(nullptr, innerOutput);
  EXPECT_EQ(clonedWrapperOutput->getNet(), innerOutput->getNet());
  EXPECT_EQ(nullptr, wrapper->getTerm(NLName(inputPort->topTermName)));
  EXPECT_EQ(nullptr, wrapper->getTerm(NLName(outputPort->topTermName)));
  EXPECT_EQ(wrapper, untouched->getModel());
}

TEST_F(DesignBoundaryTests, AliasedInputsBecomeSeparateOutputsOnTheSameNet) {
  auto* block = createScalarBlock("ALIASED_INPUTS", {"A", "B"}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* sharedInput = addPort(top, "a", SNLTerm::Direction::Input);
  auto* output = addPort(top, "y", SNLTerm::Direction::Output);
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(sharedInput);
  instance->getInstTerm(block->getScalarTerm(NLName("B")))->setNet(sharedInput);
  instance->getInstTerm(block->getScalarTerm(NLName("Y")))->setNet(output);

  BoundaryDesign boundary(top, {{"u", "u"}}, 0);
  const auto* a = findPort(boundary, 0, "A");
  const auto* b = findPort(boundary, 0, "B");
  ASSERT_NE(nullptr, a);
  ASSERT_NE(nullptr, b);
  auto* promotedA = boundary.getTop()->getScalarTerm(NLName(a->topTermName));
  auto* promotedB = boundary.getTop()->getScalarTerm(NLName(b->topTermName));
  ASSERT_NE(nullptr, promotedA);
  ASSERT_NE(nullptr, promotedB);
  EXPECT_NE(promotedA, promotedB);
  EXPECT_EQ(promotedA->getNet(), promotedB->getNet());
}

TEST_F(DesignBoundaryTests, ConnectedBoundariesRetainTheirEquality) {
  auto* block = createScalarBlock("CHAIN_BLOCK", {"A"}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* input = addPort(top, "a", SNLTerm::Direction::Input);
  auto* output = addPort(top, "y", SNLTerm::Direction::Output);
  auto* middle = SNLScalarNet::create(top, NLName("middle"));
  auto* first = SNLInstance::create(top, block, NLName("first"));
  auto* second = SNLInstance::create(top, block, NLName("second"));
  first->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(input);
  first->getInstTerm(block->getScalarTerm(NLName("Y")))->setNet(middle);
  second->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(middle);
  second->getInstTerm(block->getScalarTerm(NLName("Y")))->setNet(output);

  BoundaryPairs pairs{{"first", "first"}, {"second", "second"}};
  BoundaryDesign boundary(top, pairs, 0);
  const auto* firstOutput = findPort(boundary, 0, "Y");
  const auto* secondInput = findPort(boundary, 1, "A");
  ASSERT_NE(nullptr, firstOutput);
  ASSERT_NE(nullptr, secondInput);
  auto* promotedFirstOutput = boundary.getTop()->getScalarTerm(
      NLName(firstOutput->topTermName));
  auto* promotedSecondInput = boundary.getTop()->getScalarTerm(
      NLName(secondInput->topTermName));
  ASSERT_NE(nullptr, promotedFirstOutput);
  ASSERT_NE(nullptr, promotedSecondInput);
  ASSERT_EQ(1u, boundary.getTop()->getInstances().size());
  auto* bridge = *boundary.getTop()->getInstances().begin();
  ASSERT_TRUE(NLDB0::isAssign(bridge->getModel()));
  EXPECT_EQ(
      promotedFirstOutput->getNet(),
      bridge->getInstTerm(NLDB0::getAssignInput())->getNet());
  EXPECT_EQ(
      promotedSecondInput->getNet(),
      bridge->getInstTerm(NLDB0::getAssignOutput())->getNet());
  EXPECT_EQ(nullptr, boundary.getTop()->getInstance(NLName("first")));
  EXPECT_EQ(nullptr, boundary.getTop()->getInstance(NLName("second")));
}

TEST_F(DesignBoundaryTests, UnconnectedOutputGetsFreshBoundaryInputNet) {
  auto* block = createScalarBlock("UNUSED_OUTPUT", {"A"}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* input = addPort(top, "a", SNLTerm::Direction::Input);
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(input);

  BoundaryDesign boundary(top, {{"u", "u"}}, 0);
  const auto* outputPort = findPort(boundary, 0, "Y");
  ASSERT_NE(nullptr, outputPort);
  auto* promoted =
      boundary.getTop()->getScalarTerm(NLName(outputPort->topTermName));
  ASSERT_NE(nullptr, promoted);
  EXPECT_EQ(SNLTerm::Direction::Input, promoted->getDirection());
  ASSERT_NE(nullptr, promoted->getNet());
  ASSERT_EQ(1u, boundary.getTop()->getInstances().size());
  auto* bridge = *boundary.getTop()->getInstances().begin();
  ASSERT_TRUE(NLDB0::isAssign(bridge->getModel()));
  EXPECT_EQ(
      promoted->getNet(),
      bridge->getInstTerm(NLDB0::getAssignInput())->getNet());
  ASSERT_NE(nullptr, bridge->getInstTerm(NLDB0::getAssignOutput())->getNet());
  EXPECT_NE(
      promoted->getNet(),
      bridge->getInstTerm(NLDB0::getAssignOutput())->getNet());
}

TEST_F(DesignBoundaryTests, BuffersAConstantCheckpointOnAnInstanceFreeTop) {
  auto* block = createScalarBlock("INPUT_ONLY", {"A"}, {});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* constant = SNLScalarNet::create(top, NLName("constant"));
  constant->setType(SNLNet::Type::Assign0);
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(constant);

  BoundaryDesign boundary(top, {{"u", "u"}}, 0);
  const auto* inputPort = findPort(boundary, 0, "A");
  ASSERT_NE(nullptr, inputPort);
  auto* promoted =
      boundary.getTop()->getScalarTerm(NLName(inputPort->topTermName));
  ASSERT_NE(nullptr, promoted);
  ASSERT_NE(nullptr, promoted->getNet());
  EXPECT_FALSE(promoted->getNet()->isConstant());

  ASSERT_EQ(1u, boundary.getTop()->getInstances().size());
  auto* bridge = *boundary.getTop()->getInstances().begin();
  ASSERT_TRUE(NLDB0::isAssign(bridge->getModel()));
  auto* bridgeInput = bridge->getInstTerm(NLDB0::getAssignInput());
  auto* bridgeOutput = bridge->getInstTerm(NLDB0::getAssignOutput());
  ASSERT_NE(nullptr, bridgeInput->getNet());
  EXPECT_TRUE(bridgeInput->getNet()->isConstant());
  EXPECT_EQ(promoted->getNet(), bridgeOutput->getNet());

  // The source remains unchanged, including its direct constant connection.
  EXPECT_EQ(constant, instance->getInstTerm(
                          block->getScalarTerm(NLName("A")))->getNet());
}

TEST_F(DesignBoundaryTests, RejectsAliasedOrMultiplyDrivenOutputs) {
  auto* block = createScalarBlock("ALIASED_OUTPUTS", {"A"}, {"Y", "Z"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* input = addPort(top, "a", SNLTerm::Direction::Input);
  auto* output = addPort(top, "y", SNLTerm::Direction::Output);
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(input);
  instance->getInstTerm(block->getScalarTerm(NLName("Y")))->setNet(output);
  instance->getInstTerm(block->getScalarTerm(NLName("Z")))->setNet(output);

  EXPECT_THROW(
      BoundaryDesign(top, {{"u", "u"}}, 0), std::invalid_argument);
  EXPECT_NE(nullptr, top->getInstance(NLName("u")));
}

TEST_F(DesignBoundaryTests, RejectsInvalidDuplicateAndNestedPaths) {
  auto* block = createScalarBlock("BLOCK", {"A"}, {"Y"});
  auto* wrapper = SNLDesign::create(designs_, NLName("wrapper"));
  SNLInstance::create(wrapper, block, NLName("leaf"));
  auto* top = SNLDesign::create(designs_, NLName("top"));
  SNLInstance::create(top, wrapper, NLName("wrap"));

  EXPECT_THROW(
      BoundaryDesign(top, {{"missing", "missing"}}, 0),
      std::invalid_argument);
  EXPECT_THROW(
      BoundaryDesign(
          top, {{"wrap/leaf", "wrap/leaf"}, {"wrap/leaf", "wrap/leaf"}}, 0),
      std::invalid_argument);
  EXPECT_THROW(
      BoundaryDesign(
          top, {{"wrap", "wrap"}, {"wrap/leaf", "wrap/leaf"}}, 0),
      std::invalid_argument);
}

TEST_F(DesignBoundaryTests, RejectsUnconnectedInputAndPinlessInstance) {
  auto* block = createScalarBlock("BLOCK", {"A"}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  SNLInstance::create(top, block, NLName("u"));
  EXPECT_THROW(
      BoundaryDesign(top, {{"u", "u"}}, 0), std::invalid_argument);

  auto* empty = createScalarBlock("EMPTY", {}, {});
  auto* emptyTop = SNLDesign::create(designs_, NLName("empty_top"));
  SNLInstance::create(emptyTop, empty, NLName("u"));
  EXPECT_THROW(
      BoundaryDesign(emptyTop, {{"u", "u"}}, 0), std::invalid_argument);
}

TEST_F(DesignBoundaryTests, RejectsBoundaryInputWithoutADriver) {
  auto* block = createScalarBlock("INPUT_ONLY", {"A"}, {});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* undriven = SNLScalarNet::create(top, NLName("undriven"));
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(undriven);

  EXPECT_THROW(
      BoundaryDesign(top, {{"u", "u"}}, 0), std::invalid_argument);
  EXPECT_NE(nullptr, top->getInstance(NLName("u")));
}

TEST_F(DesignBoundaryTests, RejectsMultiplyDrivenBoundaryInput) {
  auto* block = createScalarBlock("INPUT_ONLY", {"A"}, {});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* multiplyDriven = SNLScalarNet::create(top, NLName("multiply_driven"));
  SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("i0"))
      ->setNet(multiplyDriven);
  SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("i1"))
      ->setNet(multiplyDriven);
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("A")))
      ->setNet(multiplyDriven);

  EXPECT_THROW(
      BoundaryDesign(top, {{"u", "u"}}, 0), std::invalid_argument);
  EXPECT_NE(nullptr, top->getInstance(NLName("u")));
}

TEST_F(DesignBoundaryTests, ValidatesInterfacesByPairPinAndShape) {
  BoundaryPort input;
  input.pairIndex = 0;
  input.pinName = "A";
  input.bit = 3;
  input.isInput = true;
  input.width = 2;
  input.msb = 3;
  input.lsb = 2;
  input.topTermName = "boundary_a3";

  BoundaryPort output;
  output.pairIndex = 0;
  output.pinName = "Y";
  output.bit = 0;
  output.isInput = false;
  output.topTermName = "boundary_y0";

  EXPECT_NO_THROW(validateBoundaryInterfaces({input, output}, {output, input}));

  auto wrongDirection = input;
  wrongDirection.isInput = false;
  EXPECT_THROW(
      validateBoundaryInterfaces({input}, {wrongDirection}),
      std::invalid_argument);

  auto wrongShape = input;
  wrongShape.msb = 4;
  EXPECT_THROW(
      validateBoundaryInterfaces({input}, {wrongShape}),
      std::invalid_argument);

  EXPECT_THROW(
      validateBoundaryInterfaces({input, output}, {input}),
      std::invalid_argument);
  EXPECT_THROW(
      validateBoundaryInterfaces({input, input}, {input, input}),
      std::invalid_argument);
}

// Check the diagnostic as well as the exception type so a malformed fixture
// cannot accidentally exercise a different validation failure.
void expectInvalidBoundary(SNLDesign* top,
                           const BoundaryPairs& pairs,
                           size_t side,
                           const std::string& diagnostic) {
  try {
    BoundaryDesign boundary(top, pairs, side);
    FAIL() << "Expected boundary validation failure: " << diagnostic;
  } catch (const std::invalid_argument& error) {
    EXPECT_NE(std::string::npos, std::string(error.what()).find(diagnostic))
        << error.what();
  }
}

TEST_F(DesignBoundaryTests, RejectsInvalidTopSideAndIncompletePairs) {
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* primitive = createScalarBlock("BLOCK", {}, {"Y"});
  const size_t libraryCount = db_->getGlobalLibraries().size();

  expectInvalidBoundary(nullptr, {}, 0, "top design must not be null");
  expectInvalidBoundary(top, {}, 2, "side must be 0 or 1");
  expectInvalidBoundary(primitive, {}, 0, "primitive design");
  expectInvalidBoundary(top, {{"", "u"}}, 0, "path on both sides");
  expectInvalidBoundary(top, {{"u", ""}}, 1, "path on both sides");
  // The other side is also validated even when it is not the selected side.
  expectInvalidBoundary(top, {{"u", ""}}, 0, "path on both sides");
  EXPECT_EQ(libraryCount, db_->getGlobalLibraries().size());
}

TEST_F(DesignBoundaryTests, RejectsMalformedPathsAndReportsMissingAncestor) {
  auto* block = createScalarBlock("BLOCK", {}, {"Y"});
  auto* wrapper = SNLDesign::create(designs_, NLName("wrapper"));
  SNLInstance::create(wrapper, block, NLName("leaf"));
  auto* top = SNLDesign::create(designs_, NLName("top"));
  SNLInstance::create(top, wrapper, NLName("wrap"));

  for (const std::string path : {"/wrap/leaf", "wrap/leaf/"}) {
    SCOPED_TRACE(path);
    expectInvalidBoundary(top, {{path, path}}, 0, "leading or trailing slash");
  }
  expectInvalidBoundary(
      top, {{"wrap//leaf", "wrap//leaf"}}, 0, "empty component");
  expectInvalidBoundary(
      top, {{"wrap/missing", "wrap/missing"}}, 0,
      "does not resolve at `wrap/missing`");
  expectInvalidBoundary(
      top, {{"wrap/leaf", "wrap/leaf"}, {"wrap", "wrap"}}, 0,
      "nested boundary instance paths");
  EXPECT_NE(nullptr, wrapper->getInstance(NLName("leaf")));
}

TEST_F(DesignBoundaryTests, PromotesNegativeBusIndicesWithoutNameCollisions) {
  auto* block = SNLDesign::create(
      primitives_, SNLDesign::Type::Primitive, NLName("NEGATIVE_BUS"));
  SNLBusTerm::create(block, SNLTerm::Direction::Output, -1, 1, NLName("Y"));
  auto* top = SNLDesign::create(designs_, NLName("top"));
  SNLInstance::create(top, block, NLName("u"));

  BoundaryDesign boundary(top, {{"unused_on_left", "u"}}, 1);
  ASSERT_EQ(3u, boundary.getPorts().size());
  for (int bit : {-1, 0, 1}) {
    const auto* port = findPort(boundary, 0, "Y", bit);
    ASSERT_NE(nullptr, port);
    EXPECT_EQ(3u, port->width);
    EXPECT_EQ(-1, port->msb);
    EXPECT_EQ(1, port->lsb);
    EXPECT_FALSE(port->isInput);
    auto* promoted = boundary.getTop()->getScalarTerm(NLName(port->topTermName));
    ASSERT_NE(nullptr, promoted);
    EXPECT_EQ(SNLTerm::Direction::Input, promoted->getDirection());
    EXPECT_NE(nullptr, promoted->getNet());
  }
  EXPECT_EQ("__kepler_boundary_p0_n59_bn1",
            findPort(boundary, 0, "Y", -1)->topTermName);
  EXPECT_EQ("__kepler_boundary_p0_n59_bp1",
            findPort(boundary, 0, "Y", 1)->topTermName);
  EXPECT_NE(nullptr, top->getInstance(NLName("u")));
}

TEST_F(DesignBoundaryTests, RejectsConstantDrivenBoundaryOutput) {
  auto* block = createScalarBlock("BLOCK", {}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* constant = SNLScalarNet::create(top);
  constant->setType(SNLNet::Type::Assign1);
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("Y")))->setNet(constant);

  expectInvalidBoundary(top, {{"u", "u"}}, 0, "connected to a constant net");
  EXPECT_EQ(constant, instance->getInstTerm(
                          block->getScalarTerm(NLName("Y")))->getNet());
}

TEST_F(DesignBoundaryTests, RejectsUnnamedAndUnsupportedDirectionPins) {
  auto* unnamed = SNLDesign::create(
      primitives_, SNLDesign::Type::Primitive, NLName("UNNAMED"));
  SNLScalarTerm::create(unnamed, SNLTerm::Direction::Output);
  auto* unnamedTop = SNLDesign::create(designs_, NLName("unnamed_top"));
  SNLInstance::create(unnamedTop, unnamed, NLName("u"));
  expectInvalidBoundary(unnamedTop, {{"u", "u"}}, 0, "unnamed pins");

  for (const SNLTerm::Direction direction : {SNLTerm::Direction::InOut,
                                            SNLTerm::Direction::Undefined}) {
    SCOPED_TRACE(direction.getString());
    auto* block = SNLDesign::create(primitives_, SNLDesign::Type::Primitive);
    SNLScalarTerm::create(block, direction, NLName("P"));
    auto* top = SNLDesign::create(designs_);
    SNLInstance::create(top, block, NLName("u"));
    expectInvalidBoundary(top, {{"u", "u"}}, 0, "unsupported direction");
  }
}

TEST_F(DesignBoundaryTests, RejectsScalarAndBusBundleMembers) {
  for (const bool busMember : {false, true}) {
    SCOPED_TRACE(busMember);
    auto* block = SNLDesign::create(primitives_, SNLDesign::Type::Primitive);
    auto* bundle = SNLBundleTerm::create(
        block, SNLTerm::Direction::Output, NLName("BUNDLE"));
    if (busMember) {
      SNLBusTerm::create(bundle, SNLTerm::Direction::Output, 1, 0, NLName("Y"));
    } else {
      SNLScalarTerm::create(bundle, SNLTerm::Direction::Output, NLName("Y"));
    }
    auto* top = SNLDesign::create(designs_);
    SNLInstance::create(top, block, NLName("u"));
    expectInvalidBoundary(top, {{"u", "u"}}, 0, "bundled boundary pins");
    EXPECT_NE(nullptr, top->getInstance(NLName("u")));
  }
}

TEST_F(DesignBoundaryTests, RejectsTopAndAncestorSyntheticTermCollisions) {
  auto* block = createScalarBlock("BLOCK", {}, {"Y"});
  const NLName syntheticName("__kepler_boundary_p0_n59_bp0");
  auto* directTop = SNLDesign::create(designs_, NLName("direct_top"));
  SNLInstance::create(directTop, block, NLName("u"));
  SNLScalarTerm::create(directTop, SNLTerm::Direction::Input, syntheticName);
  const size_t libraryCount = db_->getGlobalLibraries().size();
  expectInvalidBoundary(
      directTop, {{"u", "u"}}, 0, "top term collides with existing term");

  auto* wrapper = SNLDesign::create(designs_, NLName("wrapper"));
  SNLInstance::create(wrapper, block, NLName("leaf"));
  SNLScalarTerm::create(wrapper, SNLTerm::Direction::Input, syntheticName);
  auto* nestedTop = SNLDesign::create(designs_, NLName("nested_top"));
  SNLInstance::create(nestedTop, wrapper, NLName("wrap"));
  expectInvalidBoundary(
      nestedTop, {{"wrap/leaf", "wrap/leaf"}}, 0, "collides in ancestor model");

  EXPECT_EQ(libraryCount, db_->getGlobalLibraries().size());
  EXPECT_NE(nullptr, directTop->getInstance(NLName("u")));
  EXPECT_NE(nullptr, wrapper->getInstance(NLName("leaf")));
}

TEST_F(DesignBoundaryTests, MovesTransferScratchOwnershipAndRejectStaleAccess) {
  auto* block = createScalarBlock("BLOCK", {}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  SNLInstance::create(top, block, NLName("u"));
  const size_t libraryCount = db_->getGlobalLibraries().size();
  {
    BoundaryDesign source(top, {{"u", "u"}}, 0);
    auto* scratchTop = source.getTop();
    const auto portName = source.getPorts().front().topTermName;
    BoundaryDesign moved(std::move(source));
    EXPECT_EQ(nullptr, source.getTop());
    EXPECT_THROW(source.getPorts(), std::logic_error);
    EXPECT_EQ(scratchTop, moved.getTop());
    EXPECT_EQ(portName, moved.getPorts().front().topTermName);

    BoundaryDesign destination(top, {{"u", "u"}}, 0);
    EXPECT_EQ(libraryCount + 2, db_->getGlobalLibraries().size());
    destination = std::move(moved);
    EXPECT_EQ(nullptr, moved.getTop());
    EXPECT_THROW(moved.getPorts(), std::logic_error);
    EXPECT_EQ(scratchTop, destination.getTop());
    EXPECT_EQ(portName, destination.getPorts().front().topTermName);
    EXPECT_EQ(libraryCount + 1, db_->getGlobalLibraries().size());
  }
  EXPECT_EQ(libraryCount, db_->getGlobalLibraries().size());
  EXPECT_NE(nullptr, top->getInstance(NLName("u")));
}

TEST_F(DesignBoundaryTests, RestoresPreviousUniverseTopAcrossDatabases) {
  auto* otherDB = NLDB::create(universe_);
  auto* otherLibrary = NLLibrary::create(otherDB, NLName("other"));
  auto* otherTop = SNLDesign::create(otherLibrary, NLName("other_top"));
  universe_->setTopDesign(otherTop);
  auto* top = SNLDesign::create(designs_, NLName("top"));
  db_->setTopDesign(top);

  {
    BoundaryDesign boundary(top, {}, 0);
    universe_->setTopDesign(boundary.getTop());
    ASSERT_EQ(db_, universe_->getTopDB());
  }
  EXPECT_EQ(top, db_->getTopDesign());
  EXPECT_EQ(otherDB, universe_->getTopDB());
  EXPECT_EQ(otherTop, universe_->getTopDesign());
}

TEST_F(DesignBoundaryTests, RejectsMissingKeysSyntheticNamesAndRightDuplicates) {
  BoundaryPort port;
  port.pairIndex = 1;
  port.pinName = "A";
  port.bit = -1;
  port.isInput = true;
  port.width = 2;
  port.msb = 0;
  port.lsb = -1;
  port.topTermName = "boundary_a_minus1";

  auto wrongPair = port;
  wrongPair.pairIndex = 2;
  EXPECT_THROW(validateBoundaryInterfaces({port}, {wrongPair}),
               std::invalid_argument);
  auto wrongPin = port;
  wrongPin.pinName = "B";
  EXPECT_THROW(validateBoundaryInterfaces({port}, {wrongPin}),
               std::invalid_argument);
  auto wrongBit = port;
  wrongBit.bit = 0;
  EXPECT_THROW(validateBoundaryInterfaces({port}, {wrongBit}),
               std::invalid_argument);
  auto wrongName = port;
  wrongName.topTermName = "different_boundary_name";
  EXPECT_THROW(validateBoundaryInterfaces({port}, {wrongName}),
               std::invalid_argument);
  auto wrongWidth = port;
  wrongWidth.width = 3;
  EXPECT_THROW(validateBoundaryInterfaces({port}, {wrongWidth}),
               std::invalid_argument);
  auto wrongLSB = port;
  wrongLSB.lsb = -2;
  EXPECT_THROW(validateBoundaryInterfaces({port}, {wrongLSB}),
               std::invalid_argument);
  EXPECT_THROW(validateBoundaryInterfaces({port}, {port, port}),
               std::invalid_argument);
  EXPECT_NO_THROW(validateBoundaryInterfaces({}, {}));
}

}  // namespace
}  // namespace KEPLER_FORMAL
