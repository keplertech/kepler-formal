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
    naja::DNL::destroy();
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

  static const BoundaryPort* findPort(const BoundarySelection& boundary,
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

TEST_F(DesignBoundaryTests, SelectsBusBitsWithoutChangingAnyDesignObjects) {
  auto* block = SNLDesign::create(primitives_, SNLDesign::Type::Primitive);
  auto* input = SNLBusTerm::create(block, SNLTerm::Direction::Input, 0, -1, NLName("A"));
  auto* output = SNLBusTerm::create(block, SNLTerm::Direction::Output, 1, 0, NLName("Y"));
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  for (int bit : {0, -1}) {
    instance->getInstTerm(input->getBit(bit))->setNet(
        addPort(top, "a" + std::to_string(bit), SNLTerm::Direction::Input));
  }
  for (int bit : {1, 0}) {
    instance->getInstTerm(output->getBit(bit))->setNet(
        addPort(top, "y" + std::to_string(bit), SNLTerm::Direction::Output));
  }
  universe_->setTopDesign(top);
  const auto libraries = db_->getGlobalLibraries().size();
  const auto nets = top->getNets().size();
  const auto terms = top->getTerms().size();
  BoundarySelection selection(top, {{"u", "u"}}, 0);
  ASSERT_EQ(4u, selection.getPorts().size());
  const auto* port = findPort(selection, 0, "A", -1);
  ASSERT_NE(nullptr, port);
  EXPECT_EQ(2u, port->width);
  EXPECT_EQ(0, port->msb);
  EXPECT_EQ(-1, port->lsb);
  const auto* dnl = naja::DNL::get();
  const auto dnlTerms = dnl->getNBterms();
  LeafBoundary boundary(*dnl, {{"u", "u"}}, 0);
  EXPECT_EQ(2u, boundary.getInputs().size());
  EXPECT_EQ(2u, boundary.getOutputs().size());
  EXPECT_EQ(dnlTerms, dnl->getNBterms());
  EXPECT_EQ(instance, top->getInstance(NLName("u")));
  EXPECT_EQ(block, instance->getModel());
  EXPECT_EQ(libraries, db_->getGlobalLibraries().size());
  EXPECT_EQ(nets, top->getNets().size());
  EXPECT_EQ(terms, top->getTerms().size());
  EXPECT_EQ(top, universe_->getTopDesign());
}

TEST_F(DesignBoundaryTests, RejectsLeafOutputInternalWireAliasOrConstant) {
  for (bool constant : {false, true}) {
    SCOPED_TRACE(constant);
    auto* block = SNLDesign::create(designs_);
    auto* a = SNLScalarTerm::create(block, SNLTerm::Direction::Input, NLName("A"));
    auto* y = SNLScalarTerm::create(block, SNLTerm::Direction::Output, NLName("Y"));
    auto* inner = SNLScalarNet::create(block);
    y->setNet(inner);
    if (constant) inner->setType(SNLNet::Type::Assign1);
    else a->setNet(inner);
    auto* top = SNLDesign::create(designs_);
    auto* input = addPort(top, "a", SNLTerm::Direction::Input);
    auto* output = addPort(top, "y", SNLTerm::Direction::Output);
    auto* instance = SNLInstance::create(top, block, NLName("u"));
    instance->getInstTerm(a)->setNet(input);
    instance->getInstTerm(y)->setNet(output);
    universe_->setTopDesign(top);
    const auto* dnl = naja::DNL::get();
    const auto& occurrence = dnl->getTop().getChildInstance(instance);
    const auto outputID = occurrence.getTerminalFromBitTerm(y).getID();
    const auto originalIso = dnl->getDNLTerminalFromID(outputID).getIsoID();
    EXPECT_THROW(LeafBoundary(*dnl, {{"u", "u"}}, 0), std::invalid_argument);
    EXPECT_EQ(originalIso, dnl->getDNLTerminalFromID(outputID).getIsoID());
    EXPECT_EQ(inner, y->getNet());
    EXPECT_EQ(input, instance->getInstTerm(a)->getNet());
    EXPECT_EQ(output, instance->getInstTerm(y)->getNet());
    naja::DNL::destroy();
  }
}

TEST_F(DesignBoundaryTests, NestedOccurrenceDoesNotAffectSiblingUsingSameModel) {
  auto* block = createScalarBlock("BLOCK", {"A"}, {"Y"});
  auto* wrapper = SNLDesign::create(designs_, NLName("wrapper"));
  auto* a = addPort(wrapper, "A", SNLTerm::Direction::Input);
  auto* y = addPort(wrapper, "Y", SNLTerm::Direction::Output);
  auto* leaf = SNLInstance::create(wrapper, block, NLName("leaf"));
  leaf->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(a);
  leaf->getInstTerm(block->getScalarTerm(NLName("Y")))->setNet(y);
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* input = addPort(top, "a", SNLTerm::Direction::Input);
  for (const std::string name : {"first", "second"}) {
    auto* instance = SNLInstance::create(top, wrapper, NLName(name));
    instance->getInstTerm(wrapper->getScalarTerm(NLName("A")))->setNet(input);
    instance->getInstTerm(wrapper->getScalarTerm(NLName("Y")))->setNet(
        addPort(top, name, SNLTerm::Direction::Output));
  }
  universe_->setTopDesign(top);
  const auto* dnl = naja::DNL::get();
  const auto& first = dnl->getTop().getChildInstance(top->getInstance(NLName("first")));
  const auto& second = dnl->getTop().getChildInstance(top->getInstance(NLName("second")));
  LeafBoundary boundary(*dnl, {{"first/leaf", "first/leaf"}}, 0);
  EXPECT_FALSE(boundary.containsInstance(first.getID()));
  EXPECT_TRUE(boundary.containsInstance(first.getChildInstance(leaf).getID()));
  EXPECT_FALSE(boundary.containsInstance(second.getChildInstance(leaf).getID()));
  EXPECT_EQ(wrapper, top->getInstance(NLName("first"))->getModel());
  EXPECT_EQ(wrapper, top->getInstance(NLName("second"))->getModel());
  EXPECT_EQ(leaf, wrapper->getInstance(NLName("leaf")));
}

TEST_F(DesignBoundaryTests, RejectsNonLeafBoundaryWithoutChangingHierarchy) {
  auto* block = createScalarBlock("BLOCK", {"A"}, {"Y"});
  auto* wrapper = SNLDesign::create(designs_, NLName("wrapper"));
  auto* a = addPort(wrapper, "A", SNLTerm::Direction::Input);
  auto* y = addPort(wrapper, "Y", SNLTerm::Direction::Output);
  auto* leaf = SNLInstance::create(wrapper, block, NLName("leaf"));
  leaf->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(a);
  leaf->getInstTerm(block->getScalarTerm(NLName("Y")))->setNet(y);
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* instance = SNLInstance::create(top, wrapper, NLName("u"));
  instance->getInstTerm(wrapper->getScalarTerm(NLName("A")))->setNet(
      addPort(top, "a", SNLTerm::Direction::Input));
  instance->getInstTerm(wrapper->getScalarTerm(NLName("Y")))->setNet(
      addPort(top, "y", SNLTerm::Direction::Output));
  universe_->setTopDesign(top);
  try {
    BoundarySelection(top, {{"u", "u"}}, 0);
    FAIL() << "non-leaf selection must be rejected";
  } catch (const std::invalid_argument& error) {
    EXPECT_NE(std::string::npos, std::string(error.what()).find("not a leaf"));
  }
  EXPECT_EQ(wrapper, instance->getModel());
  EXPECT_EQ(leaf, wrapper->getInstance(NLName("leaf")));
  EXPECT_THROW(LeafBoundary(*naja::DNL::get(), {{"u", "u"}}, 0),
               std::invalid_argument);
}

TEST_F(DesignBoundaryTests, ConstantInputAndUnusedOutputNeedNoNewNetsOrInstances) {
  auto* block = createScalarBlock("BLOCK", {"A"}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* constant = SNLScalarNet::create(top);
  constant->setType(SNLNet::Type::Assign0);
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(constant);
  universe_->setTopDesign(top);
  const auto* dnl = naja::DNL::get();
  LeafBoundary boundary(*dnl, {{"u", "u"}}, 0);
  ASSERT_EQ(1u, boundary.getOutputs().size());
  const auto isoID = dnl->getDNLTerminalFromID(boundary.getOutputs().front()).getIsoID();
  EXPECT_TRUE(dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID).isConstant0());
  ASSERT_EQ(1u, boundary.getInputs().size());
  EXPECT_TRUE(boundary.isInput(boundary.getInputs().front()));
  EXPECT_TRUE(boundary.getPort(boundary.getOutputs().front())->isInput);
  EXPECT_EQ(nullptr, instance->getInstTerm(block->getScalarTerm(NLName("Y")))->getNet());
  EXPECT_EQ(1u, top->getInstances().size());
  EXPECT_EQ(1u, top->getNets().size());
  EXPECT_EQ(0u, top->getTerms().size());
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
      BoundarySelection(top, {{"u", "u"}}, 0), std::invalid_argument);
  EXPECT_NE(nullptr, top->getInstance(NLName("u")));
}

TEST_F(DesignBoundaryTests, RejectsInvalidDuplicateAndNestedPaths) {
  auto* block = createScalarBlock("BLOCK", {"A"}, {"Y"});
  auto* wrapper = SNLDesign::create(designs_, NLName("wrapper"));
  SNLInstance::create(wrapper, block, NLName("leaf"));
  auto* top = SNLDesign::create(designs_, NLName("top"));
  SNLInstance::create(top, wrapper, NLName("wrap"));

  EXPECT_THROW(
      BoundarySelection(top, {{"missing", "missing"}}, 0),
      std::invalid_argument);
  EXPECT_THROW(
      BoundarySelection(
          top, {{"wrap/leaf", "wrap/leaf"}, {"wrap/leaf", "wrap/leaf"}}, 0),
      std::invalid_argument);
  EXPECT_THROW(
      BoundarySelection(
          top, {{"wrap", "wrap"}, {"wrap/leaf", "wrap/leaf"}}, 0),
      std::invalid_argument);
}

TEST_F(DesignBoundaryTests, RejectsUnconnectedInputAndPinlessInstance) {
  auto* block = createScalarBlock("BLOCK", {"A"}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  SNLInstance::create(top, block, NLName("u"));
  EXPECT_THROW(
      BoundarySelection(top, {{"u", "u"}}, 0), std::invalid_argument);

  auto* empty = createScalarBlock("EMPTY", {}, {});
  auto* emptyTop = SNLDesign::create(designs_, NLName("empty_top"));
  SNLInstance::create(emptyTop, empty, NLName("u"));
  EXPECT_THROW(
      BoundarySelection(emptyTop, {{"u", "u"}}, 0), std::invalid_argument);
}

TEST_F(DesignBoundaryTests, RejectsBoundaryInputWithoutADriver) {
  auto* block = createScalarBlock("INPUT_ONLY", {"A"}, {});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  auto* undriven = SNLScalarNet::create(top, NLName("undriven"));
  auto* instance = SNLInstance::create(top, block, NLName("u"));
  instance->getInstTerm(block->getScalarTerm(NLName("A")))->setNet(undriven);

  EXPECT_THROW(
      BoundarySelection(top, {{"u", "u"}}, 0), std::invalid_argument);
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
      BoundarySelection(top, {{"u", "u"}}, 0), std::invalid_argument);
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
    BoundarySelection boundary(top, pairs, side);
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

TEST_F(DesignBoundaryTests, RejectsVirtualOutputDisplayNameCollision) {
  auto* block = createScalarBlock("BLOCK", {}, {"Y"});
  auto* top = SNLDesign::create(designs_, NLName("top"));
  SNLInstance::create(top, block, NLName("u"));
  const auto name = BoundarySelection(top, {{"u", "u"}}, 0)
                        .getPorts().front().topTermName;
  SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName(name));
  expectInvalidBoundary(top, {{"u", "u"}}, 0, "collides with existing term");
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
