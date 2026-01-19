// Copyright 2024-2025 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only
//
// Unit tests for ScopeExtraction::collectVerificationScopes()
// Adjusted to avoid duplicate SNLDesign/SNLInstance name collisions and
// double-destroy. Made tests deterministic by ensuring differences are explicit
// (net wiring or model reuse).

#include <gtest/gtest.h>
#include <set>
#include <string>
#include <utility>

#include "ScopeExtraction.h"

#include "BuildPrimaryOutputClauses.h"
#include "ConstantPropagation.h"
#include "DNL.h"
#include "MiterStrategy.h"
#include "NLLibraryTruthTables.h"
#include "NLUniverse.h"
#include "NetlistGraph.h"
#include "SNLCapnP.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLPath.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"

using namespace naja;
using namespace naja::NL;
using namespace naja::NAJA_OPT;
using namespace KEPLER_FORMAL;

// Helper to run shell commands (used by the original test for dot->svg).
static void executeCommand(const std::string& command) {
  int result = system(command.c_str());
  if (result != 0) {
    std::cerr << "Command execution failed: " << command << std::endl;
  }
}

// Test helper derived from ScopeExtraction to expose internals for assertions.
class TestScopeExtraction : public ScopeExtraction {
 public:
  TestScopeExtraction() : ScopeExtraction() {}
  void setTop0(naja::NL::SNLDesign* d) { top0_ = d; }
  void setTop1(naja::NL::SNLDesign* d) { top1_ = d; }

  const std::set<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>>&
  getDesignsToVerify() const {
    return designsToVerify_;
  }

  void runCollect() { collectVerificationScopes(); }
};

namespace {

// Utility: create a minimal library + universe and return the created library.
// Caller is responsible for calling naja::DNL::destroy() and
// NLUniverse::get()->destroy() in test teardown if needed.
static NLLibrary* createUniverseAndLibrary(NLUniverse*& outUniv) {
  outUniv = NLUniverse::create();
  NLDB* db = NLDB::create(outUniv);
  NLLibrary* lib =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  return lib;
}

// Build a simple top design with two child instances and two outputs.
// The function returns the top design pointer and also fills the instance
// pointers and term pointers so tests can wire nets differently to create
// mismatches.
//
// IMPORTANT: to avoid duplicate-name collisions in the shared NLLibrary, this
// function appends a unique suffix to created names.
struct TopDesignBundle {
  SNLDesign* top = nullptr;
  SNLScalarTerm* topOutA = nullptr;
  SNLScalarTerm* topOutB = nullptr;
  SNLDesign* childA_model = nullptr;
  SNLDesign* childB_model = nullptr;
  SNLInstance* childA_inst = nullptr;
  SNLInstance* childB_inst = nullptr;
  SNLScalarTerm* childA_out = nullptr;
  SNLScalarTerm* childB_out = nullptr;
};

static TopDesignBundle buildSimpleTop(NLLibrary* lib,
                                      const std::string& topBaseName,
                                      const std::string& childANameBase,
                                      const std::string& childBNameBase) {
  static unsigned uniqueCounter = 0;
  const unsigned id = ++uniqueCounter;
  const std::string topName = topBaseName + "_" + std::to_string(id);
  const std::string childAName = childANameBase + "_" + std::to_string(id);
  const std::string childBName = childBNameBase + "_" + std::to_string(id);

  TopDesignBundle b;
  // Create top
  b.top = SNLDesign::create(lib, SNLDesign::Type::Primitive, NLName(topName));
  // two top outputs
  b.topOutA =
      SNLScalarTerm::create(b.top, SNLTerm::Direction::Output, NLName("outA"));
  b.topOutB =
      SNLScalarTerm::create(b.top, SNLTerm::Direction::Output, NLName("outB"));

  // child A: a primitive with one output
  b.childA_model =
      SNLDesign::create(lib, SNLDesign::Type::Primitive, NLName(childAName));
  b.childA_out = SNLScalarTerm::create(
      b.childA_model, SNLTerm::Direction::Output, NLName("out"));

  // child B: another primitive with one output
  b.childB_model =
      SNLDesign::create(lib, SNLDesign::Type::Primitive, NLName(childBName));
  b.childB_out = SNLScalarTerm::create(
      b.childB_model, SNLTerm::Direction::Output, NLName("out"));

  // instantiate both in top with unique instance names
  const std::string instAName = std::string("instA_") + std::to_string(id);
  const std::string instBName = std::string("instB_") + std::to_string(id);
  b.childA_inst = SNLInstance::create(b.top, b.childA_model, NLName(instAName));
  b.childB_inst = SNLInstance::create(b.top, b.childB_model, NLName(instBName));

  return b;
}

}  // namespace

class ScopeExtractionUnitTests : public ::testing::Test {
 protected:
  void SetUp() override {
    // create universe + library for each test
    lib_ = createUniverseAndLibrary(univ_);
  }

  void TearDown() override {
    // Clean up global singletons used by the SNL framework
    NLUniverse::get()->destroy();
    univ_ = nullptr;
    lib_ = nullptr;
  }

  NLUniverse* univ_ = nullptr;
  NLLibrary* lib_ = nullptr;
};

// Test case: identical designs should not be added to designsToVerify_ at the
// top level, and the algorithm should recurse into children (which are also
// identical).
TEST_F(ScopeExtractionUnitTests, IdenticalDesigns_NoVerificationNeeded) {
  // Build two top designs with identical structure (unique names internally)
  TopDesignBundle a = buildSimpleTop(lib_, "topA", "LOGIC0", "LOGIC1");
  TopDesignBundle b = buildSimpleTop(lib_, "topB", "LOGIC0", "LOGIC1");

  // Make the child models have deterministic truth tables so library tables
  // exist
  SNLDesignModeling::setTruthTable(a.childA_model, SNLTruthTable(0, 0));
  SNLDesignModeling::setTruthTable(a.childB_model, SNLTruthTable(0, 1));
  SNLDesignModeling::setTruthTable(b.childA_model, SNLTruthTable(0, 0));
  SNLDesignModeling::setTruthTable(b.childB_model, SNLTruthTable(0, 1));
  NLLibraryTruthTables::construct(lib_);

  // Wire nets identically for both tops
  SNLNet* a_net0 = SNLScalarNet::create(a.top, NLName("net0"));
  SNLNet* a_net1 = SNLScalarNet::create(a.top, NLName("net1"));
  a.childA_inst->getInstTerm(a.childA_out)->setNet(a_net0);
  a.childB_inst->getInstTerm(a.childB_out)->setNet(a_net1);
  a.topOutA->setNet(a_net0);
  a.topOutB->setNet(a_net1);

  SNLNet* b_net0 = SNLScalarNet::create(b.top, NLName("net0"));
  SNLNet* b_net1 = SNLScalarNet::create(b.top, NLName("net1"));
  b.childA_inst->getInstTerm(b.childA_out)->setNet(b_net0);
  b.childB_inst->getInstTerm(b.childB_out)->setNet(b_net1);
  b.topOutA->setNet(b_net0);
  b.topOutB->setNet(b_net1);

  TestScopeExtraction se;
  se.setTop0(a.top);
  se.setTop1(b.top);

  se.runCollect();

  const auto& toVerify = se.getDesignsToVerify();
  EXPECT_TRUE(toVerify.empty());
}

// Test case: different number of instances at top level -> top pair should be
// added
TEST_F(ScopeExtractionUnitTests, DifferentInstanceCount_TopAddedToVerify) {
  // Build top A with two children
  TopDesignBundle a = buildSimpleTop(lib_, "topA_diff", "LOGIC0", "LOGIC1");

  // Build top B with only one child (simulate different instance count)
  TopDesignBundle b;
  static unsigned localCounter = 0;
  ++localCounter;
  const std::string topBName =
      std::string("topB_diff_") + std::to_string(localCounter);
  b.top = SNLDesign::create(lib_, SNLDesign::Type::Primitive, NLName(topBName));
  b.topOutA =
      SNLScalarTerm::create(b.top, SNLTerm::Direction::Output, NLName("outA"));
  // create only one child model and instance (unique name)
  const std::string singleChildName =
      std::string("LOGIC0_single_") + std::to_string(localCounter);
  b.childA_model = SNLDesign::create(lib_, SNLDesign::Type::Primitive,
                                     NLName(singleChildName));
  b.childA_out = SNLScalarTerm::create(
      b.childA_model, SNLTerm::Direction::Output, NLName("out"));
  const std::string instAName =
      std::string("instA_single_") + std::to_string(localCounter);
  b.childA_inst = SNLInstance::create(b.top, b.childA_model, NLName(instAName));

  // Set truth tables and construct library tables
  SNLDesignModeling::setTruthTable(a.childA_model, SNLTruthTable(0, 0));
  SNLDesignModeling::setTruthTable(a.childB_model, SNLTruthTable(0, 1));
  SNLDesignModeling::setTruthTable(b.childA_model, SNLTruthTable(0, 0));
  NLLibraryTruthTables::construct(lib_);

  // Wire nets for A
  SNLNet* a_net0 = SNLScalarNet::create(a.top, NLName("net0"));
  SNLNet* a_net1 = SNLScalarNet::create(a.top, NLName("net1"));
  a.childA_inst->getInstTerm(a.childA_out)->setNet(a_net0);
  a.childB_inst->getInstTerm(a.childB_out)->setNet(a_net1);
  a.topOutA->setNet(a_net0);
  a.topOutB->setNet(a_net1);

  // Wire nets for B (only one)
  SNLNet* b_net0 = SNLScalarNet::create(b.top, NLName("net0"));
  b.childA_inst->getInstTerm(b.childA_out)->setNet(b_net0);
  b.topOutA->setNet(b_net0);

  TestScopeExtraction se;
  se.setTop0(a.top);
  se.setTop1(b.top);

  se.runCollect();

  const auto& toVerify = se.getDesignsToVerify();
  EXPECT_FALSE(toVerify.empty());
  EXPECT_EQ(toVerify.size(), 1u);
  auto pair = *toVerify.begin();
  EXPECT_EQ(pair.first, a.top);
  EXPECT_EQ(pair.second, b.top);
}

// Test case: same instance count but different child instance ordering/IDs ->
// added to verify Make the difference explicit by wiring nets differently in
// the second top as well.
TEST_F(ScopeExtractionUnitTests, SameCountDifferentChildIDs_AddedToVerify) {
  // Build two tops with same number of instances but different child model IDs
  TopDesignBundle a = buildSimpleTop(lib_, "topA_ids", "LOGIC_A", "LOGIC_B");
  TopDesignBundle b = buildSimpleTop(lib_, "topB_ids", "LOGIC_X", "LOGIC_Y");

  // Set truth tables so nets/terms exist
  SNLDesignModeling::setTruthTable(a.childA_model, SNLTruthTable(0, 0));
  SNLDesignModeling::setTruthTable(a.childB_model, SNLTruthTable(0, 1));
  SNLDesignModeling::setTruthTable(b.childA_model, SNLTruthTable(0, 0));
  SNLDesignModeling::setTruthTable(b.childB_model, SNLTruthTable(0, 1));
  NLLibraryTruthTables::construct(lib_);

  // Wire nets for A
  SNLNet* a_net0 = SNLScalarNet::create(a.top, NLName("net0"));
  SNLNet* a_net1 = SNLScalarNet::create(a.top, NLName("net1"));
  a.childA_inst->getInstTerm(a.childA_out)->setNet(a_net0);
  a.childB_inst->getInstTerm(a.childB_out)->setNet(a_net1);
  a.topOutA->setNet(a_net0);
  a.topOutB->setNet(a_net1);

  // Wire nets for B differently (explicit difference)
  SNLNet* b_net0 = SNLScalarNet::create(b.top, NLName("net0"));
  SNLNet* b_net1 = SNLScalarNet::create(b.top, NLName("net1"));
  // swap wiring in B to ensure mismatch is detected
  b.childA_inst->getInstTerm(b.childA_out)->setNet(b_net1);
  b.childB_inst->getInstTerm(b.childB_out)->setNet(b_net0);
  b.topOutA->setNet(b_net1);
  b.topOutB->setNet(b_net0);

  TestScopeExtraction se;
  se.setTop0(a.top);
  se.setTop1(b.top);

  se.runCollect();

  const auto& toVerify = se.getDesignsToVerify();
  EXPECT_FALSE(toVerify.empty());
  EXPECT_EQ(toVerify.size(), 1u);
  auto pair = *toVerify.begin();
  EXPECT_EQ(pair.first, a.top);
  EXPECT_EQ(pair.second, b.top);
}

// Test case: same child instances but nets differ in bit terms / inst terms ->
// added to verify
TEST_F(ScopeExtractionUnitTests, SameChildrenDifferentNets_AddedToVerify) {
  // Build two tops that will share the same child models (so child model
  // pointers match) Build the first top without helper function
  TopDesignBundle a;
  const std::string aTopName = "topA_nets_reuse";
  a.top = SNLDesign::create(lib_, SNLDesign::Type::Primitive, NLName(aTopName));
  TopDesignBundle b;
  const std::string bTopName =
      std::string("topB_nets_reuse_") + std::to_string(1);
  b.top = SNLDesign::create(lib_, SNLDesign::Type::Primitive, NLName(bTopName));
  // create one net on second top
  SNLNet* b_net0 = SNLScalarNet::create(b.top, NLName("net0"));

  TestScopeExtraction se;
  se.setTop0(a.top);
  se.setTop1(b.top);

  se.runCollect();

  const auto& toVerify = se.getDesignsToVerify();
  EXPECT_FALSE(toVerify.empty());
  EXPECT_EQ(toVerify.size(), 1u);
  auto pair = *toVerify.begin();
  EXPECT_EQ(pair.first, a.top);
  EXPECT_EQ(pair.second, b.top);
}

// Required main function for Google Test
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
