// Copyright 2024-2025 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>
#include <string>

#include "gtest/gtest.h"

#include "BuildPrimaryOutputClauses.h"
#include "ConstantPropagation.h"
#include "MiterStrategy.h"
#include "NLLibraryTruthTables.h"
#include "NLUniverse.h"
#include "NetlistGraph.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLDesignModeling.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "SNLPath.h"
#include "SNLCapnP.h"
#include "DNL.h"

using namespace naja;
using namespace naja::NL;
using namespace naja::NAJA_OPT;
using namespace KEPLER_FORMAL;

namespace {

void executeCommand(const std::string& command) {
  int result = system(command.c_str());
  if (result != 0) {
    std::cerr << "Command execution failed." << std::endl;
  }
}

}  // namespace

class MiterTests : public ::testing::Test {
 protected:
  MiterTests() {
    // You can do set-up work for each test here
  }
  ~MiterTests() override {
    // You can do clean-up work that doesn't throw exceptions here
  }
  void SetUp() override {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }
  void TearDown() override {
    // Code here will be called immediately after each test (right
    // before the destructor).
    // Destroy the SNL
    naja::DNL::destroy();
    NLUniverse::get()->destroy();
    KEPLER_FORMAL::BoolExprCache::destroy();
  }
};

TEST(HelloTest, ReturnsHelloWorld) {
  EXPECT_EQ(false, false);
}

TEST_F(MiterTests, TestMiterAND) {
  // 1. Create SNL
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  // 2. Create a top model with one output
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  auto topOut2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out2"));
  // 3. create a logic_0 model
  SNLDesign* logic0 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC0"));
  // add output to logic0
  auto logic0Out =
      SNLScalarTerm::create(logic0, SNLTerm::Direction::Output, NLName("out"));
  // 4. create a logic_1 model
  SNLDesign* logic1 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC1"));
  // add output to logic0
  auto logic1Out =
      SNLScalarTerm::create(logic1, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(logic0, SNLTruthTable(0, 0));
  SNLDesignModeling::setTruthTable(logic1, SNLTruthTable(0, 1));
  NLLibraryTruthTables::construct(library);
  // 5. create a logic_0 instace in top
  SNLInstance* inst1 = SNLInstance::create(top, logic0, NLName("logic0"));
  // 6. create a logic_1 instace in top
  SNLInstance* inst2 = SNLInstance::create(top, logic1, NLName("logic1"));
  // 7. create a and model
  SNLDesign* andModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("AND"));

  // add 2 inputs and 1 output to and
  auto andIn1 =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in1"));
  auto andIn2 =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in2"));
  auto andOut = SNLScalarTerm::create(andModel, SNLTerm::Direction::Output,
                                      NLName("out"));
  // 8. create a and instance in top
  SNLInstance* inst3 = SNLInstance::create(top, andModel, NLName("and"));
  SNLInstance* inst4 = SNLInstance::create(top, andModel, NLName("and2"));
  // set truth table for and model
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8));
  // 9. connect all instances inputs
  SNLNet* net1 = SNLScalarNet::create(top, NLName("logic_0_net"));
  net1->setType(SNLNet::Type::Assign0);
  SNLNet* net2 = SNLScalarNet::create(top, NLName("logic_1_net"));
  net2->setType(SNLNet::Type::Assign1);
  SNLNet* net3 = SNLScalarNet::create(top, NLName("and_output_net"));
  SNLNet* net4 = SNLScalarNet::create(top, NLName("and2_output_net"));
  // connect logic0 to and
  inst1->getInstTerm(logic0Out)->setNet(net1);

  inst4->getInstTerm(andIn1)->setNet(net2);
  inst4->getInstTerm(andIn2)->setNet(net2);
  // connect logic1 to and
  inst2->getInstTerm(logic1Out)->setNet(net2);
  inst3->getInstTerm(andIn2)->setNet(net1);
  inst3->getInstTerm(andIn1)->setNet(net4);
  // connect the and instance output to the top output
  inst3->getInstTerm(andOut)->setNet(net3);
  topOut->setNet(net3);
  inst4->getInstTerm(andOut)->setNet(net4);
  topOut2->setNet(net4);
  // 11. create DNL
  get();
  // 12. create a constant propagation object
  {
    std::string dotFileName(
        std::string(std::string("./beforeCP") + std::string(".dot")));
    std::string svgFileName(
        std::string(std::string("./beforeCP") + std::string(".svg")));
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
                               std::string(" -o ") + svgFileName)
                       .c_str());
  }
  ConstantPropagation cp;
  // 13. collect the constants
  // cp.collectConstants();
  // 14. run the constant propagation
  {
    BuildPrimaryOutputClauses miter;
    miter.build();
    for (const auto& po : miter.getPOs()) {
      std::cout << "PO: " << po->toString() << std::endl;
    }
  }

  cp.run();
  // 15. check the output value of the top instance
  {
    std::string dotFileName(
        std::string(std::string("./afterCP") + std::string(".dot")));
    std::string svgFileName(
        std::string(std::string("./afterCP") + std::string(".svg")));
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
                               std::string(" -o ") + svgFileName)
                       .c_str());
  }
  {
    BuildPrimaryOutputClauses miter;
    miter.build();
    for (const auto& po : miter.getPOs()) {
      std::cout << "PO: " << po->toString() << std::endl;
    }
  }
  naja::DNL::destroy();
}

TEST_F(MiterTests, TestMiterANDNonConstant) {
  printf("[TEST] MiterTests.TestMiterANDNonConstant\n");
  // 1. Create NL universe and DB
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);

  // 2. Create primitives library and register truth tables
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  NLLibraryTruthTables::construct(library);

  // 3. Create top design with two inputs and two outputs
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);

  auto topOut  = SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  auto topOut2 = SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out2"));
  auto topIn1  = SNLScalarTerm::create(top, SNLTerm::Direction::Input,  NLName("In1"));
  auto topIn2  = SNLScalarTerm::create(top, SNLTerm::Direction::Input,  NLName("In2"));

  // 4. Create an AND model (primitive) and its terms
  SNLDesign* andModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("AND"));
  auto andIn1 = SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in1"));
  auto andIn2 = SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in2"));
  auto andOut = SNLScalarTerm::create(andModel, SNLTerm::Direction::Output, NLName("out"));

  // 5. Create two instances of the AND model in top
  SNLInstance* instA = SNLInstance::create(top, andModel, NLName("andA"));
  SNLInstance* instB = SNLInstance::create(top, andModel, NLName("andB"));

  // 6. Set the truth table for the AND model (2-input AND = mask 0b1000 == 8)
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8));

  // 7. Create nets
  SNLNet* netTopIn1 = SNLScalarNet::create(top, NLName("top_in1_net"));
  SNLNet* netTopIn2 = SNLScalarNet::create(top, NLName("top_in2_net"));
  SNLNet* netAndAOut = SNLScalarNet::create(top, NLName("andA_output_net"));
  SNLNet* netAndBOut = SNLScalarNet::create(top, NLName("andB_output_net"));

  // 8. Connect top-level inputs to nets
  topIn1->setNet(netTopIn1);
  topIn2->setNet(netTopIn2);

  // 9. Wire instance inputs/outputs deliberately (avoid accidental self-wiring)
  instA->getInstTerm(andIn1)->setNet(netTopIn1);
  instA->getInstTerm(andIn2)->setNet(netTopIn2);
  instA->getInstTerm(andOut)->setNet(netAndAOut);
  topOut->setNet(netAndAOut);

  instB->getInstTerm(andIn1)->setNet(netTopIn2); // both inputs tied to topIn2
  instB->getInstTerm(andIn2)->setNet(netTopIn2);
  instB->getInstTerm(andOut)->setNet(netAndBOut);
  topOut2->setNet(netAndBOut);

  // 10. Initialize DNL subsystem
  naja::DNL::get();

  // 11. Optional: dump before-CP dot for offline inspection
  {
    std::string dotFileName = "./beforeCP.dot";
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    std::cerr << "[INFO] Wrote " << dotFileName << " for inspection.\n";
  }

  // 12. Run constant propagation
  ConstantPropagation cp;
  cp.run();

  // 13. Build primary output clauses (miter)
  BuildPrimaryOutputClauses miter;
  miter.collect();
  miter.build();

  const auto& pos = miter.getPOs();
  std::cout << "[INFO] miter.getPOs().size() = " << pos.size() << std::endl;

  if (pos.empty()) {
    // When no POs are produced, write an after-CP dot and fail with diagnostics.
    std::string dotFileName = "./afterCP_debug.dot";
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    std::cerr << "[DIAGNOSTIC] BuildPrimaryOutputClauses produced zero POs. "
                 "Wrote "
              << dotFileName << " for inspection.\n";
    FAIL() << "No primary outputs generated; inspect " << dotFileName;
    // FAIL terminates the test, so no further actions here.
  }

  // 14. Print POs for debugging and make permissive assertions
  for (const auto& po : pos) {
    std::cout << "PO: " << po->toString() << std::endl;
  }

  ASSERT_GE(pos.size(), 2u);
  // Basic sanity checks: strings are non-empty
  EXPECT_FALSE(pos[0]->toString().empty());
  EXPECT_FALSE(pos[1]->toString().empty());
}


TEST_F(MiterTests, TestMiterANDNonConstantWithSequentialElements) {
  printf("[TEST] MiterTests.TestMiterANDNonConstantWithSequentialElements\n");
  // 1. Create SNL
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  // 2. Create a top model with one output
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  auto topOut2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out2"));
  auto topIn1 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("In1"));
  auto topIn2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("In2"));
  NLLibraryTruthTables::construct(library);
  // 7. create a and model
  SNLDesign* andModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("AND"));

  // add 2 inputs and 1 output to and
  auto andIn1 =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in1"));
  auto andIn2 =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in2"));
  auto andOut = SNLScalarTerm::create(andModel, SNLTerm::Direction::Output,
                                      NLName("out"));

  // Create an FF
  SNLDesign* ffModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("FF"));
  // add D, CLK and Q
  auto ffD =
      SNLScalarTerm::create(ffModel, SNLTerm::Direction::Input, NLName("D"));
  auto ffCLK =
      SNLScalarTerm::create(ffModel, SNLTerm::Direction::Input, NLName("CLK"));
  auto ffQ =
      SNLScalarTerm::create(ffModel, SNLTerm::Direction::Output, NLName("Q"));
  // Set sequential dependecies to CLK
  SNLDesignModeling::addInputsToClockArcs({ffD}, {ffCLK});
  SNLDesignModeling::addClockToOutputsArcs({ffCLK}, {ffQ});

  // Create ff instance under top
  SNLInstance* instFF = SNLInstance::create(top, ffModel, NLName("ff0"));

  // 8. create a and instance in top
  SNLInstance* inst3 = SNLInstance::create(top, andModel, NLName("and"));
  SNLInstance* inst4 = SNLInstance::create(top, andModel, NLName("and2"));
  // set truth table for and model
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8));
  // 9. connect all instances inputs
  SNLNet* net1 = SNLScalarNet::create(top, NLName("top_in1_net"));
  SNLNet* net2 = SNLScalarNet::create(top, NLName("top_in2_net"));
  SNLNet* net3 = SNLScalarNet::create(top, NLName("and_output_net"));
  SNLNet* net4 = SNLScalarNet::create(top, NLName("and2_output_net"));
  SNLNet* net5 = SNLScalarNet::create(top, NLName("ffD"));
  SNLNet* net6 = SNLScalarNet::create(top, NLName("ffCLK"));
  // connect logic0 to and
  topIn1->setNet(net1);

  inst4->getInstTerm(andIn1)->setNet(net2);
  inst4->getInstTerm(andIn2)->setNet(net2);
  // connect logic1 to and
  instFF->getInstTerm(ffQ)->setNet(net2);
  instFF->getInstTerm(ffD)->setNet(net1);
  instFF->getInstTerm(ffCLK)->setNet(net6);
  inst3->getInstTerm(andIn2)->setNet(net1);
  inst3->getInstTerm(andIn1)->setNet(net4);
  // connect the and instance output to the top output
  inst3->getInstTerm(andOut)->setNet(net3);
  topOut->setNet(net3);
  inst4->getInstTerm(andOut)->setNet(net4);
  topOut2->setNet(net4);
  topIn1->setNet(net1);
  topIn2->setNet(net6);
  // 11. create DNL
  get();
  // 12. create a constant propagation object
  {
    std::string dotFileName(
        std::string(std::string("./beforeCP") + std::string(".dot")));
    std::string svgFileName(
        std::string(std::string("./beforeCP") + std::string(".svg")));
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
                               std::string(" -o ") + svgFileName)
                       .c_str());
  }
  ConstantPropagation cp;
  // 13. collect the constants
  // cp.collectConstants();
  // 14. run the constant propagation
  {
    BuildPrimaryOutputClauses miter;
    miter.collect();
    miter.build();
    for (const auto& po : miter.getPOs()) {
      std::cout << "PO: " << po->toString() << std::endl;
    }
  }

  cp.run();
  // 15. check the output value of the top instance
  {
    std::string dotFileName(
        std::string(std::string("./afterCP") + std::string(".dot")));
    std::string svgFileName(
        std::string(std::string("./afterCP") + std::string(".svg")));
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
                               std::string(" -o ") + svgFileName)
                       .c_str());
  }
  {
    BuildPrimaryOutputClauses pc;
    pc.collect();
    pc.build();
    // print inputs
    for (naja::DNL::DNLID id : pc.getInputs()) {
      DNLTerminalFull term = naja::DNL::get()->getDNLTerminalFromID(id);
        std::cout << "Input: " << term.getSnlBitTerm()->getName().getString() << " ID=" << id << std::endl;
    }
    // print outputs
    for (naja::DNL::DNLID id : pc.getOutputs()) {
      DNLTerminalFull term = naja::DNL::get()->getDNLTerminalFromID(id);
        std::cout << "Output: " << term.getSnlBitTerm()->getName().getString() << " ID=" << id << std::endl;
    }
    for (const auto& po : pc.getPOs()) {
      std::cout << "PO: " << po->toString() << std::endl;
    }
    printf("%s\n", pc.getPOs()[0]->toString().c_str());
    //EXPECT_TRUE(miter.getPOs()[0]->toString() == std::string("((6 ∧ 6) ∧ 2)"));
    EXPECT_TRUE(pc.getPOs()[0]->toString() == std::string("2 AND 4"));
    printf("%s\n", pc.getPOs()[1]->toString().c_str());
    //EXPECT_TRUE(miter.getPOs()[1]->toString() == std::string("(6 ∧ 6)"));
    EXPECT_TRUE(pc.getPOs()[1]->toString() == std::string("4"));
    printf("%s\n", pc.getPOs()[2]->toString().c_str());
    //EXPECT_TRUE(miter.getPOs()[2]->toString() == std::string("2"));
    EXPECT_TRUE(pc.getPOs()[2]->toString() == std::string("2"));
    printf("%s\n", pc.getPOs()[3]->toString().c_str());
    //EXPECT_TRUE(miter.getPOs()[3]->toString() == std::string("3"));
    EXPECT_TRUE(pc.getPOs()[3]->toString() == std::string("3"));
  }
}

// 1. create a circuit of 2 inputs that drives and AND gate that drives top output
// 2. clone the the top and chain an inverter to the AND output
// 3. verify that the miter strategy detects the difference
TEST_F(MiterTests, TestMiterAndWithChainedInverter) {
  // 1. Create SNL
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  NLLibrary* libraryDesigns =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  // 2. Create a top model with one output
  SNLDesign* top = SNLDesign::create(libraryDesigns, SNLDesign::Type::Standard,
                                     NLName("top"));
  univ->setTopDesign(top);
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  auto topOut2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out2"));
  auto topIn1 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("In1"));
  auto topIn2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("In2"));
  // add another 2 inputs
  auto topIn3 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("In3"));
  auto topIn4 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("In4"));
  NLLibraryTruthTables::construct(library);
  // 7. create a and model
  SNLDesign* andModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("AND"));

  // add 2 inputs and 1 output to and
  auto andIn1 =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in1"));
  auto andIn2 =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in2"));
  auto andOut = SNLScalarTerm::create(andModel, SNLTerm::Direction::Output,
                                      NLName("out"));

  // set truth table for and model
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8));
  // 8. create an inverter model
  SNLDesign* inverterModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("INV"));
  // set truth table for inverter model
  auto invIn =
      SNLScalarTerm::create(inverterModel, SNLTerm::Direction::Input, NLName("in"));
  auto invOut =
      SNLScalarTerm::create(inverterModel, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(inverterModel, SNLTruthTable(1, 1));

  // create and instance in top
  SNLInstance* instAnd = SNLInstance::create(top, andModel, NLName("and"));

  // connect inputs to the and instance
  SNLNet* net1 = SNLScalarNet::create(top, NLName("top_in1_net"));
  SNLNet* net2 = SNLScalarNet::create(top, NLName("top_in2_net"));
  SNLNet* net3 = SNLScalarNet::create(top, NLName("and_output_net"));
  // connect inputs to the top instance
  topIn1->setNet(net1);
  topIn2->setNet(net2);
  // connect the and instance inputs
  instAnd->getInstTerm(andIn1)->setNet(net1);
  instAnd->getInstTerm(andIn2)->setNet(net2);
  // connect the and instance output to the top output
  instAnd->getInstTerm(andOut)->setNet(net3);
  topOut->setNet(net3);

  // add another and instance in top at the same manner
  SNLInstance* instAnd2 = SNLInstance::create(top, andModel, NLName("and2"));
  // connect the and instance inputs
  // connect inputs 2 and 3 to the top instance
  // create needed nets
  SNLNet* net4In1 = SNLScalarNet::create(top, NLName("top_in3_net"));
  SNLNet* net4In2 = SNLScalarNet::create(top, NLName("top_in4_net"));
  topIn3->setNet(net4In1);
  topIn4->setNet(net4In2);
  // connect the and instance inputs
  instAnd2->getInstTerm(andIn1)->setNet(net4In1);
  instAnd2->getInstTerm(andIn2)->setNet(net4In2);

  // connect the and instance output to the top output
  SNLNet* net4Out = SNLScalarNet::create(top, NLName("and2_output_net_out"));
  instAnd2->getInstTerm(andOut)->setNet(net4Out);
  topOut2->setNet(net4Out);


  {
    // dump top to naja_if(CapProto)
    std::filesystem::path outputPath("./top.capnp");
    SNLCapnP::dump(db, outputPath);
  }
  // Dump visual
  {
    std::string dotFileName(
        std::string(std::string("./beforeEdit") + std::string(".dot")));
    std::string svgFileName(
        std::string(std::string("./beforeEdit") + std::string(".svg")));
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
                               std::string(" -o ") + svgFileName)
                       .c_str());
  }
  // clone the top design
  SNLDesign* topClone = top->clone(NLName("topClone"));
  // create an inverter instance in the clone
  SNLInstance* instInv = SNLInstance::create(top, inverterModel, NLName("inv"));
  // connect the inverter input to the and output
  SNLNet* net4 = SNLScalarNet::create(top, NLName("and_output_net_clone"));
  instAnd->getInstTerm(andOut)->setNet(net4);
  instInv->getInstTerm(invIn)->setNet(net4);
  // connect the inverter output to the top output
  SNLNet* net5 = SNLScalarNet::create(top, NLName("top_output_net_clone"));
  instInv->getInstTerm(invOut)->setNet(net5);
  topOut->setNet(net5);

  // dump visual
  {
    std::string dotFileName(
        std::string(std::string("./afterEdit") + std::string(".dot")));
    std::string svgFileName(
        std::string(std::string("./afterEdit") + std::string(".svg")));
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
                               std::string(" -o ") + svgFileName)
                       .c_str());
  }

  // test the miter strategy
  {
    MiterStrategy MiterS(top, topClone, "CaseC");
    EXPECT_FALSE(MiterS.run());
  }
  {
    // dump top to naja_if(CapProto)
    std::filesystem::path outputPath("./topEdited1.capnp");
    SNLCapnP::dump(db, outputPath);
  }
  //Check output of binary ../../../src/bin/kepler-formal on the 2 capnp files
  executeCommand(
      std::string("../../../src/bin/kepler-formal -naja_if ./top.capnp ./topEdited1.capnp")
          .c_str());  
  // look for "DIFFERENT" in the file ./miter_log_1.txt
  // open the file  
  std::ifstream miterLogFile("./miter_log_0.txt");
  std::string line;
  bool foundDifferent = false;
  if (miterLogFile.is_open()) {
    while (getline(miterLogFile, line)) {
      if (line.find("DIFFERENT") != std::string::npos) {
        foundDifferent = true;
        break;
      }
    }
    miterLogFile.close();
  }
  EXPECT_TRUE(foundDifferent);
  // chain another inverter to the first inverter
  SNLInstance* instInv2 = SNLInstance::create(top, inverterModel, NLName("inv2"));
  // connect the second inverter input to the first inverter output
  SNLNet* net6 = SNLScalarNet::create(top, NLName("inv_output_net_clone"));
  instInv->getInstTerm(invOut)->setNet(net6);
  instInv2->getInstTerm(invIn)->setNet(net6);
  // connect the second inverter output to the top output
  SNLNet* net7 = SNLScalarNet::create(top, NLName("top_output_net_clone2"));
  instInv2->getInstTerm(invOut)->setNet(net7);
  topOut->setNet(net7);
  // test the miter strategy again
  {
    MiterStrategy MiterS(top, topClone, "CaseD");
    EXPECT_TRUE(MiterS.run());
  }
  {
    // dump top to naja_if(CapProto)
    std::filesystem::path outputPath("./topEdited2.capnp");
    SNLCapnP::dump(db, outputPath);
  }

  //Check output of binary ../../../src/bin/kepler-formal on the 2 capnp files
  executeCommand(
      std::string("../../../src/bin/kepler-formal -naja_if ./top.capnp ./topEdited2.capnp")
          .c_str()); 
  // look for "IDENTICAL" in the file ./miter_log_2.txt
  // open the file
  std::ifstream miterLogFile2("./miter_log_1.txt");
  bool foundIdentical = false;
  if (miterLogFile2.is_open()) {
    while (getline(miterLogFile2, line)) {
      if (line.find("IDENTICAL") != std::string::npos) {
        foundIdentical = true;
        break;
      }
    }
    miterLogFile2.close();
  }
  EXPECT_TRUE(foundIdentical);
  
}

// Required main function for Google Test
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// ---------------------- Tests appended for coverage (subprocess approach, tolerant) ----------------------
// Append this block at the end of the file (after main).

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>

// Path to the kepler-formal CLI binary used by the project tests.
// Adjust this path if your binary is located elsewhere.
static const char* KEPLER_BIN = "../../../src/bin/kepler-formal";

// Helper to run the CLI binary with arguments in a subprocess using std::system.
// Returns the program's exit code (child exit status) when available, otherwise EXIT_FAILURE.
static int run_kepler_cli_with_args(const std::vector<std::string>& args) {
  std::string cmd;
  cmd += KEPLER_BIN;
  for (const auto& a : args) {
    cmd += " ";
    // naive quoting: wrap in single quotes and escape any single quotes inside
    std::string quoted = "'";
    for (char c : a) {
      if (c == '\'') quoted += "'\\''";
      else quoted.push_back(c);
    }
    quoted += "'";
    cmd += quoted;
  }

  int rc = std::system(cmd.c_str());
  if (rc == -1) {
    // system failed to start the process
    return EXIT_FAILURE;
  }

#if defined(_WIN32)
  // On Windows, system returns the exit code directly
  return rc;
#else
  // On POSIX, interpret wait status
  if (WIFEXITED(rc)) {
    return WEXITSTATUS(rc);
  } else {
    // Abnormal termination (signal, etc.) -> treat as failure
    return EXIT_FAILURE;
  }
#endif
}

TEST(KeplerCliSubprocessTests, BinaryExists) {
  std::filesystem::path p(KEPLER_BIN);
  bool exists = std::filesystem::exists(p);
  if (!exists) {
    GTEST_SKIP() << "kepler-formal binary not found at " << KEPLER_BIN << "; skipping CLI subprocess tests.";
  }
  EXPECT_TRUE(std::filesystem::is_regular_file(p));
}

TEST(KeplerCliSubprocessTests, PrintUsageOnNoArgs) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  int rc = run_kepler_cli_with_args({});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, HelpFlagReturnsSuccess) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  int rc = run_kepler_cli_with_args({"--help"});
  EXPECT_EQ(rc, EXIT_SUCCESS);
  rc = run_kepler_cli_with_args({"-h"});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, MissingConfigFileArgument) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  int rc = run_kepler_cli_with_args({"--config"});
  EXPECT_NE(rc, EXIT_SUCCESS);
  rc = run_kepler_cli_with_args({"-c"});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ConfigFileNotFoundReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  std::string tmpPath = "./nonexistent_config_12345.yaml";
  int rc = run_kepler_cli_with_args({"--config", tmpPath});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ConfigUnrecognizedFormatReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  std::filesystem::path tmp = std::filesystem::temp_directory_path() / "kepler_test_bad_format.yaml";
  {
    std::ofstream ofs(tmp);
    ofs << "format: unknown_format\n";
    ofs << "input_paths:\n  - a\n  - b\n";
    ofs.close();
  }
  int rc = run_kepler_cli_with_args({"--config", tmp.string()});
  EXPECT_NE(rc, EXIT_SUCCESS);
  std::filesystem::remove(tmp);
}

TEST(KeplerCliSubprocessTests, ConfigSnlFormatLoadFailureReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  std::filesystem::path tmp = std::filesystem::temp_directory_path() / "kepler_test_snl.yaml";
  {
    std::ofstream ofs(tmp);
    ofs << "format: snl\n";
    ofs << "input_paths:\n  - /path/does/not/exist1.snl\n  - /path/does/not/exist2.snl\n";
    ofs.close();
  }
  int rc = run_kepler_cli_with_args({"--config", tmp.string()});
  // Accept any non-success result (normal nonzero exit or abnormal termination)
  EXPECT_NE(rc, EXIT_SUCCESS);
  std::filesystem::remove(tmp);
}

TEST(KeplerCliSubprocessTests, CliUnrecognizedFormatReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  int rc = run_kepler_cli_with_args({"-badformat", "a", "b"});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

// Program prints usage and returns success when argc < 4; keep that behavior expected.
TEST(KeplerCliSubprocessTests, CliNotEnoughPathsReturnsSuccess) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  int rc = run_kepler_cli_with_args({"-verilog", "only_one_path.v"});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, CliNajaIfFormatButMissingFilesReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";
  int rc = run_kepler_cli_with_args({"-naja_if", "/no/such/file1.capnp", "/no/such/file2.capnp"});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ConfigParsingViaFilesCoversYamlToVectorBehavior) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";

  // 1) Sequence of scalars -> valid config with two input_paths should proceed to further checks.
  std::filesystem::path tmpSeq = std::filesystem::temp_directory_path() / "kepler_test_seq.yaml";
  {
    std::ofstream ofs(tmpSeq);
    ofs << "format: verilog\n";
    ofs << "input_paths:\n";
    ofs << "  - fileA.v\n";
    ofs << "  - fileB.v\n";
    ofs << "liberty_files:\n";
    ofs << "  - lib1.lib\n";
    ofs.close();
  }
  {
    int rc = run_kepler_cli_with_args({"--config", tmpSeq.string()});
    EXPECT_NE(rc, EXIT_SUCCESS); // files missing or parser errors -> not success
  }
  std::filesystem::remove(tmpSeq);

  // 2) Scalar node for input_paths (invalid shape) -> should fail
  std::filesystem::path tmpScalar = std::filesystem::temp_directory_path() / "kepler_test_scalar.yaml";
  {
    std::ofstream ofs(tmpScalar);
    ofs << "format: verilog\n";
    ofs << "input_paths: \"not-a-sequence\"\n";
    ofs.close();
  }
  {
    int rc = run_kepler_cli_with_args({"--config", tmpScalar.string()});
    EXPECT_NE(rc, EXIT_SUCCESS);
  }
  std::filesystem::remove(tmpScalar);

  // 3) Null node (empty YAML) -> should fail
  std::filesystem::path tmpNull = std::filesystem::temp_directory_path() / "kepler_test_null.yaml";
  {
    std::ofstream ofs(tmpNull);
    ofs << "# empty config\n";
    ofs.close();
  }
  {
    int rc = run_kepler_cli_with_args({"--config", tmpNull.string()});
    EXPECT_NE(rc, EXIT_SUCCESS);
  }
  std::filesystem::remove(tmpNull);

  // 4) Sequence of non-scalars (maps) for input_paths -> should fail
  std::filesystem::path tmpSeqMaps = std::filesystem::temp_directory_path() / "kepler_test_seqmaps.yaml";
  {
    std::ofstream ofs(tmpSeqMaps);
    ofs << "format: verilog\n";
    ofs << "input_paths:\n";
    ofs << "  - {a: 1}\n";
    ofs << "  - {b: 2}\n";
    ofs.close();
  }
  {
    int rc = run_kepler_cli_with_args({"--config", tmpSeqMaps.string()});
    EXPECT_NE(rc, EXIT_SUCCESS);
  }
  std::filesystem::remove(tmpSeqMaps);
}

TEST_F(MiterTests, CoverDiff) {
  // 1. Create SNL
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  // 2. Create a top model with one output
  NLLibrary* libraryDesigns =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  SNLDesign* top =
      SNLDesign::create(libraryDesigns, SNLDesign::Type::Standard, NLName("top"));
  univ->setTopDesign(top);
  auto topin =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto topin2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in2"));
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  auto topOut2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out2"));
  // 3. create a logic_0 model
  SNLDesign* logic0 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC0"));
  // add output to logic0
  auto logic0Out =
      SNLScalarTerm::create(logic0, SNLTerm::Direction::Output, NLName("out"));
  // 4. create a logic_1 model
  SNLDesign* logic1 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC1"));
  // add output to logic0
  auto logic1Out =
      SNLScalarTerm::create(logic1, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(logic0, SNLTruthTable(0, 0));
  SNLDesignModeling::setTruthTable(logic1, SNLTruthTable(0, 1));
  NLLibraryTruthTables::construct(library);
  // 5. create a logic_0 instace in top
  SNLInstance* inst1 = SNLInstance::create(top, logic0, NLName("logic0"));
  // 6. create a logic_1 instace in top
  SNLInstance* inst2 = SNLInstance::create(top, logic1, NLName("logic1"));
  // 7. create a and model
  SNLDesign* seqModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("AND"));

  // add 2 inputs and 1 output to and
  auto andIn1 =
      SNLScalarTerm::create(seqModel, SNLTerm::Direction::Input, NLName("in1"));
  auto andIn2 =
      SNLScalarTerm::create(seqModel, SNLTerm::Direction::Input, NLName("in2"));
  auto andOut = SNLScalarTerm::create(seqModel, SNLTerm::Direction::Output,
                                      NLName("out"));
  // 8. create a and instance in top
  SNLInstance* inst3 = SNLInstance::create(top, seqModel, NLName("and"));
  SNLInstance* inst4 = SNLInstance::create(top, seqModel, NLName("and2"));
  // set truth table for and model
  //SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8));
  // 9. connect all instances inputs
  SNLNet* net1 = SNLScalarNet::create(top, NLName("logic_0_net"));
  net1->setType(SNLNet::Type::Assign0);
  SNLNet* net2 = SNLScalarNet::create(top, NLName("logic_1_net"));
  net2->setType(SNLNet::Type::Assign1);
  SNLNet* net3 = SNLScalarNet::create(top, NLName("and_output_net"));
  SNLNet* net4 = SNLScalarNet::create(top, NLName("and2_output_net"));
  // connect logic0 to and
  inst1->getInstTerm(logic0Out)->setNet(net1);

  inst4->getInstTerm(andIn1)->setNet(net2);
  inst4->getInstTerm(andIn2)->setNet(net2);
  // connect logic1 to and
  inst2->getInstTerm(logic1Out)->setNet(net2);
  inst3->getInstTerm(andIn2)->setNet(net1);
  inst3->getInstTerm(andIn1)->setNet(net4);
  // connect the and instance output to the top output
  inst3->getInstTerm(andOut)->setNet(net3);
  topOut->setNet(net3);
  inst4->getInstTerm(andOut)->setNet(net4);
  topOut2->setNet(net4);

  SNLDesign* topClone0 = top->clone(NLName("topClone0"));
  SNLNet* netC0a = SNLScalarNet::create(topClone0, NLName("netC0a"));
  SNLNet* netC0b = SNLScalarNet::create(topClone0, NLName("netC0b"));
  auto andC0 = topClone0->getInstance(NLName("and"));
  andC0->getInstTerm(andIn1->getID())->setNet(netC0a);
  andC0->getInstTerm(andIn2->getID())->setNet(netC0b);
  topClone0->getBitTerm(topin->getID(), 0)->setNet(netC0a);
  SNLInstance* constC0 = SNLInstance::create(topClone0, logic0, NLName("logic0C0"));
  constC0->getInstTerm(logic0Out->getID())->setNet(netC0b);

  SNLDesign* topClone1 = top->clone(NLName("topClone1"));
  SNLNet* netC1a = SNLScalarNet::create(topClone1, NLName("netC1a"));
  SNLNet* netC1b = SNLScalarNet::create(topClone1, NLName("netC1b"));
  auto andC1 = topClone1->getInstance(NLName("and"));
  andC1->getInstTerm(andIn1->getID())->setNet(netC1a);
  andC1->getInstTerm(andIn2->getID())->setNet(netC1b);
  SNLInstance* constC1 = SNLInstance::create(topClone1, logic0, NLName("logic0C1"));
  SNLDesign* inverterModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("INV"));
  // set truth table for inverter model
  auto invIn =
      SNLScalarTerm::create(inverterModel, SNLTerm::Direction::Input, NLName("in"));
  auto invOut =
      SNLScalarTerm::create(inverterModel, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(inverterModel, SNLTruthTable(1, 1));
  auto inverterC1 = SNLInstance::create(topClone1, inverterModel, NLName("inverterC1"));
  constC1->getInstTerm(logic0Out->getID())->setNet(netC1a);
  auto netC1invOut = SNLScalarNet::create(topClone1, NLName("netC1invOut"));
  inverterC1->getInstTerm(invIn->getID())->setNet(netC1b);
  inverterC1->getInstTerm(invOut->getID())->setNet(netC1invOut);
  andC1->getInstTerm(andIn1->getID())->setNet(netC1invOut);
  topClone1->getBitTerm(topin2->getID(), 0)->setNet(netC1b);
  
  // 11. create DNL
  get(); 
  naja::DNL::destroy();
  MiterStrategy MiterS(topClone0, topClone1, "CaseD");
    EXPECT_FALSE(MiterS.run());
}

TEST(KeplerCliSubprocessTests, ExampleTestRun) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) GTEST_SKIP() << "kepler-formal binary missing";

  int rc = run_kepler_cli_with_args({"-verilog", "../../../../example/tinyrocket.v", "../../../../example/tinyrocket_edited.v", 
                                         "../../../../example/NangateOpenCellLibrary_typical.lib",
                                         "../../../../example/fakeram45_64x15.lib",
                                         "../../../../example/fakeram45_64x32.lib",
                                         "../../../../example/fakeram45_1024x32.lib"});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

// End of appended tests
