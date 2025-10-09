#include <gtest/gtest.h>
#include <string>
// SPDX-FileCopyrightText: 2024 The Naja authors
// <https://github.com/najaeda/naja/blob/main/AUTHORS>
//
// SPDX-License-Identifier: Apache-2.0

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

// TEST_F(MiterTests, TestMiterANDNonConstantWithSequentialElementsFormal) {
//   // 1. Create SNL
//   NLUniverse* univ = NLUniverse::create();
//   NLDB* db = NLDB::create(univ);
//   NLLibrary* library =
//       NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
//   NLLibrary* libraryDesigns =
//       NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
//   // 2. Create a top model with one output
//   SNLDesign* top = SNLDesign::create(libraryDesigns, SNLDesign::Type::Standard,
//                                      NLName("top"));
//   univ->setTopDesign(top);
//   auto topOut =
//       SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
//   auto topOut2 =
//       SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out2"));
//   auto topIn1 =
//       SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("In1"));
//   auto topIn2 =
//       SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("In2"));
//   NLLibraryTruthTables::construct(library);
//   // 7. create a and model
//   SNLDesign* andModel =
//       SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("AND"));

//   // add 2 inputs and 1 output to and
//   auto andIn1 =
//       SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in1"));
//   auto andIn2 =
//       SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("in2"));
//   auto andOut = SNLScalarTerm::create(andModel, SNLTerm::Direction::Output,
//                                       NLName("out"));
//   // 8. create an or model
//     SNLDesign* orModel =
//         SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("OR"));
//     // add 2 inputs and 1 output to or

//     auto orIn1 =
//         SNLScalarTerm::create(orModel, SNLTerm::Direction::Input, NLName("in1"));
//     auto orIn2 =
//         SNLScalarTerm::create(orModel, SNLTerm::Direction::Input, NLName("in2"));
//     auto orOut = SNLScalarTerm::create(orModel, SNLTerm::Direction::Output,
//                                        NLName("out"));

//   // Create an FF
//   SNLDesign* ffModel =
//       SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("FF"));
//   // add D, CLK and Q
//   auto ffD =
//       SNLScalarTerm::create(ffModel, SNLTerm::Direction::Input, NLName("D"));
//   auto ffCLK =
//       SNLScalarTerm::create(ffModel, SNLTerm::Direction::Input, NLName("CLK"));
//   auto ffQ =
//       SNLScalarTerm::create(ffModel, SNLTerm::Direction::Output, NLName("Q"));
//   // Set sequential dependecies to CLK
//   SNLDesignModeling::addInputsToClockArcs({ffD}, {ffCLK});
//   SNLDesignModeling::addClockToOutputsArcs({ffCLK}, {ffQ});

//   // Create ff instance under top
//   SNLInstance* instFF = SNLInstance::create(top, ffModel, NLName("ff0"));

//   // 8. create a and instance in top
//   SNLInstance* inst3 = SNLInstance::create(top, andModel, NLName("and"));
//   SNLInstance* inst4 = SNLInstance::create(top, andModel, NLName("and2"));
//   // set truth table for and model
//   SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8));
//   SNLDesignModeling::setTruthTable(orModel, SNLTruthTable(2, 14));
//   // 9. connect all instances inputs
//   SNLNet* net1 = SNLScalarNet::create(top, NLName("top_in1_net"));
//   SNLNet* net2 = SNLScalarNet::create(top, NLName("top_in2_net"));
//   SNLNet* net3 = SNLScalarNet::create(top, NLName("and_output_net"));
//   SNLNet* net4 = SNLScalarNet::create(top, NLName("and2_output_net"));
//   SNLNet* net5 = SNLScalarNet::create(top, NLName("ffD"));
//   SNLNet* net6 = SNLScalarNet::create(top, NLName("ffCLK"));
//   // connect logic0 to and
//   topIn1->setNet(net1);

//   inst4->getInstTerm(andIn1)->setNet(net2);
//   inst4->getInstTerm(andIn2)->setNet(net2);
//   // connect logic1 to and
//   instFF->getInstTerm(ffQ)->setNet(net2);
//   instFF->getInstTerm(ffD)->setNet(net1);
//   instFF->getInstTerm(ffCLK)->setNet(net6);
//   inst3->getInstTerm(andIn2)->setNet(net1);
//   inst3->getInstTerm(andIn1)->setNet(net4);
//   // connect the and instance output to the top output
//   inst3->getInstTerm(andOut)->setNet(net3);
//   topOut->setNet(net3);
//   inst4->getInstTerm(andOut)->setNet(net4);
//   topOut2->setNet(net4);
//   topIn1->setNet(net1);
//   topIn2->setNet(net6);
//   // 11. create DNL
//   get();
//   // 12. create a constant propagation object
//   {
//     std::string dotFileName(
//         std::string(std::string("./beforeCP") + std::string(".dot")));
//     std::string svgFileName(
//         std::string(std::string("./beforeCP") + std::string(".svg")));
//     SnlVisualiser snl(top);
//     snl.process();
//     snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
//     executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
//                                std::string(" -o ") + svgFileName)
//                        .c_str());
//   }
//   // Copy current top to back it up through the clone api
//   SNLDesign* topClone = top->clone(NLName("topClone"));
//   ConstantPropagation cp;
//   // 13. collect the constants
//   // cp.collectConstants();
//   // 14. run the constant propagation
//   {
//     BuildPrimaryOutputClauses miter;
//     miter.build();
//     for (const auto& po : miter.getPOs()) {
//       std::cout << "PO: " << po->toString() << std::endl;
//     }
//   }

//   cp.run();
//   // 15. check the output value of the top instance
//   {
//     std::string dotFileName(
//         std::string(std::string("./afterCP") + std::string(".dot")));
//     std::string svgFileName(
//         std::string(std::string("./afterCP") + std::string(".svg")));
//     SnlVisualiser snl(top);
//     snl.process();
//     snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
//     executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
//                                std::string(" -o ") + svgFileName)
//                        .c_str());
//   }
//   {
//     BuildPrimaryOutputClauses miter;
//     miter.build();
//     for (const auto& po : miter.getPOs()) {
//       std::cout << "PO: " << po->toString() << std::endl;
//     }
//     printf("%s\n", miter.getPOs()[0]->toString().c_str());
//     EXPECT_TRUE(miter.getPOs()[0]->toString() == std::string("((6 ∧ 6) ∧ 2)"));
//     EXPECT_TRUE(miter.getPOs()[1]->toString() == std::string("(6 ∧ 6)"));
//     EXPECT_TRUE(miter.getPOs()[2]->toString() == std::string("2"));
//     EXPECT_TRUE(miter.getPOs()[3]->toString() == std::string("3"));
//   }
//   naja::DNL::destroy();
  
//   {
//     MiterStrategy MiterS(top, topClone, "CaseA");
//     EXPECT_TRUE(MiterS.run());
//   }
//   // changing the or model to see if the miter strategy locate it
//   inst3->setModel(orModel);
//   {
//     BuildPrimaryOutputClauses miter;
//     miter.build();
//     for (const auto& po : miter.getPOs()) {
//       std::cout << "PO: " << po->toString() << std::endl;
//     }
//     MiterStrategy MiterS(top, topClone, "CaseB");
//     EXPECT_FALSE(MiterS.run());
//   }
//   // Use the PrimaryOutputClauses to get all outputs
//         BuildPrimaryOutputClauses miter;
//         miter.build();
//         for (const auto& po : miter.getOutputs()) {
//           // Get DNLFullTerminals of po
//             naja::DNL::DNLTerminalFull dnlFullTerminal = naja::DNL::get()->getDNLTerminalFromID(po);
//             std::vector<std::string> path;
//             DNLInstanceFull currentInstance = dnlFullTerminal.getDNLInstance();
//             while (currentInstance.isTop() == false) {
//             path.push_back(currentInstance.getSNLInstance()->getName().getString());
//             currentInstance = currentInstance.getParentInstance();
//             }
//             std::reverse(path.begin(), path.end());
//             naja::NL::SNLPath snlPath(top, 
//                 path);
//             SNLNetComponentOccurrence occurrence(snlPath, dnlFullTerminal.getSnlBitTerm());
//             std::cout << "Output DNL ID: " << po << ", SNL Bit Term: "
//                         << occurrence.getString() << std::endl;
//             SNLEquipotential equipotential(occurrence);
//                           std::cout << "Equipotential: " << equipotential.getString() << std::endl;
//             std::string dotFileName(
//             std::string(std::string("./") + dnlFullTerminal.getSnlBitTerm()->getName().getString() + std::string(".dot")));
//             std::string svgFileName(
//                 std::string(std::string("./") + dnlFullTerminal.getSnlBitTerm()->getName().getString() + std::string(".svg")));
//             SnlVisualiser snl(top, true, &equipotential);
//             snl.process();
//             snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
//             executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
//                                     std::string(" -o ") + svgFileName)
//                             .c_str());

//         }
// }

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
  // {
  //   // dump top to naja_if(CapProto)
  //   std::filesystem::path outputPath("./topEdited1.capnp");
  //   SNLCapnP::dump(db, outputPath);
  // }
  // Check output of binary ../../../src/bin/kepler_formal on the 2 capnp files
  // executeCommand(
  //     std::string("../../../src/bin/kepler_formal ./top.capnp ./topEdited1.capnp")
  //         .c_str());  
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
  // {
  //   // dump top to naja_if(CapProto)
  //   std::filesystem::path outputPath("./topEdited2.capnp");
  //   SNLCapnP::dump(db, outputPath);
  // }

  // Check output of binary ../../../src/bin/kepler_formal on the 2 capnp files
  // executeCommand(
  //     std::string("../../../src/bin/kepler_formal ./top.capnp ./topEdited2.capnp")
  //         .c_str());  
}

// Required main function for Google Test
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}