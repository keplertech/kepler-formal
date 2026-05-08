// Copyright 2024-2025 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>
#include <chrono>
#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>

#include "gtest/gtest.h"

#include "BuildPrimaryOutputClauses.h"
#include "ConstantPropagation.h"
#include "MiterStrategy.h"
#include "SNLLogicCloud.h"
#include "NajaDumpableProperty.h"
#include "NLLibraryTruthTables.h"
#include "NLDB0.h"
#include "NLUniverse.h"
#include "NetlistGraph.h"
#include "../../../src/sat/SATSolverWrapper.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLDesignModeling.h"
#include "SNLBusTerm.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "SNLPath.h"
#include "SNLCapnP.h"
#include "DNL.h"
#include "Tree2BoolExpr.h"

#include "Config.h"

using namespace naja;
using namespace naja::NL;
using namespace naja::NAJA_OPT;
using namespace KEPLER_FORMAL;

// Path to the kepler-formal CLI binary used by the project tests.
// Bazel sets KEPLER_BIN env var via the test's BUILD.bazel; CMake uses
// the default relative path from the build directory.
static std::filesystem::path repoRoot() {
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path().parent_path();
}

static std::string get_kepler_bin() {
  const char* env = std::getenv("KEPLER_BIN");
  if (env && *env) {
    std::filesystem::path envPath(env);
    if (envPath.is_absolute()) {
      return envPath.string();
    }
    return std::filesystem::absolute(envPath).string();
  }
  const auto root = repoRoot();
  const std::vector<std::filesystem::path> candidates = {
      root / "build/src/bin/kepler-formal",
      root / "buildR/src/bin/kepler-formal",
      root / "buildD/src/bin/kepler-formal",
      root / "src/bin/kepler-formal",
  };
  for (const auto& candidate : candidates) {
    if (std::filesystem::exists(candidate)) {
      return candidate.string();
    }
  }
  return (root / "build/src/bin/kepler-formal").string();
}
static const std::string KEPLER_BIN_STR = get_kepler_bin();
static const char* KEPLER_BIN = KEPLER_BIN_STR.c_str();

// Prefix for test data paths. Bazel sets TEST_DATA_PREFIX="" (files are
// in runfiles at workspace-relative paths); CMake uses "../../../../".
static std::string get_test_data_prefix() {
  const char* env = std::getenv("TEST_DATA_PREFIX");
  if (env) {
    return env;
  }
  return (repoRoot().string() + "/");
}

namespace {

std::string sanitizePathComponent(std::string value) {
  for (char& ch : value) {
    if (!(std::isalnum(static_cast<unsigned char>(ch)) || ch == '_' ||
          ch == '-' || ch == '.')) {
      ch = '_';
    }
  }
  if (value.empty()) {
    value = "test";
  }
  return value;
}

std::filesystem::path makeUniqueTestTempDir() {
  static size_t counter = 0;
  const auto* testInfo = ::testing::UnitTest::GetInstance()->current_test_info();
  std::ostringstream name;
  name << "kepler_formal_";
  if (testInfo) {
    name << sanitizePathComponent(testInfo->test_suite_name()) << "_"
         << sanitizePathComponent(testInfo->name()) << "_";
  }
  name << std::chrono::steady_clock::now().time_since_epoch().count() << "_"
       << counter++;
  auto dir = std::filesystem::temp_directory_path() / name.str();
  std::filesystem::create_directories(dir);
  return dir;
}

class ScopedCurrentPath {
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path)
      : original_(std::filesystem::current_path()) {
    std::filesystem::current_path(path);
  }
  ~ScopedCurrentPath() { std::filesystem::current_path(original_); }

 private:
  std::filesystem::path original_;
};

std::vector<uint64_t> getInputFlatDependencies(const SNLDesign* design) {
  std::vector<size_t> deps;
  size_t flatPos = 0;
  for (auto term : design->getBitTerms()) {
    if (term->getDirection() != SNLTerm::Direction::Output) {
      deps.push_back(flatPos);
    }
    ++flatPos;
  }
  return NLBitDependencies::encodeBits(deps);
}

void executeCommand(const std::string& command) {
  int result = system(command.c_str());
  if (result != 0) {
    std::cerr << "Command execution failed." << std::endl;
  }
}

std::string readTextFile(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

naja::DNL::DNLID findDNLTermIDByInstanceAndTerm(const char* instanceName,
                                                const char* termName) {
  auto* dnl = naja::DNL::get();
  for (naja::DNL::DNLID id = 0; id <= dnl->getNBterms(); ++id) {
    const auto& term = dnl->getDNLTerminalFromID(id);
    if (term.isNull()) {
      continue;
    }
    auto* instance = term.getDNLInstance().getSNLInstance();
    if (instance == nullptr) {
      continue;
    }
    if (instance->getName().getString() == instanceName &&
        term.getSnlBitTerm()->getName().getString() == termName) {
      return id;
    }
  }
  return naja::DNL::DNLID_MAX;
}

std::vector<std::string> getTermLabels(const std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>& ids) {
  std::vector<std::string> labels;
  labels.reserve(ids.size());
  for (auto id : ids) {
    const auto& term = naja::DNL::get()->getDNLTerminalFromID(id);
    std::string label;
    if (auto* instance = term.getDNLInstance().getSNLInstance()) {
      label += instance->getName().getString();
      label += ".";
    }
    label += term.getSnlBitTerm()->getName().getString();
    labels.push_back(label);
  }
  std::sort(labels.begin(), labels.end());
  return labels;
}

size_t countSubstringOccurrences(const std::string& text,
                                 const std::string& needle) {
  if (needle.empty()) {
    return 0;
  }
  size_t count = 0;
  size_t pos = 0;
  while ((pos = text.find(needle, pos)) != std::string::npos) {
    ++count;
    pos += needle.size();
  }
  return count;
}

class ScopedReportSkippedPOs {
 public:
  explicit ScopedReportSkippedPOs(bool enabled)
      : original_(Config::getReportSkippedPOs()) {
    Config::setReportSkippedPOs(enabled);
  }

  ~ScopedReportSkippedPOs() { Config::setReportSkippedPOs(original_); }

 private:
  bool original_;
};

void expectGenericGateMiterEquivalent(const char* gateName,
                                      SNLTruthTable::GenericType genericType) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));

  NLDB0::GateType gateType = NLDB0::GateType::And;
  switch (genericType) {
    case SNLTruthTable::GenericType::NAND:
      gateType = NLDB0::GateType::Nand;
      break;
    case SNLTruthTable::GenericType::NOR:
      gateType = NLDB0::GateType::Nor;
      break;
    case SNLTruthTable::GenericType::XNOR:
      gateType = NLDB0::GateType::Xnor;
      break;
    default:
      throw std::runtime_error("Unsupported generic type for test");
  }

  auto gateModel = NLDB0::getOrCreateNInputGate(gateType, 2);
  ASSERT_NE(gateModel, nullptr);
  SNLBitTerm* gateOut = nullptr;
  std::vector<SNLBitTerm*> gateInputs;
  for (auto term : gateModel->getBitTerms()) {
    if (term->getDirection() == SNLTerm::Direction::Input) {
      gateInputs.push_back(term);
    } else {
      gateOut = term;
    }
  }
  ASSERT_NE(gateOut, nullptr);
  ASSERT_EQ(gateInputs.size(), 2u);

  auto buildTop = [&](const char* topName) {
    auto top =
        SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(topName));
    auto topIn0 =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
    auto topIn1 =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("b"));
    auto topOut =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"));
    auto inst = SNLInstance::create(top, gateModel, NLName("gate0"));

    auto netIn0 = SNLScalarNet::create(top, NLName("net_a"));
    auto netIn1 = SNLScalarNet::create(top, NLName("net_b"));
    auto netOut = SNLScalarNet::create(top, NLName("net_y"));

    topIn0->setNet(netIn0);
    topIn1->setNet(netIn1);
    topOut->setNet(netOut);
    inst->getInstTerm(gateInputs[0])->setNet(netIn0);
    inst->getInstTerm(gateInputs[1])->setNet(netIn1);
    inst->getInstTerm(gateOut)->setNet(netOut);
    return top;
  };

  auto top0 = buildTop("top0");
  auto top1 = buildTop("top1");
  KEPLER_FORMAL::MiterStrategy miterS(top0, top1);
  miterS.init();
  EXPECT_TRUE(miterS.run());
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
    tempDir_ = makeUniqueTestTempDir();
  }
  void TearDown() override {
    // Code here will be called immediately after each test (right
    // before the destructor).
    // Destroy the SNL
    naja::DNL::destroy();
    NLUniverse::get()->destroy();
    KEPLER_FORMAL::BoolExprCache::destroy();
    if (!tempDir_.empty()) {
      std::error_code ec;
      std::filesystem::remove_all(tempDir_, ec);
    }
  }

  std::filesystem::path testTempPath(const std::string& leaf) const {
    return tempDir_ / leaf;
  }

  std::filesystem::path tempDir_;
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
  SNLDesignModeling::setTruthTable(logic0, SNLTruthTable(0, 0, SNLTruthTable::fullDependencies(0)));
  SNLDesignModeling::setTruthTable(logic1, SNLTruthTable(0, 1, SNLTruthTable::fullDependencies(0)));
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
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8, SNLTruthTable::fullDependencies(2)));
  // 9. connect all instances inputs
  SNLNet* net1 = SNLScalarNet::create(top, NLName("logic_0_net"));
  //net1->setType(SNLNet::Type::Assign0);
  SNLNet* net2 = SNLScalarNet::create(top, NLName("logic_1_net"));
  //net2->setType(SNLNet::Type::Assign1);
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
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8, SNLTruthTable::fullDependencies(2)));

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

TEST_F(MiterTests, TestGenericNandTruthTable) {
  expectGenericGateMiterEquivalent("NAND_GENERIC", SNLTruthTable::GenericType::NAND);
}

TEST_F(MiterTests, TestGenericNorTruthTable) {
  expectGenericGateMiterEquivalent("NOR_GENERIC", SNLTruthTable::GenericType::NOR);
}

TEST_F(MiterTests, TestGenericXnorTruthTable) {
  expectGenericGateMiterEquivalent("XNOR_GENERIC", SNLTruthTable::GenericType::XNOR);
}

TEST_F(MiterTests, BuildPrimaryOutputClausesConstantTrueOutput) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);

  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  SNLDesign* logic1 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC1"));
  auto logic1Out =
      SNLScalarTerm::create(logic1, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(logic1, SNLTruthTable(0, 1, SNLTruthTable::fullDependencies(0)));
  NLLibraryTruthTables::construct(library);

  SNLInstance* inst = SNLInstance::create(top, logic1, NLName("const1"));
  SNLNet* net = SNLScalarNet::create(top, NLName("const1_net"));
  inst->getInstTerm(logic1Out)->setNet(net);
  topOut->setNet(net);

  naja::DNL::get();
  BuildPrimaryOutputClauses builder;
  builder.collect();
  builder.build();

  ASSERT_EQ(builder.getPOs().size(), 1u);
  ASSERT_NE(builder.getPOs()[0], nullptr);
  EXPECT_EQ(builder.getPOs()[0]->toString(), "1");
}

TEST_F(MiterTests, BuildPrimaryOutputClausesUsesFlatDependencyCoordinatesForPOs) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);

  auto topIn0 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
  auto topIn1 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("b"));
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"));

  SNLDesign* flatDepsModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("FLAT_DEPS_AND"));
  auto flatDepsOut = SNLScalarTerm::create(
      flatDepsModel, SNLTerm::Direction::Output, NLName("y"));
  auto flatDepsIn0 = SNLScalarTerm::create(
      flatDepsModel, SNLTerm::Direction::Input, NLName("a"));
  auto flatDepsIn1 = SNLScalarTerm::create(
      flatDepsModel, SNLTerm::Direction::Input, NLName("b"));
  SNLDesignModeling::setTruthTable(
      flatDepsModel, SNLTruthTable(2, 0b1000, getInputFlatDependencies(flatDepsModel)));

  auto inst = SNLInstance::create(top, flatDepsModel, NLName("u0"));
  auto netA = SNLScalarNet::create(top, NLName("net_a"));
  auto netB = SNLScalarNet::create(top, NLName("net_b"));
  auto netY = SNLScalarNet::create(top, NLName("net_y"));

  topIn0->setNet(netA);
  topIn1->setNet(netB);
  topOut->setNet(netY);
  inst->getInstTerm(flatDepsIn0)->setNet(netA);
  inst->getInstTerm(flatDepsIn1)->setNet(netB);
  inst->getInstTerm(flatDepsOut)->setNet(netY);

  naja::DNL::get();
  BuildPrimaryOutputClauses builder;
  builder.collect();

  ASSERT_EQ(builder.getOutputs().size(), 1u);
  const auto outputTerm =
      naja::DNL::get()->getDNLTerminalFromID(builder.getOutputs()[0]);
  EXPECT_TRUE(outputTerm.isTopPort());
  EXPECT_EQ(outputTerm.getSnlBitTerm()->getName().getString(), "y");
}

TEST_F(MiterTests, BuildPrimaryOutputClausesReportsSkippedNoDriverPO) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);

  auto topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"));

  SNLDesign* passModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("PASS0"));
  auto passIn = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("data"));
  auto unusedIn0 = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("unused0"));
  auto unusedIn1 = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("unused1"));
  auto passOut = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Output, NLName("y"));
  SNLDesignModeling::setTruthTable(
      passModel,
      SNLTruthTable(3, 0xF0,
                    NLBitDependencies::encodeBits(std::vector<size_t>{0})));

  auto inst = SNLInstance::create(top, passModel, NLName("u0"));
  auto netA = SNLScalarNet::create(top, NLName("net_a"));
  auto netY = SNLScalarNet::create(top, NLName("net_y"));
  auto floatingNet0 = SNLScalarNet::create(top, NLName("floating_net_0"));
  auto floatingNet1 = SNLScalarNet::create(top, NLName("floating_net_1"));
  auto* floatingPropA =
      naja::NajaDumpableProperty::create(floatingNet0, "report_prop_a");
  floatingPropA->addStringValue("alpha");
  auto* floatingPropB =
      naja::NajaDumpableProperty::create(floatingNet0, "report_prop_b");
  floatingPropB->addUInt64Value(7);

  topIn->setNet(netA);
  topOut->setNet(netY);
  inst->getInstTerm(passIn)->setNet(netA);
  inst->getInstTerm(unusedIn0)->setNet(floatingNet0);
  inst->getInstTerm(unusedIn1)->setNet(floatingNet1);
  inst->getInstTerm(passOut)->setNet(netY);

  ScopedCurrentPath scopedCurrentPath(tempDir_);
  ScopedReportSkippedPOs reportGuard(true);
  BuildPrimaryOutputClauses builder;
  builder.collect();
  builder.collect();

  const auto reportPath = tempDir_ / "skipped_no_driver_pos.txt";
  ASSERT_TRUE(std::filesystem::exists(reportPath));
  const auto content = readTextFile(reportPath);
  EXPECT_NE(content.find("its iso has no drivers"), std::string::npos);
  EXPECT_NE(content.find("drivers: []"), std::string::npos);
  EXPECT_NE(content.find("complex_nets"), std::string::npos);
  EXPECT_NE(content.find("floating_net_0"), std::string::npos);
  EXPECT_NE(content.find("floating_net_1"), std::string::npos);
  EXPECT_NE(content.find("report_prop_a="), std::string::npos);
  EXPECT_NE(content.find("report_prop_b="), std::string::npos);
  EXPECT_NE(content.find("See first encounter of iso="), std::string::npos);
  EXPECT_EQ(countSubstringOccurrences(content, "Skipping PO "), 4u);
}

TEST_F(MiterTests, BuildPrimaryOutputClausesReportsSkippedMultiDriverPO) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);

  auto topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"));

  SNLDesign* logic0 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC0"));
  auto logic0Out =
      SNLScalarTerm::create(logic0, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(
      logic0, SNLTruthTable(0, 0, SNLTruthTable::fullDependencies(0)));

  SNLDesign* logic1 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC1"));
  auto logic1Out =
      SNLScalarTerm::create(logic1, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(
      logic1, SNLTruthTable(0, 1, SNLTruthTable::fullDependencies(0)));

  SNLDesign* passModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("PASS1"));
  auto passIn = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("data"));
  auto unusedIn = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("unused"));
  auto passOut = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Output, NLName("y"));
  SNLDesignModeling::setTruthTable(
      passModel,
      SNLTruthTable(2, 0b1100,
                    NLBitDependencies::encodeBits(std::vector<size_t>{0})));

  auto const0 = SNLInstance::create(top, logic0, NLName("const0"));
  auto const1 = SNLInstance::create(top, logic1, NLName("const1"));
  auto inst = SNLInstance::create(top, passModel, NLName("u0"));
  auto netA = SNLScalarNet::create(top, NLName("net_a"));
  auto netY = SNLScalarNet::create(top, NLName("net_y"));
  auto sharedNet = SNLScalarNet::create(top, NLName("shared_net"));
  auto* sharedPropA =
      naja::NajaDumpableProperty::create(sharedNet, "report_prop_a");
  sharedPropA->addStringValue("beta");
  auto* sharedPropB =
      naja::NajaDumpableProperty::create(sharedNet, "report_prop_b");
  sharedPropB->addUInt64Value(11);

  topIn->setNet(netA);
  topOut->setNet(netY);
  const0->getInstTerm(logic0Out)->setNet(sharedNet);
  const1->getInstTerm(logic1Out)->setNet(sharedNet);
  inst->getInstTerm(passIn)->setNet(netA);
  inst->getInstTerm(unusedIn)->setNet(sharedNet);
  inst->getInstTerm(passOut)->setNet(netY);

  ScopedCurrentPath scopedCurrentPath(tempDir_);
  ScopedReportSkippedPOs reportGuard(true);
  BuildPrimaryOutputClauses builder;
  builder.collect();
  builder.collect();

  const auto reportPath = tempDir_ / "skipped_multi_driver_pos.txt";
  ASSERT_TRUE(std::filesystem::exists(reportPath));
  const auto content = readTextFile(reportPath);
  EXPECT_NE(content.find("its iso has multiple drivers"), std::string::npos);
  EXPECT_NE(content.find("LOGIC0"), std::string::npos);
  EXPECT_NE(content.find("LOGIC1"), std::string::npos);
  EXPECT_EQ(content.find("complex_nets"), std::string::npos);
  EXPECT_EQ(content.find("report_prop_a="), std::string::npos);
  EXPECT_EQ(content.find("report_prop_b="), std::string::npos);
  EXPECT_NE(content.find("See first encounter of iso="), std::string::npos);
  EXPECT_EQ(countSubstringOccurrences(content, "Skipping PO "), 2u);
}

TEST_F(MiterTests, BuildPrimaryOutputClausesInitializesSkippedPOReportFilesOnlyOnce) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);

  auto topInA =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
  auto topInB =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("b"));
  auto topOut0 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y0"));
  auto topOut1 =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y1"));

  SNLDesign* logic0 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC0"));
  auto logic0Out =
      SNLScalarTerm::create(logic0, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(
      logic0, SNLTruthTable(0, 0, SNLTruthTable::fullDependencies(0)));

  SNLDesign* logic1 =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("LOGIC1"));
  auto logic1Out =
      SNLScalarTerm::create(logic1, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(
      logic1, SNLTruthTable(0, 1, SNLTruthTable::fullDependencies(0)));

  SNLDesign* passModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("PASS1"));
  auto passIn = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("data"));
  auto unusedIn = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("unused"));
  auto passOut = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Output, NLName("y"));
  SNLDesignModeling::setTruthTable(
      passModel,
      SNLTruthTable(2, 0b1100,
                    NLBitDependencies::encodeBits(std::vector<size_t>{0})));

  auto const0 = SNLInstance::create(top, logic0, NLName("const0"));
  auto const1 = SNLInstance::create(top, logic1, NLName("const1"));
  auto noDriverInst = SNLInstance::create(top, passModel, NLName("u_no_driver"));
  auto multiDriverInst =
      SNLInstance::create(top, passModel, NLName("u_multi_driver"));

  auto netA = SNLScalarNet::create(top, NLName("net_a"));
  auto netB = SNLScalarNet::create(top, NLName("net_b"));
  auto netY0 = SNLScalarNet::create(top, NLName("net_y0"));
  auto netY1 = SNLScalarNet::create(top, NLName("net_y1"));
  auto floatingNet = SNLScalarNet::create(top, NLName("floating_net"));
  auto sharedNet = SNLScalarNet::create(top, NLName("shared_net"));

  topInA->setNet(netA);
  topInB->setNet(netB);
  topOut0->setNet(netY0);
  topOut1->setNet(netY1);

  const0->getInstTerm(logic0Out)->setNet(sharedNet);
  const1->getInstTerm(logic1Out)->setNet(sharedNet);

  noDriverInst->getInstTerm(passIn)->setNet(netA);
  noDriverInst->getInstTerm(unusedIn)->setNet(floatingNet);
  noDriverInst->getInstTerm(passOut)->setNet(netY0);

  multiDriverInst->getInstTerm(passIn)->setNet(netB);
  multiDriverInst->getInstTerm(unusedIn)->setNet(sharedNet);
  multiDriverInst->getInstTerm(passOut)->setNet(netY1);

  ScopedCurrentPath scopedCurrentPath(tempDir_);
  ScopedReportSkippedPOs reportGuard(true);
  BuildPrimaryOutputClauses builder;
  builder.collect();

  EXPECT_TRUE(std::filesystem::exists(tempDir_ / "skipped_no_driver_pos.txt"));
  EXPECT_TRUE(std::filesystem::exists(tempDir_ / "skipped_multi_driver_pos.txt"));
}

TEST_F(MiterTests, ReportingModePreservesCachedIsoShortcutInLogicCloud) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);

  auto topInA =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
  auto topInB =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("b"));
  auto topInC =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("c"));
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"));

  SNLDesign* andModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("AND2"));
  auto andIn0 =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("a"));
  auto andIn1 =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Input, NLName("b"));
  auto andOut =
      SNLScalarTerm::create(andModel, SNLTerm::Direction::Output, NLName("y"));
  SNLDesignModeling::setTruthTable(
      andModel, SNLTruthTable(2, 0b1000, SNLTruthTable::fullDependencies(2)));

  SNLDesign* bufModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("BUF1"));
  auto bufIn =
      SNLScalarTerm::create(bufModel, SNLTerm::Direction::Input, NLName("a"));
  auto bufOut =
      SNLScalarTerm::create(bufModel, SNLTerm::Direction::Output, NLName("y"));
  SNLDesignModeling::setTruthTable(
      bufModel, SNLTruthTable(1, 0b10, SNLTruthTable::fullDependencies(1)));

  SNLDesign* orModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("OR2"));
  auto orIn0 =
      SNLScalarTerm::create(orModel, SNLTerm::Direction::Input, NLName("a"));
  auto orIn1 =
      SNLScalarTerm::create(orModel, SNLTerm::Direction::Input, NLName("b"));
  auto orOut =
      SNLScalarTerm::create(orModel, SNLTerm::Direction::Output, NLName("y"));
  SNLDesignModeling::setTruthTable(
      orModel, SNLTruthTable(2, 0b1110, SNLTruthTable::fullDependencies(2)));

  auto andInst = SNLInstance::create(top, andModel, NLName("and0"));
  auto bufInst = SNLInstance::create(top, bufModel, NLName("buf0"));
  auto orInst = SNLInstance::create(top, orModel, NLName("or0"));

  auto netA = SNLScalarNet::create(top, NLName("net_a"));
  auto netB = SNLScalarNet::create(top, NLName("net_b"));
  auto netC = SNLScalarNet::create(top, NLName("net_c"));
  auto netAnd = SNLScalarNet::create(top, NLName("net_and"));
  auto netBuf = SNLScalarNet::create(top, NLName("net_buf"));
  auto netY = SNLScalarNet::create(top, NLName("net_y"));

  topInA->setNet(netA);
  topInB->setNet(netB);
  topInC->setNet(netC);
  topOut->setNet(netY);

  andInst->getInstTerm(andIn0)->setNet(netA);
  andInst->getInstTerm(andIn1)->setNet(netB);
  andInst->getInstTerm(andOut)->setNet(netAnd);

  bufInst->getInstTerm(bufIn)->setNet(netAnd);
  bufInst->getInstTerm(bufOut)->setNet(netBuf);

  orInst->getInstTerm(orIn0)->setNet(netBuf);
  orInst->getInstTerm(orIn1)->setNet(netC);
  orInst->getInstTerm(orOut)->setNet(netY);

  naja::DNL::destroy();
  BuildPrimaryOutputClauses builder;
  builder.collect();

  auto* dnl = naja::DNL::get();
  ASSERT_EQ(builder.getOutputs().size(), 1u);

  std::vector<bool> isPIs(dnl->getNBterms() + 1, false);
  for (auto input : builder.getInputs()) {
    isPIs[input] = true;
  }

  std::vector<bool> isPOs(dnl->getNBterms() + 1, false);
  for (auto output : builder.getOutputs()) {
    isPOs[output] = true;
  }

  const auto andOutID = findDNLTermIDByInstanceAndTerm("and0", "y");
  ASSERT_NE(andOutID, naja::DNL::DNLID_MAX);
  const auto andIsoID = dnl->getDNLTerminalFromID(andOutID).getIsoID();
  ASSERT_NE(andIsoID, naja::DNL::DNLID_MAX);

  Tree2BoolExpr::iso2boolExpr_.clear();
  Tree2BoolExpr::iso2boolExpr_.insert({andIsoID, BoolExpr::Var(999)});

  {
    ScopedReportSkippedPOs reportGuard(false);
    SNLLogicCloud cloud(builder.getOutputs()[0], isPIs, isPOs);
    cloud.compute();
    EXPECT_EQ(getTermLabels(cloud.getInputs()),
              (std::vector<std::string>{"and0.y", "c"}));
    cloud.destroy();
  }

  {
    ScopedReportSkippedPOs reportGuard(true);
    SNLLogicCloud cloud(builder.getOutputs()[0], isPIs, isPOs);
    cloud.compute();
    EXPECT_EQ(getTermLabels(cloud.getInputs()),
              (std::vector<std::string>{"and0.y", "c"}));
    cloud.destroy();
  }

  Tree2BoolExpr::iso2boolExpr_.clear();
  naja::DNL::destroy();
}

TEST_F(MiterTests, MiterStrategySummaryUsesUnnamedFallbackLabels) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryS =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("standard"));

  SNLDesign* top =
      SNLDesign::create(libraryS, SNLDesign::Type::Standard, NLName(""));
  univ->setTopDesign(top);
  auto topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName(""));
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName(""));
  auto net = SNLScalarNet::create(top, NLName("n0"));
  topIn->setNet(net);
  topOut->setNet(net);

  SNLDesign* topClone = top->clone(NLName(""));

  ScopedCurrentPath scopedCurrentPath(tempDir_);
  const auto logPath = tempDir_ / "unnamed_summary.txt";
  MiterStrategy strategy(top, topClone, logPath.string());
  strategy.init();

  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto content = readTextFile(logPath);
  EXPECT_NE(content.find("<unnamed>"), std::string::npos);
}

TEST_F(MiterTests, SNLLogicCloudReportsSkippedNoDriverRoot) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  SNLDesign* top =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("top"));
  univ->setTopDesign(top);

  auto topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"));

  SNLDesign* passModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("PASSCLOUD"));
  auto passIn = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("data"));
  auto unusedIn = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Input, NLName("unused"));
  auto passOut = SNLScalarTerm::create(
      passModel, SNLTerm::Direction::Output, NLName("y"));
  SNLDesignModeling::setTruthTable(
      passModel,
      SNLTruthTable(2, 0b1100,
                    NLBitDependencies::encodeBits(std::vector<size_t>{0})));

  auto inst = SNLInstance::create(top, passModel, NLName("u0"));
  auto netA = SNLScalarNet::create(top, NLName("net_a"));
  auto netY = SNLScalarNet::create(top, NLName("net_y"));
  auto floatingNet = SNLScalarNet::create(top, NLName("floating_net"));

  topIn->setNet(netA);
  topOut->setNet(netY);
  inst->getInstTerm(passIn)->setNet(netA);
  inst->getInstTerm(unusedIn)->setNet(floatingNet);
  inst->getInstTerm(passOut)->setNet(netY);

  ScopedCurrentPath scopedCurrentPath(tempDir_);
  ScopedReportSkippedPOs reportGuard(true);
  naja::DNL::destroy();
  BuildPrimaryOutputClauses builder;
  builder.collect();
  auto* dnl = naja::DNL::get();
  ASSERT_EQ(builder.getOutputs().size(), 1u);

  std::vector<bool> isPIs(dnl->getNBterms() + 1, false);
  for (auto input : builder.getInputs()) {
    isPIs[input] = true;
  }

  std::vector<bool> isPOs(dnl->getNBterms() + 1, false);
  for (auto output : builder.getOutputs()) {
    isPOs[output] = true;
  }

  SNLLogicCloud cloud(builder.getOutputs()[0], isPIs, isPOs);
  cloud.compute();
  EXPECT_FALSE(cloud.getTruthTable().isValid());
  SNLLogicCloud::flushSkippedPOReports();

  const auto reportPath = tempDir_ / "skipped_no_driver_pos.txt";
  ASSERT_TRUE(std::filesystem::exists(reportPath));
  const auto content = readTextFile(reportPath);
  EXPECT_NE(content.find("no drivers during cloud expansion"),
            std::string::npos);
  EXPECT_NE(content.find("current_input="), std::string::npos);
  EXPECT_NE(content.find("drivers: []"), std::string::npos);
}

TEST(MiterStandaloneTests, KissatClauseAutoExpandsTrackedVariableCount) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::KISSAT);
  solver.addClause({5});
  EXPECT_EQ(solver.newVar(), 4);
}

TEST(MiterStrategyStandaloneTests, NormalizeOutputsIgnoresOutputsOnlyPresentInSecondNetlist) {
  MiterStrategy strategy(nullptr, nullptr, "normalizeOutputs");
  using PathKey = BuildPrimaryOutputClauses::PathKey;
  using OutputMap = std::unordered_map<PathKey, naja::DNL::DNLID,
                                       BuildPrimaryOutputClauses::KeyHash>;

  const PathKey common{{1}, {10}};
  const PathKey onlySecond{{2}, {20}};

  OutputMap outputs0Map;
  OutputMap outputs1Map;
  outputs0Map.emplace(common, 100);
  outputs1Map.emplace(common, 200);
  outputs1Map.emplace(onlySecond, 300);

  std::vector<naja::DNL::DNLID> outputs0{100};
  std::vector<naja::DNL::DNLID> outputs1{200, 300};

  strategy.normalizeOutputs(outputs0, outputs1, outputs0Map, outputs1Map);

  EXPECT_EQ(outputs0, std::vector<naja::DNL::DNLID>({100}));
  EXPECT_EQ(outputs1, std::vector<naja::DNL::DNLID>({200}));
}

TEST(MiterStrategyStandaloneTests, RunCompactSnapshotsAlignsInputsOutputsAndWritesCnf) {
  using PathKey = BuildPrimaryOutputClauses::PathKey;
  auto makePathKey = [](int nameID, int objectID) -> PathKey {
    return {{static_cast<NLName::ID>(nameID)},
            {static_cast<NLID::DesignObjectID>(objectID)}};
  };

  const PathKey a = makePathKey(1, 10);
  const PathKey b = makePathKey(2, 20);
  const PathKey only0 = makePathKey(3, 30);
  const PathKey only1 = makePathKey(4, 40);
  const PathKey logicOut = makePathKey(5, 50);
  const PathKey constOut = makePathKey(6, 60);
  const PathKey drop0 = makePathKey(7, 70);
  const PathKey drop1 = makePathKey(8, 80);

  MiterStrategy::CompactSnapshot snapshot0;
  snapshot0.inputs = {a, b, only0};
  snapshot0.outputs = {logicOut, constOut, drop0};
  auto shared0 = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));
  snapshot0.POs.emplace_back(
      BoolExpr::Or(BoolExpr::Not(shared0),
                   BoolExpr::Xor(shared0, BoolExpr::Var(4))));
  snapshot0.POs.emplace_back(BoolExpr::createTrue());
  snapshot0.POs.emplace_back(BoolExpr::Var(2));

  MiterStrategy::CompactSnapshot snapshot1;
  snapshot1.inputs = {b, a, only1};
  snapshot1.outputs = {constOut, logicOut, drop1};
  auto shared1 = BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(2));
  snapshot1.POs.emplace_back(BoolExpr::createTrue());
  snapshot1.POs.emplace_back(
      BoolExpr::Or(BoolExpr::Not(shared1),
                   BoolExpr::Xor(shared1, BoolExpr::Var(4))));
  snapshot1.POs.emplace_back(BoolExpr::Var(3));

  const auto tmpDir = std::filesystem::temp_directory_path() /
                      "kepler_formal_compact_snapshot_cnf";
  std::filesystem::create_directories(tmpDir);
  const auto cnfPath = tmpDir / "compact_snapshot.cnf";
  const auto poCnfDir = tmpDir / "po_cnfs";

  MiterStrategy strategy(nullptr, nullptr,
                         (tmpDir / "compact_snapshot.log").string());
  strategy.setCnfDump(true, cnfPath.string());
  strategy.setPoCnfDump(true, poCnfDir.string());

  EXPECT_TRUE(strategy.runCompactSnapshots(snapshot0, snapshot1));
  EXPECT_TRUE(std::filesystem::exists(cnfPath));
  EXPECT_TRUE(std::filesystem::exists(poCnfDir / "top0" / "po_000000.cnf"));
  EXPECT_TRUE(std::filesystem::exists(poCnfDir / "top1" / "po_000000.cnf"));

  std::filesystem::remove_all(tmpDir);
}

TEST(MiterStrategyStandaloneTests, RunCompactSnapshotsPoCnfHandlesInvalidAndWriteFailure) {
  using PathKey = BuildPrimaryOutputClauses::PathKey;
  auto makePathKey = [](int nameID, int objectID) -> PathKey {
    return {{static_cast<NLName::ID>(nameID)},
            {static_cast<NLID::DesignObjectID>(objectID)}};
  };

  BoolExpr invalid0;
  BoolExpr invalid1;

  const PathKey in = makePathKey(1, 10);
  const PathKey invalidOut = makePathKey(2, 20);
  const PathKey validOut = makePathKey(3, 30);

  MiterStrategy::CompactSnapshot snapshot0;
  snapshot0.inputs = {in};
  snapshot0.outputs = {invalidOut, validOut};
  snapshot0.POs.emplace_back(&invalid0);
  snapshot0.POs.emplace_back(BoolExpr::Var(2));

  MiterStrategy::CompactSnapshot snapshot1;
  snapshot1.inputs = {in};
  snapshot1.outputs = {invalidOut, validOut};
  snapshot1.POs.emplace_back(&invalid1);
  snapshot1.POs.emplace_back(BoolExpr::Var(2));

  const auto tmpDir = std::filesystem::temp_directory_path() /
                      "kepler_formal_compact_snapshot_po_failures";
  const auto poCnfDir = tmpDir / "po_cnfs";
  std::filesystem::create_directories(poCnfDir / "top0" / "po_000001.cnf");
  std::filesystem::create_directories(poCnfDir / "top1" / "po_000001.cnf");

  MiterStrategy strategy(nullptr, nullptr,
                         (tmpDir / "compact_snapshot_po_failures.log").string());
  strategy.setPoCnfDump(true, poCnfDir.string());

  EXPECT_TRUE(strategy.runCompactSnapshots(snapshot0, snapshot1));

  std::filesystem::remove_all(tmpDir);
}

TEST(MiterStrategyStandaloneTests, RunCompactSnapshotsPoCnfHandlesDirectoryCreationFailure) {
  using PathKey = BuildPrimaryOutputClauses::PathKey;
  auto makePathKey = [](int nameID, int objectID) -> PathKey {
    return {{static_cast<NLName::ID>(nameID)},
            {static_cast<NLID::DesignObjectID>(objectID)}};
  };

  const PathKey in = makePathKey(1, 10);
  const PathKey out = makePathKey(2, 20);

  MiterStrategy::CompactSnapshot snapshot0;
  snapshot0.inputs = {in};
  snapshot0.outputs = {out};
  snapshot0.POs.emplace_back(BoolExpr::Var(2));

  MiterStrategy::CompactSnapshot snapshot1;
  snapshot1.inputs = {in};
  snapshot1.outputs = {out};
  snapshot1.POs.emplace_back(BoolExpr::Var(2));

  const auto tmpDir = std::filesystem::temp_directory_path() /
                      "kepler_formal_compact_snapshot_po_dir_failure";
  std::filesystem::create_directories(tmpDir);
  const auto poCnfPath = tmpDir / "po_cnfs";
  {
    std::ofstream outFile(poCnfPath);
    outFile << "not a directory\n";
  }

  MiterStrategy strategy(nullptr, nullptr,
                         (tmpDir / "compact_snapshot_po_dir_failure.log").string());
  strategy.setPoCnfDump(true, poCnfPath.string());

  EXPECT_TRUE(strategy.runCompactSnapshots(snapshot0, snapshot1));

  std::filesystem::remove_all(tmpDir);
}

TEST(MiterStrategyStandaloneTests, RunCompactSnapshotsWithNoCommonOutputsIsVacuouslyEquivalent) {
  using PathKey = BuildPrimaryOutputClauses::PathKey;
  auto makePathKey = [](int nameID, int objectID) -> PathKey {
    return {{static_cast<NLName::ID>(nameID)},
            {static_cast<NLID::DesignObjectID>(objectID)}};
  };

  MiterStrategy::CompactSnapshot snapshot0;
  snapshot0.inputs = {makePathKey(1, 10)};
  snapshot0.outputs = {makePathKey(2, 20)};
  snapshot0.POs.emplace_back(BoolExpr::Var(2));

  MiterStrategy::CompactSnapshot snapshot1;
  snapshot1.inputs = {makePathKey(1, 10)};
  snapshot1.outputs = {makePathKey(3, 30)};
  snapshot1.POs.emplace_back(BoolExpr::Var(2));

  MiterStrategy strategy(nullptr, nullptr, "compactSnapshotsNoCommonOutputs");
  EXPECT_TRUE(strategy.runCompactSnapshots(snapshot0, snapshot1));
}

TEST_F(MiterTests, CompactRunEquivalentDesignsInSeparateDBsWritesCnf) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db0 = NLDB::create(univ);
  NLDB* db1 = NLDB::create(univ);

  NLLibrary* library0 =
      NLLibrary::create(db0, NLLibrary::Type::Standard, NLName("designs0"));
  NLLibrary* library1 =
      NLLibrary::create(db1, NLLibrary::Type::Standard, NLName("designs1"));

  SNLDesign* top0 =
      SNLDesign::create(library0, SNLDesign::Type::Standard, NLName("top0"));
  SNLDesign* top1 =
      SNLDesign::create(library1, SNLDesign::Type::Standard, NLName("top1"));
  univ->setTopDesign(top0);

  auto top0In =
      SNLScalarTerm::create(top0, SNLTerm::Direction::Input, NLName("a"));
  auto top0Out =
      SNLScalarTerm::create(top0, SNLTerm::Direction::Output, NLName("y"));
  auto top1In =
      SNLScalarTerm::create(top1, SNLTerm::Direction::Input, NLName("a"));
  auto top1Out =
      SNLScalarTerm::create(top1, SNLTerm::Direction::Output, NLName("y"));

  auto net0 = SNLScalarNet::create(top0, NLName("net_a"));
  top0In->setNet(net0);
  top0Out->setNet(net0);

  auto net1 = SNLScalarNet::create(top1, NLName("net_a"));
  top1In->setNet(net1);
  top1Out->setNet(net1);

  const auto cnfPath = testTempPath("compact_run_separate_dbs.cnf");
  MiterStrategy strategy(top0, top1, testTempPath("compact_run_separate_dbs.log").string());
  strategy.init();
  strategy.setCnfDump(true, cnfPath.string());

  EXPECT_TRUE(strategy.run(true));
  EXPECT_TRUE(std::filesystem::exists(cnfPath));
}

TEST_F(MiterTests, InvertedOutputsProduceConstantTrueMiter) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  NLLibrary* designs =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));

  SNLDesign* top0 =
      SNLDesign::create(designs, SNLDesign::Type::Standard, NLName("top0"));
  SNLDesign* top1 =
      SNLDesign::create(designs, SNLDesign::Type::Standard, NLName("top1"));
  univ->setTopDesign(top0);

  auto top0In =
      SNLScalarTerm::create(top0, SNLTerm::Direction::Input, NLName("a"));
  auto top0Out =
      SNLScalarTerm::create(top0, SNLTerm::Direction::Output, NLName("y"));
  auto top1In =
      SNLScalarTerm::create(top1, SNLTerm::Direction::Input, NLName("a"));
  auto top1Out =
      SNLScalarTerm::create(top1, SNLTerm::Direction::Output, NLName("y"));

  SNLDesign* invModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("INV"));
  auto invIn =
      SNLScalarTerm::create(invModel, SNLTerm::Direction::Input, NLName("A"));
  auto invOut =
      SNLScalarTerm::create(invModel, SNLTerm::Direction::Output, NLName("ZN"));
  SNLDesignModeling::setTruthTable(invModel, SNLTruthTable(1, 1, SNLTruthTable::fullDependencies(1)));
  NLLibraryTruthTables::construct(library);

  auto top0InputNet = SNLScalarNet::create(top0, NLName("a_net"));
  top0In->setNet(top0InputNet);
  top0Out->setNet(top0InputNet);

  auto top1InputNet = SNLScalarNet::create(top1, NLName("a_net"));
  auto top1OutputNet = SNLScalarNet::create(top1, NLName("y_net"));
  top1In->setNet(top1InputNet);
  top1Out->setNet(top1OutputNet);
  SNLInstance* inv = SNLInstance::create(top1, invModel, NLName("inv0"));
  inv->getInstTerm(invIn)->setNet(top1InputNet);
  inv->getInstTerm(invOut)->setNet(top1OutputNet);

  MiterStrategy strategy(top0, top1, "constant_true_miter");
  strategy.init();
  EXPECT_FALSE(strategy.run());
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
  SNLDesignModeling::addInputsToClockArcs({ffD}, ffCLK);
  SNLDesignModeling::addClockToOutputsArcs(ffCLK, {ffQ});

  // Create ff instance under top
  SNLInstance* instFF = SNLInstance::create(top, ffModel, NLName("ff0"));

  // 8. create a and instance in top
  SNLInstance* inst3 = SNLInstance::create(top, andModel, NLName("and"));
  SNLInstance* inst4 = SNLInstance::create(top, andModel, NLName("and2"));
  // set truth table for and model
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8, SNLTruthTable::fullDependencies(2)));
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

TEST_F(MiterTests, ReducedTruthTableArityStillQueuesAllInstanceInputs) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryDesigns =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));

  SNLDesign* buf2Model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("BUF2"));
  auto bufIn1 =
      SNLScalarTerm::create(buf2Model, SNLTerm::Direction::Input, NLName("in1"));
  auto bufIn2 =
      SNLScalarTerm::create(buf2Model, SNLTerm::Direction::Input, NLName("in2"));
  auto bufOut = SNLScalarTerm::create(buf2Model, SNLTerm::Direction::Output,
                                      NLName("out"));
  SNLDesignModeling::setTruthTable(
      buf2Model,
      SNLTruthTable(1, 2, NLBitDependencies::encodeBits(std::vector<size_t>{0})));

  auto buildTop = [&](const char* topName) {
    auto top = SNLDesign::create(
        libraryDesigns, SNLDesign::Type::Standard, NLName(topName));
    univ->setTopDesign(top);

    auto topA =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
    auto topB =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("b"));
    auto topOut =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

    auto inst = SNLInstance::create(top, buf2Model, NLName("buf2"));

    auto netA = SNLScalarNet::create(top, NLName("net_a"));
    auto netB = SNLScalarNet::create(top, NLName("net_b"));
    auto netOut = SNLScalarNet::create(top, NLName("net_out"));

    topA->setNet(netA);
    topB->setNet(netB);
    topOut->setNet(netOut);

    inst->getInstTerm(bufIn1)->setNet(netA);
    inst->getInstTerm(bufIn2)->setNet(netB);
    inst->getInstTerm(bufOut)->setNet(netOut);

    return top;
  };

  auto top = buildTop("top");
  auto topClone = top->clone(NLName("topClone"));
  MiterStrategy MiterS(top, topClone, "ReducedTruthTableArity");
  MiterS.init();
  try {
    (void)MiterS.run();
    FAIL() << "Expected arity mismatch runtime_error";
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    EXPECT_NE(message.find("SNLLogicCloud arity mismatch"), std::string::npos);
    EXPECT_NE(message.find("TT arity=1"), std::string::npos);
    EXPECT_NE(message.find("model non-output term count=2"), std::string::npos);
  }
}

TEST_F(MiterTests, BuildPrimaryOutputClausesDoesNotTreatWideMuxInputsAsPOs) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryDesigns =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));

  SNLDesign* muxModel = NLDB0::getOrCreateMux2(4);
  ASSERT_NE(nullptr, muxModel);
  EXPECT_EQ(4u, SNLDesignModeling::getTruthTableCount(muxModel));

  auto top = SNLDesign::create(
      libraryDesigns, SNLDesign::Type::Standard, NLName("top"));
  univ->setTopDesign(top);

  std::array<SNLScalarTerm*, 4> topA{};
  std::array<SNLScalarTerm*, 4> topB{};
  std::array<SNLScalarTerm*, 4> topY{};
  for (size_t bit = 0; bit < 4; ++bit) {
    topA[bit] = SNLScalarTerm::create(
        top, SNLTerm::Direction::Input, NLName("a" + std::to_string(bit)));
    topB[bit] = SNLScalarTerm::create(
        top, SNLTerm::Direction::Input, NLName("b" + std::to_string(bit)));
    topY[bit] = SNLScalarTerm::create(
        top, SNLTerm::Direction::Output, NLName("y" + std::to_string(bit)));
  }
  auto topS =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("s"));

  auto muxInst = SNLInstance::create(top, muxModel, NLName("mux0"));
  auto netS = SNLScalarNet::create(top, NLName("net_s"));
  topS->setNet(netS);
  muxInst->getInstTerm(NLDB0::getMux2Select(muxModel))->setNet(netS);

  for (size_t bit = 0; bit < 4; ++bit) {
    auto* netA =
        SNLScalarNet::create(top, NLName("net_a" + std::to_string(bit)));
    auto* netB =
        SNLScalarNet::create(top, NLName("net_b" + std::to_string(bit)));
    auto* netY =
        SNLScalarNet::create(top, NLName("net_y" + std::to_string(bit)));
    topA[bit]->setNet(netA);
    topB[bit]->setNet(netB);
    topY[bit]->setNet(netY);
    muxInst->getInstTerm(NLDB0::getMux2InputA(muxModel)->getBit(bit))
        ->setNet(netA);
    muxInst->getInstTerm(NLDB0::getMux2InputB(muxModel)->getBit(bit))
        ->setNet(netB);
    muxInst->getInstTerm(NLDB0::getMux2Output(muxModel)->getBit(bit))
        ->setNet(netY);
  }

  naja::DNL::get();
  BuildPrimaryOutputClauses builder;
  builder.collect();

  ASSERT_EQ(builder.getOutputs().size(), 4u);
  std::set<std::string> outputNames;
  for (const auto outputID : builder.getOutputs()) {
    const auto& outputTerm = naja::DNL::get()->getDNLTerminalFromID(outputID);
    EXPECT_TRUE(outputTerm.isTopPort());
    outputNames.insert(outputTerm.getSnlBitTerm()->getName().getString());
  }
  EXPECT_EQ(outputNames, (std::set<std::string>{"y0", "y1", "y2", "y3"}));
}

TEST_F(MiterTests, WideMuxMiterEquivalent) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryDesigns =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));

  SNLDesign* muxModel = NLDB0::getOrCreateMux2(4);
  ASSERT_NE(nullptr, muxModel);
  EXPECT_EQ(4u, SNLDesignModeling::getTruthTableCount(muxModel));

  auto buildTop = [&](const char* topName) {
    auto top = SNLDesign::create(
        libraryDesigns, SNLDesign::Type::Standard, NLName(topName));
    univ->setTopDesign(top);

    std::array<SNLScalarTerm*, 4> topA{};
    std::array<SNLScalarTerm*, 4> topB{};
    std::array<SNLScalarTerm*, 4> topY{};
    for (size_t bit = 0; bit < 4; ++bit) {
      topA[bit] = SNLScalarTerm::create(
          top, SNLTerm::Direction::Input, NLName("a" + std::to_string(bit)));
      topB[bit] = SNLScalarTerm::create(
          top, SNLTerm::Direction::Input, NLName("b" + std::to_string(bit)));
      topY[bit] = SNLScalarTerm::create(
          top, SNLTerm::Direction::Output, NLName("y" + std::to_string(bit)));
    }
    auto topS =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("s"));

    auto muxInst = SNLInstance::create(top, muxModel, NLName("mux0"));
    auto netS = SNLScalarNet::create(top, NLName("net_s"));
    topS->setNet(netS);
    muxInst->getInstTerm(NLDB0::getMux2Select(muxModel))->setNet(netS);

    for (size_t bit = 0; bit < 4; ++bit) {
      auto* netA =
          SNLScalarNet::create(top, NLName("net_a" + std::to_string(bit)));
      auto* netB =
          SNLScalarNet::create(top, NLName("net_b" + std::to_string(bit)));
      auto* netY =
          SNLScalarNet::create(top, NLName("net_y" + std::to_string(bit)));
      topA[bit]->setNet(netA);
      topB[bit]->setNet(netB);
      topY[bit]->setNet(netY);
      muxInst->getInstTerm(NLDB0::getMux2InputA(muxModel)->getBit(bit))
          ->setNet(netA);
      muxInst->getInstTerm(NLDB0::getMux2InputB(muxModel)->getBit(bit))
          ->setNet(netB);
      muxInst->getInstTerm(NLDB0::getMux2Output(muxModel)->getBit(bit))
          ->setNet(netY);
    }

    return top;
  };

  auto top = buildTop("top");
  auto topClone = top->clone(NLName("topClone"));
  MiterStrategy MiterS(top, topClone, "WideMuxMiterEquivalent");
  MiterS.init();
  EXPECT_TRUE(MiterS.run());
}

TEST_F(MiterTests, Db0FAArityCheckFailureReproducer) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryDesigns =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));

  SNLDesign* faModel = NLDB0::getFA();
  ASSERT_NE(nullptr, faModel);

  auto buildTop = [&](const char* topName) {
    auto top = SNLDesign::create(
        libraryDesigns, SNLDesign::Type::Standard, NLName(topName));
    univ->setTopDesign(top);

    auto topA =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"));
    auto topB =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("b"));
    auto topCI =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("ci"));
    auto topS =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("s"));
    auto topCO =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("co"));

    auto inst = SNLInstance::create(top, faModel, NLName("fa0"));

    auto netA = SNLScalarNet::create(top, NLName("net_a"));
    auto netB = SNLScalarNet::create(top, NLName("net_b"));
    auto netCI = SNLScalarNet::create(top, NLName("net_ci"));
    auto netS = SNLScalarNet::create(top, NLName("net_s"));
    auto netCO = SNLScalarNet::create(top, NLName("net_co"));

    topA->setNet(netA);
    topB->setNet(netB);
    topCI->setNet(netCI);
    topS->setNet(netS);
    topCO->setNet(netCO);

    inst->getInstTerm(NLDB0::getFAInputA())->setNet(netA);
    inst->getInstTerm(NLDB0::getFAInputB())->setNet(netB);
    inst->getInstTerm(NLDB0::getFAInputCI())->setNet(netCI);
    inst->getInstTerm(NLDB0::getFAOutputS())->setNet(netS);
    inst->getInstTerm(NLDB0::getFAOutputCO())->setNet(netCO);

    return top;
  };

  auto top = buildTop("top");
  auto topClone = top->clone(NLName("topClone"));
  MiterStrategy MiterS(top, topClone, "Db0FAArityCheck");
  MiterS.init();
  EXPECT_EQ(2, SNLDesignModeling::getTruthTableCount(faModel));

  // Regression for the former workflow failure path: DB0 FA must be
  // handled through per-output truth-table counting, not the single-output
  // primitive truth-table API.
  EXPECT_TRUE(MiterS.run());
}

TEST_F(MiterTests, InternalPOAssignCycleIsSkippedAndReported) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryDesigns =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));

  auto assignModel = NLDB0::getAssign();
  ASSERT_NE(nullptr, assignModel);
  auto dffrnModel = NLDB0::getDFFRN();
  ASSERT_NE(nullptr, dffrnModel);
  auto faModel = NLDB0::getFA();
  ASSERT_NE(nullptr, faModel);
  auto muxModel = NLDB0::getMux2();
  ASSERT_NE(nullptr, muxModel);

  auto buildTop = [&](const char* topName) {
    auto top = SNLDesign::create(
        libraryDesigns, SNLDesign::Type::Standard, NLName(topName));
    univ->setTopDesign(top);

    auto topB =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("b"));
    auto topCI =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("ci"));
    auto topRN =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rn"));
    auto topC =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("c"));
    auto topMuxA =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("ma"));
    auto topMuxB =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("mb"));
    auto topMuxS =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("ms"));

    auto dffrnInst = SNLInstance::create(top, dffrnModel, NLName("ff0"));
    auto assignSeedInst =
        SNLInstance::create(top, assignModel, NLName("assign_seed"));
    auto faInst = SNLInstance::create(top, faModel, NLName("fa0"));
    auto assignSelInst =
        SNLInstance::create(top, assignModel, NLName("assign_sel"));
    auto mux0Inst = SNLInstance::create(top, muxModel, NLName("mux0"));
    auto mux1Inst = SNLInstance::create(top, muxModel, NLName("mux1"));

    auto netRN = SNLScalarNet::create(top, NLName("net_rn"));
    auto netC = SNLScalarNet::create(top, NLName("net_c"));
    auto netB = SNLScalarNet::create(top, NLName("net_b"));
    auto netCI = SNLScalarNet::create(top, NLName("net_ci"));
    auto netMuxA = SNLScalarNet::create(top, NLName("net_ma"));
    auto netMuxB = SNLScalarNet::create(top, NLName("net_mb"));
    auto netMuxS = SNLScalarNet::create(top, NLName("net_ms"));
    auto netQ = SNLScalarNet::create(top, NLName("net_q"));
    auto netSeedOut = SNLScalarNet::create(top, NLName("net_seed_out"));
    auto netSeedIn = SNLScalarNet::create(top, NLName("net_seed_in"));
    auto netFaS = SNLScalarNet::create(top, NLName("net_fa_s"));
    auto netSelOut = SNLScalarNet::create(top, NLName("net_sel_out"));
    auto netSelIn = SNLScalarNet::create(top, NLName("net_sel_in"));

    topB->setNet(netB);
    topCI->setNet(netCI);
    topRN->setNet(netRN);
    topC->setNet(netC);
    topMuxA->setNet(netMuxA);
    topMuxB->setNet(netMuxB);
    topMuxS->setNet(netMuxS);

    dffrnInst->getInstTerm(NLDB0::getDFFRNData())->setNet(netSeedOut);
    dffrnInst->getInstTerm(NLDB0::getDFFRNResetN())->setNet(netRN);
    dffrnInst->getInstTerm(NLDB0::getDFFRNClock())->setNet(netC);
    dffrnInst->getInstTerm(NLDB0::getDFFRNOutput())->setNet(netQ);

    assignSeedInst->getInstTerm(NLDB0::getAssignInput())->setNet(netSeedIn);
    assignSeedInst->getInstTerm(NLDB0::getAssignOutput())->setNet(netSeedOut);

    faInst->getInstTerm(NLDB0::getFAInputA())->setNet(netSeedOut);
    faInst->getInstTerm(NLDB0::getFAInputB())->setNet(netB);
    faInst->getInstTerm(NLDB0::getFAInputCI())->setNet(netCI);
    faInst->getInstTerm(NLDB0::getFAOutputS())->setNet(netFaS);

    assignSelInst->getInstTerm(NLDB0::getAssignInput())->setNet(netSelIn);
    assignSelInst->getInstTerm(NLDB0::getAssignOutput())->setNet(netSelOut);

    mux0Inst->getInstTerm(NLDB0::getMux2InputA()->getBit(0))->setNet(netQ);
    mux0Inst->getInstTerm(NLDB0::getMux2InputB()->getBit(0))->setNet(netFaS);
    mux0Inst->getInstTerm(NLDB0::getMux2Select())->setNet(netSelOut);
    mux0Inst->getInstTerm(NLDB0::getMux2Output()->getBit(0))->setNet(netSeedIn);

    mux1Inst->getInstTerm(NLDB0::getMux2InputA()->getBit(0))->setNet(netMuxA);
    mux1Inst->getInstTerm(NLDB0::getMux2InputB()->getBit(0))->setNet(netMuxB);
    mux1Inst->getInstTerm(NLDB0::getMux2Select())->setNet(netMuxS);
    mux1Inst->getInstTerm(NLDB0::getMux2Output()->getBit(0))->setNet(netSelIn);

    return top;
  };

  auto top = buildTop("top");
  auto topClone = top->clone(NLName("topClone"));
  ScopedCurrentPath scopedCurrentPath(tempDir_);
  const bool previousReportSkippedPOs = Config::getReportSkippedPOs();
  Config::setReportSkippedPOs(true);
  MiterStrategy MiterS(top, topClone, "InternalPOAssignCycle");
  MiterS.init();
  EXPECT_TRUE(MiterS.run());
  Config::setReportSkippedPOs(previousReportSkippedPOs);

  const auto reportPath = tempDir_ / "skipped_logical_loop_pos.txt";
  ASSERT_TRUE(std::filesystem::exists(reportPath));
  std::ifstream report(reportPath);
  ASSERT_TRUE(report.good());
  std::stringstream buffer;
  buffer << report.rdbuf();
  const std::string content = buffer.str();
  EXPECT_NE(content.find("logical loop"), std::string::npos);
  EXPECT_NE(content.find("loop_terms"), std::string::npos);
  EXPECT_NE(content.find("D0"), std::string::npos);
}

// 1. create a circuit of 2 inputs that drives and AND gate that drives top output
// 2. clone the the top and chain an inverter to the AND output
// 3. verify that the miter strategy detects the difference
TEST_F(MiterTests, TestMiterAndWithChainedInverter) {
  ScopedCurrentPath scopedCurrentPath(tempDir_);
  auto runKeplerCliWithArgs = [](const std::vector<std::string>& args) {
    std::string cmd;
    cmd += KEPLER_BIN;
    for (const auto& a : args) {
      cmd += " ";
      std::string quoted = "'";
      for (char c : a) {
        if (c == '\'') {
          quoted += "'\\''";
        } else {
          quoted.push_back(c);
        }
      }
      quoted += "'";
      cmd += quoted;
    }

    int rc = std::system(cmd.c_str());
    if (rc == -1) {
      return EXIT_FAILURE;
    }

#if defined(_WIN32)
    return rc;
#else
    if (WIFEXITED(rc)) {
      return WEXITSTATUS(rc);
    }
    return EXIT_FAILURE;
#endif
  };

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
  //NLLibraryTruthTables::construct(library);
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
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8, SNLTruthTable::fullDependencies(2)));
  // 8. create an inverter model
  SNLDesign* inverterModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("INV"));
  
   auto invIn =
      SNLScalarTerm::create(inverterModel, SNLTerm::Direction::Input, NLName("in"));
  auto invOut =
      SNLScalarTerm::create(inverterModel, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(inverterModel, SNLTruthTable(1, 1, SNLTruthTable::fullDependencies(1)));
  NLLibraryTruthTables::construct(library);
  // set truth table for inverter model
 
  
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
    std::filesystem::path outputPath("top.capnp");
    SNLCapnP::dump(db, outputPath);
  }
  // Dump visual
  {
    std::string dotFileName("beforeEdit.dot");
    std::string svgFileName("beforeEdit.svg");
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
    std::string dotFileName("afterEdit.dot");
    std::string svgFileName("afterEdit.svg");
    SnlVisualiser snl(top);
    snl.process();
    snl.getNetlistGraph().dumpDotFile(dotFileName.c_str());
    executeCommand(std::string(std::string("dot -Tsvg ") + dotFileName +
                               std::string(" -o ") + svgFileName)
                       .c_str());
  }

  // test the miter strategy
  {
    MiterStrategy MiterS(top, topClone, testTempPath("CaseC.log").string());
    MiterS.init();
    EXPECT_FALSE(MiterS.run());
  }
  {
    // dump top to naja_if(CapProto)
    std::filesystem::path outputPath("topEdited1.capnp");
    SNLCapnP::dump(db, outputPath);
  }
  const auto differentLog = testTempPath("different_miter.log");
  const auto differentCfg = testTempPath("different.yaml");
  {
    std::ofstream cfg(differentCfg);
    cfg << "format: naja_if\n";
    cfg << "input_paths:\n";
    cfg << "  - \"" << testTempPath("top.capnp").string() << "\"\n";
    cfg << "  - \"" << testTempPath("topEdited1.capnp").string() << "\"\n";
    cfg << "log_file: \"" << differentLog.string() << "\"\n";
  }
  int rc = runKeplerCliWithArgs({"--config", differentCfg.string()});
  EXPECT_EQ(rc, EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(differentLog));
  std::ifstream miterLogFile(differentLog);
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
    MiterStrategy MiterS(top, topClone, testTempPath("CaseD.log").string());
    MiterS.init();
    EXPECT_TRUE(MiterS.run());
  }
  {
    KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::SolverType::GLUCOSE);
    // print current solver type
    MiterStrategy MiterGlucose(top, topClone,
                               testTempPath("MultiDriver.log").string());
    MiterGlucose.init();
    // Expect throw in run
    EXPECT_TRUE(MiterGlucose.run());
    KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::SolverType::KISSAT);
  }
  {
    // dump top to naja_if(CapProto)
    std::filesystem::path outputPath("topEdited2.capnp");
    SNLCapnP::dump(db, outputPath);
  }

  const auto identicalLog = testTempPath("identical_miter.log");
  const auto identicalCfg = testTempPath("identical.yaml");
  {
    std::ofstream cfg(identicalCfg);
    cfg << "format: naja_if\n";
    cfg << "input_paths:\n";
    cfg << "  - \"" << testTempPath("top.capnp").string() << "\"\n";
    cfg << "  - \"" << testTempPath("topEdited2.capnp").string() << "\"\n";
    cfg << "log_file: \"" << identicalLog.string() << "\"\n";
  }
  rc = runKeplerCliWithArgs({"--config", identicalCfg.string()});
  EXPECT_EQ(rc, EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(identicalLog));
  std::ifstream miterLogFile2(identicalLog);
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

// ---------------------- Tests appended for coverage (subprocess approach, tolerant) ----------------------
// Append this block at the end of the file (after main).

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>

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
      if (c == '\'') {
        quoted += "'\\''";
      } else {
        quoted.push_back(c);
      }
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
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
  int rc = run_kepler_cli_with_args({});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, HelpFlagReturnsSuccess) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
  int rc = run_kepler_cli_with_args({"--help"});
  EXPECT_EQ(rc, EXIT_SUCCESS);
  rc = run_kepler_cli_with_args({"-h"});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, MissingConfigFileArgument) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
  int rc = run_kepler_cli_with_args({"--config"});
  EXPECT_NE(rc, EXIT_SUCCESS);
  rc = run_kepler_cli_with_args({"-c"});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ConfigFileNotFoundReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
  std::string tmpPath = "./nonexistent_config_12345.yaml";
  int rc = run_kepler_cli_with_args({"--config", tmpPath});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ConfigUnrecognizedFormatReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
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

TEST(KeplerCliSubprocessTests, ConfigUnknownKeyReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
  std::filesystem::path tmp = std::filesystem::temp_directory_path() / "kepler_test_unknown_key.yaml";
  {
    std::ofstream ofs(tmp);
    ofs << "format: verilog\n";
    ofs << "input_paths:\n  - a\n  - b\n";
    ofs << "cnf: true\n";
    ofs.close();
  }
  int rc = run_kepler_cli_with_args({"--config", tmp.string()});
  EXPECT_NE(rc, EXIT_SUCCESS);
  std::filesystem::remove(tmp);
}

TEST(KeplerCliSubprocessTests, ConfigSnlFormatLoadFailureReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
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
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
  int rc = run_kepler_cli_with_args({"-badformat", "a", "b"});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, CliNotEnoughPathsReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
  int rc = run_kepler_cli_with_args({"-verilog", "only_one_path.v"});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, CliNajaIfFormatButMissingFilesReturnsFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }
  int rc = run_kepler_cli_with_args({"-naja_if", "/no/such/file1.capnp", "/no/such/file2.capnp"});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ConfigParsingViaFilesCoversYamlToVectorBehavior) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }

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
  SNLDesignModeling::setTruthTable(logic0, SNLTruthTable(0, 0, SNLTruthTable::fullDependencies(0)));
  SNLDesignModeling::setTruthTable(logic1, SNLTruthTable(0, 1, SNLTruthTable::fullDependencies(0)));
  SNLDesign* inverterModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("INV"));
  // set truth table for inverter model
  auto invIn =
      SNLScalarTerm::create(inverterModel, SNLTerm::Direction::Input, NLName("in"));
  auto invOut =
      SNLScalarTerm::create(inverterModel, SNLTerm::Direction::Output, NLName("out"));
  SNLDesignModeling::setTruthTable(inverterModel, SNLTruthTable(1, 1, SNLTruthTable::fullDependencies(1)));
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
  //SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8, SNLTruthTable::fullDependencies(2)));
  // 9. connect all instances inputs
  SNLNet* net1 = SNLScalarNet::create(top, NLName("logic_0_net"));
  //net1->setType(SNLNet::Type::Assign0);
  SNLNet* net2 = SNLScalarNet::create(top, NLName("logic_1_net"));
  //net2->setType(SNLNet::Type::Assign1);
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
  andC0->getInstTerm(andIn1)->setNet(netC0a);
  andC0->getInstTerm(andIn2)->setNet(netC0b);
  topClone0->getBitTerm(topin->getID(), 0)->setNet(netC0a);
  SNLInstance* constC0 = SNLInstance::create(topClone0, logic0, NLName("logic0C0"));
  constC0->getInstTerm(logic0Out)->setNet(netC0b);

  SNLDesign* topClone1 = top->clone(NLName("topClone1"));
  SNLNet* netC1a = SNLScalarNet::create(topClone1, NLName("netC1a"));
  SNLNet* netC1b = SNLScalarNet::create(topClone1, NLName("netC1b"));
  auto andC1 = topClone1->getInstance(NLName("and"));
  andC1->getInstTerm(andIn1)->setNet(netC1a);
  andC1->getInstTerm(andIn2)->setNet(netC1b);
  SNLInstance* constC1 = SNLInstance::create(topClone1, logic0, NLName("logic0C1"));
  
  
  auto inverterC1 = SNLInstance::create(topClone1, inverterModel, NLName("inverterC1"));
  constC1->getInstTerm(logic0Out)->setNet(netC1a);
  auto netC1invOut = SNLScalarNet::create(topClone1, NLName("netC1invOut"));
  inverterC1->getInstTerm(invIn)->setNet(netC1b);
  inverterC1->getInstTerm(invOut)->setNet(netC1invOut);
  andC1->getInstTerm(andIn1)->setNet(netC1invOut);
  topClone1->getBitTerm(topin2->getID(), 0)->setNet(netC1b);
  
  // 11. create DNL
  get(); 
  naja::DNL::destroy();
  MiterStrategy MiterS(topClone0, topClone1, "CaseD");
  MiterS.init();
    EXPECT_FALSE(MiterS.run());
}

// Test error for multiple drivers
TEST_F(MiterTests, multiDriver) {
  // 1. Create SNL
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryS =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("Stadarts"));
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  // 2. Create a top model with one output
  SNLDesign* top =
      SNLDesign::create(libraryS, SNLDesign::Type::Standard, NLName("top"));
  univ->setTopDesign(top);
  auto topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
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
  SNLDesignModeling::setTruthTable(logic0, SNLTruthTable(0, 0, SNLTruthTable::fullDependencies(0)));
  SNLDesignModeling::setTruthTable(logic1, SNLTruthTable(0, 1, SNLTruthTable::fullDependencies(0)));
  //NLLibraryTruthTables::construct(library);
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
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8, SNLTruthTable::fullDependencies(2)));
  NLLibraryTruthTables::construct(library);
  // 9. connect all instances inputs
  SNLNet* net1 = SNLScalarNet::create(top, NLName("net1"));
  SNLNet* net2 = SNLScalarNet::create(top, NLName("net2"));
 
  // connect logic0 to and
  inst1->getInstTerm(logic0Out)->setNet(net1);
  topIn->setNet(net1);

  inst4->getInstTerm(andIn1)->setNet(net1);
  inst4->getInstTerm(andIn2)->setNet(net1);
  // connect logic1 to and
  inst2->getInstTerm(logic1Out)->setNet(net1);
  inst3->getInstTerm(andIn2)->setNet(net1);
  inst3->getInstTerm(andIn1)->setNet(net1);
  // connect the and instance output to the top output
  inst3->getInstTerm(andOut)->setNet(net1);
  topOut->setNet(net2);
  inst4->getInstTerm(andOut)->setNet(net2);
  auto topClone = top->clone(NLName("topClone"));
  // 11. create DNL
  MiterStrategy MiterS(top, topClone, "MultiDriver");
  MiterS.init();
  EXPECT_TRUE(MiterS.run());
  naja::DNL::destroy();
}

// Test error for multiple drivers
TEST_F(MiterTests, tt65In) {
  // 1. Create SNL
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryS =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("Stadarts"));
  NLLibrary* library =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("nangate45"));
  // 2. Create a top model with one output
  SNLDesign* top =
      SNLDesign::create(libraryS, SNLDesign::Type::Standard, NLName("top"));
  univ->setTopDesign(top);
  auto topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
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
  SNLDesignModeling::setTruthTable(logic0, SNLTruthTable(0, 0, SNLTruthTable::fullDependencies(0)));
  SNLDesignModeling::setTruthTable(logic1, SNLTruthTable(0, 1, SNLTruthTable::fullDependencies(0)));
  //NLLibraryTruthTables::construct(library);
  // Create a model with 65 inputs and 1 output and set the truth table so 
  // output is 1 only when all inputs are 0
  SNLDesign* tt65InModel =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("TT65IN"));
  std::vector<SNLScalarTerm*> tt65InTerms;
  std::vector<SNLScalarTerm*> topOutTerms;
  std::vector<SNLBitNet*> topOutNets;
  for (int i = 0; i < 65; ++i) {
    auto outTerm = SNLScalarTerm::create(tt65InModel, SNLTerm::Direction::Output,
                                       NLName("out" + std::to_string(i)));
    auto topOut = SNLScalarTerm::create(top, SNLTerm::Direction::Output,
                                       NLName("in" + std::to_string(i)));
    topOutTerms.push_back(topOut);

    tt65InTerms.push_back(outTerm);
    auto topOutNet = SNLScalarNet::create(top, NLName("in_net" + std::to_string(i)));
    topOutNets.push_back(topOutNet);
  }
  auto tt65In = SNLScalarTerm::create(tt65InModel, SNLTerm::Direction::Input,
                                        NLName("in"));
  auto tt65In2 = SNLScalarTerm::create(tt65InModel, SNLTerm::Direction::Input,
                                        NLName("in2"));
  // set truth tables for all 65 outputs with and function for the 2 inputs
  std::vector<SNLTruthTable> tt65InTables;
  for (int i = 0; i < 65; ++i) {
    tt65InTables.push_back(
        SNLTruthTable(2, 8, getInputFlatDependencies(tt65InModel)));
  }
  SNLDesignModeling::setTruthTables(tt65InModel,tt65InTables);
  //NLLibraryTruthTables::construct(library);
  // create the instance of the model in top
  auto tt65InInst = SNLInstance::create(top, tt65InModel, NLName("tt65in"));
  // 5. create a logic_0 instace in top
  SNLInstance* inst1 = SNLInstance::create(top, logic0, NLName("logic0"));
  // 6. create a logic_1 instace in top
  SNLInstance* inst2 = SNLInstance::create(top, logic1, NLName("logic1"));
  // create 64 nets that will be connected to the tt65In first 64 inputs
  std::vector<SNLNet*> tt65InNets;
  // connect the 65 outputs to the 65 top outputs
  for (int i = 0; i < 65; ++i) {
    tt65InInst->getInstTerm(tt65InTerms[i])->setNet(topOutNets[i]);
    topOutTerms[i]->setNet(topOutNets[i]);
  }
  // connect the last on to the top in
  
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
  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8, SNLTruthTable::fullDependencies(2)));
  //NLLibraryTruthTables::construct(library);
  // 9. connect all instances inputs
  SNLNet* net1 = SNLScalarNet::create(top, NLName("net1"));
  SNLNet* net2 = SNLScalarNet::create(top, NLName("net2"));
  auto netIn = SNLScalarNet::create(top, NLName("net_in"));
  // connect the 65th input to top in
  // connect logic0 to and
  inst1->getInstTerm(logic0Out)->setNet(net1);
  topIn->setNet(netIn);

  tt65InInst->getInstTerm(tt65In)->setNet(netIn);
  tt65InInst->getInstTerm(tt65In2)->setNet(net1);
  // connect out of tt65In to topOut2
  auto net_tt65Out = SNLScalarNet::create(top, NLName("tt65out_net"));
  topOut2->setNet(net_tt65Out);

  inst4->getInstTerm(andIn1)->setNet(net1);
  inst4->getInstTerm(andIn2)->setNet(net1);
  // connect logic1 to and
  inst2->getInstTerm(logic1Out)->setNet(net1);
  inst3->getInstTerm(andIn2)->setNet(net1);
  inst3->getInstTerm(andIn1)->setNet(net1);
  // connect the and instance output to the top output
  inst3->getInstTerm(andOut)->setNet(net1);
  topOut->setNet(net2);
  inst4->getInstTerm(andOut)->setNet(net2);
  auto topClone = top->clone(NLName("topClone"));
  // 11. create DNL
  MiterStrategy MiterS(top, topClone, "MultiDriver");
  MiterS.init();
  EXPECT_TRUE(MiterS.run());
}

TEST_F(MiterTests, ConnectedInouts) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryS =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("Stadarts"));
  // create a top model with 2 inout ports
  SNLDesign* top =
      SNLDesign::create(libraryS, SNLDesign::Type::Standard, NLName("top"));
  univ->setTopDesign(top);
  auto topInout1 =
      SNLScalarTerm::create(top, SNLTerm::Direction::InOut, NLName("inout1"));
  auto topInout2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::InOut, NLName("inout2"));
  // connect the 2 inouts together
  SNLNet* net = SNLScalarNet::create(top, NLName("net"));
  topInout1->setNet(net);
  topInout2->setNet(net);
  SNLDesign* topClone = top->clone(NLName("topClone"));
  MiterStrategy MiterS(top, topClone, "ConnectedInouts");
  MiterS.init();
  EXPECT_TRUE(MiterS.run());
}

TEST_F(MiterTests, UnconnectedTerms) {
  NLUniverse* univ = NLUniverse::create();
  NLDB* db = NLDB::create(univ);
  NLLibrary* libraryS =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("Stadarts"));
  // create a top model with 2 inout ports
  SNLDesign* top =
      SNLDesign::create(libraryS, SNLDesign::Type::Standard, NLName("top"));
  univ->setTopDesign(top);
  auto topInout1 =
      SNLScalarTerm::create(top, SNLTerm::Direction::InOut, NLName("inout1"));
  auto topInout2 =
      SNLScalarTerm::create(top, SNLTerm::Direction::InOut, NLName("inout2"));
  SNLDesign* topClone = top->clone(NLName("topClone"));
  MiterStrategy MiterS(top, topClone, "ConnectedInouts");
  MiterS.init();
  EXPECT_TRUE(MiterS.run());
}

TEST(KeplerCliSubprocessTests, ExampleTestRun) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }

  std::string config = get_test_data_prefix() + "test/strategies/miter/test_config_verilog.yaml";
  if (std::getenv("TEST_DATA_PREFIX")) {
    config = get_test_data_prefix() + "test/strategies/miter/test_config_verilog_bazel.yaml";
  }
  int rc = run_kepler_cli_with_args({"--config", config});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ExampleTestRunCommandLine) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }

  std::string pfx = get_test_data_prefix();
  int rc = run_kepler_cli_with_args({"-verilog",
                                         pfx + "example/tinyrocket.v",
                                         pfx + "example/tinyrocket_edited.v",
                                         pfx + "example/NangateOpenCellLibrary_typical.lib",
                                         pfx + "example/fakeram45_64x15.lib",
                                         pfx + "example/fakeram45_64x32.lib",
                                         pfx + "example/fakeram45_1024x32.lib"});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ExampleTestRunNajaIFWithScopeExtraction) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }

  std::string config = get_test_data_prefix() + "test/strategies/miter/test_config_naja_if_with_se.yaml";
  if (std::getenv("TEST_DATA_PREFIX")) {
    config = get_test_data_prefix() + "test/strategies/miter/test_config_naja_if_with_se_bazel.yaml";
  }
  int rc = run_kepler_cli_with_args({"--config", config});
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

// test failure with test_config_failure.yaml
TEST(KeplerCliSubprocessTests, ExampleTestRunFailure) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }

  std::string config = get_test_data_prefix() + "test/strategies/miter/test_config_failure.yaml";
  if (std::getenv("TEST_DATA_PREFIX")) {
    config = get_test_data_prefix() + "test/strategies/miter/test_config_failure_bazel.yaml";
  }
  int rc = run_kepler_cli_with_args({"--config", config});
  EXPECT_NE(rc, EXIT_SUCCESS);
}

TEST(KeplerCliSubprocessTests, ExampleRunWritesConfiguredLogFile) {
  std::filesystem::path p(KEPLER_BIN);
  if (!std::filesystem::exists(p)) {
    GTEST_SKIP() << "kepler-formal binary missing";
  }

  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_subprocess_log";
  std::filesystem::create_directories(tmpDir);
  const auto logPath = tmpDir / "configured_miter.log";
  const auto configPath = tmpDir / "config.yaml";
  const auto root = repoRoot();

  {
    std::ofstream cfg(configPath);
    cfg << "format: verilog\n";
    cfg << "input_paths:\n";
    cfg << "  - " << (root / "example/tinyrocket.v").string() << "\n";
    cfg << "  - " << (root / "example/tinyrocket_edited.v").string() << "\n";
    cfg << "liberty_files:\n";
    cfg << "  - " << (root / "example/NangateOpenCellLibrary_typical.lib").string() << "\n";
    cfg << "  - " << (root / "example/fakeram45_1024x32.lib").string() << "\n";
    cfg << "  - " << (root / "example/fakeram45_64x32.lib").string() << "\n";
    cfg << "  - " << (root / "example/fakeram45_64x15.lib").string() << "\n";
    cfg << "log_file: " << logPath.string() << "\n";
  }

  if (std::filesystem::exists(logPath)) {
    std::filesystem::remove(logPath);
  }

  int rc = run_kepler_cli_with_args({"--config", configPath.string()});
  EXPECT_EQ(rc, EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(logPath));

  std::ifstream logFile(logPath);
  std::string contents((std::istreambuf_iterator<char>(logFile)),
                       std::istreambuf_iterator<char>());
  EXPECT_NE(contents.find("DIFFERENT"), std::string::npos);

  std::filesystem::remove_all(tmpDir);
}

// Required main function for Google Test
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// End of appended tests
