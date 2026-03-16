// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <cstdlib>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "DNL.h"
#include "KeplerFormalUtils.h"
#include "SNLCapnP.h"
#include "SNLLibertyConstructor.h"
#include "SNLUtils.h"
#include "SNLVRLConstructor.h"

extern int KeplerFormalMain(int argc, char** argv);

namespace {

std::filesystem::path repoRoot() {
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path().parent_path();
}

std::filesystem::path writeTempConfig(const std::string& contents) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_cli_cfg";
  std::filesystem::create_directories(tmpDir);
  const auto cfgPath =
      tmpDir / ("config_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".yaml");
  std::ofstream cfg(cfgPath);
  cfg << contents;
  cfg.close();
  return cfgPath;
}

int runWithConfigFile(const std::filesystem::path& cfgPath) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "--config";
  std::string argv2 = cfgPath.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data()};
  int argc = 3;
  return KeplerFormalMain(argc, argv);
}

struct MultiFileVerilogFixture {
  std::filesystem::path tmpDir;
  std::filesystem::path cfgPath;
};

MultiFileVerilogFixture createVerilogPreprocessingFixture(bool enablePreprocessing) {
  MultiFileVerilogFixture fixture;
  fixture.tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_cli_preproc_v";
  std::filesystem::create_directories(fixture.tmpDir);

  const auto design0Aux = fixture.tmpDir / "design0_aux.v";
  const auto design0Top = fixture.tmpDir / "design0_top.v";
  const auto design1Aux = fixture.tmpDir / "design1_aux.v";
  const auto design1Top = fixture.tmpDir / "design1_top.v";
  fixture.cfgPath = fixture.tmpDir / "config.yaml";

  {
    std::ofstream f(design0Aux);
    f << "module AUX_MACRO (\n";
    f << "    input a,\n";
    f << "    output y\n";
    f << ");\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1Aux);
    f << "module AUX_MACRO (\n";
    f << "    input a,\n";
    f << "    output y\n";
    f << ");\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design0Top);
    f << "`define PASS_A a\n";
    f << "module top(input a, output y);\n";
    f << "  wire y_aux;\n";
    f << "  AUX_MACRO u_aux(.a(`PASS_A), .y(y_aux));\n";
    f << "  assign y = y_aux;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1Top);
    f << "`define PASS_A a\n";
    f << "module top(input a, output y);\n";
    f << "  wire y_aux;\n";
    f << "  AUX_MACRO u_aux(.a(`PASS_A), .y(y_aux));\n";
    f << "  assign y = y_aux;\n";
    f << "endmodule\n";
  }

  std::ofstream cfg(fixture.cfgPath);
  cfg << "format: verilog\n";
  cfg << "input_paths:\n";
  cfg << "  -\n";
  cfg << "    - " << design0Aux.string() << "\n";
  cfg << "    - " << design0Top.string() << "\n";
  cfg << "  -\n";
  cfg << "    - " << design1Aux.string() << "\n";
  cfg << "    - " << design1Top.string() << "\n";
  cfg << "verilog_preprocessing: " << (enablePreprocessing ? "true" : "false")
      << "\n";
  cfg.close();

  return fixture;
}

MultiFileVerilogFixture createDefaultNettypeDirectiveFixture(bool enablePreprocessing) {
  MultiFileVerilogFixture fixture;
  fixture.tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_cli_default_nettype_v";
  std::filesystem::create_directories(fixture.tmpDir);

  const auto design0 = fixture.tmpDir / "design0.v";
  const auto design1 = fixture.tmpDir / "design1.v";
  fixture.cfgPath = fixture.tmpDir / "config.yaml";

  {
    std::ofstream f(design0);
    f << "`timescale 1ns / 1ps\n";
    f << "`default_nettype none\n";
    f << "module top(input a, output y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1);
    f << "`timescale 1ns / 1ps\n";
    f << "`default_nettype none\n";
    f << "module top(input a, output y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }

  std::ofstream cfg(fixture.cfgPath);
  cfg << "format: verilog\n";
  cfg << "input_paths:\n";
  cfg << "  - " << design0.string() << "\n";
  cfg << "  - " << design1.string() << "\n";
  cfg << "verilog_preprocessing: " << (enablePreprocessing ? "true" : "false")
      << "\n";
  cfg.close();

  return fixture;
}
struct SimpleCliFixture {
  std::filesystem::path tmpDir;
  std::filesystem::path design0Path;
  std::filesystem::path design1Path;
};

struct ScopedNajaIfFixture {
  std::filesystem::path tmpDir;
  std::filesystem::path design0IfPath;
  std::filesystem::path design1IfPath;
  std::filesystem::path libertyPath;
};

struct EnvVarGuard {
  explicit EnvVarGuard(const char* name): name_(name) {
    const char* current = std::getenv(name_);
    if (current) {
      hadValue_ = true;
      oldValue_ = current;
    }
  }

  void set(const std::string& value) const {
    setenv(name_, value.c_str(), 1);
  }

  ~EnvVarGuard() {
    if (hadValue_) {
      setenv(name_, oldValue_.c_str(), 1);
    } else {
      unsetenv(name_);
    }
  }

  const char* name_;
  bool hadValue_ = false;
  std::string oldValue_;
};

MultiFileVerilogFixture createMultiFileVerilogFixture() {
  MultiFileVerilogFixture fixture;
  fixture.tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_cli_multi_v";
  std::filesystem::create_directories(fixture.tmpDir);

  const auto design0Leaf = fixture.tmpDir / "design0_leaf.v";
  const auto design0Top = fixture.tmpDir / "design0_top.v";
  const auto design1Leaf = fixture.tmpDir / "design1_leaf.v";
  const auto design1Top = fixture.tmpDir / "design1_top.v";
  fixture.cfgPath = fixture.tmpDir / "config.yaml";
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto lib0 = exampleDir / "NangateOpenCellLibrary_typical.lib";
  const auto lib1 = exampleDir / "fakeram45_1024x32.lib";
  const auto lib2 = exampleDir / "fakeram45_64x32.lib";
  const auto lib3 = exampleDir / "fakeram45_64x15.lib";

  {
    std::ofstream leaf(design0Leaf);
    leaf << "module leaf(input a, input b, output y);\n";
    leaf << "  wire n;\n";
    leaf << "  NAND2_X1 u1(.A1(a), .A2(b), .ZN(n));\n";
    leaf << "  INV_X1 u2(.A(n), .ZN(y));\n";
    leaf << "endmodule\n";
  }
  {
    std::ofstream top(design0Top);
    top << "module top(input a, input b, output y);\n";
    top << "  wire w;\n";
    top << "  leaf u1(.a(a), .b(b), .y(w));\n";
    top << "  INV_X1 u2(.A(w), .ZN(y));\n";
    top << "endmodule\n";
  }
  {
    std::ofstream leaf(design1Leaf);
    leaf << "module leaf(input a, input b, output y);\n";
    leaf << "  wire n;\n";
    leaf << "  NAND2_X1 u1(.A1(a), .A2(b), .ZN(n));\n";
    leaf << "  INV_X1 u2(.A(n), .ZN(y));\n";
    leaf << "endmodule\n";
  }
  {
    std::ofstream top(design1Top);
    top << "module top(input a, input b, output y);\n";
    top << "  wire w;\n";
    top << "  leaf u1(.a(a), .b(b), .y(w));\n";
    top << "  INV_X1 u2(.A(w), .ZN(y));\n";
    top << "endmodule\n";
  }

  std::ofstream cfg(fixture.cfgPath);
  cfg << "format: verilog\n";
  cfg << "input_paths:\n";
  cfg << "  -\n";
  cfg << "    - " << design0Leaf.string() << "\n";
  cfg << "    - " << design0Top.string() << "\n";
  cfg << "  -\n";
  cfg << "    - " << design1Leaf.string() << "\n";
  cfg << "    - " << design1Top.string() << "\n";
  cfg << "liberty_files:\n";
  cfg << "  - " << lib0.string() << "\n";
  cfg << "  - " << lib1.string() << "\n";
  cfg << "  - " << lib2.string() << "\n";
  cfg << "  - " << lib3.string() << "\n";
  cfg << "log_level: info\n";
  cfg.close();

  return fixture;
}

SimpleCliFixture createEquivalentDesignFixture(const std::string& extension,
                                               const std::string& moduleBody) {
  SimpleCliFixture fixture;
  fixture.tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_cli_simple";
  std::filesystem::create_directories(fixture.tmpDir);

  fixture.design0Path = fixture.tmpDir / ("design0." + extension);
  fixture.design1Path = fixture.tmpDir / ("design1." + extension);

  {
    std::ofstream design0(fixture.design0Path);
    design0 << moduleBody;
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << moduleBody;
  }

  return fixture;
}

void cleanupNajaTestState() {
  naja::DNL::destroy();
  if (NLUniverse::get()) {
    NLUniverse::get()->destroy();
  }
}

ScopedNajaIfFixture createEquivalentScopedNajaIfFixture() {
  ScopedNajaIfFixture fixture;
  fixture.tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_cli_scope_if";
  std::filesystem::create_directories(fixture.tmpDir);
  fixture.design0IfPath = fixture.tmpDir / "design0.capnp";
  fixture.design1IfPath = fixture.tmpDir / "design1.capnp";
  fixture.libertyPath = repoRoot() / "example" / "NangateOpenCellLibrary_typical.lib";

  const auto design0Child = fixture.tmpDir / "design0_child.v";
  const auto design0Top = fixture.tmpDir / "design0_top.v";
  const auto design1Child = fixture.tmpDir / "design1_child.v";
  const auto design1Top = fixture.tmpDir / "design1_top.v";

  {
    std::ofstream child(design0Child);
    child << "module child(input a, output y);\n";
    child << "  assign y = a;\n";
    child << "endmodule\n";
  }
  {
    std::ofstream top(design0Top);
    top << "module top(input a, output y);\n";
    top << "  child u_child(.a(a), .y(y));\n";
    top << "endmodule\n";
  }
  {
    std::ofstream child(design1Child);
    child << "module child(input a, output y);\n";
    child << "  wire n;\n";
    child << "  INV_X1 u0(.A(a), .ZN(n));\n";
    child << "  INV_X1 u1(.A(n), .ZN(y));\n";
    child << "endmodule\n";
  }
  {
    std::ofstream top(design1Top);
    top << "module top(input a, output y);\n";
    top << "  child u_child(.a(a), .y(y));\n";
    top << "endmodule\n";
  }

  const auto dumpDesign = [&](const std::vector<std::filesystem::path>& designPaths,
                              const std::filesystem::path& dumpPath) {
    cleanupNajaTestState();
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    auto* primitivesLibrary =
        NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
    SNLLibertyConstructor libertyConstructor(primitivesLibrary);
    libertyConstructor.construct(fixture.libertyPath);
    auto* designLibrary = NLLibrary::create(db, NLName("DESIGN"));
    SNLVRLConstructor constructor(designLibrary);
    constructor.construct(designPaths);
    auto* top = SNLUtils::findTop(designLibrary);
    ASSERT_NE(top, nullptr);
    db->setTopDesign(top);
    SNLCapnP::dump(db, dumpPath);
    cleanupNajaTestState();
  };

  dumpDesign({design0Child, design0Top}, fixture.design0IfPath);
  dumpDesign({design1Child, design1Top}, fixture.design1IfPath);

  return fixture;
}

}  // namespace

TEST(KeplerFormalCliTests, SanitizeFileToken) {
  EXPECT_EQ(sanitizeFileToken("scope"), "scope");
  EXPECT_EQ(sanitizeFileToken("my scope"), "my_scope");
  EXPECT_EQ(sanitizeFileToken("a/b\\c"), "a_b_c");
  EXPECT_EQ(sanitizeFileToken(""), "scope");
}


TEST(KeplerFormalCliTests, DumpCnfFromConfig) {
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto design0 = exampleDir / "tinyrocket.v";
  const auto design1 = exampleDir / "tinyrocket_edited.v";
  const auto lib0 = exampleDir / "NangateOpenCellLibrary_typical.lib";
  const auto lib1 = exampleDir / "fakeram45_1024x32.lib";
  const auto lib2 = exampleDir / "fakeram45_64x32.lib";
  const auto lib3 = exampleDir / "fakeram45_64x15.lib";

  ASSERT_TRUE(std::filesystem::exists(design0));
  ASSERT_TRUE(std::filesystem::exists(design1));
  ASSERT_TRUE(std::filesystem::exists(lib0));
  ASSERT_TRUE(std::filesystem::exists(lib1));
  ASSERT_TRUE(std::filesystem::exists(lib2));
  ASSERT_TRUE(std::filesystem::exists(lib3));

  const auto tmpDir = std::filesystem::temp_directory_path() / "kepler_formal_cli_test";
  std::filesystem::create_directories(tmpDir);
  const auto cnfPath = tmpDir / "miter_test.cnf";
  const auto cfgPath = tmpDir / "config.yaml";

  if (std::filesystem::exists(cnfPath)) {
    std::filesystem::remove(cnfPath);
  }

  std::ofstream cfg(cfgPath);
  cfg << "format: verilog\n";
  cfg << "input_paths:\n";
  cfg << "  - " << design0.string() << "\n";
  cfg << "  - " << design1.string() << "\n";
  cfg << "liberty_files:\n";
  cfg << "  - " << lib0.string() << "\n";
  cfg << "  - " << lib1.string() << "\n";
  cfg << "  - " << lib2.string() << "\n";
  cfg << "  - " << lib3.string() << "\n";
  cfg << "log_level: info\n";
  cfg << "cnf_export: true\n";
  cfg << "cnf_export_path: " << cnfPath.string() << "\n";
  cfg.close();

  std::string argv0 = "kepler-formal";
  std::string argv1 = "--config";
  std::string argv2 = cfgPath.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data()};
  int argc = 3;

  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  EXPECT_TRUE(std::filesystem::exists(cnfPath));

  std::filesystem::remove(cnfPath);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(tmpDir);
}

TEST(KeplerFormalCliTests, MultiFileVerilogConfig) {
  const auto fixture = createMultiFileVerilogFixture();
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, YamlMultiFileVerilogConfig) {
  const auto fixture = createMultiFileVerilogFixture();
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, VerilogPreprocessingEnabledParsesDirectiveInput) {
  const auto fixture = createVerilogPreprocessingFixture(true);
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, VerilogPreprocessingDisabledFailsOnDirectiveInput) {
  const auto fixture = createVerilogPreprocessingFixture(false);
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, VerilogPreprocessingEnabledDefaultNettypeIsRejected) {
  const auto fixture = createDefaultNettypeDirectiveFixture(true);
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, VerilogPreprocessingDisabledDefaultNettypeIsRejected) {
  const auto fixture = createDefaultNettypeDirectiveFixture(false);
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliVerilogPreprocessingFlagEnablesDirectiveInput) {
  const auto fixture = createVerilogPreprocessingFixture(true);
  const auto design0Aux = fixture.tmpDir / "design0_aux.v";
  const auto design0Top = fixture.tmpDir / "design0_top.v";
  const auto design1Aux = fixture.tmpDir / "design1_aux.v";
  const auto design1Top = fixture.tmpDir / "design1_top.v";

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--design1";
  std::string argv3 = design0Aux.string();
  std::string argv4 = design0Top.string();
  std::string argv5 = "--design2";
  std::string argv6 = design1Aux.string();
  std::string argv7 = design1Top.string();
  std::string argv8 = "--verilog_preprocessing";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data(),
                  argv8.data()};
  int argc = 9;

  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliWithoutVerilogPreprocessingFlagFailsOnDirectiveInput) {
  const auto fixture = createVerilogPreprocessingFixture(true);
  const auto design0Aux = fixture.tmpDir / "design0_aux.v";
  const auto design0Top = fixture.tmpDir / "design0_top.v";
  const auto design1Aux = fixture.tmpDir / "design1_aux.v";
  const auto design1Top = fixture.tmpDir / "design1_top.v";

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--design1";
  std::string argv3 = design0Aux.string();
  std::string argv4 = design0Top.string();
  std::string argv5 = "--design2";
  std::string argv6 = design1Aux.string();
  std::string argv7 = design1Top.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;

  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigMissingInputPathsFails) {
  const auto cfgPath = writeTempConfig("format: verilog\nlog_level: info\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigNestedSecondNotSequenceFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - [a.v]\n"
      "  - b.v\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigNestedEmptyDesignFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - []\n"
      "  - []\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigEmptyInputPathsFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: []\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigInputPathsNotSequenceFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: foo\n"
      "log_level: info\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigUnknownKeyFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "bogus_key: 1\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigNonScalarKeyFails) {
  const auto cfgPath = writeTempConfig(
      "? [a, b]\n"
      ": 1\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigUnknownFormatFails) {
  const auto cfgPath = writeTempConfig(
      "format: unknown_format\n"
      "input_paths: [a.v, b.v]\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigSystemVerilogAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic a, output logic y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSvAliasAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic a, output logic y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-sv";
  std::string argv2 = fixture.design0Path.string();
  std::string argv3 = fixture.design1Path.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data()};
  int argc = 4;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, FirstVerilogDesignWithoutTopFails) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_no_top_first";
  std::filesystem::create_directories(tmpDir);
  const auto design0 = tmpDir / "design0.v";
  const auto design1 = tmpDir / "design1.v";
  {
    std::ofstream f(design0);
    f << "module a(input x, output y);\n";
    f << "  assign y = x;\n";
    f << "endmodule\n";
    f << "module b(input x, output y);\n";
    f << "  assign y = x;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1);
    f << "module top(input a, output y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = design0.string();
  std::string argv3 = design1.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data()};
  int argc = 4;

  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove_all(tmpDir);
}

TEST(KeplerFormalCliTests, SecondVerilogDesignWithoutTopFails) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_no_top_second";
  std::filesystem::create_directories(tmpDir);
  const auto design0 = tmpDir / "design0.v";
  const auto design1 = tmpDir / "design1.v";
  {
    std::ofstream f(design0);
    f << "module top(input a, output y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1);
    f << "module a(input x, output y);\n";
    f << "  assign y = x;\n";
    f << "endmodule\n";
    f << "module b(input x, output y);\n";
    f << "  assign y = x;\n";
    f << "endmodule\n";
  }

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = design0.string();
  std::string argv3 = design1.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data()};
  int argc = 4;

  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove_all(tmpDir);
}

TEST(KeplerFormalCliTests, ConfigUnknownSolverFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "solver: nope\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigLogFileAccepted) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_cli_log";
  std::filesystem::create_directories(tmpDir);
  const auto logPath = tmpDir / "kepler_formal_test.log";
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "log_file: " + logPath.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(tmpDir);
}

TEST(KeplerFormalCliTests, ConfigInputPathsWrongSizeFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - only_one.v\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigInputPathsNestedWrongCountFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - [a.v]\n"
      "  - [b.v]\n"
      "  - [c.v]\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigInputPathsNestedNonScalarFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  -\n"
      "    - [nested]\n"
      "  -\n"
      "    - b.v\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, SnlMultiFileRejected) {
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  -\n"
      "    - a.if\n"
      "    - b.if\n"
      "  -\n"
      "    - c.if\n"
      "    - d.if\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, MissingFirstNajaIfFails) {
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto design1 = exampleDir / "tinyrocket_naja.if";
  ASSERT_TRUE(std::filesystem::exists(design1));

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - /definitely/missing/first.if\n"
      "  - " + design1.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, MissingSecondNajaIfFails) {
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto design0 = exampleDir / "tinyrocket_naja.if";
  ASSERT_TRUE(std::filesystem::exists(design0));

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + design0.string() + "\n"
      "  - /definitely/missing/second.if\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, CliUnknownOptionFails) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--design1";
  std::string argv3 = "a.v";
  std::string argv4 = "--design2";
  std::string argv5 = "b.v";
  std::string argv6 = "--bogus";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data()};
  int argc = 7;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliHelpPrintsUsage) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "--help";
  char* argv[] = {argv0.data(), argv1.data()};
  int argc = 2;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST(KeplerFormalCliTests, CliLibertyFlagCollectsPaths) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--design1";
  std::string argv3 = "a.v";
  std::string argv4 = "--design2";
  std::string argv5 = "b.v";
  std::string argv6 = "--liberty";
  std::string argv7 = "lib1.lib";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, PythonLibraryFilesAreLoadedByExtension) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(output y);\n"
      "  wire z;\n"
      "  LOGIC1 c0(.Z(z));\n"
      "  BUF c1(.A(z), .Z(y));\n"
      "endmodule\n");
  const auto pyPrimitives =
      repoRoot() / "thirdparty/naja/test/nl/python/pyloader/scripts/primitives1.py";
  const auto pyModuleDir =
      repoRoot() / "buildD/thirdparty/naja/src/nl/python/naja_wrapping";
  ASSERT_TRUE(std::filesystem::exists(pyPrimitives));
  ASSERT_TRUE(std::filesystem::exists(pyModuleDir / "naja.so"));
  EnvVarGuard pythonPathGuard("PYTHONPATH");
  pythonPathGuard.set(pyModuleDir.string());

  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "liberty_files:\n"
      "  - " + pyPrimitives.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, GzippedLibertyFilesAreLoadedByExtension) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto gzLiberty =
      repoRoot() / "thirdparty/naja/test/nl/formats/liberty/benchmarks/tests/small.lib.gz";
  ASSERT_TRUE(std::filesystem::exists(gzLiberty));

  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "liberty_files:\n"
      "  - " + gzLiberty.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, UnsupportedLibraryExtensionFails) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto badLibPath = fixture.tmpDir / "bad_library.txt";
  {
    std::ofstream badLib(badLibPath);
    badLib << "not a supported library file\n";
  }

  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "liberty_files:\n"
      "  - " + badLibPath.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliExplicitMissingDesignsFails) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--design1";
  std::string argv3 = "--design2";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data()};
  int argc = 4;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, ConfigDebugLogLevelAccepted) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "log_level: debug\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigUnknownLogLevelDefaultsToInfo) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "log_level: loud\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, VerilogNoLibertyCreatesDbAndFailsOnSecondParse) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_no_lib";
  std::filesystem::create_directories(tmpDir);
  const auto design0 = tmpDir / "design0.v";
  const auto design1 = tmpDir / "missing_design1.v";
  {
    std::ofstream f(design0);
    f << "module top(input a, output y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = design0.string();
  std::string argv3 = design1.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data()};
  int argc = 4;

  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove_all(tmpDir);
}

TEST(KeplerFormalCliTests, SnlScopesNoDifference) {
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto design0 = exampleDir / "tinyrocket_naja.if";
  const auto lib0 = exampleDir / "NangateOpenCellLibrary_typical.lib";
  const auto lib1 = exampleDir / "fakeram45_1024x32.lib";
  const auto lib2 = exampleDir / "fakeram45_64x32.lib";

  ASSERT_TRUE(std::filesystem::exists(design0));
  ASSERT_TRUE(std::filesystem::exists(lib0));
  ASSERT_TRUE(std::filesystem::exists(lib1));
  ASSERT_TRUE(std::filesystem::exists(lib2));

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + design0.string() + "\n"
      "  - " + design0.string() + "\n"
      "liberty_files:\n"
      "  - " + lib0.string() + "\n"
      "  - " + lib1.string() + "\n"
      "  - " + lib2.string() + "\n"
      "use_scopes: true\n");

  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, SnlScopesEquivalentEditedScopeNoDifference) {
  const auto fixture = createEquivalentScopedNajaIfFixture();
  ASSERT_TRUE(std::filesystem::exists(fixture.design0IfPath));
  ASSERT_TRUE(std::filesystem::exists(fixture.design1IfPath));
  ASSERT_TRUE(std::filesystem::exists(fixture.libertyPath));

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "liberty_files:\n"
      "  - " + fixture.libertyPath.string() + "\n"
      "use_scopes: true\n");

  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}
