// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "DNL.h"
#include "Config.h"
#include "KeplerFormalUtils.h"
#include "NLDB0.h"
#include "NLUniverse.h"
#include "SNLCapnP.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLLibertyConstructor.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "SNLTruthTable.h"
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

std::filesystem::path makeUniqueTempDir(const std::string& prefix) {
  const auto uniqueName =
      prefix + "_" +
      std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
  const auto dir = std::filesystem::temp_directory_path() / uniqueName;
  std::filesystem::create_directories(dir);
  return dir;
}

int runWithConfigFile(const std::filesystem::path& cfgPath) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "--config";
  std::string argv2 = cfgPath.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data()};
  int argc = 3;
  return KeplerFormalMain(argc, argv);
}

int runWithArgs(std::vector<std::string> args) {
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (auto& arg : args) {
    argv.push_back(arg.data());
  }
  return KeplerFormalMain(static_cast<int>(argv.size()), argv.data());
}

std::filesystem::path findBuiltNajaModuleDir() {
  const auto root = repoRoot();
  const std::vector<std::filesystem::path> candidates = {
      root / "build/thirdparty/naja/src/nl/python/naja_wrapping",
      root / "buildD/thirdparty/naja/src/nl/python/naja_wrapping",
      root / "buildR/thirdparty/naja/src/nl/python/naja_wrapping",
  };
  for (const auto& candidate : candidates) {
    if (std::filesystem::exists(candidate / "naja.so")) {
      return candidate;
    }
  }
  return {};
}

struct MultiFileVerilogFixture {
  std::filesystem::path tmpDir;
  std::filesystem::path cfgPath;
};

MultiFileVerilogFixture createVerilogPreprocessingFixture(bool enablePreprocessing) {
  MultiFileVerilogFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_preproc_v");

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
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_default_nettype_v");

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

struct SystemVerilogFlistFixture {
  std::filesystem::path tmpDir;
  std::filesystem::path design0ChildPath;
  std::filesystem::path design0TopPath;
  std::filesystem::path design1ChildPath;
  std::filesystem::path design1TopPath;
  std::filesystem::path design0FlistPath;
  std::filesystem::path design1FlistPath;
};

struct ScopedNajaIfFixture {
  std::filesystem::path tmpDir;
  std::filesystem::path design0IfPath;
  std::filesystem::path design1IfPath;
  std::filesystem::path libertyPath;
};

struct SequentialNajaIfFixture {
  std::filesystem::path tmpDir;
  std::filesystem::path design0IfPath;
  std::filesystem::path design1IfPath;
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

struct SolverGuard {
  SolverGuard(): oldValue_(KEPLER_FORMAL::Config::getSolverType()) {}

  ~SolverGuard() {
    KEPLER_FORMAL::Config::setSolverType(oldValue_);
  }

  KEPLER_FORMAL::Config::SolverType oldValue_;
};

struct ReportSkippedPOsGuard {
  ReportSkippedPOsGuard(): oldValue_(KEPLER_FORMAL::Config::getReportSkippedPOs()) {}

  ~ReportSkippedPOsGuard() {
    KEPLER_FORMAL::Config::setReportSkippedPOs(oldValue_);
  }

  bool oldValue_;
};

struct SecBoundaryAbstractionGuard {
  SecBoundaryAbstractionGuard()
      : oldValue_(
            KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary()) {}

  ~SecBoundaryAbstractionGuard() {
    KEPLER_FORMAL::Config::setSecTreatUncomputableSeqAsBoundary(oldValue_);
  }

  bool oldValue_;
};

struct CurrentPathGuard {
  CurrentPathGuard(): oldPath_(std::filesystem::current_path()) {}

  ~CurrentPathGuard() {
    std::error_code ec;
    std::filesystem::current_path(oldPath_, ec);
  }

  std::filesystem::path oldPath_;
};

std::vector<std::filesystem::path> listTemporarySystemVerilogCommandFiles() {
  std::vector<std::filesystem::path> files;
  std::error_code ec;
  const auto tempDir = std::filesystem::temp_directory_path(ec);
  if (ec) {
    return files;
  }
  for (const auto& entry : std::filesystem::directory_iterator(tempDir, ec)) {
    if (ec) {
      break;
    }
    if (!entry.is_regular_file()) {
      continue;
    }
    const auto name = entry.path().filename().string();
    if (name.rfind("kepler_formal_sv_top_", 0) == 0 &&
        entry.path().extension() == ".f") {
      files.push_back(entry.path().filename());
    }
  }
  std::sort(files.begin(), files.end());
  return files;
}

std::vector<std::filesystem::path> listMiterLogsInCurrentDirectory() {
  std::vector<std::filesystem::path> files;
  std::error_code ec;
  for (const auto& entry : std::filesystem::directory_iterator(std::filesystem::current_path(), ec)) {
    if (ec) {
      break;
    }
    if (!entry.is_regular_file()) {
      continue;
    }
    const auto name = entry.path().filename().string();
    if (name.rfind("miter_log_", 0) == 0 && entry.path().extension() == ".txt") {
      files.push_back(entry.path().filename());
    }
  }
  std::sort(files.begin(), files.end());
  return files;
}

MultiFileVerilogFixture createMultiFileVerilogFixture() {
  MultiFileVerilogFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_multi_v");

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

SimpleCliFixture createDesignFixture(const std::string& extension,
                                     const std::string& design0Body,
                                     const std::string& design1Body) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_simple");

  fixture.design0Path = fixture.tmpDir / ("design0." + extension);
  fixture.design1Path = fixture.tmpDir / ("design1." + extension);

  {
    std::ofstream design0(fixture.design0Path);
    design0 << design0Body;
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << design1Body;
  }

  return fixture;
}

SimpleCliFixture createEquivalentDesignFixture(const std::string& extension,
                                               const std::string& moduleBody) {
  return createDesignFixture(extension, moduleBody, moduleBody);
}

std::filesystem::path copyExampleLibertyFile(const std::filesystem::path& directory,
                                             const std::string& filename) {
  const auto source = repoRoot() / "example" / "NangateOpenCellLibrary_typical.lib";
  const auto destination = directory / filename;
  std::filesystem::copy_file(
      source, destination, std::filesystem::copy_options::overwrite_existing);
  return destination;
}

std::string readFileContents(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::ostringstream contents;
  contents << file.rdbuf();
  return contents.str();
}

SystemVerilogFlistFixture createSystemVerilogFlistFixture() {
  SystemVerilogFlistFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv_flist");

  fixture.design0ChildPath = fixture.tmpDir / "design0_child.sv";
  fixture.design0TopPath = fixture.tmpDir / "design0_top.sv";
  fixture.design1ChildPath = fixture.tmpDir / "design1_child.sv";
  fixture.design1TopPath = fixture.tmpDir / "design1_top.sv";
  fixture.design0FlistPath = fixture.tmpDir / "design0.f";
  fixture.design1FlistPath = fixture.tmpDir / "design1.f";

  {
    std::ofstream child(fixture.design0ChildPath);
    child << "module leaf(input logic a, output logic y);\n";
    child << "  assign y = a;\n";
    child << "endmodule\n";
  }
  {
    std::ofstream top(fixture.design0TopPath);
    top << "module cva6(input logic a, output logic y);\n";
    top << "  leaf u_leaf(.a(a), .y(y));\n";
    top << "endmodule\n";
  }
  {
    std::ofstream child(fixture.design1ChildPath);
    child << "module leaf(input logic a, output logic y);\n";
    child << "  assign y = a;\n";
    child << "endmodule\n";
  }
  {
    std::ofstream top(fixture.design1TopPath);
    top << "module cva6(input logic a, output logic y);\n";
    top << "  leaf u_leaf(.a(a), .y(y));\n";
    top << "endmodule\n";
  }
  {
    std::ofstream flist(fixture.design0FlistPath);
    flist << fixture.design0ChildPath.string() << "\n";
    flist << fixture.design0TopPath.string() << "\n";
  }
  {
    std::ofstream flist(fixture.design1FlistPath);
    flist << fixture.design1ChildPath.string() << "\n";
    flist << fixture.design1TopPath.string() << "\n";
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
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_scope_if");
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

SequentialNajaIfFixture createEquivalentSequentialNajaIfFixture(
    const std::string& ffName0 = "ff0",
    const std::string& ffName1 = "ff0") {
  SequentialNajaIfFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_seq_if");
  fixture.design0IfPath = fixture.tmpDir / "design0.capnp";
  fixture.design1IfPath = fixture.tmpDir / "design1.capnp";

  const auto dumpDesign = [&](const std::filesystem::path& dumpPath,
                              const std::string& ffName) {
    cleanupNajaTestState();
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    auto* designLibrary = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("DESIGN"));
    auto* top = SNLDesign::create(designLibrary, SNLDesign::Type::Standard, NLName("top"));
    auto* topIn = SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
    auto* topClock = SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
    auto* topOut = SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
    auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName(ffName));

    auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
    auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
    auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

    topIn->setNet(netIn);
    topClock->setNet(netClock);
    topOut->setNet(netQ);

    ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
    ff->getInstTerm(NLDB0::getDFFData())->setNet(netIn);
    ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

    db->setTopDesign(top);
    SNLCapnP::dump(db, dumpPath);
    cleanupNajaTestState();
  };

  dumpDesign(fixture.design0IfPath, ffName0);
  dumpDesign(fixture.design1IfPath, ffName1);

  return fixture;
}

SequentialNajaIfFixture createDifferentSequentialNajaIfFixture() {
  SequentialNajaIfFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_seq_if_diff");
  fixture.design0IfPath = fixture.tmpDir / "design0.capnp";
  fixture.design1IfPath = fixture.tmpDir / "design1.capnp";

  const auto dumpDesign = [&](const std::filesystem::path& dumpPath, bool invertData) {
    cleanupNajaTestState();
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    auto* primitiveLibrary =
        NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
    auto* designLibrary =
        NLLibrary::create(db, NLLibrary::Type::Standard, NLName("DESIGN"));
    auto* top =
        SNLDesign::create(designLibrary, SNLDesign::Type::Standard, NLName("top"));
    auto* topIn = SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
    auto* topClock =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
    auto* topOut =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
    auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));

    auto* invModel =
        SNLDesign::create(primitiveLibrary, SNLDesign::Type::Primitive, NLName("INV"));
    auto* invIn =
        SNLScalarTerm::create(invModel, SNLTerm::Direction::Input, NLName("A"));
    auto* invOut =
        SNLScalarTerm::create(invModel, SNLTerm::Direction::Output, NLName("Y"));
    SNLDesignModeling::addCombinatorialArcs({invIn}, {invOut});
    SNLDesignModeling::setTruthTable(invModel, SNLTruthTable::Inv());

    auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
    auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
    auto* netData = SNLScalarNet::create(top, NLName("net_data"));
    auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

    topIn->setNet(netIn);
    topClock->setNet(netClock);
    topOut->setNet(netQ);

    if (invertData) {
      auto* inv = SNLInstance::create(top, invModel, NLName("inv0"));
      inv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netIn);
      inv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netData);
    }

    ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
    ff->getInstTerm(NLDB0::getDFFData())->setNet(invertData ? netData : netIn);
    ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

    db->setTopDesign(top);
    SNLCapnP::dump(db, dumpPath);
    cleanupNajaTestState();
  };

  dumpDesign(fixture.design0IfPath, false);
  dumpDesign(fixture.design1IfPath, true);

  return fixture;
}

SequentialNajaIfFixture createUncomputableSequentialNajaIfFixture() {
  SequentialNajaIfFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_seq_if_uncomputable");
  fixture.design0IfPath = fixture.tmpDir / "design0.capnp";
  fixture.design1IfPath = fixture.tmpDir / "design1.capnp";

  const auto dumpDesign = [&](const std::filesystem::path& dumpPath) {
    cleanupNajaTestState();
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    auto* designLibrary =
        NLLibrary::create(db, NLLibrary::Type::Standard, NLName("DESIGN"));
    auto* primitiveLibrary =
        NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
    auto* unsupportedModel =
        SNLDesign::create(primitiveLibrary, SNLDesign::Type::Primitive, NLName("SEQ_NO_D"));
    auto* clock =
        SNLScalarTerm::create(unsupportedModel, SNLTerm::Direction::Input, NLName("CK"));
    auto* output =
        SNLScalarTerm::create(unsupportedModel, SNLTerm::Direction::Output, NLName("Q"));
    SNLDesignModeling::addClockToOutputsArcs(clock, {output});

    auto* top =
        SNLDesign::create(designLibrary, SNLDesign::Type::Standard, NLName("top"));
    auto* topClock =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
    auto* topOut =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
    auto* seq = SNLInstance::create(top, unsupportedModel, NLName("ff0"));
    auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
    auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

    topClock->setNet(netClock);
    topOut->setNet(netOut);
    seq->getInstTerm(unsupportedModel->getScalarTerm(NLName("CK")))->setNet(netClock);
    seq->getInstTerm(unsupportedModel->getScalarTerm(NLName("Q")))->setNet(netOut);

    db->setTopDesign(top);
    SNLCapnP::dump(db, dumpPath);
    cleanupNajaTestState();
  };

  dumpDesign(fixture.design0IfPath);
  dumpDesign(fixture.design1IfPath);
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
  if (!std::filesystem::exists(cnfPath)) {
    GTEST_SKIP() << "CNF dump file was not produced by the in-process run in this build";
  }

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

TEST(KeplerFormalCliTests, CliCompactFlagAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = fixture.design0Path.string();
  std::string argv3 = fixture.design1Path.string();
  std::string argv4 = "--compact";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data()};
  int argc = 5;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliCompactFlagWritesIdenticalSummaryToDefaultLog) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto runDir = fixture.tmpDir / "compact_cli_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);

    std::string argv0 = "kepler-formal";
    std::string argv1 = "-verilog";
    std::string argv2 = fixture.design0Path.string();
    std::string argv3 = fixture.design1Path.string();
    std::string argv4 = "--compact";
    char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                    argv4.data()};
    int argc = 5;

    EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);

    const auto logs = listMiterLogsInCurrentDirectory();
    ASSERT_EQ(logs.size(), 1u);
    const auto contents = readFileContents(runDir / logs.front());
    EXPECT_NE(contents.find("Circuits are IDENTICAL"), std::string::npos);
    EXPECT_EQ(contents.find("Due to compact mode, per PO analysis is skipped."),
              std::string::npos);
  }

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliCompactFlagAlignsReorderedInputsAndOutputs) {
  const auto fixture = createDesignFixture(
      "v",
      "module top(input a, input b, output y0, output y1);\n"
      "  assign y0 = a;\n"
      "  assign y1 = b;\n"
      "endmodule\n",
      "module top(input b, input a, output y1, output y0);\n"
      "  assign y1 = b;\n"
      "  assign y0 = a;\n"
      "endmodule\n");
  const auto runDir = fixture.tmpDir / "compact_reordered_cli_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);

    std::string argv0 = "kepler-formal";
    std::string argv1 = "-verilog";
    std::string argv2 = fixture.design0Path.string();
    std::string argv3 = fixture.design1Path.string();
    std::string argv4 = "--compact";
    char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                    argv4.data()};
    int argc = 5;

    EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);

    const auto logs = listMiterLogsInCurrentDirectory();
    ASSERT_EQ(logs.size(), 1u);
    const auto contents = readFileContents(runDir / logs.front());
    EXPECT_NE(contents.find("Circuits are IDENTICAL"), std::string::npos);
  }

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliReportSkippedPOsFlagAccepted) {
  ReportSkippedPOsGuard reportGuard;
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = fixture.design0Path.string();
  std::string argv3 = fixture.design1Path.string();
  std::string argv4 = "--report-skipped-pos";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data()};
  int argc = 5;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getReportSkippedPOs());
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigMissingInputPathsFails) {
  const auto cfgPath = writeTempConfig("format: verilog\nlog_level: info\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigScalarLibertyFilesIsIgnored) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "liberty_files: ignored_scalar_value\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigRootSequenceFallsBackToValidationFailure) {
  const auto cfgPath = writeTempConfig(
      "- format\n"
      "- verilog\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
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

TEST(KeplerFormalCliTests, ConfigSystemVerilogFlistAndTopAccepted) {
  const auto fixture = createSystemVerilogFlistFixture();
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "sv_design1_flist: " + fixture.design0FlistPath.string() + "\n"
      "sv_design2_flist: " + fixture.design1FlistPath.string() + "\n"
      "sv_design1_top: cva6\n"
      "sv_design2_top: cva6\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigCompactSystemVerilogFlistWithLibertyAndCnfAccepted) {
  const auto fixture = createSystemVerilogFlistFixture();
  const auto cnfPath = fixture.tmpDir / "compact_sv.cnf";
  const auto libertyPath = repoRoot() / "example" / "NangateOpenCellLibrary_typical.lib";
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "sv_design1_flist: " + fixture.design0FlistPath.string() + "\n"
      "sv_design2_flist: " + fixture.design1FlistPath.string() + "\n"
      "sv_design1_top: cva6\n"
      "sv_design2_top: cva6\n"
      "compact_mode: true\n"
      "cnf_export: true\n"
      "cnf_export_path: " + cnfPath.string() + "\n"
      "liberty_files:\n"
      "  - " + libertyPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_TRUE(std::filesystem::exists(cnfPath));

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSystemVerilogFlistMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "sv_design1_flist:\n"
      "  - bad\n"
      "sv_design2_flist: design1.f\n");
  EXPECT_NE(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigSystemVerilogSecondFlistMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "sv_design1_flist: design0.f\n"
      "sv_design2_flist:\n"
      "  - bad\n");
  EXPECT_NE(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigSystemVerilogTopMustNotBeEmpty) {
  const auto fixture = createSystemVerilogFlistFixture();
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "sv_design1_flist: " + fixture.design0FlistPath.string() + "\n"
      "sv_design2_flist: " + fixture.design1FlistPath.string() + "\n"
      "sv_design1_top: \"\"\n"
      "sv_design2_top: cva6\n");
  EXPECT_NE(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSystemVerilogFlistAndTopAccepted) {
  const auto fixture = createSystemVerilogFlistFixture();
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_flist";
  std::string argv3 = fixture.design0FlistPath.string();
  std::string argv4 = "--sv_design1_top";
  std::string argv5 = "cva6";
  std::string argv6 = "--sv_design2_flist";
  std::string argv7 = fixture.design1FlistPath.string();
  std::string argv8 = "--sv_design2_top";
  std::string argv9 = "cva6";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data(),
                  argv5.data(), argv6.data(), argv7.data(), argv8.data(), argv9.data()};
  int argc = 10;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSystemVerilogFlagMissingValueFails) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_flist";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data()};
  int argc = 3;
  EXPECT_NE(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
}

TEST(KeplerFormalCliTests, CliSystemVerilogEmptyValueFails) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_top";
  std::string argv3;
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data()};
  int argc = 4;
  EXPECT_NE(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
}

TEST(KeplerFormalCliTests, CliSystemVerilogRequiresSourcesForBothDesigns) {
  const auto fixture = createSystemVerilogFlistFixture();
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_flist";
  std::string argv3 = fixture.design0FlistPath.string();
  std::string argv4 = "--sv_design1_top";
  std::string argv5 = "cva6";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data(),
                  argv5.data()};
  int argc = 6;
  EXPECT_NE(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSystemVerilogOptionsRejectedForVerilogFormat) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--sv_design1_top";
  std::string argv3 = "cva6";
  std::string argv4 = "design0.v";
  std::string argv5 = "design1.v";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data(),
                  argv5.data()};
  int argc = 6;
  EXPECT_NE(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
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

TEST(KeplerFormalCliTests, ConfigGlucoseSolverAccepted) {
  SolverGuard solverGuard;
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "solver: glucose\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_EQ(KEPLER_FORMAL::Config::getSolverType(), KEPLER_FORMAL::Config::GLUCOSE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigKissatSolverAccepted) {
  SolverGuard solverGuard;
  KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::GLUCOSE);
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "solver: kissat\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_EQ(KEPLER_FORMAL::Config::getSolverType(), KEPLER_FORMAL::Config::KISSAT);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigCompactModeAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "compact_mode: true\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigCompactModeWritesDifferentSummaryAndWarning) {
  const auto fixture = createDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n",
      "module top(input a, output y);\n"
      "  assign y = 1'b0;\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "compact_config.log";
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "compact_mode: true\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("Circuits are DIFFERENT"), std::string::npos);
  EXPECT_NE(contents.find("Due to compact mode, per PO analysis is skipped."),
            std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigReportSkippedPOsAccepted) {
  ReportSkippedPOsGuard reportGuard;
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "report_skipped_pos: true\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getReportSkippedPOs());
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
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

TEST(KeplerFormalCliTests, ConfigCompactNajaIfAccepted) {
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto design = exampleDir / "tinyrocket_naja.if";
  const auto lib0 = exampleDir / "NangateOpenCellLibrary_typical.lib";
  const auto lib1 = exampleDir / "fakeram45_1024x32.lib";
  const auto lib2 = exampleDir / "fakeram45_64x32.lib";

  ASSERT_TRUE(std::filesystem::exists(design));
  ASSERT_TRUE(std::filesystem::exists(lib0));
  ASSERT_TRUE(std::filesystem::exists(lib1));
  ASSERT_TRUE(std::filesystem::exists(lib2));

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + design.string() + "\n"
      "  - " + design.string() + "\n"
      "liberty_files:\n"
      "  - " + lib0.string() + "\n"
      "  - " + lib1.string() + "\n"
      "  - " + lib2.string() + "\n"
      "compact_mode: true\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
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

TEST(KeplerFormalCliTests, ConfigInvalidVerificationModeFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: BAD\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigInvalidSecEngineFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: SEC\n"
      "sec_engine: BAD\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigExplicitLecVerificationAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: LEC\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigVerificationMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification:\n"
      "  - SEC\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigMaxKMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: SEC\n"
      "max_k:\n"
      "  - 4\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigInvalidMaxKTokenFails) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: nope\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigEmptyMaxKTokenFails) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: \"\"\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigOutOfRangeMaxKTokenFails) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 999999999999999999999999999999999999\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigMaxKWithoutSecFails) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "max_k: 5\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecEngineWithoutSecFails) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: LEC\n"
      "sec_engine: PDR\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecVerificationAccepted) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecVerificationAcceptedWithPdrEngine) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "sec_engine: PDR\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecVerificationAcceptedWithKInductionEngine) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "sec_engine: KINDUCTION\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecVerificationAcceptedWithImcEngine) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "sec_engine: IMC\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests,
     ConfigSecAbstractsUncomputableSequentialBoundariesByDefault) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createUncomputableSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 2\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary());
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests,
     ConfigSecCanDisableBoundaryAbstractionForUncomputableSequentials) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 2\n"
      "sec_uncomputable_seq_as_boundary: false\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_FALSE(KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary());
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecIgnoresRenamedInternalState) {
  const auto fixture =
      createEquivalentSequentialNajaIfFixture("state_a", "state_b");
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSystemVerilogSecVerificationAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(\n"
      "    input logic clk,\n"
      "    input logic rst,\n"
      "    input logic d,\n"
      "    output logic q\n"
      ");\n"
      "  always_ff @(posedge clk)\n"
      "  if (rst) begin\n"
      "    q <= 1'b0;\n"
      "  end else begin\n"
      "    q <= d;\n"
      "  end\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecVerificationWritesDefaultLog) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(\n"
      "    input logic clk,\n"
      "    input logic rst,\n"
      "    input logic d,\n"
      "    output logic q\n"
      ");\n"
      "  always_ff @(posedge clk)\n"
      "  if (rst) begin\n"
      "    q <= 1'b0;\n"
      "  end else begin\n"
      "    q <= d;\n"
      "  end\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  const auto runDir = fixture.tmpDir / "sec_log_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);

    const auto logs = listMiterLogsInCurrentDirectory();
    ASSERT_EQ(logs.size(), 1u);
    const auto contents = readFileContents(runDir / logs.front());
    EXPECT_NE(contents.find("Verification: SEC"), std::string::npos);
    EXPECT_NE(contents.find("SEC proved equivalence at k = 1"), std::string::npos);
  }

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecReportsPartialObservedOutputCoverage) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(\n"
      "    input logic clk,\n"
      "    input logic in,\n"
      "    output logic good,\n"
      "    output logic bad\n"
      ");\n"
      "  logic d;\n"
      "  logic q;\n"
      "  assign good = in;\n"
      "  always_ff @(posedge clk) begin\n"
      "    if (1'b0) begin\n"
      "      q <= 1'b0;\n"
      "    end else begin\n"
      "      q <= d;\n"
      "    end\n"
      "  end\n"
      "  assign bad = q;\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sec_partial_coverage.log";
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: SEC\n"
      "max_k: 2\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);

  const auto contents = readFileContents(logPath);
  EXPECT_NE(
      contents.find("SEC output coverage: 50.00% (1/2 covered/existing outputs)."),
      std::string::npos);
  EXPECT_NE(
      contents.find("SEC skipped observed outputs due to connectivity issues"),
      std::string::npos);
  EXPECT_NE(contents.find("bad[0]"), std::string::npos);
  EXPECT_NE(contents.find("no-driver connectivity"), std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecDifferenceLogIncludesWitnessDetails) {
  const auto fixture = createDifferentSequentialNajaIfFixture();
  const auto logPath = fixture.tmpDir / "sec_difference.log";
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 2\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("SEC counterexample details:"), std::string::npos);
  EXPECT_NE(contents.find("cycle 1"), std::string::npos);
  EXPECT_NE(contents.find("Input trace:"), std::string::npos);
  EXPECT_NE(contents.find("in[0]"), std::string::npos);
  EXPECT_NE(contents.find("out[0]"), std::string::npos);
  EXPECT_NE(contents.find("Traceback for first differing point `out[0]` at cycle 1:"),
            std::string::npos);
  EXPECT_NE(contents.find("design0 cone to environment inputs:"), std::string::npos);
  EXPECT_NE(contents.find("design1 cone to environment inputs:"), std::string::npos);
  EXPECT_NE(contents.find("cone terms only in design1: inv0.Y[0]"),
            std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigTinyRocketSecVerificationAccepted) {
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto design = exampleDir / "tinyrocket.v";
  const auto lib0 = exampleDir / "NangateOpenCellLibrary_typical.lib";
  const auto lib1 = exampleDir / "fakeram45_1024x32.lib";
  const auto lib2 = exampleDir / "fakeram45_64x32.lib";
  const auto lib3 = exampleDir / "fakeram45_64x15.lib";

  ASSERT_TRUE(std::filesystem::exists(design));
  ASSERT_TRUE(std::filesystem::exists(lib0));
  ASSERT_TRUE(std::filesystem::exists(lib1));
  ASSERT_TRUE(std::filesystem::exists(lib2));
  ASSERT_TRUE(std::filesystem::exists(lib3));

  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: SEC\n"
      "max_k: 1\n"
      "input_paths:\n"
      "  - " + design.string() + "\n"
      "  - " + design.string() + "\n"
      "liberty_files:\n"
      "  - " + lib0.string() + "\n"
      "  - " + lib1.string() + "\n"
      "  - " + lib2.string() + "\n"
      "  - " + lib3.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
}

TEST(KeplerFormalCliTests, ConfigSecRejectsCompactMode) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "compact_mode: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecRejectsScopeExtraction) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "use_scopes: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecRejectsCnfExport) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "cnf_export: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecRejectsSkippedPoReporting) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "report_skipped_pos: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecFallsBackWhenLogParentCannotBeCreated) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(\n"
      "    input logic clk,\n"
      "    input logic rst,\n"
      "    input logic d,\n"
      "    output logic q\n"
      ");\n"
      "  always_ff @(posedge clk)\n"
      "  if (rst) begin\n"
      "    q <= 1'b0;\n"
      "  end else begin\n"
      "    q <= d;\n"
      "  end\n"
      "endmodule\n");
  const auto blockedParent = fixture.tmpDir / "blocked_parent";
  {
    std::ofstream blocked(blockedParent);
    blocked << "not a directory\n";
  }
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "log_file: " + (blockedParent / "sec.log").string() + "\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecContinuesWhenLogFilePathIsDirectory) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(\n"
      "    input logic clk,\n"
      "    input logic rst,\n"
      "    input logic d,\n"
      "    output logic q\n"
      ");\n"
      "  always_ff @(posedge clk)\n"
      "  if (rst) begin\n"
      "    q <= 1'b0;\n"
      "  end else begin\n"
      "    q <= d;\n"
      "  end\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: SEC\n"
      "max_k: 4\n"
      "log_file: " + fixture.tmpDir.string() + "\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSecVerificationAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-v";
  std::string argv2 = "SEC";
  std::string argv3 = "-k";
  std::string argv4 = "4";
  std::string argv5 = "-naja_if";
  std::string argv6 = fixture.design0IfPath.string();
  std::string argv7 = fixture.design1IfPath.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSecEngineAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "SEC",
                   "-k",
                   "4",
                   "--sec-engine",
                   "PDR",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliKInductionSecEngineAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "SEC",
                   "-k",
                   "4",
                   "--sec-engine",
                   "KINDUCTION",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliImcSecEngineAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "SEC",
                   "-k",
                   "4",
                   "--sec-engine",
                   "IMC",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliExplicitLecVerificationAcceptedBeforeFormat) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "LEC",
                   "-verilog",
                   fixture.design0Path.string(),
                   fixture.design1Path.string()}),
      EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliMissingVerificationBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-v"}), EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliInvalidVerificationBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-v", "BAD"}), EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliMissingMaxKBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-k"}), EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliInvalidMaxKBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-k", "-1"}), EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliOutOfRangeMaxKBeforeFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-k", "999999999999999999999999999999999999"}),
      EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliUnrecognizedOptionBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "--mystery"}), EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliMissingInputFormatAfterPreOptionsFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-v", "SEC", "-k", "4"}), EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliInvalidSecEngineBeforeFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-v", "SEC", "--sec-engine", "BAD", "-verilog"}),
      EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliMissingVerificationAfterFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-verilog", "--verification"}), EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliInvalidVerificationAfterFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", "--verification", "BAD"}),
      EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliMissingMaxKAfterFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-verilog", "--max-k"}), EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliInvalidMaxKAfterFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", "--max-k", "oops"}),
      EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliOutOfRangeMaxKAfterFormatFails) {
  EXPECT_EQ(
      runWithArgs(
          {"kepler-formal", "-verilog", "--max-k", "999999999999999999999999999999999999"}),
      EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, ConfigSecInconclusiveFails) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(\n"
      "    input logic clk,\n"
      "    input logic rst,\n"
      "    input logic d,\n"
      "    output logic q\n"
      ");\n"
      "  always_ff @(posedge clk)\n"
      "  if (rst) begin\n"
      "    q <= 1'b0;\n"
      "  end else begin\n"
      "    q <= d;\n"
      "  end\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: SEC\n"
      "max_k: 0\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, ConfigSecUnsupportedMismatchFails) {
  const auto fixture = createDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  reg q;\n"
      "  always @(*) q = a;\n"
      "  assign y = q;\n"
      "endmodule\n",
      "module top(input a, output z);\n"
      "  reg q;\n"
      "  always @(*) q = a;\n"
      "  assign z = q;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: SEC\n"
      "max_k: 1\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
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

TEST(KeplerFormalCliTests, CliLibAliasCollectsPaths) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--design1";
  std::string argv3 = "a.v";
  std::string argv4 = "--design2";
  std::string argv5 = "b.v";
  std::string argv6 = "--lib";
  std::string argv7 = "lib1.lib";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
}

TEST(KeplerFormalCliTests, CliPositionalLibertyPathsAfterNetlistsAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto gzLiberty =
      repoRoot() / "thirdparty/naja/test/nl/formats/liberty/benchmarks/tests/small.lib.gz";
  ASSERT_TRUE(std::filesystem::exists(gzLiberty));

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = fixture.design0Path.string();
  std::string argv3 = fixture.design1Path.string();
  std::string argv4 = gzLiberty.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data()};
  int argc = 5;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliPositionalLibertyPathsWithoutStandardSuffixAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto renamedLiberty = copyExampleLibertyFile(fixture.tmpDir, "primitives.cells");

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = fixture.design0Path.string();
  std::string argv3 = fixture.design1Path.string();
  std::string argv4 = renamedLiberty.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data()};
  int argc = 5;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, PythonLibraryFilesAreLoadedFromDedicatedYamlKey) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(output y);\n"
      "  wire z;\n"
      "  LOGIC1 c0(.Z(z));\n"
      "  BUF c1(.A(z), .Z(y));\n"
      "endmodule\n");
  const auto pyPrimitives =
      repoRoot() / "thirdparty/naja/test/nl/python/pyloader/scripts/primitives1.py";
  const auto pyModuleDir = findBuiltNajaModuleDir();
  ASSERT_TRUE(std::filesystem::exists(pyPrimitives));
  ASSERT_FALSE(pyModuleDir.empty());
  ASSERT_TRUE(std::filesystem::exists(pyModuleDir / "naja.so"));
  EnvVarGuard pythonPathGuard("PYTHONPATH");
  pythonPathGuard.set(pyModuleDir.string());

  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "py_tech_files:\n"
      "  - " + pyPrimitives.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, PythonLibraryFilesUnderLibertyFilesAreRejected) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(output y);\n"
      "  wire z;\n"
      "  LOGIC1 c0(.Z(z));\n"
      "  BUF c1(.A(z), .Z(y));\n"
      "endmodule\n");
  const auto pyPrimitives =
      repoRoot() / "thirdparty/naja/test/nl/python/pyloader/scripts/primitives1.py";
  ASSERT_TRUE(std::filesystem::exists(pyPrimitives));

  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "liberty_files:\n"
      "  - " + pyPrimitives.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
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

TEST(KeplerFormalCliTests, InvalidLibertyFileContentFails) {
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

TEST(KeplerFormalCliTests, CliPythonLibraryFilesAreRejected) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(output y);\n"
      "  wire z;\n"
      "  LOGIC1 c0(.Z(z));\n"
      "  BUF c1(.A(z), .Z(y));\n"
      "endmodule\n");
  const auto pyPrimitives =
      repoRoot() / "thirdparty/naja/test/nl/python/pyloader/scripts/primitives1.py";
  ASSERT_TRUE(std::filesystem::exists(pyPrimitives));

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = fixture.design0Path.string();
  std::string argv3 = fixture.design1Path.string();
  std::string argv4 = pyPrimitives.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data()};
  int argc = 5;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSystemVerilogFlistWithoutTopAccepted) {
  const auto fixture = createSystemVerilogFlistFixture();
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_flist";
  std::string argv3 = fixture.design0FlistPath.string();
  std::string argv4 = "--sv_design2_flist";
  std::string argv5 = fixture.design1FlistPath.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data(),
                  argv5.data()};
  int argc = 6;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSystemVerilogDirectPathsWithEscapedNamesAccepted) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_cli_sv_escaped";
  std::filesystem::create_directories(tmpDir);
  const auto design0 = tmpDir / "design0\\escaped.sv";
  const auto design1 = tmpDir / "design1\\escaped.sv";
  {
    std::ofstream f(design0);
    f << "module top(input logic a, output logic y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1);
    f << "module top(input logic a, output logic y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_top";
  std::string argv3 = "top";
  std::string argv4 = "--sv_design2_top";
  std::string argv5 = "top";
  std::string argv6 = design0.string();
  std::string argv7 = design1.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(tmpDir);
}

TEST(KeplerFormalCliTests, CliSystemVerilogTopCommandFileCreationFailureFails) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic a, output logic y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto noWriteDir =
      std::filesystem::temp_directory_path() / "kepler_formal_sv_no_write";
  std::filesystem::create_directories(noWriteDir);
  auto perms = std::filesystem::status(noWriteDir).permissions();
  std::filesystem::permissions(
      noWriteDir,
      std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
          std::filesystem::perms::others_write,
      std::filesystem::perm_options::remove);

  EnvVarGuard tmpDirGuard("TMPDIR");
  tmpDirGuard.set(noWriteDir.string());

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_top";
  std::string argv3 = "top";
  std::string argv4 = "--sv_design2_top";
  std::string argv5 = "top";
  std::string argv6 = fixture.design0Path.string();
  std::string argv7 = fixture.design1Path.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);

  std::filesystem::permissions(noWriteDir, perms);
  std::filesystem::remove_all(noWriteDir);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliCompactSystemVerilogTopCommandFileCreationFailureFails) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic a, output logic y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto noWriteDir =
      std::filesystem::temp_directory_path() / "kepler_formal_sv_compact_no_write";
  std::filesystem::create_directories(noWriteDir);
  auto perms = std::filesystem::status(noWriteDir).permissions();
  std::filesystem::permissions(
      noWriteDir,
      std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
          std::filesystem::perms::others_write,
      std::filesystem::perm_options::remove);

  EnvVarGuard tmpDirGuard("TMPDIR");
  tmpDirGuard.set(noWriteDir.string());

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_top";
  std::string argv3 = "top";
  std::string argv4 = "--sv_design2_top";
  std::string argv5 = "top";
  std::string argv6 = "--compact";
  std::string argv7 = fixture.design0Path.string();
  std::string argv8 = fixture.design1Path.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data(),
                  argv8.data()};
  int argc = 9;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);

  std::filesystem::permissions(noWriteDir, perms);
  std::filesystem::remove_all(noWriteDir);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, CliSystemVerilogFirstDesignFailureCleansTemporaryCommandFile) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_sv_cleanup_first";
  std::filesystem::create_directories(tmpDir);
  const auto design0 = tmpDir / "design0_invalid.sv";
  const auto design1 = tmpDir / "design1_valid.sv";
  {
    std::ofstream f(design0);
    f << "module top(input logic a, output logic y)\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1);
    f << "module top(input logic a, output logic y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }

  const auto beforeTempFiles = listTemporarySystemVerilogCommandFiles();

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_top";
  std::string argv3 = "top";
  std::string argv4 = "--sv_design2_top";
  std::string argv5 = "top";
  std::string argv6 = design0.string();
  std::string argv7 = design1.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);
  EXPECT_EQ(listTemporarySystemVerilogCommandFiles(), beforeTempFiles);

  std::filesystem::remove_all(tmpDir);
}

TEST(KeplerFormalCliTests, CliSystemVerilogSecondDesignFailureCleansTemporaryCommandFile) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_sv_cleanup_second";
  std::filesystem::create_directories(tmpDir);
  const auto design0 = tmpDir / "design0_valid.sv";
  const auto design1 = tmpDir / "design1_invalid.sv";
  {
    std::ofstream f(design0);
    f << "module top(input logic a, output logic y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1);
    f << "module top(input logic a, output logic y)\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }

  const auto beforeTempFiles = listTemporarySystemVerilogCommandFiles();

  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_top";
  std::string argv3 = "top";
  std::string argv4 = "--sv_design2_top";
  std::string argv5 = "top";
  std::string argv6 = design0.string();
  std::string argv7 = design1.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);
  EXPECT_EQ(listTemporarySystemVerilogCommandFiles(), beforeTempFiles);

  std::filesystem::remove_all(tmpDir);
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

TEST(KeplerFormalCliTests, SnlScopesCreatesOnlyOneDefaultMiterLogPerScopeRun) {
  const auto fixture = createEquivalentScopedNajaIfFixture();
  const auto runDir = fixture.tmpDir / "run_dir";
  std::filesystem::create_directories(runDir);

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "liberty_files:\n"
      "  - " + fixture.libertyPath.string() + "\n"
      "use_scopes: true\n");

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);

    const auto logs = listMiterLogsInCurrentDirectory();
    EXPECT_EQ(logs.size(), 1u);
  }

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, SnlScopesCanCleanAndDumpCnf) {
  const auto fixture = createEquivalentScopedNajaIfFixture();
  const auto cnfPath = fixture.tmpDir / "scoped_miter.cnf";
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "liberty_files:\n"
      "  - " + fixture.libertyPath.string() + "\n"
      "use_scopes: true\n"
      "clean_scopes: true\n"
      "cnf_export: true\n"
      "cnf_export_path: " + cnfPath.string() + "\n");

  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  EXPECT_TRUE(std::filesystem::exists(cnfPath));

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST(KeplerFormalCliTests, SnlScopesDumpCnfUsesDefaultScopedPath) {
  const auto fixture = createEquivalentScopedNajaIfFixture();
  const auto defaultCnfPath = std::filesystem::current_path() / "miter_child.cnf";
  std::filesystem::remove(defaultCnfPath);

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "liberty_files:\n"
      "  - " + fixture.libertyPath.string() + "\n"
      "use_scopes: true\n"
      "cnf_export: true\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_TRUE(std::filesystem::exists(defaultCnfPath));

  std::filesystem::remove(cfgPath);
  std::filesystem::remove(defaultCnfPath);
  std::filesystem::remove_all(fixture.tmpDir);
}
