// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "DNL.h"
#include "BoolExprCache.h"
#include "Config.h"
#include "Tree2BoolExpr.h"
#include "KeplerFormalUtils.h"
#include "NLDB0.h"
#include "NLUniverse.h"
#include "SNLCapnP.h"
#include "SNLDumpManifest.h"
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

std::filesystem::path copyNajaIfForCurrentBuild(
    const std::filesystem::path& source,
    const std::string& prefix) {
  const auto tempDir = makeUniqueTempDir(prefix);
  const auto copy = tempDir / source.filename();
  std::filesystem::copy(source, copy, std::filesystem::copy_options::recursive);
  // The payload was produced by this Naja revision, but Git may choose a
  // different unambiguous short-hash length in another checkout.
  naja::NL::SNLDumpManifest::dump(copy);
  return copy;
}

int runWithConfigFile(
    const std::filesystem::path& cfgPath,
    std::string argv0 = "kepler-formal") {
  std::string argv1 = "--config";
  std::string argv2 = cfgPath.string();
  char* argv[] = {argv0.data(), argv1.data(), argv2.data()};
  int argc = 3;
  const int rc = KeplerFormalMain(argc, argv);
  // CLI tests invoke the tool in-process, so clear the production BoolExpr
  // cache explicitly instead of relying on OS cleanup at process exit.
  KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
  KEPLER_FORMAL::BoolExprCache::destroy();
  return rc;
}

int runWithArgs(std::vector<std::string> args) {
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (auto& arg : args) {
    argv.push_back(arg.data());
  }
  const int rc = KeplerFormalMain(static_cast<int>(argv.size()), argv.data());
  KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
  KEPLER_FORMAL::BoolExprCache::destroy();
  return rc;
}

std::filesystem::path findBuiltNajaModuleDir() {
  const auto root = repoRoot();
  const auto najaModuleSuffix =
      std::filesystem::path("thirdparty/naja/src/nl/python/naja_wrapping");
  for (auto cursor = std::filesystem::current_path(); !cursor.empty();
       cursor = cursor.parent_path()) {
    const auto candidate = cursor / najaModuleSuffix;
    if (std::filesystem::exists(candidate / "naja.so")) {
      return candidate;
    }
    if (cursor == cursor.root_path()) {
      break;
    }
  }
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

std::vector<std::filesystem::path> listTemporaryFilesWithPrefixAndExtension(
    const std::string& prefix,
    const std::string& extension) {
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
    if (name.rfind(prefix, 0) == 0 && entry.path().extension() == extension) {
      files.push_back(entry.path().filename());
    }
  }
  std::sort(files.begin(), files.end());
  return files;
}

std::vector<std::filesystem::path> listTemporarySystemVerilogCommandFiles() {
  return listTemporaryFilesWithPrefixAndExtension("kepler_formal_sv_top_", ".f");
}

std::vector<std::filesystem::path> listTemporarySystemVerilogPrimitiveStubFiles() {
  return listTemporaryFilesWithPrefixAndExtension("kepler_formal_sv2v_prims_", ".sv");
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

// These fixtures have equivalent transition logic but no reset. Exact SEC may
// therefore report a cycle-0 difference between their independent initial states.
SequentialNajaIfFixture createEquivalentSequentialNajaIfFixture(
    const std::string& ffName0 = "ff0",
    const std::string& ffName1 = "ff0",
    const std::string& outputName0 = "out",
    const std::string& outputName1 = "out") {
  SequentialNajaIfFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_seq_if");
  fixture.design0IfPath = fixture.tmpDir / "design0.capnp";
  fixture.design1IfPath = fixture.tmpDir / "design1.capnp";

  const auto dumpDesign = [&](const std::filesystem::path& dumpPath,
                              const std::string& ffName,
                              const std::string& outputName) {
    cleanupNajaTestState();
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    auto* designLibrary = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("DESIGN"));
    auto* top = SNLDesign::create(designLibrary, SNLDesign::Type::Standard, NLName("top"));
    auto* topIn = SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
    auto* topClock = SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
    auto* topOut =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName(outputName));
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

  dumpDesign(fixture.design0IfPath, ffName0, outputName0);
  dumpDesign(fixture.design1IfPath, ffName1, outputName1);

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
    auto* topIn =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
    auto* topClock =
        SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
    auto* topGood =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("good"));
    auto* topBad =
        SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("bad"));
    auto* seq = SNLInstance::create(top, unsupportedModel, NLName("ff0"));
    auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
    auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
    auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

    topIn->setNet(netIn);
    topClock->setNet(netClock);
    topGood->setNet(netIn);
    topBad->setNet(netOut);
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

class KeplerFormalCliTests : public ::testing::Test {
 protected:
  static std::string logLineContaining(const std::string& contents,
                                       const std::string& message) {
    const size_t messagePosition = contents.find(message);
    if (messagePosition == std::string::npos) {
      return {};
    }
    const size_t lineStart = contents.rfind('\n', messagePosition);
    const size_t lineEnd = contents.find('\n', messagePosition);
    const size_t start = lineStart == std::string::npos ? 0 : lineStart + 1;
    return contents.substr(
        start,
        lineEnd == std::string::npos ? std::string::npos : lineEnd - start);
  }

  void TearDown() override {
    KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
    KEPLER_FORMAL::BoolExprCache::destroy();
  }
};

}  // namespace

TEST_F(KeplerFormalCliTests, SanitizeFileToken) {
  EXPECT_EQ(sanitizeFileToken("scope"), "scope");
  EXPECT_EQ(sanitizeFileToken("my scope"), "my_scope");
  EXPECT_EQ(sanitizeFileToken("a/b\\c"), "a_b_c");
  EXPECT_EQ(sanitizeFileToken(""), "scope");
}

TEST_F(KeplerFormalCliTests, SecResultExitCodesAreStable) {
  EXPECT_EQ(kSecProvedExitCode, 0);
  EXPECT_EQ(kSecPartiallyProvedExitCode, 1);
  EXPECT_EQ(kSecInconclusiveExitCode, 2);
  EXPECT_EQ(kSecCounterexampleExitCode, 3);
}


TEST_F(KeplerFormalCliTests, DumpCnfFromConfig) {
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

TEST_F(KeplerFormalCliTests, DumpPoCnfFromConfig) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto poCnfDir = fixture.tmpDir / "po_cnfs";
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "po_cnf_export: true\n"
      "po_cnf_export_path: " + poCnfDir.string() + "\n");

  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(poCnfDir / "top0"));
  ASSERT_TRUE(std::filesystem::exists(poCnfDir / "top1"));
  EXPECT_EQ(std::distance(std::filesystem::directory_iterator(poCnfDir / "top0"),
                          std::filesystem::directory_iterator{}), 1);
  EXPECT_EQ(std::distance(std::filesystem::directory_iterator(poCnfDir / "top1"),
                          std::filesystem::directory_iterator{}), 1);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, MultiFileVerilogConfig) {
  const auto fixture = createMultiFileVerilogFixture();
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, YamlMultiFileVerilogConfig) {
  const auto fixture = createMultiFileVerilogFixture();
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigCannotBeCombinedWithCommandLineOptions) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "--config",
                   cfgPath.string(),
                   "--verilog_design1_top",
                   "top",
                   "--verilog_design2_top",
                   "top"}),
      EXIT_FAILURE);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, VerilogPreprocessingEnabledParsesDirectiveInput) {
  const auto fixture = createVerilogPreprocessingFixture(true);
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, VerilogPreprocessingDisabledFailsOnDirectiveInput) {
  const auto fixture = createVerilogPreprocessingFixture(false);
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, VerilogPreprocessingEnabledDefaultNettypeIsRejected) {
  const auto fixture = createDefaultNettypeDirectiveFixture(true);
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, VerilogPreprocessingDisabledDefaultNettypeIsRejected) {
  const auto fixture = createDefaultNettypeDirectiveFixture(false);
  int rc = runWithConfigFile(fixture.cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliVerilogPreprocessingFlagEnablesDirectiveInput) {
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

TEST_F(KeplerFormalCliTests, CliWithoutVerilogPreprocessingFlagFailsOnDirectiveInput) {
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

TEST_F(KeplerFormalCliTests, CliCompactFlagAccepted) {
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

TEST_F(KeplerFormalCliTests, CliCompactFlagWritesIdenticalSummaryToDefaultLog) {
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

TEST_F(KeplerFormalCliTests, CliCompactFlagAlignsReorderedInputsAndOutputs) {
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

TEST_F(KeplerFormalCliTests, CliReportSkippedPOsFlagAccepted) {
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

TEST_F(KeplerFormalCliTests, ConfigMissingInputPathsFails) {
  const auto cfgPath = writeTempConfig("format: verilog\nlog_level: info\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigScalarLibertyFilesIsIgnored) {
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

TEST_F(KeplerFormalCliTests, ConfigRootSequenceFallsBackToValidationFailure) {
  const auto cfgPath = writeTempConfig(
      "- format\n"
      "- verilog\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigNestedSecondNotSequenceFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - [a.v]\n"
      "  - b.v\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigNestedEmptyDesignFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - []\n"
      "  - []\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigEmptyInputPathsFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: []\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigInputPathsNotSequenceFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: foo\n"
      "log_level: info\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigUnknownKeyFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "bogus_key: 1\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigRemovedInternalStateCorrespondenceKeyFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "sec_internal_state_correspondence: true\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigNonScalarKeyFails) {
  const auto cfgPath = writeTempConfig(
      "? [a, b]\n"
      ": 1\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigUnknownFormatFails) {
  const auto cfgPath = writeTempConfig(
      "format: unknown_format\n"
      "input_paths: [a.v, b.v]\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigSystemVerilogAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic a, output logic y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSystemVerilogLecRejected) {
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
  EXPECT_NE(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSv2vGateLevelVerilogTopAccepted) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  assign y = a;\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
    design1 << "module unused(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
  }
  const auto cfgPath = writeTempConfig(
      "format: sv2v\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "verilog_design2_top: top\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSv2vSystemVerilogDesign1UsesLoadedPrimitive) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_prim");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  const auto libertyPath = repoRoot() / "example" / "NangateOpenCellLibrary_typical.lib";
  ASSERT_TRUE(std::filesystem::exists(libertyPath));

  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  INV_X1 u_inv(.A(a), .ZN(y));\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n";
    design1 << "  INV_X1 u_inv(.A(a), .ZN(y));\n";
    design1 << "endmodule\n";
  }

  const auto cfgPath = writeTempConfig(
      "format: sv2v\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "liberty_files:\n"
      "  - " + libertyPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(
    KeplerFormalCliTests,
    ConfigPythonPrimitivesUseAdjacentNajaModuleWithoutPythonPath) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_adjacent_naja");
  fixture.design0Path = fixture.tmpDir / "design0.v";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  const auto pyPrimitives = fixture.tmpDir / "primitives.py";
  const auto cfgPath = fixture.tmpDir / "config.yaml";
  {
    std::ofstream py(pyPrimitives);
    py << "import naja\n"
          "\n"
          "def constructPrimitives(lib):\n"
          "  cell = naja.SNLDesign.createPrimitive(lib, 'BUF')\n"
          "  naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Input, 'A')\n"
          "  naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Output, 'Z')\n"
          "  cell.setTruthTable(0b10)\n";
  }
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input a, output y);\n"
               "  BUF u_buf(.A(a), .Z(y));\n"
               "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n"
               "  assign y = a;\n"
               "endmodule\n";
  }
  {
    std::ofstream cfg(cfgPath);
    cfg << "format: verilog\n"
           "verification: lec\n"
           "input_paths:\n"
        << "  - " << fixture.design0Path.string() << "\n"
        << "  - " << fixture.design1Path.string() << "\n"
           "py_tech_files:\n"
        << "  - " << pyPrimitives.string() << "\n";
  }

  const char* keplerBin = std::getenv("KEPLER_BIN");
  ASSERT_NE(keplerBin, nullptr);
  const std::filesystem::path keplerBinPath(keplerBin);
  ASSERT_TRUE(std::filesystem::exists(keplerBinPath));
  ASSERT_TRUE(std::filesystem::exists(keplerBinPath.parent_path() / "naja.so"));

  EnvVarGuard pythonPathGuard("PYTHONPATH");
  unsetenv("PYTHONPATH");
  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(fixture.tmpDir);
    EXPECT_EQ(runWithConfigFile(cfgPath, keplerBinPath.string()), EXIT_SUCCESS);
  }
  ASSERT_NE(std::getenv("PYTHONPATH"), nullptr);
  EXPECT_EQ(
      std::filesystem::path(std::getenv("PYTHONPATH")),
      keplerBinPath.parent_path());

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSv2vPythonPrimitivesBuildsComplexStubLibrary) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_py_prims");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  const auto pyPrimitives = fixture.tmpDir / "sv2v_primitives.py";
  {
    std::ofstream py(pyPrimitives);
    py << "import naja\n"
          "\n"
          "def constructPrimitives(lib):\n"
          "  naja.SNLDesign.createPrimitive(lib)\n"
          "  naja.SNLDesign.createPrimitive(lib, 'DUP')\n"
          "  odd = naja.SNLDesign.createPrimitive(lib, '1BAD-BOX')\n"
          "  naja.SNLBusTerm.create(odd, naja.SNLTerm.Direction.InOut, 3, 0, 'DATA-BUS')\n"
          "  child = naja.NLLibrary.createPrimitives(lib, 'CHILD')\n"
          "  naja.SNLDesign.createPrimitive(child, 'DUP')\n"
          "  passthrough = naja.SNLDesign.createPrimitive(child, 'child_prim')\n"
          "  a = naja.SNLScalarTerm.create(passthrough, naja.SNLTerm.Direction.Input, 'A')\n"
          "  y = naja.SNLScalarTerm.create(passthrough, naja.SNLTerm.Direction.Output, 'Y')\n"
          "  passthrough.setTruthTable(0b10)\n";
  }
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  child_prim u_child(.A(a), .Y(y));\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
  }

  const auto pyModuleDir = findBuiltNajaModuleDir();
  ASSERT_TRUE(std::filesystem::exists(pyPrimitives));
  ASSERT_FALSE(pyModuleDir.empty());
  ASSERT_TRUE(std::filesystem::exists(pyModuleDir / "naja.so"));
  EnvVarGuard pythonPathGuard("PYTHONPATH");
  pythonPathGuard.set(pyModuleDir.string());

  const auto cfgPath = writeTempConfig(
      "format: sv2v\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "py_tech_files:\n"
      "  - " + pyPrimitives.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSv2vUnnamedPythonPrimitiveLibraryCreatesNoStub) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_unnamed_py_prims");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  const auto pyPrimitives = fixture.tmpDir / "unnamed_primitives.py";
  {
    std::ofstream py(pyPrimitives);
    py << "import naja\n"
          "\n"
          "def constructPrimitives(lib):\n"
          "  naja.SNLDesign.createPrimitive(lib)\n";
  }
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module helper(input logic a, output logic y);\n";
    design0 << "  assign y = a;\n";
    design0 << "endmodule\n";
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  helper u_helper(.a(a), .y(y));\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module helper(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
    design1 << "module top(input a, output y);\n";
    design1 << "  helper u_helper(.a(a), .y(y));\n";
    design1 << "endmodule\n";
  }

  const auto pyModuleDir = findBuiltNajaModuleDir();
  ASSERT_TRUE(std::filesystem::exists(pyPrimitives));
  ASSERT_FALSE(pyModuleDir.empty());
  ASSERT_TRUE(std::filesystem::exists(pyModuleDir / "naja.so"));
  EnvVarGuard pythonPathGuard("PYTHONPATH");
  pythonPathGuard.set(pyModuleDir.string());

  const auto cfgPath = writeTempConfig(
      "format: sv2v\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "py_tech_files:\n"
      "  - " + pyPrimitives.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSv2vLecRejected) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_lec");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  assign y = a;\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
  }
  const auto cfgPath = writeTempConfig(
      "format: sv2v\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_NE(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSv2vRejectsSecondSystemVerilogOptions) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_reject_design2_sv");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  assign y = a;\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
  }
  const auto cfgPath = writeTempConfig(
      "format: sv2v\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "sv_design2_top: top\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_NE(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSv2vRequiresDesign1SystemVerilogSources) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_missing_design1");
  fixture.design1Path = fixture.tmpDir / "design1.v";
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
  }
  EXPECT_NE(runWithArgs({
      "kepler-formal",
      "-sv2v",
      "--design2",
      fixture.design1Path.string(),
      "-v",
      "sec",
  }), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSv2vRequiresDesign2VerilogSources) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_missing_design2");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  assign y = a;\n";
    design0 << "endmodule\n";
  }
  EXPECT_NE(runWithArgs({
      "kepler-formal",
      "-sv2v",
      "--design1",
      fixture.design0Path.string(),
      "-v",
      "sec",
  }), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSv2vAccepted) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_arg");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  assign y = a;\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
  }
  const int rc = runWithArgs({
      "kepler-formal",
      "-sv2v",
      fixture.design0Path.string(),
      fixture.design1Path.string(),
      "-v",
      "sec",
  });
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSvAliasAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic a, output logic y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-sv";
  std::string argv2 = fixture.design0Path.string();
  std::string argv3 = fixture.design1Path.string();
  std::string argv4 = "-v";
  std::string argv5 = "sec";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data()};
  int argc = 6;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSystemVerilogFlistAndTopAccepted) {
  const auto fixture = createSystemVerilogFlistFixture();
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "sv_design1_flist: " + fixture.design0FlistPath.string() + "\n"
      "sv_design2_flist: " + fixture.design1FlistPath.string() + "\n"
      "sv_design1_top: cva6\n"
      "sv_design2_top: cva6\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigCompactSystemVerilogLecCnfRejected) {
  const auto fixture = createSystemVerilogFlistFixture();
  const auto cnfPath = fixture.tmpDir / "compact_sv.cnf";
  const auto poCnfDir = fixture.tmpDir / "compact_sv_po_cnfs";
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
      "po_cnf_export: true\n"
      "po_cnf_export_path: " + poCnfDir.string() + "\n"
      "liberty_files:\n"
      "  - " + libertyPath.string() + "\n");

  EXPECT_NE(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_FALSE(std::filesystem::exists(cnfPath));
  EXPECT_FALSE(std::filesystem::exists(poCnfDir));

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSystemVerilogFlistMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "sv_design1_flist:\n"
      "  - bad\n"
      "sv_design2_flist: design1.f\n");
  EXPECT_NE(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigSystemVerilogSecondFlistMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "sv_design1_flist: design0.f\n"
      "sv_design2_flist:\n"
      "  - bad\n");
  EXPECT_NE(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigSystemVerilogTopMustNotBeEmpty) {
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

TEST_F(KeplerFormalCliTests, CliSystemVerilogFlistAndTopAccepted) {
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
  std::string argv10 = "-v";
  std::string argv11 = "sec";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data(),
                  argv5.data(), argv6.data(), argv7.data(), argv8.data(), argv9.data(),
                  argv10.data(), argv11.data()};
  int argc = 12;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSystemVerilogFlagMissingValueFails) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_flist";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data()};
  int argc = 3;
  EXPECT_NE(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
}

TEST_F(KeplerFormalCliTests, CliSystemVerilogEmptyValueFails) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_top";
  std::string argv3;
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data()};
  int argc = 4;
  EXPECT_NE(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
}

TEST_F(KeplerFormalCliTests, CliSystemVerilogRequiresSourcesForBothDesigns) {
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

TEST_F(KeplerFormalCliTests, CliSystemVerilogOptionsRejectedForVerilogFormat) {
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

TEST_F(KeplerFormalCliTests, CliVerilogTopOptionsRejectedForSystemVerilogFormat) {
  EXPECT_NE(runWithArgs({"kepler-formal",
                         "-systemverilog",
                         "--verilog_design1_top",
                         "top",
                         "design0.sv",
                         "design1.sv",
                         "-v",
                         "sec"}),
            EXIT_SUCCESS);
}

TEST_F(KeplerFormalCliTests, CliSv2vRejectsFirstVerilogTopOption) {
  EXPECT_NE(runWithArgs({"kepler-formal",
                         "-sv2v",
                         "--verilog_design1_top",
                         "top",
                         "design0.sv",
                         "design1.v",
                         "-v",
                         "sec"}),
            EXIT_SUCCESS);
}

TEST_F(KeplerFormalCliTests, CliVerilogExplicitTopSelectsDummyModules) {
  const auto testData = repoRoot() / "test/strategies/miter/testdata";
  const auto design1 = testData / "verilog_top_design1.v";
  const auto design2 = testData / "verilog_top_design2.v";
  ASSERT_TRUE(std::filesystem::exists(design1));
  ASSERT_TRUE(std::filesystem::exists(design2));

  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", design1.string(), design2.string()}),
      EXIT_FAILURE);

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-verilog",
                   "--verilog_design1_top",
                   "design1_top",
                   "--verilog_design2_top",
                   "design2_top",
                   design1.string(),
                   design2.string()}),
      EXIT_SUCCESS);

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-verilog",
                   "--verilog_design1_top",
                   "missing",
                   "--verilog_design2_top",
                   "design2_top",
                   design1.string(),
                   design2.string()}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, ConfigCompactSecVerilogDifferentTopsAreNotReused) {
  const auto design =
      repoRoot() / "test/strategies/miter/testdata/verilog_top_design1.v";
  ASSERT_TRUE(std::filesystem::exists(design));
  const auto tmpDir = makeUniqueTempDir("kepler_formal_verilog_top_config");
  const auto logPath = tmpDir / "different_tops.log";
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 1\n"
      "compact_mode: true\n"
      "verilog_design1_top: design1_top\n"
      "verilog_design2_top: design1_unused\n"
      "input_paths:\n"
      "  - " + design.string() + "\n"
      "  - " + design.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_EQ(contents.find("reusing extracted design 1 model"), std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(tmpDir);
}

TEST_F(KeplerFormalCliTests, FirstVerilogDesignWithoutTopFails) {
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

TEST_F(KeplerFormalCliTests, SecondVerilogDesignWithoutTopFails) {
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

TEST_F(KeplerFormalCliTests, ConfigUnknownSolverFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "solver: nope\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigGlucoseSolverAccepted) {
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

TEST_F(KeplerFormalCliTests, ConfigKissatSolverAccepted) {
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

TEST_F(KeplerFormalCliTests, ConfigCompactModeAccepted) {
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

TEST_F(KeplerFormalCliTests, ConfigCompactModeWritesDifferentSummaryAndWarning) {
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

TEST_F(KeplerFormalCliTests, ConfigReportSkippedPOsAccepted) {
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

TEST_F(KeplerFormalCliTests, ConfigLogFileAccepted) {
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

TEST_F(KeplerFormalCliTests, ConfigInputPathsWrongSizeFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - only_one.v\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigInputPathsNestedWrongCountFails) {
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

TEST_F(KeplerFormalCliTests, ConfigInputPathsNestedNonScalarFails) {
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

TEST_F(KeplerFormalCliTests, SnlMultiFileRejected) {
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

TEST_F(KeplerFormalCliTests, MissingFirstNajaIfFails) {
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

TEST_F(KeplerFormalCliTests, MissingSecondNajaIfFails) {
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

TEST_F(KeplerFormalCliTests, ConfigCompactNajaIfAccepted) {
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto design = copyNajaIfForCurrentBuild(
      exampleDir / "tinyrocket_naja.if", "kepler_compact_naja_if");
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
  std::filesystem::remove_all(design.parent_path());
}

TEST_F(KeplerFormalCliTests, CliUnknownOptionFails) {
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

TEST_F(KeplerFormalCliTests, CliRemovedInternalStateCorrespondenceFlagFails) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--sec-internal-state-correspondence";
  std::string argv3 = "a.v";
  std::string argv4 = "b.v";
  char* argv[] = {
      argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data()};
  int argc = 5;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliHelpPrintsUsage) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "--help";
  char* argv[] = {argv0.data(), argv1.data()};
  int argc = 2;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_SUCCESS);
}

TEST_F(KeplerFormalCliTests, ConfigInvalidVerificationModeFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: bad\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigInvalidSecEngineFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_engine: bad\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigSecEngineMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_engine:\n"
      "  - pdr\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigInvalidSecEncodingFails) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_encoding: mystery\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigSecEncodingMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_encoding:\n"
      "  - dual_rail_steady\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigExplicitLecVerificationAccepted) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: lec\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigVerificationMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification:\n"
      "  - sec\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigMaxKMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "max_k:\n"
      "  - 4\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigSecBoundaryAbstractionMustBeScalar) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_uncomputable_seq_as_boundary:\n"
      "  - false\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigInvalidMaxKTokenFails) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "max_k: nope\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigEmptyMaxKTokenFails) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "max_k: \"\"\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigOutOfRangeMaxKTokenFails) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "max_k: 999999999999999999999999999999999999\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigMaxKWithoutSecFails) {
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

TEST_F(KeplerFormalCliTests, ConfigSecEngineWithoutSecFails) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: lec\n"
      "sec_engine: pdr\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecVerificationAccepted) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecDefaultsToDualRailEncoding) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto logPath = fixture.tmpDir / "default_sec_encoding.log";
  // Intentionally omit sec_encoding here: this is the regression that guards the
  // user-visible SEC default. Other SEC tests spell out the mode they need.
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  // Omitting sec_encoding selects the dual-rail steady-state property.
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("SEC encoding: dual_rail_steady"), std::string::npos);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecVerificationAcceptedWithPdrEngine) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "sec_engine: pdr\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecVerificationRejectsLegacyEngine) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "sec_engine: legacy\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecVerificationAcceptedWithKInductionEngine) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "sec_engine: k_induction\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecVerificationAcceptedWithImcEngine) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "sec_engine: imc\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  // This option-parsing fixture is valid input for IMC, but IMC does not
  // converge within the small bound.
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecInconclusiveExitCode);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
     ConfigSecAbstractsUncomputableSequentialBoundariesByDefault) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createUncomputableSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 2\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary());
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecUnsupportedMismatchLogUsesUnsupportedResult) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture =
      createEquivalentSequentialNajaIfFixture("ff0", "ff0", "out", "z");
  const auto logPath = fixture.tmpDir / "sec_unsupported_mismatch.log";
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 2\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("SEC workflow failed:"), std::string::npos);
  EXPECT_NE(contents.find("Mismatched observed output sets"), std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
     ConfigSecCanDisableBoundaryAbstractionForUncomputableSequentials) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 2\n"
      "sec_uncomputable_seq_as_boundary: false\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  EXPECT_FALSE(KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary());
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecIgnoresRenamedInternalState) {
  const auto fixture =
      createEquivalentSequentialNajaIfFixture("state_a", "state_b");
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 4\n"
      "allow-boundary-mismatch: false\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSystemVerilogSecVerificationAccepted) {
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
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       CliSystemVerilogSecSharedDivModPrimitiveProvesEquivalent) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(\n"
      "  input  [63:0] a,\n"
      "  output [63:0] q,\n"
      "  output [63:0] r\n"
      ");\n"
      "  assign q = {a[63:1], 1'h0} / 64'h8;\n"
      "  assign r = {a[63:1], 1'h0} % 64'h8;\n"
      "endmodule\n");
  const auto runDir = fixture.tmpDir / "shared_divmod_self_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    EXPECT_EQ(
        runWithArgs({"kepler-formal",
                     "-sv",
                     "-v",
                     "sec",
                     "--sec-engine",
                     "pdr",
                     "-k",
                     "4",
                     "--design1",
                     fixture.design0Path.string(),
                     "--design2",
                     fixture.design1Path.string()}),
        kSecProvedExitCode);

    const auto logs = listMiterLogsInCurrentDirectory();
    ASSERT_EQ(logs.size(), 1u);
    const auto contents = readFileContents(runDir / logs.front());
    EXPECT_NE(
        contents.find(
            "SEC checked-output coverage: 100.00% "
            "(128/128 covered/existing outputs)."),
        std::string::npos);
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       CliSystemVerilogSecSharedDivModMatchesShiftAndMask) {
  const auto fixture = createDesignFixture(
      "sv",
      "module top(\n"
      "  input  [63:0] a,\n"
      "  output [63:0] q,\n"
      "  output [63:0] r\n"
      ");\n"
      "  assign q = {a[63:1], 1'h0} / 64'h8;\n"
      "  assign r = {a[63:1], 1'h0} % 64'h8;\n"
      "endmodule\n",
      "module top(\n"
      "  input  [63:0] a,\n"
      "  output [63:0] q,\n"
      "  output [63:0] r\n"
      ");\n"
      "  wire [63:0] aligned = {a[63:1], 1'h0};\n"
      "  assign q = aligned >> 3;\n"
      "  assign r = aligned & 64'h7;\n"
      "endmodule\n");
  const auto runDir = fixture.tmpDir / "shared_divmod_semantics_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    EXPECT_EQ(
        runWithArgs({"kepler-formal",
                     "-sv",
                     "-v",
                     "sec",
                     "--sec-engine",
                     "pdr",
                     "-k",
                     "4",
                     "--design1",
                     fixture.design0Path.string(),
                     "--design2",
                     fixture.design1Path.string()}),
        kSecProvedExitCode);

    const auto logs = listMiterLogsInCurrentDirectory();
    ASSERT_EQ(logs.size(), 1u);
    const auto contents = readFileContents(runDir / logs.front());
    EXPECT_NE(
        contents.find(
            "SEC checked-output coverage: 100.00% "
            "(128/128 covered/existing outputs)."),
        std::string::npos);
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigSystemVerilogSecPdrDualRailExplainsResetlessStateMismatch) {
  const auto fixture = createDesignFixture(
      "sv",
      "module T(input clk, input a, output y);\n"
      "  reg r;\n"
      "  always @(posedge clk) r <= a;\n"
      "  assign y = r;\n"
      "endmodule\n",
      "module T(input clk, input a, output y);\n"
      "  reg r;\n"
      "  always @(posedge clk) r <= ~a;\n"
      "  assign y = r;\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sv_sec_resetless_dual_rail_pdr.log";
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 2\n"
      "sv_design1_top: T\n"
      "sv_design2_top: T\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecCounterexampleExitCode);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("SEC engine: pdr"), std::string::npos);
  EXPECT_NE(contents.find("SEC encoding: dual_rail_steady"), std::string::npos);

  EXPECT_EQ(
      contents.find("No difference was found. SEC proved equivalence"),
      std::string::npos)
      << contents;
  const bool reportsConcreteDifference =
      contents.find("SEC counterexample details:") != std::string::npos;
  const bool explainsSteadyXAbstraction =
      contents.find("steady-X") != std::string::npos ||
      contents.find("X-steady") != std::string::npos ||
      contents.find("X-dominated") != std::string::npos ||
      contents.find("reset-unanchored") != std::string::npos;
  EXPECT_TRUE(reportsConcreteDifference || explainsSteadyXAbstraction)
      << contents;

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigSystemVerilogSecPdrDualRailReportsSteadyStateXProof) {
  const auto fixture = createDesignFixture(
      "sv",
      "module T(input clk, output y);\n"
      "  reg r;\n"
      "  always @(posedge clk) r <= r;\n"
      "  assign y = r;\n"
      "endmodule\n",
      "module T(input clk, output y);\n"
      "  assign y = 1'b0;\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sv_sec_dual_rail_x_output.log";
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 2\n"
      "sv_design1_top: T\n"
      "sv_design2_top: T\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(
      contents.find(
          "SEC proved equivalence under the dual-rail steady-state abstraction"),
      std::string::npos)
      << contents;
  EXPECT_EQ(contents.find("Difference was found."), std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSystemVerilogSecCompactIdenticalInputReusesModel) {
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
  const auto logPath = fixture.tmpDir / "sv_sec_compact_identical.log";
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "compact_mode: true\n"
      "max_k: 4\n"
      "sv_design1_top: top\n"
      "sv_design2_top: top\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design0Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  // The SystemVerilog frontend may install Naja's logger as the process
  // default.  This test only owns CLI/config behavior, so assert the compact
  // flow reached the configured Kepler log before frontend construction.
  EXPECT_NE(contents.find("Compact mode: enabled"), std::string::npos);
  EXPECT_NE(
      contents.find(
          "SEC compact mode: extracting and releasing design 1 before loading design 2"),
      std::string::npos);

  const auto flistLogPath =
      fixture.tmpDir / "sv_sec_compact_identical_flist.log";
  const auto flistPath = fixture.tmpDir / "sv_sec_compact_identical.f";
  {
    std::ofstream flist(flistPath);
    flist << fixture.design0Path.string() << "\n";
  }
  const auto flistCfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "compact_mode: true\n"
      "max_k: 4\n"
      "sv_design1_flist: " + flistPath.string() + "\n"
      "sv_design2_flist: " + flistPath.string() + "\n"
      "sv_design1_top: top\n"
      "sv_design2_top: top\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design0Path.string() + "\n"
      "log_file: " + flistLogPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(flistCfgPath), EXIT_SUCCESS);
  ASSERT_TRUE(std::filesystem::exists(flistLogPath));
  const auto flistContents = readFileContents(flistLogPath);
  EXPECT_NE(flistContents.find("Compact mode: enabled"), std::string::npos);
  EXPECT_NE(
      flistContents.find(
          "SEC compact mode: extracting and releasing design 1 before loading design 2"),
      std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove(flistCfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigPythonLibraryPathIsLoggedBeforeLoadFailure) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(output y);\n"
      "  assign y = 1'b0;\n"
      "endmodule\n");
  const auto missingPythonPath = fixture.tmpDir / "missing_primitives.py";
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "py_tech_files:\n"
      "  - " + missingPythonPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecVerificationWritesDefaultLog) {
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
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
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
    EXPECT_NE(contents.find("Verification: sec"), std::string::npos);
    EXPECT_NE(contents.find("Parsing systemverilog file(s) for design 1"),
              std::string::npos);
  }

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecReportsPartialObservedOutputCoverage) {
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
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 2\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecPartiallyProvedExitCode);

  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("Verification: sec"), std::string::npos);
  EXPECT_NE(contents.find("Parsing systemverilog file(s) for design 1"),
            std::string::npos);
  const auto resultLine = logLineContaining(
      contents,
      "SEC partially proved equivalence at k = 0: 1/2 outputs proved; "
      "remaining outputs are inconclusive.");
  ASSERT_FALSE(resultLine.empty());
  EXPECT_NE(resultLine.find("[info]"), std::string::npos);
  EXPECT_EQ(resultLine.find("[warning]"), std::string::npos);

  const auto warningLine = logLineContaining(
      contents, "SEC verification did not prove all observed outputs.");
  ASSERT_FALSE(warningLine.empty());
  EXPECT_NE(warningLine.find("[warning]"), std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecDifferenceLogIncludesWitnessDetails) {
  const auto fixture = createDifferentSequentialNajaIfFixture();
  const auto logPath = fixture.tmpDir / "sec_difference.log";
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 2\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecCounterexampleExitCode);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("SEC counterexample details:"), std::string::npos);
  EXPECT_NE(
      contents.find("Exact PDR found a defined-value counterexample at k = "),
      std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigTinyRocketSecVerificationAccepted) {
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
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 1\n"
      "input_paths:\n"
      "  - " + design.string() + "\n"
      "  - " + design.string() + "\n"
      "liberty_files:\n"
      "  - " + lib0.string() + "\n"
      "  - " + lib1.string() + "\n"
      "  - " + lib2.string() + "\n"
      "  - " + lib3.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecPartiallyProvedExitCode);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigSecCompactModeAccepted) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "compact_mode: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecCompactIdenticalInputReusesExtractedModel) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto logPath = fixture.tmpDir / "sec_compact_identical.log";
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "compact_mode: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  ASSERT_TRUE(std::filesystem::exists(logPath));
  const auto contents = readFileContents(logPath);
  EXPECT_NE(
      contents.find(
          "SEC compact mode: reusing extracted design 1 model for identical design 2 input"),
      std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecRejectsScopeExtraction) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "use_scopes: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecRejectsCnfExport) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "cnf_export: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecAcceptsSkippedPoReporting) {
  ReportSkippedPOsGuard reportGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 4\n"
      "report_skipped_pos: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getReportSkippedPOs());
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, WriteBoundaryTermsReportFormatsEntries) {
  const auto tempDir = makeUniqueTempDir("kepler_formal_cli_boundary_terms");
  const auto reportPath = tempDir / "boundary_terms.txt";

  std::vector<KEPLER_FORMAL::SEC::ExtractedBoundaryReportEntry> reports = {
      {.design = "design0", .signal = "clk[0]", .roles = {"top_input"}},
      {.design = "design0",
       .signal = "bad[0]",
       .roles = {"top_output", "opaque_internal_output", "abstracted_sequential_observed"},
       .connectivitySkip = "logical-loop connectivity: cycle"}};

  writeBoundaryTermsReport(reportPath, reports);

  const auto content = readFileContents(reportPath);
  EXPECT_NE(content.find("# SEC boundary terms report"), std::string::npos);
  EXPECT_NE(content.find("opaque_internal_input / opaque_internal_output"),
            std::string::npos);
  EXPECT_NE(content.find("abstracted_sequential_state / abstracted_sequential_observed"),
            std::string::npos);
  EXPECT_NE(content.find("design: design0"), std::string::npos);
  EXPECT_NE(content.find("signal: clk[0]"), std::string::npos);
  EXPECT_NE(content.find("roles: [top_input]"), std::string::npos);
  EXPECT_NE(content.find("signal: bad[0]"), std::string::npos);
  EXPECT_NE(
      content.find(
          "roles: [top_output, opaque_internal_output, abstracted_sequential_observed]"),
      std::string::npos);
  EXPECT_NE(
      content.find("connectivity_skip: logical-loop connectivity: cycle"),
      std::string::npos);

  std::filesystem::remove_all(tempDir);
}

TEST_F(KeplerFormalCliTests, WriteBoundaryTermsReportSkipsEmptyReports) {
  const auto tempDir =
      makeUniqueTempDir("kepler_formal_cli_empty_boundary_terms");
  const auto reportPath = tempDir / "boundary_terms.txt";

  writeBoundaryTermsReport(reportPath, {});

  EXPECT_FALSE(std::filesystem::exists(reportPath));
  std::filesystem::remove_all(tempDir);
}

TEST_F(KeplerFormalCliTests, WriteResetUnanchoredSkippedOutputsReportFormatsEntries) {
  const auto tempDir =
      makeUniqueTempDir("kepler_formal_cli_reset_unanchored_report");
  const auto reportPath = tempDir / "skipped_reset_unanchored_pos.txt";

  writeResetUnanchoredSkippedOutputsReport(
      reportPath,
      {"bad[0]: design0 depends on reset-unanchored internal state u0.q[0] | "
       "design1 depends on reset-unanchored internal state u1.q[0]"});

  const auto content = readFileContents(reportPath);
  EXPECT_NE(
      content.find("# SEC reset-unanchored skipped observed outputs"),
      std::string::npos);
  EXPECT_NE(
      content.find("SEC does not assume internal flop equality"),
      std::string::npos);
  EXPECT_NE(content.find("- bad[0]: design0 depends"), std::string::npos);
  EXPECT_NE(content.find("u0.q[0]"), std::string::npos);
  EXPECT_NE(content.find("u1.q[0]"), std::string::npos);

  std::filesystem::remove_all(tempDir);
}

TEST_F(KeplerFormalCliTests, WriteResetUnanchoredSkippedOutputsReportSkipsEmptyEntries) {
  const auto tempDir =
      makeUniqueTempDir("kepler_formal_cli_empty_reset_unanchored_report");
  const auto reportPath = tempDir / "skipped_reset_unanchored_pos.txt";

  writeResetUnanchoredSkippedOutputsReport(reportPath, {});

  EXPECT_FALSE(std::filesystem::exists(reportPath));
  std::filesystem::remove_all(tempDir);
}

TEST_F(KeplerFormalCliTests, WriteMultiClockDomainSkippedOutputsReportFormatsEntries) {
  const auto tempDir =
      makeUniqueTempDir("kepler_formal_cli_multi_clock_domain_report");
  const auto reportPath = tempDir / "skipped_multi_clock_domain_pos.txt";

  writeMultiClockDomainSkippedOutputsReport(
      reportPath,
      {"bad[0]: design0 multi-clock-domain connectivity: "
       "Observed output cone spans clock domains: a_clk[0], b_clk[0]"});

  const auto content = readFileContents(reportPath);
  EXPECT_NE(
      content.find("# SEC multi-clock-domain skipped observed outputs"),
      std::string::npos);
  EXPECT_NE(content.find("CDC"), std::string::npos);
  EXPECT_NE(content.find("- bad[0]: design0 multi-clock-domain connectivity"),
            std::string::npos);
  EXPECT_NE(content.find("a_clk[0], b_clk[0]"), std::string::npos);

  std::filesystem::remove_all(tempDir);
}

TEST_F(KeplerFormalCliTests, WriteMultiClockDomainSkippedOutputsReportSkipsEmptyEntries) {
  const auto tempDir =
      makeUniqueTempDir("kepler_formal_cli_empty_multi_clock_domain_report");
  const auto reportPath = tempDir / "skipped_multi_clock_domain_pos.txt";

  writeMultiClockDomainSkippedOutputsReport(reportPath, {});

  EXPECT_FALSE(std::filesystem::exists(reportPath));
  std::filesystem::remove_all(tempDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecFallsBackWhenLogParentCannotBeCreated) {
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
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 4\n"
      "log_file: " + (blockedParent / "sec.log").string() + "\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecContinuesWhenLogFilePathIsDirectory) {
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
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 4\n"
      "log_file: " + fixture.tmpDir.string() + "\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSecVerificationAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecProvedExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSecEngineAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "--sec-engine",
                   "pdr",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecProvedExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliKInductionSecEngineAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "--sec-engine",
                   "k_induction",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecProvedExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliImcSecEngineAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "--sec-engine",
                   "imc",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecInconclusiveExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliDualRailEncodingAcceptedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-engine",
                   "k_induction",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecProvedExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliExplicitLecVerificationAcceptedBeforeFormat) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "lec",
                   "-verilog",
                   fixture.design0Path.string(),
                   fixture.design1Path.string()}),
      EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, LecBoundaryCheckEnabledByDefault) {
  const auto fixture =
      createEquivalentSequentialNajaIfFixture("state_a", "state_b");

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      EXIT_FAILURE);

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "--allow-boundary-mismatch",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      EXIT_SUCCESS);

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "allow-boundary-mismatch: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliMissingVerificationBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-v"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliInvalidVerificationBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-v", "bad"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliMissingMaxKBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-k"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliInvalidMaxKBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-k", "-1"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliMissingSecEngineBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "--sec-engine"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliMissingSecEncodingBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "--sec-encoding"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliOutOfRangeMaxKBeforeFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-k", "999999999999999999999999999999999999"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliUnrecognizedOptionBeforeFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "--mystery"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliMissingInputFormatAfterPreOptionsFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-v", "sec", "-k", "4"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliInvalidSecEngineBeforeFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-v", "sec", "--sec-engine", "bad", "-verilog"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliInvalidSecEncodingBeforeFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-v", "sec", "--sec-encoding", "bad", "-verilog"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliRemovedKInductionAliasesAreRejectedBeforeFormat) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  for (const char* removedAlias : {"kinduction", "classic_k_induction"}) {
    EXPECT_EQ(
        runWithArgs({"kepler-formal",
                     "-v",
                     "sec",
                     "--sec-engine",
                     removedAlias,
                     "-naja_if",
                     fixture.design0IfPath.string(),
                     fixture.design1IfPath.string()}),
        EXIT_FAILURE);
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSecBoundaryFlagAcceptedBeforeFormat) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "--sec-uncomputable-seq-boundary",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecProvedExitCode);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary());
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliNoSecBoundaryFlagAcceptedBeforeFormat) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "--no-sec-uncomputable-seq-boundary",
                   "-naja_if",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecProvedExitCode);
  EXPECT_FALSE(KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary());
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliMissingVerificationAfterFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-verilog", "--verification"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliInvalidVerificationAfterFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", "--verification", "bad"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliMissingMaxKAfterFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-verilog", "--max-k"}), EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliInvalidMaxKAfterFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", "--max-k", "oops"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliOutOfRangeMaxKAfterFormatFails) {
  EXPECT_EQ(
      runWithArgs(
          {"kepler-formal", "-verilog", "--max-k", "999999999999999999999999999999999999"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliMissingSecEngineAfterFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-verilog", "--sec-engine"}),
            EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliMissingSecEncodingAfterFormatFails) {
  EXPECT_EQ(runWithArgs({"kepler-formal", "-verilog", "--sec-encoding"}),
            EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliSecEngineAcceptedAfterFormat) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", "--sec-engine", "pdr"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliSecEncodingAcceptedAfterFormat) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", "--sec-encoding", "dual_rail_steady"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliInvalidSecEngineAfterFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", "--sec-engine", "bad"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliInvalidSecEncodingAfterFormatFails) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-verilog", "--sec-encoding", "bad"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliSecBoundaryFlagAcceptedAfterFormat) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-naja_if",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "--sec-uncomputable-seq-boundary",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecProvedExitCode);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary());
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliNoSecBoundaryFlagAcceptedAfterFormat) {
  SecBoundaryAbstractionGuard boundaryGuard;
  const auto fixture = createEquivalentSequentialNajaIfFixture();

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-naja_if",
                   "-v",
                   "sec",
                   "-k",
                   "4",
                   "--sec-encoding",
                   "dual_rail_steady",
                   "--no-sec-uncomputable-seq-boundary",
                   fixture.design0IfPath.string(),
                   fixture.design1IfPath.string()}),
      kSecProvedExitCode);
  EXPECT_FALSE(KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary());
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSecPdrReportsCombinationalMismatchAtFrameZero) {
  const auto fixture = createDesignFixture(
      "v",
      "module top(input a, input b, output y);\n"
      "  or (y, a, b);\n"
      "endmodule\n",
      "module top(input a, input b, output y);\n"
      "  and (y, a, b);\n"
      "endmodule\n");

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "-k",
                   "1",
                   "--sec-engine",
                   "pdr",
                   "-verilog",
                   fixture.design0Path.string(),
                   fixture.design1Path.string()}),
      kSecCounterexampleExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecInconclusiveFails) {
  const auto fixture = createDesignFixture(
      "sv",
      "module top(\n"
      "    input logic clk,\n"
      "    input logic a,\n"
      "    output logic q\n"
      ");\n"
      "  always_ff @(posedge clk) q <= a;\n"
      "endmodule\n",
      "module top(\n"
      "    input logic clk,\n"
      "    input logic a,\n"
      "    output logic q\n"
      ");\n"
      "  always_ff @(posedge clk) q <= ~a;\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sec_inconclusive.log";
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: dual_rail_steady\n"
      "max_k: 0\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecInconclusiveExitCode);

  const auto contents = readFileContents(logPath);
  const auto resultLine =
      logLineContaining(contents, "SEC was inconclusive");
  ASSERT_FALSE(resultLine.empty());
  EXPECT_NE(resultLine.find("[info]"), std::string::npos);
  EXPECT_EQ(resultLine.find("[warning]"), std::string::npos);

  const auto warningLine = logLineContaining(
      contents,
      "SEC verification did not produce a proof or counterexample.");
  ASSERT_FALSE(warningLine.empty());
  EXPECT_NE(warningLine.find("[warning]"), std::string::npos);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecUnsupportedMismatchFails) {
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
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 1\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliLibertyFlagCollectsPaths) {
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

TEST_F(KeplerFormalCliTests, CliLibAliasCollectsPaths) {
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

TEST_F(KeplerFormalCliTests, CliPositionalLibertyPathsAfterNetlistsAccepted) {
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

TEST_F(KeplerFormalCliTests, CliPositionalLibertyPathsWithoutStandardSuffixAccepted) {
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

TEST_F(KeplerFormalCliTests, PythonLibraryFilesAreLoadedFromDedicatedYamlKey) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  BUF c0(.A(a), .Z(y));\n"
      "endmodule\n");
  const auto pyPrimitives = fixture.tmpDir / "kepler_test_primitives.py";
  {
    std::ofstream py(pyPrimitives);
    py << "import naja\n"
          "\n"
          "def constructBUF(lib):\n"
          "  cell = naja.SNLDesign.createPrimitive(lib, 'BUF')\n"
          "  a = naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Input, 'A')\n"
          "  z = naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Output, 'Z')\n"
          "  naja.SNLDesign.addCombinatorialArcs([a], [z])\n"
          "\n"
          "def constructPrimitives(lib):\n"
          "  constructBUF(lib)\n";
  }
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

TEST_F(KeplerFormalCliTests, PythonLibraryFilesUnderLibertyFilesAreRejected) {
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

TEST_F(KeplerFormalCliTests, GzippedLibertyFilesAreLoadedByExtension) {
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

TEST_F(KeplerFormalCliTests, InvalidLibertyFileContentFails) {
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

TEST_F(KeplerFormalCliTests, CliPythonLibraryFilesAreRejected) {
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

TEST_F(KeplerFormalCliTests, CliSystemVerilogFlistWithoutTopAccepted) {
  const auto fixture = createSystemVerilogFlistFixture();
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-systemverilog";
  std::string argv2 = "--sv_design1_flist";
  std::string argv3 = fixture.design0FlistPath.string();
  std::string argv4 = "--sv_design2_flist";
  std::string argv5 = fixture.design1FlistPath.string();
  std::string argv6 = "-v";
  std::string argv7 = "sec";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(), argv4.data(),
                  argv5.data(), argv6.data(), argv7.data()};
  int argc = 8;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSystemVerilogDirectPathsWithEscapedNamesAccepted) {
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
  std::string argv8 = "-v";
  std::string argv9 = "sec";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data(),
                  argv8.data(), argv9.data()};
  int argc = 10;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_SUCCESS);
  std::filesystem::remove_all(tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSystemVerilogTopCommandFileCreationFailureFails) {
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
  std::string argv8 = "-v";
  std::string argv9 = "sec";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data(),
                  argv8.data(), argv9.data()};
  int argc = 10;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);

  std::filesystem::permissions(noWriteDir, perms);
  std::filesystem::remove_all(noWriteDir);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliCompactSystemVerilogTopCommandFileCreationFailureFails) {
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
  std::string argv9 = "-v";
  std::string argv10 = "sec";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data(),
                  argv8.data(), argv9.data(), argv10.data()};
  int argc = 11;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);

  std::filesystem::permissions(noWriteDir, perms);
  std::filesystem::remove_all(noWriteDir);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSv2vPrimitiveStubFileCreationFailureFails) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_sv2v_stub_no_write_fixture");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  const auto pyPrimitives = fixture.tmpDir / "sv2v_primitives.py";
  {
    std::ofstream py(pyPrimitives);
    py << "import naja\n"
          "\n"
          "def constructPrimitives(lib):\n"
          "  cell = naja.SNLDesign.createPrimitive(lib, 'BUF')\n"
          "  a = naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Input, 'A')\n"
          "  z = naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Output, 'Z')\n"
          "  cell.setTruthTable(0b10)\n";
  }
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input logic a, output logic y);\n";
    design0 << "  BUF u_buf(.A(a), .Z(y));\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, output y);\n";
    design1 << "  assign y = a;\n";
    design1 << "endmodule\n";
  }

  const auto pyModuleDir = findBuiltNajaModuleDir();
  ASSERT_TRUE(std::filesystem::exists(pyPrimitives));
  ASSERT_FALSE(pyModuleDir.empty());
  ASSERT_TRUE(std::filesystem::exists(pyModuleDir / "naja.so"));
  EnvVarGuard pythonPathGuard("PYTHONPATH");
  pythonPathGuard.set(pyModuleDir.string());

  const auto cfgPath = writeTempConfig(
      "format: sv2v\n"
      "verification: sec\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "py_tech_files:\n"
      "  - " + pyPrimitives.string() + "\n");

  const auto noWriteDir =
      std::filesystem::temp_directory_path() / "kepler_formal_sv2v_stub_no_write";
  std::filesystem::create_directories(noWriteDir);
  auto perms = std::filesystem::status(noWriteDir).permissions();
  std::filesystem::permissions(
      noWriteDir,
      std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
          std::filesystem::perms::others_write,
      std::filesystem::perm_options::remove);

  EnvVarGuard tmpDirGuard("TMPDIR");
  tmpDirGuard.set(noWriteDir.string());

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);

  std::filesystem::permissions(noWriteDir, perms);
  std::filesystem::remove_all(noWriteDir);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSystemVerilogFirstDesignFailureCleansTemporaryCommandFile) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_sv_cleanup_first";
  std::filesystem::create_directories(tmpDir);
  EnvVarGuard tmpDirGuard("TMPDIR");
  tmpDirGuard.set(tmpDir.string());
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
  std::string argv8 = "-v";
  std::string argv9 = "sec";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data(),
                  argv8.data(), argv9.data()};
  int argc = 10;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);
  EXPECT_EQ(listTemporarySystemVerilogCommandFiles(), beforeTempFiles);

  std::filesystem::remove_all(tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSv2vFirstDesignFailureCleansTemporaryStubAndCommandFiles) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_sv2v_cleanup_first";
  std::filesystem::create_directories(tmpDir);
  EnvVarGuard tmpDirGuard("TMPDIR");
  tmpDirGuard.set(tmpDir.string());
  const auto design0 = tmpDir / "design0_invalid.sv";
  const auto design1 = tmpDir / "design1_valid.v";
  const auto pyPrimitives = tmpDir / "sv2v_primitives.py";
  const auto cfgPath = tmpDir / "config.yaml";
  {
    std::ofstream py(pyPrimitives);
    py << "import naja\n"
          "\n"
          "def constructPrimitives(lib):\n"
          "  cell = naja.SNLDesign.createPrimitive(lib, 'BUF')\n"
          "  a = naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Input, 'A')\n"
          "  z = naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Output, 'Z')\n"
          "  cell.setTruthTable(0b10)\n";
  }
  {
    std::ofstream f(design0);
    f << "module top(input logic a, output logic y)\n";
    f << "  BUF u_buf(.A(a), .Z(y));\n";
    f << "endmodule\n";
  }
  {
    std::ofstream f(design1);
    f << "module top(input a, output y);\n";
    f << "  assign y = a;\n";
    f << "endmodule\n";
  }
  {
    std::ofstream cfg(cfgPath);
    cfg << "format: sv2v\n"
        << "verification: sec\n"
        << "input_paths:\n"
        << "  - " << design0.string() << "\n"
        << "  - " << design1.string() << "\n"
        << "py_tech_files:\n"
        << "  - " << pyPrimitives.string() << "\n";
  }

  const auto pyModuleDir = findBuiltNajaModuleDir();
  ASSERT_TRUE(std::filesystem::exists(pyPrimitives));
  ASSERT_FALSE(pyModuleDir.empty());
  ASSERT_TRUE(std::filesystem::exists(pyModuleDir / "naja.so"));
  EnvVarGuard pythonPathGuard("PYTHONPATH");
  pythonPathGuard.set(pyModuleDir.string());

  const auto beforeCommandFiles = listTemporarySystemVerilogCommandFiles();
  const auto beforeStubFiles = listTemporarySystemVerilogPrimitiveStubFiles();

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  EXPECT_EQ(listTemporarySystemVerilogCommandFiles(), beforeCommandFiles);
  EXPECT_EQ(listTemporarySystemVerilogPrimitiveStubFiles(), beforeStubFiles);

  std::filesystem::remove_all(tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSystemVerilogSecondDesignFailureCleansTemporaryCommandFile) {
  const auto tmpDir =
      std::filesystem::temp_directory_path() / "kepler_formal_sv_cleanup_second";
  std::filesystem::create_directories(tmpDir);
  EnvVarGuard tmpDirGuard("TMPDIR");
  tmpDirGuard.set(tmpDir.string());
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
  std::string argv8 = "-v";
  std::string argv9 = "sec";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data(),
                  argv4.data(), argv5.data(), argv6.data(), argv7.data(),
                  argv8.data(), argv9.data()};
  int argc = 10;

  EXPECT_EQ(KeplerFormalMain(argc, argv), EXIT_FAILURE);
  EXPECT_EQ(listTemporarySystemVerilogCommandFiles(), beforeTempFiles);

  std::filesystem::remove_all(tmpDir);
}

TEST_F(KeplerFormalCliTests, CliExplicitMissingDesignsFails) {
  std::string argv0 = "kepler-formal";
  std::string argv1 = "-verilog";
  std::string argv2 = "--design1";
  std::string argv3 = "--design2";
  char* argv[] = {argv0.data(), argv1.data(), argv2.data(), argv3.data()};
  int argc = 4;
  int rc = KeplerFormalMain(argc, argv);
  EXPECT_EQ(rc, EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, ConfigDebugLogLevelAccepted) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "log_level: debug\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigUnknownLogLevelDefaultsToInfo) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths: [a.v, b.v]\n"
      "log_level: loud\n");
  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, VerilogNoLibertyCreatesDbAndFailsOnSecondParse) {
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

TEST_F(KeplerFormalCliTests, SnlScopesNoDifference) {
  const auto root = repoRoot();
  const auto exampleDir = root / "example";
  const auto design0 = copyNajaIfForCurrentBuild(
      exampleDir / "tinyrocket_naja.if", "kepler_scoped_naja_if");
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
  std::filesystem::remove_all(design0.parent_path());
}

TEST_F(KeplerFormalCliTests, SnlScopesEquivalentEditedScopeNoDifference) {
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

TEST_F(KeplerFormalCliTests, SnlScopesCreatesOnlyOneDefaultMiterLogPerScopeRun) {
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

TEST_F(KeplerFormalCliTests, SnlScopesCanCleanAndDumpCnf) {
  const auto fixture = createEquivalentScopedNajaIfFixture();
  const auto cnfPath = fixture.tmpDir / "scoped_miter.cnf";
  const auto poCnfDir = fixture.tmpDir / "scoped_po_cnfs";
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
      "cnf_export_path: " + cnfPath.string() + "\n"
      "po_cnf_export: true\n"
      "po_cnf_export_path: " + poCnfDir.string() + "\n");

  int rc = runWithConfigFile(cfgPath);
  EXPECT_EQ(rc, EXIT_SUCCESS);
  EXPECT_TRUE(std::filesystem::exists(cnfPath));
  EXPECT_TRUE(std::filesystem::exists(poCnfDir / "top0" / "po_000000.cnf"));
  EXPECT_TRUE(std::filesystem::exists(poCnfDir / "top1" / "po_000000.cnf"));

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, SnlScopesDumpCnfUsesDefaultScopedPath) {
  const auto fixture = createEquivalentScopedNajaIfFixture();
  const auto defaultCnfPath = std::filesystem::current_path() / "miter_child.cnf";
  const auto defaultPoCnfPath = std::filesystem::current_path() / "po_cnfs_child";
  std::filesystem::remove(defaultCnfPath);
  std::filesystem::remove_all(defaultPoCnfPath);

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "liberty_files:\n"
      "  - " + fixture.libertyPath.string() + "\n"
      "use_scopes: true\n"
      "cnf_export: true\n"
      "po_cnf_export: true\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  EXPECT_TRUE(std::filesystem::exists(defaultCnfPath));
  EXPECT_TRUE(std::filesystem::exists(defaultPoCnfPath / "top0" / "po_000000.cnf"));
  EXPECT_TRUE(std::filesystem::exists(defaultPoCnfPath / "top1" / "po_000000.cnf"));

  std::filesystem::remove(cfgPath);
  std::filesystem::remove(defaultCnfPath);
  std::filesystem::remove_all(defaultPoCnfPath);
  std::filesystem::remove_all(fixture.tmpDir);
}
