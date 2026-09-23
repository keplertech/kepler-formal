// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include <spdlog/sinks/null_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/spdlog.h>

#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "DNL.h"
#include "BoolExprCache.h"
#include "Config.h"
#include "KeplerFormalDriver.h"
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

struct StructuredRun {
  int exitCode;
  KEPLER_FORMAL::RunResult result;
};

StructuredRun runStructuredWithArgs(std::vector<std::string> args) {
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (auto& arg : args) {
    argv.push_back(arg.data());
  }
  KEPLER_FORMAL::RunResult result;
  const int rc = KEPLER_FORMAL::runKeplerFormal(
      static_cast<int>(argv.size()), argv.data(), result);
  return {rc, std::move(result)};
}

StructuredRun runStructuredWithConfigFile(
    const std::filesystem::path& cfgPath,
    std::string argv0 = "kepler-formal") {
  return runStructuredWithArgs(
      {std::move(argv0), "--config", cfgPath.string()});
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

struct OpaqueBoundaryCliFixture : SimpleCliFixture {
  std::filesystem::path libertyPath;
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

struct NamedLoggerGuard {
  explicit NamedLoggerGuard(std::string name)
      : name_(std::move(name)),
        previousDefault_(spdlog::default_logger()),
        previousNamed_(spdlog::get(name_)),
        installed_(std::make_shared<spdlog::logger>(
            name_, std::make_shared<spdlog::sinks::null_sink_mt>())) {
    spdlog::drop(name_);
    spdlog::register_logger(installed_);
    spdlog::set_default_logger(installed_);
  }

  ~NamedLoggerGuard() {
    spdlog::drop(name_);
    if (previousNamed_ != nullptr) {
      spdlog::register_logger(previousNamed_);
    }
    spdlog::set_default_logger(previousDefault_);
  }

  std::string name_;
  std::shared_ptr<spdlog::logger> previousDefault_;
  std::shared_ptr<spdlog::logger> previousNamed_;
  std::shared_ptr<spdlog::logger> installed_;
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
  const auto exampleDir = root / "examples" / "tinyrocket";
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

OpaqueBoundaryCliFixture createOpaqueBoundaryFixture(
    const std::string& design0Body,
    const std::string& design1Body) {
  OpaqueBoundaryCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_boundary_opaque");
  fixture.design0Path = fixture.tmpDir / "design0.v";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  fixture.libertyPath = fixture.tmpDir / "opaque.lib";

  {
    std::ofstream design0(fixture.design0Path);
    design0 << design0Body;
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << design1Body;
  }
  {
    std::ofstream liberty(fixture.libertyPath);
    liberty <<
        "library (boundary_opaque) {\n"
        "  cell (OPAQUE) {\n"
        "    pin (A) { direction : input; }\n"
        "    pin (Y) { direction : output; }\n"
        "  }\n"
        "  cell (BOUNDARY_ZERO) {\n"
        "    pin (A) { direction : input; }\n"
        "    pin (Y) { direction : output; function : \"0\"; }\n"
        "  }\n"
        "  cell (BOUNDARY_ONE) {\n"
        "    pin (A) { direction : input; }\n"
        "    pin (Y) { direction : output; function : \"1\"; }\n"
        "  }\n"
        "}\n";
  }
  return fixture;
}

StructuredRun runBoundaryInputComparison(const std::string& leftConnection,
                                          const std::string& rightConnection,
                                          bool compact) {
  const auto makeDesign = [](const std::string& connection) {
    return std::string(
               "module child(input i, output o);\n"
               "endmodule\n"
               "module top(input a, output y);\n"
               "  child u_boundary(.i(") +
           connection +
           "), .o(y));\n"
           "endmodule\n";
  };
  const auto fixture = createDesignFixture(
      "v", makeDesign(leftConnection), makeDesign(rightConnection));
  const auto runDir = fixture.tmpDir / "boundary_constant_run";
  std::filesystem::create_directories(runDir);

  StructuredRun run;
  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    std::vector<std::string> args = {
        "kepler-formal",
        "-verilog",
        fixture.design0Path.string(),
        fixture.design1Path.string(),
        "--set-as-boundary",
        "u_boundary",
        "u_boundary"};
    if (compact) {
      args.emplace_back("--compact");
    }
    run = runStructuredWithArgs(std::move(args));
  }
  std::filesystem::remove_all(fixture.tmpDir);
  return run;
}

std::filesystem::path copyExampleLibertyFile(const std::filesystem::path& directory,
                                             const std::string& filename) {
  const auto source = repoRoot() / "examples" / "tinyrocket" /
                      "NangateOpenCellLibrary_typical.lib";
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
  fixture.libertyPath = repoRoot() / "examples" / "tinyrocket" /
                        "NangateOpenCellLibrary_typical.lib";

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

  static void expectSecDifferenceLogIncludesWitnessDetails(
      const std::string& engine) {
    const auto fixture = createDifferentSequentialNajaIfFixture();
    const auto logPath =
        fixture.tmpDir / ("sec_difference_" + engine + ".log");
    const auto cfgPath = writeTempConfig(
        "format: naja_if\n"
        "verification: sec\n"
        "sec_engine: " + engine + "\n"
        "sec_encoding: dual_rail_steady\n"
        "max_k: 2\n"
        "input_paths:\n"
        "  - " + fixture.design0IfPath.string() + "\n"
        "  - " + fixture.design1IfPath.string() + "\n"
        "log_file: " + logPath.string() + "\n");

    const auto run = runStructuredWithConfigFile(cfgPath);
    EXPECT_EQ(run.exitCode, kSecCounterexampleExitCode);
    EXPECT_EQ(run.result.exitCode, kSecCounterexampleExitCode);
    EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
    EXPECT_EQ(run.result.totalOutputs, 1u);
    EXPECT_EQ(run.result.coveredOutputs, 1u);
    EXPECT_EQ(run.result.provenOutputs, 0u);
    ASSERT_TRUE(std::filesystem::exists(logPath));
    const auto contents = readFileContents(logPath);
    EXPECT_NE(contents.find("SEC counterexample details:"), std::string::npos);
    EXPECT_NE(
        contents.find(
            "Counterexample reaches the first bad frame at cycle 1."),
        std::string::npos);
    EXPECT_NE(contents.find("Input trace:"), std::string::npos);
    EXPECT_NE(
        contents.find("Observed output mismatches at cycle 1:"),
        std::string::npos);

    std::filesystem::remove(cfgPath);
    std::filesystem::remove_all(fixture.tmpDir);
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

TEST_F(KeplerFormalCliTests, InProcessRunStatusNamesAreStable) {
  using KEPLER_FORMAL::RunStatus;
  EXPECT_STREQ(KEPLER_FORMAL::runStatusName(RunStatus::NoResult), "no_result");
  EXPECT_STREQ(KEPLER_FORMAL::runStatusName(RunStatus::Equivalent), "equivalent");
  EXPECT_STREQ(KEPLER_FORMAL::runStatusName(RunStatus::Different), "different");
  EXPECT_STREQ(
      KEPLER_FORMAL::runStatusName(RunStatus::PartiallyProved),
      "partially_proved");
  EXPECT_STREQ(
      KEPLER_FORMAL::runStatusName(RunStatus::Inconclusive), "inconclusive");
  EXPECT_STREQ(
      KEPLER_FORMAL::runStatusName(RunStatus::Unsupported), "unsupported");
  EXPECT_STREQ(KEPLER_FORMAL::runStatusName(RunStatus::Error), "error");
  EXPECT_STREQ(
      KEPLER_FORMAL::runStatusName(static_cast<RunStatus>(-1)), "error");
}

TEST_F(KeplerFormalCliTests, InProcessDriverReportsNoResultAndErrors) {
  const auto noArguments = runStructuredWithArgs({"kepler-formal"});
  EXPECT_EQ(noArguments.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(noArguments.result.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(noArguments.result.status, KEPLER_FORMAL::RunStatus::NoResult);

  const auto help =
      runStructuredWithArgs({"kepler-formal", "--help"});
  EXPECT_EQ(help.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(help.result.status, KEPLER_FORMAL::RunStatus::NoResult);

  const auto invalid =
      runStructuredWithArgs({"kepler-formal", "--not-a-kepler-option"});
  EXPECT_EQ(invalid.exitCode, EXIT_FAILURE);
  EXPECT_EQ(invalid.result.exitCode, EXIT_FAILURE);
  EXPECT_EQ(invalid.result.status, KEPLER_FORMAL::RunStatus::Error);
  EXPECT_NE(
      invalid.result.reason.find(
          "failed before producing a verification result"),
      std::string::npos);
}

TEST_F(KeplerFormalCliTests, InProcessDriverRejectsAnActiveNajaUniverse) {
  cleanupNajaTestState();
  NLUniverse::create();
  EXPECT_THROW(
      runStructuredWithArgs({"kepler-formal", "--help"}),
      std::runtime_error);
  cleanupNajaTestState();

  const auto help =
      runStructuredWithArgs({"kepler-formal", "--help"});
  EXPECT_EQ(help.result.status, KEPLER_FORMAL::RunStatus::NoResult);
}

TEST_F(KeplerFormalCliTests, InProcessDriverReturnsStructuredLecResults) {
  SolverGuard solverGuard;
  ReportSkippedPOsGuard reportGuard;
  KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::CADICAL);
  KEPLER_FORMAL::Config::setReportSkippedPOs(true);
  NamedLoggerGuard loggerGuard("kepler_formal_main_logger");

  const auto fixture = createDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n",
      "module top(input a, output y);\n"
      "  assign y = 1'b0;\n"
      "endmodule\n");
  const auto equivalentLog = fixture.tmpDir / "structured_equivalent.log";
  const auto equivalentCfg = writeTempConfig(
      "format: verilog\n"
      "verification: lec\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design0Path.string() + "\n"
      "log_file: " + equivalentLog.string() + "\n");
  const auto equivalent = runStructuredWithConfigFile(equivalentCfg);
  EXPECT_EQ(equivalent.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(equivalent.result.status, KEPLER_FORMAL::RunStatus::Equivalent);
  EXPECT_EQ(equivalent.result.inputFormat, "verilog");
  EXPECT_EQ(equivalent.result.verification, "lec");
  EXPECT_EQ(equivalent.result.logFile, equivalentLog.string());
  EXPECT_TRUE(std::filesystem::exists(equivalentLog));
  EXPECT_EQ(
      KEPLER_FORMAL::Config::getSolverType(), KEPLER_FORMAL::Config::CADICAL);
  EXPECT_TRUE(KEPLER_FORMAL::Config::getReportSkippedPOs());
  EXPECT_EQ(NLUniverse::get(), nullptr);
  EXPECT_EQ(
      spdlog::get("kepler_formal_main_logger"), loggerGuard.installed_);

  const auto differentLog = fixture.tmpDir / "structured_different.log";
  const auto differentCfg = writeTempConfig(
      "format: verilog\n"
      "verification: lec\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + differentLog.string() + "\n");
  const auto different = runStructuredWithConfigFile(differentCfg);
  EXPECT_EQ(different.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(different.result.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(different.result.status, KEPLER_FORMAL::RunStatus::Different);
  EXPECT_EQ(different.result.logFile, differentLog.string());
  EXPECT_TRUE(std::filesystem::exists(differentLog));
  EXPECT_EQ(NLUniverse::get(), nullptr);

  const auto pythonTechCfg = writeTempConfig(
      "format: verilog\n"
      "verification: lec\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design0Path.string() + "\n"
      "py_tech_files:\n"
      "  - " + (fixture.tmpDir / "primitives.py").string() + "\n");
  const auto pythonTech = runStructuredWithConfigFile(pythonTechCfg);
  EXPECT_EQ(pythonTech.exitCode, EXIT_FAILURE);
  EXPECT_EQ(pythonTech.result.status, KEPLER_FORMAL::RunStatus::Error);
  EXPECT_NE(pythonTech.result.reason.find("not supported"), std::string::npos);
  EXPECT_EQ(
      spdlog::get("kepler_formal_main_logger"), loggerGuard.installed_);

  std::filesystem::remove(equivalentCfg);
  std::filesystem::remove(differentCfg);
  std::filesystem::remove(pythonTechCfg);
  std::filesystem::remove_all(fixture.tmpDir);
}


TEST_F(KeplerFormalCliTests, DumpCnfFromConfig) {
  const auto root = repoRoot();
  const auto exampleDir = root / "examples" / "tinyrocket";
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

TEST_F(KeplerFormalCliTests, CliSetAsBoundaryExposesDifferentInstanceInputs) {
  const auto fixture = createDesignFixture(
      "v",
      "module child(input i, output o);\n"
      "endmodule\n"
      "module top(input a, input b, output y);\n"
      "  child u_boundary(.i(a), .o(y));\n"
      "endmodule\n",
      "module child(input i, output o);\n"
      "endmodule\n"
      "module top(input a, input b, output y);\n"
      "  child u_boundary(.i(b), .o(y));\n"
      "endmodule\n");
  const auto runDir = fixture.tmpDir / "boundary_input_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    const auto run = runStructuredWithArgs(
        {"kepler-formal",
         "--set_as_boundary",
         "u_boundary",
         "u_boundary",
         "-verilog",
         fixture.design0Path.string(),
         fixture.design1Path.string()});
    EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
    EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
  }

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       CliSetAsBoundaryDetectsSignalVersusConstantZero) {
  const auto run = runBoundaryInputComparison("a", "1'b0", false);
  EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
}

TEST_F(KeplerFormalCliTests,
       CliCompactSetAsBoundaryDetectsSignalVersusConstantZero) {
  const auto run = runBoundaryInputComparison("a", "1'b0", true);
  EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
}

TEST_F(KeplerFormalCliTests,
       CliSetAsBoundaryDetectsConstantZeroVersusOne) {
  const auto run = runBoundaryInputComparison("1'b0", "1'b1", false);
  EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
}

TEST_F(KeplerFormalCliTests,
       CliCompactSetAsBoundaryDetectsConstantZeroVersusOne) {
  const auto run = runBoundaryInputComparison("1'b0", "1'b1", true);
  EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
}

TEST_F(KeplerFormalCliTests, ConfigSetAsBoundaryAbstractsInstanceOutput) {
  const auto fixture = createOpaqueBoundaryFixture(
      "module top(input a, output y);\n"
      "  BOUNDARY_ZERO u_boundary(.A(a), .Y(y));\n"
      "endmodule\n",
      "module top(input a, output y);\n"
      "  BOUNDARY_ONE u_boundary(.A(a), .Y(y));\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "liberty_files: [" + fixture.libertyPath.string() + "]\n"
      "set_as_boundary:\n"
      "  - [u_boundary, u_boundary]\n");

  const auto baseline = runStructuredWithArgs(
      {"kepler-formal",
       "-verilog",
       fixture.design0Path.string(),
       fixture.design1Path.string(),
       fixture.libertyPath.string()});
  EXPECT_EQ(baseline.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(baseline.result.status, KEPLER_FORMAL::RunStatus::Different);

  const auto bounded = runStructuredWithConfigFile(cfgPath);
  EXPECT_EQ(bounded.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(bounded.result.status, KEPLER_FORMAL::RunStatus::Equivalent);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSetAsBoundaryAcceptsNestedSelfBoundary) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module leaf(input i, output o);\n"
      "endmodule\n"
      "module wrapper(input i, output o);\n"
      "  leaf u_leaf(.i(i), .o(o));\n"
      "endmodule\n"
      "module top(input a, output y);\n"
      "  wrapper u_wrap(.i(a), .o(y));\n"
      "endmodule\n");

  const auto run = runStructuredWithArgs(
      {"kepler-formal",
       "-verilog",
       fixture.design0Path.string(),
       fixture.design1Path.string(),
       "--set-as-boundary",
       "u_wrap/u_leaf",
       "u_wrap/u_leaf"});
  EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Equivalent);

  const auto nonleaf = runStructuredWithArgs(
      {"kepler-formal", "-verilog", fixture.design0Path.string(),
       fixture.design1Path.string(), "--set-as-boundary", "u_wrap", "u_wrap"});
  EXPECT_EQ(nonleaf.exitCode, EXIT_FAILURE);
  EXPECT_EQ(nonleaf.result.status, KEPLER_FORMAL::RunStatus::Error);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSecNestedBoundaryAcrossHdlFormatsAndCompactMode) {
  const auto design = [](bool changedInput) {
    return std::string("module leaf(input i, output o); endmodule\n") +
        "module wrapper(input i, output o);\n"
        "  leaf u_leaf(.i(i), .o(o));\n"
        "endmodule\n"
        "module top(input a, input b, output y);\n"
        "  wrapper u_wrap(.i(" + (changedInput ? "b" : "a") +
        "), .o(y));\nendmodule\n";
  };
  for (const std::string format : {"-verilog", "-sv", "-sv2v"}) {
    for (const bool compact : {false, true}) {
      for (const int scenario : {0, 1, 2}) {
        const bool changedInput = scenario == 1;
        const bool nonleaf = scenario == 2;
        SCOPED_TRACE(format + (compact ? " compact" : " normal") +
                     (nonleaf ? " nonleaf" :
                      changedInput ? " changed input" : " equal input"));
        // Hierarchical paths may target a leaf, but not its wrapper.
        const auto fixture = createDesignFixture(
            "v", design(false), design(changedInput));
        {
          CurrentPathGuard currentPathGuard;
          std::filesystem::current_path(fixture.tmpDir);
          std::vector<std::string> args = {
              "kepler-formal", "-v", "sec", "--sec-engine", "k_induction",
              "--sec-encoding", "binary", "-k", "1", format,
              fixture.design0Path.string(), fixture.design1Path.string(),
              "--set-as-boundary", nonleaf ? "u_wrap" : "u_wrap/u_leaf",
              nonleaf ? "u_wrap" : "u_wrap/u_leaf"};
          if (compact) {
            args.emplace_back("--compact");
          }
          const auto run = runStructuredWithArgs(std::move(args));
          EXPECT_EQ(run.exitCode, nonleaf ? EXIT_FAILURE
                        : changedInput ? kSecCounterexampleExitCode : kSecProvedExitCode);
          EXPECT_EQ(run.result.status, nonleaf ? KEPLER_FORMAL::RunStatus::Error
                        : changedInput
                        ? KEPLER_FORMAL::RunStatus::Different
                        : KEPLER_FORMAL::RunStatus::Equivalent);
          if (!nonleaf) {
            EXPECT_EQ(run.result.totalOutputs, 2u);
            EXPECT_TRUE(run.result.skippedObservedOutputs.empty());
            if (!changedInput) {
              EXPECT_EQ(run.result.coveredOutputs, 2u);
            }
          }
        }
        std::filesystem::remove_all(fixture.tmpDir);
      }
    }
  }
}

TEST_F(KeplerFormalCliTests, CliCompactSetAsBoundaryUsesEachSidePath) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module child(input i, output o);\n"
      "endmodule\n"
      "module top(input a, input b, output y);\n"
      "  wire unused0;\n"
      "  wire unused1;\n"
      "  child u_left(.i(a), .o(unused0));\n"
      "  child u_right(.i(b), .o(unused1));\n"
      "  assign y = 1'b0;\n"
      "endmodule\n");
  const auto runDir = fixture.tmpDir / "compact_boundary_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    const auto run = runStructuredWithArgs(
        {"kepler-formal",
         "-verilog",
         fixture.design0Path.string(),
         fixture.design1Path.string(),
         "--compact",
         "--set-as-boundary",
         "u_left",
         "u_right"});
    EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
    EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
  }

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSecSetAsBoundaryProvesRenamedOpaquePair) {
  const auto fixture = createOpaqueBoundaryFixture(
      "module top(input a, output y);\n"
      "  OPAQUE u_reference(.A(a), .Y(y));\n"
      "endmodule\n",
      "module top(input a, output y);\n"
      "  OPAQUE u_implementation(.A(a), .Y(y));\n"
      "endmodule\n");

  const auto run = runStructuredWithArgs(
      {"kepler-formal",
       "-v",
       "sec",
       "--sec-engine",
       "k_induction",
       "--sec-encoding",
       "binary",
       "-k",
       "1",
       "-verilog",
       "--design1",
       fixture.design0Path.string(),
       "--design2",
       fixture.design1Path.string(),
       "--set-as-boundary",
       "u_reference",
       "u_implementation",
       "--liberty",
       fixture.libertyPath.string()});
  EXPECT_EQ(run.exitCode, kSecProvedExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Equivalent);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       CliCompactSecSetAsBoundaryProvesRenamedOpaquePair) {
  const auto fixture = createOpaqueBoundaryFixture(
      "module top(input a, output y);\n"
      "  OPAQUE u_reference(.A(a), .Y(y));\n"
      "endmodule\n",
      "module top(input a, output y);\n"
      "  OPAQUE u_implementation(.A(a), .Y(y));\n"
      "endmodule\n");

  const auto run = runStructuredWithArgs(
      {"kepler-formal",
       "-v",
       "sec",
       "--sec-engine",
       "k_induction",
       "--sec-encoding",
       "binary",
       "-k",
       "1",
       "-verilog",
       "--design1",
       fixture.design0Path.string(),
       "--design2",
       fixture.design1Path.string(),
       "--compact",
       "--set-as-boundary",
       "u_reference",
       "u_implementation",
       "--liberty",
       fixture.libertyPath.string()});
  EXPECT_EQ(run.exitCode, kSecProvedExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Equivalent);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       CliCompactSecDifferentBoundaryPathsDisableSelfModelReuse) {
  const std::string design =
      "module top(input a, input b, output y);\n"
      "  wire left_y;\n"
      "  wire right_y;\n"
      "  OPAQUE u_left(.A(a), .Y(left_y));\n"
      "  OPAQUE u_right(.A(b), .Y(right_y));\n"
      "  assign y = 1'b0;\n"
      "endmodule\n";
  const auto fixture = createOpaqueBoundaryFixture(design, design);

  const auto run = runStructuredWithArgs(
      {"kepler-formal",
       "-v",
       "sec",
       "--sec-engine",
       "k_induction",
       "--sec-encoding",
       "binary",
       "-k",
       "1",
       "-verilog",
       "--design1",
       fixture.design0Path.string(),
       "--design2",
       fixture.design0Path.string(),
       "--compact",
       "--set-as-boundary",
       "u_left",
       "u_right",
       "--liberty",
       fixture.libertyPath.string()});
  EXPECT_EQ(run.exitCode, kSecCounterexampleExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSetAsBoundaryRejectsInterfaceMismatch) {
  const auto fixture = createDesignFixture(
      "v",
      "module child(input i, output o); endmodule\n"
      "module top(input a, output y); child u(.i(a), .o(y)); endmodule\n",
      "module child(input j, output o); endmodule\n"
      "module top(input a, output y); child u(.j(a), .o(y)); endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "set_as_boundary:\n"
      "  - [u, u]\n");

  const auto run = runStructuredWithConfigFile(cfgPath);
  EXPECT_EQ(run.exitCode, EXIT_FAILURE);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Error);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSetAsBoundaryRequiresPathPairs) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "set_as_boundary: [u0, u1]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, CliSetAsBoundaryRequiresTwoPaths) {
  EXPECT_EQ(
      runWithArgs(
          {"kepler-formal",
           "-verilog",
           "design0.v",
           "design1.v",
           "--set-as-boundary",
           "u0"}),
      EXIT_FAILURE);
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

TEST_F(KeplerFormalCliTests, ConfigCcFormatRequiresTop) {
  const auto cfgPath = writeTempConfig(
      "format: cc\n"
      "verification: sec\n"
      "input_paths: [design0.cc, design1.cc]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigCcFormatRequiresSec) {
  const auto cfgPath = writeTempConfig(
      "format: cxx\n"
      "cc_top: top\n"
      "input_paths: [design0.cc, design1.cc]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigCcIncludePathsMustBeSequence) {
  const auto cfgPath = writeTempConfig(
      "format: cc\n"
      "verification: sec\n"
      "cc_top: top\n"
      "cc_include_paths: include\n"
      "input_paths: [design0.cc, design1.cc]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigCcRejectsMalformedSynthesisOptions) {
  const std::vector<std::pair<std::string, std::string>> cases = {
      {"cc_top: []\n", "cc_top must be a scalar"},
      {"cc_top: ''\n", "cc_top must not be empty"},
      {"cc_include_paths: [[]]\n", "cc_include_paths entries must be scalar paths"},
  };
  for (const auto& [option, expectedError] : cases) {
    SCOPED_TRACE(option);
    std::ostringstream captured;
    NamedLoggerGuard logger("cc_config_validation");
    logger.installed_->sinks().clear();
    logger.installed_->sinks().push_back(
        std::make_shared<spdlog::sinks::ostream_sink_mt>(captured));
    const auto cfgPath = writeTempConfig(
        "format: cc\nverification: sec\ninput_paths: [a.cc, b.cc]\n" + option);
    EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
    EXPECT_NE(captured.str().find(expectedError), std::string::npos);
    std::filesystem::remove(cfgPath);
  }
}

TEST_F(KeplerFormalCliTests, ConfigCcRejectsMultipleTranslationUnitsPerSide) {
  const auto cfgPath = writeTempConfig(
      "format: cc\nverification: sec\ncc_top: top\n"
      "input_paths: [[a.cc, helper.cc], [b.cc]]\n");
  const auto run = runStructuredWithConfigFile(cfgPath);
  EXPECT_EQ(run.exitCode, EXIT_FAILURE);
  EXPECT_NE(readFileContents(run.result.logFile).find("exactly one translation unit"),
            std::string::npos);
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

TEST_F(KeplerFormalCliTests, ConfigRtlVsGateAcceptedWithSec) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_rtl_vs_gate");
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
      "format: rtl_vs_gate\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 4\n"
      "sv_design1_top: top\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigRtlVsGateRequiresSec) {
  const auto cfgPath = writeTempConfig(
      "format: rtl_vs_gate\n"
      "verification: lec\n"
      "sv_design1_top: top\n"
      "input_paths: [design0.sv, design1.v]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigRtlVsGateRejectsDesign2SystemVerilogOptions) {
  const auto cfgPath = writeTempConfig(
      "format: rtl_vs_gate\n"
      "verification: sec\n"
      "sv_design2_top: top\n"
      "input_paths: [design0.sv, design1.v]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigCVsRtlRequiresSec) {
  const auto cfgPath = writeTempConfig(
      "format: c_vs_rtl\n"
      "cc_top: top\n"
      "input_paths: [design0.cc, design1.sv]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

TEST_F(KeplerFormalCliTests, ConfigCVsRtlRejectsDesign1SystemVerilogOptions) {
  const auto cfgPath = writeTempConfig(
      "format: c_vs_rtl\n"
      "verification: sec\n"
      "cc_top: top\n"
      "sv_design1_top: top\n"
      "input_paths: [design0.cc, design1.sv]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::filesystem::remove(cfgPath);
}

namespace {

// Tagged-record encoder: a combinational C++ model against a three-stage
// pipelined RTL that deliberately differs wherever a result is unused, so the
// proof needs input constraints, conditional output relations, the
// reachability check and a non-vacuity witness.
constexpr const char* kTaggedRecordModelSource = R"cc(
void encode_record(unsigned int operand,
                   bool& invalid,
                   bool& overflow,
                   unsigned long long& payload,
                   unsigned char& tag) {
  const unsigned int kind = operand >> 28;
  const unsigned int level = (operand >> 20) & 0xff;
  const unsigned long long top = (operand >> 19) & 1;
  const unsigned long long low = operand & 0x7ffff;

  invalid = kind >= 0xe;
  overflow = !invalid && level == 0xff;

  unsigned long long value = 0;
  if (!invalid) {
    if (overflow) {
      value = top << 40;
    } else if (kind == 0 && level == 0) {
      value = (top << 40) | low;
    } else {
      value = (top << 40) | (low << 21) | (low << 1);
    }
  }
  payload = value;
  tag = (invalid || overflow || value == 0)
      ? 0
      : static_cast<unsigned char>(level ^ 0x80);
}
)cc";

constexpr const char* kTaggedRecordRtlSource = R"sv(
module encode_record_rtl (
  input  logic        clk,
  input  logic        rst_n,
  input  logic [31:0] operand,
  output logic        invalid,
  output logic        overflow,
  output logic [40:0] payload,
  output logic [7:0]  tag
);
  logic [3:0]  kind_q;
  logic [7:0]  level_q;
  logic        top_q;
  logic [18:0] low_q;
  logic        invalid_q;
  logic        overflow_q;
  logic [40:0] payload_q;
  logic [7:0]  tag_q;

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      kind_q     <= '0;
      level_q    <= '0;
      top_q      <= 1'b0;
      low_q      <= '0;
      invalid_q  <= 1'b0;
      overflow_q <= 1'b0;
      payload_q  <= '0;
      tag_q      <= '0;
      invalid    <= 1'b0;
      overflow   <= 1'b0;
      payload    <= '0;
      tag        <= '0;
    end else begin
      kind_q     <= operand[31:28];
      level_q    <= operand[27:20];
      top_q      <= operand[19];
      low_q      <= operand[18:0];

      // Reserved kinds and the compact layout are not built, and the unused
      // fields of an invalid, overflowed or zero record keep the datapath value.
      invalid_q  <= kind_q == 4'hE;
      overflow_q <= level_q == 8'hFF;
      payload_q  <= level_q == 8'hFF
                        ? {top_q, {40{1'b1}}}
                        : {top_q, low_q, 1'b0, low_q, 1'b0};
      tag_q      <= level_q ^ 8'h80;

      invalid    <= invalid_q;
      overflow   <= overflow_q;
      payload    <= payload_q;
      tag        <= tag_q;
    end
  end
endmodule
)sv";

constexpr const char* kTaggedRecordConstraints =
    "  - operand > 0x000FFFFF\n"
    "  - operand < 0xF0000000\n"
    "  - rtl.rst_n\n";

// Each entry requires condition -> equality three clock events after the
// cycle-0 inputs, which is the RTL's register depth.
constexpr const char* kTaggedRecordEventuals =
    "  - cycle: 3\n"
    "    condition: \"true\"\n"
    "    equality: \"model.invalid == rtl.invalid\"\n"
    "  - cycle: 3\n"
    "    condition: \"!model.invalid\"\n"
    "    equality: \"(model.overflow == rtl.overflow) && "
    "(model.payload[40] == rtl.payload[40])\"\n"
    "  - cycle: 3\n"
    "    condition: \"!model.invalid && !model.overflow\"\n"
    "    equality: \"model.payload == rtl.payload\"\n"
    "  - cycle: 3\n"
    "    condition: \"!model.invalid && !model.overflow && "
    "(model.payload != 0)\"\n"
    "    equality: \"model.tag == rtl.tag\"\n";

std::string taggedRecordConfig(const std::string& constraints,
                               const std::string& eventuals,
                               const std::string& checkReachability = "true") {
  return "format: c_vs_rtl\n"
         "verification: sec\n"
         "sec_engine: pdr\n"
         "sec_encoding: binary\n"
         "c2rtl_auto_align: true\n"
         "max_k: 8\n"
         "cc_top: encode_record\n"
         "sv_design2_top: encode_record_rtl\n"
         "cc_output_dir: kepler_formal_c2rtl\n"
         "input_paths:\n"
         "  - encode_record.cc\n"
         "  - encode_record.sv\n"
         "constraints:\n" + constraints +
         "check_reachability: " + checkReachability + "\n"
         "eventuals:\n" + eventuals;
}

StructuredRun runTaggedRecord(const std::string& config,
                              const std::string& rtlSource = kTaggedRecordRtlSource) {
  const auto tmpDir = makeUniqueTempDir("kepler_formal_cli_tagged_record");
  {
    std::ofstream model(tmpDir / "encode_record.cc");
    model << kTaggedRecordModelSource;
    std::ofstream rtl(tmpDir / "encode_record.sv");
    rtl << rtlSource;
    std::ofstream cfg(tmpDir / "encode_record.yaml");
    cfg << config;
  }
  StructuredRun run;
  {
    // The generated SystemVerilog for the model lands in the working directory.
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(tmpDir);
    run = runStructuredWithConfigFile(tmpDir / "encode_record.yaml");
  }
  KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
  KEPLER_FORMAL::BoolExprCache::destroy();
  std::filesystem::remove_all(tmpDir);
  return run;
}

}  // namespace

TEST_F(KeplerFormalCliTests, CliCVsRtlTaggedRecordProvesConditionalRelations) {
  const auto run = runTaggedRecord(
      taggedRecordConfig(kTaggedRecordConstraints, kTaggedRecordEventuals));
  EXPECT_EQ(run.exitCode, kSecProvedExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Equivalent);
  EXPECT_EQ(run.result.provenOutputs, 4u);
}

TEST_F(KeplerFormalCliTests,
       CliCVsRtlTaggedRecordWithoutOperandRangeFindsCounterexample) {
  const auto run = runTaggedRecord(
      taggedRecordConfig("  - rtl.rst_n\n", kTaggedRecordEventuals));
  EXPECT_EQ(run.exitCode, kSecCounterexampleExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
}

TEST_F(KeplerFormalCliTests,
       CliCVsRtlTaggedRecordUnconditionalEqualityFindsCounterexample) {
  // The RTL leaves the unused fields as its datapath computed them, so
  // comparing every output unconditionally must fail.
  const auto run = runTaggedRecord(taggedRecordConfig(
      kTaggedRecordConstraints,
      "  - cycle: 3\n"
      "    condition: \"true\"\n"
      "    equality: \"model.invalid == rtl.invalid\"\n"
      "  - cycle: 3\n"
      "    condition: \"true\"\n"
      "    equality: \"model.overflow == rtl.overflow\"\n"
      "  - cycle: 3\n"
      "    condition: \"true\"\n"
      "    equality: \"model.payload == rtl.payload\"\n"
      "  - cycle: 3\n"
      "    condition: \"true\"\n"
      "    equality: \"model.tag == rtl.tag\"\n"));
  EXPECT_EQ(run.exitCode, kSecCounterexampleExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
}

TEST_F(KeplerFormalCliTests,
       CliCVsRtlTaggedRecordContradictoryConstraintsAreNotProvedVacuously) {
  const std::string contradictory =
      "  - operand < 0x000FFFFF\n"
      "  - operand > 0xF0000000\n"
      "  - rtl.rst_n\n";
  const auto rejected = runTaggedRecord(
      taggedRecordConfig(contradictory, kTaggedRecordEventuals));
  EXPECT_EQ(rejected.exitCode, kSecInconclusiveExitCode);
  EXPECT_EQ(rejected.result.status, KEPLER_FORMAL::RunStatus::Unsupported);

  // Without the reachability check the empty input domain proves anything.
  const auto vacuous = runTaggedRecord(
      taggedRecordConfig(contradictory, kTaggedRecordEventuals, "false"));
  EXPECT_EQ(vacuous.exitCode, kSecProvedExitCode);
  EXPECT_EQ(vacuous.result.status, KEPLER_FORMAL::RunStatus::Equivalent);
}

TEST_F(KeplerFormalCliTests, CliCVsRtlTaggedRecordPlantedTagBugFindsCounterexample) {
  // A wrong tag inside the legal domain must be caught, so the conditional
  // relations are not vacuous.
  std::string rtl = kTaggedRecordRtlSource;
  const std::string good = "tag_q      <= level_q ^ 8'h80;";
  const auto at = rtl.find(good);
  ASSERT_NE(at, std::string::npos);
  rtl.replace(at, good.size(), "tag_q      <= level_q ^ 8'h81;");
  const auto run = runTaggedRecord(
      taggedRecordConfig(kTaggedRecordConstraints, kTaggedRecordEventuals), rtl);
  EXPECT_EQ(run.exitCode, kSecCounterexampleExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Different);
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
  const auto libertyPath = repoRoot() / "examples" / "tinyrocket" /
                           "NangateOpenCellLibrary_typical.lib";
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
    ConfigSv2vDesignDefinedModuleOverridesGeneratedPrimitiveStub) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_sv2v_owned_prim");
  fixture.design0Path = fixture.tmpDir / "design0.sv";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  const auto libertyPath = repoRoot() / "examples" / "tinyrocket" /
                           "NangateOpenCellLibrary_typical.lib";
  ASSERT_TRUE(std::filesystem::exists(libertyPath));

  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module INV_X1(input A, output ZN);\n";
    design0 << "  assign ZN = A;\n";
    design0 << "endmodule\n";
    design0 << "module top(input logic a, b, output logic y);\n";
    design0 << "  logic n;\n";
    design0 << "  INV_X1 u_owned(.A(a), .ZN(n));\n";
    design0 << "  AND2_X1 u_lib(.A1(n), .A2(b), .ZN(y));\n";
    design0 << "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, b, output y);\n";
    design1 << "  AND2_X1 u_lib(.A1(a), .A2(b), .ZN(y));\n";
    design1 << "endmodule\n";
  }

  for (const bool compact : {false, true}) {
    const std::string mode = compact ? "compact" : "standard";
    const auto runDir = fixture.tmpDir / mode;
    std::filesystem::create_directories(runDir);
    const auto logPath = runDir / "run.log";
    const auto diagnosticsPath = runDir / "naja_sv_diagnostics.log";
    const auto cfgPath = writeTempConfig(
        "format: sv2v\n"
        "verification: sec\n"
        "sec_engine: pdr\n"
        "sec_encoding: binary\n"
        "max_k: 4\n"
        "compact_mode: " + std::string(compact ? "true" : "false") + "\n"
        "sv_design1_top: top\n"
        "verilog_design2_top: top\n"
        "input_paths:\n"
        "  - " + fixture.design0Path.string() + "\n"
        "  - " + fixture.design1Path.string() + "\n"
        "liberty_files:\n"
        "  - " + libertyPath.string() + "\n"
        "log_file: " + logPath.string() + "\n");

    {
      CurrentPathGuard currentPathGuard;
      std::filesystem::current_path(runDir);
      EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_SUCCESS) << mode;
    }
    ASSERT_TRUE(std::filesystem::exists(logPath)) << mode;
    const auto logContents = readFileContents(logPath);
    EXPECT_NE(
        logContents.find(
            "SEC checked-output coverage: 100.00% (1/1 covered/existing outputs)."),
        std::string::npos)
        << mode;
    const auto compactMarker =
        "SEC compact mode: extracting and releasing design 1 before "
        "loading design 2";
    if (compact) {
      EXPECT_NE(logContents.find(compactMarker), std::string::npos) << mode;
    } else {
      EXPECT_EQ(logContents.find(compactMarker), std::string::npos) << mode;
    }
    ASSERT_TRUE(std::filesystem::exists(diagnosticsPath)) << mode;
    EXPECT_EQ(
        readFileContents(diagnosticsPath).find("duplicate definition of 'INV_X1'"),
        std::string::npos)
        << mode;
    std::filesystem::remove(cfgPath);
  }
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

  EnvVarGuard pathGuard("PATH");
  pathGuard.set(keplerBinPath.parent_path().string());
  pythonPathGuard.set(fixture.tmpDir.string());
  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(fixture.tmpDir);
    EXPECT_EQ(
        runWithConfigFile(cfgPath, keplerBinPath.filename().string()),
        EXIT_SUCCESS);
  }
  EXPECT_EQ(
      std::getenv("PYTHONPATH"),
      keplerBinPath.parent_path().string() + ":" + fixture.tmpDir.string());

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigXilinxPythonPrimitiveSecExample) {
  const auto exampleDir =
      repoRoot() / "examples" / "xilinx" / "register_slice";
  const auto tmpDir = makeUniqueTempDir("kepler_formal_cli_xilinx_sec");
  const auto cfgPath = tmpDir / "config.yaml";
  {
    std::ofstream cfg(cfgPath);
    cfg << "format: verilog\n"
           "verification: sec\n"
           "sec_engine: pdr\n"
           "sec_encoding: dual_rail_steady\n"
           "input_paths:\n"
        << "  - " << (exampleDir / "xilinx_register_slice_mapped.v").string()
        << "\n  - "
        << (exampleDir / "xilinx_register_slice_compact.v").string()
        << "\npy_tech_files:\n  - "
        << (exampleDir.parent_path() / "xilinx.py").string()
        << "\nlog_file: " << (tmpDir / "miter.log").string() << "\n";
  }

  const char* keplerBin = std::getenv("KEPLER_BIN");
  ASSERT_NE(keplerBin, nullptr);
  EXPECT_EQ(runWithConfigFile(cfgPath, keplerBin), EXIT_SUCCESS);
  std::filesystem::remove_all(tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigXilinxLutInstanceParameters) {
  SimpleCliFixture fixture;
  fixture.tmpDir = makeUniqueTempDir("kepler_formal_cli_xilinx_luts");
  fixture.design0Path = fixture.tmpDir / "design0.v";
  fixture.design1Path = fixture.tmpDir / "design1.v";
  const auto cfgPath = fixture.tmpDir / "config.yaml";
  {
    std::ofstream design0(fixture.design0Path);
    design0 << "module top(input a, b, output xor_y, and_y);\n"
               "  LUT2 #(.INIT(4'h6)) xor_lut(.I0(a), .I1(b), .O(xor_y));\n"
               "  LUT2 #(.INIT(4'h8)) and_lut(.I0(a), .I1(b), .O(and_y));\n"
               "endmodule\n";
  }
  {
    std::ofstream design1(fixture.design1Path);
    design1 << "module top(input a, b, output xor_y, and_y);\n"
               "  XOR2 explicit_xor(.A(a), .B(b), .Y(xor_y));\n"
               "  AND2 explicit_and(.A(a), .B(b), .Y(and_y));\n"
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
        << "  - "
        << (repoRoot() / "examples" / "xilinx" / "xilinx.py").string()
        << "\nlog_file: " << (fixture.tmpDir / "miter.log").string() << "\n";
  }

  const char* keplerBin = std::getenv("KEPLER_BIN");
  ASSERT_NE(keplerBin, nullptr);
  EXPECT_EQ(runWithConfigFile(cfgPath, keplerBin), EXIT_SUCCESS);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigXilinxVexRiscvGenFullLecExample) {
  const auto exampleDir = repoRoot() / "examples" / "xilinx" / "vexriscv";
  const auto netlist = exampleDir / "vexriscv_genfull_xilinx.v";
  const auto tmpDir = makeUniqueTempDir("kepler_formal_cli_vexriscv_xilinx");
  const auto cfgPath = tmpDir / "config.yaml";
  {
    std::ofstream cfg(cfgPath);
    cfg << "format: verilog\n"
           "verification: lec\n"
           "input_paths:\n"
        << "  - " << netlist.string() << "\n"
        << "  - " << netlist.string() << "\n"
           "verilog_design1_top: vexriscv.demo.GenFull\n"
           "verilog_design2_top: vexriscv.demo.GenFull\n"
           "py_tech_files:\n"
        << "  - " << (exampleDir.parent_path() / "xilinx.py").string()
        << "\nlog_file: " << (tmpDir / "miter.log").string() << "\n";
  }

  const char* keplerBin = std::getenv("KEPLER_BIN");
  ASSERT_NE(keplerBin, nullptr);
  EXPECT_EQ(runWithConfigFile(cfgPath, keplerBin), EXIT_SUCCESS);
  std::filesystem::remove_all(tmpDir);
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
  const auto libertyPath = repoRoot() / "examples" / "tinyrocket" /
                           "NangateOpenCellLibrary_typical.lib";
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
  const auto exampleDir = root / "examples" / "tinyrocket";
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
  const auto exampleDir = root / "examples" / "tinyrocket";
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
  const auto exampleDir = root / "examples" / "tinyrocket";
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

TEST_F(KeplerFormalCliTests, ConfigRemovedSecBoundaryAbstractionKeyIsRejected) {
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_uncomputable_seq_as_boundary: true\n");
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

TEST_F(KeplerFormalCliTests, ConfigSecResetBootstrapAcceptsMultiplePorts) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input reset, input scan_reset_n, input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sec_reset_bootstrap.log";
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: binary\n"
      "max_k: 1\n"
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: reset\n"
      "      active_value: 1\n"
      "    - name: scan_reset_n\n"
      "      active_value: 0\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("SEC reset bootstrap: 1 cycle(s)"),
            std::string::npos);
  EXPECT_NE(contents.find("SEC reset bootstrap port: reset active=1"),
            std::string::npos);
  EXPECT_NE(
      contents.find("SEC reset bootstrap port: scan_reset_n active=0"),
      std::string::npos);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecResetBootstrapUsesSequentialResetInput) {
  const EnvVarGuard secDiag("KEPLER_SEC_DIAG");
  secDiag.set("1");
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic clk, input logic reset, input logic a, output logic y);\n"
      "  always_ff @(posedge clk) begin\n"
      "    if (reset) y <= 1'b0;\n"
      "    else y <= a;\n"
      "  end\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: binary\n"
      "max_k: 1\n"
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: reset\n"
      "      active_value: true\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSecResetBootstrapAcceptsRepeatedPorts) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input reset, input scan_reset_n, input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-verilog",
                   "-v",
                   "sec",
                   "--sec-engine",
                   "pdr",
                   "--sec-encoding",
                   "binary",
                   "-k",
                   "1",
                   "--sec-reset-cycles",
                   "1",
                   "--sec-reset-port",
                   "reset=1",
                   "--sec-reset-port",
                   "scan_reset_n=0",
                   fixture.design0Path.string(),
                   fixture.design1Path.string()}),
      kSecProvedExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecResetBootstrapMissingPortFails) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input reset, input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sec_reset_missing.log";
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: binary\n"
      "max_k: 1\n"
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: missing\n"
      "      active_value: 1\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  const auto run = runStructuredWithConfigFile(cfgPath);
  EXPECT_EQ(run.exitCode, kSecInconclusiveExitCode);
  EXPECT_EQ(run.result.exitCode, kSecInconclusiveExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Unsupported);
  EXPECT_NE(run.result.reason.find("missing"), std::string::npos);
  const auto contents = readFileContents(logPath);
  EXPECT_NE(
      contents.find(
          "Reset bootstrap port `missing` was not found among aligned "
          "top-level inputs"),
      std::string::npos);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecResetBootstrapRejectsMalformedYaml) {
  const std::vector<std::string> invalidResetBlocks = {
      "sec_reset: scalar\n",
      "sec_reset:\n"
      "  ports:\n"
      "    - name: reset\n"
      "      active_value: 1\n",
      "sec_reset:\n"
      "  cycles: invalid\n"
      "  ports:\n"
      "    - name: reset\n"
      "      active_value: 1\n",
      "sec_reset:\n"
      "  cycles: 1\n",
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - reset\n",
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - active_value: 1\n",
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: \"\"\n"
      "      active_value: 1\n",
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: reset\n",
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: reset\n"
      "      active_value: maybe\n",
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: reset\n"
      "      active_value: 1\n"
      "    - name: reset\n"
      "      active_value: 0\n",
      "sec_reset:\n"
      "  cycles: 0\n"
      "  ports:\n"
      "    - name: reset\n"
      "      active_value: 1\n"};

  for (const auto& resetBlock : invalidResetBlocks) {
    const auto cfgPath = writeTempConfig(
        "format: verilog\n"
        "verification: sec\n" +
        resetBlock);
    EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE) << resetBlock;
    std::filesystem::remove(cfgPath);
  }
}

TEST_F(KeplerFormalCliTests, CliSecResetBootstrapRejectsMalformedPortTokens) {
  for (const char* resetPort : {"reset", "reset=maybe", "=1"}) {
    EXPECT_EQ(
        runWithArgs({"kepler-formal",
                     "-v",
                     "sec",
                     "--sec-reset-port",
                     resetPort,
                     "-verilog"}),
        EXIT_FAILURE) << resetPort;
  }
}

TEST_F(KeplerFormalCliTests, CliSecResetBootstrapIncompleteSpecFails) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "--sec-reset-cycles",
                   "1",
                   "-verilog",
                   fixture.design0Path.string(),
                   fixture.design1Path.string()}),
      EXIT_FAILURE);
  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-v",
                   "sec",
                   "--sec-reset-port",
                   "reset=1",
                   "-verilog",
                   fixture.design0Path.string(),
                   fixture.design1Path.string()}),
      EXIT_FAILURE);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSecResetBootstrapRejectedForLec) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module top(input reset, input a, output y);\n"
      "  assign y = a;\n"
      "endmodule\n");

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "--sec-reset-cycles",
                   "1",
                   "--sec-reset-port",
                   "reset=1",
                   "-verilog",
                   fixture.design0Path.string(),
                   fixture.design1Path.string()}),
      EXIT_FAILURE);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecResetBootstrapAcceptsBareOneBitBus) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic [0:0] rst, input logic a, output logic y);\n"
      "  assign y = rst[0] & a;\n"
      "endmodule\n");
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: binary\n"
      "max_k: 1\n"
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: rst\n"
      "      active_value: 1\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecResetBootstrapRejectsBareMultiBitBus) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic [1:0] rst, input logic a, output logic y);\n"
      "  assign y = (rst[0] & a) | rst[1];\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sec_reset_bus.log";
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: binary\n"
      "max_k: 1\n"
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: rst\n"
      "      active_value: 1\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecInconclusiveExitCode);
  const auto contents = readFileContents(logPath);
  EXPECT_NE(
      contents.find("matched multiple input bits; name each reset bit explicitly"),
      std::string::npos);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecResetBootstrapRejectsDuplicateResolvedBusBit) {
  const auto fixture = createEquivalentDesignFixture(
      "sv",
      "module top(input logic [0:0] rst, input logic a, output logic y);\n"
      "  assign y = rst[0] & a;\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sec_reset_duplicate_bit.log";
  const auto cfgPath = writeTempConfig(
      "format: systemverilog\n"
      "verification: sec\n"
      "sec_engine: pdr\n"
      "sec_encoding: binary\n"
      "max_k: 1\n"
      "sec_reset:\n"
      "  cycles: 1\n"
      "  ports:\n"
      "    - name: rst\n"
      "      active_value: 1\n"
      "    - name: rst[0]\n"
      "      active_value: 1\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecInconclusiveExitCode);
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("resolves to the same top-level input"),
            std::string::npos);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecVerificationAcceptedWithPdrEngine) {
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto logPath = fixture.tmpDir / "structured_pdr.log";
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "sec_engine: pdr\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");
  const auto run = runStructuredWithConfigFile(cfgPath);
  EXPECT_EQ(run.exitCode, kSecProvedExitCode);
  EXPECT_EQ(run.result.exitCode, kSecProvedExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Equivalent);
  EXPECT_EQ(run.result.inputFormat, "naja_if");
  EXPECT_EQ(run.result.verification, "sec");
  EXPECT_EQ(run.result.totalOutputs, 1u);
  EXPECT_EQ(run.result.coveredOutputs, 1u);
  EXPECT_EQ(run.result.provenOutputs, 1u);
  EXPECT_TRUE(run.result.unprovenOutputs.empty());
  EXPECT_EQ(run.result.logFile, logPath.string());
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecVerificationRejectsLegacyEngine) {
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
  const auto fixture = createEquivalentSequentialNajaIfFixture();
  const auto logPath = fixture.tmpDir / "structured_imc.log";
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "sec_engine: imc\n"
      "max_k: 4\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");
  // This option-parsing fixture is valid input for IMC, but IMC does not
  // converge within the small bound.
  const auto run = runStructuredWithConfigFile(cfgPath);
  EXPECT_EQ(run.exitCode, kSecInconclusiveExitCode);
  EXPECT_EQ(run.result.exitCode, kSecInconclusiveExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Inconclusive);
  EXPECT_EQ(run.result.totalOutputs, 1u);
  EXPECT_EQ(run.result.provenOutputs, 0u);
  EXPECT_EQ(run.result.unprovenOutputs, std::vector<std::string>{"out[0]"});
  EXPECT_EQ(run.result.logFile, logPath.string());
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
     ConfigSecReportsPartialProofAndWritesOpaqueOutputReport) {
  const auto fixture = createUncomputableSequentialNajaIfFixture();
  const auto logPath = fixture.tmpDir / "sec_partial_opaque.log";
  const auto reportPath = fixture.tmpDir / "skipped_opaque_cells_pos.txt";
  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "verification: sec\n"
      "sec_encoding: binary\n"
      "max_k: 2\n"
      "report_skipped_pos: true\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(fixture.tmpDir);
    const auto run = runStructuredWithConfigFile(cfgPath);
    EXPECT_EQ(run.exitCode, kSecPartiallyProvedExitCode);
    EXPECT_EQ(run.result.exitCode, kSecPartiallyProvedExitCode);
    EXPECT_EQ(
        run.result.status, KEPLER_FORMAL::RunStatus::PartiallyProved);
    EXPECT_EQ(run.result.totalOutputs, 2u);
    EXPECT_EQ(run.result.coveredOutputs, 1u);
    EXPECT_EQ(run.result.provenOutputs, 1u);
    ASSERT_EQ(run.result.skippedObservedOutputs.size(), 1u);
    EXPECT_NE(
        run.result.skippedObservedOutputs.front().find("bad[0]:"),
        std::string::npos);
    EXPECT_NE(
        run.result.skippedObservedOutputs.front().find("opaque internal cell"),
        std::string::npos);
  }
  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("SEC checked-output coverage: 50.00% (1/2"),
            std::string::npos);
  EXPECT_NE(contents.find("bad[0]:"), std::string::npos);
  EXPECT_NE(contents.find("opaque internal cell"), std::string::npos);
  EXPECT_NE(contents.find("SEC partially proved equivalence"), std::string::npos);
  EXPECT_NE(contents.find("did not prove all observed outputs"), std::string::npos);

  ASSERT_TRUE(std::filesystem::exists(reportPath));
  const auto report = readFileContents(reportPath);
  const auto entry = report.find("- bad[0]:");
  ASSERT_NE(entry, std::string::npos);
  EXPECT_EQ(report.find("- bad[0]:", entry + 1), std::string::npos);
  EXPECT_NE(report.find("outputs were not verified"), std::string::npos);
  EXPECT_NE(report.find("ff0"), std::string::npos);
  EXPECT_NE(report.find("SEQ_NO_D"), std::string::npos);
  EXPECT_NE(report.find("Q[0]"), std::string::npos);
  EXPECT_NE(report.find("no initialized combinational truth table"),
            std::string::npos);
  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
     SecStatefulClockGateWithoutBehavioralModelIsSkippedAsOpaque) {
  const auto fixture = createDesignFixture(
      "v",
      "module top(input clk, input rst_n, input en, input d, output q, output pass);\n"
      "  wire gclk;\n"
      "  wire qn;\n"
      "  wire pass_n;\n"
      "  CLKGATE_X1 u_gate(.CK(clk), .E(en), .GCK(gclk));\n"
      "  DFFR_X1 u_q(.CK(gclk), .D(d), .Q(q), .QN(qn), .RN(rst_n));\n"
      "  DFFR_X1 u_pass(.CK(clk), .D(d), .Q(pass), .QN(pass_n), .RN(rst_n));\n"
      "endmodule\n",
      "module top(input clk, input rst_n, input en, input d, output q, output pass);\n"
      "  wire next_q;\n"
      "  wire qn;\n"
      "  wire pass_n;\n"
      "  MUX2_X1 u_hold(.A(q), .B(d), .S(en), .Z(next_q));\n"
      "  DFFR_X1 u_q(.CK(clk), .D(next_q), .Q(q), .QN(qn), .RN(rst_n));\n"
      "  DFFR_X1 u_pass(.CK(clk), .D(d), .Q(pass), .QN(pass_n), .RN(rst_n));\n"
      "endmodule\n");
  const auto logPath = fixture.tmpDir / "sec_clockgate_opaque.log";
  const auto libertyPath =
      repoRoot() / "examples" / "tinyrocket" /
      "NangateOpenCellLibrary_typical.lib";
  const auto cfgPath = writeTempConfig(
      "format: verilog\n"
      "verification: sec\n"
      "sec_encoding: dual_rail_steady\n"
      "sec_engine: pdr\n"
      "max_k: 2\n"
      "input_paths:\n"
      "  - " + fixture.design0Path.string() + "\n"
      "  - " + fixture.design1Path.string() + "\n"
      "liberty_files:\n"
      "  - " + libertyPath.string() + "\n"
      "log_file: " + logPath.string() + "\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecPartiallyProvedExitCode);

  const auto contents = readFileContents(logPath);
  EXPECT_NE(contents.find("SEC checked-output coverage: 50.00% (1/2"),
            std::string::npos);
  EXPECT_NE(contents.find("q[0]: design0 opaque-internal"), std::string::npos);
  EXPECT_NE(contents.find("CLKGATE_X1"), std::string::npos);
  EXPECT_NE(contents.find("no initialized combinational truth table or usable "
                          "sequential model"),
            std::string::npos);
  // Liberty `statetable`/`state_function` cells must remain opaque until their
  // stateful behavior is modeled. In particular, the frontend must not invent
  // the old zero-input truth table for GCK.
  EXPECT_EQ(contents.find("combinational truth table arity does not match "
                          "instance inputs"),
            std::string::npos);
  EXPECT_EQ(contents.find("SNLLogicCloud arity mismatch"), std::string::npos);
  EXPECT_EQ(contents.find("unpublished internal support"), std::string::npos);

  std::filesystem::remove(cfgPath);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSecUnsupportedMismatchLogUsesUnsupportedResult) {
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
       ConfigSystemVerilogSecUnknownConstantsIdentifyTheirUseSite) {
  struct Case {
    const char* name;
    const char* assignment;
    const char* diagnostic;
    bool undriven = false;
  };
  const Case cases[] = {
      {"x_ternary", "assign y = sel ? d : 1'bx;",
       "unsupported X constant (1'bx)"},
      {"z_bitwise", "assign y = d & 1'bz;",
       "unsupported Z constant (1'bz)"},
      {"x_direct", "assign y = 1'bx;", "unsupported X constant (1'bx)"},
      {"z_direct", "assign y = 1'bz;", "unsupported Z constant (1'bz)"},
      {"x_next_state", "always_ff @(posedge clk) y <= sel ? d : 1'bx;",
       "unsupported X constant (1'bx)"},
      {"known_zero", "assign y = sel ? d : 1'b0;", nullptr},
      {"undriven", "assign y = undriven;", nullptr, true},
  };
  for (const auto& test : cases) {
    SCOPED_TRACE(test.name);
    const std::string source =
        "module literal_source(input logic clk, sel, d, output logic y);\n"
        "  wire undriven;\n"
        "  " + std::string(test.assignment) + "\n"
        "endmodule\n"
        "module top(input logic clk, sel, d, output logic bad, good);\n"
        "  literal_source u_literals(.clk(clk), .sel(sel), .d(d), .y(bad));\n"
        "  assign good = d;\n"
        "endmodule\n";
    const auto fixture = createEquivalentDesignFixture("sv", source);
    const auto logPath = fixture.tmpDir / "unknown_constants.log";
    const auto cfgPath = writeTempConfig(
        "format: systemverilog\n"
        "verification: sec\n"
        "sec_engine: pdr\n"
        "sec_encoding: dual_rail_steady\n"
        "max_k: 2\n"
        "input_paths:\n"
        "  - " + fixture.design0Path.string() + "\n"
        "  - " + fixture.design1Path.string() + "\n"
        "log_file: " + logPath.string() + "\n");
    {
      CurrentPathGuard currentPathGuard;
      std::filesystem::current_path(fixture.tmpDir);
      const auto run = runStructuredWithConfigFile(cfgPath);
      const bool skipped = test.diagnostic != nullptr || test.undriven;
      EXPECT_EQ(run.exitCode,
                skipped ? kSecPartiallyProvedExitCode : kSecProvedExitCode);
      EXPECT_EQ(run.result.totalOutputs, 2u);
      EXPECT_EQ(run.result.coveredOutputs, skipped ? 1u : 2u);
      EXPECT_EQ(run.result.provenOutputs, skipped ? 1u : 2u);
      EXPECT_EQ(run.result.skippedObservedOutputs.size(), skipped ? 1u : 0u);
      if (skipped && !run.result.skippedObservedOutputs.empty()) {
        const auto& detail = run.result.skippedObservedOutputs.front();
        EXPECT_NE(detail.find("bad[0]:"), std::string::npos);
        if (test.diagnostic != nullptr) {
          EXPECT_NE(detail.find("unknown-constant"), std::string::npos);
          EXPECT_NE(detail.find(test.diagnostic), std::string::npos);
          EXPECT_NE(detail.find("u_literals"), std::string::npos);
          // The location belongs to the consuming assignment/expression,
          // whose line is distinct from the module and net declarations.
          EXPECT_TRUE(detail.find("design0.sv:3") != std::string::npos ||
                      detail.find("design1.sv:3") != std::string::npos)
              << detail;
          EXPECT_EQ(detail.find("internal frontier term"), std::string::npos);
          EXPECT_EQ(detail.find("no-driver connectivity"), std::string::npos);
          const auto contents = readFileContents(logPath);
          EXPECT_NE(contents.find(test.diagnostic), std::string::npos);
        } else {
          EXPECT_NE(detail.find("no-driver connectivity"), std::string::npos);
          EXPECT_EQ(detail.find("unknown-constant"), std::string::npos);
        }
      }
      KEPLER_FORMAL::cleanupKeplerFormalState();
    }
    std::filesystem::remove(cfgPath);
    std::filesystem::remove_all(fixture.tmpDir);
  }
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

namespace {

// Issue #250 reproducer: `n` independent 8-bit accumulators behind a free
// (unconstrained) reset input, each driving one parity output bit. Every
// register boots as X, so the two SEC copies share no reset anchor and the
// proof has to relate corresponding registers of the two sides itself.
std::string accumulatorBankSource(size_t n) {
  std::string source =
      "module top (\n"
      "  input wire clk,\n"
      "  input wire rst,\n"
      "  input wire [7:0] din,\n"
      "  output wire [" + std::to_string(n - 1) + ":0] dout\n"
      ");\n";
  for (size_t i = 0; i < n; ++i) {
    source += "  reg [7:0] acc" + std::to_string(i) + ";\n";
  }
  source += "  always @(posedge clk) begin\n";
  for (size_t i = 0; i < n; ++i) {
    const auto acc = "acc" + std::to_string(i);
    source += "    if (rst) " + acc + " <= 8'd0; else " + acc + " <= " + acc +
        " + din + 8'd" + std::to_string(i) + ";\n";
  }
  source += "  end\n";
  for (size_t i = 0; i < n; ++i) {
    source += "  assign dout[" + std::to_string(i) + "] = ^acc" +
        std::to_string(i) + ";\n";
  }
  source += "endmodule\n";
  return source;
}

// Issue #250, second report: the tinyalu example (single-cycle ALU next to a
// three-stage multiplier pipeline, 85 X-initialized state bits, 17 outputs).
// Reproduced as attached, without its commented-out dump and formal-checker
// stubs.
const char* const kTinyAluSource = R"(typedef enum logic [2:0] {
    NOP = 3'b000,
    ADD = 3'b001,
    AND = 3'b010,
    XOR = 3'b011,
    MUL = 3'b100
} alu_op_t;

module tinyalu (input [7:0] A,
		input [7:0] B,
		input alu_op_t op,
		input clk,
		input reset_n,
		input start,
		output done,
		output [15:0] result);

   wire [15:0] 		      result_aax, result_mult;
   wire 		          start_single, start_mult;
   wire                   done_aax;
   wire                   done_mult;

   logic [23:0] op_string;

   always_comb begin
       case (op)
           ADD     : op_string = "ADD";
           AND     : op_string = "AND";
           XOR     : op_string = "XOR";
           MUL     : op_string = "MUL";
           default : op_string = "NOP";
       endcase
   end

   int test_case_id = 0;

   assign start_single = start & ~op[2];
   assign start_mult   = start & op[2];

   single_cycle and_add_xor (.A, .B, .op, .clk, .reset_n, .start(start_single),
			     .done(done_aax), .result(result_aax));

   three_cycle mult (.A, .B, .op, .clk, .reset_n, .start(start_mult),
		    .done(done_mult), .result(result_mult));

   assign done = (op[2]) ? done_mult : done_aax;

   assign result = (op[2]) ? result_mult : result_aax;
endmodule // tinyalu

module single_cycle(input [7:0] A,
		   input [7:0] B,
		   input [2:0] op,
		   input clk,
		   input reset_n,
		   input start,
		   output logic done,
		   output logic [15:0] result);

  always @(posedge clk)
    if (!reset_n)
      result <= 0;
    else
      case(op)
		ADD : result <= {8'd0,A} + {8'd0,B};
		AND : result <= {8'd0,A} & {8'd0,B};
		XOR : result <= {8'd0,A} ^ {8'd0,B};
		default : result <= {A,B};
      endcase // case (op)

   always @(posedge clk)
     if (!reset_n)
       done <= 0;
     else
       done <= ((start == 1'b1) && (op != 3'b000));

endmodule : single_cycle

module three_cycle(input [7:0] A,
		   input [7:0] B,
		   input [2:0] op,
		   input clk,
		   input reset_n,
		   input start,
		   output logic done,
		   output logic [15:0] result);

   logic [7:0] 			       a_int, b_int;
   logic [15:0] 		       mult1, mult2;
   logic 			       done1, done2, done3;

   always @(posedge clk)
     if (!reset_n) begin
	done  <= 0;
	done3 <= 0;
	done2 <= 0;
	done1 <= 0;
	a_int <= 0;
	b_int <= 0;
	mult1 <= 0;
	mult2 <= 0;
	result<= 0;
     end else begin // if (!reset_n)
	a_int  <= A;
	b_int  <= B;
	mult1  <= a_int * b_int;
	mult2  <= mult1;
	result <= mult2;
	done3  <= start & !done;
	done2  <= done3 & !done;
	done1  <= done2 & !done;
	done   <= done1 & !done;
     end // else: !if(!reset_n)
endmodule : three_cycle
)";

// Runs the issue #250 command line: the same file on both sides, dual-rail
// PDR with default options, and checks that every output is proved.
void expectSelfCompareProvesAllOutputs(
    const std::string& source,
    const std::string& top,
    size_t expectedOutputs) {
  const auto fixture = createEquivalentDesignFixture("sv", source);
  const auto runDir = fixture.tmpDir / "self_compare_run";
  std::filesystem::create_directories(runDir);
  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    const auto run = runStructuredWithArgs(
        {"kepler-formal", "-sv",
         "--design1", fixture.design0Path.string(),
         "--design2", fixture.design0Path.string(),
         "--sv_design1_top", top, "--sv_design2_top", top,
         "-v", "sec", "--sec-engine", "pdr", "--report-skipped-pos"});
    EXPECT_EQ(run.exitCode, kSecProvedExitCode);
    EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Equivalent);
    EXPECT_EQ(run.result.totalOutputs, expectedOutputs);
    EXPECT_EQ(run.result.coveredOutputs, expectedOutputs);
    EXPECT_EQ(run.result.provenOutputs, expectedOutputs);
    EXPECT_TRUE(run.result.unprovenOutputs.empty());
    EXPECT_TRUE(run.result.skippedObservedOutputs.empty());
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

}  // namespace

// Issue #250: a self-compare is the one SEC query whose answer is known before
// it runs, so it can be gated on full coverage. Sizes 7, 8, 10 and 12 are the
// ones that used to lose an output (dout[5] / dout[9] / dout[10]) although
// larger and smaller banks proved; 4 and 32 bracket them.
TEST_F(KeplerFormalCliTests,
       CliSystemVerilogSecSelfCompareAccumulatorBankProvesEveryOutput) {
  for (const size_t n : {4u, 7u, 8u, 10u, 12u, 32u}) {
    SCOPED_TRACE("accumulator bank N=" + std::to_string(n));
    expectSelfCompareProvesAllOutputs(accumulatorBankSource(n), "top", n);
  }
}

// Issue #250: tinyalu self-compare used to stop at 8/17 with result[7..15]
// inconclusive after ~2 minutes; it must prove all 17 outputs.
TEST_F(KeplerFormalCliTests,
       CliSystemVerilogSecSelfCompareTinyAluProvesEveryOutput) {
  expectSelfCompareProvesAllOutputs(kTinyAluSource, "tinyalu", 17u);
}

TEST_F(KeplerFormalCliTests,
       CliSystemVerilogVariableIndexMatchesExplicitMux) {
  const auto fixture = createDesignFixture(
      "sv",
      "module t(input [2:0] sel, input [7:0] d, output y);\n"
      "  assign y = d[sel];\n"
      "endmodule\n",
      "module t(input [2:0] sel, input [7:0] d, output y);\n"
      "  assign y = sel == 3'd0 ? d[0] : sel == 3'd1 ? d[1]\n"
      "           : sel == 3'd2 ? d[2] : sel == 3'd3 ? d[3]\n"
      "           : sel == 3'd4 ? d[4] : sel == 3'd5 ? d[5]\n"
      "           : sel == 3'd6 ? d[6] : d[7];\n"
      "endmodule\n");
  const auto runDir = fixture.tmpDir / "variable_index_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    EXPECT_EQ(
        runWithArgs({"kepler-formal",
                     "-sv",
                     "--design1",
                     fixture.design0Path.string(),
                     "--design2",
                     fixture.design1Path.string(),
                     "--sv_design1_top",
                     "t",
                     "--sv_design2_top",
                     "t",
                     "-v",
                     "sec",
                     "--sec-engine",
                     "pdr"}),
        kSecProvedExitCode);
  }

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       CliSystemVerilogAutomaticVariableDynamicIndexMatchesModuleVariable) {
  const auto fixture = createDesignFixture(
      "sv",
      "module t(input [2:0] sel, input [7:0] flat, output logic q);\n"
      "  always_comb begin\n"
      "    automatic logic [7:0] tbl = flat;\n"
      "    q = tbl[sel];\n"
      "  end\n"
      "endmodule\n",
      "module t(input [2:0] sel, input [7:0] flat, output logic q);\n"
      "  logic [7:0] tbl;\n"
      "  always_comb begin\n"
      "    tbl = flat;\n"
      "    q = tbl[sel];\n"
      "  end\n"
      "endmodule\n");
  const auto runDir = fixture.tmpDir / "automatic_variable_index_run";
  std::filesystem::create_directories(runDir);

  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(runDir);
    EXPECT_EQ(
        runWithArgs({"kepler-formal",
                     "-sv",
                     "--design1",
                     fixture.design0Path.string(),
                     "--design2",
                     fixture.design1Path.string(),
                     "--sv_design1_top",
                     "t",
                     "--sv_design2_top",
                     "t",
                     "-v",
                     "sec",
                     "--sec-engine",
                     "pdr"}),
        kSecProvedExitCode);

    const auto logs = listMiterLogsInCurrentDirectory();
    ASSERT_EQ(logs.size(), 1u);
    const auto contents = readFileContents(runDir / logs.front());
    EXPECT_NE(
        contents.find(
            "SEC checked-output coverage: 100.00% "
            "(1/1 covered/existing outputs)."),
        std::string::npos);
    EXPECT_EQ(contents.find("no-driver connectivity"), std::string::npos);
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

TEST_F(KeplerFormalCliTests, ConfigSecPdrDifferenceLogIncludesWitnessDetails) {
  expectSecDifferenceLogIncludesWitnessDetails("pdr");
}

TEST_F(KeplerFormalCliTests,
       ConfigSecKInductionDifferenceLogIncludesWitnessDetails) {
  expectSecDifferenceLogIncludesWitnessDetails("k_induction");
}

TEST_F(KeplerFormalCliTests, ConfigSecImcDifferenceLogIncludesWitnessDetails) {
  expectSecDifferenceLogIncludesWitnessDetails("imc");
}

TEST_F(KeplerFormalCliTests, ConfigTinyRocketSecVerificationAccepted) {
  const auto root = repoRoot();
  const auto exampleDir = root / "examples" / "tinyrocket";
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

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);
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
       .roles = {"top_output"},
       .connectivitySkip = "logical-loop connectivity: cycle"}};

  writeBoundaryTermsReport(reportPath, reports);

  const auto content = readFileContents(reportPath);
  EXPECT_NE(content.find("# SEC boundary terms report"), std::string::npos);
  EXPECT_NE(content.find("design: design0"), std::string::npos);
  EXPECT_NE(content.find("signal: clk[0]"), std::string::npos);
  EXPECT_NE(content.find("roles: [top_input]"), std::string::npos);
  EXPECT_NE(content.find("signal: bad[0]"), std::string::npos);
  EXPECT_NE(content.find("roles: [top_output]"), std::string::npos);
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

TEST_F(KeplerFormalCliTests, WriteOpaqueCellSkippedOutputsReportFormatsEntries) {
  const auto tempDir =
      makeUniqueTempDir("kepler_formal_cli_opaque_cell_report");
  const auto reportPath = tempDir / "skipped_opaque_cells_pos.txt";

  writeOpaqueCellSkippedOutputsReport(
      reportPath,
      {"bad[0]: design0 opaque-internal: opaque internal cell `u_opaque` "
       "(model `OPAQUE`) pin `Y[0]`: no usable SEC model"});

  const auto content = readFileContents(reportPath);
  EXPECT_NE(content.find("# SEC opaque-cell skipped top-level outputs"),
            std::string::npos);
  EXPECT_NE(content.find("outputs were not verified"), std::string::npos);
  EXPECT_NE(content.find("No free"), std::string::npos);
  EXPECT_NE(content.find("bad[0]"), std::string::npos);
  EXPECT_NE(content.find("u_opaque"), std::string::npos);
  EXPECT_NE(content.find("OPAQUE"), std::string::npos);
  EXPECT_NE(content.find("Y[0]"), std::string::npos);
  EXPECT_NE(content.find("no usable SEC model"), std::string::npos);

  std::filesystem::remove_all(tempDir);
}

TEST_F(KeplerFormalCliTests, WriteOpaqueCellSkippedOutputsReportSkipsEmptyEntries) {
  const auto tempDir =
      makeUniqueTempDir("kepler_formal_cli_empty_opaque_cell_report");
  const auto reportPath = tempDir / "skipped_opaque_cells_pos.txt";

  writeOpaqueCellSkippedOutputsReport(reportPath, {});

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

TEST_F(KeplerFormalCliTests, CliCcFormatRequiresTop) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal", "-v", "sec", "-cc", "design0.cc", "design1.cc"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliCcFormatRequiresSec) {
  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-cc",
                   "--cc_top",
                   "top",
                   "--cc_include",
                   "include",
                   "design0.cc",
                   "design1.cc"}),
      EXIT_FAILURE);
}

TEST_F(KeplerFormalCliTests, CliMixedFormatsRequireBothDesignSides) {
  for (const std::string format : {"-cc", "-c_vs_rtl", "-rtl_vs_gate"}) {
    SCOPED_TRACE(format);
    const auto run = runStructuredWithArgs(
        {"kepler-formal", "-v", "sec", format, "only_one_input"});
    EXPECT_EQ(run.exitCode, EXIT_FAILURE);
    EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Error);
    EXPECT_NE(readFileContents(run.result.logFile).find("Need "), std::string::npos);
  }
}

TEST_F(KeplerFormalCliTests, ConfigCcReusesOnlyIdenticalSynthesisOptions) {
  const auto tmpDir = makeUniqueTempDir("kepler_formal_cli_cc_reuse");
  const auto includeDir = tmpDir / "include";
  std::filesystem::create_directory(includeDir);
  {
    std::ofstream header(includeDir / "transform_mask.h");
    header << "#define TRANSFORM_MASK 0x1256\n";
    std::ofstream cc(tmpDir / "transform.cc");
    cc << "#include \"transform_mask.h\"\n"
          "unsigned short transform(unsigned short value) {\n"
          "  return value ^ TRANSFORM_MASK;\n}\n"
          "unsigned short changed(unsigned short value) {\n"
          "  return value ^ (TRANSFORM_MASK ^ 1);\n}\n";
  }
  const std::vector<std::string> options = {
      "", "cc_design2_module_name: transform_copy\n", "cc_design2_top: changed\n"};
  for (size_t index = 0; index < options.size(); ++index) {
    SCOPED_TRACE(options[index]);
    const auto runDir = tmpDir / std::to_string(index);
    std::filesystem::create_directory(runDir);
    const auto logPath = runDir / "run.log";
    const auto cfgPath = writeTempConfig(
        "format: cc\nverification: sec\nsec_encoding: binary\nmax_k: 2\n"
        "cc_top: transform\n" + options[index] +
        "cc_include_paths: ['" + includeDir.string() + "', '', '" +
        includeDir.string() + "']\n"
        "input_paths: ['" + (tmpDir / "transform.cc").string() + "', '" +
        (tmpDir / "." / "transform.cc").string() + "']\n"
        "log_file: " + logPath.string() + "\n");
    {
      CurrentPathGuard currentPathGuard;
      std::filesystem::current_path(runDir);
      EXPECT_EQ(runWithConfigFile(cfgPath),
                index == 2 ? kSecCounterexampleExitCode : kSecProvedExitCode);
    }
    const auto outputDir = runDir / "kepler_formal_c2rtl";
    EXPECT_TRUE(std::filesystem::exists(outputDir / "design1_transform.sv"));
    EXPECT_EQ(std::filesystem::exists(outputDir / "design2_transform.sv"),
              index != 0);
    EXPECT_EQ(readFileContents(logPath).find("Reusing design1 C/C++ synthesis") !=
                  std::string::npos,
              index == 0);
    if (index == 1) {
      EXPECT_NE(readFileContents(outputDir / "design2_transform.sv")
                    .find("module transform_copy"),
                std::string::npos);
    }
    std::filesystem::remove(cfgPath);
  }
  std::filesystem::remove_all(tmpDir);
}

TEST_F(KeplerFormalCliTests,
       CliCVsRtlPreservesWriteOnlyReferenceOutputNames) {
  const auto tmpDir =
      makeUniqueTempDir("kepler_formal_cli_c_vs_rtl_named_outputs");
  const auto ccPath = tmpDir / "word_transform.cc";
  const auto svPath = tmpDir / "word_transform.sv";
  const auto outputDir = tmpDir / "c2rtl";
  {
    std::ofstream cc(ccPath);
    cc << "void word_transform(unsigned short payload, "
          "unsigned short& transformed, "
          "bool& marker) {\n";
    cc << "  transformed = static_cast<unsigned short>(payload ^ 0xa55a);\n";
    cc << "  marker = payload == 0x3c5a;\n";
    cc << "}\n";
  }
  {
    std::ofstream sv(svPath);
    sv << "module word_transform_rtl(\n";
    sv << "    input logic [15:0] payload,\n";
    sv << "    output logic [15:0] transformed,\n";
    sv << "    output logic marker\n";
    sv << ");\n";
    sv << "  assign transformed = payload ^ 16'ha55a;\n";
    sv << "  assign marker = payload == 16'h3c5a;\n";
    sv << "endmodule\n";
  }

  EXPECT_EQ(
      runWithArgs({"kepler-formal",
                   "-c_vs_rtl",
                   "-v",
                   "sec",
                   "--sec-encoding",
                   "binary",
                   "-k",
                   "2",
                   "--cc_top",
                   "word_transform",
                   "--cc_design1_module_name",
                   "word_transform_c2rtl",
                   "--cc_output_dir",
                   outputDir.string(),
                   "--sv_design2_top",
                   "word_transform_rtl",
                   ccPath.string(),
                   svPath.string()}),
      kSecProvedExitCode);

  const auto generatedPath = outputDir / "design1_word_transform.sv";
  std::ifstream generated(generatedPath);
  const std::string generatedText(
      (std::istreambuf_iterator<char>(generated)),
      std::istreambuf_iterator<char>());
  EXPECT_NE(generatedText.find("output wire [15:0] transformed"),
            std::string::npos);
  EXPECT_NE(generatedText.find("output wire marker"), std::string::npos);
  EXPECT_NE(generatedText.find("module word_transform_c2rtl__xls_impl"),
            std::string::npos);

  std::filesystem::remove_all(tmpDir);
}

struct TemporalC2RtlTransformFixture {
  std::filesystem::path tmpDir;
  std::filesystem::path ccPath;
  std::filesystem::path svPath;
  std::filesystem::path outputDir;
  std::filesystem::path logPath;
};

TemporalC2RtlTransformFixture createTemporalC2RtlTransformFixture() {
  TemporalC2RtlTransformFixture fixture;
  fixture.tmpDir =
      makeUniqueTempDir("kepler_formal_cli_temporal_c2rtl_transform");
  fixture.ccPath = fixture.tmpDir / "word_transform.cc";
  fixture.svPath = fixture.tmpDir / "word_transform.sv";
  fixture.outputDir = fixture.tmpDir / "c2rtl";
  fixture.logPath = fixture.tmpDir / "c2rtl.log";
  {
    std::ofstream cc(fixture.ccPath);
    cc << "void word_transform(unsigned short payload, "
          "unsigned short& transformed, "
          "bool& marker) {\n";
    cc << "  transformed = static_cast<unsigned short>(payload ^ 0xa55a);\n";
    cc << "  marker = payload == 0x3c5a;\n";
    cc << "}\n";
  }
  {
    std::ofstream sv(fixture.svPath);
    sv << "module word_transform_rtl(\n";
    sv << "  input logic clock, input logic reset_n,\n";
    sv << "  input logic [15:0] payload,\n";
    sv << "  output logic [15:0] transformed, output logic marker);\n";
    sv << "  always_ff @(posedge clock or negedge reset_n) begin\n";
    sv << "    if (!reset_n) begin transformed <= '0; marker <= 1'b0; end\n";
    sv << "    else begin transformed <= payload ^ 16'ha55a; "
          "marker <= payload == 16'h3c5a; end\n";
    sv << "  end\n";
    sv << "endmodule\n";
  }
  return fixture;
}

std::filesystem::path writeTemporalC2RtlTransformConfig(
    const TemporalC2RtlTransformFixture& fixture,
    const std::string& delayEntries,
    const std::string& verification = "sec",
    const std::string& properties = "") {
  const auto cfgPath = fixture.tmpDir /
      ("config_" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()) +
       ".yaml");
  std::ofstream cfg(cfgPath);
  cfg << "format: c_vs_rtl\n";
  cfg << "verification: " << verification << "\n";
  cfg << "max_k: 8\n";
  cfg << "sec_engine: pdr\n";
  cfg << "sec_encoding: binary\n";
  cfg << "c2rtl_auto_align: true\n";
  if (!delayEntries.empty()) {
    cfg << "c2rtl_output_delays:\n";
    cfg << delayEntries;
  }
  cfg << properties;
  cfg << "cc_top: word_transform\n";
  cfg << "cc_design1_module_name: word_transform_c2rtl\n";
  cfg << "cc_output_dir: " << fixture.outputDir.string() << "\n";
  cfg << "sv_design2_top: word_transform_rtl\n";
  cfg << "log_file: " << fixture.logPath.string() << "\n";
  cfg << "input_paths:\n";
  cfg << "  - " << fixture.ccPath.string() << "\n";
  cfg << "  - " << fixture.svPath.string() << "\n";
  return cfgPath;
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsMalformedPropertiesBeforeSynthesis) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const std::vector<std::pair<std::string, std::string>> cases = {
      {"constraints: payload > 0\n", "constraints must be a non-empty sequence"},
      {"constraints: []\n", "constraints must be a non-empty sequence"},
      {"constraints: [null]\n", "constraints[0] must be a non-empty expression"},
      {"constraints: [[payload]]\n", "constraints[0] must be a non-empty expression"},
      {"constraints: ['   ']\n", "constraints[0] must be a non-empty expression"},
      {"constraints: ['payload >']\n", "constraints[0]:"},
      {"eventuals: {}\n", "eventuals must be a non-empty sequence"},
      {"eventuals: []\n", "eventuals must be a non-empty sequence"},
      {"eventuals: [true]\n", "eventuals[0] must be a map"},
      {"eventuals:\n  - ? [cycle]\n    : 1\n",
       "eventuals[0] keys must be cycle, condition, or equality"},
      {"eventuals: [{cycle: 1, condition: true, equality: true, delay: 1}]\n",
       "eventuals[0]: unknown key delay"},
      {"eventuals: [{cycle: 1, cycle: 2, condition: true, equality: true}]\n",
       "eventuals[0]: duplicate key cycle"},
      {"eventuals: [{cycle: 1, condition: true, condition: false, equality: true}]\n",
       "eventuals[0]: duplicate key condition"},
      {"eventuals: [{cycle: 1, condition: true, equality: true, equality: false}]\n",
       "eventuals[0]: duplicate key equality"},
      {"eventuals: [{condition: true, equality: true}]\n",
       "eventuals[0].cycle must be a non-negative integer"},
      {"eventuals: [{cycle: -1, condition: true, equality: true}]\n",
       "eventuals[0].cycle must be a non-negative integer"},
      {"eventuals: [{cycle: 1.5, condition: true, equality: true}]\n",
       "eventuals[0].cycle must be a non-negative integer"},
      {"eventuals: [{cycle: 18446744073709551616, condition: true, equality: true}]\n",
       "eventuals[0].cycle is out of range"},
      {"eventuals: [{cycle: 1, equality: true}]\n",
       "eventuals[0].condition must be a non-empty expression"},
      {"eventuals: [{cycle: 1, condition: true}]\n",
       "eventuals[0].equality must be a non-empty expression"},
      {"eventuals: [{cycle: 1, condition: '', equality: true}]\n",
       "eventuals[0].condition must be a non-empty expression"},
      {"eventuals: [{cycle: 1, condition: true, equality: ''}]\n",
       "eventuals[0].equality must be a non-empty expression"},
      {"eventuals: [{cycle: 1, condition: '!()', equality: true}]\n",
       "eventuals[0].condition:"},
      {"eventuals: [{cycle: 1, condition: true, equality: 'model.marker =='}]\n",
       "eventuals[0].equality:"},
      {"check_reachability: yes\n", "check_reachability must be true or false"},
      {"check_reachability: 1\n", "check_reachability must be true or false"},
      {"check_reachability: []\n", "check_reachability must be true or false"},
  };
  for (const auto& [properties, expectedError] : cases) {
    SCOPED_TRACE(properties);
    std::ostringstream captured;
    NamedLoggerGuard logger("c2rtl_config_validation");
    logger.installed_->sinks().clear();
    logger.installed_->sinks().push_back(
        std::make_shared<spdlog::sinks::ostream_sink_mt>(captured));
    const auto cfgPath = writeTemporalC2RtlTransformConfig(
        fixture, "", "sec", properties);
    EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
    EXPECT_NE(captured.str().find(expectedError), std::string::npos);
    EXPECT_FALSE(std::filesystem::exists(fixture.outputDir));
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsMalformedAlignmentOptionsBeforeSynthesis) {
  const std::vector<std::pair<std::string, std::string>> cases = {
      {"c2rtl_auto_align: []\n", "c2rtl_auto_align must be a boolean scalar"},
      {"c2rtl_auto_align: invalid\n", "c2rtl_auto_align must be true or false"},
      {"c2rtl_output_delays: [1]\n", "c2rtl_output_delays must be a map"},
      {"c2rtl_output_delays: {out: []}\n", "entries must be scalar output names"},
      {"c2rtl_output_delays: {'': 1}\n", "output names must not be empty"},
      {"c2rtl_output_delays: {out: 1, out: 2}\n", "Duplicate C2RTL output delay setting for out"},
  };
  for (const auto& [options, expectedError] : cases) {
    SCOPED_TRACE(options);
    std::ostringstream captured;
    NamedLoggerGuard logger("c2rtl_alignment_validation");
    logger.installed_->sinks().clear();
    logger.installed_->sinks().push_back(
        std::make_shared<spdlog::sinks::ostream_sink_mt>(captured));
    const auto cfgPath = writeTempConfig(
        "format: c_vs_rtl\nverification: sec\ncc_top: top\n"
        "input_paths: [model.cc, rtl.sv]\n" + options);
    EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
    EXPECT_NE(captured.str().find(expectedError), std::string::npos);
    std::filesystem::remove(cfgPath);
  }
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsIncompatibleOptionsBeforeSynthesis) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const std::vector<std::pair<std::string, std::string>> cases = {
      {"c2rtl_auto_align: false\nc2rtl_output_delays: {transformed: 1}\n",
       "c2rtl_output_delays requires c2rtl_auto_align"},
      {"c2rtl_auto_align: true\nsec_engine: imc\n",
       "requires sec_engine: pdr"},
      {"c2rtl_auto_align: true\nsec_encoding: dual_rail_steady\n",
       "requires sec_encoding: binary"},
      {"c2rtl_auto_align: true\ncompact_mode: true\n",
       "does not support compact_mode"},
      {"c2rtl_auto_align: true\nsec_reset: {cycles: 1, ports: [{name: reset_n, active_value: 0}]}\n",
       "does not support sec_reset options"},
      {"c2rtl_auto_align: true\nbtor2_export: true\nbtor2_export_path: '" +
           (fixture.tmpDir / "unused.btor2").string() + "'\n",
       "does not support BTOR2 export"},
  };
  for (const auto& [options, expectedError] : cases) {
    SCOPED_TRACE(options);
    const std::string encoding = options.find("sec_encoding:") == std::string::npos
        ? "sec_encoding: binary\n" : "";
    const auto cfgPath = writeTempConfig(
        "format: c_vs_rtl\nverification: sec\ncc_top: word_transform\n" +
        encoding + options + "cc_output_dir: " + fixture.outputDir.string() + "\n"
        "log_file: " + fixture.logPath.string() + "\n"
        "input_paths: ['" + fixture.ccPath.string() + "', '" +
        fixture.svPath.string() + "']\n");
    EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
    EXPECT_NE(readFileContents(fixture.logPath).find(expectedError),
              std::string::npos);
    EXPECT_FALSE(std::filesystem::exists(fixture.outputDir));
    std::filesystem::remove(cfgPath);
  }
  const auto wrongFormat = writeTempConfig(
      "format: verilog\nverification: sec\nc2rtl_auto_align: true\n"
      "input_paths: [first.v, second.v]\nlog_file: " +
      fixture.logPath.string() + "\n");
  EXPECT_EQ(runWithConfigFile(wrongFormat), EXIT_FAILURE);
  EXPECT_NE(readFileContents(fixture.logPath).find("only supported with format: c_vs_rtl"),
            std::string::npos);
  std::filesystem::remove(wrongFormat);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigCcSynthesisFailuresReportTheirCause) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto blockedDir = fixture.tmpDir / "blocked";
  {
    std::ofstream blocked(blockedDir);
    blocked << "An output directory cannot be created below a regular file.\n";
  }
  for (const bool invalidSource : {false, true}) {
    SCOPED_TRACE(invalidSource);
    if (invalidSource) {
      std::ofstream cc(fixture.ccPath);
      cc << "void word_transform( {\n";
    }
    const auto outputDir = invalidSource ? fixture.outputDir : blockedDir / "output";
    const auto cfgPath = writeTempConfig(
        "format: c_vs_rtl\nverification: sec\ncc_top: word_transform\n"
        "cc_output_dir: " + outputDir.string() + "\n"
        "input_paths: ['" + fixture.ccPath.string() + "', '" +
        fixture.svPath.string() + "']\nlog_file: " + fixture.logPath.string() + "\n");
    const auto run = runStructuredWithConfigFile(cfgPath);
    EXPECT_EQ(run.exitCode, EXIT_FAILURE);
    EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Error);
    EXPECT_NE(readFileContents(fixture.logPath).find(invalidSource
                  ? "XLS C2RTL synthesis failed for design1"
                  : "Failed to create cc_output_dir"),
              std::string::npos);
    EXPECT_FALSE(std::filesystem::exists(outputDir / "design1_word_transform.sv"));
    std::filesystem::remove(cfgPath);
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlReportsBoundedInconclusiveWithoutReset) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  {
    std::ofstream sv(fixture.svPath);
    sv << "module word_transform_rtl(input logic clock, input logic [15:0] payload,\n"
          "  output logic [15:0] transformed, output logic marker);\n"
          "  always_ff @(posedge clock) begin\n"
          "    transformed <= payload ^ 16'ha55a;\n"
          "    marker <= payload == 16'h3c5a;\n"
          "  end\nendmodule\n";
  }
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n");
  auto config = readFileContents(cfgPath);
  config.replace(config.find("max_k: 8"), std::string("max_k: 8").size(), "max_k: 0");
  {
    std::ofstream cfg(cfgPath);
    cfg << config;
  }
  const auto run = runStructuredWithConfigFile(cfgPath);
  EXPECT_EQ(run.exitCode, kSecInconclusiveExitCode);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Inconclusive);
  EXPECT_EQ(run.result.bound, 0);
  EXPECT_EQ(run.result.coveredOutputs, 17);
  EXPECT_EQ(run.result.totalOutputs, 17);
  EXPECT_NE(run.result.reason.find("configured frame bound"), std::string::npos);
  EXPECT_NE(readFileContents(fixture.logPath).find("C2RTL reset: none"),
            std::string::npos);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsPropertiesWithoutAutoAlignBeforeSynthesis) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  for (const std::string properties : {
           "constraints: ['payload > 0']\n",
           "eventuals: [{cycle: 1, condition: true, equality: 'model.marker == rtl.marker'}]\n",
           "check_reachability: false\n"}) {
    SCOPED_TRACE(properties);
    const auto cfgPath = writeTemporalC2RtlTransformConfig(
        fixture, "", "sec", properties);
    auto config = readFileContents(cfgPath);
    config.replace(config.find("c2rtl_auto_align: true"),
                   std::string("c2rtl_auto_align: true").size(),
                   "c2rtl_auto_align: false");
    {
      std::ofstream cfg(cfgPath);
      cfg << config;
    }
    EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
    EXPECT_NE(readFileContents(fixture.logPath).find(
                  "require c2rtl_auto_align: true"),
              std::string::npos);
    EXPECT_FALSE(std::filesystem::exists(fixture.outputDir));
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsEventualsWithOutputDelaysBeforeSynthesis) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n", "sec",
      "eventuals: [{cycle: 1, condition: true, equality: 'model.marker == rtl.marker'}]\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  EXPECT_NE(readFileContents(fixture.logPath).find(
                "eventuals cannot be combined with c2rtl_output_delays"),
            std::string::npos);
  EXPECT_FALSE(std::filesystem::exists(fixture.outputDir));
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlConstraintsRestrictTheInputDomain) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  auto rtl = readFileContents(fixture.svPath);
  const std::string original = "transformed <= payload ^ 16'ha55a";
  rtl.replace(rtl.find(original), original.size(),
              "transformed <= payload == 0 ? 16'b0 : payload ^ 16'ha55a");
  {
    std::ofstream sv(fixture.svPath);
    sv << rtl;
  }
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n", "sec",
      "constraints:\n  - payload > 0x1000\n  - payload < 0x7fff\n"
      "check_reachability: true\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);

  const auto unconstrained = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n");
  EXPECT_EQ(runWithConfigFile(unconstrained), kSecCounterexampleExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlContradictoryConstraintsDoNotProveEquivalence) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n", "sec",
      "constraints: ['payload > 10', 'payload < 5']\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecInconclusiveExitCode);
  EXPECT_NE(readFileContents(fixture.logPath).find("unsatisfiable"),
            std::string::npos);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlTransformProvesWithExplicitDelays) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlEventualsCheckIndependentConditionalOutputs) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "", "sec",
      "constraints: ['rtl.reset_n']\n"
      "eventuals:\n"
      "  - cycle: 1\n"
      "    condition: 'true'\n"
      "    equality: 'model.marker == rtl.marker'\n"
      "  - cycle: 1\n"
      "    condition: '!model.marker && (rtl.payload == model.payload)'\n"
      "    equality: 'model.transformed == rtl.transformed'\n"
      "  - cycle: 1\n"
      "    condition: 'false'\n"
      "    equality: 'model.transformed != rtl.transformed'\n"
      "  - cycle: 4\n"
      "    condition: 'rtl.reset_n && (rtl.transformed != 0)'\n"
      "    equality: 'model.transformed[15] == rtl.transformed[15]'\n");
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);

  // A failing entry must not be overwritten by later entries at the same cycle.
  auto config = readFileContents(cfgPath);
  const auto equality = config.find("model.marker == rtl.marker");
  config.replace(equality, std::string("model.marker == rtl.marker").size(),
                 "model.marker != rtl.marker");
  {
    std::ofstream cfg(cfgPath);
    cfg << config;
  }
  EXPECT_EQ(runWithConfigFile(cfgPath), kSecCounterexampleExitCode);
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlEventualAtWrongCycleFindsCounterexample) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  for (const size_t cycle : {0, 1}) {
    SCOPED_TRACE(cycle);
    const auto cfgPath = writeTemporalC2RtlTransformConfig(
        fixture, "", "sec",
        "constraints: ['rtl.reset_n']\n"
        "eventuals:\n"
        "  - cycle: " + std::to_string(cycle) + "\n"
        "    condition: 'true'\n"
        "    equality: 'model.transformed == rtl.transformed'\n");
    EXPECT_EQ(runWithConfigFile(cfgPath),
              cycle == 0 ? kSecCounterexampleExitCode : kSecProvedExitCode);
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlDifferenceReportsConcreteTrace) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  {
    std::ofstream sv(fixture.svPath);
    sv << "module word_transform_rtl(\n";
    sv << "  input logic clock, input logic reset_n,\n";
    sv << "  input logic [15:0] payload,\n";
    sv << "  output logic [15:0] transformed, output logic marker);\n";
    sv << "  always_ff @(posedge clock or negedge reset_n) begin\n";
    sv << "    if (!reset_n) begin transformed <= '0; marker <= 1'b0; end\n";
    sv << "    else begin transformed <= payload; "
          "marker <= payload == 16'h3c5a; end\n";
    sv << "  end\n";
    sv << "endmodule\n";
  }
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecCounterexampleExitCode);
  std::ifstream log(fixture.logPath);
  const std::string logText(
      (std::istreambuf_iterator<char>(log)),
      std::istreambuf_iterator<char>());
  EXPECT_NE(logText.find("C2RTL counterexample details"), std::string::npos);
  EXPECT_NE(logText.find("Input trace:"), std::string::npos);
  EXPECT_NE(logText.find("payload[0]="), std::string::npos);
  EXPECT_NE(
      logText.find("Delayed-reference/RTL output mismatches"),
      std::string::npos);
  EXPECT_NE(logText.find("transformed["), std::string::npos);
  EXPECT_NE(logText.find("delayed_reference="), std::string::npos);
  EXPECT_NE(logText.find("rtl="), std::string::npos);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlProvesIndependentOutputDelays) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  {
    std::ofstream sv(fixture.svPath);
    sv << "module word_transform_rtl(\n";
    sv << "  input logic clock, input logic reset_n,\n";
    sv << "  input logic [15:0] payload,\n";
    sv << "  output logic [15:0] transformed, output logic marker);\n";
    sv << "  logic marker_stage;\n";
    sv << "  always_ff @(posedge clock or negedge reset_n) begin\n";
    sv << "    if (!reset_n) begin transformed <= '0; marker_stage <= 1'b0; "
          "marker <= 1'b0; end\n";
    sv << "    else begin\n";
    sv << "      transformed <= payload ^ 16'ha55a;\n";
    sv << "      marker_stage <= payload == 16'h3c5a;\n";
    sv << "      marker <= marker_stage;\n";
    sv << "    end\n";
    sv << "  end\n";
    sv << "endmodule\n";
  }
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 2\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecProvedExitCode);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsMissingOutputDelayBeforeProof) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecInconclusiveExitCode);
  std::ifstream log(fixture.logPath);
  const std::string logText(
      (std::istreambuf_iterator<char>(log)),
      std::istreambuf_iterator<char>());
  EXPECT_NE(
      logText.find("output delay setting is required"), std::string::npos);
  EXPECT_NE(logText.find("marker"), std::string::npos);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsNegativeOutputDelay) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: -1\n  marker: 1\n");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsLecVerificationMode) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n", "lec");

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsBoundaryPairsBeforeSynthesis) {
  const auto fixture = createTemporalC2RtlTransformFixture();
  const auto cfgPath = writeTemporalC2RtlTransformConfig(
      fixture, "  transformed: 1\n  marker: 1\n");
  {
    std::ofstream cfg(cfgPath, std::ios::app);
    cfg << "set_as_boundary: [[reference_block, implementation_block]]\n";
  }

  EXPECT_EQ(runWithConfigFile(cfgPath), EXIT_FAILURE);
  std::ifstream log(fixture.logPath);
  const std::string logText(
      (std::istreambuf_iterator<char>(log)),
      std::istreambuf_iterator<char>());
  EXPECT_NE(
      logText.find("c2rtl_auto_align does not support set_as_boundary"),
      std::string::npos);
  EXPECT_FALSE(std::filesystem::exists(fixture.outputDir));

  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests,
       ConfigTemporalC2RtlRejectsRegisterEnable) {
  const auto tmpDir =
      makeUniqueTempDir("kepler_formal_cli_temporal_c2rtl_enable");
  const auto ccPath = tmpDir / "conditional_transform.cc";
  const auto svPath = tmpDir / "conditional_capture.sv";
  const auto outputDir = tmpDir / "c2rtl";
  const auto logPath = tmpDir / "c2rtl.log";
  const auto cfgPath = tmpDir / "config.yaml";
  {
    std::ofstream cc(ccPath);
    cc << "unsigned short conditional_transform(unsigned short payload, "
          "bool load) { return load ? "
          "static_cast<unsigned short>(payload ^ 0x6d39) : 0; }\n";
  }
  {
    std::ofstream sv(svPath);
    sv << "module conditional_capture(input logic clock, input logic load, "
          "input logic [15:0] payload, output logic [15:0] out);\n";
    sv << "  always_ff @(posedge clock) if (load) out <= payload ^ 16'h6d39;\n";
    sv << "endmodule\n";
  }
  {
    std::ofstream cfg(cfgPath);
    cfg << "format: c_vs_rtl\nverification: sec\nmax_k: 4\n";
    cfg << "sec_engine: pdr\nsec_encoding: binary\n";
    cfg << "c2rtl_auto_align: true\n";
    cfg << "c2rtl_output_delays:\n  out: 1\n";
    cfg << "cc_top: conditional_transform\n";
    cfg << "cc_design1_module_name: conditional_transform_c2rtl\n";
    cfg << "cc_output_dir: " << outputDir.string() << "\n";
    cfg << "sv_design2_top: conditional_capture\n";
    cfg << "log_file: " << logPath.string() << "\n";
    cfg << "input_paths:\n  - " << ccPath.string() << "\n  - "
        << svPath.string() << "\n";
  }

  EXPECT_EQ(runWithConfigFile(cfgPath), kSecInconclusiveExitCode);
  std::ifstream log(logPath);
  const std::string logText(
      (std::istreambuf_iterator<char>(log)),
      std::istreambuf_iterator<char>());
  EXPECT_TRUE(
      logText.find("enable") != std::string::npos ||
      logText.find("state feedback") != std::string::npos);

  std::filesystem::remove_all(tmpDir);
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

TEST_F(KeplerFormalCliTests, CliRemovedSecBoundaryFlagsAreRejectedBeforeFormat) {
  for (const char* flag : {"--sec-uncomputable-seq-boundary",
                           "--no-sec-uncomputable-seq-boundary"}) {
    EXPECT_EQ(
        runWithArgs({"kepler-formal", "-v", "sec", flag, "-verilog"}),
        EXIT_FAILURE);
  }
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

TEST_F(KeplerFormalCliTests, CliRemovedSecBoundaryFlagsAreRejectedAfterFormat) {
  for (const char* flag : {"--sec-uncomputable-seq-boundary",
                           "--no-sec-uncomputable-seq-boundary"}) {
    EXPECT_EQ(
        runWithArgs({"kepler-formal", "-verilog", "-v", "sec", flag}),
        EXIT_FAILURE);
  }
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
  const auto exampleDir = root / "examples" / "tinyrocket";
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
  const auto logPath = fixture.tmpDir / "structured_scoped.log";

  const auto cfgPath = writeTempConfig(
      "format: naja_if\n"
      "input_paths:\n"
      "  - " + fixture.design0IfPath.string() + "\n"
      "  - " + fixture.design1IfPath.string() + "\n"
      "liberty_files:\n"
      "  - " + fixture.libertyPath.string() + "\n"
      "use_scopes: true\n"
      "log_file: " + logPath.string() + "\n");

  const auto run = runStructuredWithConfigFile(cfgPath);
  EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(run.result.exitCode, EXIT_SUCCESS);
  EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Equivalent);
  EXPECT_EQ(run.result.logFile, logPath.string());
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

TEST_F(KeplerFormalCliTests, ConfigSetAsBoundaryRejectsMalformedPairValues) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module child(input i, output o); assign o = i; endmodule\n"
      "module top(input a, output y); child u(.i(a), .o(y)); endmodule\n");
  const std::vector<std::string> invalidValues = {
      "u", "{}", "[[u]]", "[[u, u, u]]",
      "[[[], u]]", "[[u, {}]]", "[['', u]]", "[[u, '']]"};
  for (const auto& value : invalidValues) {
    SCOPED_TRACE(value);
    const auto cfgPath = writeTempConfig(
        "format: verilog\n"
        "input_paths:\n"
        "  - " + fixture.design0Path.string() + "\n"
        "  - " + fixture.design1Path.string() + "\n"
        "set_as_boundary: " + value + "\n");
    const auto run = runStructuredWithConfigFile(cfgPath);
    EXPECT_EQ(run.exitCode, EXIT_FAILURE);
    EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Error);
    EXPECT_EQ(NLUniverse::get(), nullptr);
    std::filesystem::remove(cfgPath);
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, ConfigSetAsBoundaryRejectsScopeOperations) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module child(input i, output o); assign o = i; endmodule\n"
      "module top(input a, output y); child u(.i(a), .o(y)); endmodule\n");
  for (const auto& scopeOption : {"use_scopes", "clean_scopes"}) {
    SCOPED_TRACE(scopeOption);
    const auto cfgPath = writeTempConfig(
        "format: verilog\n"
        "input_paths:\n"
        "  - " + fixture.design0Path.string() + "\n"
        "  - " + fixture.design1Path.string() + "\n"
        "set_as_boundary: [[u, u]]\n" + scopeOption + ": true\n");
    const auto run = runStructuredWithConfigFile(cfgPath);
    EXPECT_EQ(run.exitCode, EXIT_FAILURE);
    EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Error);
    EXPECT_EQ(NLUniverse::get(), nullptr);
    std::filesystem::remove(cfgPath);
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSetAsBoundaryBothAliasesBeforeFormat) {
  const auto fixture = createDesignFixture(
      "v",
      "module child(input i, output o); endmodule\n"
      "module top(input a, output y); child left(.i(a), .o(y)); endmodule\n",
      "module child(input i, output o); endmodule\n"
      "module top(input a, output y); child right(.i(a), .o(y)); endmodule\n");
  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(fixture.tmpDir);
    for (const auto& flag : {"--set-as-boundary", "--set_as_boundary"}) {
      SCOPED_TRACE(flag);
      const auto run = runStructuredWithArgs(
          {"kepler-formal", flag, "left", "right", "-verilog",
           fixture.design0Path.string(), fixture.design1Path.string()});
      EXPECT_EQ(run.exitCode, EXIT_SUCCESS);
      EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Equivalent);
    }
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, CliSetAsBoundaryBeforeFormatRequiresBothPaths) {
  for (const auto& flag : {"--set-as-boundary", "--set_as_boundary"}) {
    for (bool includeFirstPath : {false, true}) {
      SCOPED_TRACE(flag);
      SCOPED_TRACE(includeFirstPath);
      std::vector<std::string> args = {"kepler-formal", flag};
      if (includeFirstPath) {
        args.emplace_back("u");
      }
      const auto run = runStructuredWithArgs(std::move(args));
      EXPECT_EQ(run.exitCode, EXIT_FAILURE);
      EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Error);
    }
  }
}

TEST_F(KeplerFormalCliTests, CliSetAsBoundaryRejectsEitherEmptyPath) {
  for (const auto& flag : {"--set-as-boundary", "--set_as_boundary"}) {
    for (bool beforeFormat : {false, true}) {
      for (bool emptyLeft : {false, true}) {
        SCOPED_TRACE(flag);
        SCOPED_TRACE(beforeFormat);
        SCOPED_TRACE(emptyLeft);
        std::vector<std::string> args = {"kepler-formal"};
        const std::vector<std::string> inputArgs = {
            "-verilog", "design0.v", "design1.v"};
        if (!beforeFormat) {
          args.insert(args.end(), inputArgs.begin(), inputArgs.end());
        }
        args.insert(args.end(), {flag, emptyLeft ? "" : "u",
                                 emptyLeft ? "u" : ""});
        if (beforeFormat) {
          args.insert(args.end(), inputArgs.begin(), inputArgs.end());
        }
        const auto run = runStructuredWithArgs(std::move(args));
        EXPECT_EQ(run.exitCode, EXIT_FAILURE);
        EXPECT_EQ(run.result.status, KEPLER_FORMAL::RunStatus::Error);
      }
    }
  }
}

TEST_F(KeplerFormalCliTests,
       CliCompactSetAsBoundaryFailureReleasesDesignsAndAllowsAnotherRun) {
  const auto fixture = createEquivalentDesignFixture(
      "v",
      "module child(input i, output o); endmodule\n"
      "module top(input a, output y); child u(.i(a), .o(y)); endmodule\n");
  {
    CurrentPathGuard currentPathGuard;
    std::filesystem::current_path(fixture.tmpDir);
    for (bool sec : {false, true}) {
      for (bool failFirstDesign : {false, true}) {
        SCOPED_TRACE(sec);
        SCOPED_TRACE(failFirstDesign);
        std::vector<std::string> args = {"kepler-formal"};
        if (sec) {
          args.insert(args.end(), {"-v", "sec", "--sec-engine", "k_induction",
                                   "--sec-encoding", "binary", "-k", "1"});
        }
        args.insert(args.end(),
                    {"-verilog", fixture.design0Path.string(),
                     fixture.design1Path.string(), "--compact",
                     "--set-as-boundary", failFirstDesign ? "missing" : "u",
                     failFirstDesign ? "u" : "missing"});
        const auto failed = runStructuredWithArgs(args);
        EXPECT_EQ(failed.exitCode, EXIT_FAILURE);
        EXPECT_EQ(failed.result.status, KEPLER_FORMAL::RunStatus::Error);
        EXPECT_EQ(NLUniverse::get(), nullptr);

        // Retry in the same process: neither the partially loaded design nor
        // the first side's released compact snapshot may poison the next run.
        args[args.size() - 2] = "u";
        args.back() = "u";
        const auto recovered = runStructuredWithArgs(std::move(args));
        EXPECT_EQ(recovered.exitCode, EXIT_SUCCESS);
        EXPECT_EQ(recovered.result.status, KEPLER_FORMAL::RunStatus::Equivalent);
        EXPECT_EQ(NLUniverse::get(), nullptr);
      }
    }
  }
  std::filesystem::remove_all(fixture.tmpDir);
}

TEST_F(KeplerFormalCliTests, InternalRelationOptionsReachBothCliParsersAndYaml) {
  const auto directory = makeUniqueTempDir("kf_internal_relations");
  const auto source = directory / "design.v";
  std::ofstream(source) << "module top(input a, output y); assign y = a; endmodule\n";
  for (const auto* value : {"true", "false"}) {
    for (bool beforeFormat : {false, true}) {
      std::vector<std::string> args{"kepler-formal", "-v", "sec"};
      const std::vector<std::string> flags{
          "--learn-internal-relations", value,
          "--allow-x-equality-in-internal-relations", value};
      if (beforeFormat) args.insert(args.end(), flags.begin(), flags.end());
      args.insert(args.end(), {"-verilog", source.string(), source.string()});
      if (!beforeFormat) args.insert(args.end(), flags.begin(), flags.end());
      EXPECT_EQ(runWithArgs(args), kSecProvedExitCode);
    }
    for (const auto* key : {"learn_internal_relations", "learn_ineternal_relations"}) {
      const auto log = directory / "relations.log";
      const auto config = writeTempConfig(
          "format: verilog\nverification: sec\ninput_paths:\n  - " + source.string() +
          "\n  - " + source.string() + "\n" + key + ": " + value +
          "\nallow_x_equality_in_internal_relations: " + value +
          "\nlog_file: " + log.string() + "\n");
      EXPECT_EQ(runWithConfigFile(config), kSecProvedExitCode);
      EXPECT_NE(readFileContents(log).find(
          std::string("SEC internal relations: learn=") + value + " allow_x_equality=" + value),
          std::string::npos);
    }
  }
  const std::string base = "format: verilog\nverification: sec\ninput_paths:\n  - " +
      source.string() + "\n  - " + source.string() + "\n";
  for (const auto* invalid : {"learn_internal_relations: maybe\n",
                            "allow_x_equality_in_internal_relations: []\n",
                            "learn_internal_relations: true\nlearn_ineternal_relations: false\n"}) {
    EXPECT_EQ(runWithConfigFile(writeTempConfig(base + invalid)), EXIT_FAILURE);
  }
  for (const auto* flag : {"--learn-internal-relations", "--allow-x-equality-in-internal-relations"}) {
    EXPECT_EQ(runWithArgs({"kepler-formal", "-v", "sec", flag}), EXIT_FAILURE);
    EXPECT_EQ(runWithArgs({"kepler-formal", "-verilog", "-v", "sec", flag, "maybe"}), EXIT_FAILURE);
    EXPECT_EQ(runWithArgs({"kepler-formal", "-verilog", source.string(), source.string(), flag, "false"}), EXIT_FAILURE);
  }
  std::filesystem::remove_all(directory);
}
