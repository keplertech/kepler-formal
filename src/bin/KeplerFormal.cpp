// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <optional>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <yaml-cpp/yaml.h>

#include "NajaPerf.h"

// Naja interfaces
#include "DNL.h"
#include "MiterStrategy.h"
#include "SNLCapnP.h"
#include "SNLLibertyConstructor.h"
#include "SNLVRLConstructor.h"
#include "SNLVRLDumper.h"
#include "SNLUtils.h"
#include "ScopeExtraction.h"

static void print_usage(const char* prog) {
  std::printf(
      "Usage: %s [--config <file>] | <-naja_if/-verilog> <netlist1> <netlist2> "
      "[<liberty-file>...]\n",
      prog);
}

static std::vector<std::string> yamlToVector(const YAML::Node& node) {
  std::vector<std::string> out;
  if (!node) return out;
  if (!node.IsSequence()) return out;
  for (const auto& n : node) {
    if (n.IsScalar()) out.emplace_back(n.as<std::string>());
  }
  return out;
}

int main(int argc, char** argv) {
  using namespace std::chrono;
  enum class FormatType { VERILOG, NAJA_IF };

  // Default values
  FormatType inputFormatType = FormatType::VERILOG;
  std::vector<std::string> inputPaths;
  std::vector<std::string> libertyFiles;
  std::string logLevel = "info";

  // Basic argument sanity
  if (argc < 2) {
    print_usage(argv[0]);
    return EXIT_SUCCESS;
  }

  // Check for config mode (--config or -c). If present, YAML takes precedence.
  bool usedConfig = false;

  std::string logFileName;

  bool useScopes = false;
  bool cleanScopes = false;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--config" || a == "-c") {
      if (i + 1 >= argc) {
        SPDLOG_CRITICAL("Missing config file after {}", a);
        return EXIT_FAILURE;
      }
      const std::string cfgPath = argv[i + 1];
      try {
        YAML::Node cfg = YAML::LoadFile(cfgPath);

        // format
        if (cfg["format"] && cfg["format"].IsScalar()) {
          std::string fmt = cfg["format"].as<std::string>();
          if (fmt == "naja_if" || fmt == "naja-if" || fmt == "snl")
            inputFormatType = FormatType::NAJA_IF;
          else if (fmt == "verilog" || fmt == "v")
            inputFormatType = FormatType::VERILOG;
          else {
            SPDLOG_CRITICAL("Unrecognized format in config: {}", fmt);
            return EXIT_FAILURE;
          }
        }

        // input_paths
        inputPaths = yamlToVector(cfg["input_paths"]);

        // liberty_files
        libertyFiles = yamlToVector(cfg["liberty_files"]);

        // log level
        if (cfg["log_level"] && cfg["log_level"].IsScalar()) {
          logLevel = cfg["log_level"].as<std::string>();
        }

        // Add log file name
        if (cfg["log_file"] && cfg["log_file"].IsScalar()) {
          logFileName = cfg["log_file"].as<std::string>();
        }
        
        // use_scopes
        if (cfg["use_scopes"] && cfg["use_scopes"].IsScalar()) {
          useScopes = cfg["use_scopes"].as<bool>();
        }

        // clean_scopes
        if (cfg["clean_scopes"] && cfg["clean_scopes"].IsScalar()) {
          cleanScopes = cfg["clean_scopes"].as<bool>();
        }

        usedConfig = true;
      } catch (const std::exception& e) {
        SPDLOG_CRITICAL("Failed to parse config {}: {}", cfgPath, e.what());
        return EXIT_FAILURE;
      }
      break;
    }
  }

  // If not using config, fall back to original CLI parsing
  if (!usedConfig) {
    if (argc < 4 || (std::string(argv[1]) == "--help") ||
        (std::string(argv[1]) == "-h")) {
      print_usage(argv[0]);
      return EXIT_SUCCESS;
    }

    std::string formatType = argv[1];
    if (formatType == "-naja_if" || formatType == "-naja-if") {
      inputFormatType = FormatType::NAJA_IF;
    } else if (formatType == "-verilog") {
      inputFormatType = FormatType::VERILOG;
    } else {
      SPDLOG_CRITICAL("Unrecognized input format type: {}", formatType);
      return EXIT_FAILURE;
    }

    // collect paths and liberty files from argv
    for (int i = 2; i < argc; ++i) inputPaths.emplace_back(argv[i]);

    // If user provided more than two paths, treat the rest as liberty files
    if (inputPaths.size() > 2) {
      for (size_t i = 2; i < inputPaths.size(); ++i)
        libertyFiles.push_back(inputPaths[i]);
    }
  }

  // Basic validation
  if (inputPaths.size() < 2) {
    SPDLOG_CRITICAL("Need two input netlist paths; got {}", inputPaths.size());
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  // Configure logging level
  auto console = spdlog::stdout_color_mt("console");
  if (logLevel == "debug")
    spdlog::set_level(spdlog::level::debug);
  else if (logLevel == "info")
    spdlog::set_level(spdlog::level::info);
  // else if (logLevel == "warn")
  //   spdlog::set_level(spdlog::level::warn);
  // else if (logLevel == "error")
  //   spdlog::set_level(spdlog::level::err);
  // else if (logLevel == "critical")
  //   spdlog::set_level(spdlog::level::critical);
  else
    spdlog::set_level(spdlog::level::info);

  std::printf("KEPLER FORMAL: Run.\n");
  std::printf("Input format: %s\n", (inputFormatType == FormatType::NAJA_IF) ? "SNL" : "VERILOG");
  std::printf("Netlist 1: %s\n", inputPaths[0].c_str());
  std::printf("Netlist 2: %s\n", inputPaths[1].c_str());
  if (!libertyFiles.empty()) {
    for (const auto& lf : libertyFiles) std::printf("Liberty: %s\n", lf.c_str());
  }

  // --------------------------------------------------------------------------
  // 2. Load two netlists via Cap’n Proto (or via VRL constructor)
  // --------------------------------------------------------------------------
  NLUniverse::create();
  NLDB* db0 = nullptr;
  bool primitivesAreLoaded = false;

  if (!libertyFiles.empty()) {
    db0 = NLDB::create(NLUniverse::get());
    auto primitivesLibrary =
        NLLibrary::create(db0, NLLibrary::Type::Primitives, NLName("PRIMS"));
    SNLLibertyConstructor constructor(primitivesLibrary);
    for (const auto& lf : libertyFiles) {
      std::printf("Loading liberty file: %s\n", lf.c_str());
      constructor.construct(lf.c_str());
    }
    primitivesAreLoaded = true;
  }

  if (inputFormatType == FormatType::VERILOG) {
    printf("Parsing verilog file: %s\n", inputPaths[0].c_str());
    auto designLibrary = NLLibrary::create(db0, NLName("DESIGN"));
    SNLVRLConstructor constructor(designLibrary);
    constructor.construct(inputPaths[0].c_str());
    auto top = SNLUtils::findTop(designLibrary);
    if (top) {
      db0->setTopDesign(top);
      SPDLOG_INFO("Found top design: {}", top->getString());
    } else {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("No top design was found after parsing verilog");
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  } else {  // SNL
    std::printf("Loading Naja IF: %s\n", inputPaths[0].c_str());
    naja::NL::SNLCapnP::LoadingConfiguration config;
    config.primitiveConflictPolicy_ = primitivesAreLoaded ? naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::PreferExisting :
                                                            naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::ForbidConflicts;
    db0 = SNLCapnP::load(inputPaths[0].c_str(), config);
    if (!db0) {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("Failed to load Naja IF: {}", inputPaths[0]);
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  }

  // get db0 top
  auto top0 = db0->getTopDesign();
  if (!top0) {
    // LCOV_EXCL_START
    SPDLOG_CRITICAL("Top design not set for first netlist");
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
  }
  db0->setID(2);  // Increment ID to avoid conflicts

  NLDB* db1 = nullptr;

  // Prepare second DB and primitives if needed
  if (!libertyFiles.empty()) {
    db1 = NLDB::create(NLUniverse::get());
    db1->setID(1);
    auto primitivesLibrary =
        NLLibrary::create(db1, NLLibrary::Type::Primitives, NLName("PRIMS"));
    SNLLibertyConstructor constructor(primitivesLibrary);
    for (const auto& lf : libertyFiles) {
      constructor.construct(lf.c_str());
    }
  }

  if (inputFormatType == FormatType::VERILOG) {
    printf("Parsing verilog file: %s\n", inputPaths[1].c_str());
    auto designLibrary = NLLibrary::create(db1, NLName("DESIGN"));
    SNLVRLConstructor constructor(designLibrary);
    constructor.construct(inputPaths[1].c_str());
    auto top = SNLUtils::findTop(designLibrary);
    if (top) {
      db1->setTopDesign(top);
      SPDLOG_INFO("Found top design: {}", top->getString());
    } else {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("No top design was found after parsing verilog");
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  } else {  // SNL
    std::printf("Loading Naja IF: %s\n", inputPaths[1].c_str());
    naja::NL::SNLCapnP::LoadingConfiguration config;
    config.primitiveConflictPolicy_ = primitivesAreLoaded ? naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::PreferExisting :
                                                            naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::ForbidConflicts;
    db1 = SNLCapnP::load(inputPaths[1].c_str(), config);
    if (!db1) {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("Failed to load Naja IF: {}", inputPaths[1]);
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  }

  // get db1 top
  auto top1 = db1->getTopDesign();
  if (!top1) {
    // LCOV_EXCL_START
    SPDLOG_CRITICAL("Top design not set for second netlist");
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
  }

  // --------------------------------------------------------------------------
  // 4. Hand off to the rest of the editing/analysis workflow
  // --------------------------------------------------------------------------
  if (inputFormatType == FormatType::NAJA_IF && useScopes) {
    KEPLER_FORMAL::MiterStrategy MiterS(top0, top1);
    MiterS.init();
    ScopeExtraction extractor(top0, top1);
    extractor.collectVerificationScopes();
    if (cleanScopes) {
      extractor.cleanVerificationScopes(MiterS.getPIs0(), MiterS.getPIs1());
    }
    for (auto scopes : extractor.getScopesToVerify()) {
      SPDLOG_INFO("Looking at scope: {} ",
                  scopes.first->getName().getString());
      // std::string scopeLogFile = (logFileName.empty() ? "kf_" : logFileName) + "_" +
      //                           scopes.first->getName().getString() + ".txt";
      try {
        KEPLER_FORMAL::MiterStrategy MiterScope(scopes.first, scopes.second, logFileName);
        MiterScope.init();
        if (MiterScope.run()) {
          SPDLOG_INFO("No difference was found for scope: {} , {}",
                      scopes.first->getName().getString(),
                      scopes.second->getName().getString());
        } else {
          SPDLOG_INFO("Difference was found for scope: {} , {}. Please refer to the log(miter_log_x.txt) for details.",
                      scopes.first->getName().getString(),
                      scopes.second->getName().getString());
        }
      } catch (const std::exception& e) {
        // LCOV_EXCL_START
        SPDLOG_ERROR("Workflow failed for scope: {} , {}: {}", 
                      scopes.first->getName().getString(),
                      scopes.second->getName().getString(),
                      e.what());
        return EXIT_FAILURE;
        // LCOV_EXCL_STOP
      }
    }
  } else {
    try {
      KEPLER_FORMAL::MiterStrategy MiterS(top0, top1, logFileName);
      MiterS.init();
      if (MiterS.run()) {
        SPDLOG_INFO("No difference was found.");
      } else {
        SPDLOG_INFO("Difference was found. Please refer to the log(miter_log_x.txt) for details.");
      }
    } catch (const std::exception& e) {
      // LCOV_EXCL_START
      SPDLOG_ERROR("Workflow failed: {}", e.what());
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  }

  return EXIT_SUCCESS;
}