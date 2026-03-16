// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <optional>
#include <cctype>
#include <unordered_set>
#include <filesystem>
#include <sstream>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <yaml-cpp/yaml.h>

#include "NajaPerf.h"

// Naja interfaces
#include "DNL.h"
#include "MiterStrategy.h"
#include "SNLCapnP.h"
#include "SNLLibertyConstructor.h"
#include "SNLPyLoader.h"
#include "SNLSVConstructor.h"
#include "SNLVRLConstructor.h"
#include "SNLVRLDumper.h"
#include "SNLUtils.h"
#include "ScopeExtraction.h"
#include "Config.h"

static void print_usage(const char* prog) {
  SPDLOG_INFO(
      "Usage: {} [--config <file>] | <-naja_if/-verilog/-systemverilog/-sv> "
      "<netlist1> <netlist2> [<library-file>...] | "
      "<-naja_if/-verilog/-systemverilog/-sv> --design1 <file...> --design2 "
      "<file...> [--liberty <library-file>...]",
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

static bool validateConfigKeys(const YAML::Node& cfg) {
  if (!cfg || !cfg.IsMap()) {
    return true;
  }
  static const std::unordered_set<std::string> kAllowedKeys = {
      "format",
      "input_paths",
      "liberty_files",
      "verilog_preprocessing",
      "log_level",
      "log_file",
      "use_scopes",
      "clean_scopes",
      "cnf_export",
      "cnf_export_path",
      "dump_cnf",
      "dump_cnf_path",
      "solver",
  };

  for (auto it = cfg.begin(); it != cfg.end(); ++it) {
    if (!it->first.IsScalar()) {
      SPDLOG_CRITICAL("Config key is not a scalar; invalid YAML key");
      return false;
    }
    const std::string key = it->first.as<std::string>();
    if (kAllowedKeys.find(key) == kAllowedKeys.end()) {
      SPDLOG_CRITICAL("Unknown config option: {}", key);
      return false;
    }
  }
  return true;
}

std::string sanitizeFileToken(const std::string& input) {
  std::string out;
  out.reserve(input.size());
  for (unsigned char ch : input) {
    if (std::isalnum(ch) || ch == '_' || ch == '-' || ch == '.') {
      out.push_back(static_cast<char>(ch));
    } else {
      out.push_back('_');
    }
  }
  if (out.empty()) {
    out = "scope";
  }
  return out;
}

struct DesignInputs {
  std::vector<std::string> design0;
  std::vector<std::string> design1;
};

static bool parseConfigInputPaths(const YAML::Node& node,
                                  DesignInputs& out,
                                  std::string& error) {
  out.design0.clear();
  out.design1.clear();
  if (!node) {
    error = "Missing input_paths in config";
    return false;
  }
  if (!node.IsSequence()) {
    error = "input_paths must be a sequence";
    return false;
  }
  if (node.size() == 0) {
    error = "input_paths must contain at least two entries";
    return false;
  }

  const bool firstIsSequence = node[0].IsSequence();
  if (firstIsSequence) {
    if (node.size() != 2) {
      error = "input_paths must contain exactly two design entries";
      return false;
    }
    for (size_t i = 0; i < 2; ++i) {
      if (!node[i].IsSequence()) {
        error = "input_paths must be a sequence of sequences when using multi-file designs";
        return false;
      }
      for (const auto& n : node[i]) {
        if (!n.IsScalar()) {
          error = "input_paths entries must be scalar file paths";
          return false;
        }
        if (i == 0) out.design0.emplace_back(n.as<std::string>());
        else out.design1.emplace_back(n.as<std::string>());
      }
    }
  } else {
    std::vector<std::string> flat = yamlToVector(node);
    if (flat.size() != 2) {
      error = "input_paths must contain exactly two file paths";
      return false;
    }
    out.design0.emplace_back(flat[0]);
    out.design1.emplace_back(flat[1]);
  }

  if (out.design0.empty() || out.design1.empty()) {
    error = "input_paths must contain at least one file per design";
    return false;
  }
  return true;
}

static void logDesignPaths(const char* label,
                           const std::vector<std::string>& paths) {
  if (paths.empty()) {
    SPDLOG_INFO("{}: <none>", label);
    return;
  }
  if (paths.size() == 1) {
    SPDLOG_INFO("{}: {}", label, paths.front());
    return;
  }
  std::ostringstream oss;
  oss << label << ": ";
  for (size_t i = 0; i < paths.size(); ++i) {
    if (i) oss << ", ";
    oss << paths[i];
  }
  SPDLOG_INFO("{}", oss.str());
}

static std::vector<std::filesystem::path> toPathVector(
    const std::vector<std::string>& inputs) {
  std::vector<std::filesystem::path> out;
  out.reserve(inputs.size());
  for (const auto& s : inputs) {
    out.emplace_back(s);
  }
  return out;
}

int KeplerFormalMain(int argc, char** argv) {
  using namespace std::chrono;
  enum class FormatType { VERILOG, SYSTEMVERILOG, NAJA_IF };
  const auto cleanupNajaState = []() {
    naja::DNL::destroy();
    if (NLUniverse::get()) {
      NLUniverse::get()->destroy();
    }
  };
  cleanupNajaState();
  struct CleanupGuard {
    decltype(cleanupNajaState) cleanup;
    ~CleanupGuard() { cleanup(); }
  } cleanupGuard{cleanupNajaState};

  // Default values
  FormatType inputFormatType = FormatType::VERILOG;
  DesignInputs designInputs;
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
  bool dumpCnf = false;
  bool verilogPreprocessing = false;
  std::string dumpCnfPath;

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
        if (!validateConfigKeys(cfg)) {
          return EXIT_FAILURE;
        }

        // format
        if (cfg["format"] && cfg["format"].IsScalar()) {
          std::string fmt = cfg["format"].as<std::string>();
          if (fmt == "naja_if")
            inputFormatType = FormatType::NAJA_IF;
          else if (fmt == "verilog" || fmt == "v")
            inputFormatType = FormatType::VERILOG;
          else if (fmt == "systemverilog" || fmt == "sv")
            inputFormatType = FormatType::SYSTEMVERILOG;
          else {
            SPDLOG_CRITICAL("Unrecognized format in config: {}", fmt);
            return EXIT_FAILURE;
          }
        }

        // input_paths
        {
          std::string inputError;
          if (!parseConfigInputPaths(cfg["input_paths"], designInputs, inputError)) {
            SPDLOG_CRITICAL("Invalid input_paths in config: {}", inputError);
            return EXIT_FAILURE;
          }
        }

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

        // cnf_export
        if (cfg["cnf_export"] && cfg["cnf_export"].IsScalar()) {
          dumpCnf = cfg["cnf_export"].as<bool>();
        }

        // cnf_export_path (optional)
        if (cfg["cnf_export_path"] && cfg["cnf_export_path"].IsScalar()) {
          dumpCnfPath = cfg["cnf_export_path"].as<std::string>();
        }

        // verilog_preprocessing (optional)
        if (cfg["verilog_preprocessing"] && cfg["verilog_preprocessing"].IsScalar()) {
          verilogPreprocessing = cfg["verilog_preprocessing"].as<bool>();
        }

        // solver (glucose | kissat)
        if (cfg["solver"] && cfg["solver"].IsScalar()) {
          std::string solver = cfg["solver"].as<std::string>();
          if (solver == "glucose") {
            KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::GLUCOSE);
          } else if (solver == "kissat") {
            KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::KISSAT);
          } else {
            SPDLOG_CRITICAL("Unrecognized solver in config: {}", solver);
            return EXIT_FAILURE;
          }
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
    if (formatType == "-naja_if") {
      inputFormatType = FormatType::NAJA_IF;
    } else if (formatType == "-verilog") {
      inputFormatType = FormatType::VERILOG;
    } else if (formatType == "-systemverilog" || formatType == "-sv") {
      inputFormatType = FormatType::SYSTEMVERILOG;
    } else {
      SPDLOG_CRITICAL("Unrecognized input format type: {}", formatType);
      return EXIT_FAILURE;
    }

    bool explicitDesignFlags = false;
    std::vector<std::string>* currentDesign = nullptr;
    bool currentLiberty = false;
    std::vector<std::string> inputPaths;

    for (int i = 2; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--design1") {
        explicitDesignFlags = true;
        currentDesign = &designInputs.design0;
        currentLiberty = false;
        continue;
      }
      if (arg == "--design2") {
        explicitDesignFlags = true;
        currentDesign = &designInputs.design1;
        currentLiberty = false;
        continue;
      }
      if (arg == "--liberty" || arg == "--lib") {
        explicitDesignFlags = true;
        currentDesign = nullptr;
        currentLiberty = true;
        continue;
      }
      if (arg == "--verilog_preprocessing") {
        verilogPreprocessing = true;
        continue;
      }

      if (!arg.empty() && arg[0] == '-') {
        SPDLOG_CRITICAL("Unrecognized option: {}", arg);
        return EXIT_FAILURE;
      }

      if (explicitDesignFlags) {
        if (currentLiberty) {
          libertyFiles.push_back(arg);
        } else if (currentDesign) {
          currentDesign->push_back(arg);
        } else {
          // LCOV_EXCL_START
          SPDLOG_CRITICAL("Provide --design1 or --design2 before netlist paths");
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
      } else {
        inputPaths.emplace_back(arg);
      }
    }

    if (!explicitDesignFlags) {
      if (inputPaths.size() > 2) {
        for (size_t i = 2; i < inputPaths.size(); ++i)
          libertyFiles.push_back(inputPaths[i]);
      }
      if (inputPaths.size() >= 1) designInputs.design0.emplace_back(inputPaths[0]);
      if (inputPaths.size() >= 2) designInputs.design1.emplace_back(inputPaths[1]);
    }
  }

  // Configure logging level
  auto console = spdlog::get("console");
  if (!console) {
    console = spdlog::stdout_color_mt("console");
  }
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

  SPDLOG_INFO("KEPLER FORMAL: Run.");
  std::string inputFormatName = "VERILOG";
  if (inputFormatType == FormatType::NAJA_IF) inputFormatName = "SNL";
  if (inputFormatType == FormatType::SYSTEMVERILOG) inputFormatName = "SYSTEMVERILOG";
  SPDLOG_INFO("Input format: {}", inputFormatName);
  logDesignPaths("Netlist 1", designInputs.design0);
  logDesignPaths("Netlist 2", designInputs.design1);

  // Basic validation
  if (designInputs.design0.empty() || designInputs.design1.empty()) {
    SPDLOG_CRITICAL("Need two input netlist paths (one per design)");
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }
  if (inputFormatType == FormatType::NAJA_IF &&
      (designInputs.design0.size() != 1 || designInputs.design1.size() != 1)) {
    SPDLOG_CRITICAL("SNL input only supports one file per design");
    return EXIT_FAILURE;
  }

  auto solverType = KEPLER_FORMAL::Config::getSolverType();
  SPDLOG_INFO("Solver: {}",
              solverType == KEPLER_FORMAL::Config::SolverType::KISSAT ? "KISSAT" : "GLUCOSE");
  if (!libertyFiles.empty()) {
    for (const auto& lf : libertyFiles) SPDLOG_INFO("Library: {}", lf);
  }

  // --------------------------------------------------------------------------
  // 2. Load two netlists via Cap’n Proto (or via VRL constructor)
  // --------------------------------------------------------------------------
  naja::NL::SNLDesign* top0 = nullptr;
  naja::NL::SNLDesign* top1 = nullptr;
  try {
    NLUniverse::create();
    NLDB* db0 = nullptr;
    bool primitivesAreLoaded = false;

    auto loadLibraries = [&](NLDB* db) -> bool {
      if (libertyFiles.empty()) {
        // LCOV_EXCL_START
        return false;
        // LCOV_EXCL_STOP
      }
      auto primitivesLibrary =
          NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
      for (const auto& libraryFile : libertyFiles) {
        std::filesystem::path libraryPath(libraryFile);
        const auto extension = libraryPath.extension();
        const bool isLibertyFile =
            extension == ".lib" ||
            (extension == ".gz" && libraryPath.stem().extension() == ".lib");
        SPDLOG_INFO("Loading library file: {}", libraryFile);
        if (extension == ".py") {
          SNLPyLoader::loadPrimitives(primitivesLibrary, libraryPath);
        } else if (isLibertyFile) {
          SNLLibertyConstructor constructor(primitivesLibrary);
          constructor.construct(libraryPath);
        } else {
          SPDLOG_CRITICAL("Unsupported library file extension: {}", libraryPath.string());
          return false;
        }
      }
      return true;
    };

    if (!libertyFiles.empty()) {
      db0 = NLDB::create(NLUniverse::get());
      primitivesAreLoaded = loadLibraries(db0);
      if (!primitivesAreLoaded) {
        return EXIT_FAILURE;
      }
    }

    if (inputFormatType == FormatType::VERILOG ||
        inputFormatType == FormatType::SYSTEMVERILOG) {
      if (!db0) {
        db0 = NLDB::create(NLUniverse::get());
      }
      const auto design0Paths = toPathVector(designInputs.design0);
      SPDLOG_INFO("Parsing {} file(s) for design 1",
                  inputFormatType == FormatType::SYSTEMVERILOG ? "systemverilog" : "verilog");
      auto designLibrary = NLLibrary::create(db0, NLName("DESIGN"));
      if (inputFormatType == FormatType::SYSTEMVERILOG) {
        SNLSVConstructor constructor(designLibrary);
        constructor.construct(design0Paths);
      } else {
        SNLVRLConstructor constructor(designLibrary);
        constructor.config_.preprocessEnabled_ = verilogPreprocessing;
        constructor.construct(design0Paths);
      }
      auto top = SNLUtils::findTop(designLibrary);
      if (top) {
        db0->setTopDesign(top);
        SPDLOG_INFO("Found top design: {}", top->getString());
      } else {
        SPDLOG_CRITICAL("No top design was found after parsing input");
        return EXIT_FAILURE;
      }
    } else {
      SPDLOG_INFO("Loading Naja IF: {}", designInputs.design0[0]);
      naja::NL::SNLCapnP::LoadingConfiguration config;
      config.primitiveConflictPolicy_ =
          primitivesAreLoaded
              ? naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::PreferExisting
              : naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::ForbidConflicts;
      db0 = SNLCapnP::load(designInputs.design0[0].c_str(), config);
      if (!db0) {
        // LCOV_EXCL_START
        SPDLOG_CRITICAL("Failed to load Naja IF: {}", designInputs.design0[0]);
        return EXIT_FAILURE;
        // LCOV_EXCL_STOP
      }
    }

    // get db0 top
    top0 = db0->getTopDesign();
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
      if (!loadLibraries(db1)) {
        // LCOV_EXCL_START
        return EXIT_FAILURE;
        // LCOV_EXCL_STOP
      }
    }

    if (inputFormatType == FormatType::VERILOG ||
        inputFormatType == FormatType::SYSTEMVERILOG) {
      if (!db1) {
        db1 = NLDB::create(NLUniverse::get());
      }
      const auto design1Paths = toPathVector(designInputs.design1);
      SPDLOG_INFO("Parsing {} file(s) for design 2",
                  inputFormatType == FormatType::SYSTEMVERILOG ? "systemverilog" : "verilog");
      auto designLibrary = NLLibrary::create(db1, NLName("DESIGN"));
      if (inputFormatType == FormatType::SYSTEMVERILOG) {
        SNLSVConstructor constructor(designLibrary);
        constructor.construct(design1Paths);
      } else {
        SNLVRLConstructor constructor(designLibrary);
        constructor.config_.preprocessEnabled_ = verilogPreprocessing;
        constructor.construct(design1Paths);
      }
      auto top = SNLUtils::findTop(designLibrary);
      if (top) {
        db1->setTopDesign(top);
        SPDLOG_INFO("Found top design: {}", top->getString());
      } else {
        SPDLOG_CRITICAL("No top design was found after parsing input");
        return EXIT_FAILURE;
      }
    } else {
      SPDLOG_INFO("Loading Naja IF: {}", designInputs.design1[0]);
      naja::NL::SNLCapnP::LoadingConfiguration config;
      config.primitiveConflictPolicy_ =
          primitivesAreLoaded
              ? naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::PreferExisting
              : naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::ForbidConflicts;
      db1 = SNLCapnP::load(designInputs.design1[0].c_str(), config);
      if (!db1) {
        // LCOV_EXCL_START
        SPDLOG_CRITICAL("Failed to load Naja IF: {}", designInputs.design1[0]);
        return EXIT_FAILURE;
        // LCOV_EXCL_STOP
      }
    }

    // get db1 top
    top1 = db1->getTopDesign();
    if (!top1) {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("Top design not set for second netlist");
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  } catch (const std::exception& e) {
    SPDLOG_CRITICAL("Netlist loading failed: {}", e.what());
    return EXIT_FAILURE;
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
        if (dumpCnf) {
          std::string scopeName = sanitizeFileToken(scopes.first->getName().getString());
          std::string outPath = dumpCnfPath.empty()
                                    ? ("miter_" + scopeName + ".cnf")
                                    : dumpCnfPath;
          MiterScope.setCnfDump(true, outPath);
        }
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
      if (dumpCnf) {
        const std::string outPath = dumpCnfPath.empty() ? "miter.cnf" : dumpCnfPath;
        MiterS.setCnfDump(true, outPath);
      }
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

#ifndef KEPLER_FORMAL_NO_MAIN
int main(int argc, char** argv) { return KeplerFormalMain(argc, argv); }
#endif
