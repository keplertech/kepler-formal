// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <optional>
#include <cctype>
#include <stdexcept>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
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
#include "strategy/SequentialEquivalenceStrategy.h"

static void print_usage(const char* prog) {
  SPDLOG_INFO(
      "Usage: {} [--config <file>] | <-naja_if/-verilog/-systemverilog/-sv> "
      "[-v <LEC|SEC>] [-k <max-k>] [--sec-engine <LEGACY|KINDUCTION|IMC|PDR>] <netlist1> <netlist2> [<library-file>...] | "
      "<-naja_if/-verilog/-systemverilog/-sv> --design1 <file...> --design2 "
      "<file...> [--liberty <library-file>...] [-v <LEC|SEC>] [-k <max-k>] [--sec-engine <LEGACY|KINDUCTION|IMC|PDR>] "
      "[--no-sec-uncomputable-seq-boundary] [--compact] "
      "[--report-skipped-pos] | "
      "-systemverilog/-sv [--sv_design1_flist <file>] [--sv_design1_top <name>] "
      "[--sv_design2_flist <file>] [--sv_design2_top <name>] [-v <LEC|SEC>] [-k <max-k>] [--sec-engine <LEGACY|KINDUCTION|IMC|PDR>] "
      "[--design1 <file...>] [--design2 <file...>] "
      "[--no-sec-uncomputable-seq-boundary] [--compact] "
      "[--report-skipped-pos]",
      prog);
}

static std::vector<std::string> yamlToVector(const YAML::Node& node) {
  std::vector<std::string> out;
  if (!node) {
    return out;
  }
  if (!node.IsSequence()) {
    return out;
  }
  for (const auto& n : node) {
    if (n.IsScalar()) {
      out.emplace_back(n.as<std::string>());
    }
  }
  return out;
}

static bool isPythonLoaderPath(const std::string& path) {
  return std::filesystem::path(path).extension() == ".py";
}

enum class VerificationMode {
  LEC,
  SEC,
};

static std::string toUpperCopy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
    return static_cast<char>(std::toupper(ch));
  });
  return value;
}

static bool parseVerificationModeToken(const std::string& token,
                                       VerificationMode& mode,
                                       std::string& error) {
  const std::string normalized = toUpperCopy(token);
  if (normalized == "LEC") {
    mode = VerificationMode::LEC;
    return true;
  }
  if (normalized == "SEC") {
    mode = VerificationMode::SEC;
    return true;
  }
  error = "expected LEC or SEC, got `" + token + "`";
  return false;
}

static const char* verificationModeName(VerificationMode mode) {
  switch (mode) {
    case VerificationMode::SEC:
      return "SEC";
    case VerificationMode::LEC:
    default:
      return "LEC";
  }
}

static bool parseSecEngineToken(const std::string& token,
                                KEPLER_FORMAL::SEC::SecEngine& engine,
                                std::string& error) {
  // Keep the binary-level selector intentionally small for now so the
  // user-facing SEC modes stay explicit and predictable.
  const std::string normalized = toUpperCopy(token);
  if (normalized == "LEGACY") {
    engine = KEPLER_FORMAL::SEC::SecEngine::Legacy;
    return true;
  }
  if (normalized == "KINDUCTION" || normalized == "K_INDUCTION" ||
      normalized == "CLASSIC_K_INDUCTION") {
    engine = KEPLER_FORMAL::SEC::SecEngine::KInduction;
    return true;
  }
  if (normalized == "IMC") {
    engine = KEPLER_FORMAL::SEC::SecEngine::Imc;
    return true;
  }
  if (normalized == "PDR") {
    engine = KEPLER_FORMAL::SEC::SecEngine::Pdr;
    return true;
  }
  error = "expected LEGACY, KINDUCTION, IMC, or PDR, got `" + token + "`";
  return false;
}

static const char* secEngineName(KEPLER_FORMAL::SEC::SecEngine engine) {
  switch (engine) {
    case KEPLER_FORMAL::SEC::SecEngine::KInduction:
      return "KINDUCTION";
    case KEPLER_FORMAL::SEC::SecEngine::Imc:
      return "IMC";
    case KEPLER_FORMAL::SEC::SecEngine::Pdr:
      return "PDR";
    case KEPLER_FORMAL::SEC::SecEngine::Legacy:
    default:
      return "LEGACY";
  }
}

static spdlog::level::level_enum parseLogLevel(const std::string& logLevel) {
  if (logLevel == "debug") {
    return spdlog::level::debug;
  }
  return spdlog::level::info;
}

static std::string chooseRunLogFilePath(const std::string& requestedPath) {
  if (requestedPath.empty()) {
    int logIndex = 0;
    while (true) {
      const std::string candidate = "miter_log_" + std::to_string(logIndex) + ".txt";
      std::ifstream infile(candidate);
      if (infile.good()) {
        ++logIndex;
        continue;
      }
      return candidate;
    }
  }

  std::filesystem::path path(requestedPath);
  const auto parent = path.parent_path();
  if (!parent.empty()) {
    std::error_code ec;
    std::filesystem::create_directories(parent, ec);
    if (ec) {
      std::cerr << "Warning: failed to create log directory '" << parent.string()
                << "': " << ec.message() << " (" << ec.value()
                << "). Using default SEC log path.\n";
      return chooseRunLogFilePath("");
    }
  }
  return path.string();
}

static std::string configureMainLogger(const std::string& logLevel,
                                       bool enableFileLog,
                                       const std::string& requestedLogFilePath) {
  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());

  std::string chosenLogFile;
  if (enableFileLog) {
    chosenLogFile = chooseRunLogFilePath(requestedLogFilePath);
    try {
      sinks.push_back(
          std::make_shared<spdlog::sinks::basic_file_sink_mt>(chosenLogFile, true));
    } catch (const spdlog::spdlog_ex& ex) {
      std::cerr << "Warning: failed to create SEC log file '" << chosenLogFile
                << "': " << ex.what() << ". Logging will continue on stdout only.\n";
      chosenLogFile.clear();
    }
  }

  spdlog::drop("kepler_formal_main_logger");
  auto logger = std::make_shared<spdlog::logger>(
      "kepler_formal_main_logger", sinks.begin(), sinks.end());
  const auto level = parseLogLevel(logLevel);
  logger->set_level(level);
  logger->flush_on(spdlog::level::info);
  spdlog::set_default_logger(logger);
  spdlog::set_level(level);
  return chosenLogFile;
}

static bool parseMaxKToken(const std::string& token,
                           size_t& maxK,
                           std::string& error) {
  if (token.empty()) {
    error = "max_k must not be empty";
    return false;
  }
  if (!std::all_of(token.begin(), token.end(), [](unsigned char ch) { return std::isdigit(ch); })) {
    error = "max_k must be a non-negative integer";
    return false;
  }
  try {
    maxK = static_cast<size_t>(std::stoull(token));
  } catch (const std::exception&) {
    error = "max_k is out of range";
    return false;
  }
  return true;
}

static bool validateConfigKeys(const YAML::Node& cfg) {
  if (!cfg || !cfg.IsMap()) {
    return true;
  }
  static const std::unordered_set<std::string> kAllowedKeys = {
      "format",
      "verification",
      "max_k",
      "sec_engine",
      "sec_uncomputable_seq_as_boundary",
      "input_paths",
      "liberty_files",
      "py_tech_files",
      "verilog_preprocessing",
      "log_level",
      "log_file",
      "use_scopes",
      "clean_scopes",
      "cnf_export",
      "cnf_export_path",
      "dump_cnf",
      "dump_cnf_path",
      "compact_mode",
      "report_skipped_pos",
      "solver",
      "sv_design1_flist",
      "sv_design2_flist",
      "sv_design1_top",
      "sv_design2_top",
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

struct SystemVerilogDesignOptions {
  std::optional<std::string> flist;
  std::optional<std::string> top;
};

struct SystemVerilogOptions {
  SystemVerilogDesignOptions design0;
  SystemVerilogDesignOptions design1;
};

static bool parseConfigInputPaths(const YAML::Node& node,
                                  DesignInputs& out,
                                  std::string& error) {
  out.design0.clear();
  out.design1.clear();
  if (!node) {
    // LCOV_EXCL_START
    // KeplerFormalMain only calls parseConfigInputPaths when cfg["input_paths"] exists.
    error = "Missing input_paths in config";
    return false;
    // LCOV_EXCL_STOP
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
        if (i == 0) {
          out.design0.emplace_back(n.as<std::string>());
        } else {
          out.design1.emplace_back(n.as<std::string>());
        }
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
    if (i) {
      oss << ", ";
    }
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

static bool applySystemVerilogConfigOption(const YAML::Node& cfg,
                                           const char* key,
                                           std::optional<std::string>& target,
                                           std::string& error) {
  const auto node = cfg[key];
  if (!node) {
    return true;
  }
  if (!node.IsScalar()) {
    error = std::string(key) + " must be a scalar";
    return false;
  }
  const auto value = node.as<std::string>();
  if (value.empty()) {
    error = std::string(key) + " must not be empty";
    return false;
  }
  target = value;
  return true;
}

static bool validateSystemVerilogOptions(const SystemVerilogOptions& options,
                                         std::string& error) {
  const auto validateDesign = [&](const SystemVerilogDesignOptions& designOptions,
                                  const char* designLabel) {
    if (designOptions.top && designOptions.top->empty()) {
      // LCOV_EXCL_START
      // Public config and CLI parsing already reject empty SystemVerilog values.
      error = std::string(designLabel) + " top must not be empty";
      return false;
      // LCOV_EXCL_STOP
    }
    if (designOptions.flist && designOptions.flist->empty()) {
      // LCOV_EXCL_START
      // Public config and CLI parsing already reject empty SystemVerilog values.
      error = std::string(designLabel) + " flist must not be empty";
      return false;
      // LCOV_EXCL_STOP
    }
    return true;
  };

  return validateDesign(options.design0, "design1") &&
         validateDesign(options.design1, "design2");
}

static bool hasSystemVerilogSources(const std::vector<std::string>& designInputs,
                                    const SystemVerilogDesignOptions& designOptions) {
  return !designInputs.empty() || designOptions.flist.has_value();
}

static std::vector<std::filesystem::path> buildSystemVerilogInputPaths(
    const std::vector<std::string>& designInputs,
    const SystemVerilogDesignOptions& designOptions,
    std::vector<std::filesystem::path>& temporaryFiles) {
  std::vector<std::filesystem::path> svInputPaths = toPathVector(designInputs);

  const auto quotePathForSlangCommandFile = [](const std::filesystem::path& path) {
    std::string quoted;
    quoted.reserve(path.string().size() + 2);
    quoted.push_back('"');
    for (const auto c : path.string()) {
      if (c == '\\' || c == '"') {
        quoted.push_back('\\');
      }
      quoted.push_back(c);
    }
    quoted.push_back('"');
    return quoted;
  };

  if (designOptions.top) {
    const auto temporaryTopFlistPath =
        std::filesystem::temp_directory_path() /
        std::filesystem::path(
            "kepler_formal_sv_top_" +
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
            ".f");
    std::ofstream topFlist(temporaryTopFlistPath, std::ios::out | std::ios::trunc);
    if (!topFlist) {
      throw std::runtime_error(
          "Failed to create temporary SystemVerilog command file: " +
          temporaryTopFlistPath.string());
    }
    topFlist << "--top " << *designOptions.top << "\n";
    if (designOptions.flist) {
      topFlist << "-f "
               << quotePathForSlangCommandFile(std::filesystem::path(*designOptions.flist))
               << "\n";
    }
    for (const auto& svInputPath : svInputPaths) {
      topFlist << quotePathForSlangCommandFile(svInputPath) << "\n";
    }
    topFlist.close();
    temporaryFiles.push_back(temporaryTopFlistPath);
    return {std::filesystem::path("-f"), temporaryTopFlistPath};
  }

  if (designOptions.flist) {
    svInputPaths.insert(svInputPaths.begin(), std::filesystem::path(*designOptions.flist));
    svInputPaths.insert(svInputPaths.begin(), std::filesystem::path("-f"));
  }

  return svInputPaths;
}

static KEPLER_FORMAL::MiterStrategy::CompactSnapshot captureCompactSnapshot(
    const KEPLER_FORMAL::BuildPrimaryOutputClauses& builder) {
  KEPLER_FORMAL::MiterStrategy::CompactSnapshot snapshot;
  snapshot.inputs.reserve(builder.getInputs().size());
  for (const auto input : builder.getInputs()) {
    snapshot.inputs.emplace_back(builder.getInputs2InputsIDs().at(input));
  }
  snapshot.outputs.reserve(builder.getOutputs().size());
  for (const auto output : builder.getOutputs()) {
    snapshot.outputs.emplace_back(builder.getOutputs2OutputsIDs().at(output));
  }
  for (auto* expr : builder.getPOs()) {
    snapshot.POs.push_back(expr);
  }
  return snapshot;
}

int KeplerFormalMain(int argc, char** argv) {
  using namespace std::chrono;
  enum class FormatType { VERILOG, SYSTEMVERILOG, NAJA_IF };
  constexpr size_t kDefaultSecMaxK = 32;
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
  SystemVerilogOptions systemVerilogOptions;
  std::vector<std::string> libertyFiles;
  std::vector<std::string> pythonFiles;
  std::string logLevel = "info";
  VerificationMode verificationMode = VerificationMode::LEC;
  KEPLER_FORMAL::SEC::SecEngine secEngine = KEPLER_FORMAL::SEC::SecEngine::Legacy;
  bool secEngineExplicit = false;
  size_t secMaxK = kDefaultSecMaxK;
  bool secMaxKExplicit = false;
  bool secTreatUncomputableSeqAsBoundary = true;

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
  bool compactMode = false;
  bool reportSkippedPOs = false;
  bool verilogPreprocessing = false;
  std::string dumpCnfPath;

  KEPLER_FORMAL::Config::setReportSkippedPOs(false);
  KEPLER_FORMAL::Config::setSecTreatUncomputableSeqAsBoundary(true);

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
          if (fmt == "naja_if") {
            inputFormatType = FormatType::NAJA_IF;
          } else if (fmt == "verilog" || fmt == "v") {
            inputFormatType = FormatType::VERILOG;
          } else if (fmt == "systemverilog" || fmt == "sv") {
            inputFormatType = FormatType::SYSTEMVERILOG;
          } else {
            SPDLOG_CRITICAL("Unrecognized format in config: {}", fmt);
            return EXIT_FAILURE;
          }
        }

        if (cfg["verification"]) {
          if (!cfg["verification"].IsScalar()) {
            SPDLOG_CRITICAL("verification must be a scalar");
            return EXIT_FAILURE;
          }
          std::string verificationError;
          if (!parseVerificationModeToken(
                  cfg["verification"].as<std::string>(), verificationMode, verificationError)) {
            SPDLOG_CRITICAL("Invalid verification mode in config: {}", verificationError);
            return EXIT_FAILURE;
          }
        }

        if (cfg["max_k"]) {
          if (!cfg["max_k"].IsScalar()) {
            SPDLOG_CRITICAL("max_k must be a scalar");
            return EXIT_FAILURE;
          }
          std::string maxKError;
          if (!parseMaxKToken(cfg["max_k"].as<std::string>(), secMaxK, maxKError)) {
            SPDLOG_CRITICAL("Invalid max_k in config: {}", maxKError);
            return EXIT_FAILURE;
          }
          secMaxKExplicit = true;
        }

        if (cfg["sec_engine"]) {
          if (!cfg["sec_engine"].IsScalar()) {
            SPDLOG_CRITICAL("sec_engine must be a scalar");
            return EXIT_FAILURE;
          }
          std::string secEngineError;
          if (!parseSecEngineToken(
                  cfg["sec_engine"].as<std::string>(), secEngine, secEngineError)) {
            SPDLOG_CRITICAL("Invalid sec_engine in config: {}", secEngineError);
            return EXIT_FAILURE;
          }
          secEngineExplicit = true;
        }

        if (cfg["sec_uncomputable_seq_as_boundary"]) {
          if (!cfg["sec_uncomputable_seq_as_boundary"].IsScalar()) {
            SPDLOG_CRITICAL(
                "sec_uncomputable_seq_as_boundary must be a scalar");
            return EXIT_FAILURE;
          }
          secTreatUncomputableSeqAsBoundary =
              cfg["sec_uncomputable_seq_as_boundary"].as<bool>();
        }

        // input_paths
        if (cfg["input_paths"]) {
          std::string inputError;
          if (!parseConfigInputPaths(cfg["input_paths"], designInputs, inputError)) {
            SPDLOG_CRITICAL("Invalid input_paths in config: {}", inputError);
            return EXIT_FAILURE;
          }
        }

        // liberty_files
        libertyFiles = yamlToVector(cfg["liberty_files"]);
        pythonFiles = yamlToVector(cfg["py_tech_files"]);

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

        // compact_mode
        if (cfg["compact_mode"] && cfg["compact_mode"].IsScalar()) {
          compactMode = cfg["compact_mode"].as<bool>();
        }

        // report_skipped_pos
        if (cfg["report_skipped_pos"] && cfg["report_skipped_pos"].IsScalar()) {
          reportSkippedPOs = cfg["report_skipped_pos"].as<bool>();
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

        std::string svConfigError;
        if (!applySystemVerilogConfigOption(
                cfg, "sv_design1_flist", systemVerilogOptions.design0.flist, svConfigError) ||
            !applySystemVerilogConfigOption(
                cfg, "sv_design2_flist", systemVerilogOptions.design1.flist, svConfigError) ||
            !applySystemVerilogConfigOption(
                cfg, "sv_design1_top", systemVerilogOptions.design0.top, svConfigError) ||
            !applySystemVerilogConfigOption(
                cfg, "sv_design2_top", systemVerilogOptions.design1.top, svConfigError)) {
          SPDLOG_CRITICAL("Invalid SystemVerilog config option: {}", svConfigError);
          return EXIT_FAILURE;
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
    bool formatFound = false;
    int parseStart = 1;
    while (parseStart < argc) {
      std::string arg = argv[parseStart];
      if (arg == "--help" || arg == "-h") {
        print_usage(argv[0]);
        return EXIT_SUCCESS;
      }
      if (arg == "-v" || arg == "--verification") {
        if (parseStart + 1 >= argc) {
          SPDLOG_CRITICAL("Missing verification mode after {}", arg);
          return EXIT_FAILURE;
        }
        std::string verificationError;
        if (!parseVerificationModeToken(argv[parseStart + 1], verificationMode, verificationError)) {
          SPDLOG_CRITICAL("Invalid verification mode: {}", verificationError);
          return EXIT_FAILURE;
        }
        parseStart += 2;
        continue;
      }
      if (arg == "-k" || arg == "--max-k") {
        if (parseStart + 1 >= argc) {
          SPDLOG_CRITICAL("Missing max_k after {}", arg);
          return EXIT_FAILURE;
        }
        std::string maxKError;
        if (!parseMaxKToken(argv[parseStart + 1], secMaxK, maxKError)) {
          SPDLOG_CRITICAL("Invalid max_k: {}", maxKError);
          return EXIT_FAILURE;
        }
        secMaxKExplicit = true;
        parseStart += 2;
        continue;
      }
      if (arg == "--sec-engine") {
        if (parseStart + 1 >= argc) {
          SPDLOG_CRITICAL("Missing SEC engine after {}", arg);
          return EXIT_FAILURE;
        }
        std::string secEngineError;
        if (!parseSecEngineToken(argv[parseStart + 1], secEngine, secEngineError)) {
          SPDLOG_CRITICAL("Invalid SEC engine: {}", secEngineError);
          return EXIT_FAILURE;
        }
        secEngineExplicit = true;
        parseStart += 2;
        continue;
      }
      if (arg == "--sec-uncomputable-seq-boundary") {
        secTreatUncomputableSeqAsBoundary = true;
        ++parseStart;
        continue;
      }
      if (arg == "--no-sec-uncomputable-seq-boundary") {
        secTreatUncomputableSeqAsBoundary = false;
        ++parseStart;
        continue;
      }
      if (arg == "-naja_if") {
        inputFormatType = FormatType::NAJA_IF;
        ++parseStart;
        formatFound = true;
        break;
      }
      if (arg == "-verilog") {
        inputFormatType = FormatType::VERILOG;
        ++parseStart;
        formatFound = true;
        break;
      }
      if (arg == "-systemverilog" || arg == "-sv") {
        inputFormatType = FormatType::SYSTEMVERILOG;
        ++parseStart;
        formatFound = true;
        break;
      }
      SPDLOG_CRITICAL("Unrecognized option before input format type: {}", arg);
      return EXIT_FAILURE;
    }
    if (!formatFound) {
      SPDLOG_CRITICAL("Missing input format type");
      return EXIT_FAILURE;
    }

    bool explicitDesignFlags = false;
    std::vector<std::string>* currentDesign = nullptr;
    bool currentLiberty = false;
    std::vector<std::string> inputPaths;

    for (int i = parseStart; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "-v" || arg == "--verification") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing verification mode after {}", arg);
          return EXIT_FAILURE;
        }
        std::string verificationError;
        if (!parseVerificationModeToken(argv[++i], verificationMode, verificationError)) {
          SPDLOG_CRITICAL("Invalid verification mode: {}", verificationError);
          return EXIT_FAILURE;
        }
        continue;
      }
      if (arg == "-k" || arg == "--max-k") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing max_k after {}", arg);
          return EXIT_FAILURE;
        }
        std::string maxKError;
        if (!parseMaxKToken(argv[++i], secMaxK, maxKError)) {
          SPDLOG_CRITICAL("Invalid max_k: {}", maxKError);
          return EXIT_FAILURE;
        }
        secMaxKExplicit = true;
        continue;
      }
      if (arg == "--sec-engine") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing SEC engine after {}", arg);
          return EXIT_FAILURE;
        }
        std::string secEngineError;
        if (!parseSecEngineToken(argv[++i], secEngine, secEngineError)) {
          SPDLOG_CRITICAL("Invalid SEC engine: {}", secEngineError);
          return EXIT_FAILURE;
        }
        secEngineExplicit = true;
        continue;
      }
      if (arg == "--sec-uncomputable-seq-boundary") {
        secTreatUncomputableSeqAsBoundary = true;
        continue;
      }
      if (arg == "--no-sec-uncomputable-seq-boundary") {
        secTreatUncomputableSeqAsBoundary = false;
        continue;
      }
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
      if (arg == "--compact") {
        compactMode = true;
        continue;
      }
      if (arg == "--report-skipped-pos") {
        reportSkippedPOs = true;
        continue;
      }
      if (arg == "--sv_design1_flist" || arg == "--sv_design2_flist" ||
          arg == "--sv_design1_top" || arg == "--sv_design2_top") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing value after {}", arg);
          return EXIT_FAILURE;
        }
        const std::string value = argv[++i];
        if (value.empty()) {
          SPDLOG_CRITICAL("Empty value provided for {}", arg);
          return EXIT_FAILURE;
        }
        if (arg == "--sv_design1_flist") {
          systemVerilogOptions.design0.flist = value;
        } else if (arg == "--sv_design2_flist") {
          systemVerilogOptions.design1.flist = value;
        } else if (arg == "--sv_design1_top") {
          systemVerilogOptions.design0.top = value;
        } else {
          systemVerilogOptions.design1.top = value;
        }
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
        for (size_t i = 2; i < inputPaths.size(); ++i) {
          libertyFiles.push_back(inputPaths[i]);
        }
      }
      if (inputPaths.size() >= 1) {
        designInputs.design0.emplace_back(inputPaths[0]);
      }
      if (inputPaths.size() >= 2) {
        designInputs.design1.emplace_back(inputPaths[1]);
      }
    }
  }

  const std::string runLogFilePath = configureMainLogger(
      logLevel, verificationMode == VerificationMode::SEC, logFileName);

  SPDLOG_INFO("KEPLER FORMAL: Run.");
  std::string inputFormatName = "VERILOG";
  if (inputFormatType == FormatType::NAJA_IF) {
    inputFormatName = "SNL";
  }
  if (inputFormatType == FormatType::SYSTEMVERILOG) {
    inputFormatName = "SYSTEMVERILOG";
  }
  SPDLOG_INFO("Input format: {}", inputFormatName);
  if (!runLogFilePath.empty()) {
    SPDLOG_INFO("Run log: {}", runLogFilePath);
  }
  logDesignPaths("Netlist 1", designInputs.design0);
  logDesignPaths("Netlist 2", designInputs.design1);

  // Basic validation
  if (inputFormatType == FormatType::SYSTEMVERILOG) {
    if (!hasSystemVerilogSources(designInputs.design0, systemVerilogOptions.design0) ||
        !hasSystemVerilogSources(designInputs.design1, systemVerilogOptions.design1)) {
      SPDLOG_CRITICAL(
          "Need SystemVerilog input sources for both designs (files and/or per-design flists)");
      print_usage(argv[0]);
      return EXIT_FAILURE;
    }
  } else if (designInputs.design0.empty() || designInputs.design1.empty()) {
    SPDLOG_CRITICAL("Need two input netlist paths (one per design)");
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }
  if (inputFormatType == FormatType::NAJA_IF &&
      (designInputs.design0.size() != 1 || designInputs.design1.size() != 1)) {
    SPDLOG_CRITICAL("SNL input only supports one file per design");
    return EXIT_FAILURE;
  }
  if (verificationMode == VerificationMode::LEC && secMaxKExplicit) {
    SPDLOG_CRITICAL("max_k/-k is only supported with SEC verification");
    return EXIT_FAILURE;
  }
  if (verificationMode == VerificationMode::LEC && secEngineExplicit) {
    SPDLOG_CRITICAL("sec_engine/--sec-engine is only supported with SEC verification");
    return EXIT_FAILURE;
  }
  if (verificationMode == VerificationMode::SEC) {
    if (compactMode) {
      SPDLOG_CRITICAL("SEC verification does not support compact mode");
      return EXIT_FAILURE;
    }
    if (useScopes || cleanScopes) {
      SPDLOG_CRITICAL("SEC verification does not support scope extraction/cleaning");
      return EXIT_FAILURE;
    }
    if (dumpCnf) {
      SPDLOG_CRITICAL("SEC verification does not support CNF export");
      return EXIT_FAILURE;
    }
    if (reportSkippedPOs) {
      SPDLOG_CRITICAL("SEC verification does not support skipped PO reporting");
      return EXIT_FAILURE;
    }
  }
  for (const auto& libraryFile : libertyFiles) {
    if (isPythonLoaderPath(libraryFile)) {
      SPDLOG_CRITICAL(
          "Python primitive loader {} must be provided through YAML config key "
          "py_tech_files, not liberty_files/--liberty inputs",
          libraryFile);
      return EXIT_FAILURE;
    }
  }
  if (inputFormatType != FormatType::SYSTEMVERILOG &&
      (systemVerilogOptions.design0.flist || systemVerilogOptions.design0.top ||
       systemVerilogOptions.design1.flist || systemVerilogOptions.design1.top)) {
    SPDLOG_CRITICAL("SystemVerilog design options are only valid with -systemverilog/-sv input");
    return EXIT_FAILURE;
  }
  std::string svValidationError;
  if (!validateSystemVerilogOptions(systemVerilogOptions, svValidationError)) {
    // LCOV_EXCL_START
    // validateSystemVerilogOptions only fails for empty values, which the public
    // CLI/config parsing already rejects before reaching this point.
    SPDLOG_CRITICAL("Invalid SystemVerilog options: {}", svValidationError);
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
  }

  auto solverType = KEPLER_FORMAL::Config::getSolverType();
  KEPLER_FORMAL::Config::setReportSkippedPOs(reportSkippedPOs);
  KEPLER_FORMAL::Config::setSecTreatUncomputableSeqAsBoundary(
      secTreatUncomputableSeqAsBoundary);
  SPDLOG_INFO("Solver: {}",
              solverType == KEPLER_FORMAL::Config::SolverType::KISSAT ? "KISSAT" : "GLUCOSE");
  SPDLOG_INFO("Verification: {}", verificationModeName(verificationMode));
  if (verificationMode == VerificationMode::SEC) {
    SPDLOG_INFO("SEC max_k: {}", secMaxK);
    SPDLOG_INFO("SEC engine: {}", secEngineName(secEngine));
    SPDLOG_INFO(
        "SEC uncomputable sequentials: {}",
        secTreatUncomputableSeqAsBoundary ? "boundary abstraction"
                                          : "strict failure");
  }
  SPDLOG_INFO("Compact mode: {}", compactMode ? "enabled" : "disabled");
  SPDLOG_INFO("Skipped PO reports: {}", reportSkippedPOs ? "enabled" : "disabled");
  if (!libertyFiles.empty()) {
    for (const auto& lf : libertyFiles) SPDLOG_INFO("Library: {}", lf);
  }
  if (!pythonFiles.empty()) {
    for (const auto& pf : pythonFiles) SPDLOG_INFO("Python library: {}", pf);
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
      if (libertyFiles.empty() && pythonFiles.empty()) {
        // LCOV_EXCL_START
        return false;
        // LCOV_EXCL_STOP
      }
      auto primitivesLibrary =
          NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
      for (const auto& libraryFile : libertyFiles) {
        std::filesystem::path libraryPath(libraryFile);
        SPDLOG_INFO("Loading library file: {}", libraryFile);
        SNLLibertyConstructor constructor(primitivesLibrary);
        constructor.construct(libraryPath);
      }
      for (const auto& pythonFile : pythonFiles) {
        std::filesystem::path pythonPath(pythonFile);
        SPDLOG_INFO("Loading python primitive file: {}", pythonFile);
        SNLPyLoader::loadPrimitives(primitivesLibrary, pythonPath);
      }
      return true;
    };

    auto loadOneDesign = [&](const std::vector<std::string>& designPaths,
                             const SystemVerilogDesignOptions& designOptions,
                             int designIndex,
                             int dbID) -> NLDB* {
      NLDB* db = nullptr;
      bool primitivesLoadedForDesign = false;

      if (!libertyFiles.empty() || !pythonFiles.empty()) {
        db = NLDB::create(NLUniverse::get());
        primitivesLoadedForDesign = loadLibraries(db);
        if (!primitivesLoadedForDesign) {
          throw std::runtime_error("Failed to load library files");
        }
      }

      if (inputFormatType == FormatType::VERILOG ||
          inputFormatType == FormatType::SYSTEMVERILOG) {
        if (!db) {
          db = NLDB::create(NLUniverse::get());
        }
        db->setID(dbID);
        SPDLOG_INFO("Parsing {} file(s) for design {}",
                    inputFormatType == FormatType::SYSTEMVERILOG ? "systemverilog" : "verilog",
                    designIndex + 1);
        auto designLibrary = NLLibrary::create(db, NLName("DESIGN"));
        if (inputFormatType == FormatType::SYSTEMVERILOG) {
          SNLSVConstructor constructor(designLibrary);
          std::vector<std::filesystem::path> temporaryFiles;
	          const auto svInputPaths =
	              buildSystemVerilogInputPaths(designPaths, designOptions, temporaryFiles);
	          try {
	            constructor.construct(svInputPaths);
              // LCOV_EXCL_START
	          } catch (...) {
	            for (const auto& temporaryFile : temporaryFiles) {
	              std::error_code ec;
	              std::filesystem::remove(temporaryFile, ec);
	            }
	            throw;
	          }
              // LCOV_EXCL_STOP
          for (const auto& temporaryFile : temporaryFiles) {
            std::error_code ec;
            std::filesystem::remove(temporaryFile, ec);
          }
        } else {
          SNLVRLConstructor constructor(designLibrary);
          constructor.config_.preprocessEnabled_ = verilogPreprocessing;
          constructor.construct(toPathVector(designPaths));
        }
	        auto top = SNLUtils::findTop(designLibrary);
	        if (!top) {
            // LCOV_EXCL_START
	          throw std::runtime_error("No top design was found after parsing input");
            // LCOV_EXCL_STOP
	        }
        db->setTopDesign(top);
        SPDLOG_INFO("Found top design: {}", top->getString());
      } else {
        SPDLOG_INFO("Loading Naja IF: {}", designPaths[0]);
        naja::NL::SNLCapnP::LoadingConfiguration config;
        config.primitiveConflictPolicy_ =
            primitivesLoadedForDesign
                ? naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::PreferExisting
                : naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::ForbidConflicts;
	        db = SNLCapnP::load(designPaths[0].c_str(), config);
	        if (!db) {
            // LCOV_EXCL_START
	          throw std::runtime_error("Failed to load Naja IF: " + designPaths[0]);
            // LCOV_EXCL_STOP
	        }
        db->setID(dbID);
      }

	      if (!db->getTopDesign()) {
          // LCOV_EXCL_START
	        throw std::runtime_error("Top design not set for loaded netlist");
          // LCOV_EXCL_STOP
	      }
      return db;
    };

    auto buildCompactSnapshotForTop =
        [&](naja::NL::SNLDesign* top,
            const char* designLabel) {
          KEPLER_FORMAL::BuildPrimaryOutputClauses builder;
          NLUniverse::get()->setTopDesign(top);
          naja::DNL::destroy();
          builder.collect();
          SPDLOG_INFO("Collected {} PIs for {}", builder.getInputs().size(), designLabel);
          SPDLOG_INFO("Collected {} POs for {}", builder.getOutputs().size(), designLabel);
          auto inputs = builder.getInputs();
          auto outputs = builder.getOutputs();
          builder.setInputs(inputs);
          builder.setOutputs(outputs);
          builder.build();
          return captureCompactSnapshot(builder);
        };

    if (compactMode && !useScopes) {
      NLDB* compactDb0 =
          loadOneDesign(designInputs.design0, systemVerilogOptions.design0, 0, 2);
      top0 = compactDb0->getTopDesign();
      auto snapshot0 = buildCompactSnapshotForTop(top0, "design 0");
      naja::DNL::destroy();
      compactDb0->destroy();
      top0 = nullptr;

      NLDB* compactDb1 =
          loadOneDesign(designInputs.design1, systemVerilogOptions.design1, 1, 1);
      top1 = compactDb1->getTopDesign();
      auto snapshot1 = buildCompactSnapshotForTop(top1, "design 1");
      naja::DNL::destroy();
      compactDb1->destroy();
      top1 = nullptr;

      try {
        KEPLER_FORMAL::MiterStrategy MiterS(nullptr, nullptr, logFileName);
        if (dumpCnf) {
          const std::string outPath = dumpCnfPath.empty() ? "miter.cnf" : dumpCnfPath;
          MiterS.setCnfDump(true, outPath);
        }
        if (MiterS.runCompactSnapshots(snapshot0, snapshot1)) {
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
      return EXIT_SUCCESS;
    }

    if (!libertyFiles.empty() || !pythonFiles.empty()) {
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
        std::vector<std::filesystem::path> temporaryFiles;
        const auto svInputPaths = buildSystemVerilogInputPaths(
            designInputs.design0, systemVerilogOptions.design0, temporaryFiles);
        try {
          constructor.construct(svInputPaths);
        } catch (...) {
          for (const auto& temporaryFile : temporaryFiles) {
            std::error_code ec;
            std::filesystem::remove(temporaryFile, ec);
          }
          throw;
        }
        for (const auto& temporaryFile : temporaryFiles) {
          std::error_code ec;
          std::filesystem::remove(temporaryFile, ec);
        }
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
    if (!libertyFiles.empty() || !pythonFiles.empty()) {
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
        std::vector<std::filesystem::path> temporaryFiles;
        const auto svInputPaths = buildSystemVerilogInputPaths(
            designInputs.design1, systemVerilogOptions.design1, temporaryFiles);
        try {
          constructor.construct(svInputPaths);
        } catch (...) {
          for (const auto& temporaryFile : temporaryFiles) {
            std::error_code ec;
            std::filesystem::remove(temporaryFile, ec);
          }
          throw;
        }
        for (const auto& temporaryFile : temporaryFiles) {
          std::error_code ec;
          std::filesystem::remove(temporaryFile, ec);
        }
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
  if (verificationMode == VerificationMode::SEC) {
    try {
      KEPLER_FORMAL::SEC::SequentialEquivalenceStrategy strategy(
          top0, top1, solverType, secEngine);
      const auto result = strategy.run(secMaxK);
      if (result.totalOutputs != 0) {
        SPDLOG_INFO(
            "SEC output coverage: {:.2f}% ({}/{} covered/existing outputs).",
            result.outputCoveragePercent(),
            result.coveredOutputs,
            result.totalOutputs);
      }
      if (!result.skippedObservedOutputs.empty()) {
        std::ostringstream skippedOutputs;
        for (const auto& skippedOutput : result.skippedObservedOutputs) {
          skippedOutputs << "  - " << skippedOutput << "\n";
        }
        SPDLOG_INFO(
            "SEC skipped observed outputs due to connectivity issues "
            "(no-driver or multi-driver only):\n{}",
            skippedOutputs.str());
      }
      if (!result.abstractedSequentialBoundaries.empty()) {
        std::ostringstream abstractedBoundaries;
        for (const auto& abstractedBoundary : result.abstractedSequentialBoundaries) {
          abstractedBoundaries << "  - " << abstractedBoundary << "\n";
        }
        SPDLOG_INFO(
            "SEC abstracted uncomputable sequential interfaces as boundaries:\n{}",
            abstractedBoundaries.str());
      }
      switch (result.status) {
        case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Equivalent:
          SPDLOG_INFO("No difference was found. SEC proved equivalence at k = {}.", result.bound);
          return EXIT_SUCCESS;
        case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Different:
          SPDLOG_INFO("Difference was found. SEC found a counterexample at k = {}.", result.bound);
          if (!result.reason.empty()) {
            SPDLOG_INFO("SEC counterexample details:\n{}", result.reason);
          }
          return EXIT_SUCCESS;
        case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Inconclusive:
          SPDLOG_CRITICAL("SEC was inconclusive up to max_k = {}: {}", secMaxK, result.reason);
          return EXIT_FAILURE;
        case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Unsupported:
        default:
          SPDLOG_CRITICAL("SEC cannot run on this design pair: {}", result.reason);
          return EXIT_FAILURE;
      }
    } catch (const std::exception& e) {
      SPDLOG_ERROR("SEC workflow failed: {}", e.what());
      return EXIT_FAILURE;
    }
  } else if (inputFormatType == FormatType::NAJA_IF && useScopes) {
    KEPLER_FORMAL::MiterStrategy MiterS(top0, top1);
    MiterS.init(false);
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
        if (MiterScope.run(compactMode)) {
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
      }
        // LCOV_EXCL_STOP
    }
  } else {
    try {
      KEPLER_FORMAL::MiterStrategy MiterS(top0, top1, logFileName);
      if (dumpCnf) {
        const std::string outPath = dumpCnfPath.empty() ? "miter.cnf" : dumpCnfPath;
        MiterS.setCnfDump(true, outPath);
      }
      MiterS.init();
      if (MiterS.run(compactMode)) {
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
