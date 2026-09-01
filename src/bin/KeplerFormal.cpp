// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <cctype>
#include <stdexcept>
#include <unordered_set>
#include <utility>
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
#ifndef KEPLER_FORMAL_NO_PY_TECH
#include "SNLPyLoader.h"
#endif
#include "SNLSVConstructor.h"
#include "SNLVRLConstructor.h"
#include "SNLVRLDumper.h"
#include "SNLBusTerm.h"
#include "SNLInstance.h"
#include "SNLTerm.h"
#include "SNLUtils.h"
#include "ScopeExtraction.h"
#include "Config.h"
#include "KeplerFormalDriver.h"
#include "KeplerFormalUtils.h"
#include "Tree2BoolExpr.h"
#include "model/SequentialDesignModel.h"
#include "strategy/SequentialEquivalenceStrategy.h"

#if defined(__SANITIZE_ADDRESS__)
#define KEPLER_FORMAL_ASAN_BUILD 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define KEPLER_FORMAL_ASAN_BUILD 1
#endif
#endif

static const char* kBoundaryTermsReport = "boundary_terms.txt";
static const char* kSkippedResetUnanchoredPOReport =
    "skipped_reset_unanchored_pos.txt";
static const char* kSkippedMultiClockDomainPOReport =
    "skipped_multi_clock_domain_pos.txt";
static const char* kSkippedOpaqueCellPOReport =
    "skipped_opaque_cells_pos.txt";

static void addNajaPythonPath(const char* argv0) {
  if (!argv0 || !*argv0) {
    return;
  }

  std::filesystem::path executable(argv0);
  std::error_code ec;
  if (!executable.has_parent_path()) {
    if (const char* path = std::getenv("PATH")) {
      std::istringstream paths(path);
      std::string directory;
#ifdef _WIN32
      constexpr char pathSeparator = ';';
#else
      constexpr char pathSeparator = ':';
#endif
      while (std::getline(paths, directory, pathSeparator)) {
        auto candidate = std::filesystem::path(directory) / executable;
        if (std::filesystem::exists(candidate, ec)) {
          executable = std::move(candidate);
          break;
        }
        ec.clear();
      }
    }
  }

  executable = std::filesystem::weakly_canonical(executable, ec);
  if (ec || executable.parent_path().empty()) {
    return;
  }

  std::string pythonPath = executable.parent_path().string();
  if (const char* current = std::getenv("PYTHONPATH"); current && *current) {
#ifdef _WIN32
    pythonPath += ';';
#else
    pythonPath += ':';
#endif
    pythonPath += current;
  }
#ifdef _WIN32
  if (_putenv_s("PYTHONPATH", pythonPath.c_str()) != 0) {
#else
  if (setenv("PYTHONPATH", pythonPath.c_str(), 1) != 0) {
#endif
    throw std::runtime_error("Cannot configure PYTHONPATH for Naja primitives");  // LCOV_EXCL_LINE
  }
}

// LCOV_EXCL_START
static void print_usage(const char* prog) {
  SPDLOG_INFO(
  // LCOV_EXCL_STOP
      "Usage: {} [--config <file>] | <-naja_if/-verilog/-systemverilog/-sv/-sv2v> "
      "[-v <lec|sec>] [-k <max-k>] [--sec-engine <k_induction|imc|pdr>] [--sec-encoding <binary|dual_rail_steady>] "
      "[--verilog_design1_top <name>] [--verilog_design2_top <name>] "
      "<netlist1> <netlist2> [<library-file>...] | "
      "<-naja_if/-verilog/-systemverilog/-sv/-sv2v> --design1 <file...> --design2 "
      "<file...> [--verilog_design1_top <name>] [--verilog_design2_top <name>] [--liberty <library-file>...] [-v <lec|sec>] [-k <max-k>] [--sec-engine <k_induction|imc|pdr>] [--sec-encoding <binary|dual_rail_steady>] "
      "[--allow-boundary-mismatch] [--compact] "
      "[--report-skipped-pos] | "
      "-systemverilog/-sv [--sv_design1_flist <file>] [--sv_design1_top <name>] "
      "[--sv_design2_flist <file>] [--sv_design2_top <name>] [-v <lec|sec>] [-k <max-k>] [--sec-engine <k_induction|imc|pdr>] [--sec-encoding <binary|dual_rail_steady>] "
      "[--design1 <file...>] [--design2 <file...>] "
      "[--allow-boundary-mismatch] [--compact] "
      "[--report-skipped-pos]",
      prog);
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

static std::vector<std::string> yamlToVector(const YAML::Node& node) {
  std::vector<std::string> out;
  if (!node) {
    return out;
  }
  if (!node.IsSequence()) {
    // LCOV_EXCL_START
    return out;
    // LCOV_EXCL_STOP
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
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

enum class VerificationMode {
  LEC,
  SEC,
};

static bool parseVerificationModeToken(const std::string& token,
                                       VerificationMode& mode,
                                       std::string& error) {
  if (token == "lec") {
    // LCOV_EXCL_START
    mode = VerificationMode::LEC;
    return true;
    // LCOV_EXCL_STOP
  }
  if (token == "sec") {
    mode = VerificationMode::SEC;
    return true;
  }
  // LCOV_EXCL_START
  error = "expected lec or sec, got `" + token + "`";
  return false;
  // LCOV_EXCL_STOP
}

static const char* verificationModeName(VerificationMode mode) {
  switch (mode) {
    case VerificationMode::SEC:
      return "sec";
    // LCOV_EXCL_START
    case VerificationMode::LEC:
    // LCOV_EXCL_STOP
    default:
      // LCOV_EXCL_START
      return "lec";
      // LCOV_EXCL_STOP
  }
}

static bool parseSecEngineToken(const std::string& token,
                                KEPLER_FORMAL::SEC::SecEngine& engine,
                                std::string& error) {
  // Keep the binary-level selector intentionally small for now so the
  // user-facing SEC modes stay explicit and predictable.
  if (token == "k_induction") {
    // LCOV_EXCL_START
    engine = KEPLER_FORMAL::SEC::SecEngine::KInduction;
    return true;
    // LCOV_EXCL_STOP
  }
  if (token == "imc") {
    // LCOV_EXCL_START
    engine = KEPLER_FORMAL::SEC::SecEngine::Imc;
    return true;
    // LCOV_EXCL_STOP
  }
  if (token == "pdr") {
    engine = KEPLER_FORMAL::SEC::SecEngine::Pdr;
    return true;
  }
  // LCOV_EXCL_START
  error = "expected k_induction, imc, or pdr, got `" + token + "`";
  return false;
  // LCOV_EXCL_STOP
}

static const char* secEngineName(KEPLER_FORMAL::SEC::SecEngine engine) {
  switch (engine) {
    case KEPLER_FORMAL::SEC::SecEngine::KInduction:
      // LCOV_EXCL_START
      return "k_induction";
      // LCOV_EXCL_STOP
    case KEPLER_FORMAL::SEC::SecEngine::Imc:
      // LCOV_EXCL_START
      return "imc";
      // LCOV_EXCL_STOP
    case KEPLER_FORMAL::SEC::SecEngine::Pdr:
      return "pdr";
    default:
      // LCOV_EXCL_START
      return "pdr";
      // LCOV_EXCL_STOP
  }
}

static bool parseSecEncodingToken(const std::string& token,
                                  KEPLER_FORMAL::SEC::SecEncoding& encoding,
                                  std::string& error) {
  if (token == "binary") {
    // LCOV_EXCL_START
    encoding = KEPLER_FORMAL::SEC::SecEncoding::Binary;  // LCOV_EXCL_LINE
    return true;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (token == "dual_rail_steady") {
    encoding = KEPLER_FORMAL::SEC::SecEncoding::DualRailSteady;
    return true;
  }
  // LCOV_EXCL_START
  error = "expected binary or dual_rail_steady, got `" + token + "`";
  return false;
  // LCOV_EXCL_STOP
}

static const char* secEncodingName(KEPLER_FORMAL::SEC::SecEncoding encoding) {
  switch (encoding) {
    case KEPLER_FORMAL::SEC::SecEncoding::DualRailSteady:
      return "dual_rail_steady";
    // LCOV_EXCL_START
    case KEPLER_FORMAL::SEC::SecEncoding::Binary:
    // LCOV_EXCL_STOP
    default:
      // LCOV_EXCL_START
      return "binary";
      // LCOV_EXCL_STOP
  }
}

static spdlog::level::level_enum parseLogLevel(const std::string& logLevel) {
  if (logLevel == "debug") {
    // LCOV_EXCL_START
    return spdlog::level::debug;
    // LCOV_EXCL_STOP
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

  // LCOV_EXCL_START
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
      // LCOV_EXCL_STOP
    }
  // LCOV_EXCL_START
  }
  return path.string();
  // LCOV_EXCL_STOP
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
      // LCOV_EXCL_START
      std::cerr << "Warning: failed to create SEC log file '" << chosenLogFile
                << "': " << ex.what() << ". Logging will continue on stdout only.\n";
      chosenLogFile.clear();
    }
    // LCOV_EXCL_STOP
  }

  spdlog::drop("kepler_formal_main_logger");
  auto logger = std::make_shared<spdlog::logger>(
      "kepler_formal_main_logger", sinks.begin(), sinks.end());
  const auto level = parseLogLevel(logLevel);
  logger->set_level(level);
  logger->flush_on(spdlog::level::info);
  spdlog::set_default_logger(logger);
  return chosenLogFile;
}

static bool parseNonNegativeSizeToken(const std::string& token,
                                      const char* optionName,
                                      size_t& value,
                                      std::string& error) {
  if (token.empty()) {
    error = std::string(optionName) + " must not be empty";
    return false;
  }
  if (token.find_first_not_of("0123456789") != std::string::npos) {
    error = std::string(optionName) + " must be a non-negative integer";
    return false;
  }
  try {
    value = static_cast<size_t>(std::stoull(token));
  } catch (const std::exception&) {
    error = std::string(optionName) + " is out of range";
    return false;
  }
  return true;
}

// LCOV_EXCL_START
static bool parseMaxKToken(const std::string& token,
// LCOV_EXCL_STOP
                           size_t& maxK,
                           std::string& error) {
  // LCOV_EXCL_START
  return parseNonNegativeSizeToken(token, "max_k", maxK, error);
}
// LCOV_EXCL_STOP

static bool validateConfigKeys(const YAML::Node& cfg) {
  if (!cfg || !cfg.IsMap()) {
    // LCOV_EXCL_START
    return true;
    // LCOV_EXCL_STOP
  }
  static const std::unordered_set<std::string> kAllowedKeys = {
      "format",
      "verification",
      "max_k",
      "sec_engine",
      "sec_encoding",
      "allow-boundary-mismatch",
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
      "po_cnf_export",
      "po_cnf_export_path",
      "dump_cnf",
      "dump_cnf_path",
      "compact_mode",
      "report_skipped_pos",
      "solver",
      "sv_design1_flist",
      "sv_design2_flist",
      "verilog_design1_top",
      "verilog_design2_top",
      "sv_design1_top",
      "sv_design2_top",
  };

  for (auto it = cfg.begin(); it != cfg.end(); ++it) {
    if (!it->first.IsScalar()) {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("Config key is not a scalar; invalid YAML key");
      return false;
      // LCOV_EXCL_STOP
    }
    const std::string key = it->first.as<std::string>();
    if (kAllowedKeys.find(key) == kAllowedKeys.end()) {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("Unknown config option: {}", key);
      return false;
      // LCOV_EXCL_STOP
    }
  }
  return true;
}

// LCOV_EXCL_START
std::string sanitizeFileToken(const std::string& input) {
  std::string out;
  out.reserve(input.size());
  for (unsigned char ch : input) {
    if (std::isalnum(ch) || ch == '_' || ch == '-' || ch == '.') {
      out.push_back(static_cast<char>(ch));
    } else {
      out.push_back('_');
      // LCOV_EXCL_STOP
    }
  }
  // LCOV_EXCL_START
  if (out.empty()) {
    out = "scope";
  }
  return out;
}
// LCOV_EXCL_STOP

namespace {

// LCOV_EXCL_START
std::string formatStringList(const std::vector<std::string>& values) {
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < values.size(); ++i) {
    if (i) {
      oss << ", ";
    }
    oss << values[i];
  }
  oss << "]";
  return oss.str();
}
// LCOV_EXCL_STOP

bool secInconclusiveStoppedBeforeMaxK(const std::string& reason) {
  return reason.find("budget") != std::string::npos ||
         reason.find("repair") != std::string::npos ||
         reason.find("projection") != std::string::npos ||
         reason.find("did not prove any observed output") != std::string::npos;
}

}  // namespace

// LCOV_EXCL_START
void writeBoundaryTermsReport(
// LCOV_EXCL_STOP
    const std::filesystem::path& reportPath,
    const std::vector<KEPLER_FORMAL::SEC::ExtractedBoundaryReportEntry>& reports) {
  // LCOV_EXCL_START
  if (reports.empty()) {
    return;
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  std::ofstream report(reportPath, std::ios::trunc);
  report << "# SEC boundary terms report\n";
  report << "# Categories:\n";
  report << "# - top_input / top_output: original top-level interface terms.\n\n";
  for (size_t i = 0; i < reports.size(); ++i) {
    const auto& entry = reports[i];
    report << "- design: " << entry.design << "\n";
    report << "  signal: " << entry.signal << "\n";
    report << "  roles: " << formatStringList(entry.roles) << "\n";
    if (!entry.connectivitySkip.empty()) {
      report << "  connectivity_skip: " << entry.connectivitySkip << "\n";
    }
    if (i + 1 != reports.size()) {
      report << "\n";
    }
  }
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
void writeResetUnanchoredSkippedOutputsReport(
// LCOV_EXCL_STOP
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs) {
  // LCOV_EXCL_START
  if (skippedOutputs.empty()) {
    return;
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  std::ofstream report(reportPath, std::ios::trunc);
  report << "# SEC reset-unanchored skipped observed outputs\n";
  report << "# These top outputs were removed from the proof surface because\n";
  report << "# their cones depend on internal state without an inductive\n";
  report << "# cross-design anchor. SEC does not assume internal flop equality\n";
  report << "# by name; only top-level interface signals are name-aligned.\n\n";
  for (const auto& skippedOutput : skippedOutputs) {
    report << "- " << skippedOutput << "\n";
    // LCOV_EXCL_STOP
  }
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
void writeMultiClockDomainSkippedOutputsReport(
// LCOV_EXCL_STOP
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs) {
  // LCOV_EXCL_START
  if (skippedOutputs.empty()) {
    return;
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  std::ofstream report(reportPath, std::ios::trunc);
  report << "# SEC multi-clock-domain skipped observed outputs\n";
  report << "# These top outputs were removed from the proof surface because\n";
  report << "# their cones span more than one extracted clock domain. CDC\n";
  report << "# modeling is intentionally outside this SEC pass, so the result\n";
  report << "# is reported as skipped coverage instead of assumed synchronous.\n\n";
  for (const auto& skippedOutput : skippedOutputs) {
    report << "- " << skippedOutput << "\n";
    // LCOV_EXCL_STOP
  }
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

void writeOpaqueCellSkippedOutputsReport(
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs) {
  if (skippedOutputs.empty()) {
    return;
  }

  std::ofstream report(reportPath, std::ios::trunc);
  report << "# SEC opaque-cell skipped top-level outputs\n";
  report << "# These outputs were not verified because backward cone traversal\n";
  report << "# reached an internal cell or pin without usable semantics. No free\n";
  report << "# or shared proof symbol was substituted for the opaque element.\n\n";
  for (const auto& skippedOutput : skippedOutputs) {
    report << "- " << skippedOutput << "\n";
  }
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

struct VerilogTopOptions {
  std::optional<std::string> design0;
  std::optional<std::string> design1;
};

static bool parseConfigInputPaths(const YAML::Node& node,
                                  DesignInputs& out,
                                  std::string& error) {
  out.design0.clear();
  out.design1.clear();
  if (!node) {
    // LCOV_EXCL_START
    // KeplerFormalMain only calls parseConfigInputPaths when cfg["input_paths"] exists.
    // LCOV_DISABLED_START
    error = "Missing input_paths in config";
    return false;
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  }
  if (!node.IsSequence()) {
    // LCOV_EXCL_START
    error = "input_paths must be a sequence";
    return false;
    // LCOV_EXCL_STOP
  }
  if (node.size() == 0) {
    // LCOV_EXCL_START
    error = "input_paths must contain at least two entries";
    return false;
    // LCOV_EXCL_STOP
  }

  const bool firstIsSequence = node[0].IsSequence();
  if (firstIsSequence) {
    // LCOV_EXCL_START
    if (node.size() != 2) {
      error = "input_paths must contain exactly two design entries";
      return false;
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    for (size_t i = 0; i < 2; ++i) {
      if (!node[i].IsSequence()) {
        error = "input_paths must be a sequence of sequences when using multi-file designs";
        return false;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      for (const auto& n : node[i]) {
        if (!n.IsScalar()) {
          error = "input_paths entries must be scalar file paths";
          return false;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        if (i == 0) {
          out.design0.emplace_back(n.as<std::string>());
        } else {
          out.design1.emplace_back(n.as<std::string>());
          // LCOV_EXCL_STOP
        }
      // LCOV_EXCL_START
      }
    }
  } else {
  // LCOV_EXCL_STOP
    std::vector<std::string> flat = yamlToVector(node);
    if (flat.size() != 2) {
      // LCOV_EXCL_START
      error = "input_paths must contain exactly two file paths";
      return false;
      // LCOV_EXCL_STOP
    }
    out.design0.emplace_back(flat[0]);
    out.design1.emplace_back(flat[1]);
  }

  return true;
}

static void logDesignPaths(const char* label,
                           const std::vector<std::string>& paths) {
  if (paths.empty()) {
    // LCOV_EXCL_START
    SPDLOG_INFO("{}: <none>", label);
    return;
    // LCOV_EXCL_STOP
  }
  if (paths.size() == 1) {
    SPDLOG_INFO("{}: {}", label, paths.front());
    return;
  }
  // LCOV_EXCL_START
  std::ostringstream oss;
  oss << label << ": ";
  for (size_t i = 0; i < paths.size(); ++i) {
    if (i) {
      oss << ", ";
    }
    oss << paths[i];
  }
  SPDLOG_INFO("{}", oss.str());
  // LCOV_EXCL_STOP
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

static std::string normalizeInputPathForComparison(const std::string& path) {
  std::error_code ec;
  const auto canonical = std::filesystem::weakly_canonical(path, ec);
  if (!ec) {
    return canonical.string();
  }

  // LCOV_EXCL_START
  // LCOV_DISABLED_START
  ec.clear();
  const auto absolute = std::filesystem::absolute(path, ec);
  if (!ec) {
    return absolute.lexically_normal().string();
    // LCOV_DISABLED_STOP
  }
  // LCOV_DISABLED_START
  return std::filesystem::path(path).lexically_normal().string();
  // LCOV_DISABLED_STOP
  // LCOV_EXCL_STOP
}

static std::vector<std::string> normalizeInputListForComparison(
    const std::vector<std::string>& inputs) {
  std::vector<std::string> normalized;
  normalized.reserve(inputs.size());
  for (const auto& input : inputs) {
    normalized.push_back(normalizeInputPathForComparison(input));
  }
  return normalized;
}

// LCOV_EXCL_START
static std::optional<std::string> normalizeOptionalInputPathForComparison(
// LCOV_EXCL_STOP
    const std::optional<std::string>& input) {
  // LCOV_EXCL_START
  if (!input.has_value()) {
    return std::nullopt;
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return normalizeInputPathForComparison(*input);
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
static bool sameSystemVerilogDesignOptions(
// LCOV_EXCL_STOP
    const SystemVerilogDesignOptions& lhs,
    const SystemVerilogDesignOptions& rhs) {
  // LCOV_EXCL_START
  return normalizeOptionalInputPathForComparison(lhs.flist) ==
             normalizeOptionalInputPathForComparison(rhs.flist) &&
         lhs.top == rhs.top;
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

static bool sameCompactSecDesignSpec(
    bool isSystemVerilog,
    const DesignInputs& designInputs,
    const SystemVerilogOptions& systemVerilogOptions,
    const VerilogTopOptions& verilogTopOptions) {
  if (normalizeInputListForComparison(designInputs.design0) !=
      normalizeInputListForComparison(designInputs.design1)) {
    return false;
  }
  // LCOV_EXCL_START
  if (!isSystemVerilog) {
    return verilogTopOptions.design0 == verilogTopOptions.design1;
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return sameSystemVerilogDesignOptions(
      systemVerilogOptions.design0, systemVerilogOptions.design1);
      // LCOV_EXCL_STOP
}

static naja::NL::SNLDesign* selectTopDesign(
    naja::NL::NLLibrary* library,
    const std::optional<std::string>& requestedTop,
    int designIndex) {
  if (!requestedTop) {
    return SNLUtils::findTop(library);
  }
  auto* top = library->getSNLDesign(NLName(*requestedTop));
  if (!top) {
    throw std::runtime_error(
        "Top module `" + *requestedTop + "` was not found in design " +
        std::to_string(designIndex + 1));
  }
  return top;
}

static bool applySystemVerilogConfigOption(const YAML::Node& cfg,
                                           const char* key,
                                           std::optional<std::string>& target,
                                           std::string& error) {
  const auto node = cfg[key];
  if (!node) {
    return true;
  }
  // LCOV_EXCL_START
  if (!node.IsScalar()) {
    error = std::string(key) + " must be a scalar";
    return false;
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  const auto value = node.as<std::string>();
  if (value.empty()) {
    error = std::string(key) + " must not be empty";
    return false;
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  target = value;
  return true;
  // LCOV_EXCL_STOP
}

static bool validateSystemVerilogOptions(const SystemVerilogOptions& options,
                                         std::string& error) {
  const auto validateDesign = [&](const SystemVerilogDesignOptions& designOptions,
                                  const char* designLabel) {
    if (designOptions.top && designOptions.top->empty()) {
      // LCOV_EXCL_START
      // Public config and CLI parsing already reject empty SystemVerilog values.
      // LCOV_DISABLED_START
      error = std::string(designLabel) + " top must not be empty";
      return false;
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    }
    if (designOptions.flist && designOptions.flist->empty()) {
      // LCOV_EXCL_START
      // Public config and CLI parsing already reject empty SystemVerilog values.
      // LCOV_DISABLED_START
      error = std::string(designLabel) + " flist must not be empty";
      return false;
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    }
    return true;
  };

  return validateDesign(options.design0, "design1") &&
         validateDesign(options.design1, "design2");
}

// LCOV_EXCL_START
static bool hasSystemVerilogSources(const std::vector<std::string>& designInputs,
// LCOV_EXCL_STOP
                                    const SystemVerilogDesignOptions& designOptions) {
  // LCOV_EXCL_START
  return !designInputs.empty() || designOptions.flist.has_value();
  // LCOV_EXCL_STOP
}

static bool isSimpleVerilogIdentifier(const std::string& name) {
  if (name.empty()) {
    return false;  // LCOV_EXCL_LINE
  }
  const auto isIdentifierChar = [](unsigned char c) {
    return std::isalnum(c) || c == '_' || c == '$';
  };
  const unsigned char first = static_cast<unsigned char>(name.front());
  if ((!std::isalpha(first) && first != '_') || first == '$') {
    return false;
  }
  return std::all_of(name.begin(), name.end(), [&](char c) {
    return isIdentifierChar(static_cast<unsigned char>(c));
  });
}

static std::string dumpVerilogIdentifier(const std::string& name) {
  if (isSimpleVerilogIdentifier(name)) {
    return name;
  }
  return "\\" + name + " ";
}

static const char* getVerilogDirection(naja::NL::SNLTerm::Direction direction) {
  switch (direction) {
    case naja::NL::SNLTerm::Direction::Input:
      return "input";
    case naja::NL::SNLTerm::Direction::Output:
      return "output";
    case naja::NL::SNLTerm::Direction::InOut:
    case naja::NL::SNLTerm::Direction::Undefined:
      return "inout";
  }
  return "inout";  // LCOV_EXCL_LINE
}

static void appendPrimitiveStubModules(
    naja::NL::NLLibrary* library,
    std::ofstream& out,
    std::unordered_set<std::string>& emittedNames,
    size_t& emittedCount) {
  if (library == nullptr) {
    return;  // LCOV_EXCL_LINE
  }

  for (auto* design : library->getSNLDesigns()) {
    if (design == nullptr || design->isUnnamed()) {
      continue;
    }
    const std::string modelName = design->getName().getString();
    if (modelName.empty() || !emittedNames.insert(modelName).second) {
      continue;
    }

    std::vector<naja::NL::SNLTerm*> terms;
    for (auto* term : design->getTerms()) {
      if (term != nullptr && !term->isUnnamed()) {
        terms.push_back(term);
      }
    }

    out << "(* blackbox *) module " << dumpVerilogIdentifier(modelName) << "(";
    for (size_t i = 0; i < terms.size(); ++i) {
      if (i != 0) {
        out << ", ";
      }
      out << dumpVerilogIdentifier(terms[i]->getName().getString());
    }
    out << ");\n";
    for (auto* term : terms) {
      out << "  " << getVerilogDirection(term->getDirection()) << " ";
      if (auto* busTerm = dynamic_cast<naja::NL::SNLBusTerm*>(term)) {
        out << "[" << busTerm->getMSB() << ":" << busTerm->getLSB() << "] ";
      }
      out << dumpVerilogIdentifier(term->getName().getString()) << ";\n";
    }
    out << "endmodule\n\n";
    ++emittedCount;
  }

  for (auto* childLibrary : library->getLibraries()) {
    appendPrimitiveStubModules(childLibrary, out, emittedNames, emittedCount);
  }
}

static std::filesystem::path writeSystemVerilogPrimitiveStubs(
    const std::vector<naja::NL::NLLibrary*>& primitiveLibraries,
    std::vector<std::filesystem::path>& temporaryFiles) {
  if (primitiveLibraries.empty()) {
    return {};
  }

  const auto stubPath =
      std::filesystem::temp_directory_path() /
      std::filesystem::path(
          "kepler_formal_sv2v_prims_" +
          std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
          ".sv");
  std::ofstream stubs(stubPath, std::ios::out | std::ios::trunc);
  if (!stubs) {
    throw std::runtime_error(
        "Failed to create temporary SystemVerilog primitive stub file: " +
        stubPath.string());
  }

  std::unordered_set<std::string> emittedNames;
  size_t emittedCount = 0;
  for (auto* primitiveLibrary : primitiveLibraries) {
    appendPrimitiveStubModules(primitiveLibrary, stubs, emittedNames, emittedCount);
  }
  stubs.close();

  if (emittedCount == 0) {
    std::error_code ec;
    std::filesystem::remove(stubPath, ec);
    return {};
  }

  temporaryFiles.push_back(stubPath);
  return stubPath;
}

static naja::NL::SNLDesign* findPrimitiveDesignInLibrary(
    naja::NL::NLLibrary* library,
    const naja::NL::NLName& name) {
  if (library == nullptr) {
    return nullptr;  // LCOV_EXCL_LINE
  }
  if (auto* primitive = library->getSNLDesign(name)) {
    return primitive;
  }
  for (auto* childLibrary : library->getLibraries()) {
    if (auto* primitive = findPrimitiveDesignInLibrary(childLibrary, name)) {
      return primitive;
    }
  }
  return nullptr;
}

static naja::NL::SNLDesign* findPrimitiveDesign(
    const std::vector<naja::NL::NLLibrary*>& primitiveLibraries,
    const naja::NL::NLName& name) {
  for (auto* primitiveLibrary : primitiveLibraries) {
    if (auto* primitive = findPrimitiveDesignInLibrary(primitiveLibrary, name)) {
      return primitive;
    }
  }
  return nullptr;
}

static void reconnectGeneratedPrimitiveStubs(
    naja::NL::NLLibrary* designLibrary,
    const std::vector<naja::NL::NLLibrary*>& primitiveLibraries) {
  if (designLibrary == nullptr || primitiveLibraries.empty()) {
    return;
  }

  for (auto* design : designLibrary->getSNLDesigns()) {
    if (design == nullptr || design->isPrimitive()) {
      continue;  // LCOV_EXCL_LINE
    }
    for (auto* instance : design->getInstances()) {
      auto* model = instance->getModel();
      if (model == nullptr || model->isPrimitive() || model->isUnnamed()) {
        continue;  // LCOV_EXCL_LINE
      }
      if (auto* primitive = findPrimitiveDesign(primitiveLibraries, model->getName())) {
        instance->setModel(primitive);
      }
    }
  }
}

// LCOV_EXCL_START
static std::vector<std::filesystem::path> buildSystemVerilogInputPaths(
// LCOV_EXCL_STOP
    const std::vector<std::string>& designInputs,
    const SystemVerilogDesignOptions& designOptions,
    std::vector<std::filesystem::path>& temporaryFiles,
    const std::vector<naja::NL::NLLibrary*>* primitiveLibraries = nullptr) {
  // LCOV_EXCL_START
  std::vector<std::filesystem::path> svInputPaths = toPathVector(designInputs);
  // LCOV_EXCL_STOP

  if (primitiveLibraries != nullptr) {
    const auto primitiveStubsPath =
        writeSystemVerilogPrimitiveStubs(*primitiveLibraries, temporaryFiles);
    if (!primitiveStubsPath.empty()) {
      svInputPaths.push_back(primitiveStubsPath);
    }
  }

  // LCOV_EXCL_START
  const auto quotePathForSlangCommandFile = [](const std::filesystem::path& path) {
    std::string quoted;
    quoted.reserve(path.string().size() + 2);
    quoted.push_back('"');
    for (const auto c : path.string()) {
      if (c == '\\' || c == '"') {
        quoted.push_back('\\');
      }
      quoted.push_back(c);
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    quoted.push_back('"');
    return quoted;
  };
  // LCOV_EXCL_STOP

  // LCOV_EXCL_START
  // Generated SV2V primitive stubs may be parsed with RTL that has no enclosing
  // timescale, so force a Slang command file when stubs are present.
  const bool addDefaultTimescale = primitiveLibraries != nullptr;
  if (designOptions.top || addDefaultTimescale) {
  // LCOV_EXCL_STOP
    const auto temporaryTopFlistPath =
        // LCOV_EXCL_START
        std::filesystem::temp_directory_path() /
        std::filesystem::path(
            "kepler_formal_sv_top_" +
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
            // LCOV_EXCL_STOP
            ".f");
    // LCOV_EXCL_START
    std::ofstream topFlist(temporaryTopFlistPath, std::ios::out | std::ios::trunc);
    if (!topFlist) {
      throw std::runtime_error(
          "Failed to create temporary SystemVerilog command file: " +
          temporaryTopFlistPath.string());
          // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    if (addDefaultTimescale) {
      topFlist << "--timescale 1ns/1ps\n";
    }
    if (designOptions.top) {
      topFlist << "--top " << *designOptions.top << "\n";
    }
    if (designOptions.flist) {
      topFlist << "-f "
               << quotePathForSlangCommandFile(std::filesystem::path(*designOptions.flist))
               << "\n";
    }
    for (const auto& svInputPath : svInputPaths) {
      topFlist << quotePathForSlangCommandFile(svInputPath) << "\n";
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    topFlist.close();
    temporaryFiles.push_back(temporaryTopFlistPath);
    return {std::filesystem::path("-f"), temporaryTopFlistPath};
  }
  // LCOV_EXCL_STOP

  // LCOV_EXCL_START
  if (designOptions.flist) {
    svInputPaths.insert(svInputPaths.begin(), std::filesystem::path(*designOptions.flist));
    svInputPaths.insert(svInputPaths.begin(), std::filesystem::path("-f"));
  }
  // LCOV_EXCL_STOP

  // LCOV_EXCL_START
  return svInputPaths;
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
static KEPLER_FORMAL::MiterStrategy::CompactSnapshot captureCompactSnapshot(
// LCOV_EXCL_STOP
    const KEPLER_FORMAL::BuildPrimaryOutputClauses& builder,
    std::vector<KEPLER_FORMAL::BuildPrimaryOutputClauses::PathKey> boundaryInputs) {
  // LCOV_EXCL_START
  KEPLER_FORMAL::MiterStrategy::CompactSnapshot snapshot;
  snapshot.boundaryInputs = std::move(boundaryInputs);
  snapshot.inputs.reserve(builder.getInputs().size());
  for (const auto input : builder.getInputs()) {
    snapshot.inputs.emplace_back(builder.getInputs2InputsIDs().at(input));
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  snapshot.outputs.reserve(builder.getOutputs().size());
  for (const auto output : builder.getOutputs()) {
    snapshot.outputs.emplace_back(builder.getOutputs2OutputsIDs().at(output));
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  for (auto* expr : builder.getPOs()) {
    snapshot.POs.push_back(expr);
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return snapshot;
}
// LCOV_EXCL_STOP

static int KeplerFormalMainImpl(
    int argc,
    char** argv,
    KEPLER_FORMAL::RunResult* runResult) {
  using namespace std::chrono;
  enum class FormatType { VERILOG, SYSTEMVERILOG, SV2V, NAJA_IF };
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

  if (runResult != nullptr) {
    *runResult = KEPLER_FORMAL::RunResult{};
  }

  // Default values
  FormatType inputFormatType = FormatType::VERILOG;
  DesignInputs designInputs;
  SystemVerilogOptions systemVerilogOptions;
  VerilogTopOptions verilogTopOptions;
  std::vector<std::string> libertyFiles;
  std::vector<std::string> pythonFiles;
  std::string logLevel = "info";
  VerificationMode verificationMode = VerificationMode::LEC;
  KEPLER_FORMAL::SEC::SecEngine secEngine = KEPLER_FORMAL::SEC::SecEngine::Pdr;
  KEPLER_FORMAL::SEC::SecEncoding secEncoding =
      KEPLER_FORMAL::SEC::SecEncoding::DualRailSteady;
  bool secEngineExplicit = false;
  bool secEncodingExplicit = false;
  size_t secMaxK = kDefaultSecMaxK;
  bool secMaxKExplicit = false;
  // Basic argument sanity
  if (argc < 2) {
    // LCOV_EXCL_START
    if (runResult != nullptr) {
      runResult->status = KEPLER_FORMAL::RunStatus::NoResult;
    }
    print_usage(argv[0]);  // LCOV_EXCL_LINE
    return EXIT_SUCCESS;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // Config mode (--config or -c) is exclusive with other command-line options.
  bool usedConfig = false;

  std::string logFileName;

  bool useScopes = false;
  bool cleanScopes = false;
  bool dumpCnf = false;
  bool dumpPoCnf = false;
  bool compactMode = false;
  bool allowBoundaryMismatch = false;
  bool reportSkippedPOs = false;
  bool verilogPreprocessing = false;
  std::string dumpCnfPath;
  std::string dumpPoCnfPath;

  KEPLER_FORMAL::Config::setReportSkippedPOs(false);

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--config" || a == "-c") {
      if (i + 1 >= argc) {
        // LCOV_EXCL_START
        SPDLOG_CRITICAL("Missing config file after {}", a);  // LCOV_EXCL_LINE
        return EXIT_FAILURE;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      if (argc != 3) {
        SPDLOG_CRITICAL(
            "Config mode cannot be combined with other command-line options");
        return EXIT_FAILURE;
      }
      const std::string cfgPath = argv[i + 1];
      try {
        YAML::Node cfg = YAML::LoadFile(cfgPath);
        if (!validateConfigKeys(cfg)) {
          // LCOV_EXCL_START
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }

        // format
        if (cfg["format"] && cfg["format"].IsScalar()) {
          std::string fmt = cfg["format"].as<std::string>();
          if (fmt == "naja_if") {
            // LCOV_EXCL_START
            inputFormatType = FormatType::NAJA_IF;
            // LCOV_EXCL_STOP
          } else if (fmt == "verilog" || fmt == "v") {
            inputFormatType = FormatType::VERILOG;
          } else if (fmt == "systemverilog" || fmt == "sv") {
            // LCOV_EXCL_START
            inputFormatType = FormatType::SYSTEMVERILOG;
          } else if (fmt == "sv2v") {
            inputFormatType = FormatType::SV2V;
          } else {
            SPDLOG_CRITICAL("Unrecognized format in config: {}", fmt);
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
        }

        if (cfg["verification"]) {
          if (!cfg["verification"].IsScalar()) {
            // LCOV_EXCL_START
            SPDLOG_CRITICAL("verification must be a scalar");
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
          std::string verificationError;
          if (!parseVerificationModeToken(
                  cfg["verification"].as<std::string>(), verificationMode, verificationError)) {
            // LCOV_EXCL_START
            SPDLOG_CRITICAL("Invalid verification mode in config: {}", verificationError);
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
        }

        if (cfg["max_k"]) {
          // LCOV_EXCL_START
          if (!cfg["max_k"].IsScalar()) {
            SPDLOG_CRITICAL("max_k must be a scalar");
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
          // LCOV_EXCL_START
          std::string maxKError;
          if (!parseMaxKToken(cfg["max_k"].as<std::string>(), secMaxK, maxKError)) {
            SPDLOG_CRITICAL("Invalid max_k in config: {}", maxKError);
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
          // LCOV_EXCL_START
          secMaxKExplicit = true;
        }
        // LCOV_EXCL_STOP

        if (cfg["sec_engine"]) {
          if (!cfg["sec_engine"].IsScalar()) {
            // LCOV_EXCL_START
            SPDLOG_CRITICAL("sec_engine must be a scalar");
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
          std::string secEngineError;
          if (!parseSecEngineToken(
                  cfg["sec_engine"].as<std::string>(), secEngine, secEngineError)) {
            // LCOV_EXCL_START
            SPDLOG_CRITICAL("Invalid sec_engine in config: {}", secEngineError);
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
          secEngineExplicit = true;
        }

        if (cfg["sec_encoding"]) {
          if (!cfg["sec_encoding"].IsScalar()) {
            // LCOV_EXCL_START
            SPDLOG_CRITICAL("sec_encoding must be a scalar");
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
          std::string secEncodingError;
          if (!parseSecEncodingToken(
                  cfg["sec_encoding"].as<std::string>(), secEncoding, secEncodingError)) {
            // LCOV_EXCL_START
            SPDLOG_CRITICAL("Invalid sec_encoding in config: {}", secEncodingError);
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
          secEncodingExplicit = true;  // LCOV_EXCL_LINE
        }

        // input_paths
        if (cfg["input_paths"]) {
          std::string inputError;
          if (!parseConfigInputPaths(cfg["input_paths"], designInputs, inputError)) {
            // LCOV_EXCL_START
            SPDLOG_CRITICAL("Invalid input_paths in config: {}", inputError);
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
        }

        // liberty_files
        libertyFiles = yamlToVector(cfg["liberty_files"]);
        pythonFiles = yamlToVector(cfg["py_tech_files"]);

        // log level
        if (cfg["log_level"] && cfg["log_level"].IsScalar()) {
          // LCOV_EXCL_START
          logLevel = cfg["log_level"].as<std::string>();
        }
        // LCOV_EXCL_STOP

        // Add log file name
        if (cfg["log_file"] && cfg["log_file"].IsScalar()) {
          // LCOV_EXCL_START
          logFileName = cfg["log_file"].as<std::string>();
        }
        // LCOV_EXCL_STOP
        
        // use_scopes
        if (cfg["use_scopes"] && cfg["use_scopes"].IsScalar()) {
          // LCOV_EXCL_START
          useScopes = cfg["use_scopes"].as<bool>();
        }
        // LCOV_EXCL_STOP

        // clean_scopes
        if (cfg["clean_scopes"] && cfg["clean_scopes"].IsScalar()) {
          // LCOV_EXCL_START
          cleanScopes = cfg["clean_scopes"].as<bool>();
        }
        // LCOV_EXCL_STOP

        // cnf_export
        if (cfg["cnf_export"] && cfg["cnf_export"].IsScalar()) {
          // LCOV_EXCL_START
          dumpCnf = cfg["cnf_export"].as<bool>();
        }
        // LCOV_EXCL_STOP

        // cnf_export_path (optional)
        if (cfg["cnf_export_path"] && cfg["cnf_export_path"].IsScalar()) {
          // LCOV_EXCL_START
          dumpCnfPath = cfg["cnf_export_path"].as<std::string>();
        }
        // LCOV_EXCL_STOP

        // po_cnf_export
        if (cfg["po_cnf_export"] && cfg["po_cnf_export"].IsScalar()) {
          // LCOV_EXCL_START
          dumpPoCnf = cfg["po_cnf_export"].as<bool>();
        }
        // LCOV_EXCL_STOP

        // po_cnf_export_path (optional)
        if (cfg["po_cnf_export_path"] && cfg["po_cnf_export_path"].IsScalar()) {
          // LCOV_EXCL_START
          dumpPoCnfPath = cfg["po_cnf_export_path"].as<std::string>();
        }
        // LCOV_EXCL_STOP

        // compact_mode
        if (cfg["compact_mode"] && cfg["compact_mode"].IsScalar()) {
          compactMode = cfg["compact_mode"].as<bool>();
        }

        if (cfg["allow-boundary-mismatch"] &&
            cfg["allow-boundary-mismatch"].IsScalar()) {
          allowBoundaryMismatch = cfg["allow-boundary-mismatch"].as<bool>();
        }

        // report_skipped_pos
        if (cfg["report_skipped_pos"] && cfg["report_skipped_pos"].IsScalar()) {
          // LCOV_EXCL_START
          reportSkippedPOs = cfg["report_skipped_pos"].as<bool>();
        }
        // LCOV_EXCL_STOP

        // verilog_preprocessing (optional)
        if (cfg["verilog_preprocessing"] && cfg["verilog_preprocessing"].IsScalar()) {
          // LCOV_EXCL_START
          verilogPreprocessing = cfg["verilog_preprocessing"].as<bool>();
        }
        // LCOV_EXCL_STOP

        // solver (glucose | kissat | cadical)
        if (cfg["solver"] && cfg["solver"].IsScalar()) {
          // LCOV_EXCL_START
          std::string solver = cfg["solver"].as<std::string>();
          if (solver == "glucose") {
            KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::GLUCOSE);
          } else if (solver == "kissat") {
            KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::KISSAT);
          } else if (solver == "cadical") {
            KEPLER_FORMAL::Config::setSolverType(KEPLER_FORMAL::Config::CADICAL);  // LCOV_EXCL_LINE
          } else {  // LCOV_EXCL_LINE
            SPDLOG_CRITICAL("Unrecognized solver in config: {}", solver);
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
        // LCOV_EXCL_START
        }
        // LCOV_EXCL_STOP

        std::string svConfigError;
        if (!applySystemVerilogConfigOption(
                cfg, "sv_design1_flist", systemVerilogOptions.design0.flist, svConfigError) ||
            !applySystemVerilogConfigOption(
                cfg, "sv_design2_flist", systemVerilogOptions.design1.flist, svConfigError) ||
            !applySystemVerilogConfigOption(
                cfg, "sv_design1_top", systemVerilogOptions.design0.top, svConfigError) ||
            !applySystemVerilogConfigOption(
                cfg, "sv_design2_top", systemVerilogOptions.design1.top, svConfigError) ||
            !applySystemVerilogConfigOption(
                cfg, "verilog_design1_top", verilogTopOptions.design0, svConfigError) ||
            !applySystemVerilogConfigOption(
                cfg, "verilog_design2_top", verilogTopOptions.design1, svConfigError)) {
          // LCOV_EXCL_START
          SPDLOG_CRITICAL("Invalid design config option: {}", svConfigError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }

        usedConfig = true;
      } catch (const std::exception& e) {
        // LCOV_EXCL_START
        SPDLOG_CRITICAL("Failed to parse config {}: {}", cfgPath, e.what());  // LCOV_EXCL_LINE
        return EXIT_FAILURE;  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      break;
    }
  }

  // If not using config, fall back to original CLI parsing
  if (!usedConfig) {
    // LCOV_EXCL_START
    bool formatFound = false;
    int parseStart = 1;
    while (parseStart < argc) {
      std::string arg = argv[parseStart];
      if (arg == "--help" || arg == "-h") {
        if (runResult != nullptr) {
          runResult->status = KEPLER_FORMAL::RunStatus::NoResult;
        }
        print_usage(argv[0]);
        return EXIT_SUCCESS;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "-v" || arg == "--verification") {
        if (parseStart + 1 >= argc) {
          SPDLOG_CRITICAL("Missing verification mode after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        std::string verificationError;
        if (!parseVerificationModeToken(argv[parseStart + 1], verificationMode, verificationError)) {
          SPDLOG_CRITICAL("Invalid verification mode: {}", verificationError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        parseStart += 2;
        continue;
      }
      if (arg == "-k" || arg == "--max-k") {
        if (parseStart + 1 >= argc) {
          SPDLOG_CRITICAL("Missing max_k after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        std::string maxKError;
        if (!parseMaxKToken(argv[parseStart + 1], secMaxK, maxKError)) {
          SPDLOG_CRITICAL("Invalid max_k: {}", maxKError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        secMaxKExplicit = true;
        parseStart += 2;
        continue;
      }
      if (arg == "--sec-engine") {
        if (parseStart + 1 >= argc) {
          SPDLOG_CRITICAL("Missing SEC engine after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        std::string secEngineError;
        if (!parseSecEngineToken(argv[parseStart + 1], secEngine, secEngineError)) {
          SPDLOG_CRITICAL("Invalid SEC engine: {}", secEngineError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        secEngineExplicit = true;
        parseStart += 2;
        continue;
      }
      if (arg == "--sec-encoding") {
        if (parseStart + 1 >= argc) {
          SPDLOG_CRITICAL("Missing SEC encoding after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        std::string secEncodingError;
        if (!parseSecEncodingToken(argv[parseStart + 1], secEncoding, secEncodingError)) {
          SPDLOG_CRITICAL("Invalid SEC encoding: {}", secEncodingError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        secEncodingExplicit = true;
        parseStart += 2;
        continue;
      }
      if (arg == "--allow-boundary-mismatch") {
        allowBoundaryMismatch = true;
        ++parseStart;
        continue;
      }
      // LCOV_EXCL_START
      if (arg == "-naja_if") {
        inputFormatType = FormatType::NAJA_IF;
        ++parseStart;
        formatFound = true;
        break;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "-verilog") {
        inputFormatType = FormatType::VERILOG;
        ++parseStart;
        formatFound = true;
        break;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "-systemverilog" || arg == "-sv") {
        inputFormatType = FormatType::SYSTEMVERILOG;
        ++parseStart;
        formatFound = true;
        break;
        // LCOV_EXCL_STOP
      }
      if (arg == "-sv2v") {
        inputFormatType = FormatType::SV2V;
        ++parseStart;
        formatFound = true;
        break;
      }
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("Unrecognized option before input format type: {}", arg);
      return EXIT_FAILURE;
    }
    if (!formatFound) {
      SPDLOG_CRITICAL("Missing input format type");
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }

    // LCOV_EXCL_START
    bool explicitDesignFlags = false;
    std::vector<std::string>* currentDesign = nullptr;
    bool currentLiberty = false;
    std::vector<std::string> inputPaths;
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
    for (int i = parseStart; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "-v" || arg == "--verification") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing verification mode after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        std::string verificationError;
        if (!parseVerificationModeToken(argv[++i], verificationMode, verificationError)) {
          SPDLOG_CRITICAL("Invalid verification mode: {}", verificationError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        continue;
      }
      if (arg == "-k" || arg == "--max-k") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing max_k after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        std::string maxKError;
        if (!parseMaxKToken(argv[++i], secMaxK, maxKError)) {
          SPDLOG_CRITICAL("Invalid max_k: {}", maxKError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        secMaxKExplicit = true;
        continue;
      }
      if (arg == "--sec-engine") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing SEC engine after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        std::string secEngineError;
        if (!parseSecEngineToken(argv[++i], secEngine, secEngineError)) {
          SPDLOG_CRITICAL("Invalid SEC engine: {}", secEngineError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        secEngineExplicit = true;
        continue;
      }
      if (arg == "--sec-encoding") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing SEC encoding after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        std::string secEncodingError;
        if (!parseSecEncodingToken(argv[++i], secEncoding, secEncodingError)) {
          SPDLOG_CRITICAL("Invalid SEC encoding: {}", secEncodingError);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        secEncodingExplicit = true;
        continue;
      }
      if (arg == "--allow-boundary-mismatch") {
        allowBoundaryMismatch = true;
        continue;
      }
      // LCOV_EXCL_START
      if (arg == "--design1") {
        explicitDesignFlags = true;
        currentDesign = &designInputs.design0;
        currentLiberty = false;
        continue;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "--design2") {
        explicitDesignFlags = true;
        currentDesign = &designInputs.design1;
        currentLiberty = false;
        continue;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "--liberty" || arg == "--lib") {
        explicitDesignFlags = true;
        currentDesign = nullptr;
        currentLiberty = true;
        continue;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "--verilog_preprocessing") {
        verilogPreprocessing = true;
        continue;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "--compact") {
        compactMode = true;
        continue;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "--report-skipped-pos") {
        reportSkippedPOs = true;
        continue;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "--sv_design1_flist" || arg == "--sv_design2_flist" ||
          arg == "--verilog_design1_top" || arg == "--verilog_design2_top" ||
          arg == "--sv_design1_top" || arg == "--sv_design2_top") {
        if (i + 1 >= argc) {
          SPDLOG_CRITICAL("Missing value after {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        const std::string value = argv[++i];
        if (value.empty()) {
          SPDLOG_CRITICAL("Empty value provided for {}", arg);
          return EXIT_FAILURE;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        if (arg == "--sv_design1_flist") {
          systemVerilogOptions.design0.flist = value;
        } else if (arg == "--sv_design2_flist") {
          systemVerilogOptions.design1.flist = value;
        } else if (arg == "--sv_design1_top") {
          systemVerilogOptions.design0.top = value;
        } else if (arg == "--sv_design2_top") {
          systemVerilogOptions.design1.top = value;
        } else if (arg == "--verilog_design1_top") {
          verilogTopOptions.design0 = value;
        } else {
          verilogTopOptions.design1 = value;
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        continue;
      }
      // LCOV_EXCL_STOP

      // LCOV_EXCL_START
      if (!arg.empty() && arg[0] == '-') {
        SPDLOG_CRITICAL("Unrecognized option: {}", arg);
        return EXIT_FAILURE;
        // LCOV_EXCL_STOP
      }

      // LCOV_EXCL_START
      if (explicitDesignFlags) {
        if (currentLiberty) {
          libertyFiles.push_back(arg);
        } else if (currentDesign) {
          currentDesign->push_back(arg);
        } else {
        // LCOV_EXCL_STOP
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          SPDLOG_CRITICAL("Provide --design1 or --design2 before netlist paths");
          return EXIT_FAILURE;
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
      // LCOV_EXCL_START
      } else {
        inputPaths.emplace_back(arg);
        // LCOV_EXCL_STOP
      }
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
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
  // LCOV_EXCL_STOP

  const std::string runLogFilePath = configureMainLogger(
      logLevel, verificationMode == VerificationMode::SEC, logFileName);

  SPDLOG_INFO("KEPLER FORMAL: Run.");
  std::string inputFormatName = "VERILOG";
  std::string inputFormatToken = "verilog";
  if (inputFormatType == FormatType::NAJA_IF) {
    // LCOV_EXCL_START
    inputFormatName = "SNL";
    inputFormatToken = "naja_if";
  }
  // LCOV_EXCL_STOP
  if (inputFormatType == FormatType::SYSTEMVERILOG) {
    // LCOV_EXCL_START
    inputFormatName = "SYSTEMVERILOG";
    inputFormatToken = "systemverilog";
  }
  // LCOV_EXCL_STOP
  if (inputFormatType == FormatType::SV2V) {
    inputFormatName = "SV2V";
    inputFormatToken = "sv2v";
  }
  if (runResult != nullptr) {
    runResult->inputFormat = inputFormatToken;
    runResult->verification = verificationModeName(verificationMode);
    runResult->logFile = runLogFilePath;
  }
  SPDLOG_INFO("Input format: {}", inputFormatName);
  if (!runLogFilePath.empty()) {
    SPDLOG_INFO("Run log: {}", runLogFilePath);
  }
  logDesignPaths("Netlist 1", designInputs.design0);
  logDesignPaths("Netlist 2", designInputs.design1);

  // Basic validation
  if (inputFormatType == FormatType::SYSTEMVERILOG) {
    // LCOV_EXCL_START
    if (!hasSystemVerilogSources(designInputs.design0, systemVerilogOptions.design0) ||
        !hasSystemVerilogSources(designInputs.design1, systemVerilogOptions.design1)) {
      SPDLOG_CRITICAL(
      // LCOV_EXCL_STOP
          "Need SystemVerilog input sources for both designs (files and/or per-design flists)");
      // LCOV_EXCL_START
      print_usage(argv[0]);
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  } else if (inputFormatType == FormatType::SV2V) {
    if (!hasSystemVerilogSources(designInputs.design0, systemVerilogOptions.design0)) {
      SPDLOG_CRITICAL(
          "Need SystemVerilog input sources for design 1 (files and/or flist)");
      print_usage(argv[0]);
      return EXIT_FAILURE;
    }
    if (designInputs.design1.empty()) {
      SPDLOG_CRITICAL("Need Verilog input sources for design 2");
      print_usage(argv[0]);
      return EXIT_FAILURE;
    }
  } else if (designInputs.design0.empty() || designInputs.design1.empty()) {
    // LCOV_EXCL_START
    SPDLOG_CRITICAL("Need two input netlist paths (one per design)");
    print_usage(argv[0]);
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
  }
  if (inputFormatType == FormatType::NAJA_IF &&
      // LCOV_EXCL_START
      (designInputs.design0.size() != 1 || designInputs.design1.size() != 1)) {
    SPDLOG_CRITICAL("SNL input only supports one file per design");
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
  }
  if ((inputFormatType == FormatType::SYSTEMVERILOG ||
       inputFormatType == FormatType::SV2V) &&
      verificationMode != VerificationMode::SEC) {
    SPDLOG_CRITICAL(
        "SystemVerilog input formats require SEC verification (-v sec or verification: sec)");
    return EXIT_FAILURE;
  }
  if (verificationMode == VerificationMode::LEC && secMaxKExplicit) {
    // LCOV_EXCL_START
    SPDLOG_CRITICAL("max_k/-k is only supported with SEC verification");
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
  }
  if (verificationMode == VerificationMode::LEC && secEngineExplicit) {
    // LCOV_EXCL_START
    SPDLOG_CRITICAL("sec_engine/--sec-engine is only supported with SEC verification");
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
  }
  if (verificationMode == VerificationMode::LEC && secEncodingExplicit) {
    // LCOV_EXCL_START
    SPDLOG_CRITICAL("sec_encoding/--sec-encoding is only supported with SEC verification");  // LCOV_EXCL_LINE
    return EXIT_FAILURE;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (verificationMode == VerificationMode::SEC) {
    if (useScopes || cleanScopes) {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("SEC verification does not support scope extraction/cleaning");
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
    if (dumpCnf) {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL("SEC verification does not support CNF export");
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  }
  for (const auto& libraryFile : libertyFiles) {
    if (isPythonLoaderPath(libraryFile)) {
      // LCOV_EXCL_START
      SPDLOG_CRITICAL(
      // LCOV_EXCL_STOP
          "Python primitive loader {} must be provided through YAML config key "
          "py_tech_files, not liberty_files/--liberty inputs",
          libraryFile);
      // LCOV_EXCL_START
      return EXIT_FAILURE;
      // LCOV_EXCL_STOP
    }
  }
  if (inputFormatType != FormatType::SYSTEMVERILOG &&
      inputFormatType != FormatType::SV2V &&
      (systemVerilogOptions.design0.flist || systemVerilogOptions.design0.top ||
       systemVerilogOptions.design1.flist || systemVerilogOptions.design1.top)) {
    // LCOV_EXCL_START
    SPDLOG_CRITICAL(
        "SystemVerilog design options are only valid with -systemverilog/-sv/-sv2v input");
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
  }
  if (inputFormatType == FormatType::SV2V &&
      (systemVerilogOptions.design1.flist || systemVerilogOptions.design1.top)) {
    SPDLOG_CRITICAL(
        "sv2v format only accepts SystemVerilog options for design 1; "
        "design 2 is parsed as Verilog");
    return EXIT_FAILURE;
  }
  if (inputFormatType != FormatType::VERILOG &&
      inputFormatType != FormatType::SV2V &&
      (verilogTopOptions.design0 || verilogTopOptions.design1)) {
    SPDLOG_CRITICAL(
        "Verilog top options are only valid with -verilog/-sv2v input");
    return EXIT_FAILURE;
  }
  if (inputFormatType == FormatType::SV2V && verilogTopOptions.design0) {
    SPDLOG_CRITICAL(
        "sv2v format only accepts a Verilog top option for design 2; "
        "design 1 is parsed as SystemVerilog");
    return EXIT_FAILURE;
  }
  std::string svValidationError;
  if (!validateSystemVerilogOptions(systemVerilogOptions, svValidationError)) {
    // LCOV_EXCL_START
    // validateSystemVerilogOptions only fails for empty values, which the public
    // CLI/config parsing already rejects before reaching this point.
    // LCOV_DISABLED_START
    SPDLOG_CRITICAL("Invalid SystemVerilog options: {}", svValidationError);
    return EXIT_FAILURE;
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  }

  auto solverType = KEPLER_FORMAL::Config::getSolverType();
  KEPLER_FORMAL::Config::setReportSkippedPOs(reportSkippedPOs);
  const char* solverName =
      solverType == KEPLER_FORMAL::Config::SolverType::KISSAT
          ? "KISSAT"
          // LCOV_EXCL_START
          : (solverType == KEPLER_FORMAL::Config::SolverType::CADICAL
          // LCOV_EXCL_STOP
                 ? "CADICAL"
                 : "GLUCOSE");
  SPDLOG_INFO("Solver: {}", solverName);
  SPDLOG_INFO("Verification: {}", verificationModeName(verificationMode));
  if (verificationMode == VerificationMode::SEC) {
    SPDLOG_INFO("SEC max_k: {}", secMaxK);
    SPDLOG_INFO("SEC engine: {}", secEngineName(secEngine));
    SPDLOG_INFO("SEC encoding: {}", secEncodingName(secEncoding));
  }
  SPDLOG_INFO("Compact mode: {}", compactMode ? "enabled" : "disabled");
  SPDLOG_INFO("Skipped PO reports: {}", reportSkippedPOs ? "enabled" : "disabled");
  if (!libertyFiles.empty()) {
    for (const auto& lf : libertyFiles) SPDLOG_INFO("Library: {}", lf);
  }
  if (!pythonFiles.empty()) {
    if (runResult != nullptr) {
      runResult->reason =
          "py_tech_files are not supported by the in-process Python API";
      SPDLOG_CRITICAL("{}", runResult->reason);
      return EXIT_FAILURE;
    }
#ifdef KEPLER_FORMAL_NO_PY_TECH
    SPDLOG_CRITICAL(
        "py_tech_files are not available in this Kepler Formal build");
    return EXIT_FAILURE;
#else
    // LCOV_EXCL_START
    addNajaPythonPath(argv[0]);
    for (const auto& pf : pythonFiles) SPDLOG_INFO("Python library: {}", pf);
#endif
  }
  // LCOV_EXCL_STOP

  auto emitSecResult =
      [&](const KEPLER_FORMAL::SEC::SequentialEquivalenceResult& result) {
        if (runResult != nullptr) {
          runResult->bound = result.bound;
          runResult->reason = result.reason;
          runResult->coveredOutputs = result.coveredOutputs;
          runResult->totalOutputs = result.totalOutputs;
          runResult->skippedObservedOutputs = result.skippedObservedOutputs;
          if (result.proofProgress.has_value()) {
            runResult->provenOutputs = result.proofProgress->provenOutputs;
            runResult->unprovenOutputs.clear();
            runResult->unprovenOutputs.reserve(
                result.proofProgress->unprovenOutputs.size());
            for (const auto& output : result.proofProgress->unprovenOutputs) {
              runResult->unprovenOutputs.push_back(output.name);
            }
          } else if (
              result.status ==
                  KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Equivalent ||
              result.status == KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::PartiallyProved) {
            runResult->provenOutputs = result.coveredOutputs;
          }
        }
        // Naja creates its logger lazily and may replace spdlog's default
        // logger while loading SystemVerilog. Restore the run logger before
        // reporting the result so the requested SEC log remains complete.
        if (auto mainLogger = spdlog::get("kepler_formal_main_logger")) {
          spdlog::set_default_logger(mainLogger);
        }
        if (result.totalOutputs != 0) {
          SPDLOG_INFO(
              "SEC checked-output coverage: {:.2f}% ({}/{} covered/existing outputs).",
              result.outputCoveragePercent(),
              result.coveredOutputs,
              result.totalOutputs);
        }
        if (result.proofProgress.has_value()) {
          const auto& progress = *result.proofProgress;
          SPDLOG_INFO(
              "SEC {} proven outputs: {}/{}",
              progress.engineLabel,
              progress.provenOutputs,
              progress.totalOutputs);
          for (const auto& output : progress.unprovenOutputs) {
            SPDLOG_INFO(  // LCOV_EXCL_LINE
                "SEC {} not proven output[{}]={}",
                progress.engineLabel,
                output.index,
                output.name);
          }
        }
        if (!result.skippedObservedOutputs.empty()) {
          // LCOV_EXCL_START
          std::ostringstream skippedOutputs;
          for (const auto& skippedOutput : result.skippedObservedOutputs) {
            skippedOutputs << "  - " << skippedOutput << "\n";
            // LCOV_EXCL_STOP
          }
          // LCOV_EXCL_START
          SPDLOG_INFO(
          // LCOV_EXCL_STOP
              "SEC skipped observed outputs due to extraction or coverage "
              "limitations:\n{}",
              skippedOutputs.str());
        // LCOV_EXCL_START
        }
        // LCOV_EXCL_STOP
        // LCOV_EXCL_STOP
        if (reportSkippedPOs) {
          // LCOV_EXCL_START
          writeBoundaryTermsReport(
              kBoundaryTermsReport, result.extractedBoundaryReports);
          writeResetUnanchoredSkippedOutputsReport(
              kSkippedResetUnanchoredPOReport,
              result.resetUnanchoredSkippedOutputs);
          writeMultiClockDomainSkippedOutputsReport(
              kSkippedMultiClockDomainPOReport,
              result.multiClockDomainSkippedOutputs);
          writeOpaqueCellSkippedOutputsReport(
              kSkippedOpaqueCellPOReport,
              result.opaqueCellSkippedOutputs);
        }
        // LCOV_EXCL_STOP
        switch (result.status) {
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Equivalent:
            if (runResult != nullptr) {
              runResult->status = KEPLER_FORMAL::RunStatus::Equivalent;
            }
            if (secEncoding ==
                KEPLER_FORMAL::SEC::SecEncoding::DualRailSteady) {
              SPDLOG_INFO(
                  "No binary-defined difference was found. SEC proved "
                  "equivalence under the dual-rail steady-state abstraction "
                  "at k = {}.",
                  result.bound);
            } else {
              SPDLOG_INFO(
                  "No difference was found. SEC proved equivalence at k = {}.",
                  result.bound);
            }
            return kSecProvedExitCode;
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::PartiallyProved: {
            if (runResult != nullptr) {
              runResult->status = KEPLER_FORMAL::RunStatus::PartiallyProved;
            }
            const size_t provedOutputs = result.proofProgress.has_value()
                                             ? result.proofProgress->provenOutputs
                                             : result.coveredOutputs;
            const size_t totalOutputs = result.totalOutputs;
            SPDLOG_INFO(
                "SEC partially proved equivalence at k = {}: {}/{} outputs "
                "proved; remaining outputs are inconclusive.",
                result.bound,
                provedOutputs,
                totalOutputs);
            SPDLOG_WARN(
                "SEC verification did not prove all observed outputs.");
            if (!result.reason.empty()) {
              SPDLOG_INFO("SEC partial-proof details: {}", result.reason);
            }
            return kSecPartiallyProvedExitCode;
          }
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Different:
            if (runResult != nullptr) {
              runResult->status = KEPLER_FORMAL::RunStatus::Different;
            }
            // LCOV_EXCL_START
            SPDLOG_INFO(
            // LCOV_EXCL_STOP
                "Difference was found. SEC found a counterexample at k = {}.",
                result.bound);
            // LCOV_EXCL_START
            if (!result.reason.empty()) {
              SPDLOG_INFO("SEC counterexample details:\n{}", result.reason);
            }
            return kSecCounterexampleExitCode;
            // LCOV_EXCL_STOP
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Inconclusive:
            if (runResult != nullptr) {
              runResult->status = KEPLER_FORMAL::RunStatus::Inconclusive;
            }
            // LCOV_EXCL_START
            if (secInconclusiveStoppedBeforeMaxK(result.reason)) {
              SPDLOG_INFO(
                  "SEC was inconclusive before completing max_k = {}: {}",
                  secMaxK,
                  result.reason);
            } else {
              SPDLOG_INFO(
                  "SEC was inconclusive up to max_k = {}: {}",
                  secMaxK,
                  result.reason);
            }
            SPDLOG_WARN(
                "SEC verification did not produce a proof or counterexample.");
            return kSecInconclusiveExitCode;
            // LCOV_EXCL_STOP
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Unsupported:
            if (runResult != nullptr) {
              runResult->status = KEPLER_FORMAL::RunStatus::Unsupported;
            }
          // LCOV_DISABLED_STOP
          default:
            // LCOV_EXCL_START
            SPDLOG_CRITICAL(
            // LCOV_EXCL_STOP
                "SEC cannot run on this design pair: {}", result.reason);
            // LCOV_EXCL_START
            return kSecInconclusiveExitCode;
            // LCOV_EXCL_STOP
        }
      };

  // --------------------------------------------------------------------------
  // 2. Load two netlists via Cap’n Proto (or via VRL constructor)
  // --------------------------------------------------------------------------
  naja::NL::SNLDesign* top0 = nullptr;
  naja::NL::SNLDesign* top1 = nullptr;
  try {
    NLUniverse::create();
    NLDB* db0 = nullptr;
    bool primitivesAreLoaded = false;

    auto loadLibraries = [&](NLDB* db,
                             std::vector<NLLibrary*>* primitiveLibraries = nullptr) -> bool {
      if (libertyFiles.empty() && pythonFiles.empty()) {
        // LCOV_EXCL_START
        // LCOV_DISABLED_START
        return false;
        // LCOV_DISABLED_STOP
        // LCOV_EXCL_STOP
      }
      auto primitivesLibrary =
          NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
      if (primitiveLibraries != nullptr) {
        primitiveLibraries->push_back(primitivesLibrary);
      }
      for (const auto& libraryFile : libertyFiles) {
        std::filesystem::path libraryPath(libraryFile);
        SPDLOG_INFO("Loading library file: {}", libraryFile);
        SNLLibertyConstructor constructor(primitivesLibrary);
        constructor.construct(libraryPath);
      }
#ifndef KEPLER_FORMAL_NO_PY_TECH
      for (const auto& pythonFile : pythonFiles) {
        // LCOV_EXCL_START
        std::filesystem::path pythonPath(pythonFile);
        SPDLOG_INFO("Loading python primitive file: {}", pythonFile);
        SNLPyLoader::loadPrimitives(primitivesLibrary, pythonPath);
      }
#endif
      // LCOV_EXCL_STOP
      return true;
    };

    const auto isHdlFormat = [&]() {
      return inputFormatType == FormatType::VERILOG ||
             inputFormatType == FormatType::SYSTEMVERILOG ||
             inputFormatType == FormatType::SV2V;
    };

    const auto designUsesSystemVerilog = [&](int designIndex) {
      return inputFormatType == FormatType::SYSTEMVERILOG ||
             (inputFormatType == FormatType::SV2V && designIndex == 0);
    };

    auto loadOneDesign = [&](const std::vector<std::string>& designPaths,
                             const SystemVerilogDesignOptions& designOptions,
                             int designIndex,
                             int dbID) -> NLDB* {
      NLDB* db = nullptr;
      bool primitivesLoadedForDesign = false;
      std::vector<NLLibrary*> primitiveLibraries;

      if (!libertyFiles.empty() || !pythonFiles.empty()) {
        db = NLDB::create(NLUniverse::get());
        primitivesLoadedForDesign = loadLibraries(db, &primitiveLibraries);
        if (!primitivesLoadedForDesign) {
          // LCOV_EXCL_START
          throw std::runtime_error("Failed to load library files");
          // LCOV_EXCL_STOP
        }
      }

      if (isHdlFormat()) {
        if (!db) {
          // LCOV_EXCL_START
          db = NLDB::create(NLUniverse::get());
        }
        // LCOV_EXCL_STOP
        db->setID(dbID);
        const bool useSystemVerilog = designUsesSystemVerilog(designIndex);
        SPDLOG_INFO("Parsing {} file(s) for design {}",
                    useSystemVerilog ? "systemverilog" : "verilog",
                    designIndex + 1);
        auto designLibrary = NLLibrary::create(db, NLName("DESIGN"));
        if (useSystemVerilog) {
          // LCOV_EXCL_START
          SNLSVConstructor constructor(designLibrary);
          std::vector<std::filesystem::path> temporaryFiles;
          // LCOV_EXCL_STOP
          const auto* sv2vPrimitiveLibraries =
              (inputFormatType == FormatType::SV2V && designIndex == 0)
                  ? &primitiveLibraries
                  : nullptr;
          const auto svInputPaths =
              // LCOV_EXCL_START
              buildSystemVerilogInputPaths(
                  designPaths, designOptions, temporaryFiles, sv2vPrimitiveLibraries);
              // LCOV_EXCL_STOP
          try {
            // LCOV_EXCL_START
            constructor.construct(svInputPaths);
            if (sv2vPrimitiveLibraries != nullptr) {
              reconnectGeneratedPrimitiveStubs(designLibrary, *sv2vPrimitiveLibraries);
            }
            // LCOV_EXCL_STOP
            // LCOV_EXCL_START
            // LCOV_DISABLED_START
          } catch (...) {
            for (const auto& temporaryFile : temporaryFiles) {
              std::error_code ec;
              std::filesystem::remove(temporaryFile, ec);
              // LCOV_DISABLED_STOP
            }
            // LCOV_DISABLED_START
            throw;
          }
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
          // LCOV_EXCL_START
          for (const auto& temporaryFile : temporaryFiles) {
            std::error_code ec;
            std::filesystem::remove(temporaryFile, ec);
            // LCOV_EXCL_STOP
          }
        // LCOV_EXCL_START
        } else {
        // LCOV_EXCL_STOP
          SNLVRLConstructor constructor(designLibrary);
          constructor.config_.preprocessEnabled_ = verilogPreprocessing;
          constructor.construct(toPathVector(designPaths));
        }
        auto top = useSystemVerilog
                       ? SNLUtils::findTop(designLibrary)
                       : selectTopDesign(
                             designLibrary,
                             designIndex == 0 ? verilogTopOptions.design0
                                              : verilogTopOptions.design1,
                             designIndex);
        if (!top) {
            // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("No top design was found after parsing input");
          // LCOV_DISABLED_STOP
            // LCOV_EXCL_STOP
        }
        db->setTopDesign(top);
        SPDLOG_INFO("Found top design: {}", top->getString());
      } else {
        // LCOV_EXCL_START
        SPDLOG_INFO("Loading Naja IF: {}", designPaths[0]);
        naja::NL::SNLCapnP::LoadingConfiguration config;
        config.primitiveConflictPolicy_ =
            primitivesLoadedForDesign
            // LCOV_EXCL_STOP
                ? naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::PreferExisting
                : naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::ForbidConflicts;
	        // LCOV_EXCL_START
	        db = SNLCapnP::load(designPaths[0].c_str(), config);
	        if (!db) {
	        // LCOV_EXCL_STOP
            // LCOV_EXCL_START
	          // LCOV_DISABLED_START
	          throw std::runtime_error("Failed to load Naja IF: " + designPaths[0]);
	          // LCOV_DISABLED_STOP
            // LCOV_EXCL_STOP
	        }
        // LCOV_EXCL_START
        db->setID(dbID);
        // LCOV_EXCL_STOP
      }

	      if (!db->getTopDesign()) {
          // LCOV_EXCL_START
	        // LCOV_DISABLED_START
	        throw std::runtime_error("Top design not set for loaded netlist");
	        // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
	      }
      return db;
    // LCOV_EXCL_START
    };
    // LCOV_EXCL_STOP

    auto buildCompactSnapshotForTop =
        // LCOV_EXCL_START
        [&](naja::NL::SNLDesign* top,
        // LCOV_EXCL_STOP
            const char* designLabel) {
          // LCOV_EXCL_START
          KEPLER_FORMAL::BuildPrimaryOutputClauses builder;
          NLUniverse::get()->setTopDesign(top);
          naja::DNL::destroy();
          builder.collect();
          SPDLOG_INFO("Collected {} PIs for {}", builder.getInputs().size(), designLabel);
          SPDLOG_INFO("Collected {} POs for {}", builder.getOutputs().size(), designLabel);
          std::vector<KEPLER_FORMAL::BuildPrimaryOutputClauses::PathKey>
              boundaryInputs;
          if (!allowBoundaryMismatch) {
            boundaryInputs = builder.getLecBoundaryInputs();
          }
          auto inputs = builder.getInputs();
          auto outputs = builder.getOutputs();
          builder.setInputs(inputs);
          builder.setOutputs(outputs);
          builder.build();
          return captureCompactSnapshot(builder, std::move(boundaryInputs));
        };
        // LCOV_EXCL_STOP

    if (compactMode && verificationMode == VerificationMode::LEC && !useScopes) {
      // LCOV_EXCL_START
      NLDB* compactDb0 =
          loadOneDesign(designInputs.design0, systemVerilogOptions.design0, 0, 2);
      top0 = compactDb0->getTopDesign();
      auto snapshot0 = buildCompactSnapshotForTop(top0, "design 0");
      naja::DNL::destroy();
      compactDb0->destroy();
      top0 = nullptr;
      // LCOV_EXCL_STOP

      // LCOV_EXCL_START
      NLDB* compactDb1 =
          loadOneDesign(designInputs.design1, systemVerilogOptions.design1, 1, 1);
      top1 = compactDb1->getTopDesign();
      auto snapshot1 = buildCompactSnapshotForTop(top1, "design 1");
      naja::DNL::destroy();
      compactDb1->destroy();
      top1 = nullptr;
      // LCOV_EXCL_STOP

      try {
        // LCOV_EXCL_START
        KEPLER_FORMAL::MiterStrategy MiterS(nullptr, nullptr, logFileName);
        MiterS.setAllowBoundaryMismatch(allowBoundaryMismatch);
        if (dumpCnf) {
          const std::string outPath = dumpCnfPath.empty() ? "miter.cnf" : dumpCnfPath;
          MiterS.setCnfDump(true, outPath);
        }
        if (dumpPoCnf) {
          const std::string outPath = dumpPoCnfPath.empty() ? "po_cnfs" : dumpPoCnfPath;
          MiterS.setPoCnfDump(true, outPath);
        }
        const bool equivalent = MiterS.runCompactSnapshots(snapshot0, snapshot1);
        if (runResult != nullptr) {
          runResult->logFile =
              KEPLER_FORMAL::MiterStrategy::getActualLogFileName();
          runResult->status = equivalent ? KEPLER_FORMAL::RunStatus::Equivalent
                                         : KEPLER_FORMAL::RunStatus::Different;
        }
        if (equivalent) {
          SPDLOG_INFO("No difference was found.");
        } else {
          SPDLOG_INFO("Difference was found. Please refer to the log(miter_log_x.txt) for details.");
          // LCOV_EXCL_STOP
        }
	      // LCOV_EXCL_START
	      // LCOV_DISABLED_START
	      } catch (const std::exception& e) {
	        SPDLOG_ERROR("Workflow failed: {}", e.what());
	        return EXIT_FAILURE;
      }
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
      // LCOV_EXCL_START
      return EXIT_SUCCESS;
    }
    // LCOV_EXCL_STOP

    if (compactMode && verificationMode == VerificationMode::SEC) {
      auto releaseCompactDb = [](NLDB* db) {
        naja::DNL::destroy();
        if (auto* universe = NLUniverse::get()) {
          // Sequential extraction restores the current top when possible. In
          // compact mode we deliberately drop that elaborated DB before loading
          // the next design, so clear the universe top DB first to avoid a
          // dangling top pointer.
          universe->setTopDB(nullptr);
        }
        if (db != nullptr) {
          db->destroy();
        }
      };

      auto extractCompactSecModel =
          [&](const std::vector<std::string>& designPaths,
              const SystemVerilogDesignOptions& designOptions,
              int designIndex,
              int dbID,
              const char* designLabel) {
            NLDB* db =
                loadOneDesign(designPaths, designOptions, designIndex, dbID);
            try {
              auto* top = db->getTopDesign();
              if (top == nullptr) {
                // LCOV_EXCL_START
                // LCOV_DISABLED_START
                throw std::runtime_error(
                    std::string("Top design not set for ") + designLabel);
                    // LCOV_DISABLED_STOP
                // LCOV_EXCL_STOP
              }
              auto model = KEPLER_FORMAL::SEC::SequentialDesignModel::extract(top);
              releaseCompactDb(db);
              return model;
            } catch (...) {
              // LCOV_EXCL_START
              // LCOV_DISABLED_START
              releaseCompactDb(db);
              throw;
              // LCOV_DISABLED_STOP
              // LCOV_EXCL_STOP
            // LCOV_EXCL_START
            }  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
          };

      try {
        SPDLOG_INFO(
            "SEC compact mode: extracting and releasing design 1 before "
            "loading design 2");
        const auto model0 = extractCompactSecModel(
            designInputs.design0,
            systemVerilogOptions.design0,
            0,
            2,
            "design 1");
        if (inputFormatType != FormatType::SV2V &&
            sameCompactSecDesignSpec(
                inputFormatType == FormatType::SYSTEMVERILOG,
                designInputs,
                systemVerilogOptions,
                verilogTopOptions)) {
          // CVA6-style smoke runs often compare a design against itself. In
          // compact SEC, extracting that identical second side would require
          // holding the already extracted value model while elaborating the
          // same large SystemVerilog design again. Reusing the immutable model
          // keeps the formal problem identical and removes that avoidable peak.
          // LCOV_EXCL_START
          SPDLOG_INFO(
          // LCOV_EXCL_STOP
              "SEC compact mode: reusing extracted design 1 model for "
              "identical design 2 input");
          // LCOV_EXCL_START
          KEPLER_FORMAL::SEC::SequentialEquivalenceStrategy strategy(
              nullptr,
              nullptr,
              solverType,
              secEngine,
              secEncoding);
          return emitSecResult(
              strategy.runExtractedModels(model0, model0, secMaxK));
              // LCOV_EXCL_STOP
        }
        SPDLOG_INFO(
            "SEC compact mode: extracting and releasing design 2 before "
            "starting proof");
        const auto model1 = extractCompactSecModel(
            designInputs.design1,
            systemVerilogOptions.design1,
            1,
            1,
            "design 2");

        KEPLER_FORMAL::SEC::SequentialEquivalenceStrategy strategy(
            nullptr,
            nullptr,
            solverType,
            secEngine,
            secEncoding);
        return emitSecResult(
            strategy.runExtractedModels(model0, model1, secMaxK));
      // LCOV_EXCL_START
      } catch (const std::exception& e) {
        // LCOV_DISABLED_START
        SPDLOG_ERROR("SEC compact workflow failed: {}", e.what());
        return EXIT_FAILURE;
      }
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
    std::vector<NLLibrary*> db0PrimitiveLibraries;
    if (!libertyFiles.empty() || !pythonFiles.empty()) {
      db0 = NLDB::create(NLUniverse::get());
      primitivesAreLoaded = loadLibraries(db0, &db0PrimitiveLibraries);
      if (!primitivesAreLoaded) {
        return EXIT_FAILURE;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
    if (isHdlFormat()) {
      if (!db0) {
        db0 = NLDB::create(NLUniverse::get());
      }
      const auto design0Paths = toPathVector(designInputs.design0);
      const bool design0UsesSystemVerilog = designUsesSystemVerilog(0);
      SPDLOG_INFO("Parsing {} file(s) for design 1",
      // LCOV_EXCL_STOP
                  design0UsesSystemVerilog ? "systemverilog" : "verilog");
      // LCOV_EXCL_START
      auto designLibrary = NLLibrary::create(db0, NLName("DESIGN"));
      if (design0UsesSystemVerilog) {
        SNLSVConstructor constructor(designLibrary);
        std::vector<std::filesystem::path> temporaryFiles;
        const auto* sv2vPrimitiveLibraries =
            inputFormatType == FormatType::SV2V ? &db0PrimitiveLibraries : nullptr;
        const auto svInputPaths = buildSystemVerilogInputPaths(
            designInputs.design0,
            systemVerilogOptions.design0,
            temporaryFiles,
            sv2vPrimitiveLibraries);
            // LCOV_EXCL_STOP
        try {
          // LCOV_EXCL_START
          constructor.construct(svInputPaths);
          if (sv2vPrimitiveLibraries != nullptr) {
            reconnectGeneratedPrimitiveStubs(designLibrary, *sv2vPrimitiveLibraries);
          }
        } catch (...) {
          for (const auto& temporaryFile : temporaryFiles) {
            std::error_code ec;
            std::filesystem::remove(temporaryFile, ec);
            // LCOV_EXCL_STOP
          }
          // LCOV_EXCL_START
          throw;
        }
        for (const auto& temporaryFile : temporaryFiles) {
          std::error_code ec;
          std::filesystem::remove(temporaryFile, ec);
          // LCOV_EXCL_STOP
        }
      // LCOV_EXCL_START
      } else {
        SNLVRLConstructor constructor(designLibrary);
        constructor.config_.preprocessEnabled_ = verilogPreprocessing;
        constructor.construct(design0Paths);
      }
      auto top = design0UsesSystemVerilog
                     ? SNLUtils::findTop(designLibrary)
                     : selectTopDesign(designLibrary, verilogTopOptions.design0, 0);
      if (top) {
        db0->setTopDesign(top);
        SPDLOG_INFO("Found top design: {}", top->getString());
      } else {
        SPDLOG_CRITICAL("No top design was found after parsing input");
        return EXIT_FAILURE;
        // LCOV_EXCL_STOP
      }
    // LCOV_EXCL_START
    } else {
      SPDLOG_INFO("Loading Naja IF: {}", designInputs.design0[0]);
      naja::NL::SNLCapnP::LoadingConfiguration config;
      config.primitiveConflictPolicy_ =
          primitivesAreLoaded
          // LCOV_EXCL_STOP
              ? naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::PreferExisting
              : naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::ForbidConflicts;
      // LCOV_EXCL_START
      db0 = SNLCapnP::load(designInputs.design0[0].c_str(), config);
      if (!db0) {
      // LCOV_EXCL_STOP
        // LCOV_EXCL_START
        // LCOV_DISABLED_START
        SPDLOG_CRITICAL("Failed to load Naja IF: {}", designInputs.design0[0]);
        return EXIT_FAILURE;
        // LCOV_DISABLED_STOP
        // LCOV_EXCL_STOP
      }
    }

    // get db0 top
    // LCOV_EXCL_START
    top0 = db0->getTopDesign();
    if (!top0) {
    // LCOV_EXCL_STOP
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      SPDLOG_CRITICAL("Top design not set for first netlist");
      return EXIT_FAILURE;
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    db0->setID(2);  // Increment ID to avoid conflicts
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
    NLDB* db1 = nullptr;
    // LCOV_EXCL_STOP

    // Prepare second DB and primitives if needed
    // LCOV_EXCL_START
    if (!libertyFiles.empty() || !pythonFiles.empty()) {
      db1 = NLDB::create(NLUniverse::get());
      db1->setID(1);
      if (!loadLibraries(db1)) {
      // LCOV_EXCL_STOP
        // LCOV_EXCL_START
        // LCOV_DISABLED_START
        return EXIT_FAILURE;
        // LCOV_DISABLED_STOP
        // LCOV_EXCL_STOP
      }
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
    if (isHdlFormat()) {
      if (!db1) {
        db1 = NLDB::create(NLUniverse::get());
      }
      const auto design1Paths = toPathVector(designInputs.design1);
      const bool design1UsesSystemVerilog = designUsesSystemVerilog(1);
      SPDLOG_INFO("Parsing {} file(s) for design 2",
      // LCOV_EXCL_STOP
                  design1UsesSystemVerilog ? "systemverilog" : "verilog");
      // LCOV_EXCL_START
      auto designLibrary = NLLibrary::create(db1, NLName("DESIGN"));
      if (design1UsesSystemVerilog) {
        SNLSVConstructor constructor(designLibrary);
        std::vector<std::filesystem::path> temporaryFiles;
        const auto svInputPaths = buildSystemVerilogInputPaths(
            designInputs.design1, systemVerilogOptions.design1, temporaryFiles);
            // LCOV_EXCL_STOP
        try {
          // LCOV_EXCL_START
          constructor.construct(svInputPaths);
        } catch (...) {
          for (const auto& temporaryFile : temporaryFiles) {
            std::error_code ec;
            std::filesystem::remove(temporaryFile, ec);
            // LCOV_EXCL_STOP
          }
          // LCOV_EXCL_START
          throw;
        }
        for (const auto& temporaryFile : temporaryFiles) {
          std::error_code ec;
          std::filesystem::remove(temporaryFile, ec);
          // LCOV_EXCL_STOP
        }
      // LCOV_EXCL_START
      } else {
        SNLVRLConstructor constructor(designLibrary);
        constructor.config_.preprocessEnabled_ = verilogPreprocessing;
        constructor.construct(design1Paths);
      }
      auto top = design1UsesSystemVerilog
                     ? SNLUtils::findTop(designLibrary)
                     : selectTopDesign(designLibrary, verilogTopOptions.design1, 1);
      if (top) {
        db1->setTopDesign(top);
        SPDLOG_INFO("Found top design: {}", top->getString());
      } else {
        SPDLOG_CRITICAL("No top design was found after parsing input");
        return EXIT_FAILURE;
        // LCOV_EXCL_STOP
      }
    // LCOV_EXCL_START
    } else {
      SPDLOG_INFO("Loading Naja IF: {}", designInputs.design1[0]);
      naja::NL::SNLCapnP::LoadingConfiguration config;
      config.primitiveConflictPolicy_ =
          primitivesAreLoaded
          // LCOV_EXCL_STOP
              ? naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::PreferExisting
              : naja::NL::SNLCapnP::LoadingConfiguration::PrimitiveConflictPolicy::ForbidConflicts;
      // LCOV_EXCL_START
      db1 = SNLCapnP::load(designInputs.design1[0].c_str(), config);
      if (!db1) {
      // LCOV_EXCL_STOP
        // LCOV_EXCL_START
        // LCOV_DISABLED_START
        SPDLOG_CRITICAL("Failed to load Naja IF: {}", designInputs.design1[0]);
        return EXIT_FAILURE;
        // LCOV_DISABLED_STOP
        // LCOV_EXCL_STOP
      }
    }

    // get db1 top
    // LCOV_EXCL_START
    top1 = db1->getTopDesign();
    if (!top1) {
    // LCOV_EXCL_STOP
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      SPDLOG_CRITICAL("Top design not set for second netlist");
      return EXIT_FAILURE;
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    }
  // LCOV_EXCL_START
  } catch (const std::exception& e) {
    SPDLOG_CRITICAL("Netlist loading failed: {}", e.what());
    return EXIT_FAILURE;
  }
  // LCOV_EXCL_STOP

  // --------------------------------------------------------------------------
  // 4. Hand off to the rest of the editing/analysis workflow
  // --------------------------------------------------------------------------
  // LCOV_EXCL_START
  if (verificationMode == VerificationMode::SEC) {
  // LCOV_EXCL_STOP
    try {
      // LCOV_EXCL_START
      KEPLER_FORMAL::SEC::SequentialEquivalenceStrategy strategy(
          top0,
          top1,
          solverType,
          secEngine,
          secEncoding);
      return emitSecResult(strategy.run(secMaxK));
      // LCOV_EXCL_STOP
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    } catch (const std::exception& e) {
      SPDLOG_ERROR("SEC workflow failed: {}", e.what());
      return EXIT_FAILURE;
    }
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  // LCOV_EXCL_START
  } else if (inputFormatType == FormatType::NAJA_IF && useScopes) {
    KEPLER_FORMAL::MiterStrategy MiterS(top0, top1);
    MiterS.setAllowBoundaryMismatch(true);
    MiterS.init(false);
    ScopeExtraction extractor(top0, top1);
    extractor.collectVerificationScopes();
    if (cleanScopes) {
      extractor.cleanVerificationScopes(MiterS.getPIs0(), MiterS.getPIs1());
    }
    bool allScopesEquivalent = true;
    for (auto scopes : extractor.getScopesToVerify()) {
      SPDLOG_INFO("Looking at scope: {} ",
      // LCOV_EXCL_STOP
                  scopes.first->getName().getString());
      // std::string scopeLogFile = (logFileName.empty() ? "kf_" : logFileName) + "_" +
      //                           scopes.first->getName().getString() + ".txt";
      try {
        // LCOV_EXCL_START
        KEPLER_FORMAL::MiterStrategy MiterScope(scopes.first, scopes.second, logFileName);
        MiterScope.setAllowBoundaryMismatch(allowBoundaryMismatch);
        if (dumpCnf) {
          std::string scopeName = sanitizeFileToken(scopes.first->getName().getString());
          std::string outPath = dumpCnfPath.empty()
                                    ? ("miter_" + scopeName + ".cnf")
                                    : dumpCnfPath;
          MiterScope.setCnfDump(true, outPath);
        }
        if (dumpPoCnf) {
          std::string scopeName = sanitizeFileToken(scopes.first->getName().getString());
          std::string outPath = dumpPoCnfPath.empty()
                                    ? ("po_cnfs_" + scopeName)
                                    : dumpPoCnfPath;
          MiterScope.setPoCnfDump(true, outPath);
        }
        MiterScope.init();
        if (MiterScope.run(compactMode)) {
          SPDLOG_INFO("No difference was found for scope: {} , {}",
          // LCOV_EXCL_STOP
                      scopes.first->getName().getString(),
                      scopes.second->getName().getString());
        // LCOV_EXCL_START
        } else {
          allScopesEquivalent = false;
          SPDLOG_INFO("Difference was found for scope: {} , {}. Please refer to the log(miter_log_x.txt) for details.",  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
                      scopes.first->getName().getString(),
                      scopes.second->getName().getString());
        }
      // LCOV_EXCL_START
      } catch (const std::exception& e) {
      // LCOV_EXCL_STOP
        // LCOV_EXCL_START
        // LCOV_DISABLED_START
        SPDLOG_ERROR("Workflow failed for scope: {} , {}: {}", 
        // LCOV_DISABLED_STOP
                      scopes.first->getName().getString(),
                      scopes.second->getName().getString(),
                      e.what());
        // LCOV_DISABLED_START
        return EXIT_FAILURE;
      }
      // LCOV_DISABLED_STOP
        // LCOV_EXCL_STOP
    }
    if (runResult != nullptr) {
      runResult->logFile =
          KEPLER_FORMAL::MiterStrategy::getActualLogFileName();
      runResult->status = allScopesEquivalent
                              ? KEPLER_FORMAL::RunStatus::Equivalent
                              : KEPLER_FORMAL::RunStatus::Different;
    }
  // LCOV_EXCL_START
  } else {
  // LCOV_EXCL_STOP
    try {
      // LCOV_EXCL_START
      KEPLER_FORMAL::MiterStrategy MiterS(top0, top1, logFileName);
      MiterS.setAllowBoundaryMismatch(allowBoundaryMismatch);
      if (dumpCnf) {
        const std::string outPath = dumpCnfPath.empty() ? "miter.cnf" : dumpCnfPath;
        MiterS.setCnfDump(true, outPath);
      }
      if (dumpPoCnf) {
        const std::string outPath = dumpPoCnfPath.empty() ? "po_cnfs" : dumpPoCnfPath;
        MiterS.setPoCnfDump(true, outPath);
      }
      MiterS.init();
      const bool equivalent = MiterS.run(compactMode);
      if (runResult != nullptr) {
        runResult->logFile =
            KEPLER_FORMAL::MiterStrategy::getActualLogFileName();
        runResult->status = equivalent ? KEPLER_FORMAL::RunStatus::Equivalent
                                       : KEPLER_FORMAL::RunStatus::Different;
      }
      if (equivalent) {
        SPDLOG_INFO("No difference was found.");
      } else {
        SPDLOG_INFO("Difference was found. Please refer to the log(miter_log_x.txt) for details.");
        // LCOV_EXCL_STOP
      }
    // LCOV_EXCL_START
    } catch (const std::exception& e) {
    // LCOV_EXCL_STOP
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      SPDLOG_ERROR("Workflow failed: {}", e.what());
      return EXIT_FAILURE;
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  return EXIT_SUCCESS;
  // LCOV_EXCL_STOP
}

namespace KEPLER_FORMAL {

const char* runStatusName(RunStatus status) {
  switch (status) {
    case RunStatus::NoResult:
      return "no_result";
    case RunStatus::Equivalent:
      return "equivalent";
    case RunStatus::Different:
      return "different";
    case RunStatus::PartiallyProved:
      return "partially_proved";
    case RunStatus::Inconclusive:
      return "inconclusive";
    case RunStatus::Unsupported:
      return "unsupported";
    case RunStatus::Error:
    default:
      return "error";
  }
}

void cleanupKeplerFormalState() {
  MiterStrategy::cleanupProcessState();
  Tree2BoolExpr::iso2boolExpr_.clear();
  BoolExprCache::destroy();
}

int runKeplerFormal(int argc, char** argv, RunResult& result) {
  static std::mutex runMutex;
  static thread_local bool runInProgress = false;
  if (runInProgress) {
    throw std::runtime_error("Kepler Formal in-process calls are not reentrant");
  }
  struct ReentrancyGuard {
    bool& inProgress;
    explicit ReentrancyGuard(bool& value) : inProgress(value) { inProgress = true; }
    ~ReentrancyGuard() { inProgress = false; }
  } reentrancyGuard{runInProgress};
  std::lock_guard<std::mutex> runLock(runMutex);

  if (NLUniverse::get() != nullptr) {
    throw std::runtime_error(
        "Kepler Formal cannot start while another Naja universe is active");
  }

  const auto previousSolver = Config::getSolverType();
  const bool previousReportSkippedPOs = Config::getReportSkippedPOs();
  const auto previousDefaultLogger = spdlog::default_logger();
  const auto previousNamedLogger = spdlog::get("kepler_formal_main_logger");
  const auto previousMiterLogger = spdlog::get("miter_logger");
  const auto previousMiterFallbackLogger =
      spdlog::get("miter_logger_fallback");
  struct RunStateGuard {
    Config::SolverType solver;
    bool reportSkippedPOs;
    std::shared_ptr<spdlog::logger> defaultLogger;
    std::shared_ptr<spdlog::logger> namedLogger;
    std::shared_ptr<spdlog::logger> miterLogger;
    std::shared_ptr<spdlog::logger> miterFallbackLogger;
    ~RunStateGuard() {
      cleanupKeplerFormalState();
      Config::setSolverType(solver);
      Config::setReportSkippedPOs(reportSkippedPOs);
      const auto restoreNamedLogger = [](
                                          const char* name,
                                          const std::shared_ptr<spdlog::logger>&
                                              previousLogger) {
        if (spdlog::get(name) != previousLogger) {
          spdlog::drop(name);
          if (previousLogger != nullptr) {
            spdlog::register_logger(previousLogger);
          }
        }
      };
      restoreNamedLogger("kepler_formal_main_logger", namedLogger);
      restoreNamedLogger("miter_logger", miterLogger);
      restoreNamedLogger("miter_logger_fallback", miterFallbackLogger);
      spdlog::set_default_logger(defaultLogger);
    }
  } stateGuard{
      previousSolver,
      previousReportSkippedPOs,
      previousDefaultLogger,
      previousNamedLogger,
      previousMiterLogger,
      previousMiterFallbackLogger};

  Config::setSolverType(Config::SolverType::KISSAT);
  Config::setReportSkippedPOs(false);
  const int rc = KeplerFormalMainImpl(argc, argv, &result);
  result.exitCode = rc;
  if (rc != EXIT_SUCCESS && result.status == RunStatus::Error &&
      result.reason.empty()) {
    result.reason =
        "Kepler Formal failed before producing a verification result; "
        "see the native log output for details";
  }
  return rc;
}

}  // namespace KEPLER_FORMAL

int KeplerFormalMain(int argc, char** argv) {
  return KeplerFormalMainImpl(argc, argv, nullptr);
}
