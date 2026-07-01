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
#include "CFrontend.h"
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

// LCOV_EXCL_START
static void print_usage(const char* prog) {
  SPDLOG_INFO(
  // LCOV_EXCL_STOP
      "Usage: {} [--config <file>] | <-naja_if/-verilog/-systemverilog/-sv> "
      "[-v <lec|sec>] [-k <max-k>] [--sec-engine <legacy|k_induction|imc|pdr>] [--sec-encoding <binary|dual_rail_steady>] <netlist1> <netlist2> [<library-file>...] | "
      "<-naja_if/-verilog/-systemverilog/-sv> --design1 <file...> --design2 "
      "<file...> [--liberty <library-file>...] [-v <lec|sec>] [-k <max-k>] [--sec-engine <legacy|k_induction|imc|pdr>] [--sec-encoding <binary|dual_rail_steady>] "
      "[--sec-internal-state-correspondence] [--no-sec-uncomputable-seq-boundary] [--compact] "
      "[--report-skipped-pos] | "
      "-systemverilog/-sv [--sv_design1_flist <file>] [--sv_design1_top <name>] "
      "[--sv_design2_flist <file>] [--sv_design2_top <name>] [-v <lec|sec>] [-k <max-k>] [--sec-engine <legacy|k_induction|imc|pdr>] [--sec-encoding <binary|dual_rail_steady>] "
      "[--design1 <file...>] [--design2 <file...>] "
      "[--sec-internal-state-correspondence] [--no-sec-uncomputable-seq-boundary] [--compact] "
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
  if (token == "legacy") {
    // LCOV_EXCL_START
    engine = KEPLER_FORMAL::SEC::SecEngine::Legacy;
    return true;
    // LCOV_EXCL_STOP
  }
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
  error = "expected legacy, k_induction, imc, or pdr, got `" + token + "`";
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
    // LCOV_EXCL_START
    case KEPLER_FORMAL::SEC::SecEngine::Legacy:
    // LCOV_EXCL_STOP
    default:
      // LCOV_EXCL_START
      return "legacy";
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
  spdlog::set_level(level);
  return chosenLogFile;
}

// LCOV_EXCL_START
static bool parseMaxKToken(const std::string& token,
// LCOV_EXCL_STOP
                           size_t& maxK,
                           std::string& error) {
  // LCOV_EXCL_START
  if (token.empty()) {
    error = "max_k must not be empty";
    return false;
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  if (!std::all_of(token.begin(), token.end(), [](unsigned char ch) { return std::isdigit(ch); })) {
    error = "max_k must be a non-negative integer";
    return false;
    // LCOV_EXCL_STOP
  }
  try {
    // LCOV_EXCL_START
    maxK = static_cast<size_t>(std::stoull(token));
  } catch (const std::exception&) {
    error = "max_k is out of range";
    return false;
  }
  return true;
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
      "sec_uncomputable_seq_as_boundary",
      "sec_steady_frontier_guard",
      "sec_internal_state_correspondence",
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
      "sv_design1_top",
      "sv_design2_top",
      "design1",
      "design2",
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
  report << "# - top_input / top_output: original top-level interface terms.\n";
  report << "# - opaque_internal_input / opaque_internal_output: internal leaf cut points\n";
  report << "#   that SEC could not reconstruct combinationally and did not model as sequential.\n";
  report << "# - abstracted_sequential_state / abstracted_sequential_observed: interface terms\n";
  report << "#   exposed when an uncomputable sequential instance is abstracted as a SEC boundary.\n\n";
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

  if (out.design0.empty() || out.design1.empty()) {
    // LCOV_EXCL_START
    error = "input_paths must contain at least one file per design";
    return false;
    // LCOV_EXCL_STOP
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
    const SystemVerilogOptions& systemVerilogOptions) {
  if (normalizeInputListForComparison(designInputs.design0) !=
      normalizeInputListForComparison(designInputs.design1)) {
    return false;
  }
  // LCOV_EXCL_START
  if (!isSystemVerilog) {
    return true;
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return sameSystemVerilogDesignOptions(
      systemVerilogOptions.design0, systemVerilogOptions.design1);
      // LCOV_EXCL_STOP
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

// LCOV_EXCL_START
static std::vector<std::filesystem::path> buildSystemVerilogInputPaths(
// LCOV_EXCL_STOP
    const std::vector<std::string>& designInputs,
    const SystemVerilogDesignOptions& designOptions,
    std::vector<std::filesystem::path>& temporaryFiles) {
  // LCOV_EXCL_START
  std::vector<std::filesystem::path> svInputPaths = toPathVector(designInputs);
  // LCOV_EXCL_STOP

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
  if (designOptions.top) {
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
    topFlist << "--top " << *designOptions.top << "\n";
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
    const KEPLER_FORMAL::BuildPrimaryOutputClauses& builder) {
  // LCOV_EXCL_START
  KEPLER_FORMAL::MiterStrategy::CompactSnapshot snapshot;
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

int KeplerFormalMain(int argc, char** argv) {
  using namespace std::chrono;
  enum class FormatType { VERILOG, SYSTEMVERILOG, NAJA_IF, C };
  struct CDesignOptions {
    std::optional<std::string> top;
    std::optional<std::string> clock;
    std::optional<std::string> reset;
    std::vector<std::string> includePaths;
    std::optional<std::string> workDir;
    bool keepArtifacts = false;
  };
  struct DesignFormatOptions {
    FormatType format = FormatType::VERILOG;
    CDesignOptions c;
  };
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
  DesignFormatOptions designFormatOptions[2];
  bool usedDesignSections = false;
  std::vector<std::string> libertyFiles;
  std::vector<std::string> pythonFiles;
  std::string logLevel = "info";
  VerificationMode verificationMode = VerificationMode::LEC;
  KEPLER_FORMAL::SEC::SecEngine secEngine = KEPLER_FORMAL::SEC::SecEngine::Legacy;
  KEPLER_FORMAL::SEC::SecEncoding secEncoding =
      KEPLER_FORMAL::SEC::SecEncoding::DualRailSteady;
  bool secEngineExplicit = false;
  bool secEncodingExplicit = false;
  size_t secMaxK = kDefaultSecMaxK;
  bool secMaxKExplicit = false;
  bool secTreatUncomputableSeqAsBoundary = true;

  // Basic argument sanity
  if (argc < 2) {
    // LCOV_EXCL_START
    print_usage(argv[0]);  // LCOV_EXCL_LINE
    return EXIT_SUCCESS;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // Check for config mode (--config or -c). If present, YAML takes precedence.
  bool usedConfig = false;

  std::string logFileName;

  bool useScopes = false;
  bool cleanScopes = false;
  bool dumpCnf = false;
  bool dumpPoCnf = false;
  bool compactMode = false;
  bool reportSkippedPOs = false;
  bool verilogPreprocessing = false;
  bool secSteadyFrontierGuard = false;
  bool secInternalStateCorrespondence = false;
  std::string dumpCnfPath;
  std::string dumpPoCnfPath;

  KEPLER_FORMAL::Config::setReportSkippedPOs(false);
  KEPLER_FORMAL::Config::setSecTreatUncomputableSeqAsBoundary(true);
  KEPLER_FORMAL::Config::setSecSteadyFrontierGuard(false);
  KEPLER_FORMAL::Config::setSecInternalStateCorrespondence(false);

  const auto parseFormatToken = [](const std::string& fmt,
                                   FormatType& format) {
    if (fmt == "naja_if") {
      format = FormatType::NAJA_IF;
      return true;
    }
    if (fmt == "verilog" || fmt == "v") {
      format = FormatType::VERILOG;
      return true;
    }
    if (fmt == "systemverilog" || fmt == "sv") {
      format = FormatType::SYSTEMVERILOG;
      return true;
    }
    if (fmt == "c") {
      format = FormatType::C;
      return true;
    }
    return false;
  };

  const auto formatTypeName = [](FormatType format) {
    switch (format) {
      case FormatType::NAJA_IF:
        return "SNL";
      case FormatType::SYSTEMVERILOG:
        return "SYSTEMVERILOG";
      case FormatType::C:
        return "C";
      case FormatType::VERILOG:
      default:
        return "VERILOG";
    }
  };

  const auto yamlScalarOrSequenceToVector =
      [](const YAML::Node& node, std::vector<std::string>& out,
         std::string& error, const char* key) {
        out.clear();
        if (!node) {
          return true;
        }
        if (node.IsScalar()) {
          out.emplace_back(node.as<std::string>());
          return true;
        }
        if (!node.IsSequence()) {
          error = std::string(key) + " must be a scalar or sequence";
          return false;
        }
        for (const auto& entry : node) {
          if (!entry.IsScalar()) {
            error = std::string(key) + " entries must be scalar file paths";
            return false;
          }
          out.emplace_back(entry.as<std::string>());
        }
        return true;
      };

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--config" || a == "-c") {
      if (i + 1 >= argc) {
        // LCOV_EXCL_START
        SPDLOG_CRITICAL("Missing config file after {}", a);  // LCOV_EXCL_LINE
        return EXIT_FAILURE;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
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
          if (!parseFormatToken(fmt, inputFormatType)) {
            SPDLOG_CRITICAL("Unrecognized format in config: {}", fmt);
            return EXIT_FAILURE;
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

        if (cfg["sec_uncomputable_seq_as_boundary"]) {
          // LCOV_EXCL_START
          if (!cfg["sec_uncomputable_seq_as_boundary"].IsScalar()) {
            SPDLOG_CRITICAL(
            // LCOV_EXCL_STOP
                "sec_uncomputable_seq_as_boundary must be a scalar");
            // LCOV_EXCL_START
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          }
          // LCOV_EXCL_START
          secTreatUncomputableSeqAsBoundary =
              cfg["sec_uncomputable_seq_as_boundary"].as<bool>();
        }
        // LCOV_EXCL_STOP

        if (cfg["sec_steady_frontier_guard"]) {
          if (!cfg["sec_steady_frontier_guard"].IsScalar()) {
            SPDLOG_CRITICAL("sec_steady_frontier_guard must be a scalar");
            return EXIT_FAILURE;
          }
          secSteadyFrontierGuard = cfg["sec_steady_frontier_guard"].as<bool>();
        }

        if (cfg["sec_internal_state_correspondence"]) {
          if (!cfg["sec_internal_state_correspondence"].IsScalar()) {
            SPDLOG_CRITICAL("sec_internal_state_correspondence must be a scalar");
            return EXIT_FAILURE;
          }
          secInternalStateCorrespondence =
              cfg["sec_internal_state_correspondence"].as<bool>();
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

        if (cfg["design1"] || cfg["design2"]) {
          if (!cfg["design1"] || !cfg["design2"]) {
            SPDLOG_CRITICAL("design1 and design2 must both be present when using per-design config");
            return EXIT_FAILURE;
          }
          if (!cfg["design1"].IsMap() || !cfg["design2"].IsMap()) {
            SPDLOG_CRITICAL("design1 and design2 must be YAML maps");
            return EXIT_FAILURE;
          }
          usedDesignSections = true;

          const auto parseDesignSection =
              [&](const YAML::Node& designNode, int designIndex,
                  const char* designKey) {
                static const std::unordered_set<std::string> kAllowedDesignKeys = {
                    "format",
                    "input_paths",
                    "top",
                    "sv_flist",
                    "flist",
                    "clock",
                    "reset",
                    "include_paths",
                    "work_dir",
                    "keep_generated",
                    "keep_c2rtl_work",
                };

                for (auto it = designNode.begin(); it != designNode.end(); ++it) {
                  if (!it->first.IsScalar()) {
                    SPDLOG_CRITICAL("{} contains a non-scalar key", designKey);
                    return false;
                  }
                  const auto key = it->first.as<std::string>();
                  if (kAllowedDesignKeys.find(key) == kAllowedDesignKeys.end()) {
                    SPDLOG_CRITICAL("Unknown {} option: {}", designKey, key);
                    return false;
                  }
                }

                auto& designFormat = designFormatOptions[designIndex];
                auto& designInput =
                    designIndex == 0 ? designInputs.design0 : designInputs.design1;
                auto& svOptions = designIndex == 0
                                      ? systemVerilogOptions.design0
                                      : systemVerilogOptions.design1;

                if (!designNode["format"] || !designNode["format"].IsScalar()) {
                  SPDLOG_CRITICAL("{} requires scalar format", designKey);
                  return false;
                }
                const auto fmt = designNode["format"].as<std::string>();
                if (!parseFormatToken(fmt, designFormat.format)) {
                  SPDLOG_CRITICAL("Unrecognized {} format: {}", designKey, fmt);
                  return false;
                }
                if (designFormat.format == FormatType::NAJA_IF) {
                  SPDLOG_CRITICAL("{} format naja_if is not supported in per-design config yet", designKey);
                  return false;
                }

                std::string vectorError;
                if (!yamlScalarOrSequenceToVector(
                        designNode["input_paths"], designInput, vectorError, "input_paths")) {
                  SPDLOG_CRITICAL("Invalid {} input_paths: {}", designKey, vectorError);
                  return false;
                }

                const auto parseScalar = [&](const char* key,
                                             std::optional<std::string>& target) {
                  const auto node = designNode[key];
                  if (!node) {
                    return true;
                  }
                  if (!node.IsScalar()) {
                    SPDLOG_CRITICAL("{} {} must be a scalar", designKey, key);
                    return false;
                  }
                  const auto value = node.as<std::string>();
                  if (value.empty()) {
                    SPDLOG_CRITICAL("{} {} must not be empty", designKey, key);
                    return false;
                  }
                  target = value;
                  return true;
                };

                if (!parseScalar("top", svOptions.top) ||
                    !parseScalar("top", designFormat.c.top)) {
                  return false;
                }
                if (!parseScalar("sv_flist", svOptions.flist) ||
                    !parseScalar("flist", svOptions.flist) ||
                    !parseScalar("clock", designFormat.c.clock) ||
                    !parseScalar("reset", designFormat.c.reset) ||
                    !parseScalar("work_dir", designFormat.c.workDir)) {
                  return false;
                }

                if (designNode["include_paths"]) {
                  if (!yamlScalarOrSequenceToVector(
                          designNode["include_paths"],
                          designFormat.c.includePaths,
                          vectorError,
                          "include_paths")) {
                    SPDLOG_CRITICAL("Invalid {} include_paths: {}", designKey, vectorError);
                    return false;
                  }
                }
                if (designNode["keep_generated"]) {
                  if (!designNode["keep_generated"].IsScalar()) {
                    SPDLOG_CRITICAL("{} keep_generated must be a scalar", designKey);
                    return false;
                  }
                  designFormat.c.keepArtifacts =
                      designNode["keep_generated"].as<bool>();
                }
                if (designNode["keep_c2rtl_work"]) {
                  if (!designNode["keep_c2rtl_work"].IsScalar()) {
                    SPDLOG_CRITICAL("{} keep_c2rtl_work must be a scalar", designKey);
                    return false;
                  }
                  designFormat.c.keepArtifacts =
                      designNode["keep_c2rtl_work"].as<bool>();
                }

                if (designFormat.format == FormatType::C && designInput.empty()) {
                  SPDLOG_CRITICAL("{} format c requires input_paths", designKey);
                  return false;
                }
                if ((designFormat.format == FormatType::VERILOG ||
                     designFormat.format == FormatType::SYSTEMVERILOG) &&
                    designInput.empty() && !svOptions.flist) {
                  SPDLOG_CRITICAL(
                      "{} requires input_paths or flist/sv_flist", designKey);
                  return false;
                }
                return true;
              };

          if (!parseDesignSection(cfg["design1"], 0, "design1") ||
              !parseDesignSection(cfg["design2"], 1, "design2")) {
            return EXIT_FAILURE;
          }

          const bool anyC =
              designFormatOptions[0].format == FormatType::C ||
              designFormatOptions[1].format == FormatType::C;
          const bool anySystemVerilog =
              designFormatOptions[0].format == FormatType::SYSTEMVERILOG ||
              designFormatOptions[1].format == FormatType::SYSTEMVERILOG;
          if (anyC || anySystemVerilog) {
            inputFormatType = FormatType::SYSTEMVERILOG;
          } else {
            inputFormatType = FormatType::VERILOG;
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
                cfg, "sv_design2_top", systemVerilogOptions.design1.top, svConfigError)) {
          // LCOV_EXCL_START
          SPDLOG_CRITICAL("Invalid SystemVerilog config option: {}", svConfigError);
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
      if (arg == "--sec-uncomputable-seq-boundary") {
        secTreatUncomputableSeqAsBoundary = true;
        ++parseStart;
        continue;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "--no-sec-uncomputable-seq-boundary") {
        secTreatUncomputableSeqAsBoundary = false;
        ++parseStart;
        continue;
        // LCOV_EXCL_STOP
      }
      if (arg == "--sec-internal-state-correspondence") {
        secInternalStateCorrespondence = true;
        ++parseStart;
        continue;
      }
      if (arg == "--no-sec-internal-state-correspondence") {
        secInternalStateCorrespondence = false;
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
      if (arg == "--sec-uncomputable-seq-boundary") {
        secTreatUncomputableSeqAsBoundary = true;
        continue;
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (arg == "--no-sec-uncomputable-seq-boundary") {
        secTreatUncomputableSeqAsBoundary = false;
        continue;
        // LCOV_EXCL_STOP
      }
      if (arg == "--sec-internal-state-correspondence") {
        secInternalStateCorrespondence = true;
        continue;
      }
      if (arg == "--no-sec-internal-state-correspondence") {
        secInternalStateCorrespondence = false;
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
        } else {
          systemVerilogOptions.design1.top = value;
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
  std::string inputFormatName = formatTypeName(inputFormatType);
  if (usedDesignSections) {
    inputFormatName = std::string(formatTypeName(designFormatOptions[0].format)) +
                      "/" + formatTypeName(designFormatOptions[1].format);
  }
  SPDLOG_INFO("Input format: {}", inputFormatName);
  if (!runLogFilePath.empty()) {
    SPDLOG_INFO("Run log: {}", runLogFilePath);
  }
  logDesignPaths("Netlist 1", designInputs.design0);
  logDesignPaths("Netlist 2", designInputs.design1);

  std::vector<KEPLER_FORMAL::C2RTL::CFrontendResult> cFrontendResults;
  struct CFrontendCleanupGuard {
    std::vector<KEPLER_FORMAL::C2RTL::CFrontendResult>* results;
    ~CFrontendCleanupGuard() {
      if (results == nullptr) {
        return;
      }
      for (const auto& result : *results) {
        if (result.keepArtifacts) {
          continue;
        }
        std::error_code ec;
        std::filesystem::remove_all(result.workDir, ec);
      }
    }
  } cFrontendCleanupGuard{&cFrontendResults};
  const auto materializeCDesign =
      [&](int designIndex,
          std::vector<std::string>& designPaths,
          SystemVerilogDesignOptions& svOptions) {
        const auto& cOptions = designFormatOptions[designIndex].c;
        KEPLER_FORMAL::C2RTL::CFrontendOptions options;
        options.designIndex = designIndex;
        options.top = cOptions.top;
        options.clock = cOptions.clock;
        options.reset = cOptions.reset;
        options.keepArtifacts = cOptions.keepArtifacts;
        if (cOptions.workDir) {
          options.workDir = std::filesystem::path(*cOptions.workDir);
        }
        for (const auto& inputPath : designPaths) {
          options.inputPaths.emplace_back(inputPath);
        }
        for (const auto& includePath : cOptions.includePaths) {
          options.includePaths.emplace_back(includePath);
        }

        KEPLER_FORMAL::C2RTL::CFrontend frontend;
        auto result = frontend.translateToSystemVerilog(options);
        SPDLOG_INFO(
            "C frontend generated SystemVerilog for design {}: {}",
            designIndex + 1,
            result.generatedSystemVerilog.string());
        designPaths.clear();
        designPaths.emplace_back(result.generatedSystemVerilog.string());
        svOptions.top = cOptions.top;
        cFrontendResults.push_back(std::move(result));
      };

  if (inputFormatType == FormatType::C && !usedDesignSections) {
    SPDLOG_CRITICAL("format: c requires per-design YAML sections design1 and design2");
    return EXIT_FAILURE;
  }

  if (usedDesignSections) {
    try {
      if (designFormatOptions[0].format == FormatType::C) {
        materializeCDesign(
            0, designInputs.design0, systemVerilogOptions.design0);
        designFormatOptions[0].format = FormatType::SYSTEMVERILOG;
      }
      if (designFormatOptions[1].format == FormatType::C) {
        materializeCDesign(
            1, designInputs.design1, systemVerilogOptions.design1);
        designFormatOptions[1].format = FormatType::SYSTEMVERILOG;
      }
    } catch (const std::exception& e) {
      SPDLOG_CRITICAL("C frontend failed: {}", e.what());
      return EXIT_FAILURE;
    }
    inputFormatType = FormatType::SYSTEMVERILOG;
  }

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
      (systemVerilogOptions.design0.flist || systemVerilogOptions.design0.top ||
       systemVerilogOptions.design1.flist || systemVerilogOptions.design1.top)) {
    // LCOV_EXCL_START
    SPDLOG_CRITICAL("SystemVerilog design options are only valid with -systemverilog/-sv input");
    return EXIT_FAILURE;
    // LCOV_EXCL_STOP
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
  KEPLER_FORMAL::Config::setSecTreatUncomputableSeqAsBoundary(
      secTreatUncomputableSeqAsBoundary);
  KEPLER_FORMAL::Config::setSecSteadyFrontierGuard(secSteadyFrontierGuard);
  KEPLER_FORMAL::Config::setSecInternalStateCorrespondence(
      secInternalStateCorrespondence);
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
    SPDLOG_INFO(
        "SEC uncomputable sequentials: {}",
        secTreatUncomputableSeqAsBoundary ? "boundary abstraction"
                                          : "strict failure");
    SPDLOG_INFO(
        "SEC internal state correspondence: {}",
        secInternalStateCorrespondence ? "enabled" : "disabled");
  }
  SPDLOG_INFO("Compact mode: {}", compactMode ? "enabled" : "disabled");
  SPDLOG_INFO("Skipped PO reports: {}", reportSkippedPOs ? "enabled" : "disabled");
  if (!libertyFiles.empty()) {
    for (const auto& lf : libertyFiles) SPDLOG_INFO("Library: {}", lf);
  }
  if (!pythonFiles.empty()) {
    // LCOV_EXCL_START
    for (const auto& pf : pythonFiles) SPDLOG_INFO("Python library: {}", pf);
  }
  // LCOV_EXCL_STOP

  auto emitSecResult =
      [&](const KEPLER_FORMAL::SEC::SequentialEquivalenceResult& result) {
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
            SPDLOG_INFO(
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
        // LCOV_EXCL_START
        if (!result.abstractedSequentialBoundaries.empty()) {
          // LCOV_DISABLED_START
          std::ostringstream abstractedBoundaries;
          for (const auto& abstractedBoundary :
               result.abstractedSequentialBoundaries) {
            abstractedBoundaries << "  - " << abstractedBoundary << "\n";
            // LCOV_DISABLED_STOP
          }
          // LCOV_DISABLED_START
          SPDLOG_INFO(
          // LCOV_DISABLED_STOP
              "SEC abstracted uncomputable sequential interfaces as "
              "boundaries:\n{}",
              abstractedBoundaries.str());
        // LCOV_DISABLED_START
        }
        // LCOV_DISABLED_STOP
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
        }
        // LCOV_EXCL_STOP
        switch (result.status) {
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Equivalent:
            SPDLOG_INFO(
                "No difference was found. SEC proved equivalence at k = {}.",
                result.bound);
            return EXIT_SUCCESS;
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Different:
            // LCOV_EXCL_START
            SPDLOG_INFO(
            // LCOV_EXCL_STOP
                "Difference was found. SEC found a counterexample at k = {}.",
                result.bound);
            // LCOV_EXCL_START
            if (!result.reason.empty()) {
              SPDLOG_INFO("SEC counterexample details:\n{}", result.reason);
            }
            return EXIT_SUCCESS;
            // LCOV_EXCL_STOP
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Inconclusive:
            // LCOV_EXCL_START
            if (secInconclusiveStoppedBeforeMaxK(result.reason)) {
              SPDLOG_CRITICAL(
                  "SEC was inconclusive before completing max_k = {}: {}",
                  secMaxK,
                  result.reason);
            } else {
              SPDLOG_CRITICAL(
                  "SEC was inconclusive up to max_k = {}: {}",
                  secMaxK,
                  result.reason);
            }
            return EXIT_FAILURE;
            // LCOV_EXCL_STOP
          case KEPLER_FORMAL::SEC::SequentialEquivalenceStatus::Unsupported:
          // LCOV_EXCL_STOP
          default:
            // LCOV_EXCL_START
            SPDLOG_CRITICAL(
            // LCOV_EXCL_STOP
                "SEC cannot run on this design pair: {}", result.reason);
            // LCOV_EXCL_START
            return EXIT_FAILURE;
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

    auto loadLibraries = [&](NLDB* db) -> bool {
      if (libertyFiles.empty() && pythonFiles.empty()) {
        // LCOV_EXCL_START
        // LCOV_DISABLED_START
        return false;
        // LCOV_DISABLED_STOP
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
        // LCOV_EXCL_START
        std::filesystem::path pythonPath(pythonFile);
        SPDLOG_INFO("Loading python primitive file: {}", pythonFile);
        SNLPyLoader::loadPrimitives(primitivesLibrary, pythonPath);
      }
      // LCOV_EXCL_STOP
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
          // LCOV_EXCL_START
          throw std::runtime_error("Failed to load library files");
          // LCOV_EXCL_STOP
        }
      }

      if (inputFormatType == FormatType::VERILOG ||
          // LCOV_EXCL_START
          inputFormatType == FormatType::SYSTEMVERILOG) {
          // LCOV_EXCL_STOP
        if (!db) {
          // LCOV_EXCL_START
          db = NLDB::create(NLUniverse::get());
        }
        // LCOV_EXCL_STOP
        db->setID(dbID);
        SPDLOG_INFO("Parsing {} file(s) for design {}",
                    inputFormatType == FormatType::SYSTEMVERILOG ? "systemverilog" : "verilog",
                    designIndex + 1);
        auto designLibrary = NLLibrary::create(db, NLName("DESIGN"));
        if (inputFormatType == FormatType::SYSTEMVERILOG) {
          // LCOV_EXCL_START
          SNLSVConstructor constructor(designLibrary);
          std::vector<std::filesystem::path> temporaryFiles;
          // LCOV_EXCL_STOP
	          const auto svInputPaths =
	              // LCOV_EXCL_START
	              buildSystemVerilogInputPaths(designPaths, designOptions, temporaryFiles);
	              // LCOV_EXCL_STOP
	          try {
	            // LCOV_EXCL_START
	            constructor.construct(svInputPaths);
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
	        auto top = SNLUtils::findTop(designLibrary);
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
          auto inputs = builder.getInputs();
          auto outputs = builder.getOutputs();
          builder.setInputs(inputs);
          builder.setOutputs(outputs);
          builder.build();
          return captureCompactSnapshot(builder);
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
        if (dumpCnf) {
          const std::string outPath = dumpCnfPath.empty() ? "miter.cnf" : dumpCnfPath;
          MiterS.setCnfDump(true, outPath);
        }
        if (dumpPoCnf) {
          const std::string outPath = dumpPoCnfPath.empty() ? "po_cnfs" : dumpPoCnfPath;
          MiterS.setPoCnfDump(true, outPath);
        }
        if (MiterS.runCompactSnapshots(snapshot0, snapshot1)) {
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
        if (sameCompactSecDesignSpec(
                inputFormatType == FormatType::SYSTEMVERILOG,
                designInputs,
                systemVerilogOptions)) {
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
              nullptr, nullptr, solverType, secEngine, secEncoding);
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
            nullptr, nullptr, solverType, secEngine, secEncoding);
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
    if (!libertyFiles.empty() || !pythonFiles.empty()) {
      db0 = NLDB::create(NLUniverse::get());
      primitivesAreLoaded = loadLibraries(db0);
      if (!primitivesAreLoaded) {
        return EXIT_FAILURE;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
    if (inputFormatType == FormatType::VERILOG ||
        inputFormatType == FormatType::SYSTEMVERILOG) {
      if (!db0) {
        db0 = NLDB::create(NLUniverse::get());
      }
      const auto design0Paths = toPathVector(designInputs.design0);
      SPDLOG_INFO("Parsing {} file(s) for design 1",
      // LCOV_EXCL_STOP
                  inputFormatType == FormatType::SYSTEMVERILOG ? "systemverilog" : "verilog");
      // LCOV_EXCL_START
      auto designLibrary = NLLibrary::create(db0, NLName("DESIGN"));
      if (inputFormatType == FormatType::SYSTEMVERILOG) {
        SNLSVConstructor constructor(designLibrary);
        std::vector<std::filesystem::path> temporaryFiles;
        const auto svInputPaths = buildSystemVerilogInputPaths(
            designInputs.design0, systemVerilogOptions.design0, temporaryFiles);
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
        constructor.construct(design0Paths);
      }
      auto top = SNLUtils::findTop(designLibrary);
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
    if (inputFormatType == FormatType::VERILOG ||
        inputFormatType == FormatType::SYSTEMVERILOG) {
      if (!db1) {
        db1 = NLDB::create(NLUniverse::get());
      }
      const auto design1Paths = toPathVector(designInputs.design1);
      SPDLOG_INFO("Parsing {} file(s) for design 2",
      // LCOV_EXCL_STOP
                  inputFormatType == FormatType::SYSTEMVERILOG ? "systemverilog" : "verilog");
      // LCOV_EXCL_START
      auto designLibrary = NLLibrary::create(db1, NLName("DESIGN"));
      if (inputFormatType == FormatType::SYSTEMVERILOG) {
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
      auto top = SNLUtils::findTop(designLibrary);
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
          top0, top1, solverType, secEngine, secEncoding);
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
    MiterS.init(false);
    ScopeExtraction extractor(top0, top1);
    extractor.collectVerificationScopes();
    if (cleanScopes) {
      extractor.cleanVerificationScopes(MiterS.getPIs0(), MiterS.getPIs1());
    }
    for (auto scopes : extractor.getScopesToVerify()) {
      SPDLOG_INFO("Looking at scope: {} ",
      // LCOV_EXCL_STOP
                  scopes.first->getName().getString());
      // std::string scopeLogFile = (logFileName.empty() ? "kf_" : logFileName) + "_" +
      //                           scopes.first->getName().getString() + ".txt";
      try {
        // LCOV_EXCL_START
        KEPLER_FORMAL::MiterStrategy MiterScope(scopes.first, scopes.second, logFileName);
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
  // LCOV_EXCL_START
  } else {
  // LCOV_EXCL_STOP
    try {
      // LCOV_EXCL_START
      KEPLER_FORMAL::MiterStrategy MiterS(top0, top1, logFileName);
      if (dumpCnf) {
        const std::string outPath = dumpCnfPath.empty() ? "miter.cnf" : dumpCnfPath;
        MiterS.setCnfDump(true, outPath);
      }
      if (dumpPoCnf) {
        const std::string outPath = dumpPoCnfPath.empty() ? "po_cnfs" : dumpPoCnfPath;
        MiterS.setPoCnfDump(true, outPath);
      }
      MiterS.init();
      if (MiterS.run(compactMode)) {
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

#ifndef KEPLER_FORMAL_NO_MAIN
int main(int argc, char** argv) {
  const int rc = KeplerFormalMain(argc, argv);
#if defined(KEPLER_FORMAL_ASAN_BUILD)
  // Subprocess sanitizer tests check leaks before process teardown.
  KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
  KEPLER_FORMAL::BoolExprCache::destroy();
#endif
  return rc;
}
#endif
