// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "CFrontend.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace KEPLER_FORMAL::C2RTL {
namespace {

std::string shellQuote(const std::string& text) {
  std::string quoted;
  quoted.reserve(text.size() + 2);
  quoted.push_back('\'');
  for (char c : text) {
    if (c == '\'') {
      quoted += "'\\''";
    } else {
      quoted.push_back(c);
    }
  }
  quoted.push_back('\'');
  return quoted;
}

std::string jsonQuote(const std::string& text) {
  std::string quoted;
  quoted.reserve(text.size() + 2);
  quoted.push_back('"');
  for (char c : text) {
    switch (c) {
      case '\\':
        quoted += "\\\\";
        break;
      case '"':
        quoted += "\\\"";
        break;
      case '\n':
        quoted += "\\n";
        break;
      case '\r':
        quoted += "\\r";
        break;
      case '\t':
        quoted += "\\t";
        break;
      default:
        quoted.push_back(c);
        break;
    }
  }
  quoted.push_back('"');
  return quoted;
}

std::filesystem::path defaultWorkDir(int designIndex) {
  const auto stamp =
      std::chrono::steady_clock::now().time_since_epoch().count();
  return std::filesystem::temp_directory_path() /
         ("kepler_c2rtl_design" + std::to_string(designIndex + 1) + "_" +
          std::to_string(stamp));
}

std::string translatorExecutable() {
  if (const char* env = std::getenv("KEPLER_C_FRONTEND_TRANSLATOR")) {
    if (*env != '\0') {
      return env;
    }
  }
  if (const char* env = std::getenv("KEPLER_C2RTL_TRANSLATOR")) {
    if (*env != '\0') {
      return env;
    }
  }
  return "metron";
}

std::set<std::filesystem::path> collectIncludePaths(
    const std::vector<std::filesystem::path>& inputPaths,
    const std::vector<std::filesystem::path>& explicitIncludePaths) {
  std::set<std::filesystem::path> includes;
  for (const auto& inputPath : inputPaths) {
    const auto parent = inputPath.parent_path();
    if (!parent.empty()) {
      includes.insert(parent);
    }
  }
  for (const auto& includePath : explicitIncludePaths) {
    includes.insert(includePath);
  }
  return includes;
}

void writeManifest(const CFrontendOptions& options,
                   const CFrontendResult& result,
                   const std::string& translator,
                   const std::string& command) {
  std::ofstream manifest(result.manifest, std::ios::out | std::ios::trunc);
  if (!manifest) {
    throw std::runtime_error(
        "C frontend failed to create manifest `" + result.manifest.string() + "`");
  }

  manifest << "{\n";
  manifest << "  \"frontend\": \"kepler-c\",\n";
  manifest << "  \"translator\": " << jsonQuote(translator) << ",\n";
  manifest << "  \"command\": " << jsonQuote(command) << ",\n";
  manifest << "  \"design_index\": " << options.designIndex << ",\n";
  manifest << "  \"top\": "
           << (options.top ? jsonQuote(*options.top) : std::string("null"))
           << ",\n";
  manifest << "  \"clock\": "
           << (options.clock ? jsonQuote(*options.clock) : std::string("null"))
           << ",\n";
  manifest << "  \"reset\": "
           << (options.reset ? jsonQuote(*options.reset) : std::string("null"))
           << ",\n";
  manifest << "  \"inputs\": [";
  for (size_t i = 0; i < options.inputPaths.size(); ++i) {
    if (i != 0) {
      manifest << ", ";
    }
    manifest << jsonQuote(options.inputPaths[i].string());
  }
  manifest << "],\n";
  manifest << "  \"generated_systemverilog\": "
           << jsonQuote(result.generatedSystemVerilog.string()) << ",\n";
  manifest << "  \"stdout_log\": " << jsonQuote(result.stdoutLog.string()) << ",\n";
  manifest << "  \"stderr_log\": " << jsonQuote(result.stderrLog.string()) << "\n";
  manifest << "}\n";
}

}  // namespace

CFrontendResult CFrontend::translateToSystemVerilog(
    const CFrontendOptions& options) const {
  if (options.inputPaths.empty()) {
    throw std::runtime_error("C frontend requires at least one C input path");
  }
  if (!options.top || options.top->empty()) {
    throw std::runtime_error("C frontend requires a top name");
  }

  std::error_code ec;
  for (const auto& inputPath : options.inputPaths) {
    if (!std::filesystem::exists(inputPath, ec)) {
      throw std::runtime_error(
          "C frontend input does not exist: `" + inputPath.string() + "`");
    }
  }

  CFrontendResult result;
  result.workDir = options.workDir.value_or(defaultWorkDir(options.designIndex));
  result.keepArtifacts = options.keepArtifacts;
  std::filesystem::create_directories(result.workDir);

  const auto outputBase =
      "design" + std::to_string(options.designIndex + 1) + "_" + *options.top;
  result.generatedSystemVerilog = result.workDir / (outputBase + "_from_c.sv");
  result.stdoutLog = result.workDir / "c_frontend_stdout.log";
  result.stderrLog = result.workDir / "c_frontend_stderr.log";
  result.manifest = result.workDir / "input_manifest.json";

  const auto translator = translatorExecutable();
  const auto includes = collectIncludePaths(options.inputPaths, options.includePaths);

  std::ostringstream command;
  command << shellQuote(translator)
          << " --quiet"
          << " --convert " << shellQuote(options.inputPaths.front().string())
          << " --output " << shellQuote(result.generatedSystemVerilog.string());
  for (const auto& includePath : includes) {
    command << " --include " << shellQuote(includePath.string());
  }
  command << " > " << shellQuote(result.stdoutLog.string())
          << " 2> " << shellQuote(result.stderrLog.string());

  const auto commandString = command.str();
  writeManifest(options, result, translator, commandString);

  const int rc = std::system(commandString.c_str());
  if (rc != 0) {
    throw std::runtime_error(
        "C frontend failed to translate `" + options.inputPaths.front().string() +
        "` to SystemVerilog; see `" + result.stderrLog.string() + "`");
  }

  ec.clear();
  if (!std::filesystem::exists(result.generatedSystemVerilog, ec) ||
      std::filesystem::file_size(result.generatedSystemVerilog, ec) == 0) {
    throw std::runtime_error(
        "C frontend translator did not produce SystemVerilog output `" +
        result.generatedSystemVerilog.string() + "`");
  }

  return result;
}

}  // namespace KEPLER_FORMAL::C2RTL
