// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace KEPLER_FORMAL::C2RTL {

struct CFrontendOptions {
  std::string frontend = "auto";
  std::vector<std::filesystem::path> inputPaths;
  std::vector<std::filesystem::path> includePaths;
  std::optional<std::filesystem::path> blockProto;
  std::optional<std::string> moduleName;
  std::optional<std::string> top;
  std::optional<std::string> clock;
  std::optional<std::string> reset;
  std::optional<std::filesystem::path> workDir;
  bool keepArtifacts = false;
  int designIndex = 0;
};

struct CFrontendResult {
  std::filesystem::path workDir;
  std::filesystem::path generatedSystemVerilog;
  std::filesystem::path stdoutLog;
  std::filesystem::path stderrLog;
  std::filesystem::path manifest;
  bool keepArtifacts = false;
};

class CFrontend {
 public:
  CFrontendResult translateToSystemVerilog(const CFrontendOptions& options) const;
};

}  // namespace KEPLER_FORMAL::C2RTL
