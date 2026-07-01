// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace KEPLER_FORMAL::C2RTL {

struct XlsTranslationOptions {
  std::filesystem::path inputPath;
  std::filesystem::path outputPath;
  std::filesystem::path stdoutLog;
  std::filesystem::path stderrLog;
  std::vector<std::filesystem::path> includePaths;
  std::optional<std::filesystem::path> blockProto;
  std::optional<std::string> moduleName;
  std::string top;
};

void translateWithXls(const XlsTranslationOptions& options);

}  // namespace KEPLER_FORMAL::C2RTL
