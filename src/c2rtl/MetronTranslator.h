// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <filesystem>
#include <vector>

namespace KEPLER_FORMAL::C2RTL {

struct MetronTranslationOptions {
  std::filesystem::path inputPath;
  std::filesystem::path outputPath;
  std::filesystem::path stdoutLog;
  std::filesystem::path stderrLog;
  std::vector<std::filesystem::path> includePaths;
};

void translateWithMetron(const MetronTranslationOptions& options);

}  // namespace KEPLER_FORMAL::C2RTL
