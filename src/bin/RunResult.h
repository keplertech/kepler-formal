// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace KEPLER_FORMAL {

enum class RunStatus {
  NoResult,
  Error,
  Equivalent,
  Different,
  PartiallyProved,
  Inconclusive,
  Unsupported,
  Exported,
};

struct RunResult {
  RunStatus status = RunStatus::Error;
  int exitCode = 1;
  std::string inputFormat;
  std::string verification;
  std::string logFile;
  size_t bound = 0;
  std::string reason;
  size_t coveredOutputs = 0;
  size_t totalOutputs = 0;
  size_t provenOutputs = 0;
  std::vector<std::string> unprovenOutputs;
  std::vector<std::string> skippedObservedOutputs;
};

const char *runStatusName(RunStatus status);

}  // namespace KEPLER_FORMAL
