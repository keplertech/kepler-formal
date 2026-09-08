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

// Runs one complete verification in the current process. Calls are not
// reentrant: the implementation owns Naja's process-global universe while it
// runs. The Python binding serializes access to this entry point.
int runKeplerFormal(int argc, char **argv, RunResult &result);

// Releases process-global expression caches retained after a run. The native
// Python binding calls this after every invocation because, unlike the command
// line program, its process remains alive.
void cleanupKeplerFormalState();

} // namespace KEPLER_FORMAL

// Compatibility entry point used by the command-line executable and existing
// in-process CLI tests.
int KeplerFormalMain(int argc, char **argv);
