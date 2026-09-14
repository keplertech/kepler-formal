// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <filesystem>

#include <string>
#include <vector>

namespace naja::NL {
class NLLibrary;
}

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

// Frontend-specific technology loading. The shared run API invokes prepare
// once, after configuration validation, before creating any design library.
class PrimitiveLibraryLoader {
 public:
  virtual ~PrimitiveLibraryLoader() = default;
  virtual void prepare(const char* executable) const = 0;
  virtual void load(naja::NL::NLLibrary* library,
                    const std::filesystem::path& path) const = 0;
};

// Shared workflow used by both adapters: argument/YAML parsing, file loading,
// LEC/SEC, exports and results. It owns and releases the run's Naja designs.
// The standalone adapter uses this directly; embedding callers use the guarded
// API below to preserve their process state and reject a foreign universe.
int runKeplerFormalWorkflow(int argc, char** argv, RunResult& result,
                            const PrimitiveLibraryLoader& primitiveLoader);

// Serialize an in-process run and restore its logger/solver/cache state.
// A live external Naja universe is rejected without being destroyed.
int runKeplerFormal(int argc, char** argv, RunResult& result,
                    const PrimitiveLibraryLoader& primitiveLoader);

// Compatibility entry point for existing in-process callers. It uses the
// Python driver's technology policy and is supplied by that driver library.
int runKeplerFormal(int argc, char** argv, RunResult& result);

// Releases process-global expression caches retained after a run. The native
// Python binding calls this after every invocation because, unlike the command
// line program, its process remains alive.
void cleanupKeplerFormalState();

} // namespace KEPLER_FORMAL

// Compatibility entry point used by the command-line executable and existing
// in-process CLI tests.
int KeplerFormalMain(int argc, char **argv);
