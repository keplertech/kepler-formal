// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <filesystem>
#include "RunResult.h"

namespace naja::NL {
class NLLibrary;
}

namespace KEPLER_FORMAL {

// Frontend-specific technology loading. The shared run API invokes prepare
// once, after configuration validation, before creating any design library.
class PrimitiveLibraryLoader {
 public:
  virtual ~PrimitiveLibraryLoader() = default;
  virtual void prepare(const char* executable) const = 0;
  virtual void load(naja::NL::NLLibrary* library,
                    const std::filesystem::path& path) const = 0;
};

// File workflow: argument/YAML parsing, file loading,
// LEC/SEC, exports and results. It owns and releases the run's Naja designs.
// The standalone adapter uses this directly; embedding callers use the guarded
// API below to preserve their process state and reject a foreign universe.
int runKeplerFormalWorkflow(int argc, char** argv, RunResult& result,
                            const PrimitiveLibraryLoader& primitiveLoader);

// Serialize an in-process run and restore its logger/solver/cache state.
// A live external Naja universe is rejected without being destroyed.
int runKeplerFormal(int argc, char** argv, RunResult& result,
                    const PrimitiveLibraryLoader& primitiveLoader);

// Compatibility entry point for C++ hosts that do not embed Python primitive
// loading. Python bindings use verifyBorrowedDesigns instead of this file API.
int runKeplerFormal(int argc, char** argv, RunResult& result);

// Releases process-global expression caches retained after a file-based run.
void cleanupKeplerFormalState();

} // namespace KEPLER_FORMAL

// Compatibility entry point used by the command-line executable and existing
// in-process CLI tests.
int KeplerFormalMain(int argc, char **argv);
