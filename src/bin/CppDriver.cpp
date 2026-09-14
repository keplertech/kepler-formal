// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "KeplerFormalDriver.h"
#include "SNLPyLoader.h"

#include <cstdlib>
#include <sstream>
#include <stdexcept>

namespace {

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

class StandalonePrimitiveLoader final : public KEPLER_FORMAL::PrimitiveLibraryLoader {
 public:
  void prepare(const char* executable) const override {
    addNajaPythonPath(executable);
  }
  void load(naja::NL::NLLibrary* library,
            const std::filesystem::path& path) const override {
    naja::NL::SNLPyLoader::loadPrimitives(library, path);
  }
};

}  // namespace

int KeplerFormalMain(int argc, char** argv) {
  StandalonePrimitiveLoader primitives;
  KEPLER_FORMAL::RunResult result;
  return KEPLER_FORMAL::runKeplerFormalWorkflow(argc, argv, result, primitives);
}
