// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "KeplerFormalDriver.h"
#include "SNLPyLoader.h"

#include <cstdlib>
#include <memory>
#include <sstream>
#include <stdexcept>

#ifdef KEPLER_PYTHON3_ROOTPATH
#include "rules_cc/cc/runfiles/runfiles.h"
#endif

namespace {

#ifdef KEPLER_PYTHON3_ROOTPATH
// Bazel builds embed the toolchain's hermetic Python, whose compiled-in
// prefix does not exist at run time. Point it at the runtime that ships in
// this binary's runfiles, unless the caller chose a PYTHONHOME (the release
// tarball's wrapper script does).
static void setBundledPythonHome(const char* argv0) {
  if (std::getenv("PYTHONHOME")) {
    return;
  }
  std::string error;
  std::unique_ptr<rules_cc::cc::runfiles::Runfiles> runfiles(
      rules_cc::cc::runfiles::Runfiles::Create(argv0, BAZEL_CURRENT_REPOSITORY, &error));
  if (!runfiles) {
    return;
  }
  // $(PYTHON3_ROOTPATH) of an external repository is "../<repo>/bin/python3";
  // its runfiles location drops the "../".
  std::string rootpath(KEPLER_PYTHON3_ROOTPATH);
  const std::string rlocation =
      rootpath.rfind("../", 0) == 0 ? rootpath.substr(3) : "_main/" + rootpath;
  const std::filesystem::path interpreter(runfiles->Rlocation(rlocation));
  std::error_code ec;
  if (interpreter.empty() || !std::filesystem::exists(interpreter, ec)) {
    return;
  }
  const auto home = interpreter.parent_path().parent_path();
  if (setenv("PYTHONHOME", home.string().c_str(), 1) != 0) {
    throw std::runtime_error("Cannot configure PYTHONHOME for Naja primitives");  // LCOV_EXCL_LINE
  }
}
#endif

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
#ifdef KEPLER_PYTHON3_ROOTPATH
    setBundledPythonHome(executable);
#endif
    addNajaPythonPath(executable);
  }
  void load(naja::NL::NLLibrary* library,
            const std::filesystem::path& path) const override {
    naja::NL::SNLPyLoader::loadPrimitives(library, path);
  }
};

class InProcessPrimitiveLoader final : public KEPLER_FORMAL::PrimitiveLibraryLoader {
 public:
  void prepare(const char*) const override {
    throw std::runtime_error(
        "py_tech_files are not supported by the in-process file API");
  }
  void load(naja::NL::NLLibrary*, const std::filesystem::path&) const override {
    // Shared validation calls prepare before any loading begins.
    throw std::logic_error("In-process Python primitive loading was not validated");
  }
};

}  // namespace

namespace KEPLER_FORMAL {

int runKeplerFormal(int argc, char** argv, RunResult& result) {
  InProcessPrimitiveLoader primitives;
  return runKeplerFormal(argc, argv, result, primitives);
}

}  // namespace KEPLER_FORMAL

int KeplerFormalMain(int argc, char** argv) {
  StandalonePrimitiveLoader primitives;
  KEPLER_FORMAL::RunResult result;
  return KEPLER_FORMAL::runKeplerFormalWorkflow(argc, argv, result, primitives);
}
