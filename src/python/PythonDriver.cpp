// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "PythonDriver.h"

#include <stdexcept>

namespace KEPLER_FORMAL {
namespace {

class InProcessPrimitiveLoader final : public PrimitiveLibraryLoader {
 public:
  void prepare(const char*) const override {
    throw std::runtime_error(
        "py_tech_files are not supported by the in-process Python API");
  }
  void load(naja::NL::NLLibrary*, const std::filesystem::path&) const override {
    // Shared validation calls prepare before any loading begins.
    throw std::logic_error("In-process Python primitive loading was not validated");
  }
};

}  // namespace

int runPythonVerification(int argc, char** argv, RunResult& result) {
  InProcessPrimitiveLoader primitives;
  return runKeplerFormal(argc, argv, result, primitives);
}

int runKeplerFormal(int argc, char** argv, RunResult& result) {
  return runPythonVerification(argc, argv, result);
}

}  // namespace KEPLER_FORMAL
