// Copyright 2024-2025 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <chrono>
#include <cstdlib>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>
#include "NajaPerf.h"

// Naja interfaces
#include "DNL.h"
#include "MiterStrategy.h"
#include "SNLCapnP.h"
#include "SNLLibertyConstructor.h"
#include "SNLVRLConstructor.h"
#include "SNLVRLDumper.h"
#include "SNLUtils.h"

int main(int argc, char** argv) {
  using namespace std::chrono;
  enum class FormatType { VERILOG, SNL };

  FormatType inputFormatType = FormatType::VERILOG;

  // Help print when --help or -h is provided
  if (argc < 3 || (std::string(argv[1]) == "--help") ||
      (std::string(argv[1]) == "-h")) {
    printf(
        "Usage: kepler_formal <naja-if-dir-1> <naja-if-dir-2> "
        "[<liberty-file>...]\n");
    return EXIT_SUCCESS;
  }

  // --------------------------------------------------------------------------
  // 1. Parse command‐line arguments into inputPaths (requires exactly 2 paths)
  // --------------------------------------------------------------------------
  // if (argc != 3) {
  //   SPDLOG_CRITICAL("Usage: {} <naja-if-dir-1> <naja-if-dir-2>", argv[0]);
  //   return EXIT_FAILURE;
  // }
  printf("KEPLER FORMAL: Run.\n");

  size_t inputPathsIndex = 2;
  size_t inputLibraryIndex = 4;

  std::string formatType = argv[1];

  if (formatType == "-snl") {
    inputFormatType = FormatType::SNL;
  } else if (formatType == "-verilog") {
    inputFormatType = FormatType::VERILOG;
  } else {
    SPDLOG_CRITICAL("Unrecognized input format type: {}", formatType);
    return EXIT_FAILURE;
  }

  std::vector<std::string> inputPaths;
  for (int i = inputPathsIndex; i < argc; ++i) {
    inputPaths.emplace_back(argv[i]);
  }
  printf("number of library files: %zu\n", inputPaths.size() - 2);

  // --------------------------------------------------------------------------
  // 2. Load two netlists via Cap’n Proto
  // --------------------------------------------------------------------------
  // naja::NajaPerf::Scope scope("Parsing SNL format");
  // const auto t0 = steady_clock::now();
  // Load liberty
  NLUniverse::create();
  NLDB* db0 = nullptr;
  bool primitivesAreLoaded = false;
  if (inputPaths.size() > 2) {
    db0 = NLDB::create(NLUniverse::get());
    auto primitivesLibrary =
        NLLibrary::create(db0, NLLibrary::Type::Primitives, NLName("PRIMS"));
    SNLLibertyConstructor constructor(primitivesLibrary);
    for (size_t i = inputLibraryIndex; i < argc; ++i) {
      printf("Loading liberty file: %s\n", argv[i]);
      constructor.construct(argv[i]);
    }
    primitivesAreLoaded = true;
  }
  if (inputFormatType == FormatType::VERILOG) {
    auto designLibrary = NLLibrary::create(db0, NLName("DESIGN"));
    SNLVRLConstructor constructor(designLibrary);
    constructor.construct(argv[inputPathsIndex]);
    auto top = SNLUtils::findTop(designLibrary);
    if (top) {
      db0->setTopDesign(top);
      SPDLOG_INFO("Found top design: " + top->getString());
    } else {
      SPDLOG_ERROR("No top design was found after parsing verilog");
    }
  } else if (inputFormatType == FormatType::SNL) {
    printf("Loading SNL file: %s\n", argv[inputPathsIndex]);
    db0 = SNLCapnP::load(argv[inputPathsIndex], primitivesAreLoaded);
  } else {
    SPDLOG_CRITICAL("Unrecognized input format type: {}", formatType);
    return EXIT_FAILURE;
  }

  // get db0 top
  auto top0 = db0->getTopDesign();
  db0->setID(2);  // Increment ID to avoid conflicts
  NLDB* db1 = nullptr;

  // Increment ID to avoid conflicts
  if (inputPaths.size() > 2) {
    db1 = NLDB::create(NLUniverse::get());
    db1->setID(1);
    auto primitivesLibrary =
        NLLibrary::create(db1, NLLibrary::Type::Primitives, NLName("PRIMS"));
    SNLLibertyConstructor constructor(primitivesLibrary);
    for (size_t i = inputLibraryIndex; i < argc; ++i) {
      constructor.construct(argv[i]);
    }
  }
  if (inputFormatType == FormatType::VERILOG) {
    auto designLibrary = NLLibrary::create(db1, NLName("DESIGN"));
    SNLVRLConstructor constructor(designLibrary);
    constructor.construct(argv[inputPathsIndex + 1]);
    auto top = SNLUtils::findTop(designLibrary);
    if (top) {
      db1->setTopDesign(top);
      SPDLOG_INFO("Found top design: " + top->getString());
    } else {
      SPDLOG_ERROR("No top design was found after parsing verilog");
    }
  } else if (inputFormatType == FormatType::SNL) {
    printf("Loading SNL file: %s\n", argv[inputPathsIndex + 1]);
    db1 = SNLCapnP::load(argv[inputPathsIndex + 1], primitivesAreLoaded);
  } else {
    SPDLOG_CRITICAL("Unrecognized input format type: {}", formatType);
    return EXIT_FAILURE;
  }
  // get db1 top
  auto top1 = db1->getTopDesign();

  // --------------------------------------------------------------------------
  // 4. Hand off to the rest of the editing/analysis workflow
  // --------------------------------------------------------------------------
  /*try {
    KEPLER_FORMAL::MiterStrategy MiterS(top0, top1);
    if (MiterS.run()) {
      SPDLOG_INFO("Miter strategy succeeded: outputs are identical.");
    } else {
      SPDLOG_INFO("Miter strategy failed: outputs differ.");
    }
  } catch (const std::exception& e) {
    SPDLOG_ERROR("Workflow failed: {}", e.what());
    return EXIT_FAILURE;
  }*/
  KEPLER_FORMAL::MiterStrategy MiterS(top0, top1);
  MiterS.run();

  return EXIT_SUCCESS;
}
