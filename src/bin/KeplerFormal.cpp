// File: src/apps/naja_edit/NajaEdit.cpp

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

int main(int argc, char** argv) {
  printf("KEPLER FORMAL: Run.\n");
  using namespace std::chrono;

  // --------------------------------------------------------------------------
  // 1. Parse command‐line arguments into inputPaths (requires exactly 2 paths)
  // --------------------------------------------------------------------------
  // if (argc != 3) {
  //   SPDLOG_CRITICAL("Usage: {} <naja-if-dir-1> <naja-if-dir-2>", argv[0]);
  //   return EXIT_FAILURE;
  // }

  std::vector<std::string> inputPaths;
  for (int i = 1; i < argc; ++i) {
    inputPaths.emplace_back(argv[i]);
  }

  // --------------------------------------------------------------------------
  // 2. Load two netlists via Cap’n Proto
  // --------------------------------------------------------------------------
  //naja::NajaPerf::Scope scope("Parsing SNL format");
  //const auto t0 = steady_clock::now();
  // Load liberty
  NLUniverse::create();
  auto db0 = NLDB::create(NLUniverse::get());
  if (inputPaths.size() > 2) {
    auto primitivesLibrary = NLLibrary::create(db0, NLLibrary::Type::Primitives,
                                               NLName("PRIMS"));
    SNLLibertyConstructor constructor(primitivesLibrary);
    for (size_t i = 2; i < inputPaths.size(); ++i) {
      constructor.construct(inputPaths[i]);
    }
  }
  db0 = SNLCapnP::load(inputPaths[0], true);
  // get db0 top
  auto top0 = db0->getTopDesign();
  db0->setID(2);  // Increment ID to avoid conflicts
  auto db1 = NLDB::create(NLUniverse::get());
  
  db1->setID(1);  // Increment ID to avoid conflicts
  if (inputPaths.size() > 2) {
    auto primitivesLibrary = NLLibrary::create(db1, NLLibrary::Type::Primitives,
                                               NLName("PRIMS"));
    SNLLibertyConstructor constructor(primitivesLibrary);
    for (size_t i = 2; i < inputPaths.size(); ++i) {
      constructor.construct(inputPaths[i]);
    }
  }
  db1 = SNLCapnP::load(inputPaths[1], true);
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
