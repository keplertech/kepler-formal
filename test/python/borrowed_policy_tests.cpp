// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "KeplerBorrowedDesigns.h"
#include "NLUniverse.h"
#include "NLLibrary.h"
#include "SNLDesign.h"
#include "SNLInstance.h"
#include "SNLInstTerm.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "latch/LatchSupportOptions.h"

namespace {
using namespace KEPLER_FORMAL;
using namespace naja::NL;

void check(bool value, const char* detail) {
  if (!value) throw std::runtime_error(detail);
}

void runTests(const std::filesystem::path& directory) {
  auto* universe = NLUniverse::create();
  auto* database = NLDB::create(universe);
  auto* library = NLLibrary::create(database, NLName("designs"));
  auto* top = SNLDesign::create(library, NLName("top"));
  auto* net = SNLScalarNet::create(top, NLName("net"));
  SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a"))->setNet(net);
  SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"))->setNet(net);
  BorrowedDesignOptions options;
  options.mode = BorrowedVerificationMode::SEC;
  options.secEncoding = SEC::SecEncoding::Binary;
  options.secEngine = SEC::SecEngine::KInduction;
  options.logFile = (directory / "run.log").string();
  RunResult result;

  // An incomplete ambient event contract would make any extraction unsupported
  // if borrowed verification accidentally inherited it.
  SEC::LATCH::SupportOptions ambient;
  ambient.enabled = true;
  ambient.workers = 3;
  SEC::LATCH::ScopedSupportOptions scope(ambient);
  check(verifyBorrowedDesigns(top, top, options, result) == 0 &&
            result.status == RunStatus::Equivalent,
        "borrowed verification inherited ambient event semantics");
  check(SEC::LATCH::supportOptions().enabled &&
            !SEC::LATCH::supportOptions().initialInputs.has_value() &&
            SEC::LATCH::supportOptions().workers == 3,
        "borrowed verification failed to restore ambient event contract");

  auto* primitives = NLLibrary::create(database, NLLibrary::Type::Primitives, NLName("prims"));
  auto* unknown = SNLDesign::create(primitives, SNLDesign::Type::Primitive, NLName("OPAQUE"));
  auto* d = SNLScalarTerm::create(unknown, SNLTerm::Direction::Input, NLName("D"));
  auto* q = SNLScalarTerm::create(unknown, SNLTerm::Direction::Output, NLName("Q"));
  auto* instance = SNLInstance::create(top, unknown, NLName("unused"));
  instance->getInstTerm(d)->setNet(net);
  instance->getInstTerm(q)->setNet(SNLScalarNet::create(top, NLName("unused_net")));
  Config::setErrorOnOpaque(true);
  check(verifyBorrowedDesigns(top, top, options, result) == 0 &&
            result.status == RunStatus::Equivalent,
        "borrowed default inherited caller's strict opaque policy");
  check(Config::getErrorOnOpaque(), "borrowed default did not restore strict opaque policy");
  Config::setErrorOnOpaque(false);
  options.errorOnOpaque = true;
  check(verifyBorrowedDesigns(top, top, options, result) != 0 &&
            result.status == RunStatus::Unsupported &&
            result.reason.find("unused") != std::string::npos,
        "borrowed strict policy ignored disconnected opacity");
  check(!Config::getErrorOnOpaque(), "borrowed strict policy leaked to caller");
  check(SEC::LATCH::supportOptions().enabled,
        "unsupported borrowed run failed to restore ambient event scope");
  options.mode = BorrowedVerificationMode::LEC;
  check(verifyBorrowedDesigns(top, top, options, result) != 0 &&
            result.status == RunStatus::Error,
        "borrowed LEC accepted SEC-only strict opaque policy");
  check(!Config::getErrorOnOpaque() && SEC::LATCH::supportOptions().enabled,
        "error path leaked borrowed policy state");
  universe->destroy();
}
}  // namespace

int main() {
  const auto directory = std::filesystem::temp_directory_path() /
      ("kepler_borrowed_policy_" + std::to_string(
          std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directories(directory);
  try {
    runTests(directory);
    std::filesystem::remove_all(directory);
    std::cout << "Borrowed policy isolation tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    if (auto* universe = naja::NL::NLUniverse::get()) universe->destroy();
    std::filesystem::remove_all(directory);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
