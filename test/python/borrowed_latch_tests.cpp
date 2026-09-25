// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <chrono>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "KeplerBorrowedDesigns.h"
#include "DNL.h"
#include "NLUniverse.h"
#include "NLLibrary.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLInstTerm.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"

namespace {
using namespace KEPLER_FORMAL;
using namespace naja::NL;
size_t checks = 0;
void check(bool value, const std::string& detail) {
  ++checks;
  if (!value) throw std::runtime_error(detail);
}

SNLDesign* latchModel(NLLibrary* primitives, const char* name, bool invert) {
  auto* cell = SNLDesign::create(primitives, SNLDesign::Type::Primitive, NLName(name));
  auto* data = SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("D"));
  auto* enable = SNLScalarTerm::create(cell, SNLTerm::Direction::Input, NLName("E"));
  auto* output = SNLScalarTerm::create(cell, SNLTerm::Direction::Output, NLName("Q"));
  SNLDesignModeling::SequentialModel model;
  model.kind = SNLDesignModeling::SequentialModel::Kind::Latch;
  model.clockedOn.root = model.clockedOn.addTerm(enable);
  SNLDesignModeling::SequentialState storage;
  storage.nextState.root = storage.nextState.addTerm(data);
  model.states.push_back(storage);
  SNLDesignModeling::BooleanExpression visible;
  visible.root = visible.addState(0);
  if (invert) visible.root = visible.addOperation(
      SNLDesignModeling::BooleanExpression::Operator::Not, {visible.root});
  model.outputs.push_back({output, visible});
  SNLDesignModeling::setSequentialModel(cell, model);
  return cell;
}

SNLDesign* top(NLLibrary* library, const char* name, SNLDesign* cell) {
  auto* design = SNLDesign::create(library, NLName(name));
  auto* instance = SNLInstance::create(design, cell, NLName("latch"));
  for (auto* pin : cell->getScalarTerms()) {
    auto* net = SNLScalarNet::create(design, pin->getName());
    SNLScalarTerm::create(design, pin->getDirection(), pin->getName())->setNet(net);
    instance->getInstTerm(pin)->setNet(net);
  }
  return design;
}

void run(const std::filesystem::path& directory) {
  auto* universe = NLUniverse::create();
  auto* database = NLDB::create(universe);
  auto* designs = NLLibrary::create(database, NLName("designs"));
  auto* primitives = NLLibrary::create(database, NLLibrary::Type::Primitives, NLName("primitives"));
  auto* positive = latchModel(primitives, "explicit_latch", false);
  auto* negative = latchModel(primitives, "explicit_inverted_latch", true);
  auto* first = top(designs, "first", positive);
  auto* second = top(designs, "second", positive);
  auto* different = top(designs, "different", negative);
  universe->setTopDesign(first);
  auto* callerGraph = naja::DNL::get();
  const auto firstReference = first->getReference();
  const auto secondReference = second->getReference();
  SEC::LATCH::SupportOptions ambient;
  ambient.enabled = true;
  ambient.workers = 7;
  SEC::LATCH::ScopedSupportOptions ambientScope(ambient);

  BorrowedDesignOptions options;
  options.mode = BorrowedVerificationMode::SEC;
  options.logFile = (directory / "latches.log").string();
  RunResult result;
  const auto unchanged = [&] {
    check(universe->getTopDesign() == first && naja::DNL::get() == callerGraph,
          "borrowed run changed the caller's selected design or DNL");
    check(universe->getSNLDesign(firstReference) == first &&
          universe->getSNLDesign(secondReference) == second &&
          first->getInstances().size() == 1 && second->getInstances().size() == 1,
          "borrowed run modified or deleted source designs");
    check(SEC::LATCH::supportOptions().enabled && SEC::LATCH::supportOptions().workers == 7 &&
          !SEC::LATCH::supportOptions().initialInputs,
          "borrowed run leaked its event contract");
  };
  verifyBorrowedDesigns(first, second, options, result);
  check(result.status == RunStatus::Unsupported && result.coveredOutputs == 0,
        "default did not preserve opaque-latch behavior");
  unchanged();

  check(options.latchSupport.enabled, "latch support is not enabled by default");
  options.latchSupport.enabled = false;
  verifyBorrowedDesigns(first, second, options, result);
  check(result.status == RunStatus::Unsupported && result.coveredOutputs == 0,
        "explicit false did not preserve opaque-latch behavior");
  unchanged();
  options.latchSupport = {};
  options.latchSupport.inputChanges = LatchInputChanges::Single;
  options.latchSupport.initialInputs = false;
  options.latchSupport.initialStorage = false;
  options.latchSupport.workers = 2;
  for (auto engine : {SEC::SecEngine::Pdr, SEC::SecEngine::KInduction, SEC::SecEngine::Imc}) {
    for (auto encoding : {SEC::SecEncoding::Binary, SEC::SecEncoding::DualRailSteady}) {
      options.secEngine = engine;
      options.secEncoding = encoding;
      check(verifyBorrowedDesigns(first, second, options, result) == 0 &&
            result.status == RunStatus::Equivalent && result.coveredOutputs == 1 && result.totalOutputs == 1,
            "explicit borrowed latch contract failed self equivalence: " + result.reason);
      unchanged();
    }
  }
  check(verifyBorrowedDesigns(first, different, options, result) != 0 &&
        result.status == RunStatus::Different, "different latch outputs incorrectly proved equivalent");
  unchanged();
  options.latchSupport.inputChanges = LatchInputChanges::Any;
  verifyBorrowedDesigns(first, second, options, result);
  check(result.coveredOutputs == 0, "order-dependent any-change latch should remain opaque");
  unchanged();
  options.latchSupport.inputChanges = LatchInputChanges::Single;
  options.latchSupport.initialStorage.reset();
  check(verifyBorrowedDesigns(first, second, options, result) != 0 &&
        result.status == RunStatus::Error && result.reason.find("requires explicit") != std::string::npos,
        "incomplete event contract accepted");
  unchanged();
  options.latchSupport.initialStorage = false;
  options.latchSupport.enabled = false;
  check(verifyBorrowedDesigns(first, second, options, result) != 0 &&
        result.status == RunStatus::Error && result.reason.find("requires latch_support") != std::string::npos,
        "tuning enabled the master gate");
  unchanged();
  options.latchSupport.enabled = true;
  options.mode = BorrowedVerificationMode::LEC;
  check(verifyBorrowedDesigns(first, second, options, result) != 0 && result.status == RunStatus::Error,
        "LEC accepted latch event semantics");
  unchanged();
  options.mode = BorrowedVerificationMode::SEC;
  options.setAsBoundary = {{"latch", "latch"}};
  check(verifyBorrowedDesigns(first, second, options, result) != 0 && result.status == RunStatus::Error &&
        result.reason.find("complete top interface") != std::string::npos, "event mode accepted leaf boundaries");
  unchanged();
  options.setAsBoundary.clear();
  options.latchSupport = {};
  verifyBorrowedDesigns(first, second, options, result);
  check(result.coveredOutputs == 0, "enabled call contaminated the following default call");
  unchanged();
  naja::DNL::destroy();
  universe->destroy();
}
}  // namespace

int main() {
  const auto directory = std::filesystem::temp_directory_path() /
      ("kepler_borrowed_latches_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directories(directory);
  try {
    run(directory);
    std::filesystem::remove_all(directory);
    std::cout << "Borrowed latch integration: " << checks << " checks passed\n";
    return 0;
  } catch (const std::exception& error) {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) universe->destroy();
    std::filesystem::remove_all(directory);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
