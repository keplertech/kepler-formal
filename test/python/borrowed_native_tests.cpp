// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "KeplerBorrowedDesigns.h"
#include "BoolExpr.h"
#include "DNL.h"
#include "NLDB.h"
#include "NLDB0.h"
#include "NLLibrary.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLInstParameter.h"
#include "SNLInstTerm.h"
#include "SNLParameter.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "Tree2BoolExpr.h"
#include <spdlog/spdlog.h>

using namespace KEPLER_FORMAL;
using namespace naja::NL;

namespace {

void check(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

SNLDesign* makeWire(NLLibrary* library, const std::string& name,
                    const char* inputName = "I") {
  auto* design = SNLDesign::create(library, NLName(name));
  auto* input = SNLScalarTerm::create(design, SNLTerm::Direction::Input,
                                     NLName(inputName));
  auto* output = SNLScalarTerm::create(design, SNLTerm::Direction::Output,
                                      NLName("O"));
  auto* net = SNLScalarNet::create(design, NLName("wire"));
  input->setNet(net);
  output->setNet(net);
  return design;
}

SNLDesign* makeHierarchy(NLLibrary* library, const std::string& name) {
  auto* model = makeWire(library, name + "_model");
  auto* design = SNLDesign::create(library, NLName(name));
  auto* input = SNLScalarTerm::create(design, SNLTerm::Direction::Input, NLName("I"));
  auto* output = SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName("O"));
  auto* instance = SNLInstance::create(design, model, NLName("child"));
  auto* inputNet = SNLScalarNet::create(design, NLName("input"));
  auto* outputNet = SNLScalarNet::create(design, NLName("output"));
  input->setNet(inputNet);
  output->setNet(outputNet);
  instance->getInstTerm(model->getScalarTerm(NLName("I")))->setNet(inputNet);
  instance->getInstTerm(model->getScalarTerm(NLName("O")))->setNet(outputNet);
  return design;
}

SNLDesign* makeConstant(NLLibrary* library, const std::string& name, bool value) {
  auto* design = SNLDesign::create(library, NLName(name));
  auto* output = SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName("O"));
  auto* net = SNLScalarNet::create(design, NLName("constant"));
  net->setType(value ? SNLNet::Type::Assign1 : SNLNet::Type::Assign0);
  output->setNet(net);
  return design;
}

std::pair<SNLDesign*, SNLInstParameter*> makeMutableLut(
    NLLibrary* library, NLLibrary* primitives, const std::string& name) {
  auto* model = SNLDesign::create(primitives, SNLDesign::Type::Primitive,
                                  NLName(name + "_lut"));
  auto* modelInput = SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("A"));
  auto* modelOutput = SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Y"));
  auto* init = SNLParameter::create(model, NLName("INIT"),
                                   SNLParameter::Type::Binary, "2'h2");
  SNLDesignModeling::setTruthTableFromParameter(model, modelOutput, {modelInput}, init);
  auto* design = SNLDesign::create(library, NLName(name));
  auto* input = SNLScalarTerm::create(design, SNLTerm::Direction::Input, NLName("I"));
  auto* output = SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName("O"));
  auto* otherOutput = SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName("O2"));
  auto* instance = SNLInstance::create(design, model, NLName("lut"));
  auto* parameter = SNLInstParameter::create(instance, init, "2'h2");
  auto* inputNet = SNLScalarNet::create(design, NLName("input"));
  auto* outputNet = SNLScalarNet::create(design, NLName("output"));
  input->setNet(inputNet);
  output->setNet(outputNet);
  otherOutput->setNet(outputNet);
  instance->getInstTerm(modelInput)->setNet(inputNet);
  instance->getInstTerm(modelOutput)->setNet(outputNet);
  return {design, parameter};
}

SNLDesign* makeUnresetFlop(NLLibrary* library, const std::string& name,
                          bool addCombinationalOutput) {
  auto* design = SNLDesign::create(library, NLName(name));
  auto* input = SNLScalarTerm::create(design, SNLTerm::Direction::Input, NLName("I"));
  auto* clock = SNLScalarTerm::create(design, SNLTerm::Direction::Input, NLName("CLK"));
  auto* output = SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName("Q"));
  auto* instance = SNLInstance::create(design, NLDB0::getDFF(), NLName("ff"));
  auto* inputNet = SNLScalarNet::create(design, NLName("input"));
  auto* clockNet = SNLScalarNet::create(design, NLName("clock"));
  auto* outputNet = SNLScalarNet::create(design, NLName("output"));
  input->setNet(inputNet);
  clock->setNet(clockNet);
  output->setNet(outputNet);
  instance->getInstTerm(NLDB0::getDFFData())->setNet(inputNet);
  instance->getInstTerm(NLDB0::getDFFClock())->setNet(clockNet);
  instance->getInstTerm(NLDB0::getDFFOutput())->setNet(outputNet);
  if (addCombinationalOutput) {
    SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName("O"))->setNet(inputNet);
  }
  return design;
}

void runTests() {
  auto* universe = NLUniverse::create();
  auto* db0 = NLDB::create(universe);
  auto* db1 = NLDB::create(universe);
  auto* library0 = NLLibrary::create(db0, NLName("first"));
  auto* library1 = NLLibrary::create(db1, NLName("second"));
  auto* primitives0 = NLLibrary::create(db0, NLLibrary::Type::Primitives, NLName("primitives0"));
  auto* primitives1 = NLLibrary::create(db1, NLLibrary::Type::Primitives, NLName("primitives1"));
  auto* first = makeHierarchy(library0, "first");
  auto* second = makeHierarchy(library1, "second");
  auto* anchor = makeWire(library0, "anchor");
  auto* mismatch = makeWire(library1, "mismatch", "OTHER");
  auto* zero = makeConstant(library0, "zero", false);
  auto* otherZero = makeConstant(library1, "other_zero", false);
  auto mutableFirst = makeMutableLut(library0, primitives0, "mutable_first");
  auto mutableSecond = makeMutableLut(library1, primitives1, "mutable_second");
  auto* partialFirst = makeUnresetFlop(library0, "partial_first", true);
  auto* partialSecond = makeUnresetFlop(library1, "partial_second", true);
  auto* inconclusiveFirst = makeUnresetFlop(library0, "inconclusive_first", false);
  auto* inconclusiveSecond = makeUnresetFlop(library1, "inconclusive_second", false);
  universe->setTopDesign(anchor);
  db1->setTopDesign(nullptr);
  auto* savedDnl = naja::DNL::get();

  auto* savedExpression = BoolExpr::Var(987654);
  Tree2BoolExpr::iso2boolExpr_.insert({123456, savedExpression});
  Config::setSolverType(Config::SolverType::CADICAL);
  Config::setReportSkippedPOs(true);
  auto savedLogger = spdlog::default_logger();

  // DNL assigns these metadata fields while traversing the borrowed graphs.
  auto* firstInput = first->getScalarTerm(NLName("I"));
  auto* child = first->getInstance(NLName("child"));
  auto* modelInput = child->getModel()->getScalarTerm(NLName("I"));
  firstInput->setOrderID(100);
  child->setOrderID(101);
  modelInput->setOrderID(102);
  const auto firstRevision = first->getRevisionCount();
  const auto secondRevision = second->getRevisionCount();

  auto checkState = [&]() {
    check(NLUniverse::get() == universe, "universe identity changed");
    check(universe->getTopDB() == db0 && universe->getTopDesign() == anchor,
          "caller top selection changed");
    check(db0->getTopDesign() == anchor && db1->getTopDesign() == nullptr,
          "a per-DB top selection changed");
    check(naja::DNL::isCreated() && naja::DNL::get() == savedDnl,
          "caller DNL was replaced or destroyed");
    check(firstInput->getOrderID() == 100 && child->getOrderID() == 101 &&
              modelInput->getOrderID() == 102, "caller ordering metadata changed");
    check(first->getRevisionCount() == firstRevision &&
              second->getRevisionCount() == secondRevision, "design was modified");
    check(BoolExpr::Var(987654) == savedExpression, "caller expression cache changed");
    auto found = Tree2BoolExpr::iso2boolExpr_.find(123456);
    check(found != Tree2BoolExpr::iso2boolExpr_.end() && found->second == savedExpression,
          "caller expression lookup map changed");
    check(Config::getSolverType() == Config::SolverType::CADICAL &&
              Config::getReportSkippedPOs(), "caller configuration changed");
    check(spdlog::default_logger() == savedLogger, "caller logger changed");
    check(Config::getVerificationGeneration() == 0, "caller cache generation changed");
  };

  BorrowedDesignOptions options;
  options.logFile = "borrowed_native_lec.log";
  RunResult result;
  for (int repeat = 0; repeat < 2; ++repeat) {
    check(verifyBorrowedDesigns(first, second, options, result) == 0 &&
              result.status == RunStatus::Equivalent,
          "hierarchical LEC failed: " + result.reason);
    checkState();
    mutableSecond.second->setValue("2'h1");
    check(verifyBorrowedDesigns(mutableFirst.first, mutableSecond.first, options, result) == 0 &&
              result.status == RunStatus::Different,
          "different LEC verdict/exit code changed: " + result.reason);
    mutableSecond.second->setValue("2'h2");
    checkState();
  }
  check(verifyBorrowedDesigns(first, mismatch, options, result) == 1 &&
            result.status == RunStatus::Error &&
            result.reason.find("boundary mismatch") != std::string::npos,
        "boundary mismatch did not return a structured error: " + result.reason);
  checkState();

  options.mode = BorrowedVerificationMode::SEC;
  options.secEncoding = SEC::SecEncoding::Binary;
  options.maxK = 2;
  options.logFile = "borrowed_native_sec.log";
  check(verifyBorrowedDesigns(mutableFirst.first, mutableSecond.first, options, result) == 0 &&
            result.status == RunStatus::Equivalent && result.coveredOutputs == 2 &&
            result.totalOutputs == 2 && result.provenOutputs == 2,
        "equivalent SEC verdict/counts changed: " + result.reason);
  checkState();
  mutableSecond.second->setValue("2'h1");
  check(verifyBorrowedDesigns(mutableFirst.first, mutableSecond.first, options, result) == 3 &&
            result.status == RunStatus::Different,
        "different SEC verdict/exit code changed: " + result.reason);
  mutableSecond.second->setValue("2'h2");
  checkState();

  check(verifyBorrowedDesigns(zero, otherZero, options, result) == 2 &&
            result.status == RunStatus::Unsupported && !result.reason.empty(),
        "unsupported SEC extraction was not returned as a verdict");
  checkState();

  options.secEngine = SEC::SecEngine::KInduction;
  options.secEncoding = SEC::SecEncoding::Binary;
  options.maxK = 1;
  check(verifyBorrowedDesigns(partialFirst, partialSecond, options, result) == 1 &&
            result.status == RunStatus::PartiallyProved &&
            result.provenOutputs == 1 && result.coveredOutputs == 1 &&
            result.totalOutputs == 2 && !result.skippedObservedOutputs.empty(),
        "partial SEC verdict/counts changed: " + result.reason);
  checkState();
  options.secEncoding = SEC::SecEncoding::DualRailSteady;
  options.maxK = 0;
  check(verifyBorrowedDesigns(inconclusiveFirst, inconclusiveSecond, options, result) == 2 &&
            result.status == RunStatus::Inconclusive && result.provenOutputs == 0,
        "inconclusive SEC verdict/exit code changed: " + result.reason);
  checkState();

  // Editing an instance parameter keeps its pointer and topology stable, but
  // changes its primitive truth table. Neither worker-local extraction caches
  // nor the multi-output KI base-prefix cache may reuse the previous proof.
  for (auto mode : {BorrowedVerificationMode::LEC, BorrowedVerificationMode::SEC}) {
    options.mode = mode;
    options.secEngine = SEC::SecEngine::KInduction;
    options.secEncoding = SEC::SecEncoding::DualRailSteady;
    options.maxK = 1;
    for (bool invert : {false, true, false, true}) {
      mutableSecond.second->setValue(invert ? "2'h1" : "2'h2");
      const auto expected = invert ? RunStatus::Different : RunStatus::Equivalent;
      verifyBorrowedDesigns(mutableFirst.first, mutableSecond.first, options, result);
      check(result.status == expected,
            "edited live LUT reused a stale proof: " + result.reason);
      checkState();
    }
  }

  // A universe with no top and no DNL is valid; restoring it must not call
  // NLUniverse::setTopDesign(nullptr), which dereferences its argument.
  naja::DNL::destroy();
  db0->setTopDesign(nullptr);
  db1->setTopDesign(nullptr);
  universe->setTopDB(nullptr);
  options.mode = BorrowedVerificationMode::LEC;
  check(verifyBorrowedDesigns(first, second, options, result) == 0 &&
            result.status == RunStatus::Equivalent, "no-top verification failed");
  check(universe->getTopDB() == nullptr && universe->getTopDesign() == nullptr &&
            db0->getTopDesign() == nullptr && db1->getTopDesign() == nullptr &&
            !naja::DNL::isCreated(), "absent caller top/DNL was not restored");
  check(verifyBorrowedDesigns(first, mismatch, options, result) == 1 &&
            result.status == RunStatus::Error, "no-top boundary error was lost");
  check(universe->getTopDB() == nullptr && !naja::DNL::isCreated(),
        "exception failed to restore absent top/DNL");

  Tree2BoolExpr::iso2boolExpr_.clear();
  BoolExprCache::destroy();
  universe->destroy();
}

}  // namespace

int main() {
  try {
    runTests();
    std::cout << "Borrowed-design native tests passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Borrowed-design native test failed: " << error.what() << '\n';
    naja::DNL::destroy();
    if (NLUniverse::get()) {
      NLUniverse::get()->destroy();
    }
    return 1;
  }
}
