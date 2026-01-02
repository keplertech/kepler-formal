// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "SNLLogicCloud.h"
#include <tbb/tbb_allocator.h>
#include <cassert>
#include "SNLDesignModeling.h"
#include "tbb/concurrent_vector.h"
#include "tbb/enumerable_thread_specific.h"

typedef std::pair<
    std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>,
    size_t>
    IterationInputsETSPair;
tbb::enumerable_thread_specific<IterationInputsETSPair>
    currentIterationInputsETS;
tbb::enumerable_thread_specific<IterationInputsETSPair> newIterationInputsETS;

tbb::concurrent_vector<IterationInputsETSPair*>
    currentIterationInputsETSvector =
        tbb::concurrent_vector<IterationInputsETSPair*>(40, nullptr);

tbb::concurrent_vector<IterationInputsETSPair*> newIterationInputsETSvector =
    tbb::concurrent_vector<IterationInputsETSPair*>(40, nullptr);

void initCurrentIterationInputsETS() {
  // LCOV_EXCL_START
  if (currentIterationInputsETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
    for (size_t i = currentIterationInputsETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
      currentIterationInputsETSvector.push_back(nullptr);
    }
  }
  // LCOV_EXCL_STOP
  if (currentIterationInputsETSvector
          [tbb::this_task_arena::current_thread_index()] == nullptr) {
    currentIterationInputsETSvector
        [tbb::this_task_arena::current_thread_index()] =
            &currentIterationInputsETS.local();
  }
}

IterationInputsETSPair& getCurrentIterationInputsETS() {
  return *currentIterationInputsETSvector
      [tbb::this_task_arena::current_thread_index()];
}

void initNewIterationInputsETS() {
  // LCOV_EXCL_START
  if (newIterationInputsETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
    for (size_t i = newIterationInputsETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
      newIterationInputsETSvector.push_back(nullptr);
    }
  }
  // LCOV_EXCL_STOP
  if (newIterationInputsETSvector
          [tbb::this_task_arena::current_thread_index()] == nullptr) {
    newIterationInputsETSvector[tbb::this_task_arena::current_thread_index()] =
        &newIterationInputsETS.local();
  }
}

IterationInputsETSPair& getNewIterationInputsETS() {
  return *newIterationInputsETSvector
      [tbb::this_task_arena::current_thread_index()];
}

void clearCurrentIterationInputsETS() {
  auto& currentIterationInputs = getCurrentIterationInputsETS();
  currentIterationInputs.second = 0;
}

void pushBackCurrentIterationInputsETS(naja::DNL::DNLID input) {
  auto& currentIterationInputs = getCurrentIterationInputsETS();
  auto& vec = currentIterationInputs.first;
  auto& sz = currentIterationInputs.second;
  if (vec.size() > sz) {
    vec[sz] = input;
    sz++;
    return;
  }
  vec.push_back(input);
  sz++;
}

size_t sizeOfCurrentIterationInputsETS() {
  return getCurrentIterationInputsETS().second;
}

void copyCurrentIterationInputsETS(std::vector<naja::DNL::DNLID>& res) {
  res.clear();
  auto& current = getCurrentIterationInputsETS();
  for (size_t i = 0; i < current.second; i++) {
    res.push_back(current.first[i]);
  }
}

void clearNewIterationInputsETS() {
  auto& newIterationInputs = getNewIterationInputsETS();
  newIterationInputs.second = 0;
}

void pushBackNewIterationInputsETS(naja::DNL::DNLID input) {
  auto& newIterationInputs = getNewIterationInputsETS();
  auto& vec = newIterationInputs.first;
  auto& sz = newIterationInputs.second;
  if (vec.size() > sz) {
    vec[sz] = input;
    sz++;
    return;
  }
  vec.push_back(input);
  sz++;
}

bool emptyNewIterationInputsETS() {
  return getNewIterationInputsETS().second == 0;
}

size_t sizeOfNewIterationInputsETS() {
  return getNewIterationInputsETS().second;
}

void copyNewIterationInputsETStoCurrent() {
  clearCurrentIterationInputsETS();
  auto& newIterationInputs = getNewIterationInputsETS();
  for (size_t i = 0; i < newIterationInputs.second; i++) {
    pushBackCurrentIterationInputsETS(newIterationInputs.first[i]);
  }
  assert(sizeOfCurrentIterationInputsETS() == sizeOfNewIterationInputsETS());
}

// #define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace KEPLER_FORMAL;
using namespace naja::DNL;

bool SNLLogicCloud::isInput(naja::DNL::DNLID termID) {
  return PIs_[termID];
}

bool SNLLogicCloud::isOutput(naja::DNL::DNLID termID) {
  return POs_[termID];
}

void SNLLogicCloud::compute() {
  // std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>
  // newIterationInputs;
  initNewIterationInputsETS();
  initCurrentIterationInputsETS();
  clearNewIterationInputsETS();
  clearCurrentIterationInputsETS();
  DEBUG_LOG("---- Begin!!\n");
  if (dnl_.getDNLTerminalFromID(seedOutputTerm_).isTopPort() ||
      isOutput(seedOutputTerm_)) {
    auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(seedOutputTerm_).getIsoID());
    // LCOV_EXCL_START
    if (iso.getDrivers().size() > 1) {
      
      for (auto driver : iso.getDrivers()) {
        DEBUG_LOG("Driver: %s\n", dnl_.getDNLTerminalFromID(driver)
                                      .getSnlBitTerm()
                                      ->getName()
                                      .getString()
                                      .c_str());
      }
      throw std::runtime_error("Seed output term is not a single driver");
    } else if (iso.getDrivers().empty()) {
      std::string termName = dnl_.getDNLTerminalFromID(seedOutputTerm_)
                                 .getSnlBitTerm()
                                 ->getName()
                                 .getString();
      std::string error =
          "Seed output term '" + termName + "' has no drivers";
      throw std::runtime_error(error);
    }
    // LCOV_EXCL_STOP
    auto driver = iso.getDrivers().front();
    auto inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
    if (isInput(driver)) {
      pushBackCurrentIterationInputsETS(driver);
      table_ = SNLTruthTableTree(inst.getID(), driver,
                                 SNLTruthTableTree::Node::Type::P);
      return;
    }
    DEBUG_LOG("Instance name: %s\n",
              inst.getSNLInstance()->getName().getString().c_str());
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Output) {
        pushBackNewIterationInputsETS(termID);
        DEBUG_LOG("Add input with id: %zu\n", termID);
      }
    }
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    table_ = SNLTruthTableTree(inst.getID(), driver);
    auto* model = const_cast<SNLDesign*>(inst.getSNLModel());
    assert(SNLDesignModeling::getTruthTable(model, 
                dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID())
            .isInitialized() &&
        "Truth table is not initialized");
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  } else {
    auto inst = dnl_.getDNLInstanceFromID(seedOutputTerm_);
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Output) {
        // newIterationInputs.push_back(termID);
        pushBackNewIterationInputsETS(termID);
        DEBUG_LOG("Add input with id: %zu\n", termID);
      }
    }
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    table_ = SNLTruthTableTree(inst.getID(), seedOutputTerm_);
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  }

  if (emptyNewIterationInputsETS()) {
    DEBUG_LOG("No inputs found for seed output term %zu\n", seedOutputTerm_);
    return;
  }

  bool reachedPIs = true;
  size_t size = sizeOfNewIterationInputsETS();
  for (size_t i = 0; i < size; i++) {
    if (!isInput(
            getNewIterationInputsETS().first
                [i]) /* && !isOutput(getNewIterationInputsETS().first[i])*/) {
      reachedPIs = false;
      break;
    }
  }

  std::set<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>> handledTerms;
  size_t iter = 0;

  while (!reachedPIs) {
    DEBUG_LOG("---iter %lu---\n", iter);
    DEBUG_LOG("Current iteration inputs size: %zu\n",
              sizeOfNewIterationInputsETS());
    copyNewIterationInputsETStoCurrent();

    clearNewIterationInputsETS();
    DEBUG_LOG("table size: %zu, currentIterationInputs_ size: %zu\n",
              table_.size(), sizeOfCurrentIterationInputsETS());

    std::vector<
        std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
        tbb::tbb_allocator<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>>>
        inputsToMerge;

    size_t sizeOfCurrentInputs = sizeOfCurrentIterationInputsETS();
    for (size_t i = 0; i < sizeOfCurrentInputs; i++) {
      auto input = getCurrentIterationInputsETS().first[i];
      if (isInput(input) /*|| isOutput(input)*/) {
        pushBackNewIterationInputsETS(input);
        DEBUG_LOG("Adding input id: %zu %s\n", input,
                  dnl_.getDNLTerminalFromID(input)
                      .getSnlBitTerm()
                      ->getName()
                      .getString()
                      .c_str());
        inputsToMerge.push_back(
            {naja::DNL::DNLID_MAX, input});  // Placeholder for PI/PO
        continue;
      }

      auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          dnl_.getDNLTerminalFromID(input).getIsoID());
      DEBUG_LOG("number of drivers: %zu\n", iso.getDrivers().size());

      for (auto driver : iso.getDrivers()) {
        DEBUG_LOG("Driver: %s\n", dnl_.getDNLTerminalFromID(driver)
                                      .getSnlBitTerm()
                                      ->getName()
                                      .getString()
                                      .c_str());
      }

      if (iso.getDrivers().size() >= 1) {
        assert(iso.getDrivers().size() <= 1 &&
               "Iso have more than one driver, not supported");
      } else if (iso.getDrivers().empty()) {
        assert(iso.getDrivers().size() == 1 &&
               "Iso have no drivers and more than one reader, not supported");
      }

      auto driver = iso.getDrivers().front();
      if (isInput(driver) /* || isOutput(driver)*/) {
        pushBackNewIterationInputsETS(driver);
        DEBUG_LOG(
            "- %lu After analyzing input %s(%lu), addings driver %s(%lu) is a "
            "primary input\n",
            iter,
            dnl_.getDNLTerminalFromID(input)
                .getSnlBitTerm()
                ->getName()
                .getString()
                .c_str(),
            input,
            dnl_.getDNLTerminalFromID(driver)
                .getSnlBitTerm()
                ->getName()
                .getString()
                .c_str(),
            driver);
        inputsToMerge.push_back(
            {naja::DNL::DNLID_MAX, driver});  // Placeholder for PI/PO
        continue;
      }

      auto inst = dnl_.getDNLInstanceFromID(
          dnl_.getDNLTerminalFromID(driver).getDNLInstance().getID());
      auto* model = const_cast<SNLDesign*>(inst.getSNLModel());
      // if (!model
      //          ->getTruthTable(dnl_.getDNLTerminalFromID(driver)
      //                              .getSnlBitTerm()
      //                              ->getOrderID())
      //          .isInitialized()) {
      //   DEBUG_LOG(
      //       "#####Truth table for instance %s is not initialized\n",
      //       inst.getSNLInstance()->getModel()->getName().getString().c_str());
      //   auto* model = const_cast<SNLDesign*>(inst.getSNLModel());
      //   assert(model
      //              ->getTruthTable(dnl_.getDNLTerminalFromID(driver)
      //                                  .getSnlBitTerm()
      //                                  ->getOrderID())
      //              .isInitialized() &&
      //          "Truth table for instance is not initialized");
      // }

      DEBUG_LOG("Adding driver id: %zu %s(%s)\n", driver,
                dnl_.getDNLTerminalFromID(driver)
                    .getSnlBitTerm()
                    ->getName()
                    .getString()
                    .c_str(),
                dnl_.getDNLTerminalFromID(driver)
                    .getSnlBitTerm()
                    ->getDesign()
                    ->getName()
                    .getString()
                    .c_str());
      inputsToMerge.push_back({inst.getID(), driver});

      for (DNLID termID = inst.getTermIndexes().first;
           termID <= inst.getTermIndexes().second; termID++) {
        const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Output) {
          if (handledTerms.find({driver, termID}) != handledTerms.end()) {
            DEBUG_LOG(
                "#### iter %lu 1 Term (%zu) %s of inst %s already handled, "
                "skipping\n",
                iter, input,
                naja::DNL::get()
                    ->getDNLTerminalFromID(input)
                    .getSnlBitTerm()
                    ->getName()
                    .getString()
                    .c_str(),
                naja::DNL::get()
                    ->getDNLTerminalFromID(input)
                    .getDNLInstance()
                    .getSNLModel()
                    ->getName()
                    .getString()
                    .c_str());
            continue;
          }
          handledTerms.insert({driver, termID});
          pushBackNewIterationInputsETS(termID);
        }
      }
    }

    if (inputsToMerge.empty()) {
      break;
    }

    DEBUG_LOG("--- Merging truth tables with %zu inputs\n",
              inputsToMerge.size());
    table_.concatFull(inputsToMerge);
    reachedPIs = true;
    size_t sizeOfNewInputs = sizeOfNewIterationInputsETS();
    for (size_t i = 0; i < sizeOfNewInputs; i++) {
      if (!isInput(getNewIterationInputsETS().first[i])) {
        reachedPIs = false;
        break;
      }
    }
    DEBUG_LOG("--- End of iteration %zu\n", iter);
    iter++;
  }

  copyNewIterationInputsETStoCurrent();
  copyCurrentIterationInputsETS(currentIterationInputs_);
  assert(currentIterationInputs_.size() == sizeOfCurrentIterationInputsETS());
  for (auto input : currentIterationInputs_) {
    assert(isInput(input));
  }
}
