// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "SNLLogicCloud.h"
#include <tbb/tbb_allocator.h>
#include <cassert>
#include "SNLDesignModeling.h"
#include "tbb/concurrent_vector.h"
#include "tbb/enumerable_thread_specific.h"
#include "SNLPath.h"
#include "Tree2BoolExpr.h"


// #define DEBUG_CHECKS
// #define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace KEPLER_FORMAL;
using namespace naja::DNL;
typedef std::pair<
    std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>,
    size_t>
    IterationInputsETSPair;

thread_local IterationInputsETSPair currentIterationInputsETS;

IterationInputsETSPair& getCurrentIterationInputsETS() {
  return currentIterationInputsETS;
}

thread_local IterationInputsETSPair newIterationInputsETS;

IterationInputsETSPair& getNewIterationInputsETS() {
  return newIterationInputsETS;
}

void clearCurrentIterationInputsETS() {
  auto& currentIterationInputs = getCurrentIterationInputsETS();
  currentIterationInputs.first.clear();
}

void pushBackCurrentIterationInputsETS(naja::DNL::DNLID input) {
  auto& currentIterationInputs = getCurrentIterationInputsETS();
  currentIterationInputs.first.emplace_back(input);
}

size_t sizeOfCurrentIterationInputsETS() {
  return getCurrentIterationInputsETS().first.size();
}

void copyCurrentIterationInputsETS(std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>& res) {
  res.clear();
  auto& current = getCurrentIterationInputsETS();
  res = std::move(current.first);
}

void clearNewIterationInputsETS() {
  auto& newIterationInputs = getNewIterationInputsETS();
  newIterationInputs.first.clear();
}

void pushBackNewIterationInputsETS(naja::DNL::DNLID input) {
  getNewIterationInputsETS().first.emplace_back(input);
}

bool emptyNewIterationInputsETS() {
  return getNewIterationInputsETS().first.empty();
}

size_t sizeOfNewIterationInputsETS() {
  return getNewIterationInputsETS().first.size();
}

void copyNewIterationInputsETStoCurrent() {
  auto& newIterationInputs = getNewIterationInputsETS();
  auto& currentIterationInputs = getCurrentIterationInputsETS();
  #ifdef DEBUG_CHECKS
  size_t newSize = newIterationInputs.first.size();
  #endif
  currentIterationInputs = std::move(newIterationInputs);
  #ifdef DEBUG_CHECKS
  assert(currentIterationInputs.first.size() == newSize &&
         "copyNewIterationInputsETStoCurrent: size mismatch after copy");
  #endif
}

thread_local std::pair<
    std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
                tbb::tbb_allocator<std::pair<naja::DNL::DNLID,
                                            naja::DNL::DNLID>>>,
    size_t>
    inputsToMergeETS;

std::pair<std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
                      tbb::tbb_allocator<std::pair<naja::DNL::DNLID,
                                                  naja::DNL::DNLID>>>, size_t>&
getInputsToMergeETS() {
  return inputsToMergeETS;
}

void clearInputsToMergeETS() {
  auto& inputsToMerge = getInputsToMergeETS();
  inputsToMerge.first.clear();
}

void pushBackInputsToMergeETS(
    const std::pair<naja::DNL::DNLID, naja::DNL::DNLID>& input) {
  getInputsToMergeETS().first.emplace_back(input);
}

size_t sizeOfInputsToMergeETS() {
  return getInputsToMergeETS().first.size();
}

// 2 level vector visited terms pair - 1st: termID, 2nd: termID
typedef std::vector<
    std::unordered_set<naja::DNL::DNLID, 
                                                 std::hash<naja::DNL::DNLID>,
                                                 std::equal_to<naja::DNL::DNLID>,
                                                 tbb::tbb_allocator<naja::DNL::DNLID>>,
    tbb::tbb_allocator<std::unordered_set<naja::DNL::DNLID, 
                                                 std::hash<naja::DNL::DNLID>,
                                                 std::equal_to<naja::DNL::DNLID>,
                                                 tbb::tbb_allocator<naja::DNL::DNLID>>>>
    VisitedTermsPairsVec;

thread_local VisitedTermsPairsVec visitedTermsPairsETS;

struct PairHash {
  size_t operator()(const std::pair<naja::DNL::DNLID,naja::DNL::DNLID>& p) const noexcept {
    // 64-bit combine; tweak for your DNLID type
    uint64_t a = static_cast<uint64_t>(p.first);
    uint64_t b = static_cast<uint64_t>(p.second);
    return (a * 11400714819323198485ull) ^ (b + 0x9e3779b97f4a7c15ull + (a<<6) + (a>>2));
    }
  };
  struct PairEq {
    bool operator()(const std::pair<naja::DNL::DNLID,naja::DNL::DNLID>& x, const std::pair<naja::DNL::DNLID,naja::DNL::DNLID>& y) const noexcept {
      return x.first == y.first && x.second == y.second;
    }
  };
   using HandledSet = std::unordered_set<
     std::pair<naja::DNL::DNLID,naja::DNL::DNLID>,
     PairHash, PairEq,
     tbb::tbb_allocator<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>>>;

thread_local HandledSet visitedTermsPairsETSSet;

void clearVisitedTermsPairsETS() {
  visitedTermsPairsETSSet.clear();
}

thread_local std::pair<naja::DNL::DNLID, naja::DNL::DNLID> tempPairETS;

bool isPairVisitedETS(naja::DNL::DNLID termA,
                              naja::DNL::DNLID termB) {
  tempPairETS.first = termA;
  tempPairETS.second = termB;
  if (!(visitedTermsPairsETSSet.insert(tempPairETS)).second) {
    return true;
  }
  return false;
}

bool SNLLogicCloud::isInput(naja::DNL::DNLID termID) {
  return PIs_[termID];
}

bool SNLLogicCloud::isOutput(naja::DNL::DNLID termID) {
  return POs_[termID];
}

void SNLLogicCloud::compute() {
  clearNewIterationInputsETS();
  clearCurrentIterationInputsETS();
  DEBUG_LOG("---- Begin!!\n");
  if (dnl_.getDNLTerminalFromID(seedOutputTerm_).isTopPort() ||
      isOutput(seedOutputTerm_)) {
    const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(seedOutputTerm_).getIsoID());
    // LCOV_EXCL_START
    if (iso.getDrivers().size() > 1) {
      
      for (const auto& driver : iso.getDrivers()) {
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
    const auto& driver = iso.getDrivers().front();
    auto& inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
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
    auto* model = inst.getSNLModel();
    assert(SNLDesignModeling::getTruthTable(model, 
                dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID())
            .isInitialized() &&
        "Truth table is not initialized");
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  } else {
    const auto& inst = dnl_.getDNLInstanceFromID(seedOutputTerm_);
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Output) {
        // newIterationInputs.emplace_back(termID);
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
  } // allocator for buckets

  // HandledSet handledTerms;
  // handledTerms.reserve(naja::DNL::get()->getDNLTerms().size() / 4);
  clearVisitedTermsPairsETS();
  size_t iter = 0;

  while (!reachedPIs) {
    DEBUG_LOG("---iter %lu---\n", iter);
    DEBUG_LOG("Current iteration inputs size: %zu\n",
              sizeOfNewIterationInputsETS());
    copyNewIterationInputsETStoCurrent();

    clearNewIterationInputsETS();
    DEBUG_LOG("table size: %zu, currentIterationInputs_ size: %zu\n",
              table_.size(), sizeOfCurrentIterationInputsETS());
    clearInputsToMergeETS();
    size_t sizeOfCurrentInputs = sizeOfCurrentIterationInputsETS();
    for (size_t i = 0; i < sizeOfCurrentInputs; i++) {
      const auto& input = getCurrentIterationInputsETS().first[i];
      if (isInput(input) /*|| isOutput(input)*/) {
        pushBackNewIterationInputsETS(input);
        DEBUG_LOG("Adding input id: %zu %s\n", input,
                  dnl_.getDNLTerminalFromID(input)
                      .getSnlBitTerm()
                      ->getName()
                      .getString()
                      .c_str());
        pushBackInputsToMergeETS(
            {naja::DNL::DNLID_MAX, input});  // Placeholder for PI/PO
        continue;
      }

      const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          dnl_.getDNLTerminalFromID(input).getIsoID());
      DEBUG_LOG("number of drivers: %zu\n", iso.getDrivers().size());

      for (const auto& driver : iso.getDrivers()) {
        DEBUG_LOG("Driver: %s\n", dnl_.getDNLTerminalFromID(driver)
                                      .getSnlBitTerm()
                                      ->getName()
                                      .getString()
                                      .c_str());
      }

      if (iso.getDrivers().size() >= 1) {
        // proper error with names of all the drivers
        // throw an error and separate names by comma
        if (iso.getDrivers().size() > 1) {
          std::vector<std::string> namesOfDrivers;
          for (auto dnlid : iso.getDrivers()) {
            auto driver = dnl_.getDNLTerminalFromID(dnlid);
            auto path = driver.getDNLInstance().getPath().getPathNames();
            std::string fullName;
            for (size_t i = 0; i < path.size(); i++) {
              fullName += path[i].getString();
              if (i != path.size() - 1) {
                // LCOV_EXCL_START
                fullName += ".";
                // LCOV_EXCL_STOP
              }
            }
            // add terminal name and bit
            std::string termName =
                driver.getSnlBitTerm()->getName().getString();
            fullName += "." + termName;
            // add bit
            fullName += std::to_string(driver.getSnlBitTerm()->getBit());
            namesOfDrivers.push_back(fullName);
          }
          std::string error = "Iso has multiple drivers: ";
          for (size_t i = 0; i < namesOfDrivers.size(); i++) {
            error += namesOfDrivers[i];
            if (i != namesOfDrivers.size() - 1) {
              error += ", ";
            }
          }
          throw std::runtime_error(error);
        }
      } else if (iso.getDrivers().empty()) {
        assert(iso.getDrivers().size() == 1 &&
               "Iso have no drivers and more than one reader, not supported");
      }
      const auto& driver = iso.getDrivers().front();
      
      if (isInput(driver) /* || isOutput(driver)*/
        || (Tree2BoolExpr::iso2boolExpr_.find(iso.getIsoID()) !=
            Tree2BoolExpr::iso2boolExpr_.end() && iter > 0)) {
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
        pushBackInputsToMergeETS(
            {naja::DNL::DNLID_MAX, driver});  // Placeholder for PI/PO
        continue;
      }

      const auto& inst = dnl_.getDNLInstanceFromID(
          dnl_.getDNLTerminalFromID(driver).getDNLInstance().getID());

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
      pushBackInputsToMergeETS({inst.getID(), driver});

      for (DNLID termID = inst.getTermIndexes().first;
           termID <= inst.getTermIndexes().second; termID++) {
        const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Output) {
          if (isPairVisitedETS(driver, termID)) {
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
          pushBackNewIterationInputsETS(termID);
        }
      }
    }

    if (sizeOfInputsToMergeETS() == 0) {
      break;
    }

    DEBUG_LOG("--- Merging truth tables with %zu inputs\n",
              sizeOfInputsToMergeETS());
    table_.concatFull(getInputsToMergeETS().first,
                      sizeOfInputsToMergeETS());
    reachedPIs = true;
    size_t sizeOfNewInputs = sizeOfNewIterationInputsETS();
    for (size_t i = 0; i < sizeOfNewInputs; i++) {
      auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          dnl_.getDNLTerminalFromID(
              getNewIterationInputsETS().first[i])
              .getIsoID());
      if (!isInput(getNewIterationInputsETS().first[i]) &&
        (Tree2BoolExpr::iso2boolExpr_.find(
             dnl_.getDNLTerminalFromID(
                 getNewIterationInputsETS().first[i])
                 .getIsoID()) == Tree2BoolExpr::iso2boolExpr_.end() ||
         iso.getDrivers().front() != getNewIterationInputsETS().first[i])) {
        reachedPIs = false;
        break;
      }
    }
    DEBUG_LOG("--- End of iteration %zu\n", iter);
    iter++;
  }

  copyNewIterationInputsETStoCurrent();
  #ifdef DEBUG_CHECKS
  size_t finalSize = sizeOfCurrentIterationInputsETS();
  #endif
  copyCurrentIterationInputsETS(currentIterationInputs_);
  #ifdef DEBUG_CHECKS
  assert(finalSize == currentIterationInputs_.size() &&
         "compute: size mismatch after final copy");
  //assert(currentIterationInputs_.size() == sizeOfCurrentIterationInputsETS());
  for (const auto& input : currentIterationInputs_) {
    auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(input).getIsoID());
    assert(isInput(input) || (Tree2BoolExpr::iso2boolExpr_.find(
            dnl_.getDNLTerminalFromID(input)
                .getIsoID()) != Tree2BoolExpr::iso2boolExpr_.end() &&
            iso.getDrivers().front() == input));
  }
  #endif
}
