#include "SNLLogicCloud.h"
#include "SNLDesignTruthTable.h"
#include "SNLTruthTableMerger.h"
#include <cassert>

//#define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace KEPLER_FORMAL;
using namespace naja::DNL;

bool SNLLogicCloud::isInput(naja::DNL::DNLID termID) {
  return std::find(PIs_.begin(), PIs_.end(), termID) != PIs_.end();
}

bool SNLLogicCloud::isOutput(naja::DNL::DNLID termID) {
  return std::find(POs_.begin(), POs_.end(), termID) != POs_.end();
}

void SNLLogicCloud::compute() {
  std::vector<naja::DNL::DNLID> newIterationInputs;

  if (dnl_.getDNLTerminalFromID(seedOutputTerm_).isTopPort() || isOutput(seedOutputTerm_)) {
    auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(seedOutputTerm_).getIsoID());
    if (iso.getDrivers().size() != 1) {
      throw std::runtime_error("Seed output term is not a single driver");
    }
    auto driver = iso.getDrivers().front();
    auto inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
    if (isInput(driver)) {
      currentIterationInputs_.push_back(driver);
      SNLTruthTable tt(1, 2);
      table_ = tt;
      return;
    }
    DEBUG_LOG("Instance name: %s\n",
              inst.getSNLInstance()->getName().getString().c_str());
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
        newIterationInputs.push_back(termID);
      }
    }
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    table_ = SNLDesignTruthTable::getTruthTable(inst.getSNLModel());
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  } else {
    auto inst = dnl_.getDNLInstanceFromID(seedOutputTerm_);
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
        newIterationInputs.push_back(termID);
      }
    }
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    table_ = SNLDesignTruthTable::getTruthTable(inst.getSNLModel());
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  }

  if (newIterationInputs.empty()) {
    DEBUG_LOG("No inputs found for seed output term %zu\n", seedOutputTerm_);
    return;
  }

  bool reachedPIs = true;
  for (auto input : newIterationInputs) {
    if (!isInput(input) && !isOutput(input)) {
      reachedPIs = false;
      break;
    }
  }

  for (auto input : newIterationInputs) {
    DEBUG_LOG("newIterationInputs Input: %s(%s)\n",
              dnl_.getDNLTerminalFromID(input).getSnlBitTerm()->getName().getString().c_str(),
              dnl_.getDNLTerminalFromID(input).getSnlBitTerm()->getDesign()->getName().getString().c_str());
  }

  while (!reachedPIs) {
    DEBUG_LOG("---iter---\n");
    DEBUG_LOG("Current iteration inputs: %lu\n", newIterationInputs.size());
    currentIterationInputs_ = newIterationInputs;
    newIterationInputs.clear();
    DEBUG_LOG("Truth table: %s\n", table_.getString().c_str());
    assert(currentIterationInputs_.size() == table_.size());

    std::vector<naja::NL::SNLTruthTable> inputsToMerge;
    for (auto input : currentIterationInputs_) {
      if (isInput(input) || isOutput(input)) {
        SNLTruthTable tt(1, 2);
        newIterationInputs.push_back(input);
        DEBUG_LOG("Adding input: %s\n",
                  dnl_.getDNLTerminalFromID(input).getSnlBitTerm()->getName().getString().c_str());
        inputsToMerge.push_back(tt);
        continue;
      }

      auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          dnl_.getDNLTerminalFromID(input).getIsoID());
      DEBUG_LOG("number of drivers: %zu\n", iso.getDrivers().size());

      for (auto driver : iso.getDrivers()) {
        DEBUG_LOG("Driver: %s\n",
                  dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getName().getString().c_str());
      }

      if (iso.getDrivers().size() >= 1) {
        assert(iso.getDrivers().size() <= 1 &&
               "Iso have more than one driver, not supported");
      } else if (iso.getDrivers().empty()) {
        assert(iso.getDrivers().size() == 1 &&
               "Iso have no drivers and more than one reader, not supported");
      }

      auto driver = iso.getDrivers().front();
      if (isInput(driver) || isOutput(driver)) {
        SNLTruthTable tt(1, 2);
        newIterationInputs.push_back(driver);
        DEBUG_LOG("Adding top input: %s\n",
                  dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getName().getString().c_str());
        inputsToMerge.push_back(tt);
        continue;
      }

      auto inst = dnl_.getDNLInstanceFromID(
          dnl_.getDNLTerminalFromID(driver).getDNLInstance().getID());
      assert(SNLDesignTruthTable::getTruthTable(inst.getSNLModel()).isInitialized() &&
             "Truth table for instance is not initialized");

      inputsToMerge.push_back(SNLDesignTruthTable::getTruthTable(inst.getSNLModel()));

      for (DNLID termID = inst.getTermIndexes().first;
           termID <= inst.getTermIndexes().second; termID++) {
        const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
        if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
          DEBUG_LOG("Adding input: %s(%s)\n",
                    term.getSnlBitTerm()->getName().getString().c_str(),
                    term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
          newIterationInputs.push_back(termID);
        }
      }
    }

    if (inputsToMerge.empty()) {
      break;
    }

    DEBUG_LOG("Merging truth tables with %zu inputs\n", inputsToMerge.size());
    DEBUG_LOG("Truth table %s\n", table_.getString().c_str());

    SNLTruthTableMerger merger(inputsToMerge, table_);
    merger.computeMerged();
    table_ = merger.getMergedTable();

    reachedPIs = true;
    for (auto input : newIterationInputs) {
      if (!isInput(input) && !isOutput(input)) {
        reachedPIs = false;
        break;
      }
    }
  }

  currentIterationInputs_ = newIterationInputs;
  for (auto input : currentIterationInputs_) {
    DEBUG_LOG("Input: %s\n",
              dnl_.getDNLTerminalFromID(input).getSnlBitTerm()->getName().getString().c_str());
  }
}
