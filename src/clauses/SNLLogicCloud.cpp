#include "SNLLogicCloud.h"
#include "SNLDesignTruthTable.h"
#include "SNLTruthTableMerger.h"
// For std::assert
#include <cassert>

using namespace KEPLER_FORMAL;
using namespace naja::DNL;

bool SNLLogicCloud::isInput(naja::DNL::DNLID termID) {
  return std::find(PIs_.begin(), PIs_.end(), termID) != PIs_.end();
}

void SNLLogicCloud::compute() {
  std::vector<naja::DNL::DNLID> newIterationInputs;
  // Collect the inputs of the leaf that his output is the seedOutputTerm_
  // and push them into newIterationInputs
  if (dnl_.getDNLTerminalFromID(seedOutputTerm_).isTopPort()) {
    // Get equiposequential of the seed output term
    auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(seedOutputTerm_).getIsoID());
    if (iso.getDrivers().size() != 1) {
      throw std::runtime_error("Seed output term is not a single driver");
    }
    auto driver = iso.getDrivers().front();
    auto inst = dnl_.getDNLInstanceFromID(driver);
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      if (isInput(termID)) {
        newIterationInputs.push_back(termID);
      }
    }
    table_ = SNLDesignTruthTable::getTruthTable(inst.getSNLModel());
  } else {
    auto inst = dnl_.getDNLInstanceFromID(seedOutputTerm_);
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      if (isInput(termID)) {
        newIterationInputs.push_back(termID);
      }
    }
    table_ = SNLDesignTruthTable::getTruthTable(inst.getSNLModel());
  }
  // If no inputs, return
  if (newIterationInputs.empty()) {
    return;
  }

  // Expand the truth table until reaching a terminal which will be postive for
  // isInput
  bool reachedPIs = true;
  for (auto input : newIterationInputs) {
    if (!isInput(input)) {
      reachedPIs = false;
      break;
    }
  }
  while (!reachedPIs) {
    inputs_ = newIterationInputs;
    std::vector<const naja::NL::SNLTruthTable*> inputsToMerge;
    for (auto input : inputs_) {
      if (!isInput(input)) {
        // Create a feed trhough snl truth table with mask 10 and 1 input
        SNLTruthTable tt(2, 1);
        newIterationInputs.push_back(input);
        continue;
      }
      auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          dnl_.getDNLTerminalFromID(seedOutputTerm_).getIsoID());
      if (iso.getDrivers().size() != 1) {
        throw std::runtime_error("Seed output term is not a single driver");
      }
      auto driver = iso.getDrivers().front();
      auto inst = dnl_.getDNLInstanceFromID(driver);
      for (DNLID termID = inst.getTermIndexes().first;
           termID <= inst.getTermIndexes().second; termID++) {
        if (isInput(termID)) {
          newIterationInputs.push_back(termID);
        }
      }
    }
    if (inputsToMerge.empty()) {
      break;
    }
    // Merge the truth tables
    auto base = SNLDesignTruthTable::getTruthTable(
        dnl_.getDNLTerminalFromID(seedOutputTerm_)
            .getDNLInstance()
            .getSNLModel());
    if (!base.isInitialized()) {
      throw std::runtime_error("Base truth table is null");
    }
    SNLTruthTableMerger merger(inputsToMerge, table_);
    merger.computeMerged();
    table_ = merger.getMergedTable();
    for (auto input : newIterationInputs) {
      if (!isInput(input)) {
        reachedPIs = false;
        break;
      }
    }
  }
}