#include "SNLLogicCloud.h"
#include "SNLDesignModeling.h"
#include "SNLTruthTableMerger.h"
#include <cassert>
#include <tbb/tbb_allocator.h>
#include "tbb/enumerable_thread_specific.h"
#include "tbb/concurrent_vector.h"

typedef std::pair<std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>, size_t> IterationInputsETSPair;
tbb::enumerable_thread_specific< IterationInputsETSPair> currentIterationInputsETS;
tbb::enumerable_thread_specific< IterationInputsETSPair> newIterationInputsETS;

tbb::concurrent_vector<IterationInputsETSPair*> currentIterationInputsETSvector = tbb::concurrent_vector<IterationInputsETSPair*>(40, nullptr);

tbb::concurrent_vector<IterationInputsETSPair*> newIterationInputsETSvector = tbb::concurrent_vector<IterationInputsETSPair*>(40, nullptr);

void initCurrentIterationInputsETS() {
  // if (currentIterationInputsETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
  //   for (size_t i = currentIterationInputsETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
  //     currentIterationInputsETSvector.push_back(nullptr);
  //   }
  // }
  assert(tbb::this_task_arena::current_thread_index() < currentIterationInputsETSvector.size());
  if (currentIterationInputsETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
    currentIterationInputsETSvector[tbb::this_task_arena::current_thread_index()] = &currentIterationInputsETS.local();
  }
}

IterationInputsETSPair& getCurrentIterationInputsETS() {
  //initCurrentIterationInputsETS();
  //assert(tbb::this_task_arena::current_thread_index() < currentIterationInputsETSvector.size());
  return *currentIterationInputsETSvector[tbb::this_task_arena::current_thread_index()];
}

void initNewIterationInputsETS() {
  // if (newIterationInputsETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
  //   for (size_t i = newIterationInputsETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
  //     newIterationInputsETSvector.push_back(nullptr);
  //   }
  // }
  assert(tbb::this_task_arena::current_thread_index() < newIterationInputsETSvector.size());
  if (newIterationInputsETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
    newIterationInputsETSvector[tbb::this_task_arena::current_thread_index()] = &newIterationInputsETS.local();
  }
}

IterationInputsETSPair& getNewIterationInputsETS() {
  //initNewIterationInputsETS();
  //assert(tbb::this_task_arena::current_thread_index() < newIterationInputsETSvector.size());
  return *newIterationInputsETSvector[tbb::this_task_arena::current_thread_index()];
}

void clearCurrentIterationInputsETS() {
  getCurrentIterationInputsETS().second = 0;
}

void pushBackCurrentIterationInputsETS(naja::DNL::DNLID input) {
  auto & vec = getCurrentIterationInputsETS().first;
  auto & sz = getCurrentIterationInputsETS().second;
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
  auto & current = getCurrentIterationInputsETS();
  for (size_t i = 0; i < current.second; i++) {
    res.push_back(current.first[i]);
  }
}

void clearNewIterationInputsETS() {
  getNewIterationInputsETS().second = 0;
}

void pushBackNewIterationInputsETS(naja::DNL::DNLID input) {
  auto & vec = getNewIterationInputsETS().first;
  auto & sz = getNewIterationInputsETS().second;
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

//#define DEBUG_PRINTS

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
  //std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>> newIterationInputs;
  initNewIterationInputsETS();
  initCurrentIterationInputsETS();
  clearNewIterationInputsETS();
  clearCurrentIterationInputsETS();
  if (dnl_.getDNLTerminalFromID(seedOutputTerm_).isTopPort() || isOutput(seedOutputTerm_)) {
    auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(seedOutputTerm_).getIsoID());
    if (iso.getDrivers().size() != 1) {
      throw std::runtime_error("Seed output term is not a single driver");
    }
    auto driver = iso.getDrivers().front();
    auto inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
    if (isInput(driver)) {
      //currentIterationInputs_.push_back(driver);
      pushBackCurrentIterationInputsETS(driver);
      table_ = SNLTruthTableTree(inst.getID(), driver, SNLTruthTableTree::Node::Type::P);
      //DEBUG_LOG("Driver %s is a primary input, returning\n",
      //       dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getName().getString().c_str());
      return;
    }
    DEBUG_LOG("Instance name: %s\n",
              inst.getSNLInstance()->getName().getString().c_str());
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
        //newIterationInputs.push_back(termID);
        pushBackNewIterationInputsETS(termID);
        DEBUG_LOG("Add input with id: %zu\n", termID);
      }
    }
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    table_ = SNLTruthTableTree(inst.getID(), driver);
    auto* model = const_cast<SNLDesign*>(inst.getSNLModel());
    assert(model->getTruthTable(
        dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID()).isInitialized() &&
           "Truth table is not initialized");
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  } else {
    auto inst = dnl_.getDNLInstanceFromID(seedOutputTerm_);
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
        //newIterationInputs.push_back(termID);
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
  // for (auto input : newIterationInputs) {
  //   if (!isInput(input)/* && !isOutput(input)*/) {
  //     reachedPIs = false;
  //     break;
  //   }
  // }
  size_t size = sizeOfNewIterationInputsETS();
  for (size_t i = 0; i < size; i++) {
    if (!isInput(getNewIterationInputsETS().first[i])/* && !isOutput(getNewIterationInputsETS().first[i])*/) {
      reachedPIs = false;
      break;
    }
  }

  

  while (!reachedPIs) {
    //DEBUG_LOG("size of truth table tree: %zu\n", table_.size());
    DEBUG_LOG("---iter---\n");
    //DEBUG_LOG("Current iteration inputs: %lu\n", newIterationInputs.size());
    DEBUG_LOG("Current iteration inputs size: %zu\n", sizeOfNewIterationInputsETS());
    //DEBUG_LOG("term %lu: newIterationInputs size: %zu\n", seedOutputTerm_, newIterationInputs.size());
    // for (auto input : newIterationInputs) {
    //   DEBUG_LOG("newIterationInputs Input: %s(%s)\n",
    //             dnl_.getDNLTerminalFromID(input).getSnlBitTerm()->getName().getString().c_str(),
    //             dnl_.getDNLTerminalFromID(input).getSnlBitTerm()->getDesign()->getName().getString().c_str());
    // }
    //currentIterationInputs_ = newIterationInputs;
    copyNewIterationInputsETStoCurrent();
    // for (auto input : currentIterationInputs_) {
    //   DEBUG_LOG("Input: %s\n",
    //             dnl_.getDNLTerminalFromID(input).getSnlBitTerm()->getName().getString().c_str());
    // }
  #ifdef DEBUG_PRINTS
    for (size_t i = 0; i < sizeOfNewIterationInputsETS(); i++) {
      DEBUG_LOG("Input: %s\n",
                dnl_.getDNLTerminalFromID(currentIterationInputsETS.local().first[i]).getSnlBitTerm()->getName().getString().c_str());
    }
  #endif
    //newIterationInputs.clear();
    clearNewIterationInputsETS();
    //DEBUG_LOG("Truth table: %s\n", table_.getString().c_str());
    //DEBUG_LOG("Truth table size: %zu\n", table_.size());
    //DEBUG_LOG("Current iteration inputs size: %zu\n", currentIterationInputs_.size());
    //DEBUG_LOG("table size: %zu, currentIterationInputs_ size: %zu\n", table_.size(), currentIterationInputs_.size());
    DEBUG_LOG("table size: %zu, currentIterationInputs_ size: %zu\n", table_.size(), sizeOfCurrentIterationInputsETS());
    //assert(currentIterationInputs_.size() == table_.size());
    assert(sizeOfCurrentIterationInputsETS() == table_.size());

    std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
            tbb::tbb_allocator<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>>> inputsToMerge;

    //for (auto input : currentIterationInputs_) {
    size_t sizeOfCurrentInputs = sizeOfCurrentIterationInputsETS();
    for (size_t i = 0; i < sizeOfCurrentInputs; i++) {
      auto input = getCurrentIterationInputsETS().first[i];
      if (isInput(input)/*|| isOutput(input)*/) {
        //SNLTruthTable tt(1, 2); // uncommented
        //newIterationInputs.push_back(input);
        pushBackNewIterationInputsETS(input);
        DEBUG_LOG("Add input with id: %zu\n", input);
        DEBUG_LOG("Adding input: %s\n",
                  dnl_.getDNLTerminalFromID(input).getSnlBitTerm()->getName().getString().c_str());
        inputsToMerge.push_back({naja::DNL::DNLID_MAX, input}); // Placeholder for PI/PO
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
      if (isInput(driver)/* || isOutput(driver)*/) {
        //SNLTruthTable tt(1, 2);
        //newIterationInputs.push_back(driver);
        pushBackNewIterationInputsETS(driver);
        DEBUG_LOG("Add input with id: %zu\n", driver);
        DEBUG_LOG("Adding top input: %s\n",
                  dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getName().getString().c_str());
        inputsToMerge.push_back({naja::DNL::DNLID_MAX, driver}); // Placeholder for PI/PO
        continue;
      }

      auto inst = dnl_.getDNLInstanceFromID(
          dnl_.getDNLTerminalFromID(driver).getDNLInstance().getID());
      auto* model = const_cast<SNLDesign*>(inst.getSNLModel());
      if (!model->getTruthTable(
        dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID()).isInitialized())
      {
        //DEBUG_LOG("#####Truth table for instance %s is not initialized\n",
        //          inst.getSNLModel()->getName().getString().c_str());
        DEBUG_LOG("#####Truth table for instance %s is not initialized\n",
                  inst.getSNLInstance()->getModel()->getName().getString().c_str());
        auto* model = const_cast<SNLDesign*>(inst.getSNLModel());
        assert(model->getTruthTable(
        dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID()).isInitialized() &&
             "Truth table for instance is not initialized");
      }
      
      DEBUG_LOG("Instance name: %s\n",
                inst.getSNLInstance()->getName().getString().c_str());
      inputsToMerge.push_back({inst.getID(), driver});

      for (DNLID termID = inst.getTermIndexes().first;
           termID <= inst.getTermIndexes().second; termID++) {
        const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
        if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
          DEBUG_LOG("Add input with id: %zu\n", termID);
          DEBUG_LOG("Adding input: %s(%s)\n",
                    term.getSnlBitTerm()->getName().getString().c_str(),
                    term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
          //newIterationInputs.push_back(termID);
          pushBackNewIterationInputsETS(termID);
          //inputsToMerge.push_back({term.getDNLInstance().getID(), termID});
        }
      }
    }

    if (inputsToMerge.empty()) {
      break;
    }

    DEBUG_LOG("--- Merging truth tables with %zu inputs\n", inputsToMerge.size());
    //DEBUG_LOG("Truth table %s\n", table_.getString().c_str());
    table_.concatFull(inputsToMerge);
    reachedPIs = true;
    // for (auto input : newIterationInputs) {
    //   if (!isInput(input)/*&& !isOutput(input)*/) {
    //     reachedPIs = false;
    //     break;
    //   }
    // }
    size_t sizeOfNewInputs = sizeOfNewIterationInputsETS();
    for (size_t i = 0; i < sizeOfNewInputs; i++) {
      if (!isInput(getNewIterationInputsETS().first[i])) {
        reachedPIs = false;
        break;
      }
    }
  }
  copyNewIterationInputsETStoCurrent();
  copyCurrentIterationInputsETS(currentIterationInputs_);
  assert(currentIterationInputs_.size() == sizeOfCurrentIterationInputsETS());
  // Assert all currentIterationInputs_ are PIs
  for (auto input : currentIterationInputs_) {
    assert(isInput(input));
  }
  if (getAllInputs().size() != currentIterationInputs_.size()) {
    DEBUG_LOG("Number of inputs in the truth table: %zu, number of current iteration inputs: %zu\n",
           getAllInputs().size(), currentIterationInputs_.size());
    assert(false && "Number of inputs in the truth table does not match the number of current iteration inputs");
  }
}
