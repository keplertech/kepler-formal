// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "BuildPrimaryOutputClauses.h"
#include "DNL.h"
#include "SNLDesignModeling.h"
#include "SNLLogicCloud.h"
#include "Tree2BoolExpr.h"
#include "SNLPath.h"
#include <thread>
#include <tbb/global_control.h>

//#define DEBUG_PRINTS
//#define DEBUG_CHECKS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace KEPLER_FORMAL;
using namespace naja::DNL;
using namespace naja::NL;

namespace {

BuildPrimaryOutputClauses::PathNameIDs getPathNameIDs(const SNLPath& path) {
  BuildPrimaryOutputClauses::PathNameIDs ids;
  const auto pathNames = path.getPathNames();
  ids.reserve(pathNames.size());
  for (const auto& name : pathNames) {
    ids.push_back(name.getID());
  }
  return ids;
}

BuildPrimaryOutputClauses::PathKey getTerminalPathKey(const DNLTerminalFull& terminal) {
  auto pathIDs = getPathNameIDs(terminal.getDNLInstance().getPath());
  pathIDs.push_back(terminal.getSnlBitTerm()->getName().getID());
  BuildPrimaryOutputClauses::PathObjectIDs objectIDs = {
      static_cast<NLID::DesignObjectID>(terminal.getSnlBitTerm()->getBit())};
  return {std::move(pathIDs), std::move(objectIDs)};
}

}  // namespace

std::vector<DNLID> BuildPrimaryOutputClauses::collectInputs() {
  std::vector<DNLID> inputs;
  auto dnl = get();
  DNLInstanceFull top = dnl->getTop();

  for (DNLID termId = top.getTermIndexes().first;
       termId != DNLID_MAX && termId <= top.getTermIndexes().second; termId++) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
    if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
      DEBUG_LOG("Collecting input %s\n",
                term.getSnlBitTerm()->getName().getString().c_str());
      assert(termId < naja::DNL::get()->getDNLTerms().size());
      inputs.emplace_back(termId);
    }
  }

  for (DNLID leaf : dnl->getLeaves()) {
    auto iter = modelCache_.find(dnl->getDNLInstanceFromID(leaf).getSNLModel());
    const DNLInstanceFull& instance = dnl->getDNLInstanceFromID(leaf);
    if ((iter != modelCache_.end()) && iter->second.analyzedPIs) {
      const auto& cache = iter->second;
      for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
         const auto& term = dnl->getDNLTerminalFromID(termId);
          if (cache.PIs.find(term.getSnlBitTerm()) != cache.PIs.end()) {
            inputs.emplace_back(termId);
          }
      }
      continue;
    }
    modelCache_[instance.getSNLModel()].analyzedPIs = true;
    size_t numberOfInputs = 0, numberOfOutputs = 0;
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output)
        numberOfInputs++;
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Input)
        numberOfOutputs++;
    }

    if (numberOfInputs == 0 && numberOfOutputs > 0) {
      // no inputs primtive, assuming generator so adding as PI
      for (DNLID termId = instance.getTermIndexes().first;
           termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
           termId++) {
        const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Input) {
          assert(termId < naja::DNL::get()->getDNLTerms().size());
          inputs.emplace_back(termId);
          modelCache_[instance.getSNLModel()].PIs.insert(
              term.getSnlBitTerm());
          DEBUG_LOG(
              "Collecting input %s of model %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
        }
      }
      continue;
    }

    bool isSequential = false;
    std::vector<SNLBitTerm*> seqBitTerms;
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      auto related =
          SNLDesignModeling::getClockRelatedOutputs(term.getSnlBitTerm());
      if (!related.empty()) {
        isSequential = true;
        for (auto bitTerm : related) {
          seqBitTerms.emplace_back(bitTerm);
        }
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Input) {
          assert(termId < naja::DNL::get()->getDNLTerms().size());
          inputs.emplace_back(termId);
          modelCache_[instance.getSNLModel()].PIs.insert(
              term.getSnlBitTerm());
          DEBUG_LOG(
              "Collecting seq input %s of model %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
        }
      }
    }
    if (!isSequential) {
      for (DNLID termId = instance.getTermIndexes().first;
           termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
           termId++) {
        const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Input) {
          auto deps =
              SNLDesignModeling::getCombinatorialInputs(term.getSnlBitTerm());
          const auto tt = SNLDesignModeling::getTruthTable(term.getSnlBitTerm()->getDesign(), 
              term.getSnlBitTerm()->getOrderID());
          if (!tt.isInitialized()) {
            assert(termId < naja::DNL::get()->getDNLTerms().size());
            inputs.emplace_back(termId);
            modelCache_[instance.getSNLModel()].PIs.insert(
                term.getSnlBitTerm());
            DEBUG_LOG("Collecting input %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(),
                      term.getSnlBitTerm()
                          ->getDesign()
                          ->getName()
                          .getString()
                          .c_str());
          }
          
          if (tt.all0() ||
              tt.all1()) {
            assert(termId < naja::DNL::get()->getDNLTerms().size());
            inputs.emplace_back(termId);
            modelCache_[instance.getSNLModel()].PIs.insert(
                term.getSnlBitTerm());
            DEBUG_LOG("Collecting constant input %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(),
                      term.getSnlBitTerm()
                          ->getDesign()
                          ->getName()
                          .getString()
                          .c_str());
          }
        }
      }
      continue;
    }
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (term.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Input) {
        if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                      term.getSnlBitTerm()) != seqBitTerms.end()) {
          assert(termId < naja::DNL::get()->getDNLTerms().size());
          inputs.emplace_back(termId);
          modelCache_[instance.getSNLModel()].PIs.insert(
              term.getSnlBitTerm());
          DEBUG_LOG(
              "Collecting seq input %s of model %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
        }
      }
    }
  }
  std::set<DNLID> inputSet(inputs.begin(), inputs.end());
  inputs.clear();
  inputs.assign(inputSet.begin(), inputSet.end());
  DEBUG_LOG("Collected %zu inputs\n", inputs.size());
  return inputs;
}

std::vector<DNLID> BuildPrimaryOutputClauses::collectOutputs() {
  std::vector<DNLID> outputs;
  std::set<DNLID> outputsSet;
  auto dnl = get();
  DNLInstanceFull top = dnl->getTop();

  for (DNLID termId = top.getTermIndexes().first;
       termId != DNLID_MAX && termId <= top.getTermIndexes().second; termId++) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
    if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Input) {
      outputsSet.insert(termId);
      DEBUG_LOG(
          "Collecting top output %s of model %s\n",
          term.getSnlBitTerm()->getName().getString().c_str(),
          term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
    }
  }
  for (DNLID leaf : dnl->getLeaves()) {
    const DNLInstanceFull& instance = dnl->getDNLInstanceFromID(leaf);
    auto iter = modelCache_.find(instance.getSNLModel());
    if ((iter != modelCache_.end()) && iter->second.analyzedPOs) {
      const auto& cache = iter->second;
      for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
         const auto& term = dnl->getDNLTerminalFromID(termId);
          if (cache.POs.find(term.getSnlBitTerm()) != cache.POs.end()) {
            outputsSet.insert(termId);
          }
      }
      continue;
    }
    modelCache_[instance.getSNLModel()].analyzedPOs = true;
    bool isSequential = false;
    std::vector<SNLBitTerm*> seqBitTerms;

    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      auto related =
          SNLDesignModeling::getClockRelatedInputs(term.getSnlBitTerm());
      if (!related.empty()) {
        isSequential = true;
        for (auto bitTerm : related) {
          seqBitTerms.emplace_back(bitTerm);
        }
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Output) {
          outputsSet.insert(termId);
          modelCache_[instance.getSNLModel()].POs.insert(
              term.getSnlBitTerm());
          DEBUG_LOG(
              "Collecting seq output %s of model %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
        }
      }
    }

    if (!isSequential) {
      uint64_t inputNum = 0;
      for (DNLID termId = instance.getTermIndexes().first;
           termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
           termId++) {
        const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Output) {
          uint64_t orderID = inputNum;
          inputNum++;
          auto deps =
              SNLDesignModeling::getCombinatorialOutputs(term.getSnlBitTerm());
          // Collect all tt on the model
          std::vector<SNLTruthTable> tts;
          for (DNLID tId = instance.getTermIndexes().first;
               tId != DNLID_MAX && tId <= instance.getTermIndexes().second;
               tId++) {
            const DNLTerminalFull& tTerm = dnl->getDNLTerminalFromID(tId);
            // If direction is input, skip
            if (tTerm.getSnlBitTerm()->getDirection() ==
                SNLBitTerm::Direction::Input) {
              continue;
            }
            const auto tt = SNLDesignModeling::getTruthTable(tTerm.getSnlBitTerm()->getDesign(), 
              tTerm.getSnlBitTerm()->getOrderID());
            if (tt.isInitialized()) {
              tts.emplace_back(tt);
              // print deps
              for (const auto& d : tt.getDependencies()) {
                DEBUG_LOG("TT deps: %llu\n", d);
              }
            } else if (tt.all0() || tt.all1()) {
              tts.emplace_back(tt);
            }
          }
          bool inTermInTTDeps = false;
          for (const auto tt : tts) {
            const auto ttDeps =
                tt.getDependencies();  // expect std::vector<uint64_t>
            // WRONG!!! We need the input's "orderID" not general terminal orderID
            //uint64_t orderID =
            //    term.getSnlBitTerm()->getOrderID();  // assume 0-based

            // If orderID is 1-based, uncomment:
            // if (orderID > 0) --orderID;

            for (size_t index = 0; index < ttDeps.size(); ++index) {
              uint64_t d = ttDeps[index];

              uint64_t blockMin = index * 64ULL;     // inclusive
              uint64_t blockMax = blockMin + 64ULL;  // exclusive

              // If orderID is before this block, nothing more to find
              if (orderID < blockMin) {
                break;
              }

              // Skip if orderID beyond this block
              if (orderID >= blockMax) {
                continue;
              }

              uint64_t localBit = orderID - blockMin;  // 0..63

              // std::printf("TT deps: %llu\n", static_cast<unsigned long
              // long>(d)); std::printf("d = %llu, orderID = %llu, index = %zu,
              // localBit = %llu\n",
              //  static_cast<unsigned long long>(d),
              //  static_cast<unsigned long long>(orderID),
              //  index,
              //  static_cast<unsigned long long>(localBit));

              // Defensive: check localBit < 64
              if (localBit >= 64) {
                // std::fprintf(stderr, "localBit out of range: %llu\n",
                //             static_cast<unsigned long long>(localBit));
                // LCOV_EXCL_START
                continue;
                // LCOV_EXCL_STOP
              }

              // Correct shift using 1ULL and parentheses
              assert(localBit < 64);
              uint64_t mask = (1ULL << localBit);
              // std::printf("mask = 0x%llx\n", static_cast<unsigned long
              // long>(mask));

              if ((d & mask) != 0ULL) {
                inTermInTTDeps = true;
              }

              // we handled the block which contains orderID; stop scanning
              // ttDeps
              break;
            }

            if (inTermInTTDeps) {
              // Found — act and break outer loop if that's desired
              // std::puts("Found matching bit in this TT deps");
              break;
            }
          }
          if (/*deps.empty() &&*/ !inTermInTTDeps) {
            outputsSet.insert(termId);
            modelCache_[instance.getSNLModel()].POs.insert(
              term.getSnlBitTerm());
            DEBUG_LOG("Collecting output %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(),
                      term.getSnlBitTerm()
                          ->getDesign()
                          ->getName()
                          .getString()
                          .c_str());
          }
        }
      }
      continue;
    }
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (term.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Output) {
        if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                      term.getSnlBitTerm()) != seqBitTerms.end())
          outputsSet.insert(termId);
          modelCache_[instance.getSNLModel()].POs.insert(
              term.getSnlBitTerm());
        DEBUG_LOG(
            "Collecting seq output %s of model %s\n",
            term.getSnlBitTerm()->getName().getString().c_str(),
            term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
      }
    }
  }
  outputs.clear();
  // keep only terminals who are connected to nets
  for (const auto& out : outputsSet) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(out);
    if (term.getIsoID() == DNLID_MAX) {
      DEBUG_LOG("Skipping output %s of model %s as it is not connected to any net\n",
                term.getSnlBitTerm()->getName().getString().c_str(),
                term.getSnlBitTerm()
                    ->getDesign()
                    ->getName()
                    .getString()
                    .c_str());
      continue;
    }
    if (term.getIsoID() != DNLID_MAX && 
      dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID()).getDrivers().empty()) {
      DEBUG_LOG("Skipping output %s of model %s as it is not connected to any net\n",
                term.getSnlBitTerm()->getName().getString().c_str(),
                term.getSnlBitTerm()
                    ->getDesign()
                    ->getName()
                    .getString()
                    .c_str());
      continue;
    }
    if (term.getIsoID() != DNLID_MAX && 
      dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID()).getDrivers().size() > 1) {
      DEBUG_LOG("Skipping output %s of model %s as it is driven by multiple drivers\n",
                term.getSnlBitTerm()->getName().getString().c_str(),
                term.getSnlBitTerm()
                    ->getDesign()
                    ->getName()
                    .getString()
                    .c_str());
      continue;
    }
    outputs.emplace_back(out);
  }
  //outputs.assign(outputsSet.begin(), outputsSet.end());
  return outputs;
}

void BuildPrimaryOutputClauses::collect() {
  inputs_ = collectInputs();
  //sortInputs(); <- cannot sort inputs as it has to respect the inputs vector order
  for (const auto& input : inputs_) {
    PathKey key = getTerminalPathKey(naja::DNL::get()->getDNLTerminalFromID(input));
    inputsMap_[std::move(key)]  =
            input;
  }
  outputs_ = collectOutputs();
  //sortOutputs(); <- cannot sort as it needs to keep the order for POs_
  for (const auto& output : outputs_) {
    PathKey key = getTerminalPathKey(naja::DNL::get()->getDNLTerminalFromID(output));
    outputsMap_[std::move(key)]  =
            output;
    DEBUG_LOG("Output collected: %s\n", naja::DNL::get()
                                         ->getDNLTerminalFromID(output)
                                         .getSnlBitTerm()
                                         ->getName()
                                         .getString()
                                         .c_str());
  }
  POs_.resize(outputs_.size());
}

void BuildPrimaryOutputClauses::initVarNames() {
  termDNLID2varID_.resize(naja::DNL::get()->getDNLTerms().size(), (size_t)-1);
  for (size_t i = 0; i < inputs_.size(); ++i) {
    // Get Truth Table for terminal
    const DNLTerminalFull& tTerm = naja::DNL::get()->getDNLTerminalFromID(inputs_[i]);
    // If direction is input, skip
    if (!tTerm.isTopPort()) {
      const auto tt = SNLDesignModeling::getTruthTable(tTerm.getSnlBitTerm()->getDesign(), 
      tTerm.getSnlBitTerm()->getOrderID());
      if (tt.isInitialized()) {
        if (tt.all0()) {
          termDNLID2varID_[inputs_[i]] = 0;
          continue;
        } else if (tt.all1()) {
          termDNLID2varID_[inputs_[i]] = 1;
          continue;
        }
      }
    }
    termDNLID2varID_[inputs_[i]] =
        i + 2;  // +2 to avoid 0 and 1 which are reserved for constants
  }
  for (DNLID constIsoID : naja::DNL::get()->getDNLIsoDB().getConstant0Isos()) {
    const auto& constIso = naja::DNL::get()->getDNLIsoDB().getIsoFromIsoIDconst(constIsoID);
    for (auto termID : constIso.getReaders()) {
        termDNLID2varID_[termID] = 0;
    }
    for (auto termID : constIso.getDrivers()) {
        termDNLID2varID_[termID] = 0;
    }
  }
  for (DNLID constIsoID : naja::DNL::get()->getDNLIsoDB().getConstant1Isos()) {
    const auto& constIso = naja::DNL::get()->getDNLIsoDB().getIsoFromIsoIDconst(constIsoID);
    for (auto termID : constIso.getReaders()) {
        termDNLID2varID_[termID] = 1;
    }
    for (auto termID : constIso.getDrivers()) {
        termDNLID2varID_[termID] = 1;
    }
  }
}

void BuildPrimaryOutputClauses::build() {
  //printf("Building primary output clauses\n");
  naja::DNL::get();
  POs_.clear();
  POs_ = tbb::concurrent_vector<BoolExpr*>(outputs_.size());
  initVarNames();
  // Init var names(counting on the fact that normalization happened before)

  // inputs_ = collectInputs();
  // sortInputs();
  // outputs_ = collectOutputs();
  // sortOutputs();
  size_t processedOutputs = 0;
  // tbb::task_arena arena(20);
  //  init arena with automatic number of threads
  // unsigned hw = std::thread::hardware_concurrency(); 
  // if (hw == 0) hw = 1; // fallback 
  tbb::task_arena arena(20);
  IsPIs_ = std::vector<bool>(naja::DNL::get()->getNBterms(), false);
  for (auto pi : inputs_) {
    if (pi >= IsPIs_.size()) {
      // LCOV_EXCL_START
      std::string error = "PI " + std::to_string(pi) + " is out of range";
      throw std::runtime_error(error);
      // LCOV_EXCL_STOP
    }
    IsPIs_[pi] = true;
  }
  IsPOs_ = std::vector<bool>(naja::DNL::get()->getNBterms(), false);
  for (auto po : outputs_) {
    if (po >= IsPOs_.size()) {
      // LCOV_EXCL_START
      std::string error = "PO " + std::to_string(po) + " is out of range";
      throw std::runtime_error(error);
      // LCOV_EXCL_STOP
    }
    IsPOs_[po] = true;
  }
  auto processOutput = [&](size_t i) {
    DNLID out = outputs_[i];
    #ifdef DEBUG_PRINTS
    printf("Procssing output %zu/%zu: %s\n", ++processedOutputs,
           outputs_.size(),
           get()
               ->getDNLTerminalFromID(out)
               .getSnlBitTerm()
               ->getName()
               .getString()
               .c_str());
    #endif

    DNLID isoID = get()->getDNLTerminalFromID(out).getIsoID();
    DEBUG_LOG("isoID: %zu\n", isoID);
    if (Tree2BoolExpr::iso2boolExpr_.find(isoID) != Tree2BoolExpr::iso2boolExpr_.end()) {
      POs_[i] = Tree2BoolExpr::iso2boolExpr_[isoID];
      #ifdef DEBUG_CHECKS
      assert(POs_[i] != nullptr);
      #endif
      #ifdef DEBUG_PRINTS
      printf("Reusing iso output %s for output %s\n",
             POs_[i]->toString().c_str(),
             get()
                 ->getDNLTerminalFromID(out)
                 .getSnlBitTerm()
                 ->getName()
                 .getString()
                 .c_str());
      #endif
      return;
    }
    
    SNLLogicCloud cloud(out, IsPIs_, IsPOs_);
    #ifdef DEBUG_CHECKS
    auto startComp = std::chrono::steady_clock::now();
    #endif
    cloud.compute();
    #ifdef DEBUG_CHECKS
    auto endComp = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_comp = endComp - startComp;
    printf("Computation time for %lu: %f seconds\n", i, elapsed_seconds_comp.count());
    #endif
    // //cloud.SNLDesignModeling::getTruthTable().print();
    // std::vector<DNLID> test1;
    // std::vector<DNLID> test2;
    // for (auto in : cloud.getAllInputs()) {
    //   printf("Input in tree cloud: %lu\n", in);
    //   // if (in >= cloud.getInputs().size()) {
    //   //   printf("size of inputs in cloud: %lu\n",
    //   cloud.getInputs().size());
    //   //   //assert(false && "Input in cloud is out of range");
    //   // }
    //  test1.emplace_back(in);
    // }
    // for (auto in : cloud.getInputs()) {
    //   printf("Input in cloud: %lu\n", in);
    //   test2.emplace_back(in);
    // }
    // std::sort(test1.begin(), test1.end());
    // std::sort(test2.begin(), test2.end());
    // assert(test1 == test2);
    // std::vector<std::string> varNames;
    /*for (auto input : cloud.getInputs()) {
      DNLTerminalFull term = get()->getDNLTerminalFromID(input);
      if (term.getSnlTerm() != nullptr) {
        auto net = term.getSnlTerm()->getNet();
        if (net != nullptr) {
          if (net->isConstant0()) {
            varNames.emplace_back("0");
            continue;
          } else if (net->isConstant1()) {
            varNames.emplace_back("1");
            continue;
          }
        }
        auto model = const_cast<SNLDesign*>(
            term.getSnlBitTerm()->getDesign());
        auto tt = model->SNLDesignModeling::getTruthTable(term.getSnlBitTerm()->getOrderID());
        if (tt.isInitialized()) {
          if (tt.all0()) {
            varNames.emplace_back("0");
            continue;
          } else if (tt.all1()) {
            varNames.emplace_back("1");
            continue;
          }
        }
      }
      // find the index of input in inputs_
      auto it = std::find(inputs_.begin(), inputs_.end(), input);
      // printf("Input: %s\n",
      //
    get()->getDNLTerminalFromID(input).getSnlBitTerm()->getName().getString().c_str());
      // printf("Model: %s\n",
      //
    get()->getDNLTerminalFromID(input).getSnlBitTerm()->getDesign()->getName().getString().c_str());
      assert(it != inputs_.end());
      size_t index = std::distance(inputs_.begin(), it);
      varNames.emplace_back(std::to_string(index + 2)); // +2 to avoid 0 and 1
    which are reserved for constants
    }*/
#ifdef DEBUG_CHECKS
    assert(cloud.SNLDesignModeling::getTruthTable().isInitialized());
#endif
    // DEBUG_LOG("Truth Table: %s\n",
    //           cloud.SNLDesignModeling::getTruthTable().print().c_str());
    /*std::shared_ptr<BoolExpr> expr = Tree2BoolExpr::convert(
        cloud.SNLDesignModeling::getTruthTable(), varNames);*/
    // BoolExpr::getMutex().lock();
    //  if (POs_.size() - 1 < i) {
    //    for (size_t j = POs_.size(); j <= i; ++j) {
    //      POs_.emplace_back(nullptr);
    //    }
    //  }
    assert(POs_.size() - 1 >= i);
    // add run time counter here
    #ifdef DEBUG_CHECKS
    auto startFin = std::chrono::steady_clock::now();    
    #endif
    cloud.getTruthTable().finalize();
    #ifdef DEBUG_CHECKS
    auto endFin = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_fin = endFin - startFin;
    printf("Finalization time for %lu: %f seconds\n", i, elapsed_seconds_fin.count());
    #endif
    #ifdef DEBUG_CHECKS
    auto startConv = std::chrono::steady_clock::now();
    #endif
    POs_[i] = Tree2BoolExpr::convert(cloud.getTruthTable(), termDNLID2varID_);
    #ifdef DEBUG_CHECKS
    auto endConv = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_conv = endConv - startConv;
    printf("Conversion time for %lu: %f seconds\n", i, elapsed_seconds_conv.count());
    #endif
    cloud.destroy();
    // BoolExpr::getMutex().unlock();
    // printf("size of expr: %lu\n", POs_.back()->size());
    Tree2BoolExpr::iso2boolExpr_[isoID] = POs_[i];
  };
  Tree2BoolExpr::iso2boolExpr_.clear();
  if (getenv("KEPLER_NO_MT")) {
    for (size_t i = 0; i < outputs_.size(); ++i) {
      processOutput(i);
    }
  } else {
    // compute grain safely
    size_t n = outputs_.size();
    size_t default_grain = 1000;
    size_t computed = (n >= 1000) ? (n / 1000) : 1; // never zero
    size_t grain = std::max<size_t>(computed, default_grain); // or clamp as you prefer

    tbb::parallel_for(
      tbb::blocked_range<DNLID>(0, n, grain),
      [&](const tbb::blocked_range<DNLID>& r) {
        for (DNLID i = r.begin(); i < r.end(); ++i) {
          processOutput(i);
        }
      },
      tbb::static_partitioner()
    );
  }
  destroy();  // Clean up DNL instance
}

void BuildPrimaryOutputClauses::setInputs2InputsIDs() {
  //printf("Setting inputs to input IDs mapping\n");
  inputs2inputsIDs_.clear();
  for (const auto& input : inputs_) {
    if (get()->getDNLTerminalFromID(input).isNull()) {
      throw std::runtime_error("Input terminal is null");
    }
    const DNLInstanceFull& currentInstance =
        get()->getDNLTerminalFromID(input).getDNLInstance();
   
    //std::vector<NLID::DesignObjectID> termIDs;
    //termIDs.emplace_back(
    //     get()->getDNLTerminalFromID(input).getSnlBitTerm()->getID());
    //termIDs.emplace_back(
    //    get()->getDNLTerminalFromID(input).getSnlBitTerm()->getBit());
    PathKey& pair = inputs2inputsIDs_[input];
    pair.first = getPathNameIDs(currentInstance.getPath());
    pair.first.emplace_back(
        get()->getDNLTerminalFromID(input).getSnlBitTerm()->getName().getID());
    pair.second.emplace_back(
        get()->getDNLTerminalFromID(input).getSnlBitTerm()->getBit());
  }
}

void BuildPrimaryOutputClauses::setOutputs2OutputsIDs() {
  //printf("Setting outputs to output IDs mapping\n");
  outputs2outputsIDs_.clear();
  for (const auto& output : outputs_) {
    //std::vector<NLID::DesignObjectID> termIDs;
    const DNLInstanceFull& currentInstance =
        get()->getDNLTerminalFromID(output).getDNLInstance();
    //termIDs.emplace_back(
    //     get()->getDNLTerminalFromID(output).getSnlBitTerm()->getID());
    //termIDs
    PathKey& pair = outputs2outputsIDs_[output];
    pair.first = getPathNameIDs(currentInstance.getPath());
    pair.first.emplace_back(
        get()->getDNLTerminalFromID(output).getSnlBitTerm()->getName().getID());
    pair.second.emplace_back(
        get()->getDNLTerminalFromID(output).getSnlBitTerm()->getBit());
  }
}

// Sort functions are retierd for now as they break the mapping between the 2 circuits, normalize is used instead

// void BuildPrimaryOutputClauses::sortInputs() {
//   // Sort based on inputs2inputsIDs_ content
//   std::sort(inputs_.begin(), inputs_.end(),
//             [this](const DNLID& a, const DNLID& b) {
//               return inputs2inputsIDs_[a].first < inputs2inputsIDs_[b].first && 
//                       inputs2inputsIDs_[a].second < inputs2inputsIDs_[b].second;
//             });
// }

// void BuildPrimaryOutputClauses::sortOutputs() {
//   // Sort based on outputs2outputsIDs_ content
//   std::sort(
//       outputs_.begin(), outputs_.end(), [this](const DNLID& a, const DNLID& b) {
//         return outputs2outputsIDs_[a].first < outputs2outputsIDs_[b].first && 
//                outputs2outputsIDs_[a].second < outputs2outputsIDs_[b].second;
//       });
// }

// const naja::NL::SNLTruthTable& BuildPrimaryOutputClauses::getTruthTable(naja::NL::SNLDesign* design, size_t orderID) {
//   auto designID = design->getID();
//   auto iter = ttCache_.find({designID, orderID});
//   if (iter != ttCache_.end()) {
//     return iter->second;
//   }
//   const auto tt = SNLDesignModeling::getTruthTable(design, orderID);
//   ttCache_[{designID, orderID}] = tt;
//   return tt;
// }
