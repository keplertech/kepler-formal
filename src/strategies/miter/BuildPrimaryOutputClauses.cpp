#include "BuildPrimaryOutputClauses.h"
#include "DNL.h"
#include "SNLDesignModeling.h"
#include "SNLLogicCloud.h"
#include "SNLTruthTable2BoolExpr.h"
#include "Tree2BoolExpr.h"

//#define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...)  printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace KEPLER_FORMAL;
using namespace naja::DNL;
using namespace naja::NL;

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
      inputs.push_back(termId);
    }
  }

  for (DNLID leaf : dnl->getLeaves()) {
    DNLInstanceFull instance = dnl->getDNLInstanceFromID(leaf);
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

    if (numberOfInputs == 0 && numberOfOutputs > 1) {
      for (DNLID termId = instance.getTermIndexes().first;
           termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
           termId++) {
        const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Input)
          inputs.push_back(termId);
          DEBUG_LOG("Collecting input %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(), 
                      term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
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
          seqBitTerms.push_back(bitTerm);
        }
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Input) {
          inputs.push_back(termId);
          DEBUG_LOG("Collecting seq input %s of model %s\n",
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
          auto deps = SNLDesignModeling::getCombinatorialInputs(
              term.getSnlBitTerm());
          if (/*deps.empty() ||*/ !(term.getSnlBitTerm()->getDesign()->getTruthTable(term.getSnlBitTerm()->getOrderID()).isInitialized())) {
            inputs.push_back(termId);
            DEBUG_LOG("Collecting input %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(), 
                      term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
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
                      term.getSnlBitTerm()) != seqBitTerms.end())
          inputs.push_back(termId);
          DEBUG_LOG("Collecting seq input %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(), 
                      term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
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
      DEBUG_LOG("Collecting top output %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(), 
                      term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
    }
  }
  for (DNLID leaf : dnl->getLeaves()) {
    DNLInstanceFull instance = dnl->getDNLInstanceFromID(leaf);
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
          seqBitTerms.push_back(bitTerm);
        }
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Output) {
          outputsSet.insert(termId);
          DEBUG_LOG("Collecting seq output %s of model %s\n",
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
        if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
          auto deps = SNLDesignModeling::getCombinatorialOutputs(
              term.getSnlBitTerm());
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
            const auto& tt = tTerm.getSnlBitTerm()->getDesign()->getTruthTable(tTerm.getSnlBitTerm()->getOrderID());
            if (tt.isInitialized()) {
              tts.push_back(tt);
              // print deps
              for (const auto& d : tt.getDependencies()) {
                DEBUG_LOG("TT deps: %llu\n",
                          d);
              }
            }
          }
          bool inTermInTTDeps = false;
         for (const auto& tt : tts) {
              const auto& ttDeps = tt.getDependencies();            // expect std::vector<uint64_t>
              uint64_t orderID = term.getSnlBitTerm()->getOrderID(); // assume 0-based

              // If orderID is 1-based, uncomment:
              // if (orderID > 0) --orderID;

              for (size_t index = 0; index < ttDeps.size(); ++index) {
                  uint64_t d = ttDeps[index];

                  uint64_t blockMin = index * 64ULL;                // inclusive
                  uint64_t blockMax = blockMin + 64ULL;             // exclusive

                  // If orderID is before this block, nothing more to find
                  if (orderID < blockMin) {
                      break;
                  }

                  // Skip if orderID beyond this block
                  if (orderID >= blockMax) {
                      continue;
                  }

                  uint64_t localBit = orderID - blockMin;          // 0..63

                  //std::printf("TT deps: %llu\n", static_cast<unsigned long long>(d));
                  //std::printf("d = %llu, orderID = %llu, index = %zu, localBit = %llu\n",
                              // static_cast<unsigned long long>(d),
                              // static_cast<unsigned long long>(orderID),
                              // index,
                              // static_cast<unsigned long long>(localBit));

                  // Defensive: check localBit < 64
                  if (localBit >= 64) {
                      //std::fprintf(stderr, "localBit out of range: %llu\n",
                      //            static_cast<unsigned long long>(localBit));
                      break;
                  }

                  // Correct shift using 1ULL and parentheses
                  uint64_t mask = (1ULL << localBit);
                  //std::printf("mask = 0x%llx\n", static_cast<unsigned long long>(mask));

                  if ((d & mask) != 0ULL) {
                      inTermInTTDeps = true;
                  }

                  // we handled the block which contains orderID; stop scanning ttDeps
                  break;
              }

              if (inTermInTTDeps) {
                  // Found — act and break outer loop if that's desired
                  //std::puts("Found matching bit in this TT deps");
                  break;
              }
          }
          if (/*deps.empty() &&*/ !inTermInTTDeps) {
            outputsSet.insert(termId);
            DEBUG_LOG("Collecting output %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(), 
                      term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
          }
        }
      }
      continue;
    }
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
        if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                      term.getSnlBitTerm()) != seqBitTerms.end())
          outputsSet.insert(termId);
          DEBUG_LOG("Collecting seq output %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(), 
                      term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
      }
    }
  }
  outputs.clear();
  outputs.assign(outputsSet.begin(), outputsSet.end());
  return outputs;
}

void BuildPrimaryOutputClauses::collect() {
  inputs_ = collectInputs();
  sortInputs();
  for (const auto& input : inputs_) {
    inputsMap_[naja::DNL::get()->getDNLTerminalFromID(input).getFullPathIDs()] = input;
  }
  outputs_ = collectOutputs();
  sortOutputs();
  for (const auto& output : outputs_) {
    outputsMap_[naja::DNL::get()->getDNLTerminalFromID(output).getFullPathIDs()] = output;
    printf("Output collected: %s\n",
          naja::DNL::get()->getDNLTerminalFromID(output).getSnlBitTerm()->getName().getString().c_str());
  }
  POs_.resize(outputs_.size());
}

void BuildPrimaryOutputClauses::initVarNames() {
  termDNLID2varID_.resize(naja::DNL::get()->getDNLTerms().size(), (size_t) -1);
  for (size_t i = 0; i < inputs_.size(); ++i) {
    termDNLID2varID_[inputs_[i]] = i + 2;  // +2 to avoid 0 and 1 which are reserved for constants
  }
}

void BuildPrimaryOutputClauses::build() {
  naja::DNL::get();
  POs_.clear();
  POs_ = tbb::concurrent_vector<std::shared_ptr<BoolExpr>>(outputs_.size());
  initVarNames();
  // Init var names(counting on the fact that normalization happened before)

  // inputs_ = collectInputs();
  // sortInputs();
  // outputs_ = collectOutputs();
  // sortOutputs();
  size_t processedOutputs = 0;
  //tbb::task_arena arena(20);
  // init arena with automatic number of threads
  tbb::task_arena arena(20);
  auto processOutput = [&](size_t i) {
    DNLID out = outputs_[i];
    printf("Procssing output %zu/%zu: %s\n",
                               ++processedOutputs, outputs_.size(),
                               get()
                                   ->getDNLTerminalFromID(out)
                                   .getSnlBitTerm()
                                   ->getName()
                                   .getString()
                                   .c_str());
                        
                        SNLLogicCloud cloud(out, inputs_, outputs_);
                        cloud.compute();
                        // //cloud.getTruthTable().print();
                        // std::vector<DNLID> test1;
                        // std::vector<DNLID> test2;
                        // for (auto in : cloud.getAllInputs()) {
                        //   printf("Input in tree cloud: %lu\n", in);
                        //   // if (in >= cloud.getInputs().size()) {
                        //   //   printf("size of inputs in cloud: %lu\n", cloud.getInputs().size());
                        //   //   //assert(false && "Input in cloud is out of range");
                        //   // }
                        //  test1.push_back(in);
                        // }
                        // for (auto in : cloud.getInputs()) {
                        //   printf("Input in cloud: %lu\n", in);
                        //   test2.push_back(in);
                        // }
                        // std::sort(test1.begin(), test1.end());
                        // std::sort(test2.begin(), test2.end());
                        // assert(test1 == test2);
                        //std::vector<std::string> varNames;
                        /*for (auto input : cloud.getInputs()) {
                          DNLTerminalFull term = get()->getDNLTerminalFromID(input);
                          if (term.getSnlTerm() != nullptr) {
                            auto net = term.getSnlTerm()->getNet();
                            if (net != nullptr) {
                              if (net->isConstant0()) {
                                varNames.push_back("0");
                                continue;
                              } else if (net->isConstant1()) {
                                varNames.push_back("1");
                                continue;
                              }
                            }
                            auto model = const_cast<SNLDesign*>(
                                term.getSnlBitTerm()->getDesign());
                            auto tt = model->getTruthTable(term.getSnlBitTerm()->getOrderID());
                            if (tt.isInitialized()) {
                              if (tt.all0()) {
                                varNames.push_back("0");
                                continue;
                              } else if (tt.all1()) {
                                varNames.push_back("1");
                                continue;
                              }
                            }
                          }
                          // find the index of input in inputs_
                          auto it = std::find(inputs_.begin(), inputs_.end(), input);
                          // printf("Input: %s\n",
                          //       get()->getDNLTerminalFromID(input).getSnlBitTerm()->getName().getString().c_str());
                          // printf("Model: %s\n",
                          //       get()->getDNLTerminalFromID(input).getSnlBitTerm()->getDesign()->getName().getString().c_str());
                          assert(it != inputs_.end());
                          size_t index = std::distance(inputs_.begin(), it);
                          varNames.push_back(std::to_string(index + 2)); // +2 to avoid 0 and 1 which are reserved for constants
                        }*/

                        assert(cloud.getTruthTable().isInitialized());
                        //DEBUG_LOG("Truth Table: %s\n",
                        //          cloud.getTruthTable().print().c_str());
                        /*std::shared_ptr<BoolExpr> expr = Tree2BoolExpr::convert(
                            cloud.getTruthTable(), varNames);*/
                        //BoolExpr::getMutex().lock();
                        // if (POs_.size() - 1 < i) {
                        //   for (size_t j = POs_.size(); j <= i; ++j) {
                        //     POs_.push_back(nullptr);
                        //   }
                        // }
                        assert(POs_.size()  - 1 >= i);
                        POs_[i] = Tree2BoolExpr::convert(
                            cloud.getTruthTable(), termDNLID2varID_);
                        //BoolExpr::getMutex().unlock();
                        //printf("size of expr: %lu\n", POs_.back()->size());
                      };
  
  if (getenv("KEPLER_NO_MT")) {
    for (size_t i = 0; i < outputs_.size(); ++i) {
      processOutput(i);
    }
  } else {
    tbb::parallel_for(tbb::blocked_range<DNLID>(0, outputs_.size()),
                      [&](const tbb::blocked_range<DNLID>& r) {
                        for (DNLID i = r.begin(); i < r.end(); ++i) {
                          processOutput(i);
                      }});
  }
  destroy();  // Clean up DNL instance
}

void BuildPrimaryOutputClauses::setInputs2InputsIDs() {
  inputs2inputsIDs_.clear();
  for (const auto& input : inputs_) {
    std::vector<NLID::DesignObjectID> path;
    DNLInstanceFull currentInstance =
        get()->getDNLTerminalFromID(input).getDNLInstance();
    while (currentInstance.isTop() == false) {
      path.push_back(currentInstance.getSNLInstance()->getID());
      currentInstance = currentInstance.getParentInstance();
    }
    std::reverse(path.begin(), path.end());
    std::vector<NLID::DesignObjectID> termIDs;
    termIDs.push_back(
        get()->getDNLTerminalFromID(input).getSnlBitTerm()->getID());
    termIDs.push_back(
        get()->getDNLTerminalFromID(input).getSnlBitTerm()->getBit());
    inputs2inputsIDs_[input] = std::make_pair(path, termIDs);
  }
}

void BuildPrimaryOutputClauses::setOutputs2OutputsIDs() {
  outputs2outputsIDs_.clear();
  for (const auto& output : outputs_) {
    std::vector<NLID::DesignObjectID> path;
    DNLInstanceFull currentInstance =
        get()->getDNLTerminalFromID(output).getDNLInstance();
    while (currentInstance.isTop() == false) {
      path.push_back(currentInstance.getSNLInstance()->getID());
      currentInstance = currentInstance.getParentInstance();
    }
    std::reverse(path.begin(), path.end());
    std::vector<NLID::DesignObjectID> termIDs;
    termIDs.push_back(
        get()->getDNLTerminalFromID(output).getSnlBitTerm()->getID());
    termIDs.push_back(
        get()->getDNLTerminalFromID(output).getSnlBitTerm()->getBit());
    outputs2outputsIDs_[output] = std::make_pair(path, termIDs);
  }
}

void BuildPrimaryOutputClauses::sortInputs() {
  // Sort based on inputs2inputsIDs_ content
  std::sort(inputs_.begin(), inputs_.end(),
            [this](const DNLID& a, const DNLID& b) {
              return inputs2inputsIDs_[a].first < inputs2inputsIDs_[b].first;
            });
}

void BuildPrimaryOutputClauses::sortOutputs() {
  // Sort based on outputs2outputsIDs_ content
  std::sort(
      outputs_.begin(), outputs_.end(), [this](const DNLID& a, const DNLID& b) {
        return outputs2outputsIDs_[a].first < outputs2outputsIDs_[b].first;
      });
}