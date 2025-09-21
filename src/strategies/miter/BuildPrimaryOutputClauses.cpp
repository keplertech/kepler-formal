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
        inputs.push_back(termId);
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
          if (deps.empty()) {
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
      if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                    term.getSnlBitTerm()) != seqBitTerms.end())
        inputs.push_back(termId);
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
  auto dnl = get();
  DNLInstanceFull top = dnl->getTop();

  for (DNLID termId = top.getTermIndexes().first;
       termId != DNLID_MAX && termId <= top.getTermIndexes().second; termId++) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
    if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Input)
      outputs.push_back(termId);
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
        outputs.push_back(termId);
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
          if (deps.empty()) {
            outputs.push_back(termId);
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
      if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                    term.getSnlBitTerm()) != seqBitTerms.end())
        outputs.push_back(termId);
    }
  }

  std::set<DNLID> outputSet(outputs.begin(), outputs.end());
  outputs.clear();
  outputs.assign(outputSet.begin(), outputSet.end());
  return outputs;
}

void BuildPrimaryOutputClauses::build() {
  POs_.clear();
  inputs_ = collectInputs();
  sortInputs();
  outputs_ = collectOutputs();
  sortOutputs();
  size_t processedOutputs = 0;
  tbb::task_arena arena(20);
  auto processOutput = [&](DNLID out) {
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

                        std::vector<std::string> varNames;
                        for (auto input : cloud.getInputs()) {
                          varNames.push_back(std::to_string(input));
                        }

                        assert(cloud.getTruthTable().isInitialized());
                        //DEBUG_LOG("Truth Table: %s\n",
                        //          cloud.getTruthTable().print().c_str());
                        /*std::shared_ptr<BoolExpr> expr = Tree2BoolExpr::convert(
                            cloud.getTruthTable(), varNames);*/
                        
                        POs_.push_back(Tree2BoolExpr::convert(
                            cloud.getTruthTable(), varNames));
                        //printf("size of expr: %lu\n", POs_.back()->size());
                      };
  if (getenv("KEPLER_NO_MT")) {
    for (DNLID i = 0; i < outputs_.size(); ++i) {
      auto out = outputs_[i];
      processOutput(out);
    }
  } else {
    tbb::parallel_for(tbb::blocked_range<DNLID>(0, outputs_.size()),
                      [&](const tbb::blocked_range<DNLID>& r) {
                        for (DNLID i = r.begin(); i < r.end(); ++i) {
                          auto out = outputs_[i];
                          processOutput(out);
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