#include "MiterStrategy.h"
#include "BoolExpr.h"
#include "DNL.h"
#include "SNLDesignModeling.h"

using namespace KEPLER_FORMAL;
using namespace naja::DNL;
using namespace naja::NL;

std::vector<DNLID> MiterStrategy::collectInputs() {
  std::vector<DNLID> inputs;
  auto dnl = get();
  // Add top ports
  DNLInstanceFull top = dnl->getTop();
  for (DNLID termId = top.getTermIndexes().first;
       termId != DNLID_MAX and termId <= top.getTermIndexes().second;
       termId++) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
    if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
      inputs.push_back(termId);
    }
  }
  // Add Sequential inputs
  for (DNLID leaf : dnl->getLeaves()) {
    DNLInstanceFull instance = dnl->getDNLInstanceFromID(leaf);
    bool isSequential = false;
    std::vector<SNLBitTerm*> seqBitTerms;
    for (DNLID termId = top.getTermIndexes().first;
         termId != DNLID_MAX and termId <= top.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (SNLDesignModeling::getClockRelatedOutputs(term.getSnlBitTerm())
              .size() > 0) {
        isSequential = true;
        for (auto bitTerm :
             SNLDesignModeling::getClockRelatedOutputs(term.getSnlBitTerm())) {
          seqBitTerms.push_back(bitTerm);
        }
        inputs.push_back(termId);
      }
    }
    if (!isSequential) {
      continue;
    }
    for (DNLID termId = top.getTermIndexes().first;
         termId != DNLID_MAX and termId <= top.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                    term.getSnlBitTerm()) != seqBitTerms.end()) {
        inputs.push_back(termId);
      }
    }
  }
  return inputs;
}

std::vector<DNLID> MiterStrategy::collectOutputs() {
  std::vector<DNLID> outputs;
  auto dnl = get();
  // Add top ports
  DNLInstanceFull top = dnl->getTop();
  for (DNLID termId = top.getTermIndexes().first;
       termId != DNLID_MAX and termId <= top.getTermIndexes().second;
       termId++) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
    if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Input) {
      outputs.push_back(termId);
    }
  }
  // Add Sequential outputs
  for (DNLID leaf : dnl->getLeaves()) {
    DNLInstanceFull instance = dnl->getDNLInstanceFromID(leaf);
    bool isSequential = false;
    std::vector<SNLBitTerm*> seqBitTerms;
    for (DNLID termId = top.getTermIndexes().first;
         termId != DNLID_MAX and termId <= top.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (SNLDesignModeling::getClockRelatedInputs(term.getSnlBitTerm())
              .size() > 0) {
        isSequential = true;
        for (auto bitTerm :
             SNLDesignModeling::getClockRelatedInputs(term.getSnlBitTerm())) {
          seqBitTerms.push_back(bitTerm);
        }
        outputs.push_back(termId);
      }
    }
    if (!isSequential) {
      continue;
    }
    for (DNLID termId = top.getTermIndexes().first;
         termId != DNLID_MAX and termId <= top.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                    term.getSnlBitTerm()) != seqBitTerms.end()) {
        outputs.push_back(termId);
      }
    }
  }
  return outputs;
}

void MiterStrategy::build() {
  // Collect inputs
}
