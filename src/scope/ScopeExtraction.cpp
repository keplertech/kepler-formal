#include "ScopeExtraction.h"
#include <tbb/concurrent_unordered_set.h>
#include <algorithm>
#include <stack>
#include "DNL.h"
#include "NLUniverse.h"
#include "RemoveLoadlessLogic.h"
#include "SNLBitNet.h"
#include "SNLBitTerm.h"
#include "SNLBusNetBit.h"
#include "SNLDesign.h"
#include "SNLInstParameter.h"
#include "SNLInstTerm.h"
#include "SNLInstance.h"
#include "SNLLogicCone.h"
#include "SNLParameter.h"
#include <tbb/enumerable_thread_specific.h>

// #define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace naja::NL;
using namespace naja::DNL;

namespace {

// Name-based key of one bit of a term or net. Bits of a bus share the bus
// name and differ by index; scalars use bit 0.
std::string bitKey(const NLName& name, NLID::Bit bit) {
  return name.getString() + "[" + std::to_string(bit) + "]";
}

NLID::Bit netBit(const SNLBitNet* net) {
  if (auto* busBit = dynamic_cast<const SNLBusNetBit*>(net)) {
    return busBit->getBit();
  }
  return 0;
}

std::string termKey(const SNLBitTerm* term) {
  // A term without a name (anonymous primitive ports) is identified by its
  // position, which the lowering assigns deterministically.
  if (term->getName().empty()) {
    return "#" + std::to_string(term->getID()) + "[" +
           std::to_string(term->getBit()) + "]";
  }
  return bitKey(term->getName(), term->getBit());
}

template <class Container>
bool sameSorted(Container& left, Container& right) {
  if (left.size() != right.size()) {
    return false;
  }
  std::sort(left.begin(), left.end());
  std::sort(right.begin(), right.end());
  return left == right;
}

}  // namespace

bool ScopeExtraction::sameInterface(SNLDesign* design0, SNLDesign* design1) {
  std::vector<std::string> keys0;
  std::vector<std::string> keys1;
  for (auto* term : design0->getBitTerms()) {
    keys0.push_back(termKey(term) + ":" +
                    std::to_string(static_cast<int>(term->getDirection())));
  }
  for (auto* term : design1->getBitTerms()) {
    keys1.push_back(termKey(term) + ":" +
                    std::to_string(static_cast<int>(term->getDirection())));
  }
  return sameSorted(keys0, keys1);
}

bool ScopeExtraction::sameParameters(SNLDesign* design0, SNLDesign* design1) {
  std::vector<std::string> keys0;
  std::vector<std::string> keys1;
  for (auto* parameter : design0->getParameters()) {
    keys0.push_back(parameter->getName().getString() + "=" +
                    (parameter->hasDefaultValue() ? parameter->getValue() : ""));
  }
  for (auto* parameter : design1->getParameters()) {
    keys1.push_back(parameter->getName().getString() + "=" +
                    (parameter->hasDefaultValue() ? parameter->getValue() : ""));
  }
  return sameSorted(keys0, keys1);
}

bool ScopeExtraction::sameInstanceParameters(SNLInstance* instance0,
                                             SNLInstance* instance1) {
  std::vector<std::string> keys0;
  std::vector<std::string> keys1;
  for (auto* parameter : instance0->getInstParameters()) {
    keys0.push_back(parameter->getName().getString() + "=" + parameter->getValue());
  }
  for (auto* parameter : instance1->getInstParameters()) {
    keys1.push_back(parameter->getName().getString() + "=" + parameter->getValue());
  }
  return sameSorted(keys0, keys1);
}

bool ScopeExtraction::sameLeafModels(SNLDesign* model0, SNLDesign* model1) {
  if (model0 == model1) {
    return true;
  }
  if (!model0->isLeaf() || !model1->isLeaf()) {
    return false;
  }
  if (model0->getName() != model1->getName()) {
    return false;
  }
  if (model0->getType().getString() != model1->getType().getString()) {
    return false;
  }
  return sameInterface(model0, model1) && sameParameters(model0, model1);
}

bool ScopeExtraction::pairInstances(SNLDesign* design0,
                                    SNLDesign* design1,
                                    InstancePairs& pairs) {
  pairs.clear();
  if (design0->getInstances().size() != design1->getInstances().size()) {
    return false;
  }
  std::map<NLName, SNLInstance*> named0;
  std::map<NLName, SNLInstance*> named1;
  std::vector<SNLInstance*> unnamed0;
  std::vector<SNLInstance*> unnamed1;
  for (auto* instance : design0->getInstances()) {
    if (instance->getName().empty()) {
      unnamed0.push_back(instance);
    } else {
      named0.emplace(instance->getName(), instance);
    }
  }
  for (auto* instance : design1->getInstances()) {
    if (instance->getName().empty()) {
      unnamed1.push_back(instance);
    } else {
      named1.emplace(instance->getName(), instance);
    }
  }
  if (named0.size() != named1.size() || unnamed0.size() != unnamed1.size()) {
    return false;
  }
  for (const auto& [name, instance0] : named0) {
    auto it = named1.find(name);
    if (it == named1.end()) {
      return false;
    }
    pairs.emplace_back(instance0, it->second);
  }
  for (size_t i = 0; i < unnamed0.size(); ++i) {
    pairs.emplace_back(unnamed0[i], unnamed1[i]);
  }
  return true;
}

bool ScopeExtraction::sameNets(
    SNLDesign* design0,
    SNLDesign* design1,
    const std::map<SNLInstance*, SNLInstance*>& instanceMap) {
  if (design0->getBitNets().size() != design1->getBitNets().size()) {
    return false;
  }
  std::map<SNLInstance*, SNLInstance*> reverseMap;
  for (const auto& [instance0, instance1] : instanceMap) {
    reverseMap.emplace(instance1, instance0);
  }
  // Pair nets like instances: by name where the RTL named them, in creation
  // order where the lowering created them.
  std::map<std::string, SNLBitNet*> named0;
  std::map<std::string, SNLBitNet*> named1;
  std::vector<SNLBitNet*> unnamed0;
  std::vector<SNLBitNet*> unnamed1;
  for (auto* net : design0->getBitNets()) {
    if (net->getName().empty()) {
      unnamed0.push_back(net);
    } else {
      named0.emplace(bitKey(net->getName(), netBit(net)), net);
    }
  }
  for (auto* net : design1->getBitNets()) {
    if (net->getName().empty()) {
      unnamed1.push_back(net);
    } else {
      named1.emplace(bitKey(net->getName(), netBit(net)), net);
    }
  }
  if (named0.size() != named1.size() || unnamed0.size() != unnamed1.size()) {
    return false;
  }
  std::vector<std::pair<SNLBitNet*, SNLBitNet*>> pairs;
  for (const auto& [key, net0] : named0) {
    auto it = named1.find(key);
    if (it == named1.end()) {
      return false;
    }
    pairs.emplace_back(net0, it->second);
  }
  for (size_t i = 0; i < unnamed0.size(); ++i) {
    pairs.emplace_back(unnamed0[i], unnamed1[i]);
  }
  for (const auto& [net0, net1] : pairs) {
    if (!(net0->getType() == net1->getType())) {
      return false;
    }
    // Connections to the module ports.
    std::vector<std::string> ports0;
    std::vector<std::string> ports1;
    for (auto* term : net0->getBitTerms()) {
      ports0.push_back(termKey(term));
    }
    for (auto* term : net1->getBitTerms()) {
      ports1.push_back(termKey(term));
    }
    if (!sameSorted(ports0, ports1)) {
      return false;
    }
    // Connections to instance pins, expressed on design 0's instances so both
    // sides use one vocabulary.
    std::vector<std::string> pins0;
    std::vector<std::string> pins1;
    for (auto* instTerm : net0->getInstTerms()) {
      pins0.push_back(std::to_string(instTerm->getInstance()->getID()) + ":" +
                      termKey(instTerm->getBitTerm()));
    }
    for (auto* instTerm : net1->getInstTerms()) {
      auto it = reverseMap.find(instTerm->getInstance());
      if (it == reverseMap.end()) {
        return false;
      }
      pins1.push_back(std::to_string(it->second->getID()) + ":" +
                      termKey(instTerm->getBitTerm()));
    }
    if (!sameSorted(pins0, pins1)) {
      return false;
    }
  }
  return true;
}

const ScopeExtraction::Comparison& ScopeExtraction::compareModels(
    SNLDesign* design0, SNLDesign* design1) {
  const ModelPair key(design0, design1);
  auto found = comparisons_.find(key);
  if (found != comparisons_.end()) {
    return found->second;
  }
  // Reserve the slot first so a malformed cyclic hierarchy cannot recurse
  // forever; the hierarchy is a DAG so the slot is final after this call.
  comparisons_[key];
  Comparison result;
  if (design0 == design1) {
    result.ownIdentical = true;
    result.identical = true;
  } else if (design0->isLeaf() || design1->isLeaf()) {
    result.ownIdentical = sameLeafModels(design0, design1);
    result.identical = result.ownIdentical;
  } else {
    DEBUG_LOG(" - Comparing models %s and %s\n",
              design0->getName().getString().c_str(),
              design1->getName().getString().c_str());
    bool own = sameInterface(design0, design1) && sameParameters(design0, design1);
    bool childrenIdentical = true;
    InstancePairs pairs;
    if (own && !pairInstances(design0, design1, pairs)) {
      DEBUG_LOG(" - Instances do not pair.\n");
      own = false;
    }
    std::map<SNLInstance*, SNLInstance*> instanceMap;
    if (own) {
      for (const auto& [instance0, instance1] : pairs) {
        instanceMap.emplace(instance0, instance1);
        if (!sameInstanceParameters(instance0, instance1)) {
          own = false;
          break;
        }
        auto* model0 = instance0->getModel();
        auto* model1 = instance1->getModel();
        if (model0->isLeaf() || model1->isLeaf()) {
          // Leaf models are part of this module's own implementation.
          if (!sameLeafModels(model0, model1)) {
            own = false;
            break;
          }
          continue;
        }
        // A child whose interface changed changes this module's connections
        // too, and the pair could not be verified as a scope of its own.
        if (!sameInterface(model0, model1)) {
          own = false;
          break;
        }
        // Bottom-up: the child pair is decided before this one.
        if (!compareModels(model0, model1).identical) {
          childrenIdentical = false;
        }
      }
    }
    if (own && !sameNets(design0, design1, instanceMap)) {
      DEBUG_LOG(" - Nets differ.\n");
      own = false;
    }
    result.ownIdentical = own;
    result.identical = own && childrenIdentical;
  }
  if (result.ownIdentical) {
    ++identicalOwnPairs_;
  }
  Comparison& slot = comparisons_[key];
  slot = result;
  return slot;
}

void ScopeExtraction::collectScopesBelow(SNLDesign* design0,
                                         SNLDesign* design1,
                                         std::set<ModelPair>& visited) {
  if (!visited.insert(ModelPair(design0, design1)).second) {
    return;
  }
  const Comparison& comparison = compareModels(design0, design1);
  if (!comparison.ownIdentical) {
    // This module's own implementation differs: verify it as a whole, its
    // differing descendants included.
    designsToVerify_.insert(ModelPair(design0, design1));
    return;
  }
  if (comparison.identical) {
    return;
  }
  // Own wiring identical, some descendant differs: descend to find it.
  InstancePairs pairs;
  pairInstances(design0, design1, pairs);
  for (const auto& [instance0, instance1] : pairs) {
    auto* model0 = instance0->getModel();
    auto* model1 = instance1->getModel();
    if (model0->isLeaf() || model1->isLeaf()) {
      continue;
    }
    collectScopesBelow(model0, model1, visited);
  }
}

void ScopeExtraction::collectVerificationScopes() {
  comparisons_.clear();
  designsToVerify_.clear();
  identicalOwnPairs_ = 0;
  if (top0_ == nullptr || top1_ == nullptr) {
    return;
  }
  // Bottom-up comparison: every model pair is decided after its children,
  // and each pair is compared once however many instances share it.
  compareModels(top0_, top1_);
  std::set<ModelPair> visited;
  collectScopesBelow(top0_, top1_, visited);
}

void ScopeExtraction::cleanVerificationScopes(
    const std::vector<naja::DNL::DNLID>& pis0,
    const std::vector<naja::DNL::DNLID>& pis1) {
  std::vector<naja::NL::SNLDesign*> tops;
  tops.push_back(top0_);
  tops.push_back(top1_);
  std::vector<std::vector<naja::DNL::DNLID>> pis;
  pis.push_back(pis0);
  pis.push_back(pis1);
  for (size_t i = 0; i < tops.size(); ++i) {
    naja::DNL::destroy();
    // Setting top to top0
    auto top = tops[i];
    naja::NL::NLUniverse::get()->setTopDesign(top);
    // Collecting all scopes to verify under top0
    std::set<const naja::NL::SNLDesign*> scopes;
    for (const auto& scope : designsToVerify_) {
      naja::NL::SNLDesign* design = i == 0 ? scope.first : scope.second;
      scopes.insert(design);
    }
    // Collecting all leaves who are under model is in scopes
    std::set<naja::DNL::DNLID> scopeLeaves;
    auto dnl = naja::DNL::get();
    tbb::concurrent_unordered_set<naja::DNL::DNLID> isosToKeep;
    for (const auto& leaf : dnl->getLeaves()) {
      const naja::NL::SNLDesign* model =
          dnl->getDNLInstanceFromID(leaf).getSNLModel();
      naja::DNL::DNLInstanceFull currentInstance =
          dnl->getDNLInstanceFromID(leaf);
      while (currentInstance.isTop() == false) {
        if (scopes.find(model) != scopes.end()) {
          break;
        }
        currentInstance = currentInstance.getParentInstance();
        model = currentInstance.getSNLModel();
        // keep all isos on terminals of model
      }
      if (scopes.find(model) != scopes.end()) {
        scopeLeaves.insert(leaf);
        for (DNLID termId = currentInstance.getTermIndexes().first;
             termId != DNLID_MAX && termId <= currentInstance.getTermIndexes().second;
             termId++) {
          const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
          DNLID isoID = term.getIsoID();
          if (isoID != DNLID_MAX) {
            isosToKeep.insert(isoID);
          }
        }
      }
    }
    std::set<DNLID> readers;
    for (const auto& leaf : dnl->getLeaves()) {
      // For each terminal collect the readers of it's iso
      DNLInstanceFull instance = dnl->getDNLInstanceFromID(leaf);
      for (DNLID termId = instance.getTermIndexes().first;
           termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
           termId++) {
        const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
        DNLID isoID = term.getIsoID();
        if (isoID != DNLID_MAX) {
          auto iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
          for (const auto& reader : iso.getReaders()) {
            DNLInstanceFull readerInstance =
                dnl->getDNLTerminalFromID(reader).getDNLInstance();
            DNLID readerLeafID = readerInstance.getID();
            if (scopeLeaves.find(readerLeafID) != scopeLeaves.end()) {
              readers.insert(reader);
            }
          }
        }
      }
    }
    // After this loop the dnl variable contain deleted DNL as cone destroy it in constructor for safety
    std::vector<DNLID> readersVector;
    readersVector.assign(readers.begin(), readers.end());
    KEPLER_FORMAL::SNLLogicCone cone(readersVector, pis[i]);
    cone.run();
    for (const auto& isoID : cone.getConeIsoIDs()) {
      isosToKeep.insert(isoID);
    }
    naja::NAJA_OPT::LoadlessLogicRemover remover;
    auto loadlessInstances = remover.getLoadlessInstances(*naja::DNL::get(), isosToKeep);
    remover.removeLoadlessInstances(top, loadlessInstances);
  }
}
