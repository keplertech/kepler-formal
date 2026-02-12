#include "ScopeExtraction.h"
#include <tbb/concurrent_unordered_set.h>
#include <stack>
#include "DNL.h"
#include "NLUniverse.h"
#include "RemoveLoadlessLogic.h"
#include "SNLDesign.h"
#include "SNLLogicCone.h"
#include <tbb/enumerable_thread_specific.h>
#include "SNLBitNet.h"

// #define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace naja::NL;
using namespace naja::DNL;

void ScopeExtraction::collectVerificationScopes() {
  // Find leaf models(model that contain leaf only)
  // DFS
  std::stack<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>> stack;
  std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*> toCompare;
  toCompare.first = top0_;
  toCompare.second = top1_;
  stack.push(toCompare);
  while (!stack.empty()) {
    std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*> toCompareNew =
        stack.top();
    stack.pop();
    auto design0 = toCompareNew.first;
    auto design1 = toCompareNew.second;
    bool modelsAreEqual = true;
    // First check, same number of instances
    if (design0->getInstances().size() == design1->getInstances().size()) {
      // We have same number of instances moving on to next checks
      // 1. Check same child instances
      DEBUG_LOG(" - Comparing models %s and %s\n",
                design0->getName().getString().c_str(),
                design1->getName().getString().c_str());
      std::vector<naja::NL::NLID::DesignObjectID> childInstances0;
      for (auto instance : design0->getInstances()) {
        childInstances0.push_back(instance->getID());
      }
      std::vector<naja::NL::NLID::DesignObjectID> childInstances1;
      for (auto instance : design1->getInstances()) {
        childInstances1.push_back(instance->getID());
      }
      if (childInstances0 != childInstances1) {
        DEBUG_LOG(" - Child instances are different.\n");
        modelsAreEqual = false;
      }
      if (modelsAreEqual) {
        DEBUG_LOG(" - Child instances are the same -> Comparing nets.\n");
        // Same child instances, now we will check all nets
        std::vector<naja::NL::SNLBitNet*> nets0;
        for (auto net : design0->getBitNets()) {
          nets0.push_back(net);
        }
        std::vector<naja::NL::SNLBitNet*> nets1;
        for (auto net : design1->getBitNets()) {
          nets1.push_back(net);
        }
        if (nets0.size() != nets1.size()) {
          DEBUG_LOG(" - Different number of nets.\n");
          modelsAreEqual = false;
        }
        if (modelsAreEqual) {
          for (size_t i = 0; i < nets0.size(); ++i) {
            DEBUG_LOG(" - Comparing net %s and %s\n",
                      nets0[i]->getName().getString().c_str(),
                      nets1[i]->getName().getString().c_str());
            DEBUG_LOG("Comparing bit terms.\n");
            auto net0 = nets0[i];
            auto net1 = nets1[i];
            // Check drivers
            std::vector<SNLBitTerm*> bitTerms0;
            for (auto bitterm : net0->getBitTerms()) {
              bitTerms0.push_back(bitterm);
            }
            std::vector<SNLBitTerm*> bitTerms1;
            for (auto bitterm : net1->getBitTerms()) {
              bitTerms1.push_back(bitterm);
            }
            if (bitTerms0.size() != bitTerms1.size()) {
              modelsAreEqual = false;
              break;
            }
            for (size_t j = 0; j < bitTerms0.size(); ++j) {
              auto term0 = bitTerms0[j];
              auto term1 = bitTerms1[j];
              if (term0->getID() != term1->getID()) {
                modelsAreEqual = false;
                break;
              }
            }
            if (modelsAreEqual) {
              DEBUG_LOG("Comparing inst terms.\n");
              // bit terms are same, now check ins terms
              std::vector<SNLInstTerm*> instTerms0;
              for (auto instterm : net0->getInstTerms()) {
                instTerms0.push_back(instterm);
              }
              std::vector<SNLInstTerm*> instTerms1;
              for (auto instterm : net1->getInstTerms()) {
                instTerms1.push_back(instterm);
              }
              if (instTerms0.size() != instTerms1.size()) {
                modelsAreEqual = false;
                break;
              }
              for (size_t j = 0; j < instTerms0.size(); ++j) {
                auto term0 = instTerms0[j];
                auto term1 = instTerms1[j];
                if (term0->getInstance()->getID() !=
                        term1->getInstance()->getID() ||
                    term0->getBitTerm()->getID() !=
                        term1->getBitTerm()->getID()) {
                  modelsAreEqual = false;
                  break;
                }
              }
            }
          }
        }
      }
    } else {
      DEBUG_LOG(" - Different number of instances.\n");
      modelsAreEqual = false;
    }
    if (!modelsAreEqual) {
      designsToVerify_.insert(toCompareNew);
    } else {
      // Add checks for children instances in next interation
      std::vector<naja::NL::SNLInstance*> instancesVector0;
      for (auto instance : design0->getInstances()) {
        if (instance->isPrimitive()) {
          continue;
        }
        instancesVector0.push_back(instance);
      }
      std::vector<naja::NL::SNLInstance*> instancesVector1;
      for (auto instance : design1->getInstances()) {
        if (instance->isPrimitive()) {
          continue;
        }
        instancesVector1.push_back(instance);
      }
      for (size_t i = 0; i < instancesVector1.size(); ++i) {
        std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*> pair;
        pair.first = instancesVector0[i]->getModel();
        pair.second = instancesVector1[i]->getModel();
        stack.push(pair);
      }
    }
  }
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
    size_t readerCount = 0;
    std::vector<DNLID> readersVector;
    readersVector.assign(readers.begin(), readers.end());   
    KEPLER_FORMAL::SNLLogicCone cone(readersVector, pis[i]);
    cone.run();
    for (const auto& isoID : cone.getConeIsoIDs()) {
      isosToKeep.insert(isoID);
    }
    // if (getenv("KEPLER_NO_MT") || true) {
    //   for (const auto& reader : readersVector) {
    //     printf("Processing reader %zu / %zu\r", readerCount + 1,
    //            readersVector.size());
    //     readerCount++;
    //     KEPLER_FORMAL::SNLLogicCone cone(reader, pis[i]);
    //     cone.initConeIsos(isosToKeep);
    //     cone.run();
    //     printf("inserting isos %zu\n", cone.getConeIsoIDs().size());
    //     // for (const auto& isoID : cone.getConeIsoIDs()) {
    //     //   isosToKeep.insert(isoID);
    //     // }
    //     printf("inserted\n");
    //   }
    // } else {
    //   tbb::enumerable_thread_specific<std::set<DNLID>> localIsosToKeep;
    //   tbb::parallel_for(tbb::blocked_range<DNLID>(0, readersVector.size()),
    //                     [&](const tbb::blocked_range<DNLID>& r) {
    //                       for (DNLID j = r.begin(); j < r.end(); ++j) {
    //                         KEPLER_FORMAL::SNLLogicCone cone(readersVector[j], pis[i]);
    //                         //cone.initConeIsos(isosToKeep);
    //                         cone.run();
    //                         auto& localIsosToKeepSet = localIsosToKeep.local();
    //                         localIsosToKeepSet.insert(
    //                             cone.getConeIsoIDs().begin(),
    //                             cone.getConeIsoIDs().end());
    //                       }
    //                     });
    //   for (const auto& localSet : localIsosToKeep) {
    //     for (const auto& isoID : localSet) {
    //       isosToKeep.insert(isoID);
    //     }
    //   }
    // }
    naja::NAJA_OPT::LoadlessLogicRemover remover;
    auto loadlessInstances = remover.getLoadlessInstances(*naja::DNL::get(), isosToKeep);
    remover.removeLoadlessInstances(top, loadlessInstances);
  }
}
