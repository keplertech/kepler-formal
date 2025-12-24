#include "ScopeExtraction.h"
#include <stack>

void ScopeExtraction::collectVerificationScopes() {
  // Find leaf models(model that contain leaf only)
  // DFS
  std::stack<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>> stack;
  std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*> toCompare;
  toCompare.first = top0_;
  toCompare.second = top1_;
  auto currentEqualDesigns = equalDesigns_;
  while ((currentEqualDesigns.size() > equalDesigns_.size()) ||
         equalDesigns_.empty()) {
    stack.push(toCompare);
    currentEqualDesigns = equalDesigns_;
    designsToVerify_.clear();
    while (!stack.empty()) {
      std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*> toCompareNew =
          stack.top();
      stack.pop();
      auto design0 = toCompareNew.first;
      auto design1 = toCompareNew.second;
      if (design0->getInstances().size() == design1->getInstances().size()) {
        bool modelsAreEqual = true;
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
        if (instancesVector1.size() != instancesVector0.size()) {
          // Different number of hierarchical instances, for sure different
          modelsAreEqual = false;
        }
        if (modelsAreEqual) {
          for (size_t i = 0; i < instancesVector1.size(); ++i) {
            if ((instancesVector0[i]->getModel()->getID() !=
                 instancesVector1[i]->getModel()->getID()) ||
                (instancesVector0[i]->getID() !=
                 instancesVector1[i]->getID())) {
              // Different instance ID or model ID of an instance the same
              // position means that the models differ
              modelsAreEqual = false;
            } else {
              std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*> pair;
              pair.first = instancesVector0[i]->getModel();
              pair.second = instancesVector1[i]->getModel();
              if (equalDesigns_.find(pair) == equalDesigns_.end()) {
                modelsAreEqual = false;
                // Only in one case we need to go down the hierarchy
                // When the 2 models have the same interfaces but were not
                // yes checked
                stack.push(pair);
              }
            }
          }
        }
        if (modelsAreEqual) {
          // Hierarchical instances passed all the checks, now we need to
          // compare the primitive instances
          auto primitiveInstances0 = design0->getPrimitiveInstances();
          std::vector<naja::NL::SNLInstance*> primitivesVector0;
          for (auto instance : primitiveInstances0) {
            primitivesVector0.push_back(instance);
          }
          auto primitiveInstances1 = design1->getPrimitiveInstances();
          std::vector<naja::NL::SNLInstance*> primitivesVector1;
          for (auto instance : primitiveInstances1) {
            primitivesVector1.push_back(instance);
          }
          bool primitivesAreEqual = true;
          for (size_t i = 0; i < primitivesVector0.size(); ++i) {
            if ((primitivesVector0[i]->getModel()->getID() !=
                 primitivesVector1[i]->getModel()->getID()) ||
                (primitivesVector0[i]->getID() !=
                 primitivesVector1[i]->getID())) {
              modelsAreEqual = false;
              break;
            }
          }
        }
        if (modelsAreEqual) {
          equalDesigns_.insert(toCompareNew);
        } else {
          // models are not equal, need to compare with formal verification
          designsToVerify_.insert(toCompareNew);
        }
      }
    }
  }
}
