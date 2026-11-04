// src/scope/ScopeExtraction.h
// SPDX-License-Identifier: Apache-2.0
//
// Hierarchical scope extraction: pairs the two design hierarchies bottom-up
// and reports the topmost module pairs whose own implementation differs.
// Those pairs are the scopes a verification run has to compare; every other
// module pair is structurally identical and needs no proof.
//
// Small, test-friendly adjustment:
//  - add a protected default constructor so test subclasses can default-construct
//  - make internal members protected so test subclasses can set/inspect them
// The public API is preserved (existing constructor and getScopesToVerify()).
#ifndef SCOPE_EXTRACTION_H
#define SCOPE_EXTRACTION_H
#include "DNL.h"
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
class ScopeExtraction {
 public:
    // Existing public constructor (preserve original API)
    // LCOV_EXCL_START
    ScopeExtraction(naja::NL::SNLDesign* top0, naja::NL::SNLDesign* top1)
        : top0_(top0), top1_(top1) {}
        // LCOV_EXCL_STOP
    // LCOV_EXCL_START
    virtual ~ScopeExtraction() = default;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    // Compares the two hierarchies bottom-up and fills the scopes to verify:
    // the topmost model pairs whose own implementation differs. A pair whose
    // own wiring is identical but holds a differing child is not a scope; the
    // child is. A pair whose interface differs makes its parent the scope,
    // because the parent's connections to it differ and the pair could not be
    // compared on its own.
    void collectVerificationScopes();
    // LCOV_EXCL_START
    std::set<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>> getScopesToVerify() const {
        return designsToVerify_;
        // LCOV_EXCL_STOP
    }
    // Number of compared model pairs whose own implementation is identical
    // (interface, parameters, instances, leaf models and nets), regardless of
    // their children. Available after collectVerificationScopes().
    size_t getIdenticalModelPairCount() const { return identicalOwnPairs_; }
    // Number of model pairs compared in total.
    size_t getComparedModelPairCount() const { return comparisons_.size(); }
    // Provide extracted verification scopes for debugging purposes by:
    // 1 Collecting scopes via collectVerificationScopes()
    // 2 Collecting logic cones for all differ elements in the scopes with SNLLogicCone
    // 4 Deleting all elements not in the collected cones from the cloned tops
    void cleanVerificationScopes(const std::vector<naja::DNL::DNLID>& pis0, const std::vector<naja::DNL::DNLID>& pis1);
 protected:
    // Allow derived test helpers to default-construct and inspect internals.
    // Protected so only subclasses (e.g., test helpers) can use it.
    ScopeExtraction() = default;
    // Internal state made protected so test subclasses can set and inspect them.
    naja::NL::SNLDesign* top0_ = nullptr;
    naja::NL::SNLDesign* top1_ = nullptr;
    std::set<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>> designsToVerify_;

 private:
    struct Comparison {
        bool ownIdentical = false;  // interface, parameters, instances, nets
        bool identical = false;     // ownIdentical and every paired child identical
    };
    using ModelPair = std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>;
    using InstancePairs =
        std::vector<std::pair<naja::NL::SNLInstance*, naja::NL::SNLInstance*>>;
    // Bottom-up, memoized structural comparison of one model pair.
    const Comparison& compareModels(naja::NL::SNLDesign* design0,
                                    naja::NL::SNLDesign* design1);
    // Pairs the instances of two models: named instances by name, unnamed ones
    // (lowered gates and flops) in creation order. Returns false when no such
    // bijection exists.
    static bool pairInstances(naja::NL::SNLDesign* design0,
                              naja::NL::SNLDesign* design1,
                              InstancePairs& pairs);
    static bool sameInterface(naja::NL::SNLDesign* design0, naja::NL::SNLDesign* design1);
    static bool sameParameters(naja::NL::SNLDesign* design0, naja::NL::SNLDesign* design1);
    static bool sameInstanceParameters(naja::NL::SNLInstance* instance0,
                                       naja::NL::SNLInstance* instance1);
    // Leaf (primitive or black-box) models compare by name, type, interface
    // and parameters, never by model ID: shared primitive models are created
    // lazily, in first-use order across both designs.
    static bool sameLeafModels(naja::NL::SNLDesign* model0, naja::NL::SNLDesign* model1);
    static bool sameNets(
        naja::NL::SNLDesign* design0,
        naja::NL::SNLDesign* design1,
        const std::map<naja::NL::SNLInstance*, naja::NL::SNLInstance*>& instanceMap);
    // Collects the topmost own-differing pairs under a pair.
    void collectScopesBelow(naja::NL::SNLDesign* design0, naja::NL::SNLDesign* design1,
                            std::set<ModelPair>& visited);
    std::map<ModelPair, Comparison> comparisons_;
    size_t identicalOwnPairs_ = 0;
};
#endif // SCOPE_EXTRACTION_H
