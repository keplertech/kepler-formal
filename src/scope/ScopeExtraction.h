// src/scope/ScopeExtraction.h
// SPDX-License-Identifier: GPL-3.0-only
//
// Small, test-friendly adjustment:
//  - add a protected default constructor so test subclasses can default-construct
//  - make internal members protected so test subclasses can set/inspect them
// The public API is preserved (existing constructor and getScopesToVerify()).

#ifndef SCOPE_EXTRACTION_H
#define SCOPE_EXTRACTION_H

#include "DNL.h"
#include <set>
#include <utility>

class ScopeExtraction {
 public:
    // Existing public constructor (preserve original API)
    ScopeExtraction(naja::NL::SNLDesign* top0, naja::NL::SNLDesign* top1)
        : top0_(top0), top1_(top1) {}

    virtual ~ScopeExtraction() = default;

    void collectVerificationScopes();

    std::set<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>> getScopesToVerify() const {
        return designsToVerify_;
    }
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
};

#endif // SCOPE_EXTRACTION_H
