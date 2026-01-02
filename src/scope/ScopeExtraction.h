#include "DNL.h"

class ScopeExtraction {
    public:
        ScopeExtraction(naja::NL::SNLDesign* top0, naja::NL::SNLDesign* top1) : top0_(top0), top1_(top1) {}
        void collectVerificationScopes();
        std::set<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>> getScopesToVerify() const {
            return designsToVerify_;
        }
    private:
        naja::NL::SNLDesign* top0_ = nullptr;
        naja::NL::SNLDesign* top1_ = nullptr;
        std::set<std::pair<naja::NL::SNLDesign*, naja::NL::SNLDesign*>> designsToVerify_;
};