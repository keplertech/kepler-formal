#include <vector>
#include "BoolExpr.h"
#include "DNL.h"
#include <tbb/concurrent_vector.h>

#pragma once

namespace naja {
namespace NL {
class SNLDesign;
}
}  // namespace naja

namespace KEPLER_FORMAL {

class MiterStrategy {
 public:
  MiterStrategy(naja::NL::SNLDesign* top0, naja::NL::SNLDesign* top1, const std::string& prefix = "")
      : top0_(top0), top1_(top1), prefix_(prefix) {}

  bool run();

  void normalizeInputs(std::vector<naja::DNL::DNLID>& inputs0,
                       std::vector<naja::DNL::DNLID>& inputs1,
                        const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& inputs0Map,
                        const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& inputs1Map);

  void normalizeOutputs(std::vector<naja::DNL::DNLID>& outputs0,
                        std::vector<naja::DNL::DNLID>& outputs1,
                        const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& outputs0Map,
                        const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& outputs1Map);
  

 private:
  BoolExpr* buildMiter(
      const tbb::concurrent_vector<BoolExpr*>& A,
      const tbb::concurrent_vector<BoolExpr*>& B) const;

  naja::NL::SNLDesign* top0_ = nullptr;
  naja::NL::SNLDesign* top1_ = nullptr;
  tbb::concurrent_vector<BoolExpr> POs0_;
  tbb::concurrent_vector<BoolExpr> POs1_;
  std::vector<naja::DNL::DNLID> failedPOs_;
  BoolExpr miterClause_;
  std::string prefix_;
  naja::NL::SNLDesign* topInit_ = nullptr;
};

}  // namespace KEPLER_FORMAL