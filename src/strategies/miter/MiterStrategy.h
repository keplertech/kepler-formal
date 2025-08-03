#include <vector>
#include "BoolExpr.h"
#include "DNL.h"

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

 private:
  std::shared_ptr<BoolExpr> buildMiter(
      const std::vector<std::shared_ptr<BoolExpr>>& A,
      const std::vector<std::shared_ptr<BoolExpr>>& B) const;

  naja::NL::SNLDesign* top0_ = nullptr;
  naja::NL::SNLDesign* top1_ = nullptr;
  std::vector<BoolExpr> POs0_;
  std::vector<BoolExpr> POs1_;
  std::vector<naja::DNL::DNLID> failedPOs_;
  BoolExpr miterClause_;
  std::string prefix_;
  naja::NL::SNLDesign* topInit_ = nullptr;
};

}  // namespace KEPLER_FORMAL