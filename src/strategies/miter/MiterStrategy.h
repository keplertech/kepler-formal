// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

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
  MiterStrategy(naja::NL::SNLDesign* top0, naja::NL::SNLDesign* top1, const std::string& logFileName = "", const std::string& prefix = "");

  bool run();

  void normalizeInputs(std::vector<naja::DNL::DNLID>& inputs0,
                       std::vector<naja::DNL::DNLID>& inputs1,
                        const std::map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID>& inputs0Map,
                        const std::map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID>& inputs1Map);

  void normalizeOutputs(std::vector<naja::DNL::DNLID>& outputs0,
                        std::vector<naja::DNL::DNLID>& outputs1,
                        const std::map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID>& outputs0Map,
                        const std::map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID>& outputs1Map);
  
  static std::string logFileName_;
 private:
  std::shared_ptr<BoolExpr> buildMiter(
      const tbb::concurrent_vector<std::shared_ptr<BoolExpr>>& A,
      const tbb::concurrent_vector<std::shared_ptr<BoolExpr>>& B) const;
  
  static naja::NL::SNLDesign* top0_;
  static naja::NL::SNLDesign* top1_;
  tbb::concurrent_vector<BoolExpr> POs0_;
  tbb::concurrent_vector<BoolExpr> POs1_;
  std::vector<naja::DNL::DNLID> failedPOs_;
  BoolExpr miterClause_;
  std::string prefix_;
  naja::NL::SNLDesign* topInit_ = nullptr;
  std::vector<naja::DNL::DNLFull> dnls_;
};

}  // namespace KEPLER_FORMAL