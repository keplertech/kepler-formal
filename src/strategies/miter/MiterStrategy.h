// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <vector>
#include "BoolExpr.h"
#include "DNL.h"
#include <tbb/concurrent_vector.h>
#include "BuildPrimaryOutputClauses.h"

#pragma once

namespace naja {
namespace NL {
class SNLDesign;
}
}  // namespace naja

namespace KEPLER_FORMAL {

class MiterStrategy {
 public:
  struct CompactSnapshot {
    std::vector<BuildPrimaryOutputClauses::PathKey> inputs;
    std::vector<BuildPrimaryOutputClauses::PathKey> outputs;
    tbb::concurrent_vector<BoolExpr*> POs;
  };

  MiterStrategy(naja::NL::SNLDesign* top0, naja::NL::SNLDesign* top1, const std::string& logFileName = "", const std::string& prefix = "");

  void init(bool enableLogging = true);

  bool run(bool compact = false);
  bool runCompactSnapshots(const CompactSnapshot& snapshot0,
                           const CompactSnapshot& snapshot1);

  void setCnfDump(bool enabled, const std::string& path = "");
  void setPoCnfDump(bool enabled, const std::string& path = "");

  size_t normalizeInputs(std::vector<naja::DNL::DNLID>& inputs0,
                       std::vector<naja::DNL::DNLID>& inputs1,
                        const std::unordered_map<KEPLER_FORMAL::BuildPrimaryOutputClauses::PathKey, naja::DNL::DNLID, KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash>& inputs0Map,
                        const std::unordered_map<KEPLER_FORMAL::BuildPrimaryOutputClauses::PathKey, naja::DNL::DNLID, KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash>& inputs1Map);

  void normalizeOutputs(std::vector<naja::DNL::DNLID>& outputs0,
                        std::vector<naja::DNL::DNLID>& outputs1,
                        const std::unordered_map<KEPLER_FORMAL::BuildPrimaryOutputClauses::PathKey, naja::DNL::DNLID, KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash>& outputs0Map,
                        const std::unordered_map<KEPLER_FORMAL::BuildPrimaryOutputClauses::PathKey, naja::DNL::DNLID, KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash>& outputs1Map);
  
  static std::string logFileName_;
  const std::vector<naja::DNL::DNLID>& getPIs0() const { return PIs0_; }
  const std::vector<naja::DNL::DNLID>& getPIs1() const { return PIs1_; }
 private:
  bool runCompactPOs(const tbb::concurrent_vector<BoolExpr*>& POs0,
                     const tbb::concurrent_vector<BoolExpr*>& POs1);
  BoolExpr* buildMiter(
      const tbb::concurrent_vector<BoolExpr*>& A,
      const tbb::concurrent_vector<BoolExpr*>& B) const;
  BuildPrimaryOutputClauses builder0_;
  BuildPrimaryOutputClauses builder1_;
  static naja::NL::SNLDesign* top0_;
  static naja::NL::SNLDesign* top1_;
  std::vector<naja::DNL::DNLID> PIs0_;
  std::vector<naja::DNL::DNLID> PIs1_;
  tbb::concurrent_vector<BoolExpr> POs0_;
  tbb::concurrent_vector<BoolExpr> POs1_;
  std::vector<naja::DNL::DNLID> failedPOs_;
  BoolExpr miterClause_;
  std::string prefix_;
  naja::NL::SNLDesign* topInit_ = nullptr;
  std::vector<naja::DNL::DNLFull> dnls_;
  size_t lastCommonVarID_ = 1;
  bool dumpCnf_ = false;
  std::string dumpCnfPath_;
  bool dumpPoCnf_ = false;
  std::string dumpPoCnfPath_;
};

}  // namespace KEPLER_FORMAL
