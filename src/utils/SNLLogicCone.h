// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include "DNL.h"
#include <tbb/concurrent_unordered_set.h>

namespace naja {
namespace NL {
class SNLEquipotential;
}
}  // namespace naja

namespace KEPLER_FORMAL {

class SNLLogicCone {
 public:
  SNLLogicCone(naja::DNL::DNLID seedOutputTerm,
               std::vector<naja::DNL::DNLID> pis)
      : PIs_(pis) {
    //naja::DNL::destroy();
    seedOutputTerms_.push_back(seedOutputTerm);
    dnl_ = naja::DNL::get();
  }
  SNLLogicCone(std::vector<naja::DNL::DNLID> seedOutputTerms,
               std::vector<naja::DNL::DNLID> pis)
      : PIs_(pis) {
    //naja::DNL::destroy();
    for (const auto& term : seedOutputTerms) {
      seedOutputTerms_.push_back(term);
    }
    dnl_ = naja::DNL::get();
  }
  SNLLogicCone(naja::DNL::DNLID seedOutputTerm,
               std::vector<naja::DNL::DNLID> pis,
               naja::DNL::DNLFull* dnl)
      : PIs_(pis) {
    //naja::DNL::destroy();
    seedOutputTerms_.push_back(seedOutputTerm);
    dnl_ = dnl;
  }
  void run();
  std::vector<naja::NL::SNLEquipotential> getEquipotentials() const;

  const tbb::concurrent_unordered_set<naja::DNL::DNLID>& getConeIsoIDs() const {
    return coneIsos_;
  }

  void initConeIsos(const tbb::concurrent_unordered_set<naja::DNL::DNLID>& isoIDs) {
    coneIsos_ = isoIDs;
  }

 private:
  std::vector<naja::DNL::DNLID> seedOutputTerms_;
  tbb::concurrent_unordered_set<naja::DNL::DNLID> coneIsos_;
  std::vector<naja::DNL::DNLID> PIs_;
  naja::DNL::DNLFull* dnl_;
};

}  // namespace KEPLER_FORMAL
