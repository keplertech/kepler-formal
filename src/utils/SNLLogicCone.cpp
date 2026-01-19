// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "SNLLogicCone.h"
#include "SNLEquipotential.h"

using namespace KEPLER_FORMAL;
using namespace naja::DNL;

void SNLLogicCone::run() {
  std::vector<naja::DNL::DNLID> currentIterationDrivers;
  std::vector<naja::DNL::DNLID> newIterationIsos;
  for (const auto& seedOutputTerm : seedOutputTerms_) {
    newIterationIsos.push_back(
        dnl_->getDNLTerminalFromID(seedOutputTerm).getIsoID());
  } 
  while (!newIterationIsos.empty()) {
    currentIterationDrivers.clear();
    for (const auto& isoID : newIterationIsos) {
      if (isoID != naja::DNL::DNLID_MAX) {
        coneIsos_.insert(isoID);
        for (auto driver :
             dnl_->getDNLIsoDB().getIsoFromIsoIDconst(isoID).getDrivers()) {
          collectedTerms_.insert(driver);
          currentIterationDrivers.push_back(driver);
        }
      }
    }
    newIterationIsos.clear();
    for (auto driver : currentIterationDrivers) {
      if (isPIs_[driver]) {
        continue;  // Skip PIs and loops(?)
      }
      collectedTerms_.insert(driver);
      const DNLInstanceFull& inst =
          dnl_->getDNLTerminalFromID(driver).getDNLInstance();
      for (DNLID termID = inst.getTermIndexes().first;
           termID <= inst.getTermIndexes().second && termID != DNLID_MAX;
           termID++) {
        const DNLTerminalFull& term = dnl_->getDNLTerminalFromID(termID);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Output && term.getIsoID() != naja::DNL::DNLID_MAX) {
          if (coneIsos_.find(term.getIsoID())  ==
              coneIsos_.end()) {
            newIterationIsos.push_back(term.getIsoID());
            collectedTerms_.insert(termID);
          }
        }
      }
    }
  }
}

std::vector<naja::NL::SNLEquipotential> SNLLogicCone::getEquipotentials()
    const {
  std::vector<naja::NL::SNLEquipotential> equipotentials;
  for (const auto& isoID : coneIsos_) {
    equipotentials.push_back(
        dnl_->getDNLTerminalFromID(
                dnl_->getDNLIsoDB().getIsoFromIsoIDconst(isoID).getDrivers()[0])
            .getEquipotential());
  }
  return equipotentials;
}