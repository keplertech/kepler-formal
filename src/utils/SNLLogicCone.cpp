#include "SNLLogicCone.h"
#include "SNLEquipotential.h"

using namespace KEPLER_FORMAL;
using namespace naja::DNL;

void SNLLogicCone::run() {
  std::vector<naja::DNL::DNLID> currentIterationDrivers;
  std::vector<naja::DNL::DNLID> newIterationIsos;
  newIterationIsos.push_back(
      dnl_->getDNLTerminalFromID(seedOutputTerm_).getIsoID());
  while (!newIterationIsos.empty()) {
    currentIterationDrivers.clear();
    for (const auto& isoID : newIterationIsos) {
      if (isoID != naja::DNL::DNLID_MAX) {
        coneIsos_.push_back(isoID);
        for (auto driver :
             dnl_->getDNLIsoDB().getIsoFromIsoIDconst(isoID).getDrivers()) {
          currentIterationDrivers.push_back(driver);
        }
      }
    }
    newIterationIsos.clear();
    for (auto driver : currentIterationDrivers) {
      if (std::find(PIs_.begin(), PIs_.end(), driver) != PIs_.end()) {
        continue;  // Skip PIs and loops(?)
      }
      DNLInstanceFull inst = dnl_->getDNLTerminalFromID(driver).getDNLInstance();
      for (DNLID termID = inst.getTermIndexes().first;
           termID <= inst.getTermIndexes().second && termID != DNLID_MAX; termID++) {
        const DNLTerminalFull& term = dnl_->getDNLTerminalFromID(termID);
        if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
          if (std::find(coneIsos_.begin(), coneIsos_.end(), term.getIsoID()) ==
                coneIsos_.end()) {
                newIterationIsos.push_back(term.getIsoID());
            }
        }
      }      
    }
  }
}

std::vector<naja::NL::SNLEquipotential> SNLLogicCone::getEquipotentials() const {
  std::vector<naja::NL::SNLEquipotential> equipotentials;
  for (const auto& isoID : coneIsos_) {
    equipotentials.push_back(
       dnl_->getDNLTerminalFromID(dnl_->getDNLIsoDB().getIsoFromIsoIDconst(isoID).getDrivers()[0]).getEquipotential());
  }
  return equipotentials;
}