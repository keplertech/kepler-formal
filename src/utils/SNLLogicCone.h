#pragma once

#include "DNL.h"

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
      : seedOutputTerm_(seedOutputTerm), PIs_(pis) {
    naja::DNL::destroy();
    dnl_ = naja::DNL::get();
  }
  SNLLogicCone(naja::DNL::DNLID seedOutputTerm,
               std::vector<naja::DNL::DNLID> pis,
               naja::DNL::DNLFull* dnl)
      : seedOutputTerm_(seedOutputTerm), PIs_(pis) {
    naja::DNL::destroy();
    dnl_ = dnl;
  }
  void run();
  std::vector<naja::NL::SNLEquipotential> getEquipotentials() const;

 private:
  naja::DNL::DNLID seedOutputTerm_;
  std::vector<naja::DNL::DNLID> coneIsos_;
  std::vector<naja::DNL::DNLID> PIs_;
  naja::DNL::DNLFull* dnl_;
};

}  // namespace KEPLER_FORMAL
