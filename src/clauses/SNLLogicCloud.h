#include "DNL.h"
#include "SNLTruthTable.h"

namespace KEPLER_FORMAL {

class SNLLogicCloud {
 public:
  SNLLogicCloud(naja::DNL::DNLID seedOutputTerm, 
                const std::vector<naja::DNL::DNLID>& PIs)
      : seedOutputTerm_(seedOutputTerm), dnl_(*naja::DNL::get()), 
      PIs_(PIs) {}
  void compute();
  bool isInput(naja::DNL::DNLID inputTerm);
  const SNLTruthTable& getTruthTable() const { return table_; }
  const std::vector<naja::DNL::DNLID>& getInputs() const { return inputs_; }

 private:
  naja::DNL::DNLID seedOutputTerm_;
  std::vector<naja::DNL::DNLID> inputs_;
  naja::NL::SNLTruthTable table_;
  const naja::DNL::DNLFull& dnl_;
  const std::vector<naja::DNL::DNLID>& PIs_;
};

}  // namespace KEPLER_FORMAL