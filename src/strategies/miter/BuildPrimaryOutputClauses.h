#include <vector>
#include "BoolExpr.h"
#include "DNL.h"
#include <tbb/concurrent_vector.h>

#pragma once

namespace KEPLER_FORMAL {

class BuildPrimaryOutputClauses {
 public:
  BuildPrimaryOutputClauses() = default;

  void build();

  const tbb::concurrent_vector<std::shared_ptr<BoolExpr>>& getPOs() const { return POs_; }
  const std::vector<naja::DNL::DNLID>& getInputs() const { return inputs_; }
  const std::vector<naja::DNL::DNLID>& getOutputs() const { return outputs_; }
  const std::map<naja::DNL::DNLID, std::pair<std::vector<naja::NL::NLID::DesignID>, std::vector<naja::NL::NLID::DesignID>>>& getInputs2InputsIDs() const { return inputs2inputsIDs_; }
  const std::map<naja::DNL::DNLID, std::pair<std::vector<naja::NL::NLID::DesignID>, std::vector<naja::NL::NLID::DesignID>>>& getOutputs2OutputsIDs() const { return outputs2outputsIDs_; }
  
 private:

  std::vector<naja::DNL::DNLID> collectInputs();
  void setInputs2InputsIDs();
  void sortInputs();
  std::vector<naja::DNL::DNLID> collectOutputs();
  void setOutputs2OutputsIDs();
  void sortOutputs();
  
  tbb::concurrent_vector<std::shared_ptr<BoolExpr>> POs_;
  std::vector<naja::DNL::DNLID> inputs_;
  std::vector<naja::DNL::DNLID> outputs_;
  std::map<naja::DNL::DNLID, std::pair<std::vector<naja::NL::NLID::DesignID>, std::vector<naja::NL::NLID::DesignID>>> inputs2inputsIDs_;
  std::map<naja::DNL::DNLID, std::pair<std::vector<naja::NL::NLID::DesignID>, std::vector<naja::NL::NLID::DesignID>>> outputs2outputsIDs_;
};

}  // namespace KEPLER_FORMAL