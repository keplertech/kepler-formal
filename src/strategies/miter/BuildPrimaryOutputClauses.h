#include <vector>
#include "BoolExpr.h"
#include "DNL.h"
#include <tbb/concurrent_vector.h>

#pragma once

namespace KEPLER_FORMAL {

class BuildPrimaryOutputClauses {
 public:
  BuildPrimaryOutputClauses() = default;
  void collect();
  void build();

  const tbb::concurrent_vector<BoolExpr*>& getPOs() const { return POs_; }
  const std::vector<naja::DNL::DNLID>& getInputs() const { return inputs_; }
  const std::vector<naja::DNL::DNLID>& getOutputs() const { return outputs_; }
  const std::map<naja::DNL::DNLID, std::pair<std::vector<naja::NL::NLID::DesignID>, std::vector<naja::NL::NLID::DesignID>>>& getInputs2InputsIDs() const { return inputs2inputsIDs_; }
  const std::map<naja::DNL::DNLID, std::pair<std::vector<naja::NL::NLID::DesignID>, std::vector<naja::NL::NLID::DesignID>>>& getOutputs2OutputsIDs() const { return outputs2outputsIDs_; }
  void setInputs(const std::vector<naja::DNL::DNLID>& inputs) { inputs_ = inputs; sortInputs(); setInputs2InputsIDs(); }
  void setOutputs(const std::vector<naja::DNL::DNLID>& outputs) { outputs_ = outputs; sortOutputs(); setOutputs2OutputsIDs(); }
  const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& getInputsMap() const { return inputsMap_; }
  const std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID>& getOutputsMap() const { return outputsMap_; }
 private:

  std::vector<naja::DNL::DNLID> collectInputs();
  void setInputs2InputsIDs();
  void sortInputs();
  std::vector<naja::DNL::DNLID> collectOutputs();
  void setOutputs2OutputsIDs();
  void sortOutputs();
  void initVarNames();
  
  tbb::concurrent_vector<BoolExpr*> POs_;
  std::vector<naja::DNL::DNLID> inputs_;
  std::vector<naja::DNL::DNLID> outputs_;
  std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID> inputsMap_;
  std::map<std::vector<NLID::DesignObjectID>, naja::DNL::DNLID> outputsMap_;
  std::map<naja::DNL::DNLID, std::pair<std::vector<naja::NL::NLID::DesignID>, std::vector<naja::NL::NLID::DesignID>>> inputs2inputsIDs_;
  std::map<naja::DNL::DNLID, std::pair<std::vector<naja::NL::NLID::DesignID>, std::vector<naja::NL::NLID::DesignID>>> outputs2outputsIDs_;
  std::vector<size_t> termDNLID2varID_; // Only for PIs
};

}  // namespace KEPLER_FORMAL