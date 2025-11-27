// Copyright 2024-2025 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <tbb/concurrent_vector.h>
#include <vector>
#include "BoolExpr.h"
#include "DNL.h"

#pragma once

namespace KEPLER_FORMAL {

class BuildPrimaryOutputClauses {
 public:
  BuildPrimaryOutputClauses() = default;
  void collect();
  void build();

  const tbb::concurrent_vector<std::shared_ptr<BoolExpr>>& getPOs() const {
    return POs_;
  }
  const std::vector<naja::DNL::DNLID>& getInputs() const { return inputs_; }
  const std::vector<naja::DNL::DNLID>& getOutputs() const { return outputs_; }
  const std::map<naja::DNL::DNLID,
                 std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>>&
  getInputs2InputsIDs() const {
    return inputs2inputsIDs_;
  }
  const std::map<naja::DNL::DNLID,
                 std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>>&
  getOutputs2OutputsIDs() const {
    return outputs2outputsIDs_;
  }
  void setInputs(const std::vector<naja::DNL::DNLID>& inputs) {
    inputs_ = inputs; /*sortInputs();*/
    setInputs2InputsIDs();
  }
  void setOutputs(const std::vector<naja::DNL::DNLID>& outputs) {
    outputs_ = outputs; /*sortOutputs();*/
    setOutputs2OutputsIDs();
  }
  const std::map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID>&
  getInputsMap() const {
    return inputsMap_;
  }
  const std::map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID>&
  getOutputsMap() const {
    return outputsMap_;
  }
  naja::DNL::DNLID getDNLIDforOutput(size_t index) const {
    return outputs_[index];
  }

 private:
  std::vector<naja::DNL::DNLID> collectInputs();
  void setInputs2InputsIDs();
  void sortInputs();
  std::vector<naja::DNL::DNLID> collectOutputs();
  void setOutputs2OutputsIDs();
  void sortOutputs();
  void initVarNames();

  tbb::concurrent_vector<std::shared_ptr<BoolExpr>> POs_;
  std::vector<naja::DNL::DNLID> inputs_;
  std::vector<naja::DNL::DNLID> outputs_;
  std::map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID> inputsMap_;
  std::map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID> outputsMap_;
  std::map<naja::DNL::DNLID,
           std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>>
      inputs2inputsIDs_;
  std::map<naja::DNL::DNLID,
           std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>>
      outputs2outputsIDs_;
  std::vector<size_t> termDNLID2varID_;  // Only for PIs
};

}  // namespace KEPLER_FORMAL