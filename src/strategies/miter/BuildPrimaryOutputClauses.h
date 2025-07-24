#include <vector>
#include "BoolExpr.h"
#include "DNL.h"

#pragma once

namespace KEPLER_FORMAL {

class BuildPrimaryOutputClauses {
 public:
  BuildPrimaryOutputClauses() = default;

  void build();

  const std::vector<std::shared_ptr<BoolExpr>>& getPOs() const { return POs_; }
  const std::vector<naja::DNL::DNLID>& getInputs() const { return inputs_; }
  const std::vector<naja::DNL::DNLID>& getOutputs() const { return outputs_; }

 private:

  std::vector<naja::DNL::DNLID> collectInputs();
  std::vector<naja::DNL::DNLID> collectOutputs();

  std::vector<std::shared_ptr<BoolExpr>> POs_;
  std::vector<naja::DNL::DNLID> inputs_;
  std::vector<naja::DNL::DNLID> outputs_;
};

}  // namespace KEPLER_FORMAL