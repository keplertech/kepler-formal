// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <tbb/concurrent_vector.h>
#include <vector>
#include "BoolExpr.h"
#include "DNL.h"
#include <tbb/concurrent_unordered_map.h>
#include "SNLTruthTable.h"
#include <unordered_map>
#include <unordered_set>

#pragma once

namespace KEPLER_FORMAL {

class BuildPrimaryOutputClauses {
 public:

   struct KeyHash {
    size_t operator()(const std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>& k) const {
      size_t res = 0;
      for (const auto& name : k.first) {
        res ^= std::hash<std::string>()(name.getString()) + 0x9e3779b9 + (res << 6) + (res >> 2);
      }
      for (const auto& id : k.second) {
        res ^= std::hash<NLID::DesignObjectID>()(id) + 0x9e3779b9 + (res << 6) + (res >> 2);
      }
      return res;
    }
  };

  BuildPrimaryOutputClauses() = default;
  void collect();
  void build();

  const tbb::concurrent_vector<BoolExpr*>& getPOs() const {
    return POs_;
  }
  const std::vector<naja::DNL::DNLID>& getInputs() const { return inputs_; }
  const std::vector<naja::DNL::DNLID>& getOutputs() const { return outputs_; }
  const std::unordered_map<naja::DNL::DNLID,
                 std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>>&
  getInputs2InputsIDs() const {
    return inputs2inputsIDs_;
  }
  const std::unordered_map<naja::DNL::DNLID,
                 std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>>&
  getOutputs2OutputsIDs() const {
    return outputs2outputsIDs_;
  }
  void setInputs(const std::vector<naja::DNL::DNLID>& inputs) {
    inputs_ = std::move(inputs); /*sortInputs();*/
    setInputs2InputsIDs();
  }
  void setOutputs(const std::vector<naja::DNL::DNLID>& outputs) {
    outputs_ = std::move(outputs); /*sortOutputs();*/
    setOutputs2OutputsIDs();
  }
  const std::unordered_map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID, KeyHash>&
  getInputsMap() const {
    return inputsMap_;
  }
  const std::unordered_map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID, KeyHash>&
  getOutputsMap() const {
    return outputsMap_;
  }
  naja::DNL::DNLID getDNLIDforOutput(size_t index) const {
    return outputs_[index];
  }
  naja::DNL::DNLID getDNLIDforInput(size_t index) const {
    return inputs_[index];
  }

  const std::vector<size_t>& getTermDNLID2VarID() const {
    return termDNLID2varID_;
  }
  void setLastCommonID(size_t id) { lastCommonID = id; }

 private:

  //const naja::NL::SNLTruthTable& getTruthTable(naja::NL::SNLDesign* design, size_t orderID);
  
  std::vector<naja::DNL::DNLID> collectInputs();
  void setInputs2InputsIDs();
  // void sortInputs();
  std::vector<naja::DNL::DNLID> collectOutputs();
  void setOutputs2OutputsIDs();
  // void sortOutputs();
  void initVarNames();
  
  tbb::concurrent_vector<BoolExpr*> POs_;
  std::vector<naja::DNL::DNLID> inputs_;
  std::vector<bool> IsPIs_;
  std::vector<naja::DNL::DNLID> outputs_;
  std::vector<bool> IsPOs_;

  std::unordered_map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID, KeyHash> inputsMap_;
  std::unordered_map<std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>, naja::DNL::DNLID, KeyHash> outputsMap_;
  std::unordered_map<naja::DNL::DNLID,
           std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>>
      inputs2inputsIDs_;
  std::unordered_map<naja::DNL::DNLID,
           std::pair<std::vector<NLName>, std::vector<NLID::DesignObjectID>>>
      outputs2outputsIDs_;
  std::vector<size_t> termDNLID2varID_;  // Only for PIs
  size_t lastCommonID = 1;

  struct hash {
    size_t operator()(const std::pair<unsigned int, unsigned long>& p) const noexcept {
      // A simple, solid 64-bit mix
      uint64_t h = (uint64_t(p.first) << 32) ^ uint64_t(p.second);
      h ^= (h >> 33);
      h *= 0xff51afd7ed558ccdULL;
      h ^= (h >> 33);
      h *= 0xc4ceb9fe1a85ec53ULL;
      h ^= (h >> 33);
      return size_t(h);
    }
  };

  std::unordered_map<std::pair<naja::NL::NLID::DesignObjectID, size_t>, naja::NL::SNLTruthTable, hash> ttCache_;

  struct CachedModel {
    CachedModel() : analyzedPIs(false), analyzedPOs(false) {}
    std::unordered_set<naja::NL::SNLBitTerm*> PIs;
    bool analyzedPIs = false;
    std::unordered_set<naja::NL::SNLBitTerm*> POs;
    bool analyzedPOs = false;
  };

  std::unordered_map<const naja::NL::SNLDesign*, CachedModel> modelCache_;
};

}  // namespace KEPLER_FORMAL