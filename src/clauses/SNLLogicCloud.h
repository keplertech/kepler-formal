// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "DNL.h"
#include "SNLTruthTableTree.h"

namespace KEPLER_FORMAL {

class SNLLogicCloud {
 public:
  enum class SkipReason {
    None,
    MultiDriver,
    NoDriver,
    LogicalLoop,
  };

  SNLLogicCloud(naja::DNL::DNLID seedOutputTerm,
                const std::vector<bool>& PIs,
                const std::vector<bool>& POs)
      : seedOutputTerm_(seedOutputTerm), dnl_(*naja::DNL::get()),
        PIs_(PIs), POs_(POs) {
  }
  void compute();
  static void resetAnalysisCaches();
  static void flushSkippedPOReports();
  bool isInput(naja::DNL::DNLID inputTerm);
  bool isOutput(naja::DNL::DNLID inputTerm);
  SkipReason getSkipReason() const { return skipReason_; }
  const std::string& getSkipReasonText() const { return skipReasonText_; }
  SNLTruthTableTree& getTruthTable() { return table_; }
  const std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>& getInputs() const {
    return currentIterationInputs_;
  }
  // Get all inputs from the tree SNLTruthTableTree directly
  std::vector<naja::DNL::DNLID> getAllInputs() const {
    std::vector<naja::DNL::DNLID> allInputs;
    std::vector<std::shared_ptr<SNLTruthTableTree::Node>> stk;
    stk.push_back(table_.getRoot());
    while (!stk.empty()) {
      auto f = stk.back();
      stk.pop_back();
      // printf("Node type: %d\n", (int)f->type);
      if (f->type == SNLTruthTableTree::Node::Type::P) {
        allInputs.push_back(f->data.termid);
      } else if (f->type == SNLTruthTableTree::Node::Type::Table ||
                 f->type == SNLTruthTableTree::Node::Type::Input) {
        for (auto& c : f->childrenIds)
          stk.push_back(table_.nodeFromId(c));
      }
    }
    return allInputs;
  }
  void destroy() { table_.destroy(); }

 private:
  naja::DNL::DNLID seedOutputTerm_;
  std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>> currentIterationInputs_;
  SNLTruthTableTree table_;
  const naja::DNL::DNLFull& dnl_;
  const std::vector<bool>& PIs_;
  const std::vector<bool>& POs_;
  SkipReason skipReason_ = SkipReason::None;
  std::string skipReasonText_;
};

}  // namespace KEPLER_FORMAL
