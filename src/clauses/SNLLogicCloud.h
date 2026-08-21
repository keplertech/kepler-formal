// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "DNL.h"
#include "SNLTruthTableTree.h"
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace KEPLER_FORMAL {

class SNLLogicCloud {
 public:
  enum class SkipReason {
    None,
    MultiDriver,
    NoDriver,
    LogicalLoop,
    OpaqueInternal,
  };

  SNLLogicCloud(naja::DNL::DNLID seedOutputTerm,
                const std::vector<bool>& PIs,
                const std::vector<bool>& POs,
                bool stopAtOpaqueInternalOutputs = false)
      : seedOutputTerm_(seedOutputTerm), dnl_(*naja::DNL::get()),
        PIs_(PIs), POs_(POs),
        stopAtOpaqueInternalOutputs_(stopAtOpaqueInternalOutputs) {
  }
  void compute();
  static void flushSkippedPOReports();
  bool isInput(naja::DNL::DNLID inputTerm);
  bool isOutput(naja::DNL::DNLID inputTerm);
  SkipReason getSkipReason() const { return skipReason_; }
  const std::string& getSkipReasonText() const { return skipReasonText_; }
  naja::DNL::DNLID getOpaqueInternalTerm() const {
    return opaqueInternalTerm_;
  }
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
  using TermIDVector =
      std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>;

  naja::DNL::DNLID getIsoIDCached(
      naja::DNL::DNLID termID,
      const std::shared_ptr<const std::vector<naja::DNL::DNLID>>&
          termIsoIDs) const;
  std::string formatTermName(naja::DNL::DNLID termID) const;
  void appendTermList(std::ostream& out,
                      const TermIDVector& termIDs,
                      size_t limit = 24) const;
  void appendInstNonOutputs(std::ostream& out,
                            naja::DNL::DNLID instID,
                            size_t limit = 24) const;
  void appendMergeListDetailed(std::ostream& out,
                               size_t limit = 24) const;
  std::string buildIterationSnapshot(size_t iter) const;

  naja::DNL::DNLID resolveInstanceInputTerm(
      const naja::DNL::DNLInstanceFull& inst,
      size_t flatTermID,
      naja::DNL::DNLID driver,
      const char* role) const;
  size_t getRelevantInstanceInputCount(naja::DNL::DNLID driver) const;
  void appendRelevantInstanceInputs(
      naja::DNL::DNLID driver,
      TermIDVector& relevantTerms) const;
  TermIDVector collectRelevantInstanceInputs(naja::DNL::DNLID driver) const;
  void throwIfTruthTableArityMismatch(naja::DNL::DNLID driver) const;
  void throwIfFrontierMismatch(
      size_t iter,
      const std::vector<std::string>& frontierHistory) const;
  void throwIfNextFrontierMismatch(
      size_t iter,
      const std::vector<std::string>& frontierHistory) const;
  naja::DNL::DNLID resolveTransparentLoopTarget(
      naja::DNL::DNLID termID,
      const std::shared_ptr<const std::vector<naja::DNL::DNLID>>&
          termIsoIDs) const;
  bool rejectOpaqueInternalOutput(naja::DNL::DNLID termID);

  naja::DNL::DNLID seedOutputTerm_;
  TermIDVector currentIterationInputs_;
  SNLTruthTableTree table_;
  const naja::DNL::DNLFull& dnl_;
  const std::vector<bool>& PIs_;
  const std::vector<bool>& POs_;
  bool stopAtOpaqueInternalOutputs_ = false;
  SkipReason skipReason_ = SkipReason::None;
  std::string skipReasonText_;
  naja::DNL::DNLID opaqueInternalTerm_ = naja::DNL::DNLID_MAX;
};

}  // namespace KEPLER_FORMAL
