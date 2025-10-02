#include "DNL.h"
#include "SNLTruthTableTree.h"

namespace KEPLER_FORMAL {

class SNLLogicCloud {
 public:
  SNLLogicCloud(naja::DNL::DNLID seedOutputTerm, 
                const std::vector<naja::DNL::DNLID>& PIs, 
                const std::vector<naja::DNL::DNLID>& POs)
      : seedOutputTerm_(seedOutputTerm), dnl_(*naja::DNL::get()) {
        PIs_ = std::vector<bool>(naja::DNL::get()->getNBterms(), false);
        for (auto pi : PIs) {
          PIs_[pi] = true;
        }
        POs_ = std::vector<bool>(naja::DNL::get()->getNBterms(), false);
        for (auto po : POs) {
          POs_[po] = true;
        }
      }
  void compute();
  bool isInput(naja::DNL::DNLID inputTerm);
  bool isOutput(naja::DNL::DNLID inputTerm);
  const SNLTruthTableTree& getTruthTable() const { return table_; }
  const std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>& getInputs() const { return currentIterationInputs_; }
  // Get all inputs from the tree SNLTruthTableTree directly
  std::vector<naja::DNL::DNLID> getAllInputs() const {
    std::vector<naja::DNL::DNLID> allInputs;
    std::vector<const SNLTruthTableTree::Node*> stk;
    stk.push_back(table_.getRoot());
    while(!stk.empty()) {
      auto f = stk.back(); stk.pop_back();
      //printf("Node type: %d\n", (int)f->type);
      if (f->type == SNLTruthTableTree::Node::Type::P) {
        allInputs.push_back(f->termid);
      }
      else if (f->type == SNLTruthTableTree::Node::Type::Table || 
               f->type == SNLTruthTableTree::Node::Type::Input) {
        for (auto& c : f->children)
          stk.push_back(c.get());
      }
    }
    return allInputs;
  }
  

 private:
  naja::DNL::DNLID seedOutputTerm_;
  std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>> currentIterationInputs_;
  SNLTruthTableTree table_;
  const naja::DNL::DNLFull& dnl_;
  std::vector<bool> PIs_;
  std::vector<bool> POs_;
};

}  // namespace KEPLER_FORMAL