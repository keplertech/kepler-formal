#include <vector>
#include "BoolExpr.h"

namespace naja {
namespace NL {
class SNLTruthTable;
}
}  // namespace naja

namespace KEPLER_FORMAL {

class SNLTruthTableMerger {
 public:

  SNLTruthTableMerger(const std::vector<naja::NL::SNLTruthTable>& inputsToMerge,
                   const naja::NL::SNLTruthTable& base) : inputsToMerge_(inputsToMerge), base_(base) {}

  void computeMerged() {
    mergedTable_ = mergeTruthTables(inputsToMerge_, base_);
  }

  const naja::NL::SNLTruthTable& getMergedTable() const {
    return mergedTable_;
  }

 private:
  naja::NL::SNLTruthTable mergeTruthTables(const std::vector<naja::NL::SNLTruthTable>& inputsToMerge,
                   const naja::NL::SNLTruthTable& base);

  const std::vector<naja::NL::SNLTruthTable>& inputsToMerge_;
  const naja::NL::SNLTruthTable& base_;
  naja::NL::SNLTruthTable mergedTable_;
};

}  // namespace KEPLER_FORMAL