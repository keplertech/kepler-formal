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

  SNLTruthTableMerger(const std::vector<const naja::NL::SNLTruthTable*>& inputsToMerge,
                   const naja::NL::SNLTruthTable& base) : inputsToMerge_(inputsToMerge), base_(base) {}

  naja::NL::SNLTruthTable mergeTruthTables(const std::vector<const naja::NL::SNLTruthTable*>&,
                   const naja::NL::SNLTruthTable& base);

 private:
  const std::vector<const naja::NL::SNLTruthTable*>& inputsToMerge_;
  const naja::NL::SNLTruthTable& base_;
  BoolExpr boolExpr_;
};

}  // namespace KEPLER_FORMAL