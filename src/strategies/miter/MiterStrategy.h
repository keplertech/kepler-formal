#include <vector>
#include "DNL.h"

namespace KEPLER_FORMAL {
class BoolExpr;

class MiterStrategy {
 public:
  MiterStrategy() = default;

  void build();

  

 private:
 
  std::vector<naja::DNL::DNLID> collectInputs();
  std::vector<naja::DNL::DNLID> collectOutputs();

  std::vector<BoolExpr> POs_;
  std::vector<naja::DNL::DNLID> inputs_;
  std::vector<naja::DNL::DNLID> outputs_;
};

}  // namespace KEPLER_FORMAL