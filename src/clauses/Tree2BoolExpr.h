#pragma once

#include <string>
#include <vector>
#include <memory>

#include "BoolExpr.h"
#include "SNLTruthTableTree.h"

namespace KEPLER_FORMAL {

/// Convert a truth-table tree directly into a BoolExpr
class Tree2BoolExpr {
    public:
    static std::shared_ptr<BoolExpr>
    convert(const SNLTruthTableTree& tree,
                const std::vector<std::string>&     varNames);
};

} // namespace KEPLER_FORMAL
