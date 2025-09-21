#ifndef TRUTHTABLETOBOOLEXPR_H
#define TRUTHTABLETOBOOLEXPR_H

#include "SNLTruthTableTree.h"
#include "BoolExpr.h"
#include <memory>
#include <string>
#include <vector>

namespace KEPLER_FORMAL {

/// Converts a (small) merged truth‐table tree into a BoolExpr (sum-of-products).
/// Enumerates all 2^n rows, so only use when n is modest (e.g. ≤ 10).
struct TruthTableToBoolExpr {
    /// @param tree     The merged‐tree to convert
    /// @param varNames Name for each of the external inputs (size == tree.size())
    static std::shared_ptr<BoolExpr> convert(
        const SNLTruthTableTree& tree,
        const std::vector<std::string>& varNames);
};

} // namespace KEPLER_FORMAL

#endif // TRUTHTABLETOBOOLEXPR_H
