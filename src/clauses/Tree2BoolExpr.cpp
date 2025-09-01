#include "Tree2BoolExpr.h"
#include "SNLTruthTableTree.h"   // for ITTNode, InputNode, TableNode
#include "SNLTruthTable.h"       // for SNLTruthTable::bits()
#include "BoolExpr.h"

#include <unordered_map>
#include <stdexcept>
#include <cstdint>
#include <vector>

using namespace naja::NL;
using namespace KEPLER_FORMAL;

// Fold a list of literals into a single AND
static std::shared_ptr<BoolExpr>
mkAnd(const std::vector<std::shared_ptr<BoolExpr>>& lits)
{
    if (lits.empty()) 
        return BoolExpr::createTrue();
    auto cur = lits[0];
    for (size_t i = 1; i < lits.size(); ++i)
        cur = BoolExpr::And(cur, lits[i]);
    return cur;
}

// Fold a list of terms into a single OR
static std::shared_ptr<BoolExpr>
mkOr(const std::vector<std::shared_ptr<BoolExpr>>& terms)
{
    if (terms.empty()) 
        return BoolExpr::createFalse();
    auto cur = terms[0];
    for (size_t i = 1; i < terms.size(); ++i)
        cur = BoolExpr::Or(cur, terms[i]);
    return cur;
}

namespace {
using NodeKey = const ITTNode*;

/// Recursively build a BoolExpr for the subtree at `node`
std::shared_ptr<BoolExpr>
buildExpr(NodeKey node,
          const std::vector<std::string>&                    varNames,
          std::unordered_map<NodeKey,std::shared_ptr<BoolExpr>>& memo)
{
    // 1) Memo‐reuse
    auto it = memo.find(node);
    if (it != memo.end())
        return it->second;

    std::shared_ptr<BoolExpr> result;

    // 2) InputNode -> variable
    if (auto in = dynamic_cast<const InputNode*>(node)) {
        //printf("%p: InputNode %zu\n", (void*)node, in->inputIndex);
        //printf("varNames.size() = %zu\n", varNames.size());
        //printf("InputNode %zu -> var %s\n", in->inputIndex, varNames[in->inputIndex].c_str());
        result = BoolExpr::Var(varNames[in->inputIndex]);
    }
    // 3) TableNode -> build local DNF
    else if (auto tn = dynamic_cast<const TableNode*>(node)) {
        const SNLTruthTable& tbl = tn->table;
        uint32_t k    = tbl.size();
        uint64_t rows = uint64_t{1} << k;

        // a) constant‐table shortcuts
        if (tbl.all0()) {
            result = BoolExpr::createFalse();
        }
        else if (tbl.all1()) {
            result = BoolExpr::createTrue();
        }
        else {
            // b) recurse on children
            std::vector<std::shared_ptr<BoolExpr>> childF(k);
            for (uint32_t i = 0; i < k; ++i)
                childF[i] = buildExpr(tn->children[i].get(), varNames, memo);

            // c) collect terms for each true‐row
            std::vector<std::shared_ptr<BoolExpr>> terms;
            terms.reserve(rows);
            for (uint64_t m = 0; m < rows; ++m) {
                if (!tbl.bits().bit(m))  // skip false outputs
                    continue;

                // build the conjunction of literals for row m
                std::vector<std::shared_ptr<BoolExpr>> lits;
                lits.reserve(k);
                for (uint32_t j = 0; j < k; ++j) {
                    bool bitIs1 = ((m >> j) & 1) != 0;
                    if (bitIs1)
                        lits.push_back(childF[j]);
                    else
                        lits.push_back(BoolExpr::Not(childF[j]));
                }
                terms.push_back(mkAnd(lits));
            }
            result = mkOr(terms);
        }
    }
    else {
        throw std::logic_error(
          "Tree2BoolExpr: encountered unknown ITTNode subtype"
        );
    }

    // 4) memoize & return
    memo[node] = result;
    return result;
}
} // anonymous namespace

std::shared_ptr<BoolExpr>
KEPLER_FORMAL::Tree2BoolExpr::convert(
    const SNLTruthTableTree&        tree,
    const std::vector<std::string>& varNames)
{
    if (varNames.size() != tree.size())
        throw std::invalid_argument(
          "Tree2BoolExpr: varNames.size() != tree.size()"
        );

    std::unordered_map<NodeKey,std::shared_ptr<BoolExpr>> memo;
    return buildExpr(tree.getRoot(), varNames, memo);
}
