// Tree2BoolExpr.cpp

#include "Tree2BoolExpr.h"
#include "SNLTruthTableTree.h"     // for SNLTruthTableTree::Node
#include "SNLTruthTable.h"
#include "BoolExpr.h"

#include <unordered_map>
#include <stdexcept>
#include <cstdint>
#include <vector>
#include <utility>   // for std::pair
#include <bitset>    // needed for bit‐manipulation

using namespace naja::NL;
using namespace KEPLER_FORMAL;

size_t toSizeT(const std::string& s) {
    if (s.empty()) {
        throw std::invalid_argument("toSizeT: input string is empty");
    }

    size_t result = 0;
    const size_t max = std::numeric_limits<size_t>::max();

    for (char c : s) {
        if (!std::isdigit(static_cast<unsigned char>(c))) {
            throw std::invalid_argument(
                "toSizeT: invalid character '" + std::string(1, c) + "' in input"
            );
        }
        size_t digit = static_cast<size_t>(c - '0');

        // Check for overflow: result * 10 + digit > max
        if (result > (max - digit) / 10) {
            throw std::out_of_range("toSizeT: value out of range for size_t");
        }

        result = result * 10 + digit;
    }

    return result;
}

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

std::shared_ptr<BoolExpr>
Tree2BoolExpr::convert(
    const SNLTruthTableTree&        tree,
    const std::vector<std::string>& varNames)
{
    if (varNames.size() != tree.size()) {
        throw std::invalid_argument{
            "Tree2BoolExpr: varNames.size() != tree.size()"
        };
    }

    const auto* root = tree.getRoot();
    if (!root) 
        return nullptr;

    // 1) find maxID
    size_t maxID = 0;
    {
        std::vector<const SNLTruthTableTree::Node*> dfs{root};
        while (!dfs.empty()) {
            auto* n = dfs.back(); dfs.pop_back();
            maxID = std::max(maxID, n->nodeID);
            if (n->type == SNLTruthTableTree::Node::Type::Table ||
                n->type == SNLTruthTableTree::Node::Type::P)
            {
                for (auto& c : n->children)
                    dfs.push_back(c.get());
            }
        }
    }

    // 2) memo table
    std::vector<std::shared_ptr<BoolExpr>> memo(maxID + 1);

    // 3) post-order build
    using Frame = std::pair<const SNLTruthTableTree::Node*, bool>;
    std::vector<Frame> stack;
    stack.reserve(maxID + 1);
    stack.emplace_back(root, false);

    while (!stack.empty()) {
        auto [node, visited] = stack.back();
        stack.pop_back();
        size_t id = node->nodeID;

        if (!visited) {
            if (memo[id]) 
                continue;

            if (node->type == SNLTruthTableTree::Node::Type::Table ||
                node->type == SNLTruthTableTree::Node::Type::P)
            {
                stack.emplace_back(node, true);
                for (auto& c : node->children)
                    stack.emplace_back(c.get(), false);
            }
            else {
                // input leaf
                size_t vidx = node->inputIndex;
                size_t vid  = toSizeT(varNames[vidx]);
                memo[id]    = BoolExpr::Var(vid);
            }
        }
        else {
            // post-visit for Table / P
            auto tbl  = node->getTruthTable();
            uint32_t k = tbl.size();
            uint64_t rows = uint64_t{1} << k;

            if (tbl.all0()) {
                memo[id] = BoolExpr::createFalse();
            }
            else if (tbl.all1()) {
                memo[id] = BoolExpr::createTrue();
            }
            else {
                // gather children
                std::vector<std::shared_ptr<BoolExpr>> childF(k);
                for (uint32_t i = 0; i < k; ++i) {
                    auto cid = node->children[i]->nodeID;
                    childF[i] = memo[cid];
                }

                // find which inputs actually matter
                std::vector<bool> relevant(k, false);
                for (uint32_t j = 0; j < k; ++j) {
                    for (uint64_t m = 0; m < rows; ++m) {
                        bool b0 = tbl.bits().bit(m);
                        bool b1 = tbl.bits().bit(m ^ (uint64_t{1} << j));
                        if (b0 != b1) { relevant[j] = true; break; }
                    }
                }

                // collect the indices of relevant vars
                std::vector<uint32_t> relIdx;
                for (uint32_t j = 0; j < k; ++j)
                    if (relevant[j]) relIdx.push_back(j);

                // if nothing matters, fall back to constant‐false
                if (relIdx.empty()) {
                    memo[id] = BoolExpr::createFalse();
                }
                else {
                    // build the DNF terms
                    std::vector<std::shared_ptr<BoolExpr>> terms;
                    terms.reserve(rows);

                    for (uint64_t m = 0; m < rows; ++m) {
                        if (!tbl.bits().bit(m)) continue;

                        std::shared_ptr<BoolExpr> term;
                        bool firstLit = true;

                        for (auto j : relIdx) {
                            bool bit1 = ((m >> j) & 1) != 0;
                            auto lit  = bit1 
                                      ? childF[j] 
                                      : BoolExpr::Not(childF[j]);

                            if (firstLit) {
                                term = lit;
                                firstLit = false;
                            } else {
                                term = BoolExpr::And(term, lit);
                            }
                        }

                        // only push if we actually got a literal
                        if (term)
                            terms.push_back(std::move(term));
                    }

                    // guard against an empty terms list
                    if (terms.empty()) {
                        memo[id] = BoolExpr::createFalse();
                    }
                    else {
                        // fold into OR
                        auto expr = terms[0];
                        for (size_t t = 1; t < terms.size(); ++t)
                            expr = BoolExpr::Or(expr, terms[t]);
                        memo[id] = expr;
                    }
                }
            }
        }
    }

    // 4) return root
    return memo[root->nodeID];
}

