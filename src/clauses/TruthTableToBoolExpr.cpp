#include "TruthTableToBoolExpr.h"
#include <cstdio>
#include <stdexcept>

using namespace KEPLER_FORMAL;
using namespace naja::NL;

BoolExpr* TruthTableToBoolExpr::convert(
    const SNLTruthTableTree& tree,
    const std::vector<std::string>& varNames)
{
    size_t n = tree.size();
    if (varNames.size() < n) {
        std::fprintf(stderr,
                     "VarNames vector size (%zu) is less than tree arity (%zu).\n",
                     varNames.size(), n);
        throw std::invalid_argument(
            "varNames vector size must match truth-table tree arity");
    }

    // handle the zero-input degenerate case
    if (n == 0) {
        bool out0 = tree.eval(std::vector<bool>{});
        auto v = BoolExpr::Var(out0 ? "TRUE" : "FALSE");
        return out0 ? v : BoolExpr::Not(v);
    }

    uint64_t rows = 1ULL << n;
    std::vector<bool> inputPattern(n);

    std::vector<BoolExpr*> cubes;
    cubes.reserve(rows / 2); // guess half the rows may be 1

    // build one cube for each row where tree.eval(...) == true
    for (uint64_t row = 0; row < rows; ++row) {
        // unpack row bits into inputPattern
        for (size_t b = 0; b < n; ++b) {
            inputPattern[b] = ((row >> b) & 1) != 0;
        }

        if (!tree.eval(inputPattern))
            continue;

        // form the AND-cube for this row
        BoolExpr* cube;
        for (size_t b = 0; b < n; ++b) {
            auto var = BoolExpr::Var(varNames[b]);
            auto lit = inputPattern[b] ? var : BoolExpr::Not(var);
            cube = cube ? BoolExpr::And(cube, lit) : lit;
        }
        cubes.push_back(cube);
    }

    // if no minterms, return constant FALSE
    if (cubes.empty()) {
        auto v = BoolExpr::Var("FALSE");
        return BoolExpr::Not(v);
    }

    // OR together all the cubes
    auto expr = cubes.front();
    for (size_t i = 1; i < cubes.size(); ++i) {
        expr = BoolExpr::Or(expr, cubes[i]);
    }

    return expr;
}
