// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <unordered_map>

#include "BoolExpr.h"

namespace KEPLER_FORMAL::SEC {

BoolExpr* remapBoolExprVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, size_t>& varMap,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo);

BoolExpr* remapBoolExprVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, size_t>& varMap);

BoolExpr* substituteBoolExprVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, bool>& assignments,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo);

BoolExpr* substituteBoolExprVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, bool>& assignments);

BoolExpr* substituteBoolExprSubexpressions(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& replacements,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo);

BoolExpr* substituteBoolExprSubexpressions(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& replacements);

inline BoolExpr* makeEqualityExpr(BoolExpr* lhs, BoolExpr* rhs) {
  return BoolExpr::Not(BoolExpr::Xor(lhs, rhs));
}

}  // namespace KEPLER_FORMAL::SEC
