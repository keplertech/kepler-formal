// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>
#include <unordered_map>

#include "BoolExpr.h"
#include "../../config/Config.h"

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

bool isBoolFormulaSatisfiable(
    BoolExpr* formula,
    KEPLER_FORMAL::Config::SolverType solverType);

bool boolFormulaImplies(
    BoolExpr* assumptions,
    BoolExpr* conclusion,
    KEPLER_FORMAL::Config::SolverType solverType);

std::optional<bool> boolFormulaImpliesWithConflictLimit(
    BoolExpr* assumptions,
    BoolExpr* conclusion,
    KEPLER_FORMAL::Config::SolverType solverType,
    unsigned conflictLimit);

inline BoolExpr* makeEqualityExpr(BoolExpr* lhs, BoolExpr* rhs) {
  return BoolExpr::Not(BoolExpr::Xor(lhs, rhs));
}

}  // namespace KEPLER_FORMAL::SEC
