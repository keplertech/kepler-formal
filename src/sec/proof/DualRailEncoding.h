// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <unordered_map>

#include "BoolExpr.h"

namespace KEPLER_FORMAL::SEC {

struct DualRailBoolExpr {
  BoolExpr* mayBeOne = nullptr;
  BoolExpr* mayBeZero = nullptr;
};

class DualRailVariableMapper {
 public:
  // LCOV_EXCL_START
  virtual ~DualRailVariableMapper() = default;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  virtual DualRailBoolExpr mapVariable(size_t symbol) = 0;
};

// Lift a Boolean expression into the paper-style two-rail X model.  This is an
// SEC encoding utility: BoolExpr remains a plain Boolean DAG, while this helper
// builds the "may be 1" and "may be 0" formulas consumed by SEC proof engines.
DualRailBoolExpr buildDualRailBoolExpr(
    BoolExpr* root,
    DualRailVariableMapper& mapper,
    std::unordered_map<BoolExpr*, DualRailBoolExpr>& memo);

}  // namespace KEPLER_FORMAL::SEC
