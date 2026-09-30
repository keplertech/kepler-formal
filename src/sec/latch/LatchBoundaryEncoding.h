// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "BoolExpr.h"
#include "latch/LatchSettlingCompiler.h"

namespace KEPLER_FORMAL::SEC::LATCH {
struct BoundaryEncoding {
  std::vector<BoolExpr*> nextState;
  // Value at the completed observation of the incoming external transaction.
  std::vector<BoolExpr*> observedNets;
};

size_t boundaryEncodingBits(size_t states);

// The finite table is the exact certified macro relation, encoded injectively.
// State IDs retain the entire quiescent state (table is its reconstruction).
// In single-change mode selectors name global external pins; other selector
// codes mean a stutter for this component. No input/scheduler path is assumed
// away. The shared decoder generates precisely the declared environment.
BoundaryEncoding encodeBoundaryTable(
    const TransitionTable& table, const Network& network,
    const std::vector<BoolExpr*>& state,
    const std::vector<BoolExpr*>& inputs,
    const std::vector<BoolExpr*>& selector = {},
    BoolExpr* eventValue = nullptr,
    const std::vector<size_t>& globalInputIndices = {});
}  // namespace KEPLER_FORMAL::SEC::LATCH
