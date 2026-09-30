// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "latch/LatchSymbolicModel.h"
#include "latch/LatchSettlingCompiler.h"

namespace KEPLER_FORMAL::SEC::LATCH {
// The finite compiler's state ID encodes the complete boundary. Decode a
// remembered level, not the post-transaction observation, for reset composition.
inline BoolExpr* finiteInputHistory(const TransitionTable& table,
                                    const SymbolicBits& state, size_t net) {
  auto* result = BoolExpr::createFalse();
  for (size_t row = 0; row < table.boundaries.size(); ++row) {
    if (!table.boundaries[row].current.at(net)) continue;
    auto* match = BoolExpr::createTrue();
    for (size_t bit = 0; bit < state.size(); ++bit)
      match = BoolExpr::And(match, (row >> bit) & 1 ? state[bit] : BoolExpr::Not(state[bit]));
    result = BoolExpr::Or(result, match);
  }
  return result;
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
