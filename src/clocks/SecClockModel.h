// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <optional>
#include <unordered_map>

#include "BoolExpr.h"
#include "common/SignalKey.h"

namespace KEPLER_FORMAL::SEC {

enum class ClockPhase {
  Pos,
  Neg,
};

struct ClockEvent {  // LCOV_EXCL_LINE
  SignalKey domain;
  ClockPhase phase = ClockPhase::Pos;
  // nullptr means the clock event is unconditionally active.
  BoolExpr* enable = nullptr;
};

struct ClockCarrierClass {  // LCOV_EXCL_LINE
  size_t varID = 0;
  SignalKey domain;
  ClockPhase phase = ClockPhase::Pos;
};

ClockPhase invertClockPhase(ClockPhase phase);
const char* clockPhaseName(ClockPhase phase);
BoolExpr* clockEventEnableOrTrue(const ClockEvent& event);
bool clockEventIsUngated(const ClockEvent& event);

std::optional<ClockEvent> classifyClockEventExpression(
    BoolExpr* expr,
    const std::unordered_map<size_t, ClockEvent>& carrierEvents);

BoolExpr* substituteBoolExprVariableExpressions(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& replacements);

// Same substitution with a memo the caller shares across roots, so logic
// shared by many expressions is rewritten once.
BoolExpr* substituteBoolExprVariableExpressions(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& replacements,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo);

}  // namespace KEPLER_FORMAL::SEC
