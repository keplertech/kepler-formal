// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "phase/PeriodicSignalAnalysis.h"

namespace KEPLER_FORMAL::SEC::PHASE {

enum class TransparencyKind : uint8_t {
  Closed,
  Transparent,
  Conditional,
};

struct PhaseTransparency {
  TransparencyKind kind = TransparencyKind::Conditional;
  // For Conditional phases, this is the transparency predicate after all
  // proven generator values for the phase have been substituted.  Variable
  // indexes retain their meaning from the NormalizedTransitionSystem.
  ExprID residualExpression = InvalidExpr;
};

struct LatchPhaseProfile {
  size_t latchIndex = 0;
  std::vector<PhaseTransparency> phases;
};

struct LatchPhaseProfileOptions {
  // 0 means the TBB default, 1 forces serial execution.
  size_t maxConcurrency = 0;
};

struct LatchPhaseProfileResult {
  // Owns every residual expression referenced by profiles.
  ExpressionArena residualExpressions;
  std::vector<LatchPhaseProfile> profiles;
  bool complete = true;
  std::string diagnostic;
};

// Interpret the paper's periodic state-signal generators for Kepler latches.
// This is intentionally separate from generator discovery: latches with
// `clock & gate1` and `clock & gate2` may share carrier phases while retaining
// distinct residual local gates.
LatchPhaseProfileResult buildLatchPhaseProfiles(
    const NormalizedTransitionSystem& system,
    const PeriodicSignalAnalysisResult& periodic,
    const LatchPhaseProfileOptions& options = {});

const char* toString(TransparencyKind kind);

}  // namespace KEPLER_FORMAL::SEC::PHASE
