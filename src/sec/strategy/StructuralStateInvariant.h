// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <unordered_map>

#include "BoolExpr.h"
#include "common/AlignedSignals.h"
#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC {

using LocalToAbstractVarMap = std::unordered_map<size_t, size_t>;

// Rewrites each design into a shared abstract symbol space where matched SEC
// inputs and already-correlated state bits use the same variable IDs. This is
// the common base used by the structural matcher and by later proof
// strengtheners that need naming-independent comparisons.
std::pair<LocalToAbstractVarMap, LocalToAbstractVarMap> buildAbstractTransitionMaps(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedStates);

// Checks whether two BoolExpr transition fragments are structurally identical
// after remapping each design into the shared abstract symbol space.
bool areEquivalentUnderAbstractMaps(
    BoolExpr* expr0,
    BoolExpr* expr1,
    const LocalToAbstractVarMap& abstractMap0,
    const LocalToAbstractVarMap& abstractMap1);

// Matches state bits by a fixed point over their transition structure instead
// of by display names. This keeps the invariant reusable for renamed designs.
AlignedSignals inferStructurallyEquivalentStatePairs(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs);

}  // namespace KEPLER_FORMAL::SEC
