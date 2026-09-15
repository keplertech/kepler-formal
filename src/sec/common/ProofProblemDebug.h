// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

// Debug-only textual rendering helpers for the normalized SEC proof problem.
// Production code uses these behind explicit trace flags; readability is more
// important than compactness because the output is meant to explain a proof.
std::string formatKInductionProblemForDebug(const KInductionProblem& problem);

}  // namespace KEPLER_FORMAL::SEC
