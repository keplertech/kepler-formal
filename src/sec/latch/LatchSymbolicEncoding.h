// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "latch/LatchBoundaryEncoding.h"
#include "latch/LatchSymbolicCompiler.h"

namespace KEPLER_FORMAL::SEC::LATCH {
// Map certified local symbols to SEC variables. In single-event mode raw new
// input levels are decoded from one globally aligned selector/value pair.
BoundaryEncoding encodeSymbolicMacro(const SymbolicMacro& macro, const Network& network,
    const SymbolicBits& state, const SymbolicBits& inputs, bool singleInputChange,
    const SymbolicBits& selector = {}, BoolExpr* eventValue = nullptr,
    const std::vector<size_t>& globalInputIndices = {});
}  // namespace KEPLER_FORMAL::SEC::LATCH
