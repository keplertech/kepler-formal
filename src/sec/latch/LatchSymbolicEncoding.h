// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "latch/LatchBoundaryEncoding.h"
#include "latch/LatchSymbolicCompiler.h"

namespace KEPLER_FORMAL::SEC::LATCH {
// Fixed-once origins are substituted in initialParameterSymbols order. Shared
// input levels must be aligned by the caller; storage origins are not implicitly
// coupled across designs. No parameter appears in the ordinary transition map.
SymbolicBits encodeSymbolicInitialState(const SymbolicMacro& macro,
                                      const SymbolicBits& parameters);
BoolExpr* encodeSymbolicInitialRelation(const SymbolicMacro& macro,
    const SymbolicBits& state, const SymbolicBits& parameters);

// Map certified local symbols to SEC variables. In single-event mode raw new
// input levels are decoded from one globally aligned selector/value pair.
BoundaryEncoding encodeSymbolicMacro(const SymbolicMacro& macro, const Network& network,
    const SymbolicBits& state, const SymbolicBits& inputs, bool singleInputChange,
    const SymbolicBits& selector = {}, BoolExpr* eventValue = nullptr,
    const std::vector<size_t>& globalInputIndices = {});
}  // namespace KEPLER_FORMAL::SEC::LATCH
