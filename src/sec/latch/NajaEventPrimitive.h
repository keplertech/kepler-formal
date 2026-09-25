// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <map>
#include "latch/LatchEventModel.h"

namespace naja::NL { class SNLInstance; class SNLBitTerm; }
namespace KEPLER_FORMAL::SEC::LATCH {
// No cell-name heuristics: only explicit sequential expressions/truth tables.
Primitive makeNajaEventPrimitive(
    naja::NL::SNLInstance* instance, std::string path,
    const std::map<const naja::NL::SNLBitTerm*, size_t>& nets);
}  // namespace KEPLER_FORMAL::SEC::LATCH
