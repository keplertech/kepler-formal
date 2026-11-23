// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <filesystem>
#include <set>
#include <string>

namespace naja::NL {
class NLLibrary;
class NLUniverse;
}  // namespace naja::NL

namespace KEPLER_FORMAL {

// Integrated clock-gating cells (ICGs).
//
// Liberty describes an ICG with a `statetable` for its enable latch and a
// `state_function` for its gated clock. Naja keeps no function for such a
// cell, so every flop clocked through one is opaque to SEC. The functions
// below give the cell the structure SEC already reduces: a transparent latch
// on the enable and a gate on the clock. The latch abstraction turns the
// latch into a wire and the clock model turns the gate into a flop enable
// (docs/sec-clock-handling.md).

// Reads the ICG cells of `libertyPath`. For each one that `primitives` holds,
// whose statetable is a single latch open while its clock is low and whose
// outputs are gated copies of that clock, adds to `primitives` the latch and
// the gate the cell decomposes into. Returns the number of cells modeled. A
// cell that does not fit stays as Naja loaded it, opaque: a negative-edge
// clock gate, an inverted gated clock, an edge-triggered statetable, a
// compressed Liberty file.
size_t addClockGateModels(naja::NL::NLLibrary* primitives,
                          const std::filesystem::path& libertyPath);

// Replaces instances of modeled ICG cells, in any database of `universe`, by
// one instance of the latch and one of the gate. An instance is replaced
// only when its gated clock reaches nothing but clock pins of flops that
// capture on its rising edge, through buffers and inverters of the same
// module; any other instance stays opaque. So does an instance named in
// `keptInstanceNames`: one the caller refers to by name, such as a
// verification boundary. Returns the number replaced.
size_t expandClockGates(naja::NL::NLUniverse* universe,
                        const std::set<std::string>& keptInstanceNames = {});

}  // namespace KEPLER_FORMAL
