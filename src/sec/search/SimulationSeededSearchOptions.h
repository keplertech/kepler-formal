// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>

namespace KEPLER_FORMAL::SEC {

// Simulation-seeded bounded search for sequential counterexamples.
//
// A long deterministic start-up sequence (an invalidation walk, a memory
// initialization loop) pushes every reachable difference hundreds of cycles
// away from the initial state, where neither bounded unrolling nor
// property-directed reachability gets in reasonable time. Three-valued
// simulation with held random inputs is cheap at that depth, so the search
// first simulates the proof transition system for `prefixCycles` cycles,
// keeps the end state that determines the most state bits, and only then
// runs a bounded search for a binary output difference from that state.
// Any trace it finds is a concrete input sequence of the original problem;
// finding nothing proves nothing, and the usual proof flow follows.
struct SimulationSeededSearchOptions {
  // Simulated prefix length. Zero disables the search.
  size_t prefixCycles = 0;
  // Number of held-input simulation trials; the best end state is kept.
  size_t trials = 32;
  // Bounded search depth after the simulated prefix.
  size_t searchDepth = 40;
  // Cycles the first random input vector is held before the second one
  // takes over for the rest of the prefix (a reset-then-run pattern).
  size_t holdCycles = 2;
  // Random seed for the held input vectors.
  uint64_t seed = 1;
  // Property check schedule of the bounded search: 0 checks once, at the
  // full depth, with a one-shot solve; n > 0 also checks the first two frames
  // and then every n frames with an incremental solver. Every intermediate
  // check is an unsatisfiability proof the search does not need, so the
  // default is the single final check.
  size_t checkEvery = 0;

  bool enabled() const {
    return prefixCycles != 0 && trials != 0 && searchDepth != 0;
  }
};

}  // namespace KEPLER_FORMAL::SEC
