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
  // Target enlargement (Baumgartner, Kuehlmann and Abraham, CAV 2002): widen
  // the target by this many cycles before the forward search. Structurally,
  // the last enlargeSteps frames unroll only the target's cone; with
  // enlargeCubes > 0 and no structured memories, up to that many generalized
  // preimage cubes per step are enumerated and added to the target. Zero
  // disables enlargement.
  size_t enlargeSteps = 0;
  size_t enlargeCubes = 0;
  // Waypoint chaining (lighthouses: Ganai, Aziz and Kuehlmann, DAC 1999;
  // Ho et al., ICCAD 2000): instead of one deep search, up to waypointHops
  // short searches of waypointDepth frames, each toward the real target or a
  // waypoint, chained from the concrete state the previous one reached. A
  // waypoint is a register in the target's backward cone, within
  // waypointLayers register steps of it at first, that has not changed value
  // so far; the goal is that it changes (or becomes known). Zero hops
  // disables chaining.
  size_t waypointHops = 0;
  size_t waypointDepth = 8;
  size_t waypointLayers = 1;
  // Solver conflict budget of one hop; a hop that exceeds it counts as
  // unreachable and the chain tries the next mode. Zero means no budget.
  size_t waypointConflicts = 300000;
  // Abstraction-guided search: the abstract model keeps the registers within
  // abstractionLayers steps of the divergence predicate, its trace pins the
  // boundary registers in the concrete search, and the registers behind a
  // refuted pin join the model, for up to abstractionIterations rounds. A
  // concrete check stops after abstractionConflicts solver conflicts (0: no
  // budget). Zero layers disables it.
  size_t abstractionLayers = 0;
  size_t abstractionIterations = 16;
  size_t abstractionConflicts = 1000000;
  // Start the abstract model from the smallest design instance that encloses
  // the differing registers (and more than them) instead of the predicate's
  // register support.
  bool abstractionByInstance = false;
  // Seed screening: a bounded search of probeDepth frames, under
  // probeConflicts solver conflicts, from every simulated end state in rank
  // order before the deep search; the first that diverges becomes the start
  // state. 0 disables it.
  size_t probeDepth = 0;
  size_t probeConflicts = 200000;

  bool enabled() const {
    return prefixCycles != 0 && trials != 0 && searchDepth != 0;
  }
};

}  // namespace KEPLER_FORMAL::SEC
