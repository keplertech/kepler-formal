# Generalized phase analysis

This directory separates the paper-aligned analysis from Kepler's netlist
frontend and latch interpretation.

1. `NormalizedTransitionSystem` copies a Naja design into an owned,
   implicitly-clocked transition relation. A latch state is normalized to
   `enable ? data : old_state`, with modeled clear and preset behavior outside
   that update. Flip-flop next-state expressions retain Kepler's existing
   implicit-edge-step interpretation. Reset-state values are explicit input to
   normalization; clear or preset pins are not treated as initial values.
2. `TernarySimulation` starts from that reset vector, drives unconstrained
   inputs with `X`, and simulates until the state vector reaches a cycle. Time
   steps are sequential, while independent next-state roots are evaluated in
   parallel.
3. `PeriodicSignalAnalysis` implements the periodic-generator portion of
   Bjesse and Kukula, *Automatic Generalized Phase Abstraction for Formal
   Verification* (ICCAD 2005): Boolean periodic state signals receive minimum
   generator words and a small global unfolding factor `N` is selected. A
   signal that is `X` anywhere in the reset stem or cycle is conservatively not
   classified as a generator.
4. `LatchPhaseProfile` is Kepler-specific interpretation. It substitutes the
   proven generator value at each phase into every latch transparency
   predicate. This preserves residual local gates, so `clock & gate1` and
   `clock & gate2` can share a carrier phase without being conflated.

The actual `N`-step phase-abstraction transform described later in the paper is
not performed here. This module only produces the reset trace, generator
evidence, selected phase count, and latch profiles required by that future
transform.

## Concurrency and lifetime

Naja's `NLUniverse`, DNL, and the temporary `BoolExpr` caches are process-global.
Calls through this frontend are therefore serialized; callers must also avoid
running unrelated DNL clients concurrently with collection. Before releasing
DNL, the frontend copies all expressions into a dense owned DAG. Simulation,
generator discovery, and phase classification operate only on that immutable
representation and are reentrant. Parallel stages write to preallocated
index-addressed slots and merge in stable order. `maxConcurrency = 1` or
`KEPLER_NO_MT` selects the serial path.

## Conservative limits

The frontend reports an incomplete model for unsupported sequential cells,
structured memories, opaque combinational roots, combinational loops,
multi-driver roots, or unknown simultaneous clear/preset behavior. Unknown
initial state is valid, but normally prevents that state from becoming a
periodic generator. No temporal ordering, reachability-based enable relation,
or pseudo-clock is inferred from latch enable syntax.
