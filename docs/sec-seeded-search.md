# SEC Simulation-Seeded Search

Simulation-seeded search is an optional counterexample search that runs on the
freshly built SEC problem, before internal-relation learning and the proof
engines. It targets designs whose reachable differences sit
behind a long deterministic start-up sequence, such as a cache that walks
every tag entry after reset before it accepts lookups, or a block that clears
a memory before use. Bounded unrolling and property-directed reachability
both have to pass through every one of those cycles; three-valued simulation
does not.

The search is a bug hunter, not a proof. When it finds a trace, SEC reports
`Different` with that trace. When it finds nothing, the normal proof flow
runs unchanged.

## Method

1. **Simulated prefix.** Both extracted designs are simulated for
   `prefixCycles` cycles in the same three-valued value domain the dual-rail
   proof uses (two rails per bit: may-be-1 and may-be-0). Every gate is lifted
   exactly as the proof encoding lifts it, so a simulated state is a state of
   the proof transition system. Inputs are held random vectors: one vector for
   the first `holdCycles` cycles and a second one for the rest, which covers
   the common reset-then-run pattern. Reset ports named by the reset bootstrap
   spec follow the spec instead: active for its cycles, then released.
2. **End-state selection.** Each of `trials` prefixes ends in one state. The
   search keeps the run with the highest toggle coverage, the number of
   registers that took two different known values during the prefix; among
   equals it prefers the state that determines the most state bits, then the
   run in which the state kept changing. A design held in reset toggles
   nothing, one that a held input keeps restarting its initialization toggles
   only that loop's counters, and one that completed its start-up and went on
   to work toggles the most, so that run wins without any name-based knowledge
   of the design. The ranking is a heuristic: a control input held at an
   unlucky value can leave the best-ranked run in an initialization loop. With
   `--sec-seed-probe-depth n`, every end state is first probed in rank order
   with a bounded search of `n` frames under a conflict budget, and the first
   state that diverges becomes the start state, its probe trace the suffix
   (short searches from several simulated seeds: Ho et al., ICCAD 2000).
3. **Divergence predicate.** Registers of the two designs are paired by name
   (output net name first, because RTL signal names survive a local edit, then
   display name, whose instance numbering shifts where a module changed);
   names are only hints. For every pair,
   design 2's next-state formula is renamed to design 1's symbols and compared
   with design 1's formula on the hash-consed DAG. Pairs with identical
   formulas can never take different next values from equal states; pairs
   whose formulas differ contribute "binary and different" terms, as do the
   aligned outputs whose cones differ. Their disjunction is the divergence
   predicate: on one design's state it says that the lockstep run of both
   designs splits at the next cycle. Pairs whose renamed formula still reads an
   unpaired design-2 register are left out and counted in the diagnostics.
4. **Bounded search on one copy.** From the selected state, a bounded search
   of `searchDepth` cycles over design 1 alone looks for the divergence
   predicate with an incremental SAT solver. A single copy keeps every
   bounded check an ordinary bounded-model-checking query; unrolling both
   designs would turn each check into a combinational equivalence proof of
   the two copies. Structured memories are not unrolled cell by cell: a read
   at a frame returns the data of the latest earlier write to the same
   address, under its enables and bit mask, and otherwise the selected
   state's contents at that address. Resets of the memory are ordered after
   the writes of their cycle, as in the extracted model. The predicate is
   checked once, at the full depth, as a disjunction over all frames and with
   a single solver call: every intermediate check would be an
   unsatisfiability proof the search does not need.
5. **From divergence to the outputs.** The found inputs are replayed on both
   exact models; if no observed output differs yet, random continuations of
   up to `searchDepth` cycles are simulated, and failing that a short bounded
   search over both copies starts from the diverged pair of states. When the
   selected state itself already holds different values in paired registers,
   the two-copy search runs directly from it.
6. **Confirmation.** Every reported trace, prefix and suffix, has been
   replayed on the exact extracted models with the output difference
   observed in that replay, so the result never depends on the memory
   encoding of step 4.

## Flags

| CLI flag | YAML key | Default | Effect |
| --- | --- | --- | --- |
| `--sec-seed-cycles <n>` | `sec_seed_cycles: <n>` | `0` (off) | Simulated prefix length. Enables the search when positive. |
| `--sec-seed-trials <n>` | `sec_seed_trials: <n>` | `32` | Number of held-input prefixes to simulate. |
| `--sec-seed-depth <n>` | `sec_seed_depth: <n>` | `40` | Bounded search depth after the prefix. |
| `--sec-seed-enlarge <n>` | `sec_seed_enlarge: <n>` | `0` (off) | Target enlargement by `n` cycles (see below). |
| `--sec-seed-enlarge-cubes <n>` | `sec_seed_enlarge_cubes: <n>` | `0` (off) | With enlargement, up to `n` generalized preimage cubes per step. |
| `--sec-seed-waypoint-hops <n>` | `sec_seed_waypoint_hops: <n>` | `0` (off) | Waypoint chaining: up to `n` short searches (see below). |
| `--sec-seed-waypoint-depth <n>` | `sec_seed_waypoint_depth: <n>` | `8` | Frames per hop. |
| `--sec-seed-waypoint-layers <n>` | `sec_seed_waypoint_layers: <n>` | `1` | Register distance from the target within which waypoints are taken first. |
| `--sec-seed-waypoint-conflicts <n>` | `sec_seed_waypoint_conflicts: <n>` | `300000` | Solver conflict budget per hop (`0`: none). |
| `--sec-seed-abstraction-layers <n>` | `sec_seed_abstraction_layers: <n>` | `0` (off) | Abstraction-guided search from the registers within `n` steps of the predicate (see below). |
| `--sec-seed-abstraction-iterations <n>` | `sec_seed_abstraction_iterations: <n>` | `16` | Refinement rounds. |
| `--sec-seed-abstraction-conflicts <n>` | `sec_seed_abstraction_conflicts: <n>` | `1000000` | Conflict budget of one concretization check (`0`: none). |
| `--sec-seed-probe-depth <n>` | `sec_seed_probe_depth: <n>` | `0` (off) | Seed screening: bounded search of `n` frames from every simulated end state, in rank order, before the deep search (see below). |
| `--sec-seed-probe-conflicts <n>` | `sec_seed_probe_conflicts: <n>` | `200000` | Solver conflict budget of one probe (`0`: none). |

The search needs the default `dual_rail_steady` encoding. It is skipped under
the binary encoding and whenever a formula reads memory cells outside a read
port.

```sh
kepler-formal -sv -v sec \
  --sv_design1_flist before.f --sv_design1_top top \
  --sv_design2_flist after.f  --sv_design2_top top \
  --sec-seed-cycles 300 --sec-seed-depth 40
```

The reported trace lists the held input vectors once per stretch of cycles,
followed by the per-cycle inputs of the bounded suffix and the differing
outputs.

## Held inputs

Single-bit inputs that feed the next-state function of at least half of the
registers (memory cells and read outputs not counted) are global controls such
as resets and scan or test enables. A bounded search
that may drive them turns every step into a reset, so during the suffix
searches they keep the value the selected prefix ran with. Bits of a bus are
never held: at module scope a memory response bus reaches most registers, and
holding it would fix the very values the search has to choose. Reset ports
named by the reset bootstrap spec are governed by that spec instead.

## Target enlargement

Target enlargement (Baumgartner, Kuehlmann and Abraham, CAV 2002) widens the
target before the forward search so that fewer frames have to be searched.
`--sec-seed-enlarge n` enables two forms:

- **Structural.** The forward search unrolls `depth` full frames and then `n`
  more frames in which only the backward cone of the target is encoded (the
  registers the target depends on at that distance, and the memories they
  read). A check at the last frame therefore covers every state within `n`
  cycles of the target while the remainder of the design is unrolled only
  `depth` times.
- **State sets.** With `--sec-seed-enlarge-cubes m`, up to `m` preimage cubes
  per step are enumerated by SAT: each is a conjunction of register literals
  from which the target is reached in one step under some input, generalized
  by assumption-based lifting, then blocked and the next one found. Cubes of a
  step are the targets of the next. The cubes are added to the target at every
  frame. Memory cells are not search variables, so a problem with structured
  memories keeps the structural form only, and a step whose cube budget runs
  out is kept (its cubes are real preimage states) but not expanded further.

A hit on an enlarged target is completed by an exact search of `n` frames
from the state reached, so the reported trace always ends in the original
target, and it is replayed like any other trace.

## Waypoint chaining

A long input sequence is hard to find in one query because the solver has no
sense of progress toward the target. Waypoint chaining (lighthouses: Ganai,
Aziz and Kuehlmann, DAC 1999; Ho et al., ICCAD 2000; distance guidance:
Kuehlmann, McMillan and Brayton, ICCAD 1999; Shyam and Bertacco, DATE 2006)
replaces the one deep search on the single-copy divergence problem by up to
`hops` short searches of `depth` frames each, chained from the concrete state
the previous one reached.

A target state is sampled as the start of a short free-start bounded search
(the selected state's memory contents, free registers) that ends in the
divergence predicate, so the sampled registers are consistent with reaching
it; registers the sample leaves unknown carry no goal. The guide is the
backward cone of the predicate's register part (output terms are left out,
because output cones span the whole design), and distance is disagreement
with the sample. Each hop's goal is the divergence predicate itself or one
more cone register, within `layers` register steps at first, taking its
target value while the registers the chain itself aligned keep theirs; if
that is unreachable within the hop's conflict budget the hop is retried
without the second condition, then as a plain
change of any cone register that has not changed yet, then with a wider
layer, then with a new sample (up to four). When a hop reaches a
waypoint, its inputs are replayed on both exact models and the next hop
starts from the state reached. When a hop reaches the divergence predicate,
the chained inputs form the suffix, which is replayed to the outputs as
usual. When the chain stops short, the plain bounded search continues from
the chain's last state, and the chain's inputs are prepended to its trace.

## Abstraction-guided search

The deep search is hard because the core has to be driven to produce the
addresses and control writes the changed logic needs; the changed logic
itself is small. Abstraction-guided search (localization: Kurshan 1994;
counterexample-guided refinement: Clarke, Grumberg, Jha, Lu and Veith,
CAV 2000; guidance by abstract traces: Nanshi and Somenzi, DAC 2006) works on
the single-copy divergence problem in rounds:

1. **Abstract model.** The guide is the origin of the difference: the two
   next-state formulas of each differing register pair are walked in lockstep
   and the deepest mismatching subformulas (the changed logic itself) form a
   "binary and different" predicate whose support is that logic's fan-in. The
   walk descends through a node pair only while the operand correspondence
   is unambiguous (one operand identical on both sides, in either order);
   where both operands differ, that node pair is the origin, so the
   predicate never compares unrelated subformulas. The
   registers it reads (`layers` 1), or within `layers - 1` further steps, keep
   their next-state functions; every other register is free after the
   selected state. Memories are kept with their read outputs; their cells are
   never boundary registers. With `--sec-seed-abstraction-instance true` the
   first model is instead the smallest design instance that encloses the
   differing registers (widened until it holds more than them): the component
   keeps its logic and its environment is free, the usual compositional
   abstraction (Clarke, Long and McMillan, LICS 1989; Grumberg and Long,
   TOPLAS 1994). This suits a change whose registers read, in one step, far
   more of the design than the change itself involves.
2. **Abstract trace.** A bounded search of `depth` frames on that model finds
   a divergence quickly; if it finds none, there is none within that depth.
3. **Concretization.** The values the abstract trace gave the boundary
   registers (the free registers the kept part reads) are first lifted to the
   ones the abstract model needs: with every boundary value and every input of
   the trace assumed, the negated predicate at the last frame is
   unsatisfiable, and the failed assumptions are the relevant pins (minimal
   assignments by unsatisfiable core: Ravi and Somenzi, TACAS 2004; Een,
   Mishchenko and Brayton, FMCAD 2011). The remaining pins are assumed frame
   by frame in a bounded search on the full design, which then only has to
   find inputs that reproduce them. A trace found here is a real one.
4. **Refinement.** When the pinned trace is infeasible, the registers behind
   the lifted pins join the abstract model (the abstract trace relied on each
   of them; the failed subset alone is often a single pin and would refine one
   register per round), and the round repeats; an inconclusive check widens
   the model by one register step instead.

## Diagnostics

With `KEPLER_SEC_DIAG=1` or `KEPLER_SEC_SUMMARY_STATS=1` the search prints
`SEC seeded search:` lines: simulator size, the determined state bits and
active cycles of each trial, the selected prefix, and per-frame solver status
and time of the bounded search.

## Background

The combination of simulation to reach a deep state with a bounded formal
search from that state goes back to Yuan, Shen, Abraham and Aziz (CAV 1997)
and Ho, Shiple, Harer, Kukula, Damiano, Bertacco, Taylor and Long
(ICCAD 2000). Three-valued simulation of gate networks follows Bryant
(IEEE Trans. Computers 1987) and Seger and Bryant (FMSD 1995); bounded model
checking follows Biere, Cimatti, Clarke and Zhu (TACAS 1999). Reading a memory
through its write history is the read-over-write model of Burch and Dill
(CAV 1994), in the three-valued form of Velev and Bryant (TACAS 1998).
