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
   search keeps the state that determines the most state bits; among equals it
   prefers the run in which the state kept changing. A design held in reset or
   stuck before its initialization leaves memories and data paths unknown, so
   the run that completed its start-up sequence wins without any name-based
   knowledge of the design.
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
