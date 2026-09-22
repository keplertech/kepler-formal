# SEC Internal Relation Learning

Before the output proof, SEC tries to prove that pairs of registers always hold
the same value. Every proved pair is handed to the proof engines as an
invariant, which turns a proof over two unrelated state spaces into a proof
over one. The pass is optional (`learn_internal_relations`, default `true`) and
lives in `src/sec/proof/InternalRelations.cpp`.

A learned relation is never trusted because of how it was guessed. It is used
only after a SAT proof over the original transitions and the boot state.

## 1. Guessing candidate pairs

`learnInternalStateRelations` pairs two state bits when either holds:

- **Same full name.** Both have the same display name: the full instance path
  from the top, the pin name and the bit, for example `core.csr.reg_mie.Q[0]`.
  Memory cells are named from the instance path plus cell and bit index. The
  strings must match exactly, so a flow that renames or flattens instances
  produces no candidate for the affected registers.
- **Same driver.** Inside one design, both have the identical next-state
  expression.

Candidates are not guessed from simulation.

In the dual-rail encoding a candidate equates both rails of the two registers.
With `allow_x_equality_in_internal_relations: false` it also requires both
registers to be binary-defined.

## 2. Base case

A candidate is kept only if its equalities already hold in the initial state
assignments, and both registers have a transition.

## 3. Size gate

The logic behind all remaining candidates is counted, each shared node once.
Above about 8 million nodes the learner is skipped and the output proof runs
exactly as it does without learning. Counting stops at the limit, so a large
design is never built in memory just to be measured. tinyrocket has about 0.9
million nodes; nangate45_black_parrot, with 666,543 candidates, exceeds the
limit and would otherwise need over 13 GiB and tens of minutes. The gate is an
engineering limit, not a technique from the papers.

## 4. Inductive step

`proveInternalRelations` looks for the largest subset of candidates that is
jointly 1-step inductive: if every pair is equal now, every pair is equal after
one transition.

- **Speculative reduction.** The hypotheses are applied by literal
  substitution: both registers of a pair share one current-frame literal. The
  two sides' transitions are then encoded over the same literals, so identical
  logic collapses structurally and needs no search. (Mony et al., DAC 2005;
  Mishchenko et al., ICCAD 2008, section 3.2.)
- **Partitioning.** One-step register correspondence needs a single time frame,
  so the candidates are split into partitions bounded by solver variables.
  Every hypothesis is merged in every partition and each candidate is proved in
  exactly one, so splitting loses no relation. A large design is split rather
  than skipped. (Mishchenko et al., section 3.3.)
- **Variables on first use.** A partition reads a small part of the design, so
  its solver creates a variable only when the encoded logic first mentions a
  symbol, not one per symbol per frame. This is an implementation choice, not a
  technique from the papers. It lowers memory and encode time per partition.
- **Query.** Each partition asks whether some candidate in it can differ in the
  next frame.
  - UNSAT: all of its candidates hold under the hypotheses.
  - SAT: the counterexample is replayed (below).
  - Undecided within budget: each pair of that partition is asked separately
    on the same solver with its own budget, and only the pairs that stay
    undecided are dropped. (Mony et al., sections 2 and 4.1.)

## 5. Refinement by simulation

A counterexample is replayed on the original transitions as one of 64 parallel
patterns; the other 63 are random states that also satisfy the hypotheses. The
replay runs for up to 16 steps, and every candidate seen differing on a valid
pattern is dropped. One counterexample therefore refines all candidates, not
only the pairs the solver model happens to separate. (Mony et al., section
3.2.)

Simulation only drops candidates. It never proves one.

## 6. Fixed point

Dropping a hypothesis weakens every other proof, so the rounds repeat with a
fresh encoding. The survivors are returned only after a whole round drops
nothing; if the round limit is reached first, nothing is returned.

## What it does and does not help

- Two designs with the same register names and near-identical logic behind
  them are proved quickly, because the merged transitions collapse.
- Designs that really differ keep only the registers the difference cannot
  reach; the output check still reports the counterexample.
- Equivalent designs whose logic was rebuilt (for example synthesis netlist
  versus final netlist) prove few pairs: register equalities alone are often
  not inductive there. Mishchenko et al. address this with signal
  correspondence, which also relates internal nodes. That is not implemented.
