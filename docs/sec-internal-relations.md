# SEC Internal Relation Learning

Before the output proof, SEC tries to prove that pairs of registers always hold
the same value. Every proved pair is handed to the proof engines as an
invariant, which turns a proof over two unrelated state spaces into a proof
over one. The pass is optional (`learn_internal_relations`, default `true`) and
lives in `src/sec/proof/InternalRelations.cpp`.

A learned relation is never trusted because of how it was guessed. It is used
only after a proof over the original transitions and the boot state.

## 1. Guessing candidate pairs

`learnInternalStateRelations` pairs two state bits when either holds:

- **Same full name.** Both have the same display name: the full instance path
  from the top, the pin name and the bit, for example `core.csr.reg_mie.Q[0]`.
  Memory cells are named from the instance path plus cell and bit index. The
  strings must match exactly, so a flow that renames or flattens instances
  produces no candidate for the affected registers.
- **Same output net.** Both state bits' output pins sit on nets with the same
  name: the parent instance path, the net name and the bit. Netlist flows keep
  net names where they rename or re-cluster instances, for example multibit
  flop clustering turns `u0.w[0][13]$_DFF_P_/QN` into `_tray_size2_1098/QN0`
  while the net `_00641_` stays. This guess is applied in addition to the pin
  name, never instead of it. A state bit is its output pin, so two pins on one
  net are proposed as equal; the tray-versus-flop polarity of the Liberty state
  variable does not enter.
- **Same driver.** Inside one design, both have the identical next-state
  expression.

Candidates are not guessed from simulation. After the first proof the
candidates are extended by propagation, see section 7.

In the dual-rail encoding a candidate equates both rails of the two registers.
With `allow_x_equality_in_internal_relations: false` it also requires both
registers to be binary-defined.

## 2. Base case

A candidate is kept only if its equalities already hold in the initial state
assignments, and both registers have a transition.

## 3. Size gate

The gate applies only to the candidates that are sent to the solver (section
4). Their logic is counted, each shared node once. Above about 8 million nodes
those candidates are given up for the round; candidates proved structurally are
not affected. Counting stops at the limit, so a large design is never built in
memory just to be measured. The gate is an engineering limit, not a technique
from the papers.

## 4. Inductive step

`proveInternalRelations` looks for the largest subset of candidates that is
jointly 1-step inductive: if every pair is equal now, every pair is equal after
one transition.

- **Structural proof.** Each transition is rewritten under the hypotheses: a
  register merged into another reads that register's leaf. The rewrite is
  numbered in a temporary table rather than built: two transitions get one
  number exactly when they are the same formula after the rewrite, and that
  candidate is proved without the solver. The table is freed when the round
  ends. Expressions themselves are never freed, so building the rewrite of a
  whole design would keep a second copy of its logic for the rest of the run
  (about 3.5 GB on nangate45_black_parrot). Short AND/OR
  chains are rebuilt from their sorted operands, because physical optimization
  swaps the pins of symmetric gates and that re-brackets the chain. Sorting the
  operands of an AND or an OR cannot change the function, so a pin swap on a
  gate that is not symmetric, a mux for example, still compares as different.
  The step is skipped when X equality is disallowed, since definedness is not
  structural. On nangate45_black_parrot all 666,417 candidates are proved this
  way in one round. This is an implementation choice, not a technique from the
  papers.
- **Speculative reduction.** The candidates that remain go to the solver. The
  hypotheses are applied by literal substitution: both registers of a pair
  share one current-frame literal. For these candidates only, the rewritten
  transitions are built as expressions and the solver reads them, so logic the
  two sides have in common is encoded once and needs no search.
  (Mony et al., DAC 2005; Mishchenko et al., ICCAD 2008, section 3.2.)
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

The engines receive the survivors as a list of pairs. The exact IMC engine also
reads their conjunction as one formula, so that formula is built only for the
problems that engine handles: binary encoding, or dual rail with at most 12
state bits. Larger dual-rail problems are proved per output from the pairs.

With `KEPLER_SEC_DIAG` set, each round prints how many candidates were active,
proved structurally and dropped, and each solver partition prints its status
and the work it used.

## 7. Propagation to a fixed point

A register one stage behind proved pairs is not a candidate by name or driver:
its next-state formula reads its own design's registers, so it never hashes
equal to its counterpart. Once the proved pairs are merged, that is, one
register of each pair is read in place of the other, the two formulas become
one, and the pair can be proposed. This is the structural-similarity iteration
of van Eijk (TCAD 2000) and the merge-and-rehash of AIG sweeping (Kuehlmann,
Paruthi, Krohm, Ganai, TCAD 2002).

After the guessed candidates are proved (possibly none of them),
`learnInternalStateRelations` repeats:

1. Build the merge map from the proved pairs (a union of the equalities, so
   chains of pairs share one representative). In the first round nothing is
   merged yet, so that round is the cross-design structural guess: registers
   whose transitions already read the same inputs.
2. Number every state bit's transition under that merge with the same table
   section 4 uses for structural proofs, one number per rail. The table is
   freed when the round ends, so no rewritten logic is kept.
3. Put state bits with equal numbers in one bucket, across both designs. In
   the dual-rail encoding the two rails are ordered by number, so a register
   whose rails match another's crosswise lands in the same bucket with the
   opposite orientation: the complement of a value is its rails exchanged (X
   stays X), and such a pair is proposed with the rails swapped. These are the
   complemented edges of Kuehlmann et al.; in the binary encoding a complement
   has no rail to swap and is not proposed.
4. Propose every unrelated member of a bucket against its first member, then
   re-prove the whole candidate set with section 4. The result is still one
   jointly inductive certificate; nothing proved in an earlier round is kept
   without being re-proved with the rest.
5. Stop when a round proposes nothing new, when the proof does not grow, or
   after 32 rounds.

Each round reaches one register stage deeper, so a pipeline of depth d needs
about d rounds. Registers in a feedback loop are reached as a group: all
members of a cycle are proposed in the round their inputs are merged, and the
joint induction proves them together.

The `SEC summary` diagnostics line reports the candidates by source
(`name_candidates`, `net_name_candidates`, `driver_candidates`,
`propagated_candidates`) and the number of propagation rounds.

## What it does and does not help

- Two designs with the same register names and near-identical logic behind
  them are proved quickly, because the merged transitions collapse.
- Designs that really differ keep only the registers the difference cannot
  reach; the output check still reports the counterexample.
- Flows that re-cluster or rename registers, for example multibit flop
  clustering, keep their net names, so the net-name guess and the propagation
  rounds recover the pairs the pin names lose. On the asap7 aes-mbff
  synthesis-versus-final pair (562 registers, 26 of them multibit trays holding
  64 bits), pin names alone proposed 498 pairs, 368 of them were refuted
  because their cones read the unpaired tray bits, 130 survived, and PDR was
  still inconclusive after 27 minutes. With the net-name guess the 64 tray bits
  are proposed as well, all 562 pairs are proved structurally in one round
  (1,124 rail equalities), and PDR proves the design at k = 1 in 3 seconds.
  The plain aes pair is unchanged: 530 pairs by pin name, k = 1, 2 seconds.
- Equivalent designs whose logic was rebuilt (for example synthesis netlist
  versus final netlist) can prove few pairs: one pair the solver leaves
  undecided is dropped, and so is every pair that reads it. sky130hs_riscv32i
  keeps 64 of 2,112 rail equalities this way, while the synthesis-versus-final
  nangate45_black_parrot pair proves all of them structurally. Mishchenko et
  al. address the general case with signal correspondence from random
  simulation, which also relates internal nodes. That is not implemented.
