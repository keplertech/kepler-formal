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
  flop clustering moves a flop's output pin onto a tray cell with a new
  instance name and pin name while the net it drives keeps its name. This
  guess is applied in addition to the pin
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
  (several gigabytes on a large design). AND/OR chains are rebuilt from their
  sorted operands, the canonical form of structural hashing (Kuehlmann and
  Krohm, "Equivalence checking using cuts and heaps", DAC 1997), because
  physical optimization swaps the pins of symmetric gates and that re-brackets
  the chain. Wide reduction trees are chains of hundreds of operands, so the
  rebuild applies to chains of any practical width. Sorting the operands of an
  AND or an OR cannot change the function, so a pin swap on a gate that is not
  symmetric, a mux for example, still compares as different. The step is
  skipped when X equality is disallowed, since definedness is not structural.
  On a large design with several hundred thousand candidates, all of them are
  proved this way in one round.
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

## 8. Merging the proved relations into the proof

Proved relations used to reach the engines only as frame constraints: two
copies of the same cone, one per design, plus the facts that their registers
agree. Proving that both copies agree at the outputs was then a SAT search over
the whole cone, which is slow exactly when the cone is large and identical.

After learning, every proved same-frame equality merges its two state symbols
into one representative, and every transition and output formula is rewritten
under that merge. The expression factory hash-conses, so logic that is
identical over related registers becomes one node and the output comparison is
a syntactic equality. Logic that differs structurally but is functionally equal
is proved once per node, not once per output cone that reads it. Both registers
remain state; only reads of the merged one change. This is the register merge
of van Eijk (TCAD 2000) followed by structural hashing, the SAT-sweeping step of
Kuehlmann et al. (TCAD 2002).

Lazily remapped problems are materialized first, since the merge changes every
formula's support. Designs above a size limit keep the lazy store and the
constraint-only behavior.

## 9. Structural discharge of outputs

The pairs are a certificate that, in every reachable state, each register of
design 1 holds the value of its design-0 partner. Van Eijk (TCAD 2000) uses
such a certificate twice: as hypotheses for the engine, and to merge the
product machine, reading one register of each pair in place of the other so
that logic the two designs share collapses into one copy. The first use is
sections 1 to 6, the merge of the proof problem is section 8. This section
applies the merge to the observed outputs.

After the survivors are added to the problem, every certified same-frame
equality (the proved pairs and the definitional complement pairs) is put into
one equality class with a single representative. The two cones of each
observed output are then numbered with the table of section 4 under that
merge, so that they get one number exactly when they are the same formula once
every merged register is read through its representative. In the dual-rail
encoding both rails must match. An output whose cones match is proved: its
obligation becomes the constant true, the property, bad and induction formulas
are rebuilt from the remaining outputs, and it is flagged in
`structurallyProvedOutputs`. No engine query is spent on it. An output whose
cones differ anywhere keeps its obligation unchanged.

This is what makes a design compared with itself, or with a netlist that kept
the structure of its arithmetic, independent of the engine budget: once the
registers are paired, the outputs are the same formula and nothing is left to
search. A cone that synthesis restructured is not caught by the hash; it goes
to the engine as before.

With `KEPLER_SEC_DIAG` set, the pass prints the number of outputs discharged
and the number of symbols merged.

## 10. Binary round

The rail check compares the ternary values that the per-gate dual-rail
lifting computes, and that lifting is a shape-dependent over-approximation of
the exact unknown value: two shapes of one Boolean function can give X and a
defined value for the same unknown inputs. A register pair whose logic the
flow restructured can therefore fail the rail check although its two
next-state functions are equal.

Candidates the rail check leaves unproved are retried in binary semantics on
the next-state functions the lazy store keeps before lifting: a register reads
by its may-be-one rail, an input by its shared symbol. The base case is the
common initial state (both registers hold the same value, or both are
unknown), the step is the joint induction of section 4 with every rail-equal
pair merged as well, and the query is the combinational equivalence a LEC
solves, so it is allowed a hundred times the rail check's budget. A pair
proved this way is equal from every common initial state (van Eijk, TCAD
2000; Savoj et al., FMCAD 2010, Theorem 1).

A proved pair is not merged on the rails, where one of its registers can be
X while the other is known. It gets the strongest rule that is proved on the
rails, one register at least as known as the other or the two never showing
opposite known values, and the rule is added to the engines' queries as
clauses. An output whose two functions agree in binary semantics once the
proved pairs are read as one is discharged like a structurally identical one
(section 9). Each step, and the article behind it, is in
[sec-dual-rail-binary-passes.md](sec-dual-rail-binary-passes.md). Every
other pair and output stays with the dual-rail engines unchanged.

The round is on by default and switched by `binary_relation_round`
(`--binary-relation-round <bool>`). It needs learning, the dual-rail encoding
and an engine that keeps the lazy transition store (k-induction or PDR). The
summary line `SEC summary: binary_relation_round candidates=... binary_pairs=... covering_clauses=... compatible_clauses=... discharged_outputs=...`
reports the pairs it retried, the pairs it proved, the clauses of the rules
proved for them and the outputs it discharged. The same switch covers the run of section 11 on the binary
functions. The common initial state matters only for state that is unknown
when the comparison starts: a run whose reset leaves such state says so in the
log, and says that these passes are on (see
[sec-reset-bootstrap.md](sec-reset-bootstrap.md)).

## 11. Signal correspondence

Sections 1 to 10 relate registers. That is not enough when one design keeps in
a register a value the other recomputes, for example a register copy of a
memory word that is read in place of the word while it is valid. The copy has
no partner. The two memories are equal only while the outputs are, because
what is written depends on the output, and the outputs only while the
memories are. No register pair is inductive by itself, and an engine has to
find the missing fact one state at a time.

Outputs that sections 9 and 10 leave open are therefore tried by signal
correspondence (van Eijk, TCAD 2000), decided with a SAT solver instead of
BDDs (Bjesse and Claessen, FMCAD 2000). `proveSignalCorrespondence` in
`src/sec/proof/SignalCorrespondence.cpp` works on the problem the engines
receive, the rules of section 10 included:

1. **Signals.** The set F of the paper: every register, both cones of every
   open output (each rail in the dual-rail encoding) and every register's
   next-state function. Any set is sound; this one is small.
2. **Start state.** The boot state. With a reset bootstrap it is the state the
   reset window leaves, found by propagating constants through the window
   with the reset ports active and every other input unknown. A register
   that stays unknown is left free, so every state the window can leave is
   covered. After the window the reset ports are held inactive, as in the
   engines.
3. **First partition.** Sequential simulation with random inputs from the
   start state groups the signals that never differ. Each signal is
   normalised by its value in the first pattern, so signals that always hold
   opposite values share a class. In the dual-rail encoding this is how
   "this register is binary-defined" is stated: its two rails are
   complements. A run goes on while the partition still refines, because a
   memory behind an address counter is only reached after many frames.
4. **Start (eq. 4 of the paper).** Members of a class must be equal in the
   start state for every input.
5. **Step (eq. 5).** With every class assumed in one time frame, the members
   of each class must be equal in the next. The assumed equalities are plain
   constraints of the current frame, the state compatibility condition of
   the paper; nothing is merged before it is proved. No register holds the
   empty rail value, the reset ports rest at their inactive value and the
   equalities certified by sections 1 to 10 hold: these are constraints of
   both frames.
6. **Refinement.** A counterexample is replayed on the logic and splits every
   class it separates. The step is repeated with a fresh solver until a round
   splits nothing. The classes are then a signal correspondence relation:
   their members are equal in every state reachable from the start state.

Classes are asked in groups of 64, one query per group, so a round that changes
nothing costs a query per group. A group the solver cannot decide in one query
is halved, down to single classes, which are asked member by member. A class
proved in a round may be used by the later queries of that round. A signal the
solver leaves undecided is removed from its class, never assumed.

A relation can fail to be inductive over one time frame although it holds. A
design that reads a value in one phase and writes it back in a later one, a
program counter read nibble by nibble and incremented three phases later for
example, only shows in the reading phase whether a register copy of that value
is right; in the other phases a state with a wrong copy passes for a state of
the relation. The step is then taken over several frames (Sheeran, Singh and
Stalmarck, FMCAD 2000; Bjesse and Claessen combine it with signal
correspondence): the classes must hold in the first k frames from the start
state, and, assumed in k consecutive frames, in the next. With k at least the
length of the operation, the phase that exposes a wrong copy is always among
the assumed frames. The pass tries k = 1, 2, 4 and 8 while outputs are open.

If outputs are still open and the logic is small enough, the same is repeated
with F set to every node of the next-state and output logic.

An output whose two cones already take different values in a state the start
state reaches, in simulation or in the check of the first frames, is not tried
again: no larger set and no deeper induction proves it. A design that really
differs therefore costs one run before the engine looks for the
counterexample.

### On the binary functions

On the rails the relation inherits the weakness section 10 describes for
register pairs: the lifting is shape dependent, so two shapes of one function
can hold X and a defined value in the same state. One rewritten adder is
enough. Every register behind it then differs between the designs in where it
is unknown, few classes are left, and an output that is X on one side and
defined on the other is in no class at any induction depth, although nothing
is wrong with it: the obligation of an output is only that two defined values
agree.

The pass therefore runs first on the binary functions the lazy store keeps
before lifting, the functions of the binary round: a register reads by its
may-be-one rail, an input by its shared symbol, an output by the function of
each design. The run differs from the one above in three points:

1. **Start state.** The boot state, taken through the reset window. A
   register is known there when its rails hold a binary value, and the two
   registers of a candidate pair of section 1 that are both unknown hold one
   value, the common initial state of section 10. The window is simulated,
   and unrolled in the solver in front of the first compared frame, with the
   reset ports active and every other input free. The start state is then
   exactly what the window leaves, including what it loads from the inputs.
2. **Signals.** Registers and outputs only, half as many as on the rails, and
   with induction over 1, 2, 4, 8 and 16 frames: two machine cycles of eight
   clock phases are what a value stored by one instruction and read by the
   next one spans.
3. **Result.** Outputs only. The register classes hold from the common boot
   states and not from every boot state, so they are not handed to the
   engines.

An output proved this way is discharged like one proved on the rails. The
argument is the one of ternary simulation (Seger and Bryant, Formal Methods in
System Design 6(2), 1995): a rail value that is defined is the value the
output takes for every binary completion of the unknown boot state. The run
proves the two binary outputs equal for every completion in which the paired
registers agree, and at least one such completion exists. Two outputs that
are both defined on the rails therefore hold the value of that completion in
both designs, and agree. Any choice of pairs is sound for this reason; the
pairs of section 1 are the ones that make equivalent designs equal.

The rail runs follow for the outputs that are left: those that are equal
only because both sides are unknown, which no binary function shows.

The result is used twice:

- An output whose two cones are in one class with the same polarity, on both
  rails in the dual-rail encoding, is discharged like a structurally
  identical one (section 9).
- Registers of one class are added to the same-frame equalities, and are then
  merged like the pairs of sections 1 to 10. A pair is added only as a
  complete value (both rails, plain or crosswise for a complement) and only
  if the two registers are also known and equal in every frame before the
  start state, because the engines read these equalities in the reset window
  too.

The pass is bounded by engineering limits (`SignalCorrespondenceLimits`), not
techniques from the papers:
16384 state symbols, about two million nodes counting every unrolled frame
(262144 for all signals), 64 rounds, 8 induction frames on the rails and 16
on the binary functions, and a solver budget of a couple of minutes shared by
the runs of one kind, twice the ticks for the binary runs. A query is limited
by its conflicts and decisions, not by its ticks, because propagation alone
costs ticks in proportion to the design. When a limit is reached nothing is
reported and the engines run as before.

It is on by default and switched by `signal_correspondence`
(`--signal-correspondence <bool>`). It needs learning and the dual-rail
encoding. Its classes may relate signals that are unknown, X on both sides, so
it is also left out with `allow_x_equality_in_internal_relations: false`. The
runs on the binary functions need what the binary round needs as well:
`binary_relation_round` and an engine that keeps the lazy transition store.
The line `SEC summary: signal_correspondence scope=... frames=... completed=... signals=... classes=... rounds=... solver_calls=... proved_outputs=... state_equalities=... seconds=... semantics=... stopped_by=...`
reports each run: `semantics` is `binary` or `rails`, and `stopped_by` names
the limit that ended a run with `completed=0` (`model` for a size limit or an
operator the pass does not read, `rounds`, or the part of the solver budget
that ran out: `conflicts`, `decisions` or `ticks`).

## What it does and does not help

- Two designs with the same register names and near-identical logic behind
  them are proved quickly, because the merged transitions collapse.
- Designs that really differ keep only the registers the difference cannot
  reach; the output check still reports the counterexample.
- Logic the flow restructured, such as a resynthesized arithmetic cone or a
  mux rewritten as and-or logic, is proved by the binary round, where the
  rail check of the first rounds is shape sensitive.
- Flows that re-cluster or rename registers, for example multibit flop
  clustering, keep their net names, so the net-name guess and the propagation
  rounds recover the pairs the pin names lose. On a synthesis-versus-final
  pair where about a tenth of the register bits were clustered into multibit
  trays, pin names alone proposed most pairs, but the majority were refuted
  because their cones read the unpaired tray bits, and PDR was still
  inconclusive after tens of minutes. With the net-name guess the tray bits
  are proposed as well, every pair is proved structurally in one round, and
  PDR proves the design at k = 1 in seconds. A pair without clustering is
  unchanged.
- Equivalent designs whose logic was rebuilt (for example synthesis netlist
  versus final netlist) can prove few pairs: one pair the solver leaves
  undecided is dropped, and so is every pair that reads it. One such pair
  keeps only a few percent of its rail equalities this way, while another
  proves all of them structurally. Section 11 relates internal signals found
  by random simulation, but only within its size limits; a design above them
  keeps the behavior of sections 8 to 10, which merge proved registers,
  recognize outputs whose cones become identical and prove registers and
  outputs whose functions are equal.
- A design that keeps a register copy of a value the other recomputes, such
  as a shadow of a memory word, has outputs and registers that only hold
  together. Section 11 proves them in one induction from the state after
  reset. Without a reset bootstrap the start state is the boot state, where
  every register may still be unknown, and less is proved.
- The same with logic rebuilt beside it, a processor core whose program
  counter gained a shadow and whose adder became carry-lookahead for example,
  differs on the rails in where values are unknown, and the rail runs of
  section 11 prove nothing at any depth. The binary runs prove it, with
  enough frames for the longest path from a value stored to the same value
  read.
