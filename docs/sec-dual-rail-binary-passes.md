# Binary Passes of Dual-Rail SEC

Dual-rail SEC proves equivalence on rail values, where a register can be X.
Two passes work on the 0/1 functions of the designs instead: the binary round
and signal correspondence on the binary functions. This document says what
each part of them does, what it hands to the dual-rail proof, and which
published article each step comes from. The code is in
`src/sec/proof/InternalRelations.cpp` and
`src/sec/proof/SignalCorrespondence.cpp`.

Both passes are on by default and switched together by
`binary_relation_round` (`--binary-relation-round <bool>`). They need
internal relation learning, the dual-rail encoding and an engine that keeps
the lazy transition store (k-induction or PDR).

References are listed in full at the end and cited by a short name.

## 1. The model the passes serve

Every register is two rails, may-be-one and may-be-zero; X is both rails set
and the pair with no rail set is not a value. The proof starts from the boot
state, X where a design gives no value, and two outputs differ only in a
frame where both are known and unequal. This is the dual-rail model and the
steady-state equivalence of [KH03], sections 4 and 6.1.

Rail values are computed gate by gate. Two facts about that computation are
used below:

- **A known value is certain.** If a signal is 0 or 1 on the rails, it has
  that value for every 0/1 completion of the unknown registers. This is the
  monotonicity of three-valued simulation: a refined input gives a refined
  output ([HCC01], section 2.6, after [BBS95]).
- **An X can be false.** A signal can be X on the rails although it has one
  value for every completion, because the computation does not see that two
  X come from the same register. [HCC01] (appendix) separates such an X from
  a true one, and [vE00] notes that combinational synthesis does not preserve
  three-valued equivalence for this reason.

The second fact is why two registers that are equal in 0/1 can differ on the
rails, one X and the other known, and why the dual-rail engines alone can
fail on logic that synthesis reshaped.

## 2. Binary round: register pairs equal in 0/1

The candidate pairs of [sec-internal-relations.md](sec-internal-relations.md)
(section 1) that the rail check leaves unproved are tried again on the 0/1
next-state functions, which the lazy store keeps before lifting.

The proof is an induction over the pairs together ([vE00]): with every pair
assumed equal in one frame, each pair must be equal in the next, and a pair
that fails is dropped until nothing fails. The base case is the boot state: a
pair is kept when both registers have the same boot value or both are
unknown. Two unknown registers count as equal at the start in the
flip-flop-pair induction of [HCC01] (section 5.1, Lemma 5) too; [SBMB10]
(Theorem 1) states what such a proof gives, equality from every common
initial state. Section 6 says what rests on this start.

A pair proved here is used in two ways, sections 3 and 4.

## 3. What such a pair gives on the rails

A pair that is equal in 0/1 but not on the rails is **never merged**. In
three-valued logic a merge is only correct for signals that are equal X
included; merging two signals of which one can be X while the other is known
excludes that combination and "could lead to incorrect verification results"
([HCC01], section 5). [WLC21] (section V) gives the same rule for
equivalence checking with X values: such a pair cannot be merged unless the
relation holds in both directions, and is added as clauses instead.

An example. Register `r` is 1 after reset and 0 afterwards in both designs.
Design 1 computes that 0 from a register `m` that is never reset, in a shape
where `m` cancels out, so its rails show X; design 2 writes the 0 directly.
The two `r` are equal in 0/1. If design 2's `r` were given the logic of
design 1's, it would show X as well, and an output that is 1 in design 1 and
`r` in design 2 would compare 1 against X: a real difference would be lost.

For each such pair `proveRailRules` asks, on the rails, for the strongest
rule that holds:

| Rule | Meaning | Clauses, for registers L and M | Source |
|---|---|---|---|
| Covering | M is at least as known as L: whenever L is known, M has the same value | `M.one -> L.one`, `M.zero -> L.zero` | [HCC01] section 2.4; [WLC21] Definition 1 |
| Compatible | L and M never show opposite known values | `L.zero or M.one`, `L.one or M.zero` | [KH03] section 6.1, applied to the pair |

Covering is tried in both directions. A pair that is equal on the rails, X
included, needs none of this: the rail learner has already proved and merged
it (sections 4 to 8 of
[sec-internal-relations.md](sec-internal-relations.md)).

The rules are proved on the rails, not derived from the 0/1 proof: by
induction from the boot state, all rules together, with the certified rail
equalities as hypotheses ([vE00]; the proof obligation is changed from
equality to the rule, as [KH03] section 6.1 does for outputs). A rule that
fails, or that the solver budget cannot decide, is dropped, and the rest
counts only once a round drops nothing. A pair with no rule left gets
nothing.

A proved rule is a constraint ([HCC01], section 5): its two clauses are
added to every frame of the k-induction and PDR queries and of the
signal-correspondence run on the rails. Each design keeps its own logic. A
clause is left out where a query does not hold both symbols; that only loses
pruning, because a rule removes states no run reaches.

## 4. Outputs closed in 0/1

An output whose two 0/1 functions are equal once the registers of the proved
pairs are read as one is closed without the dual-rail engine
(`dischargeBinaryEqualOutputs`).

Deciding a three-valued verdict by looking for a 0-against-1 difference from
known states is the method of [HCC01] (Lemmas 2 and 4). The step to the
rails is the first fact of section 1: if both outputs were known and unequal
on the rails, they would be unequal for every completion, which the 0/1
proof excludes. The verdict of section 1 asks no more than that ([KH03],
section 6.1).

## 5. Signal correspondence on the binary functions

Outputs still open are tried by signal correspondence ([vE00]), decided with
a SAT solver ([BC00]) and with induction over several frames ([SSS00]), on
the 0/1 functions. It relates signals and not only registers, which is what
a design needs that keeps in a register a value the other design recomputes.
The run starts at the boot state, goes through the reset window, and reports
outputs only; its register classes are not handed to the engines. The steps
are in section 11 of
[sec-internal-relations.md](sec-internal-relations.md). An output it proves
is closed for the reason of section 4.

## 6. The common start value of unknown registers

Sections 2 and 5 start the two registers of a candidate pair on one value
when both are unknown at boot. What rests on it:

- **The rules of section 3 do not.** They are proved on the rails, where
  every unknown register starts as X on its own.
- **The closed outputs of sections 4 and 5 do**, in this sense: the 0/1
  proof covers the completions in which the paired registers agree. The
  verdict is still about the rails. A known rail value holds for every
  completion (section 1), one of the covered ones included, so two outputs
  that are both known on the rails are equal.

The start value only matters for state that is still unknown when the
comparison starts. A run whose reset leaves such state says so in the log,
and says that these passes are on; see
[sec-reset-bootstrap.md](sec-reset-bootstrap.md). With
`binary_relation_round: false` none of this document runs.

## 7. Diagnostics

With `KEPLER_SEC_DIAG=1` or `KEPLER_SEC_SUMMARY_STATS=1`:

```text
SEC summary: binary_relation_round candidates=... binary_pairs=... covering_clauses=... compatible_clauses=... discharged_outputs=.../...
```

`candidates` are the pairs the rail check left unproved, `binary_pairs` the
ones proved in 0/1, the two clause counts what section 3 proved for them
(two clauses per rule), and `discharged_outputs` the outputs of section 4.
Signal correspondence prints its own line per run, with `semantics=binary`.

## References

- [BBS95] D. Brand, R. A. Bergamaschi, L. Stok, "Be careful with don't
  cares", ICCAD 1995. Cited through [HCC01].
- [BC00] P. Bjesse, K. Claessen, "SAT-based verification without state space
  traversal", FMCAD 2000.
- [HCC01] S.-Y. Huang, K.-T. Cheng, K.-C. Chen, "Verifying sequential
  equivalence using ATPG techniques", ACM Transactions on Design Automation
  of Electronic Systems 6(2), 2001.
- [KH03] Z. Khasidashvili, Z. Hanna, "SAT-based methods for sequential
  hardware equivalence verification without synchronization", Electronic
  Notes in Theoretical Computer Science 89, 2003.
- [SBMB10] H. Savoj, D. Berthelot, A. Mishchenko, R. Brayton, "Combinational
  techniques for sequential equivalence checking", FMCAD 2010.
- [SSS00] M. Sheeran, S. Singh, G. Stålmarck, "Checking safety properties
  using induction and a SAT-solver", FMCAD 2000.
- [vE00] C. A. J. van Eijk, "Sequential equivalence checking based on
  structural similarities", IEEE Transactions on Computer-Aided Design
  19(7), 2000.
- [WLC21] Y.-N. Wang, Y.-R. Luo, P.-C. Chien, P.-L. Wang, H.-R. Wang, W.-H.
  Lin, J.-H. R. Jiang, C.-Y. R. Huang, "Compatible equivalence checking of
  X-valued circuits", ICCAD 2021.
