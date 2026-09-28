# Latch Support Algorithm

This is a Boolean event model for latch behavior and sequential equivalence.
It describes a conditional method, not a guarantee that every latch circuit
settles or can be reduced. The [behavior guide](sec-latch-implementation.md)
contains the primitive diagrams and event examples.

## 1. Approach

1. Represent each latch as remembered storage plus a transparent mux.
2. Propagate changes through connected gates and storage elements.
3. Identify feedback loops and event-dependent regions.
4. Prove that every allowed propagation episode settles within a finite bound.
5. Prove its final retained state is independent of internal ordering choices.
6. Unfold those internal rounds into one boundary-to-boundary transition.
7. Compare both designs under the same external stimulus.

Phase reduction is optional and comes afterward. The cited research supports
individual steps; no single article proves this complete combination.

## 2. Latch and event behavior

For an active-high latch with data D, enable E, remembered bit H, and output Q:

```text
Q = E ? D : H
next(H) = Q
```

H is abstract storage, not a flip-flop connected to the hardware clock.
The output follows data while open and retains history while closed [R1].
Reset and preset follow the element's declared priority.

The block is evaluated as one primitive. Combining its register and mux does
not change these equations. Splitting them into independently scheduled
elements can change event timing and requires a separate equivalence argument.

An external transaction changes permitted inputs at a settled boundary.
Their new values remain fixed until internal propagation finishes. Multiple
inputs may change together unless explicitly restricted. A formal round is
not a hardware clock tick or a chosen physical sampling interval.

Each round:

1. Freeze current signals, previous signals, and stored values.
2. Evaluate elements activated by changed inputs.
3. Gates evaluate their Boolean functions. Storage elements process every
   permitted order of their changed input pins, retaining storage between visits.
4. Latches apply transparency; flip-flops capture only on the appropriate
   modeled edge, subject to asynchronous controls.
5. Publish each element's final output/storage tuple together with the others.
6. Activate receivers of changed signals and repeat.

Intermediate storage updates within one element's pin ordering are retained;
only its final outputs are published for that round. Changes between rounds
remain visible to downstream elements. This primitive granularity is part of
the model, not a claim to reproduce arbitrary physical glitches [R2].

```mermaid
flowchart LR
    A[External transaction] --> E[Evaluate activated elements]
    E --> C[Commit updates together]
    C -->|Signal changes| E
    C -->|No pending work or error| B[Settled observation]
```

A simultaneous data change and latch closure can retain either old or new
data, depending on event order. Do not silently choose the favorable result.

## 3. Starting state and reset

Unspecified inputs and storage are arbitrary Boolean values chosen once, not
repeatedly resampled or forced to zero. They are not literal unknown-valued
signals. Matching external stimulus is shared across the designs; internal
storage is independent unless an explicit relation constrains it.

Initially project remembered storage to its outputs, then force an evaluation.
Latches apply transparency and asynchronous controls; flip-flops retain state
unless asynchronous controls apply. No clock edge is invented merely because
a clock starts high. Subsequent generated edges are processed normally.

Auxiliary initial values on other internal signals must not determine the
retained result. For each genuine starting input/storage choice, prove settling
and independence from those auxiliary values and permitted event orders.
This Boolean startup convention is an adaptation, not a universal hardware rule.

Reset is applied only when requested. Reset cycles mean hardware clock cycles,
not propagation rounds. A cycle schedule must represent assertion, clock edges,
intervening settling, and release. See the
[reset-cycle behavior](sec-latch-implementation.md#reset-cycle-adapter).

## 4. Feedback loops and regions

A loop is a directed path returning to itself through latches and logic.
Find these paths using both data and control dependencies.

- An open latch feeding itself unchanged keeps its previous value. The equation
  Q = Q alone loses that history.
- An open latch with inverting feedback can oscillate.
- A closed latch can break a loop. To classify a region as inactive, prove
  that every feedback cycle is broken under the allowed conditions.
- A snapshot with no active loop does not prove settling if controls can
  change during propagation.

A scheduling region is not necessarily a feedback loop: it may contain several
loops and connecting logic [R5, R9]. Flip-flop data paths can supply boundaries
when capture cannot occur inside an episode; generated clocks and asynchronous
controls cannot be cut without justification.

If one region briefly raises another latch's enable, sending only the final
low level loses the capture. Preserve that boundary event sequence, enlarge
the region, or prove those intermediate events irrelevant.

Merging all connected elements avoids such interfaces but can create enormous
regions. A failure then affects every observation retained in that region.
Finer regions require valid event interfaces; connectivity alone does not prove
they may settle independently.

## 5. Proving and unfolding settling

Retain the complete state needed for future behavior: storage, current and
previous signals, pending activity, and errors. A state is settled only when
no modeled work remains.

Before accepting a bound, establish:

- **Valid starts:** startup includes every admitted initial condition and
  establishes a boundary invariant B.
- **Complete entries:** every B-state and allowed transaction has an admission
  successor, and Entry includes every admission outcome. Establish this
  separately from the bounded query below.
- **Total behavior:** every allowed transaction and internal state has a
  successor. Undefined behavior becomes an explicit persistent error, not
  a missing execution.
- **Stable padding:** a settled state repeats unchanged for the rest of its
  episode. Errors persist and never count as settled.
- **Closure:** completed episodes starting in B return to B.

For candidate bound K, ask whether any admitted entry and permitted internal
ordering can remain unsettled after K rounds:

```text
Entry(x0, u)
AND T(x0, x1, u) AND ... AND T(x[K-1], xK, u)
AND NOT Settled(xK)
```

Here u is held external stimulus; T is one complete internal round.
Proving this formula impossible establishes the bound [R3, R4].
Keeping errors alive prevents short failing executions from disappearing before K.

A finite state space alone does not establish convergence. A reachable cycle
outside settled states permits endless propagation. If no such cycle exists,
the longest nonsettled path gives a bound; exhaustive exploration can still
be impractical. Alternatively, test candidate bounds or prove a decreasing
ranking measure.

A failed candidate may need more rounds; resource exhaustion proves neither
convergence nor oscillation. An apparent failure from an over-approximation
must be shown reachable before calling it a circuit defect.

### Unique retained result

Settling alone is insufficient. Run two copies from the same pre-transaction
state and stimulus, with independent admission and internal choices. Require:

```text
B(q) AND Allowed(u) AND Episode(q, u, a) AND Episode(q, u, b)
IMPLIES a = b
```

Compare the complete retained state, not only current outputs: a later input
may expose hidden differences. Startup needs the analogous check with genuine
initial values shared and auxiliary choices independent.

If uniqueness fails or remains unproved, this deterministic approach does not
support the affected behavior. It must not align favorable choices between designs.

### Build once, reuse each transition

After proving the bound and uniqueness, replace one propagation episode with
K chained copies of its internal-round relation. Early completion uses unchanged
padding. Build this chain once and reuse it for successive external transactions.

Why this is exact under the stated conditions: every permitted episode finishes
by K and can be padded without changing its result; every unfolded execution
therefore corresponds to a permitted completed episode. Closure extends that
argument across transaction sequences [R3, R4, R9].

The bound concerns complete episodes. Combining independently guessed local
bounds by a maximum or sum is not a proof of the whole episode.

## 6. Local scheduling and parallel evaluation

A safe scheduling reduction preserves the reference round number:

1. Each region reads only complete inputs and history from the same round.
2. Each predecessor supplies an update or an explicit unchanged indication;
   silence is not evidence that nothing changed.
3. Regions publish updates for the next round, preserving boundary transitions,
   shared choices, and their effects on clocks, enables, and resets.
4. Declare global settling only when all pending work is complete.

Independent evaluations may run together. Evaluation order must not become
circuit behavior. By induction, equal complete inputs at one round produce
the same permitted updates at the next.

Never mix an old value on one input with a new value from another round:
doing so can invent a pulse at a reconvergent gate. Likewise, one producer's
choice must remain the same at all its receivers.

This allows local scheduling without requiring each region to settle privately.
Skipping rounds, publishing only endpoints, or reordering dependent events
needs stronger preservation arguments [R6–R8].

## 7. Equivalence and unsupported behavior

Give both designs the same allowed external transactions and specified
initial/reset relation. Compare corresponding settled observations, even when
the designs need different numbers of internal rounds.

Two obligations remain separate:

- Every admitted episode completes.
- Required outputs agree at the corresponding boundaries.

Agreement only when both sides finish can pass vacuously if one never finishes.
Do not discard that execution. Without reset, do not assume independent
uninitialized storage happens to agree.

Unsupported or unproved behavior remains opaque, with affected observations
excluded and explicitly reported—not proved equivalent or replaced by shared
arbitrary values. Dependencies include controls as well as data. An optional
strict policy stops on any opacity; it is off by default.

This method does not cover arbitrary delays, metastability, unrestricted input
changes during settling, undefined controls, or general unknown/multi-driver
semantics. Large bounds may be impractical even when behavior is valid.

## 8. Optional phase reduction

Only after defining correct transitions, look for deterministic periodic state
signals and use them to reduce repeated phases [R10]. Preserve residual gating,
capture effects, and the observation boundary.

This is not one phase per latch, and internal round numbers are not automatically
hardware clock phases. Failure to find a period leaves the original model intact.

## References

- **R1.** Håkan Hjort, [On Applying Model Checking in Formal Verification](https://fmcad.org/FMCAD22/presentations/00%20-%20tutorials/02_hjort.pdf),
  FMCAD 2022 tutorial, slides 41 and 44–50: latch register/mux construction and
  feedback hazards; no general safe sampling interval or settling bound.
- **R2.** Raffelsieper, Roorda, Mousavi,
  [Model Checking Verilog Descriptions of Cell Libraries](https://doi.org/10.1109/ACSD.2009.18),
  ACSD 2009. Accessible treatment:
  [Cell Libraries and Verification, Chapter 3](https://pure.tue.nl/ws/files/3499974/717717.pdf#page=24),
  2011, especially pp. 20–27: pin history, evaluation, and updates.
  Published fixed-order/input restrictions and unknown-value conventions are
  not adopted implicitly here; the experiments do not establish whole-design scalability.
- **R3.** Claessen and Sörensson,
  [A Liveness Checking Algorithm that Counts](https://www.cs.utexas.edu/~hunt/fmcad/FMCAD12/fmcad2012.pdf#page=59),
  FMCAD 2012, Section III-A: finite-state eventuality bounds.
  Applying the result to complete latch episodes is our specialization.
- **R4.** Schuppan and Biere,
  [Efficient Reduction of Finite State Model Checking to Reachability Analysis](https://www.schuppan.de/viktor/VSchuppanABiere-STTT-2004.pdf),
  2004: liveness-to-safety reasoning; no guarantee that a particular loop settles.
- **R5.** DeVane,
  [Efficient Circuit Partitioning to Extend Cycle Simulation Beyond Synchronous Circuits](https://cecs.uci.edu/~papers/compendium94-03/papers/1997/iccad97/pdffiles/03a_1.pdf),
  ICCAD 1997, Sections 3.5–4: trigger-based regions, latches, generated clocks,
  and asynchronous controls. The principal algorithm restricts combinational feedback.
- **R6.** Lang and Mateescu,
  [Partial Order Reductions using Compositional Confluence Detection](https://cadp.inria.fr/publications/Lang-Mateescu-09.html),
  FM 2009: conditional composition and reordering results, not permission to
  serialize arbitrary latch regions.
- **R7.** Neele, Valmari, Willemse,
  [A Detailed Account of the Inconsistent Labelling Problem of Stutter-Preserving Partial-Order Reduction](https://lmcs.episciences.org/7709/pdf),
  2021, Section 5: conditions for preserving observations under reduction.
- **R8.** McDonald and Bryant,
  [Symbolic Timing Simulation Using Cluster Scheduling](https://www.cs.cmu.edu/~bryant/pubdir/dac00a.pdf),
  DAC 2000, Section 3.2: local symbolic event queues with ordering safeguards;
  not a proof of this zero-delay latch model.
- **R9.** Alur and Henzinger,
  [Reactive Modules](https://www.cis.upenn.edu/~alur/FMSD99.pdf),
  1999, Section 6.2, pp. 37–38, Figure 13: transparent-latch feedback and
  internal-round abstraction; not universal convergence of arbitrary loops.
- **R10.** Bjesse and Kukula,
  [Automatic Generalized Phase Abstraction for Formal Verification](http://www.perbjesse.com/iccad05.pdf),
  ICCAD 2005: periodic reduction after correct latch/clock modeling and loop
  resolution, not a replacement for those foundations.
