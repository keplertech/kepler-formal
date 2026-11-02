# SEC Sequential Models

This document describes which Naja sequential models Kepler Formal can use
when extracting a transition system for Sequential Equivalence Checking (SEC).

## Naja Model Contents

A Naja sequential model records:

- whether the element is a flip-flop or latch
- the clock or enable expression in `clockedOn`
- each state's next-state, clear, and preset expressions
- the relationship between physical output pins and modeled state

Models can come from Naja DB0 primitives, supported frontend lowering,
explicit Python primitive descriptions, or Liberty construction. A cell can
still lack a usable model when its sequential behavior is absent, incomplete,
or not representable by the current constructor.

Combinational primitives, including the word-level arithmetic and read-only
memories the SystemVerilog frontend produces, are expanded by the cone builder
as described in [sec-primitive-expansion.md](sec-primitive-expansion.md).

## Flip-Flops

SEC consumes flip-flop models when their state expressions and physical output
mapping are valid. The clock expression is classified using the clock model
described in [sec-clock-handling.md](sec-clock-handling.md), and the extracted
state is used by all three SEC engines.

Invalid or incomplete state fragments are opaque at output-terminal
granularity. Other independently modeled outputs of the same instance remain
eligible unless their own dependency cones reach an opaque terminal.

## Latches

The full description of latch handling, with the decision rules, reporting
and limits, is [sec-latch-support.md](sec-latch-support.md). This section is
the short form.

A latch is level-sensitive. For an active-high latch, its behavior is:

```text
next_q = enable ? data : q
```

This is different from a flip-flop's edge-triggered update. Naja represents
that distinction with `SequentialModel::Kind::Latch`, and its built-in
`naja_dlatch` has such a model.

Kepler SEC does not implement generic level-sensitive transition semantics.
Instead it applies phase abstraction (Hasteer, Mathur, Banerjee, "Efficient
equivalence checking of multi-phase designs using phase abstraction and
retiming", ACM TODAES 1998) to latches whose enable resolves to one phase of a
clock through the clock carrier:

- Enable active while the clock is high: the latch closes at the falling edge
  and is kept as negative-phase state with `next = data_in`. A gated enable,
  `clk & g`, becomes the usual clock-enable form `next = g ? data_in : q`.
- Enable active while the clock is low: the latch is open when the
  positive-edge state samples, so it is a wire. Its state variable is replaced
  by its `data_in` formula everywhere and removed from the state. This is the
  integrated clock-gate case: the gated clock's enable becomes the latch data
  input. A gated enable, `!clk & g`, as in a latch-based memory with a
  registered write enable, keeps the latch closed while `g` is zero: it is
  kept as positive-phase state with `next = g ? data_in : q`, and its
  positive-edge readers see that value, the data-enabled latch form below.
  Negative-phase readers see the value it held through the high phase, which
  is its next value.

The published method assumes alternating phases, so SEC does not reduce a latch
whose data input reads another latch; such latches and latches with a clear or
preset that can fire stay opaque and are reported. A latch whose enable is not
a clock phase takes the data-enabled latch model below. A clear or preset whose
cone folds to constant false, for example a pin tied to its inactive level,
never fires and does not disqualify the latch. SEC does not infer latch
behavior from cell or pin names. The abstraction is on by default and
`--latch-abstraction false`
(`latch_abstraction: false`) restores the previous behavior where every latch
output is opaque.

### Data-enabled latches

A latch whose enable is data logic rather than a clock phase, typically
inferred from an incomplete assignment in RTL, is modeled as state on the
global step with next value `enable ? data_in : q`, and its readers see that
value, so the latch is transparent within the step. This is the synchronous
latch model, the same lowering open-source formal flows apply to latches
before a proof. A latch reading another data-enabled latch takes that latch's
current-step value, so cascades ripple through within the step. Such a latch
stays opaque when its output reaches a clock, when a clear or preset can fire,
or when it sits in a loop of latches that can be open together. The model is on
by default and `--data-latch-model false` (`data_latch_model: false`) leaves these
latches opaque.

## Read-Only Memories

A memory primitive with read ports but no write port holds no state. SEC keeps
it out of the structured-memory model and expands each read-data bit as a
constant table lookup over the address bits, as described in
[sec-primitive-expansion.md](sec-primitive-expansion.md). Memories with write
ports keep the structured-memory model.

## X Constant Literals

An RTL `x` literal, typically a `case` default such as `default: y = 'x;`, is
kept by the frontends as an X-constant net. SEC models each such net as a
state bit with no initial value whose next state is itself: under the
dual-rail encoding it is X forever, so a binary value on the other side is not
a binary-defined difference. This is the usual don't-care reading of such
literals. Under the binary encoding the bit behaves like any uninitialized
register, and dependent outputs are skipped. Z literals remain unsupported and
skip the cones that read them.

## Opaque Outputs

Opacity is strict and local to an output terminal. During backward cone
construction, reaching an opaque internal output stops construction of that
requested top-level output. SEC does not create a free, shared, or substitute
symbol for the opaque value and does not compare the affected output.

Unaffected top-level outputs remain eligible. The main result reports partial
checked-output coverage, and `--report-skipped-pos` writes details to
`skipped_opaque_cells_pos.txt`, including the top-level output, opaque cell and
pin, and reason.

The complete SEC reporting and flag behavior is documented in
[sec-flags-spec.md](sec-flags-spec.md).
