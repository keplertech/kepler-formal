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

## Flip-Flops

SEC consumes flip-flop models when their state expressions and physical output
mapping are valid. The clock expression is classified using the clock model
described in [sec-clock-handling.md](sec-clock-handling.md), and the extracted
state is used by all three SEC engines.

Invalid or incomplete state fragments are opaque at output-terminal
granularity. Other independently modeled outputs of the same instance remain
eligible unless their own dependency cones reach an opaque terminal.

## Latches

A latch is level-sensitive. For an active-high latch, its behavior is:

```text
next_q = enable ? data : q
```

This is different from a flip-flop's edge-triggered update. Naja represents
that distinction with `SequentialModel::Kind::Latch`, and its built-in
`naja_dlatch` has such a model.

Kepler SEC does not currently implement generic level-sensitive transition
semantics. It therefore does not consume a Naja latch model as ordinary SEC
state. Every latch output is opaque, including latches used in clock-gating
structures, and any requested top-level output whose cone reaches it is
skipped. SEC does not infer latch behavior from cell or pin names.

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
