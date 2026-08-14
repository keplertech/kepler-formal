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
state. A generic latch output is opaque and any requested top-level output
whose cone reaches it is skipped.

## Clock-Gate Latch Rewrite

SEC has one narrow latch-specific rewrite for recognized clock-gating
structures. This is not a general latch transition model.

A candidate must have one latch-like output, one data input, and a gate input.
If the primitive has a Naja sequential model, its kind must explicitly be
`Latch`; a `FlipFlop` model never enters this rewrite. The gate must
structurally trace to a pure clock carrier, and both the data and gate
dependency cones must be fully modelable. When all checks pass, SEC
substitutes the latch output with its data expression while reconstructing the
clock-gate enable. The latch output is then removed from the extracted SEC
state and is not published as an environment input.

This rewrite captures the clock-gating abstraction used by SEC: the latch
holds an enable stable around the active clock edge, while the consuming
flip-flop remains the actual modeled state. If candidate recognition or
dependency validation fails, the latch output remains opaque.

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
