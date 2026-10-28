# SEC Latch Support

This document describes how Sequential Equivalence Checking (SEC) handles
latches: which latches it can model, the form each one takes in the extracted
transition system, what stays opaque, and how that is reported. The short
reference version lives in [sec-sequential-models.md](sec-sequential-models.md);
the clock terms used below (domain, phase, carrier, clock event) are defined in
[sec-clock-handling.md](sec-clock-handling.md).

## 1. Summary

SEC does not simulate level-sensitive timing. Each latch is classified once,
during extraction, into one of four forms, or left opaque:

| Enable of the latch | Form in the transition system |
| --- | --- |
| A clock phase, open while the clock is high (`clk`, or `clk & g`) | Negative-phase state. The latch closes at the falling edge. |
| A clock phase, open while the clock is low (`!clk`) | A wire. The latch is open when the positive-edge state samples. |
| A gated low phase (`!clk & g`) | Positive-phase state that is transparent within the step. |
| Data logic, not a clock phase | Data-enabled state that is transparent within the step. |

The first three forms are the *latch abstraction*; the fourth is the
*data-enabled latch model*. Both are on by default:

| CLI | YAML | Default | Effect when off |
| --- | --- | --- | --- |
| `--latch-abstraction <bool>` | `latch_abstraction: <bool>` | `true` | Every latch output is opaque, the behavior before latch support. This also disables the data-enabled model, because latches never reach transition construction. |
| `--data-latch-model <bool>` | `data_latch_model: <bool>` | `true` | Latches whose enable is not a clock phase are opaque; the abstraction still applies to the others. |

The run log states both settings at startup
(`SEC latch abstraction: enabled` and `SEC data-enabled latch model: enabled`).

## 2. Input: the latch model

Naja records a latch as a sequential model of kind `Latch`. Its `clockedOn`
expression is the level-sensitive enable, `nextState` is the data, and the
model may carry `clear` and `preset` expressions and the mapping from the
physical output pins, true and complemented, to the state. Models come from
Liberty `latch` groups, from the built-in latch primitive the SystemVerilog
frontend infers for incomplete assignments, or from Python primitive
descriptions. SEC never infers latch behavior from cell or pin names: a cell
without a latch model is an ordinary opaque cell.

For an active-high latch the level-sensitive behavior is

```text
q_next = enable ? data : q
```

The forms below are the ways SEC approximates that behavior at the step
boundaries where it compares the two designs.

## 3. Classification

For each latch, transition construction proceeds in this order.

1. **Connectivity.** A pin with no net at all is recorded as a connectivity
   skip with the reason `its pin is not connected`; the latch state and the
   state depending on it are skipped and the rest of the design is still
   proved. A pin on a net without a driver is a no-driver skip, as for any
   cell.
2. **Enable classification.** The enable cone is built and classified by the
   clock model. When the enable is a clock carrier, its inversion, or a clock
   gated by data (`clk & g`, `!clk & g`), the result is a clock event: a
   domain, a phase, and the gate `g` when present. Otherwise there is no clock
   event and the enable is data logic.
3. **Clear and preset.** The clear and preset cones are built. A cone that
   folds to constant false, for example a pin tied to its inactive level,
   never fires and is ignored. A clear or preset that can fire makes the latch
   opaque with the reason `latch with clear or preset is not supported`.
4. **Form.** The table in section 1 decides the form. Two checks can still
   make the latch opaque:
   - On the abstraction path, a latch whose data input reads another latch is
     opaque: `latch data input reads another latch; cascaded latches are not
     reduced`. The published method assumes alternating phases and does not
     reduce such chains.
   - On the transparent paths, a latch whose output reaches a clock is opaque
     (`latch drives a clock` or `data-enabled latch drives a clock`), because
     clock events must resolve to a top-level clock.
   When the data-enabled model is off, a latch without a clock event is opaque
   with the reason `latch enable does not resolve to one phase of a clock`.

A latch with a complemented output pin gets the complement of its next-value
expression for that pin. When the latch becomes a wire, both pins are replaced.

## 4. Latch abstraction

The abstraction follows phase abstraction (Hasteer, Mathur, Banerjee,
"Efficient equivalence checking of multi-phase designs using phase abstraction
and retiming", ACM TODAES 1998) inside the two-phase composition that the clock
model already uses for negative-edge flip-flops: positive-edge state samples
first, and negative-edge state in the same domain then sees the positive-edge
next-state values.

- **Open while the clock is high.** The latch holds the value it captured
  when the clock fell, so it is negative-phase state with `next = data`. A
  gated enable `clk & g` is folded by the clock model into the usual
  clock-enable form, `next = g ? data : q`.
- **Open while the clock is low.** When the positive-edge state samples, the
  latch is open: its output is its data input. The state variable is replaced
  by the data formula in every next-state and output expression and removed
  from the state. This is the integrated clock-gate case: the gate's enable
  becomes the latch data input, and the gated clock resolves through it.
- **Open while the clock is low with a gated enable** `!clk & g`, as in a
  latch-based storage array with a registered write enable. When the
  positive-edge state samples, the latch holds `g ? data : q`, the same
  step-transparent form as a data-enabled latch, so it is kept as
  positive-phase state with that next value. Positive-edge readers see the
  next value. Negative-phase readers keep the stored value: the latch was held
  through the high phase, and the phase composition already supplies its next
  value to them.

The method assumes that a latch on one phase feeds logic read on the other
phase. A latch reading a latch breaks that assumption, so such chains are not
reduced and stay opaque.

## 5. Data-enabled latch model

A latch whose enable is data logic, typically inferred from an incomplete
assignment in RTL, is state on the global step with next value
`enable ? data : q`. Its readers see that next value rather than the stored
one, so within a step the latch is transparent when open and holds when
closed. This is the synchronous latch lowering that open-source formal flows
apply to latches before a proof.

- **Cascades.** A data-enabled latch reading another one sees that one's
  value for the step. The next-value expressions are resolved in dependency
  order, so a cascade ripples through within the step.
- **Loops.** Latches that can be open together in a loop have no order and
  stay opaque: `data-enabled latch in a loop of latches open together`.
- **Readers.** The latch variable is replaced by the next-value expression in
  every positive-phase next-state expression and in every observed output.
  Negative-phase readers are rewritten too, except for the gated low-phase
  latches of section 4, which they read as stored. One substitution memo is
  shared across all readers, so logic common to many expressions is rewritten
  once. After the rewrite, complemented-state relations that named a
  data-enabled latch are dropped, since the complemented pin now carries its
  own expression.

## 6. Reset and initial state

Kept latch state, whether negative-phase, positive-phase or data-enabled, is
ordinary state: reset bootstrap and initial-state handling apply to it exactly
as to flip-flops (see [sec-reset-bootstrap.md](sec-reset-bootstrap.md)). A
latch reduced to a wire has no state and no initial value. Asynchronous clear
and preset are not modeled; a latch whose clear or preset can fire is opaque.

## 7. Reporting and diagnostics

An opaque latch is handled like any opaque cell: each of its physical outputs
is marked opaque with the reason above, top-level outputs whose cones reach it
are skipped, and `--report-skipped-pos` lists them in
`skipped_opaque_cells_pos.txt` with the output, cell, pin and reason. The run
log summarizes the skipped outputs. With `--latch-abstraction false` the reason
is `Naja latch sequential models are not supported by SEC`.

With `KEPLER_SEC_DIAG=1` two lines report the rewrites:

```text
SEC diag: latch abstraction wires=<n> stale_references=<m>
SEC diag: data-enabled latches=<n> rewritten_readers=<m>
```

`wires` and `data-enabled latches` count the latches of each kind;
`rewritten_readers` counts the expressions that changed; `stale_references`
counts remaining uses of a replaced variable and is expected to be zero.

## 8. Limits

- No generic level-sensitive semantics: a latch is modeled only in the forms
  of section 1. Behavior between step boundaries, including the order in which
  several transparent latches settle, is not modeled.
- Chains of clock-phase latches are not reduced.
- A latch whose output reaches a clock, a latch whose clear or preset can
  fire, and a loop of data-enabled latches stay opaque.
- Phase is relative to one clock domain; outputs whose cones span several
  domains follow the rules in [sec-clock-handling.md](sec-clock-handling.md).

## References

- Hasteer, Mathur, Banerjee, "Efficient equivalence checking of multi-phase
  designs using phase abstraction and retiming", ACM Transactions on Design
  Automation of Electronic Systems, 1998.
