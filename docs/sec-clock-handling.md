# SEC Clock Handling

This document describes the SEC clock modeling behavior used during
gate-level and RTL-level sequential equivalence checking.

## Scope

This page covers how Sequential Equivalence Checking (SEC) interprets clocks
while extracting the transition system. Combinational LEC/CNF generation does
not use this clock model.

The SEC clock model is enabled by default for all SEC engines:

- `k_induction`
- `imc`
- `pdr`

There is currently no user flag to opt in or opt out. The model is part of the
normal SEC extraction path.

## Clock Concepts

SEC uses three related concepts during extraction:

| Concept | Meaning |
| --- | --- |
| Clock domain | The abstract source clock event, usually a top-level clock input such as `clk`, `clk_i`, or `wb_clk_i`. |
| Clock phase | Whether the state samples on the positive or negative edge of that domain. |
| Clock enable | Non-carrier gating logic that controls whether a state update happens on that event. |

Pure routed clock signals are treated as carriers of an abstract event. A
carrier can preserve phase or invert phase, but it is still part of the same
clock domain when it is structurally derived from the same source clock.

## Carrier Classification

SEC seeds clock domains from top-level clock carrier inputs and then classifies
clock expressions seen on sequential clock pins.

A top-level input is a clock carrier when its fanout, through combinational
cells, reaches only sequential clock pins. Loads without a function, such as
antenna diodes or clock-load cells with unconnected outputs, are ignored. The
port name does not matter; the names `CLK`, `CLOCK`, `CK` and `*_CLK` remain a
fallback for inputs whose fanout cannot be classified. Routed clock-tree names
(`clkbuf*`, `clknet*`) only extend a clock domain that was found this way; they
never create one on their own.

The classifier handles:

- Direct clock use, such as `clk`.
- Phase inversion, such as `!clk` or an inverter in the clock tree.
- Simple routed clock-tree buffers.
- Gated clock expressions, such as `clk & en` or equivalent Boolean forms.

For each recognized clock-pin expression, SEC extracts a `ClockEvent`:

```text
domain = source clock carrier
phase  = posedge or negedge
enable = optional non-clock gate
```

If an expression is a pure ungated carrier, it can become part of the carrier
set used by downstream state/input matching. If an expression is gated, the
carrier part remains the event and the non-clock part becomes enable logic.

Expressions with dynamic or ambiguous polarity are not treated as pure carriers.
For example, a clock expression whose edge can behave as both rising and falling
depending on data logic is not collapsed into one abstract carrier event.

## Positive And Negative Edge State

SEC no longer assumes every clock carrier represents the same positive-edge
tick.

For every modeled state bit, extraction records the clock event associated with
its clock pin:

- positive-edge primitive flops become positive-phase state
- negative-edge primitive flops become negative-phase state
- inverted clock-tree carriers flip the phase
- inverted clock pins such as `CKN` or `CLKN` flip the phase

Within a single clock domain, SEC samples the macro-cycle just before the
positive edge:

1. Positive-edge state updates first.
2. Negative-edge state in the same domain then sees the positive-edge next-state
   value.

This supports common same-domain communication patterns where a positive-edge
flop feeds a negative-edge flop through combinational logic.

SEC does not compose state updates across different clock domains. Different
domains are treated as separate abstract events, and crossings are not assumed
to be synchronous.

## Gated Clocks

SEC tries to avoid modeling a gated clock as an unrelated clock domain.

When the clock expression can be classified as:

```text
clock_event & enable
```

the state transition becomes:

```text
next_state = enable ? data_next : current_state
```

This means clock gating is modeled as state enable behavior rather than as a
new independent clock when the gate is combinational and fully modelable.
A clock-gate latch in the gating cone is reduced by the latch abstraction in
[sec-latch-support.md](sec-latch-support.md): open while the clock is
low, it is a wire when the state samples, so the gate's enable becomes the latch
data input. Latches the abstraction does not cover stay opaque.

### Integrated Clock-Gating Cells

A library clock gate (ICG) is one cell with the latch inside. Liberty marks it
with `clock_gating_integrated_cell` and describes it with a `statetable` for
the latch and a `state_function` for the gated clock:

```text
statetable ("CLK ENA SE", "IQ") { table : "L L L : - : L, ..., H - - : - : N"; }
pin (GCLK) { direction : output; state_function : "CLK & IQ"; }
```

Naja keeps no function for such a cell. SEC reads the statetable itself and,
when it is one latch that is open while its clock pin is low, replaces every
instance of the cell by that latch and a gate carrying the `state_function`.
The rules above then apply: the latch becomes a wire, the gate becomes a clock
enable, and a flop behind the cell takes the form

```text
next_state = (ENA | SE) ? data_next : current_state
```

so a clock-gated netlist compares against RTL that writes the enable as a mux
on the data. The run log reports the cells modeled per Liberty file and the
number of instances expanded.

The latch is a wire only for state that samples on the rising edge of the
gated clock. An instance is therefore replaced only when its gated clock
reaches nothing but clock pins of flops that capture on that edge, through
buffers and inverters of the same module: a falling-edge flop behind an
inverter qualifies, a falling-edge flop on the gated clock itself does not.
Any other instance stays opaque, as does a cell the model does not cover:

- a clock gate for falling-edge flops (its latch is open while the clock is
  high);
- an output that is not a gated copy of the clock: one that falls when the
  clock rises, such as an inverted gated clock, or one that follows the latch
  while the clock is low;
- a statetable with edge-triggered rows, several internal nodes, or a next
  state that depends on the current state;
- a gated clock that leaves its module, or reaches a latch, a memory or other
  logic;
- an instance named in `--set-as-boundary`, which stays the cell it is;
- a cell from a compressed Liberty file (`.gz`, `.zip`).

The log names each instance or cell left opaque and why. Outputs behind them
are skipped and reported like those behind any opaque cell.

`--clock-gate-model false` (`clock_gate_model: false`) turns the model off and
leaves every such cell opaque. Combinational LEC does not use it: a
clock-gated flop and an enable-mux flop are different there. The model is
applied by the command-line driver while it loads the Liberty files and the
designs; designs handed over through the Python API are verified as they were
loaded.

## Complex Clock Trees

Routed clock trees can contain many internal clock nets and buffer cells. SEC
tries to classify those as carriers of the same abstract clock event instead of
publishing every clock-tree leaf as a separate input.

The current handling is structural first:

- Top-level clock inputs seed carrier events.
- Internal clock-tree expressions are classified by Boolean relation to known
  carriers.
- Inverters flip phase.
- Buffers preserve phase.
- Gated expressions keep the same domain/phase and expose an enable.

Name-based clock-tree hints are intentionally limited and design-local. They are
used only to recover local routed clock branches inside the same extracted
design, not to relate internal names across the two compared designs.

## Multiple Clock Domains

SEC can extract state from multiple clock domains in one design. However, this
is not full clock-domain-crossing verification.

For same-domain, mixed-phase state, SEC composes positive-edge and negative-edge
updates as described above. For different domains, SEC does not invent an
ordering or timing relation.

During extraction, SEC records a domain for each modeled state bit when the
clock event is known. It then checks whether:

- a state next-state cone depends on state from another clock domain
- an observed output cone depends on state from more than one clock domain

When an observed output cone spans more than one extracted clock domain, SEC
marks that output as skipped with the `multi-clock-domain` origin. This is a
coverage decision: SEC avoids silently proving such an output under an
unstated synchronous-clock assumption.

## CDC Boundary

SEC clock handling is not CDC analysis.

The current model does not prove synchronizer correctness, metastability
containment, handshakes, asynchronous FIFO protocols, or clock-ratio
relationships. If a design has CDC logic, SEC may still prove outputs whose
cones remain inside one abstract clock domain, and it may skip outputs whose
cones span domains. That skip should be read as "outside this SEC timing
model", not as a CDC pass or fail.

## Reports And Coverage

When SEC skips observed outputs due to clock-domain mixing, the main log reports
partial coverage:

```text
SEC checked-output coverage: <percent> (<covered>/<existing> covered/existing outputs)
```

The skipped reason is summarized as `multi-clock-domain`.

When `--report-skipped-pos` or `report_skipped_pos: true` is enabled, SEC may
also write:

```text
skipped_multi_clock_domain_pos.txt
```

This file lists observed outputs skipped because their cones span more than one
extracted clock domain. The report is emitted only when there are entries.

Other skip reasons can hide a multi-clock-domain concern. For example, if an
output is already skipped because it depends on reset-unanchored state, that
reset-unanchored reason may be the reported coverage blocker.

## Practical Reading Of Results

Use this interpretation when reviewing SEC output:

| Observation | Interpretation |
| --- | --- |
| Full coverage on a single-clock design | SEC modeled the clocked state under one abstract event. |
| Full coverage on a multi-clock design | Covered observed outputs did not require an unsupported cross-domain timing assumption under the extracted model. |
| Partial coverage with `multi-clock-domain` skips | Some observed outputs crossed extracted clock domains and were intentionally skipped. |
| Partial coverage with reset/unconnected/loop skips | Those issues blocked coverage first; the run does not prove whether skipped outputs are also CDC-sensitive. |
