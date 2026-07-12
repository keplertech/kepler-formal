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
new independent clock. Clock-gate latch data that has been folded during SEC
extraction is substituted before the clock event is classified, so common
integrated clock-gating structures can still expose the intended enable.

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
