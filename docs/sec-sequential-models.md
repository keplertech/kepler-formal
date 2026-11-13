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

## Memories

A memory with write ports is a structured memory: a Naja memory primitive
inferred from RTL, or a library cell that carries a memory interface. SEC has
two models for it.

### Bit-level model

Every memory bit is one state bit with no initial value. A write port updates
a cell when its address matches, all its enables are 1 and, for a masked port,
the mask bit is 1. A read port selects the addressed cell. The model is exact
but grows with the memory: a 16,384 x 32 memory is 524,288 state bits per
design, and proofs stop finishing well before that size. SEC logs a warning
for each memory above 65,536 bits that it models this way.

The cells update on every step, like the registers of the model. For the
clock domain checks they belong to the domain of the memory's clock; see
[sec-clock-handling.md](sec-clock-handling.md#multiple-clock-domains).

In the `binary` encoding the unwritten cells of the two designs are
independent free values, so a design compared with itself differs as soon as
it reads a word it never wrote. SEC does not report that as a difference. When
a binary counterexample is found on a model with bit-level memories, the
dual-rail model, where an unwritten cell is X, is asked for a binary-defined
difference: within the bound of that counterexample first, then up to
`max_k`. If one exists, its trace is the counterexample reported. If none
exists, the run ends inconclusive, with the first trace and the advice to use
`dual_rail_steady`. The search up to `max_k` runs on a solver budget, so on a
large memory it can end without an answer; the message then says so.

### Memory abstraction

A memory found in both designs under the same instance path, with the same
primitive, the same instance parameters, the same connected pins and the same
clock, is not modeled. Its read data becomes one input shared by the two designs, and its
input pins become proof obligations next to the outputs of the designs. Two
identical memories that hold the same contents and are driven alike return
the same data, so proving the pins equal on every cycle justifies sharing the
read data. The proof no longer depends on the memory size.

The pins are compared as follows:

- Write enables, read addresses, clocks and resets are compared on every
  cycle.
- The address and mask of a write port are compared only while all the
  enables of that port are 1, and a data bit only while its mask bit is 1 as
  well: the same condition under which the bit-level model writes a cell. A
  write address that differs while nothing is written does not matter.
- In `dual_rail_steady` a pin differs when it is 0 in one design and 1 in the
  other; an enable that is X counts as possibly 1. Nothing is assumed about
  the memory contents: unwritten cells are X in the bit-level model, and equal
  contents are one of the cases X stands for.
- In `binary` the two memories are taken to hold the same contents when the
  comparison starts. This replaces the independent free values of the
  bit-level model.

With a reset sequence the engines compare from the first cycle after it, so
SEC separately checks the cycles of the sequence itself: on the dual-rail
model, in either encoding, no pin may differ while the reset ports are
asserted. Two designs that write a memory differently during reset are
therefore not abstracted. The `binary` encoding assumes its property on the
first cycle it compares when some state has no initial value: the first cycle
after the sequence, or the very first cycle without one. For the memory pins
SEC checks that cycle instead, and requires every pin that matters to be
defined and equal. It also checks that no output of the designs has a
binary-defined difference on that cycle, which the assumption would hide and
the bit-level model reports.

The abstraction only ever proves. If a pin is not proved equal, if a check
above fails, or if the proof finds a counterexample, which could rely on read
data no memory would return, SEC logs the reason and repeats the run with the
bit-level model. The verdict is then the bit-level one.

A memory keeps its cells from the start when:

- the other design has no memory under the same path, or it differs in
  primitive, parameters, connected pins or clock;
- an instance on its path is unnamed, so the two designs cannot align it;
- a pin breaks a rule of a boundary leaf, for example a read-data net with a
  second driver;
- a write enable is unconnected, undriven or depends on an X literal. The
  bit-level model reads such an enable as a port that never writes, while a
  compared pin would treat it as an unknown value;
- its read data is sampled on the falling clock edge, or drives a clock. The
  shared read data is one value per cycle, the one the rising edge samples; a
  falling-edge register reads after the rising-edge registers have moved;
- the instance is selected with `--set-as-boundary`, which keeps its own
  meaning.

A read whose address depends on its own data is a combinational loop and is
reported as unsupported in both models. A register that reads a memory
addressed from, or written on, another clock domain is skipped as
`multi-clock-domain` in both models.

The two models differ on one point. An X literal in the cone of a memory
address, data or mask pin leaves a bit-level memory unmodeled, and the outputs
that read it are skipped. Compared at its pins, the memory sees such a pin as
X like any other X literal: under `dual_rail_steady` it does not differ from a
binary value on the other side, and the outputs are compared.

Proved memory pins are not counted as outputs. The log names the memories
compared this way and, for each memory that keeps its cells, the reason;
`boundary_terms.txt` lists the pins of the compared ones under the roles
`memory_interface` and `memory_read_data`. The abstraction is on by default;
`--memory-abstraction false` (`memory_abstraction: false`) models every memory
bit by bit. A BTOR2 export always writes the bit-level model.

## Read-Only Memories

A memory primitive with read ports but no write port holds no state. SEC keeps
it out of the structured-memory model and expands each read-data bit as a
constant table lookup over the address bits, as described in
[sec-primitive-expansion.md](sec-primitive-expansion.md). Memories with write
ports are handled as described above.

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
