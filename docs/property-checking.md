# Property Checking

Kepler Formal can prove the `assert`, `assume` and `cover` statements written
in a SystemVerilog design. This mode checks one design against its own
properties; it does not compare two designs.

```text
SystemVerilog + assertions → frontend keeps each assertion and the signals it reads
                           → the assertion is translated into logic
                           → one monitor per property → PDR, k-induction or IMC
                           → proved / failed with a trace / covered / ...
```

## CLI

```sh
kepler-formal -sv -v fpv \
  --sec-reset-cycles 2 --sec-reset-port rst_n=0 \
  design_pkg.sv design.sv
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `-v fpv`, `--verification fpv` | | Selects property checking. Only SystemVerilog input is supported. |
| `--top <name>` | the single top module | Top module. |
| `-f <file>` | | SystemVerilog file list, passed to the frontend. |
| `-k <n>`, `--max-k <n>` | `32` | Proof and search bound of the engine. |
| `--sec-engine <engine>` | `pdr` | `pdr`, `k_induction` or `imc`. |
| `--sec-reset-cycles <n>` | none | Number of cycles the reset ports are held active. |
| `--sec-reset-port <name=0\|1>` | none | A reset input and its active value. Repeatable. |

The reset flags are the ones of SEC, with the same syntax. Other flags are not
accepted with `-v fpv`, and the mode is not available yet from a YAML
configuration or from the Python API.

Without reset flags, `imc` finds failures but leaves every other property
`inconclusive`; use `pdr` or `k_induction` there.

## Supported properties

A property is supported when it is a safety property on one clock. Anything
else is reported as `unsupported` with the reason; it is never dropped and
never counted as proved.

| Part | Supported |
| --- | --- |
| Statements | `assert property`, `assume property` and `cover property` written as module items. The label is the property name; without a label the name is the kind and the line, for example `assert_line42`. Action blocks are ignored. |
| Clock | `@(posedge clk)` on the design clock, written in the property or taken from `default clocking`. |
| Disable | `disable iff (expr)` and `default disable iff`. |
| Expressions | Packed signals of the module, constants and parameters; logical, bitwise and reduction operators; `==`, `!=`, `<`, `<=`, `>`, `>=`; `+`, `-`; shifts; `?:`; `inside`; concatenation and replication; bit selects, part selects with constant bounds and packed structure fields; `$past(x)`, `$past(x, n)`, `$rose`, `$fell`, `$stable`, `$changed`, `$onehot`, `$onehot0`, `$countones`, `$isunknown`, `$sampled`, `$signed`, `$unsigned`. Widths and signs follow the SystemVerilog rules. |
| Sequences on the left of `\|->` or `\|=>` | Boolean, `##N`, `##[M:N]`, `##[M:$]`, `[*N]`, `[*M:N]`, `[*M:$]`, `[->N]`, `[->M:N]`, `or`. |
| Sequences on the right | Boolean, `##N`, `[*N]`, and one final delay window `##[M:N] expr`. |
| Properties | Boolean, `seq \|-> prop`, `seq \|=> prop`, `not seq`, `and`, `if`/`else`, `nexttime`, `nexttime [N]`, `always`, weak `until` and `until_with` on Boolean operands. Named `property` and `sequence` declarations with arguments. |
| Covers | Sequences. |

Not supported yet: liveness (`s_eventually`, `s_until`, `strong`), negedge
clocks, clocks with `iff`, more than one clock in a property, local
variables, `intersect`, `within`, `throughout`, `first_match`, nonconsecutive
repetition, a delay range that is not the last element of a consequent,
`accept_on` and `reject_on`, covers of implications, and in expressions:
multiplication, division, modulo and power of signals, part selects with a
variable bound, an index that can fall outside its vector, and constants with
unknown bits.

The frontend stops the whole run, with the source position, on a form it
cannot keep: a concurrent assertion inside procedural code, `restrict
property`, `cover sequence`, a call of a user function, and a reference to
anything but a packed signal of the same module (a hierarchical reference, an
unpacked array). Immediate assertions are ignored.

## Model

- **One step is one rising edge of the design clock.** Signals in a property
  have the value they hold just before that edge.
- **Reset.** With `--sec-reset-cycles N`, the listed ports are held active for
  the first `N` steps and inactive afterwards. Properties are checked from the
  first step after that window. During the window, properties and assumptions
  make no claim and their pending attempts are dropped.
- **Initial state.** A register with an initial value starts there. Every other
  register starts from an arbitrary value, before the reset window. Without a
  reset, properties are therefore checked from every possible state.
- **Nothing is assumed about any step.** In particular, the first step after
  reset is checked like every other one.
- **Values are two-valued.** An unknown value is an arbitrary 0 or 1,
  `!==`/`===` behave as `!=`/`==`, and `$isunknown` is always 0.
- **History.** `$past` and the other sampled-value functions keep their
  history during reset. Before the first step the history is arbitrary.
- **Assumptions** hold in every step of a trace, up to and including the step
  where an assertion fails or a cover is hit. A trace that breaks an assumption
  at any earlier step is not a counterexample.
- **`disable iff`** removes the claim in that step and drops the attempts in
  progress.

## Results

Each property gets one line:

```text
Property a_range [assert] counter.sv:9: proved (k = 3)
Property a_next_bad [assert] counter.sv:17: failed at cycle 2
Input trace:
  reset cycle 0: en[0]=1, initial cnt[0]=0, initial cnt[1]=1, initial cnt[2]=1
  cycle 0: en[0]=1
  cycle 1: en[0]=1
  cycle 2: en[0]=0
Property c_limit [cover] counter.sv:19: covered at cycle 5
Property u_fifo/a_live [assert] fifo.sv:40: unsupported: `s_eventually` is a liveness obligation, which is not supported yet
```

| Status | Meaning |
| --- | --- |
| `proved` | The assertion holds in every reachable step. |
| `failed` | A trace violates the assertion. |
| `covered` | A trace reaches the cover. |
| `unreachable` | The cover is proved impossible. |
| `assumed` | The assumption was applied to every other property. |
| `inconclusive` | Neither a proof nor a trace within the bound. |
| `unsupported` | The property is outside the supported subset. |

Cycle numbers count from the first cycle after reset. A trace lists the input
values of every cycle, the reset cycles first, and in its first line the
initial value of every register that had none. Reset ports and clocks are not
listed.

If an assumption is unsupported it is not applied. Proofs stay valid, since
they then cover more behaviours; a `failed` or `covered` trace carries a note
that it may break that assumption.

| Exit code | Meaning |
| --- | --- |
| `0` | Every assertion is proved and every cover is covered. |
| `1` | No failure; some properties are proved or covered, others inconclusive or unsupported. |
| `2` | Nothing was established: no properties, only inconclusive or unsupported ones, or the design could not be modelled. |
| `3` | An assertion failed or a cover is unreachable. |

## How it works

1. The SystemVerilog frontend is run with `preserveAssertions` and
   `keepASTLink`. Each assertion becomes one input-only
   `naja_formal_property` instance. Its `DATA` pins are connected to the
   signals the assertion reads, and it records the source text and where each
   signal sits in `DATA`. The frontend builds no logic for the assertion, and
   the instance drives nothing, so the design is unchanged.
2. These instances are selected as verification boundaries, which makes the
   extracted model publish the logic of every `DATA` bit.
3. `SvaPropertyTranslator` reads the elaborated assertion the frontend kept
   for the instance. It compiles the clock, the disable condition and every
   Boolean operand into logic over the `DATA` bits, and writes the temporal
   structure as a prefix expression, for example
   `(implies_next (leaf 0) (leaf 1))`. Each `$past` adds history registers.
4. `PropertyCompiler` turns that expression into a 1-bit violation (or match)
   signal plus monitor state. A sequence on the left of an implication is
   matched with one state bit per delay step, shared by all attempts. A
   sequence on the right is followed step by step; a final delay window fails
   at its deadline when the operand was false throughout.
5. The reset sequence and the free initial values are added to the model as
   ordinary logic, so every state bit has an initial value.
6. Each assertion and each cover is proved separately by the selected engine,
   as "this signal is always 1".

Set `KEPLER_FPV_DIAG=1` to print the transcribed expression and the monitor
size of every property.
