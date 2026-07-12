# SEC Flag Specification

This document describes the user-facing surface for Sequential Equivalence
Checking (SEC). SEC supports both gate-level sequential netlists and RTL-level
Verilog/SystemVerilog designs through the same extracted transition-system
model.

## Scope

This page covers the flags, YAML keys, reports, and debug environment variables
that affect Sequential Equivalence Checking (SEC). General input-format,
library, logging, solver, CNF export, and LEC flags remain documented in
[flags-spec.md](flags-spec.md).
SEC clock extraction and multi-clock-domain coverage handling are documented in
[sec-clock-handling.md](sec-clock-handling.md).

Supported SEC flows:

| Flow | Typical inputs |
| --- | --- |
| Gate-level SEC | Sequential gate-level Verilog/SystemVerilog netlists with Liberty/Python primitive libraries as needed. |
| RTL-level SEC | RTL Verilog/SystemVerilog sources, including SystemVerilog flists with explicit tops. |

## CLI Shape

SEC is selected with `-v sec` or `--verification sec`. LEC remains the default
verification mode when no SEC selector is provided.

```sh
kepler-formal -verilog \
  -v sec \
  -k 4 \
  --sec-engine pdr \
  --sec-encoding dual_rail_steady \
  --sec-uncomputable-seq-boundary \
  --report-skipped-pos \
  design0.v design1.v library.lib
```

The SEC selector flags are accepted before or after the input format flag:

```sh
kepler-formal -v sec --sec-engine k_induction --sec-encoding binary -k 8 -verilog design0.v design1.v
kepler-formal -verilog design0.v design1.v -v sec --sec-engine imc --sec-encoding dual_rail_steady -k 8
```

For SystemVerilog, SEC uses the normal SystemVerilog source options:

```sh
kepler-formal -sv \
  --sv_design1_flist design0.f \
  --sv_design1_top top0 \
  --sv_design2_flist design1.f \
  --sv_design2_top top1 \
  -v sec \
  --sec-engine pdr \
  --sec-encoding dual_rail_steady \
  -k 32
```

For RTL-vs-gate runs where design 1 is SystemVerilog and design 2 is Verilog,
use `-sv2v`:

```sh
kepler-formal -sv2v \
  --design1 rtl_pkg.sv rtl_top.sv \
  --design2 gate_top.v \
  -v sec \
  --sec-engine pdr \
  --sec-encoding dual_rail_steady
```

## YAML Shape

When `--config` or `-c` is present, YAML config mode takes precedence over the
rest of the command line.

```yaml
format: systemverilog
verification: sec
sec_engine: pdr
sec_encoding: dual_rail_steady
max_k: 32
sec_uncomputable_seq_as_boundary: true
compact_mode: true
report_skipped_pos: true
solver: kissat
sv_design1_flist: design0.f
sv_design1_top: cva6
sv_design2_flist: design1.f
sv_design2_top: cva6
liberty_files:
  - NangateOpenCellLibrary_typical.lib
```

## SEC Flags And Keys

| CLI flag | YAML key | Default | Values | Effect |
| --- | --- | --- | --- | --- |
| `-v sec`, `--verification sec` | `verification: sec` | `lec` | `lec`, `sec` | Selects SEC instead of combinational LEC. Values are lowercase. |
| `-k <n>`, `--max-k <n>` | `max_k: <n>` | `32` | Non-negative integer | Sets the maximum SEC proof/search bound used by the selected engine. `0` is valid and only permits zero-bound checks. |
| `--sec-engine <engine>` | `sec_engine: <engine>` | `pdr` | `k_induction`, `imc`, `pdr` | Selects the top-level SEC proof engine. Engine names are lowercase. |
| `--sec-encoding <mode>` | `sec_encoding: <mode>` | `dual_rail_steady` | `binary`, `dual_rail_steady` | Selects how SEC models unknown or reset-unanchored state values. Omit the key/flag to use the dual-rail default. |
| `--sec-uncomputable-seq-boundary` | `sec_uncomputable_seq_as_boundary: true` | `true` | boolean | Abstracts unsupported sequential instances as SEC boundaries instead of failing immediately. |
| `--no-sec-uncomputable-seq-boundary` | `sec_uncomputable_seq_as_boundary: false` | `true` | boolean | Uses strict mode: unsupported sequential interfaces cause SEC to fail as unsupported. |
| `--compact` | `compact_mode: true` | `false` | boolean | Enables compact SEC extraction: design 1 is extracted and released before design 2 is loaded; identical SEC inputs can reuse the extracted design 1 model. |
| `--report-skipped-pos` | `report_skipped_pos: true` | `false` | boolean | Enables skipped-output reporting and writes SEC boundary reporting when entries exist. |

Accepted values for `sec_engine`:

| Engine | Accepted value |
| --- | --- |
| `k_induction` | `k_induction` |
| `imc` | `imc` |
| `pdr` | `pdr` |

Accepted values for `sec_encoding`:

| Encoding | Accepted value | Meaning |
| --- | --- | --- |
| Dual-rail steady-state mode | `dual_rail_steady` | Default SEC mode. Tracks value/knownness rails so outputs driven by resetless or reset-unanchored state can be represented without assuming cross-design internal flop equality. |
| Binary mode | `binary` | Historical concrete 0/1 SEC mode. Outputs whose cones depend on reset-unanchored internal state can be skipped because SEC cannot assume equality between internal flops in different designs. |

There is no `default` token for `sec_encoding`. To use the default, omit
`--sec-encoding` or omit `sec_encoding` from YAML. Tests, scripts, and regression
flows that require stable behavior should always spell out either `binary` or
`dual_rail_steady` explicitly.

## Engine Semantics

| Engine | Current behavior |
| --- | --- |
| `k_induction` | Explicit classic k-induction flow: bounded base-case search followed by induction-step proof over the extracted SEC transition system. |
| `imc` | Interpolation-Based Model Checking flow over the same extracted SEC problem. It uses the shared base-case search and exact interpolant strengthening where applicable. |
| `pdr` | Property Directed Reachability flow over the extracted SEC transition system. It first accepts immediate zero-bound k-induction results, then runs PDR frames up to `max_k`. |

All engines use the same extracted SEC model: aligned environment inputs,
state bits, observed outputs, next-state formulas, initial-state information,
and complemented-state relations.

SEC only aligns environment inputs and observed top outputs across the two
designs. It must not use matching internal instance, net, or flop names as a
cross-design equivalence assumption. Internal names may still be used inside a
single extracted design for diagnostics, state updates, and local recovery
heuristics.

## Bounds And Results

`max_k` is parsed as a non-negative integer.

SEC result handling is currently:

| Result | Exit code | Meaning |
| --- | --- | --- |
| Equivalent | `0` | SEC proved equivalence at the reported bound. |
| Different | `0` | SEC found a concrete counterexample at the reported bound. |
| Inconclusive | non-zero | The selected engine reached `max_k` without a proof or counterexample. |
| Unsupported | non-zero | The extracted model was incomplete or unsupported for SEC. |

The log always prints:

```text
SEC max_k: <n>
SEC engine: <engine>
SEC encoding: binary|dual_rail_steady
SEC uncomputable sequentials: boundary abstraction|strict failure
Compact mode: enabled|disabled
Skipped PO reports: enabled|disabled
```

## Boundary And Skipped Reports

When `--report-skipped-pos` or `report_skipped_pos: true` is enabled, SEC may
write the following files in the current working directory:

| File | Producer | Contents |
| --- | --- | --- |
| `boundary_terms.txt` | SEC | Extracted SEC boundary surface. Includes top inputs, top outputs, opaque internal cut points, abstracted sequential state terms, abstracted sequential observed terms, and connectivity-skip annotations when present. |
| `skipped_no_driver_pos.txt` | shared cone builder | Outputs skipped because the relevant iso has no driver. |
| `skipped_multi_driver_pos.txt` | shared cone builder | Outputs skipped because the relevant iso has multiple drivers. |
| `skipped_logical_loop_pos.txt` | shared cone builder | Outputs skipped because the relevant cone contains a logical loop. |
| `skipped_reset_unanchored_pos.txt` | SEC | Outputs skipped in binary SEC because their cones depend on reset-unanchored internal state. |
| `skipped_multi_clock_domain_pos.txt` | SEC clock model | Outputs skipped because the observed output cone spans multiple extracted clock domains. |

`boundary_terms.txt` starts with a category legend. Current categories are:

| Role | Meaning |
| --- | --- |
| `top_input` | Original top-level input term. |
| `top_output` | Original top-level output term. |
| `opaque_internal_input` | Internal cut-point input that SEC could not reconstruct combinationally and did not model as sequential. |
| `opaque_internal_output` | Internal cut-point output paired with an opaque internal boundary. |
| `abstracted_sequential_state` | State-facing term exposed when an uncomputable sequential instance is abstracted as a SEC boundary. |
| `abstracted_sequential_observed` | Observed-output-facing term exposed when an uncomputable sequential instance is abstracted as a SEC boundary. |

Connectivity skipped outputs are also summarized in the main run log, including
no-driver, multi-driver, logical-loop, reset-unanchored, and multi-clock-domain
skips.

## Compact SEC

Compact SEC is intended for large designs where holding both elaborated designs
at the same time can be too expensive.

Current behavior:

1. Load and extract design 1.
2. Release design 1's database.
3. If design 2 has the same SEC input specification, reuse the immutable design
   1 extracted model for design 2.
4. Otherwise load and extract design 2, release it, then run the proof over the
   two extracted models.

Compact SEC can reduce peak memory. The tradeoff is that counterexample cone
traceback may be less detailed because the elaborated design database can be
released before witness formatting.

## Diagnostic Environment Variables

These are debug-only controls, not stable public flags.

| Environment variable | Effect |
| --- | --- |
| `KEPLER_SEC_DIAG=1` | Prints SEC extraction, alignment, boundary, proof-problem, and engine-progress diagnostics to stderr/log output. |
| `KEPLER_SEC_PDR_TRACE=1` | Prints PDR-specific problem and frame traces. Most useful together with `--sec-engine pdr`. |

## Related Non-SEC Inputs

SEC still depends on the normal front-end and library flags:

| CLI/YAML | Why it matters for SEC |
| --- | --- |
| `-verilog`, `format: verilog` | Verilog SEC input mode. |
| `-sv`, `-systemverilog`, `format: systemverilog` | SystemVerilog SEC input mode. |
| `-sv2v`, `format: sv2v` | Mixed SEC input mode: design 1 is SystemVerilog RTL and design 2 is Verilog gate-level netlist. |
| `--design1`, `--design2`, `input_paths` | Source lists for the two compared designs. |
| `--sv_design1_flist`, `--sv_design2_flist` | Per-design SystemVerilog file lists. In `sv2v` mode, only `--sv_design1_flist` is accepted. |
| `--sv_design1_top`, `--sv_design2_top` | Per-design SystemVerilog top names. In `sv2v` mode, only `--sv_design1_top` is accepted. |
| `--liberty`, `--lib`, `liberty_files` | Liberty primitives, including structured memory information used during SEC extraction. |
| `py_tech_files` | Python primitive loaders. Must be provided through YAML, not through `--liberty`. |
| `solver: kissat|glucose` | SAT solver used by the selected SEC engine. |
| `log_file` | Run log path. SEC defaults to `miter_log_<n>.txt` when no path is provided. |

## Known Construction Notes

- `pdr` is the default SEC engine. Tests, scripts, and regression flows can
  still choose an explicit engine (`k_induction`, `imc`, or `pdr`) when
  comparing behavior.
- `dual_rail_steady` is the default SEC encoding. New tests and regressions should
  still set `sec_encoding` explicitly when they depend on a specific mode.
- `report_skipped_pos` still has historical LEC naming even though SEC also uses
  it to gate `boundary_terms.txt`.
- The exact contents of SEC diagnostics and boundary reports are expected to
  evolve as memory modeling and boundary extraction mature.
