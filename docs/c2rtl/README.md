# C-to-RTL SEC Frontend

Status: experimental.

The C frontend lets a SEC run compare a C model against RTL by translating the C
side to temporary SystemVerilog and then using the normal SystemVerilog SEC
flow.

```text
C source
  -> Kepler C frontend
  -> internal C-to-RTL translator
  -> generated SystemVerilog
  -> existing SystemVerilog loader
  -> existing SEC engine
```

## YAML Shape

Use per-design YAML sections when one side is C:

```yaml
verification: sec
sec_engine: pdr
sec_encoding: dual_rail_steady
max_k: 32

design1:
  format: c
  input_paths:
    - model/foo.c
  top: foo
  clock: clk
  reset: rst_n

design2:
  format: systemverilog
  input_paths:
    - rtl/foo.sv
  top: foo
```

`format: c` is user-facing. The translator backend is an implementation detail.
For local development, the executable can be overridden with
`KEPLER_C_FRONTEND_TRANSLATOR` or `KEPLER_C2RTL_TRANSLATOR`; otherwise
`kepler-formal` looks for `metron` on `PATH`.

## Artifacts

By default, generated frontend artifacts are created under a temporary directory
and removed at the end of the run. Set `keep_generated: true` to preserve them.

Each preserved C frontend work directory contains:

- `designN_<top>_from_c.sv`: generated SystemVerilog passed to the RTL loader.
- `input_manifest.json`: source paths, top, clock/reset metadata, command line,
  and log paths.
- `c_frontend_stdout.log`
- `c_frontend_stderr.log`

## Current Contract

This is not arbitrary software C. The public contract is a synthesizable Kepler C
subset with explicit hardware semantics. The top name is required because the
generated RTL is compared through the same top input/output contract as normal
SystemVerilog SEC.
