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
`kepler-formal` selects the internal translator, performs the C-to-RTL
conversion in-process, and then loads the generated RTL through the existing RTL
flow. No frontend command is exposed to users.

Current backends:

- `metron`: working backend for the existing Metron-style benchmark subset.
- `xls`: selected by `frontend: xls` or automatically for XLS/xlscc-looking C++;
  `thirdparty/xls` is vendored from `https://github.com/keplertech/xls` and is
  linked through the in-process XLS wrapper when the XLS frontend is enabled.

For XLS block/channel examples, a hidden per-design `block_proto` key can point
at the xlscc HLS block textproto:

```yaml
design1:
  format: c
  frontend: xls
  input_paths:
    - source.cc
  top: Mux3_f
  block_proto: mux3.textproto
  module_name: mux3_comb
```

`module_name` is also used as the generated RTL top when it is present. This is
useful for XLS block examples where the C top function and emitted Verilog
module name naturally differ.

## Artifacts

By default, generated frontend artifacts are created under a temporary directory
and removed at the end of the run. Set `keep_generated: true` to preserve them.

Each preserved C frontend work directory contains:

- `designN_<top>_from_c.sv`: generated SystemVerilog passed to the RTL loader.
- `input_manifest.json`: source paths, top, clock/reset metadata, internal
  translation command description,
  and log paths.
- `c_frontend_stdout.log`
- `c_frontend_stderr.log`

## Current Contract

This is not arbitrary software C. The public contract is a synthesizable Kepler C
subset with explicit hardware semantics. The top name is required because the
generated RTL is compared through the same top input/output contract as normal
SystemVerilog SEC.

The intended XLS path is:

```text
C/C++ source
  -> in-process xlscc::Translator
  -> XLS IR package
  -> in-process XLS optimization
  -> in-process XLS Verilog codegen
  -> generated SystemVerilog artifact
  -> existing Kepler RTL loader and LEC/SEC flow
```

## Enabling XLS

The default Kepler build keeps the XLS frontend disabled so ordinary RTL and
Metron C2RTL builds do not need the full XLS dependency stack. To enable the
embedded XLS path and let CMake build the fork-local XLS wrapper with Bazel:

```sh
cmake -S . -B build-xls \
  -DKEPLER_BUILD_XLS_FRONTEND=ON
cmake --build build-xls --target kepler-formal
```

If the XLS wrapper was already built separately, point CMake at it instead:

```sh
cmake -S . -B build-xls \
  -DKEPLER_ENABLE_XLS_FRONTEND=ON \
  -DKEPLER_XLS_C2RTL_LIBRARY=thirdparty/xls/bazel-bin/xls/contrib/kepler/libkepler_xls_c2rtl.dylib
```

Use `libkepler_xls_c2rtl.so` for Linux builds.
