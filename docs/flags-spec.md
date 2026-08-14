# Full Flag Specification

This document captures the current command-line and YAML configuration
surface accepted by `kepler-formal`.

The stable top-level overview remains in [README.md](../README.md).
The SEC-specific flag surface is documented separately in
[sec-flags-spec.md](sec-flags-spec.md).

## Verification modes

`kepler-formal` supports full equivalence checking flows:

| Mode | Typical use | Inputs |
| --- | --- | --- |
| `lec` | Gate-level combinational equivalence checking | Verilog or Naja IF netlists plus Liberty/Python primitive libraries as needed. |
| `sec` | Gate-level sequential equivalence checking | Sequential Verilog/SystemVerilog netlists plus Liberty/Python primitive libraries as needed. |
| `sec` | RTL-level sequential equivalence checking | RTL Verilog/SystemVerilog sources, including SystemVerilog flists with explicit tops. |
| `sec` | SystemVerilog-to-Verilog RTL-vs-gate checking (`sv2v`) | SystemVerilog design 1 and Verilog design 2, plus Liberty/Python primitive libraries as needed. |

LEC is the default. Select SEC with `-v sec`, `--verification sec`, or
`verification: sec` in YAML.

## Binary flags

| Flag | Meaning |
| --- | --- |
| `--verification <lec\|sec>`, `-v <lec\|sec>` | Select LEC or SEC. Defaults to `lec`. |
| `--max-k <n>`, `-k <n>` | Set the SEC proof/search bound. Defaults to `32`; SEC only. |
| `--sec-engine <k_induction\|imc\|pdr>` | Select the SEC engine. Defaults to `pdr`; SEC only. |
| `--sec-encoding <binary\|dual_rail_steady>` | Select the SEC encoding. Defaults to `dual_rail_steady`; SEC only. |
| `--sec-uncomputable-seq-boundary` | Abstract unsupported sequential instances as SEC boundaries. This is the default. |
| `--no-sec-uncomputable-seq-boundary` | Fail SEC when an unsupported sequential instance is encountered. |
| `--allow-boundary-mismatch` | Allow LEC to continue when top-level inputs or sequential-element outputs do not match by name. Without this flag, a mismatch stops the run before SAT solving. LEC only. |
| `-verilog` | Use Verilog Format. |
| `-naja_if` | Use naja-if format. |
| `-systemverilog`, `-sv` | Use SystemVerilog format for both designs. Requires SEC verification. |
| `-sv2v` | Use mixed SystemVerilog-to-Verilog format for SEC RTL-vs-gate comparison: design 1 is parsed as SystemVerilog, design 2 is parsed as Verilog. |
| `-cc`, `-cxx`, `-cpp` | Synthesize one C/C++ translation unit per design to SystemVerilog through XLS, then run SEC on the generated RTL. |
| `--help`, `-h` | Print usage and exit. |
| `--config <file>`, `-c <file>` | Load a YAML config file. If present anywhere on the CLI, YAML parsing takes precedence over the rest of the arguments. |
| `--design1 <file...>` | Explicit source list for design 1 in multi-file Verilog mode. |
| `--design2 <file...>` | Explicit source list for design 2 in multi-file Verilog mode. |
| `--liberty <file...>`, `--lib <file...>` | Liberty library files. |
| `--verilog_preprocessing` | Enable preprocessing for Verilog inputs. |
| `--cc_top <function>` | C/C++ top function for both designs. |
| `--cc_design1_top <function>` | C/C++ top function for design 1 when it differs from `--cc_top`. |
| `--cc_design2_top <function>` | C/C++ top function for design 2 when it differs from `--cc_top`. |
| `--cc_module_name <name>` | Generated SystemVerilog module name for both designs. |
| `--cc_design1_module_name <name>` | Generated SystemVerilog module name for design 1. |
| `--cc_design2_module_name <name>` | Generated SystemVerilog module name for design 2. |
| `--cc_include <dir>`, `--cc_include_path <dir>`, `-I<dir>` | Add an include directory for XLS C/C++ synthesis. May be repeated. |
| `--cc_output_dir <dir>` | Directory for generated SystemVerilog. Defaults to `./kepler_formal_c2rtl`. |
| `--sv_design1_flist <file>`, `--sv_design2_flist <file>` | Per-design SystemVerilog file lists. Only design 1 is valid in `sv2v` mode. |
| `--sv_design1_top <top>`, `--sv_design2_top <top>` | Per-design SystemVerilog top modules. Only design 1 is valid in `sv2v` mode. |
| `--verilog_design1_top <top>`, `--verilog_design2_top <top>` | Per-design Verilog top modules. Only design 2 is valid in `sv2v` mode. |
| `--compact` | Reduce peak memory. In SEC, extract and release design 1 before loading design 2. |
| `--report-skipped-pos` | Emit skipped-PO reports in the current working directory. |

## YAML config flags

| Key | Type | Meaning |
| --- | --- | --- |
| `format` | string | Input format: `verilog`, `v`, `naja_if`, `systemverilog`, `sv`, `sv2v`, `cc`, `c`, `cxx`, `cpp`, or `c2rtl`. If omitted, the implementation defaults to `verilog`. |
| `verification` | string | `lec` or `sec`. Defaults to `lec`. |
| `max_k` | integer | SEC proof/search bound. Defaults to `32`. |
| `sec_engine` | string | `k_induction`, `imc`, or `pdr`. Defaults to `pdr`. |
| `sec_encoding` | string | `binary` or `dual_rail_steady`. Defaults to `dual_rail_steady`. |
| `sec_uncomputable_seq_as_boundary` | bool | Abstract unsupported sequential instances as SEC boundaries. Defaults to `true`. |
| `allow-boundary-mismatch` | bool | Allow an LEC boundary mismatch. Defaults to `false`; ignored for SEC. |
| `input_paths` | list | Required for normal runs. Accepts either `[design0, design1]` or `[[design0_file...], [design1_file...]]`. The nested form is for multi-file Verilog. |
| `liberty_files` | list[string] | Liberty libraries loaded through `SNLLibertyConstructor`. |
| `py_tech_files` | list[string] | Python primitive loaders loaded through `SNLPyLoader`. |
| `verilog_preprocessing` | bool | Enable preprocessing for Verilog inputs. |
| `cc_top` | string | C/C++ top function for both designs. Required for `cc` unless both design-specific tops are set. |
| `cc_design1_top` | string | C/C++ top function for design 1. |
| `cc_design2_top` | string | C/C++ top function for design 2. |
| `cc_module_name` | string | Generated SystemVerilog module name for both designs. Defaults to the resolved C/C++ top name. |
| `cc_design1_module_name` | string | Generated SystemVerilog module name for design 1. |
| `cc_design2_module_name` | string | Generated SystemVerilog module name for design 2. |
| `cc_block_proto_path` | string | Optional XLSCC block proto path for both designs. |
| `cc_design1_block_proto_path` | string | Optional XLSCC block proto path for design 1. |
| `cc_design2_block_proto_path` | string | Optional XLSCC block proto path for design 2. |
| `cc_include_paths` | list[string] | Include directories passed to XLS C/C++ synthesis. |
| `cc_output_dir` | string | Directory for generated SystemVerilog. Defaults to `./kepler_formal_c2rtl`. |
| `log_level` | string | `debug` and `info` are handled explicitly. Other values currently fall back to `info`. |
| `log_file` | string | Output path for the miter log file. If omitted, the tool writes `miter_log_<n>.txt` in the current working directory. |
| `use_scopes` | bool | Enable scoped verification for `naja_if` inputs. |
| `clean_scopes` | bool | Clean extracted scopes before scoped verification. Used with `use_scopes`. |
| `cnf_export` | bool | Enable CNF export. |
| `cnf_export_path` | string | Output path for CNF export. Defaults to `miter.cnf`, or `miter_<scope>.cnf` in scoped `naja_if` mode. |
| `po_cnf_export` | bool | Enable per-primary-output CNF export for each compared top. |
| `po_cnf_export_path` | string | Output directory for per-PO CNF export. Defaults to `po_cnfs`, or `po_cnfs_<scope>` in scoped `naja_if` mode. |
| `compact_mode` | bool | Same behavior as `--compact`. |
| `report_skipped_pos` | bool | Same behavior as `--report-skipped-pos`. |
| `sv_design1_flist`, `sv_design2_flist` | string | Per-design SystemVerilog file lists. Only design 1 is valid in `sv2v` mode. |
| `sv_design1_top`, `sv_design2_top` | string | Per-design SystemVerilog top modules. Only design 1 is valid in `sv2v` mode. |
| `verilog_design1_top`, `verilog_design2_top` | string | Per-design Verilog top modules. Only design 2 is valid in `sv2v` mode. |
| `solver` | string | SAT solver selection: `kissat`, `glucose`, or `cadical`. Defaults to `kissat`. |

Example:

```yaml
format: verilog
verification: lec
input_paths:
  - [design0_part1.v, design0_part2.v]
  - [design1_part1.v, design1_part2.v]
liberty_files:
  - library_file0.lib
  - library_file1.lib
py_tech_files:
  - primitives.py
verilog_preprocessing: true
solver: kissat
compact_mode: true
report_skipped_pos: true
cnf_export: true
cnf_export_path: ./miter.cnf
po_cnf_export: true
po_cnf_export_path: ./po_cnfs
```

Example C/C++ C2RTL config:

```yaml
format: cc
verification: sec
cc_top: top_function
cc_include_paths:
  - include
cc_output_dir: build/c2rtl
input_paths:
  - design0.cc
  - design1.cc
```

The `kepler-formal` CMake build compiles the KF C2RTL bridge as the normal
`kepler_xls_c2rtl` library target and links it into the binary:

```bash
cmake --build build --target kepler-formal
```

SEC `sv2v` example:

```yaml
format: sv2v
verification: sec
max_k: 32
sec_engine: pdr
sec_encoding: dual_rail_steady
sec_uncomputable_seq_as_boundary: true
input_paths:
  - [rtl_pkg.sv, rtl_top.sv]
  - [gate_top.v]
liberty_files:
  - stdcells.lib
solver: kissat
compact_mode: true
report_skipped_pos: true
```
