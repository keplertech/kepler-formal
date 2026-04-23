# Full Flag Specification

This document captures the current command-line and YAML configuration
surface accepted by `kepler-formal`.

The stable top-level overview remains in [README.md](../README.md).

## Binary flags

### Stable formats

- `-verilog`
- `-naja_if`

### Stable invocation forms

```bash
# YAML config
build/src/bin/kepler-formal --config <file.yaml>

# Single file per design
build/src/bin/kepler-formal <-verilog/-naja_if> [options] <design1> <design2> [<library-file>...]

# Multi-file Verilog
build/src/bin/kepler-formal -verilog [options] --design1 <file...> --design2 <file...> \
  [--liberty <library-file>...]
```

### Stable binary flags

| Flag | Meaning |
| --- | --- |
| `--help`, `-h` | Print usage and exit. |
| `--config <file>`, `-c <file>` | Load a YAML config file. If present anywhere on the CLI, YAML parsing takes precedence over the rest of the arguments. |
| `--design1 <file...>` | Explicit source list for design 1 in multi-file Verilog mode. |
| `--design2 <file...>` | Explicit source list for design 2 in multi-file Verilog mode. |
| `--liberty <file...>`, `--lib <file...>` | Liberty library files. |
| `--verilog_preprocessing` | Enable preprocessing for Verilog inputs. |
| `--compact` | Per-PO analysis is skipped in case the design is different. |
| `--report-skipped-pos` | Emit skipped-PO reports in the current working directory. |

### Stable binary behavior

- `-naja_if` requires exactly one input file per design.
- In the single-file positional form, extra positional arguments after the
  first two design inputs are treated as Liberty files.
- `--design1` and `--design2` switch the parser into explicit multi-file
  mode.
- `--liberty` and `--lib` load files through `SNLLibertyConstructor`.
- Python loader files (`.py`) are rejected on the `--liberty` / `--lib`
  path. Use YAML `py_tech_files` instead.
- `--report-skipped-pos` may produce `skipped_multi_driver_pos.txt`,
  `skipped_no_driver_pos.txt`, and `skipped_logical_loop_pos.txt`.

## YAML config flags

### Stable format values

- `verilog`
- `v`
- `naja_if`

### Stable config keys

| Key | Type | Meaning |
| --- | --- | --- |
| `format` | string | Input format. If omitted, the implementation defaults to `verilog`. |
| `input_paths` | list | Required for normal runs. Accepts either `[design0, design1]` or `[[design0_file...], [design1_file...]]`. The nested form is for multi-file Verilog. |
| `liberty_files` | list[string] | Liberty libraries loaded through `SNLLibertyConstructor`. |
| `py_tech_files` | list[string] | Python primitive loaders loaded through `SNLPyLoader`. |
| `verilog_preprocessing` | bool | Enable preprocessing for Verilog inputs. |
| `log_level` | string | `debug` and `info` are handled explicitly. Other values currently fall back to `info`. |
| `log_file` | string | Output path for the miter log file. If omitted, the tool writes `miter_log_<n>.txt` in the current working directory. |
| `use_scopes` | bool | Enable scoped verification for `naja_if` inputs. |
| `clean_scopes` | bool | Clean extracted scopes before scoped verification. Used with `use_scopes`. |
| `cnf_export` | bool | Enable CNF export. |
| `cnf_export_path` | string | Output path for CNF export. Defaults to `miter.cnf`, or `miter_<scope>.cnf` in scoped `naja_if` mode. |
| `compact_mode` | bool | Same behavior as `--compact`. |
| `report_skipped_pos` | bool | Same behavior as `--report-skipped-pos`. |
| `solver` | string | SAT solver selection. Supported values: `kissat`, `glucose`. If omitted, the implementation defaults to `kissat`. |

### Stable YAML behavior

- Unknown YAML keys are rejected.
- `input_paths` must contain exactly two design entries.
- In nested-list form, each design entry must contain at least one file.
- `use_scopes` and `clean_scopes` are meaningful only for `naja_if`.
- `liberty_files` entries are loaded as Liberty without suffix checks.
- `py_tech_files` is YAML-only.

Example:

```yaml
format: verilog
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
```

## Under development

SystemVerilog support is accepted by the parser today but remains under
development and should not be treated as a stable public interface.

### SystemVerilog format values

- CLI: `-systemverilog`, `-sv`
- YAML `format`: `systemverilog`, `sv`

### SystemVerilog binary forms

```bash
# Single file per design
build/src/bin/kepler-formal <-systemverilog/-sv> [options] <design1> <design2> [<library-file>...]

# Multi-file source lists
build/src/bin/kepler-formal <-systemverilog/-sv> [options] --design1 <file...> --design2 <file...> \
  [--liberty <library-file>...]

# Flist mode with explicit tops
build/src/bin/kepler-formal -systemverilog \
  --sv_design1_flist <file> --sv_design1_top <name> \
  --sv_design2_flist <file> --sv_design2_top <name> \
  [--design1 <file...>] [--design2 <file...>] [--liberty <library-file>...]
```

### SystemVerilog binary flags

| Flag | Meaning |
| --- | --- |
| `--sv_design1_flist <file>` | Per-design flist / command file for design 1. |
| `--sv_design2_flist <file>` | Per-design flist / command file for design 2. |
| `--sv_design1_top <name>` | Explicit top name for design 1. |
| `--sv_design2_top <name>` | Explicit top name for design 2. |

### SystemVerilog YAML keys

| Key | Type | Meaning |
| --- | --- | --- |
| `sv_design1_flist` | string | Per-design flist / command file for design 1. |
| `sv_design2_flist` | string | Per-design flist / command file for design 2. |
| `sv_design1_top` | string | Explicit top name for design 1. |
| `sv_design2_top` | string | Explicit top name for design 2. |

### SystemVerilog behavior notes

- These flags and keys are valid only with SystemVerilog input mode.
- Each design must provide at least one source, a flist, or both.
- Empty `sv_design*_flist` and `sv_design*_top` values are rejected.
- When a top name is provided, the implementation builds a temporary
  command file and invokes the parser through `-f`.
