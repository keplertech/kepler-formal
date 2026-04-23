# Full Flag Specification

This document captures the current command-line and YAML configuration
surface accepted by `kepler-formal`.

The stable top-level overview remains in [README.md](../README.md).

## Binary flags

| Flag | Meaning |
| --- | --- |
| `-verilog` | Use Verilog Format. |
| `-naja_if` | Use naja-if format. |
| `--help`, `-h` | Print usage and exit. |
| `--config <file>`, `-c <file>` | Load a YAML config file. If present anywhere on the CLI, YAML parsing takes precedence over the rest of the arguments. |
| `--design1 <file...>` | Explicit source list for design 1 in multi-file Verilog mode. |
| `--design2 <file...>` | Explicit source list for design 2 in multi-file Verilog mode. |
| `--liberty <file...>`, `--lib <file...>` | Liberty library files. |
| `--verilog_preprocessing` | Enable preprocessing for Verilog inputs. |
| `--compact` | Per-PO analysis is skipped in case the design is different. |
| `--report-skipped-pos` | Emit skipped-PO reports in the current working directory. |

## YAML config flags

| Key | Type | Meaning |
| --- | --- | --- |
| `format` | string | Input format[verilog, naja_if]. If omitted, the implementation defaults to `verilog`. |
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
