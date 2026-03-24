# SystemVerilog Support

This document tracks the current experimental SystemVerilog flow in `kepler-formal`.

Status:

- ongoing work
- experimental
- subject to CLI and YAML format changes
- not part of the stable interface documented in the top-level [README](../../README.md)

## Current scope

The SystemVerilog path is intended for development and validation work while the interface and behavior are still being refined.

Current entry points include:

- CLI:
  - `-systemverilog`
  - `-sv`
- YAML `format`:
  - `systemverilog`
  - `sv`

## Current CLI forms

```bash
# Classic (single file per design)
build/src/bin/kepler-formal <-systemverilog/-sv> [--verilog_preprocessing] <netlist1> <netlist2> [<library-file>...]

# Multi-file SystemVerilog designs
build/src/bin/kepler-formal <-systemverilog/-sv> [--verilog_preprocessing] --design1 <file...> --design2 <file...> \
  [--liberty <library-file>...]

# slang flists with explicit tops
build/src/bin/kepler-formal -systemverilog \
  --sv_design1_flist <file> --sv_design1_top <name> \
  --sv_design2_flist <file> --sv_design2_top <name>
```

`--verilog_preprocessing` is also accepted as `--verilog-preprocessing`.

## Flist mode

For SystemVerilog designs that are already driven by a slang command file or flist, use:

- `sv_design1_flist`
- `sv_design2_flist`
- `sv_design1_top`
- `sv_design2_top`

This mode is flist-based. Do not combine it with:

- `input_paths`
- `--design1`
- `--design2`

## YAML examples

Multi-file SystemVerilog example:

```yaml
format: systemverilog
input_paths:
  - [design0_pkg.sv, design0_top.sv]
  - [design1_pkg.sv, design1_top.sv]
liberty_files:
  - stdcells.lib.gz
  - primitives.py
```

Flist example:

```yaml
format: systemverilog
sv_design1_flist: /path/to/design1.f
sv_design1_top: top1
sv_design2_flist: /path/to/design2.f
sv_design2_top: top2
liberty_files:
  - stdcells.lib.gz
```

## Notes

- Expect behavior, validation rules, and documentation to evolve.
- If you are documenting stable usage for new users, prefer the top-level [README](../../README.md).
