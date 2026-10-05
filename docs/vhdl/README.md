# VHDL Support

This document tracks the VHDL flow in `kepler-formal`.

Status:

- experimental: the VHDL frontend in Naja is in beta and accepts a restricted
  RTL subset
- supported for RTL-level SEC, with both designs given as VHDL
- the `vhdl` input mode requires SEC verification

## Usage

CLI flag: `-vhdl`. YAML: `format: vhdl` (or `vhd`).

```bash
# Single file per design
build/src/bin/kepler-formal -vhdl -v sec <design1.vhd> <design2.vhd>

# Multi-file VHDL designs
build/src/bin/kepler-formal -vhdl -v sec \
  --design1 <file...> --design2 <file...> \
  [--vhdl_design1_top <top>] [--vhdl_design2_top <top>]
```

```yaml
format: vhdl
verification: sec
vhdl_design1_top: top
vhdl_design2_top: top
input_paths:
  - [design1/pkg.vhd, design1/leaf.vhd, design1/top.vhd]
  - [design2/top.vhd]
```

## File order

Files are loaded one at a time, in the order given, and each file can use the
units of the files before it. List them in compile order: packages and
instantiated entities first, the top-level unit last. A file that instantiates
an entity from a later file is rejected with a `missing entity` error.

## Top selection

`vhdl_design1_top` and `vhdl_design2_top` name the top entity of each design.
The top is elaborated after the last file is loaded, so it may be declared in
any of the files.

Without a top option, the top is the design built from the last file that
produces one; package-only files produce none. When that file declares several
entities, the top is the one no other entity instantiates. Name the top
explicitly whenever the last file is not the top-level unit.

## Supported language subset

Naja lowers VHDL to the same primitives as its SystemVerilog frontend, as
two-state hardware. The subset covers `bit`, `bit_vector` and imported
`std_logic` types, concurrent and conditional assignments, positive-edge
processes with synchronous reset and enable, constrained arrays with static
indexing, and direct entity instantiation. See
`thirdparty/naja/src/vhdl/README.md` for the current limits.

A construct outside the subset stops the run with an error that names the file
and line. It is never approximated.

## Notes

- Warnings from the VHDL frontend are printed on the console. No separate
  diagnostics report is written.
- Registers without a reset bootstrap need the default `dual_rail_steady`
  encoding. With `sec_encoding: binary`, give a reset bootstrap through
  `sec_reset`; see [SEC reset bootstrap](../sec-reset-bootstrap.md).
- Mixed comparisons, such as VHDL against Verilog or SystemVerilog, are not
  available.
