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

All files of a design are registered before anything is elaborated, so they
may be listed in any order. An entity instantiated from another file resolves
from the complete unit index, and a child entity takes the generic actuals of
its parent's `generic map` rather than its own defaults, so a child whose
defaults alone would not elaborate (a null array, say) is fine. Files that are
listed twice are read once.

## Top selection

`vhdl_design1_top` and `vhdl_design2_top` name the top entity of each design.
The top is elaborated once all files are known, so it may be declared in any
of them.

Without a top option, the top is the one entity no other entity instantiates
across the design's files. Name the top explicitly when several entities are
left uninstantiated.

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
