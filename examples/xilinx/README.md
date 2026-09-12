# Xilinx Examples

The Xilinx examples use the shared [`xilinx.py`](xilinx.py) formal primitive
models through the YAML `py_tech_files` option.

| Directory | Contents |
| --- | --- |
| [`register_slice`](register_slice) | Small mapped-versus-compact SEC equivalence example. |
| [`vexriscv`](vexriscv) | Large VexRiscv GenFull LEC self-check and intentional-difference examples. |

`xilinx.py` models the combinational, parameterized LUT, sequential, carry,
DSP, distributed RAM, and block RAM primitives used by these netlists.
