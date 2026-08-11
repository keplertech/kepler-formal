# Examples

The examples are grouped by design and technology so each directory can be run
independently.

| Directory | Contents |
| --- | --- |
| [`tinyrocket`](tinyrocket) | Verilog and Naja IF LEC/SEC examples, Liberty libraries, and input-generation scripts. |
| [`xilinx/register_slice`](xilinx/register_slice) | Small mapped-versus-compact SEC example using Python primitive models. |
| [`xilinx/vexriscv`](xilinx/vexriscv) | Large VexRiscv GenFull Xilinx LEC examples with equivalent and different designs. |
| [`xilinx/xilinx.py`](xilinx/xilinx.py) | Shared formal models for the Xilinx primitives used by both Xilinx examples. |

## TinyRocket

Run these commands from `examples/tinyrocket`:

```bash
# Verilog LEC
../../build/src/bin/kepler-formal -verilog \
  tinyrocket.v tinyrocket_edited.v \
  NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib \
  fakeram45_64x32.lib fakeram45_64x15.lib

# Naja IF LEC
../../build/src/bin/kepler-formal -naja_if \
  tinyrocket_naja.if tinyrocket_naja_edited.if \
  NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib \
  fakeram45_64x32.lib fakeram45_64x15.lib

# YAML examples
../../build/src/bin/kepler-formal --config test_config_verilog.yaml
../../build/src/bin/kepler-formal \
  --config test_config_verilog_tinyrocket_sec.yaml
```

The `to_naja_if.py`, `edit.py`, and `extract_tinyrocket_csrfile_sec.py`
scripts regenerate derived TinyRocket inputs used by tests and regressions.

## Xilinx Python Primitives

Run the register-slice SEC example from `examples/xilinx/register_slice`:

```bash
../../../build/src/bin/kepler-formal \
  --config test_config_verilog_xilinx_sec.yaml
```

Run the VexRiscv examples from `examples/xilinx/vexriscv`:

```bash
# Equivalent self-check
../../../build/src/bin/kepler-formal \
  --config test_config_verilog_vexriscv_genfull_xilinx_lec.yaml

# Intentional difference
../../../build/src/bin/kepler-formal \
  --config test_config_verilog_vexriscv_genfull_xilinx_lec_different.yaml
```

The shared `xilinx.py` file provides exact combinational truth tables,
instance-parameterized `LUT1` through `LUT6` models, and sequential models. The
VexRiscv self-check exercises a 4,995-cell generated Xilinx netlist containing
carry, DSP, distributed RAM, and block RAM macros.
