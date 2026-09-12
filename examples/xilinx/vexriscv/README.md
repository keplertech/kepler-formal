# Xilinx VexRiscv LEC

Run these commands from `examples/xilinx/vexriscv`.

The self-check compares the 4,995-cell VexRiscv GenFull netlist with itself and
is expected to report no difference.

```bash
../../../build/src/bin/kepler-formal \
  --config test_config_verilog_vexriscv_genfull_xilinx_lec.yaml
```

The negative case compares it with an intentionally changed netlist and is
expected to report a difference.

```bash
../../../build/src/bin/kepler-formal \
  --config test_config_verilog_vexriscv_genfull_xilinx_lec_different.yaml
```
