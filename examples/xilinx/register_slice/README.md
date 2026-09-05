# Xilinx Register-Slice SEC

This compares a Xilinx-mapped register slice with a compact equivalent model.
SEC is expected to prove equivalence.

Run from `examples/xilinx/register_slice`:

```bash
../../../build/src/bin/kepler-formal \
  --config test_config_verilog_xilinx_sec.yaml
```
