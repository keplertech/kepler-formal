# Examples

## Instructions

```bash
cd example
pip install najaeda
python edit.py
# For naja_if
../build/src/bin/kepler-formal -naja_if tinyrocket_naja.if tinyrocket_naja_edited.if \
  NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib fakeram45_64x32.lib fakeram45_64x15.lib
# For verilog
../build/src/bin/kepler-formal -verilog tinyrocket.v tinyrocket_edited.v \
  NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib fakeram45_64x32.lib fakeram45_64x15.lib
# For systemverilog
../build/src/bin/kepler-formal -systemverilog design0.sv design1.sv \
  --liberty stdcells.lib.gz
# Short alias for systemverilog
../build/src/bin/kepler-formal -sv design0.sv design1.sv \
  --liberty stdcells.lib.gz
# For multi-file Verilog designs
../build/src/bin/kepler-formal -verilog --design1 <file...> --design2 <file...> \
  --liberty NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib fakeram45_64x32.lib fakeram45_64x15.lib
# For multi-file SystemVerilog designs
../build/src/bin/kepler-formal -sv --design1 <file...> --design2 <file...> \
  --liberty stdcells.lib.gz
# Through config file
../build/src/bin/kepler-formal --config test_config_naja_if.yaml
../build/src/bin/kepler-formal --config test_config_verilog.yaml
../build/src/bin/kepler-formal --config test_config_verilog_sec.yaml
../build/src/bin/kepler-formal --config test_config_verilog_tinyrocket_sec.yaml
python edit_tinyrocket_output_buffers.py
../build/src/bin/kepler-formal --config test_config_verilog_tinyrocket_buffered_sec.yaml
python extract_tinyrocket_csrfile_sec.py
../build/src/bin/kepler-formal --config test_config_verilog_tinyrocket_csrfile_sec.yaml
```

`test_config_verilog_sec.yaml` is a small SEC example that uses the local
`seq_design0.sv` / `seq_design1.sv` pair.

`test_config_verilog_tinyrocket_sec.yaml` keeps the original full-TinyRocket
SEC self-compare setup for experimentation with the same Liberty files as the
LEC example.

`edit_tinyrocket_output_buffers.py` generates
`tinyrocket_output_buffered.v`, which adds one combinational buffer on every
top-level TinyRocket output bit. The matching
`test_config_verilog_tinyrocket_buffered_sec.yaml` keeps behavior identical
while making the SEC pair structurally different.

`test_config_verilog_tinyrocket_csrfile_sec.yaml` extracts `CSRFile` from
`tinyrocket.v` with Najaeda and then runs a TinyRocket-derived SEC example that
first fails at `k = 3` on the existing `io_decode_0_fp_illegal` output.

The current SEC flow compares observed outputs only. Internal register naming
does not need to match between the two designs.

### YAML input_paths notes

- Flat list (single file per design):
  `input_paths: [design0.v, design1.v]`
- Nested list (multi-file per design):
  `input_paths: [[design0_a.v, design0_b.v], [design1_a.v, design1_b.v]]`

### Notes

- Supported YAML formats: `verilog`, `systemverilog`, `sv`, `naja_if`
- Supported library files on the existing `--liberty` / `liberty_files` path:
  - `.lib`
  - `.lib.gz`
- Python tech loaders are YAML-only under `py_tech_files`
