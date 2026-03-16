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
  --liberty stdcells.lib.gz primitives.py
# Short alias for systemverilog
../build/src/bin/kepler-formal -sv design0.sv design1.sv \
  --liberty stdcells.lib.gz primitives.py
# For multi-file Verilog designs
../build/src/bin/kepler-formal -verilog --design1 <file...> --design2 <file...> \
  --liberty NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib fakeram45_64x32.lib fakeram45_64x15.lib
# For multi-file SystemVerilog designs
../build/src/bin/kepler-formal -sv --design1 <file...> --design2 <file...> \
  --liberty stdcells.lib.gz primitives.py
# Through config file
../build/src/bin/kepler-formal --config test_config_naja_if.yaml
../build/src/bin/kepler-formal --config test_config_verilog.yaml
```

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
  - `.py`
