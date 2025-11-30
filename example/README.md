# Examples

## Instructions

```bash
cd example
pip install najaeda
python edit.py
# For naja_if
../build/src/bin/kepler_formal -naja_if tinyrocket_naja.if tinyrocket_naja_edited.if NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib fakeram45_64x32.lib
# For verilog
../build/src/bin/kepler_formal -verilog tinyrocket_pre_edited.v tinyrocket_edited.v NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib /example/fakeram45_64x32.lib
# Through config file
../build/src/bin/kepler_formal --config test_config.yaml 
```
