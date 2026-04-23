# Examples

## Instructions

```bash

# Through binary flags

## Verilog example
../build/src/bin/kepler-formal -verilog tinyrocket.v tinyrocket_edited.v \
  NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib fakeram45_64x32.lib fakeram45_64x15.lib

## naja_if example
../build/src/bin/kepler-formal -naja_if tinyrocket_naja.if tinyrocket_naja_edited.if \
  NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib fakeram45_64x32.lib fakeram45_64x15.lib

# Through config file

../build/src/bin/kepler-formal --config test_config_verilog.yaml

```
