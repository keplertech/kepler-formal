# TinyRocket

Run these commands from `examples/tinyrocket`.

The Python generators require `najaeda` from the same Naja revision as the
kepler-formal build. Install it once with `python3 -m pip install ../../thirdparty/naja`.

## LEC

The LEC examples compare the original TinyRocket netlist with an intentionally
edited version and are expected to find a difference.

```bash
# Verilog LEC
../../build/src/bin/kepler-formal -verilog \
  tinyrocket.v tinyrocket_edited.v \
  NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib \
  fakeram45_64x32.lib fakeram45_64x15.lib

# Naja IF LEC
python3 to_naja_if.py
python3 edit.py
../../build/src/bin/kepler-formal -naja_if \
  tinyrocket_naja.if tinyrocket_naja_edited.if \
  NangateOpenCellLibrary_typical.lib fakeram45_1024x32.lib \
  fakeram45_64x32.lib fakeram45_64x15.lib
```

## SEC Self-Check

This compares `tinyrocket.v` with itself. It exercises the full sequential
model and reset-unanchored state handling, and must not find a counterexample.

```bash
../../build/src/bin/kepler-formal \
  --config test_config_verilog_tinyrocket_sec.yaml
```

## SEC CSRFile Difference

This isolates `CSRFile` and compares it with a three-register delayed edit. The
test is expected to find a counterexample within the configured four frames.

```bash
python3 extract_tinyrocket_csrfile_sec.py
../../build/src/bin/kepler-formal \
  --config test_config_verilog_tinyrocket_csrfile_sec.yaml
```
