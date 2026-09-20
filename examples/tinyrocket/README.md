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

## Python: NajaEDA loads, Kepler verifies

Use the matching standalone `najaeda` and `kepler_formal` packages built from
this checkout. The [source regression runner](../../docs/python-regression.md)
builds both locally from pinned sources with CMake, without wheels or publishing,
and tests the two commands below.

```bash
# Original versus edited TinyRocket: expected DIFFERENT.
python verify_python.py

# Load the original twice into separate databases: expected EQUIVALENT.
python verify_python.py tinyrocket.v tinyrocket.v
```

The script imports `najaeda` directly, loads the Liberty libraries and each
Verilog netlist into its own Naja database, then passes both live designs to
`kepler_formal.verify_designs()`. It prints the package locations, the verdict,
and the design names after verification to show that they remain available.
Default input and library paths are relative to the script, so it can also be
run from another directory.

The log defaults to `tinyrocket_python.log` in the working directory; override
it with `--log-file PATH`. Exit codes are 0 for equivalent, 1 for different,
and 2 for an error or another verdict. A difference is expected for the edited
example and is not an installation failure.

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
