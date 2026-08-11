# Custom Python Primitives

Kepler-formal can load technology primitive models written with the Naja Python
API. This is useful when a technology library needs formal models that are not
available in Liberty, including parameterized truth tables and sequential-cell
models.

## Configure Primitive Files

List each primitive file under `py_tech_files` in the YAML configuration:

```yaml
format: verilog
verification: lec
input_paths:
  - design1.v
  - design2.v
py_tech_files:
  - ./my_primitives.py
```

Each file imports `naja` and defines `constructPrimitives(lib)`. Kepler-formal
passes the shared primitives library to this function:

```python
import naja


def constructPrimitives(lib):
    inv = naja.SNLDesign.createPrimitive(lib, "INV")
    naja.SNLScalarTerm.create(
        inv, naja.SNLTerm.Direction.Input, "I"
    )
    naja.SNLScalarTerm.create(
        inv, naja.SNLTerm.Direction.Output, "O"
    )
    inv.setTruthTable(0b01)
```

Files are loaded in their listed order. See the
[Xilinx primitive model](../examples/xilinx/xilinx.py) and its
[register-slice](../examples/xilinx/register_slice) and
[VexRiscv](../examples/xilinx/vexriscv) examples for combinational,
parameterized, and sequential models.

## How `naja.so` Is Found

Python primitive files require the `naja` extension module. When
`py_tech_files` is present, kepler-formal:

1. Resolves the running executable, searching `PATH` when it was invoked by
   name.
2. Prepends the executable directory to `PYTHONPATH`.
3. Preserves every directory already present in `PYTHONPATH` after that entry.

The CMake build and install place `naja.so` next to `kepler-formal`, so the
normal layout works without manual environment configuration. An adjacent
module takes precedence over entries supplied by the user.

## Moving the Binary

Do not copy only the `kepler-formal` executable. Copy or deploy its complete
binary directory so that the compatible `naja.so` remains beside it:

```text
bin/
  kepler-formal
  naja.so
```

The module should come from the same build or distribution as the executable.
Using a module from an incompatible Naja build can cause import or ABI errors.

## Set `PYTHONPATH` Manually

If `naja.so` cannot be kept beside the executable, set `PYTHONPATH` to the
directory that contains it. Specify the directory, not the `.so` file:

```bash
export PYTHONPATH="/opt/kepler/python${PYTHONPATH:+:$PYTHONPATH}"
kepler-formal --config config.yaml
```

Kepler-formal keeps this value when it adds its own executable directory. If no
compatible `naja` module is available through either location, loading the
Python primitive file fails with a Python import error. Runs without
`py_tech_files` do not require the Python extension module.
