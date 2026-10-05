# Python Design Input

A design can be given as a Python script instead of a netlist file. The script
builds the design with the Naja Python API, the same API that
[Python primitive files](python-primitives.md) use for technology cells. The
mode is for netlists produced or transformed by a program: every module,
instance, net and connection is created by calling the API.

## Usage

CLI flag: `-python` (or `-py`). YAML: `format: python` (or `py`). The mode
works for LEC and SEC.

```bash
build/src/bin/kepler-formal -python -v sec \
  --design1 <script...> --design2 <script...> \
  --liberty <cells.lib> \
  [--python_design1_top <top>] [--python_design2_top <top>]
```

```yaml
format: python
verification: lec
input_paths:
  - [design1/cells.py, design1/top.py]
  - [design2/top.py]
py_tech_files:
  - ./my_primitives.py
python_design1_top: top
python_design2_top: top
```

## Script contract

Each script defines `constructLibrary(lib)`. Kepler-formal calls it with the
design library of one side, after the technology primitives are loaded. The
scripts of a design run in the listed order against the same library, so a
design may be split over several files.

```python
import naja


def constructLibrary(lib):
    inv = primitive(lib, "INV")
    top = naja.SNLDesign.create(lib, "top")
    a = naja.SNLScalarTerm.create(top, naja.SNLTerm.Direction.Input, "a")
    y = naja.SNLScalarTerm.create(top, naja.SNLTerm.Direction.Output, "y")
    net_a = naja.SNLScalarNet.create(top, "a")
    net_y = naja.SNLScalarNet.create(top, "y")
    a.setNet(net_a)
    y.setNet(net_y)
    gate = naja.SNLInstance.create(top, inv, "u_inv")
    gate.getInstTerm(inv.getScalarTerm("A")).setNet(net_a)
    gate.getInstTerm(inv.getScalarTerm("Y")).setNet(net_y)
```

Primitive cells come from `py_tech_files` or `--liberty`, as in every other
mode. They live in the primitive libraries of the design's database. A Liberty
library keeps the name declared in the file, so look a cell up in every
primitive library rather than by library name:

```python
def primitive(lib, name):
    for primitives in lib.getDB().getPrimitiveLibraries():
        model = primitives.getSNLDesign(name)
        if model is not None:
            return model
    raise RuntimeError("primitive %s was not loaded" % name)
```

## Top selection

`python_design1_top` and `python_design2_top` name the top module of each
design. Without them the top is the module that no other module instantiates.
When several modules qualify, name the top explicitly; the run stops with
`No top design was found` otherwise. A named top that the script did not
create stops the run with an error.

## Notes

- The script's `print` output is only shown when written to `sys.stderr` or
  flushed; a run that stops on an error does not flush Python's `stdout`.
- `naja.so` must be importable, as for Python primitive files; see
  [How `naja.so` is found](python-primitives.md#how-najaso-is-found).
- The in-process `kepler_formal` Python API does not accept this format. It
  borrows live NajaEDA designs instead; see [Python API](python-api.md).
