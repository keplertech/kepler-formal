# Python API

Kepler Formal provides a native CPython package named `kepler_formal`. The
package calls the same C++ LEC and SEC engine as the `kepler-formal` executable
and returns an owning, structured result after each run.

This is a direct, in-process binding. It does not start the command-line
executable, use a subprocess, or communicate through MCP or another service.
It can verify source files or borrow live designs from NajaEDA without
serializing, cloning, or rebuilding them. The separately installed `najaeda`
package provides the single native netlist runtime shared by both APIs.

## Build and install

This development change requires the matching NajaEDA shared-runtime SDK,
currently version `0.7.24.dev0` in `thirdparty/naja`. That SDK has not been
published. Build both packages from this recursive checkout in one virtual
environment, with the native build dependencies installed:

```bash
python -m pip install 'scikit-build-core>=0.11.3,<0.12' build wheel
python -m pip install --no-build-isolation ./thirdparty/naja
python -m pip install --no-build-isolation .
```

The wheel CI builds and installs a separate local provider wheel for every
Python/platform combination. Publishing Kepler is blocked until the SDK is
released and its released version is pinned in `pyproject.toml` and
`ci/shared_naja_wheels.py`. An older NajaEDA wheel without this SDK cannot be
used for direct object sharing. Once that prerequisite is satisfied, ordinary
isolated `pip install .` builds can resolve the provider from the package index.

`BUILD_KEPLER_PYTHON=ON` is a Python-only CMake build. Build the standalone
executable separately with `BUILD_KEPLER_PYTHON=OFF`; it continues to use the
vendored Naja and does not depend on an installed NajaEDA package.

The build uses `scikit-build-core`, following NajaEDA's package layout. Native
build dependencies are the same as for the CMake build. The wheel matrix covers
Linux x86_64/aarch64, macOS arm64, and Windows AMD64; see the exact Python
versions in [wheel coverage](python-release.md#wheel-coverage).

Windows source builds require clang-cl with the MSVC SDK environment, CMake,
Ninja, and the vcpkg dependencies installed by `ci/windows_setup.ps1`. They
also require pregenerated Verilog parser sources and
`PREGENERATED_PARSER_SOURCES=ON`; the wheel workflow generates these sources
on Linux, following NajaEDA's Windows build approach. Ordinary MSVC `cl.exe`
is not supported for KF's bundled solver sources.

Maintainers can publish tested wheels using the manual
[Python release workflow](python-release.md).

## Shared NajaEDA runtime

`najaeda` is a runtime dependency of `kepler_formal` and is imported before
Kepler's native extension. It initializes the Naja runtime and publishes a
versioned native API that Kepler validates at import and before each live-design
call. A mismatched build, ABI, or runtime identity fails explicitly instead of
passing objects across an unsafe binary boundary.

For compatibility, `kepler_formal.najaeda` and all of its submodules are
aliases to the original package:

```python
import najaeda
import kepler_formal
from najaeda import netlist

assert kepler_formal.najaeda is najaeda
assert kepler_formal.najaeda.netlist is netlist
```

There is one package, one native runtime, and one live universe. Either import
spelling refers to the same Python module objects; new code should generally
use the original `najaeda` namespace.

## Compare live NajaEDA designs

Capture a raw `najaeda.naja.SNLDesign` or a high-level
`najaeda.netlist.Instance` with `from_najaeda()`, then pass the handles to
`verify_designs()`:

```python
import najaeda
from najaeda import netlist
from kepler_formal import from_najaeda, verify_designs

netlist.reset()
netlist.load_verilog("reference.v")
universe = najaeda.naja.NLUniverse.get()
reference_design = universe.getTopDesign()
implementation_design = reference_design.clone("implementation")

universe.setTopDesign(reference_design)
reference = from_najaeda(netlist.get_top())

# Select and edit the distinct implementation design in the same universe.
universe.setTopDesign(implementation_design)
implementation = from_najaeda(netlist.get_top())
# ...edit implementation through NajaEDA...

result = verify_designs(reference, implementation)
print(result.status)
```

`from_najaeda(Instance)` resolves the instance's current model immediately.
Changing NajaEDA's selected top later does not retarget the returned
`NativeDesign`. The handle retains the original Python object and its resolved
raw `SNLDesign`; its `source` and `najaeda_design` properties expose those two
objects. It retains the Python wrappers, not ownership of the native design or
universe. It is not a netlist snapshot: edits to the captured native design
remain visible to later calls, and explicit Naja destruction invalidates it.

`verify_designs()` also accepts raw `SNLDesign` objects and captures them for
the call. It deliberately rejects a high-level `Instance`; call
`from_najaeda(instance)` while the intended model is current so the selection
cannot change implicitly. If the design or its universe is destroyed, using
the handle raises `ReferenceError`.

The live API accepts the common `VerificationOptions` fields `mode`, `solver`,
`max_k`, `sec_engine`, `sec_encoding`, `allow_boundary_mismatch`,
`report_skipped_outputs`, `log_file`, and `log_level`. File-loading fields do
not apply. Non-default `input_format`, nonempty `libraries`,
`verilog_preprocessing=True`, and `compact=True` are rejected before native
verification starts. Load primitives and libraries into the NajaEDA universe
before capturing the designs.

The call is synchronous, serialized, and holds Python's GIL. Do not mutate,
delete, or reset the shared NajaEDA universe from another native thread while
verification is running. Kepler temporarily selects and analyzes the borrowed
designs, then restores the caller's universe top selections, DNL, ordering
metadata, expression caches, configuration, and logger references.

The file-oriented `verify()`, `run_config()`, and `run_cli()` entry points
still own their load/run lifecycle and reject an already-live Naja universe.
Finish with `netlist.reset()` before calling them, or use `verify_designs()`
when the editor universe must remain live.

## Compare two designs

```python
from kepler_formal import (
    Design,
    InputFormat,
    SecEngine,
    VerificationMode,
    VerificationOptions,
    VerificationStatus,
    verify,
)

result = verify(
    Design("reference.v", top="top"),
    Design("implementation.v", top="top"),
    options=VerificationOptions(
        input_format=InputFormat.VERILOG,
        mode=VerificationMode.SEC,
        libraries=("cells.lib",),
        sec_engine=SecEngine.PDR,
        max_k=32,
        log_file="verification.log",
    ),
)

if result.status is VerificationStatus.EQUIVALENT:
    print("proved equivalent")
else:
    print(result.status.value, result.reason)
```

`verify()` accepts a `Design`, one path, or a sequence of paths for each side.
It creates a temporary JSON configuration and synchronously runs the native
engine. Paths passed through `Design` and `VerificationOptions` are expanded
and made absolute relative to the process's current working directory.

## Design inputs and flists

`Design` has three fields:

- `files`: one path-like value or a sequence of source paths. It defaults to an
  empty sequence so that a SystemVerilog flist can be the only input.
- `top`: an optional top-module name.
- `flist`: an optional path to a SystemVerilog file list.

Each design must provide at least one source in `files`, a permitted `flist`,
or both. The exact rules depend on `input_format`:

| Input format | Design 1 | Design 2 | Verification modes |
| --- | --- | --- | --- |
| `verilog` | Verilog file(s), optional `top`, no flist | Verilog file(s), optional `top`, no flist | LEC or SEC |
| `systemverilog` | SystemVerilog file(s) and/or flist, optional `top` | SystemVerilog file(s) and/or flist, optional `top` | SEC only |
| `sv2v` | SystemVerilog file(s) and/or flist, optional `top` | Verilog file(s), optional `top`, no flist | SEC only |
| `naja_if` | Exactly one Naja IF snapshot, no `top` or flist | Exactly one Naja IF snapshot, no `top` or flist | LEC or SEC |

For example, a flist-only SystemVerilog comparison is:

```python
result = verify(
    Design(flist="reference.f", top="top"),
    Design(flist="implementation.f", top="top"),
    options=VerificationOptions(
        input_format=InputFormat.SYSTEMVERILOG,
        mode=VerificationMode.SEC,
    ),
)
```

The native SystemVerilog loader interprets the contents of each flist. The
Python layer resolves the flist path itself but does not rewrite paths inside
the flist.

## Verification options

Enum fields accept either the exported enum member or its exact string value.

| Field | Default | Meaning |
| --- | --- | --- |
| `input_format` | `InputFormat.VERILOG` | `verilog`, `systemverilog`, `sv2v`, or `naja_if` |
| `mode` | `VerificationMode.LEC` | `lec` or `sec` |
| `libraries` | empty | One Liberty path or a sequence of Liberty paths |
| `solver` | `Solver.KISSAT` | `kissat`, `cadical`, or `glucose` |
| `max_k` | native default (32 for SEC) | Non-negative SEC bound |
| `sec_engine` | native default (`pdr`) | `pdr`, `k_induction`, or `imc` |
| `sec_encoding` | native default (`dual_rail_steady`) | `dual_rail_steady` or `binary` |
| `verilog_preprocessing` | `False` | Enable Verilog preprocessing |
| `compact` | `False` | Enable compact verification mode |
| `allow_boundary_mismatch` | `False` | Permit supported extracted-boundary mismatches |
| `report_skipped_outputs` | `False` | Ask the native engine to write detailed skipped-output reports |
| `log_file` | native-generated path | Requested native log path |
| `log_level` | `info` | Native log level (`debug` enables debug logging) |

`max_k`, `sec_engine`, and `sec_encoding` are SEC-only and are rejected when
`mode` is LEC. SystemVerilog-based formats are also SEC-only in the high-level
API.

## Other entry points

An existing YAML or JSON configuration can be used directly:

```python
from kepler_formal import run_config

result = run_config("verify.yml")
```

`run_config()` passes the configuration to the same in-process engine. Unlike
`verify()`, it does not translate the configuration into high-level Python
options, so the native configuration parser is authoritative.

The lower-level `run_cli()` accepts the same argument sequence as the native
executable, excluding `argv[0]`:

```python
from kepler_formal import run_cli

result = run_cli(("-verilog", "reference.v", "implementation.v"))
```

`run_cli()` is still an in-process call; “CLI” describes only its argument
shape. Pass a sequence, not one string or path.

## Statuses and errors

Always use `result.status` for the semantic outcome:

| Status | Meaning |
| --- | --- |
| `NO_RESULT` | The invocation intentionally did not attempt verification, for example an empty argument list or `--help`. |
| `EQUIVALENT` | LEC found no difference, or SEC completed a proof of equivalence under the selected model and encoding. |
| `DIFFERENT` | LEC found a difference, or SEC found a counterexample. |
| `PARTIALLY_PROVED` | SEC proved some observed outputs, but not all of them. |
| `INCONCLUSIVE` | SEC completed without either a full proof or a counterexample, commonly because a bound or engine limit was reached. |
| `UNSUPPORTED` | The selected SEC workflow cannot analyze the design pair. |
| `ERROR` | Argument parsing, configuration, loading, or another operational step failed before a semantic verdict was produced. |

`DIFFERENT`, `PARTIALLY_PROVED`, `INCONCLUSIVE`, `UNSUPPORTED`, `NO_RESULT`,
and ordinary native `ERROR` outcomes are returned as values, not raised as
Python exceptions. For an `ERROR`, inspect `reason`, the native log output, and
`exit_code`; some early failures can provide only a general reason.

Python argument validation still raises `TypeError` or `ValueError`. A native
safety violation (for example, an active external Naja universe) or an
unexpected native exception raises `RuntimeError` because no normal result can
be produced.

Do not infer equivalence from `exit_code == 0`. The value preserves the native
program's historical exit convention, and LEC uses zero for both equivalent
and different designs. SEC currently uses zero for a proof, one for a partial
proof, two for inconclusive/unsupported, and three for a counterexample, but
`status` is the stable, mode-independent interpretation.

## Result fields

`VerificationResult` is a frozen, value-only dataclass:

| Field or property | Meaning |
| --- | --- |
| `status` | A `VerificationStatus` semantic outcome |
| `exit_code` | The native return code; do not use it alone as the verdict |
| `input_format` | Parsed input-format string, or `None` if parsing ended before it was available |
| `verification` | Parsed `lec` or `sec` mode, or `None` if unavailable |
| `log_file` | The actual log path selected by the engine, or `None` if no log was created |
| `bound` | The SEC bound associated with the result; zero for LEC or when unavailable |
| `reason` | Engine detail or an operational error explanation, when available |
| `covered_outputs` | Outputs included in the extracted SEC comparison |
| `total_outputs` | Existing observed outputs considered during SEC extraction |
| `proven_outputs` | Outputs reported proved by the SEC engine |
| `unproven_outputs` | Names of not-yet-proved outputs when the engine reports per-output proof progress |
| `skipped_observed_outputs` | Observed-output names excluded because of extraction or coverage limitations |
| `equivalent` | `True` only when `status is EQUIVALENT` |
| `conclusive` | `True` only for `EQUIVALENT` or `DIFFERENT` |
| `coverage_percent` | `100 * covered_outputs / total_outputs`, or `None` when `total_outputs` is zero |

The output counters are primarily meaningful for SEC. LEC and early failures
normally leave them at zero and the name tuples empty.

### Extraction coverage is not proof progress

`covered_outputs / total_outputs` measures whether observed outputs survived
extraction and were included in the SEC problem. Consequently,
`coverage_percent` is **not** “percent proved.” An output can be covered but
remain unproved.

`proven_outputs` and `unproven_outputs` describe proof progress. Some engine
paths report only an aggregate count; for equivalent or partial results the
binding then uses `covered_outputs` as the best available proved count. The
engine may not provide corresponding names, so an empty `unproven_outputs`
tuple does not by itself mean that every output was proved. Check `status`
first. `skipped_observed_outputs` belongs to extraction coverage, not to the
set of covered-but-unproved outputs.

## Lifetime, global state, and concurrency

Results contain only Python strings, integers, tuples, and enums, so they remain
valid after either kind of call. `NativeDesign` is different: it retains live
NajaEDA wrappers and is valid only while the captured design remains in the
active shared universe. The binding restores the Kepler solver/report settings
and spdlog logger references that it changes.

Kepler and Naja still use global design, solver, and logging state. The binding
therefore has these constraints:

- Verification calls are synchronous, serialized by a process-wide mutex, and
  not reentrant.
- The binding intentionally keeps Python's GIL for the entire native run.
  Other Python threads cannot execute Python code until verification returns.
  On free-threaded CPython, importing the binding normally enables the GIL.
  If it is forcibly disabled with `PYTHON_GIL=0` or `-X gil=0`, verification
  raises `RuntimeError`; use `PYTHON_GIL=1` or `-X gil=1` instead.
- The API does not currently provide an in-process timeout or cancellation
  hook. A caller that needs hard cancellation or crash isolation should place
  the Python call in a separately managed process.
- `verify_designs()` uses the shared live NajaEDA universe. Do not concurrently
  read, mutate, delete, or reset that universe from native threads. Destroying
  a captured design or calling `netlist.reset()` invalidates its handles.
- File entry points reject a live universe because their loader lifecycle must
  own the runtime. Reset it first, or use `verify_designs()`.
- `Design` and result fields remain value/file oriented. Live objects enter
  only through `from_najaeda()` or the raw-`SNLDesign` convenience accepted by
  `verify_designs()`.
- Kepler temporarily installs a process-global spdlog default/named logger and
  restores the previous loggers after the run. The mutex protects Kepler calls,
  but it cannot protect unrelated native threads. A host with C++ threads that
  concurrently use spdlog's global default logger must coordinate those
  threads or run verification in an isolated process.

## Python technology files

`py_tech_files` are not supported by the in-process verification driver. The
driver deliberately excludes Naja's CPython technology loader. This
restriction applies to file-oriented Kepler configuration. NajaEDA can still
load primitives into its shared live universe for `verify_designs()`. Liberty
libraries are supported by the file API through `VerificationOptions.libraries`.

If a YAML/JSON configuration supplied to `run_config()` contains
`py_tech_files`, the call returns `VerificationStatus.ERROR` with a nonzero
`exit_code` and an explanatory `reason`. Use the standalone `kepler-formal`
executable for a workflow that currently requires Python technology files.
