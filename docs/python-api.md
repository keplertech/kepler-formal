# Python API

Kepler Formal provides a native CPython package named `kepler_formal`. The
package calls the same C++ LEC and SEC engine as the `kepler-formal` executable
and returns an owning, structured result after each run.

This is a direct, in-process binding. It does not start the command-line
executable, use a subprocess, or communicate through MCP or another service.
NajaEDA loads or creates the netlists; Kepler borrows those live designs for
verification without serializing, cloning, or rebuilding them. The separately
installed `najaeda` package provides the single native netlist runtime shared
by both APIs.

## Build and install

For local regression without wheels or publishing, use the
[source regression runner](python-regression.md). It compiles both Python
packages from this checkout and tests their shared runtime.

The default development build uses the matching NajaEDA shared-runtime SDK,
version `0.7.24.dev0` in `thirdparty/naja`. Build both packages from this
recursive checkout in one virtual environment, with the native build
dependencies installed:

```bash
python -m pip install 'scikit-build-core>=0.11.3,<0.12' build wheel
python -m pip install --no-build-isolation ./thirdparty/naja
python -m pip install --no-build-isolation .
```

Default wheel CI runs on relevant pull requests, pushes to `main`, `v*` tags, and manual
runs. It builds and installs a separate local provider wheel for every
Python/platform combination. Only a manual `publish` request adds a second,
parallel set of jobs against published NajaEDA `0.7.24` wheels through KF's
version-specific compatibility adapter. Only published-provider jobs must
pass for publication; development-provider jobs remain independent regression
checks. Only wheels tested against the published provider are uploaded.
The adapter obtains matching release headers, links the installed native
libraries, and checks their identity before sharing designs. Those wheels
require `najaeda==0.7.24`, which pip installs normally. There is no separate
provider checkbox; see [release instructions](python-release.md).

`BUILD_KEPLER_PYTHON=ON` is a Python-only CMake build. Build the standalone
executable separately with `BUILD_KEPLER_PYTHON=OFF`; it continues to use the
vendored Naja and does not depend on an installed NajaEDA package.

The Python binding calls the shared verification engine with existing NajaEDA
designs. Source loading and configuration parsing belong to the standalone
executable; the Python extension does not link that file driver.

The build uses `scikit-build-core`, following NajaEDA's package layout. Native
build dependencies are the same as for the CMake build. The wheel matrix covers
Linux x86_64/aarch64, macOS arm64, and Windows AMD64; see the exact Python
versions in [wheel coverage](python-release.md#wheel-coverage).

Windows source builds require clang-cl with the MSVC SDK environment, CMake,
Ninja, and the vcpkg dependencies installed by `ci/windows_setup.ps1`. Building
the Windows NajaEDA provider also requires pregenerated Verilog parser sources
and `PREGENERATED_PARSER_SOURCES=ON`; the wheel workflow generates these sources
on Linux. Ordinary MSVC `cl.exe` is not supported for KF's bundled solver
sources.

Maintainers can publish tested wheels using the manual
[Python release workflow](python-release.md).

## Shared NajaEDA runtime

`najaeda` is a runtime dependency of `kepler_formal` and is imported before
Kepler's native extension. Development builds use its versioned native API to
check build and runtime identity. The published-provider adapter instead checks
the pinned release's native-file fingerprints, exported Python types, and live
universe identity. Both paths reject mismatches before accepting live designs.

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

## Load in NajaEDA, verify in Kepler

Use NajaEDA to load each design into the same live universe, then pass the raw
`SNLDesign` objects to `verify_designs()`:

```python
from najaeda import naja
from kepler_formal import VerificationOptions, verify_designs

universe = naja.NLUniverse.create()
reference_db = naja.NLDB.create(universe)
reference_db.loadVerilog(["reference.v"])
reference = reference_db.getTopDesign()

implementation_db = naja.NLDB.create(universe)
implementation_db.loadVerilog(["implementation.v"])
implementation = implementation_db.getTopDesign()

result = verify_designs(
    reference,
    implementation,
    options=VerificationOptions(log_file="verification.log"),
)
print(result.status, result.reason)
# Both designs remain available for edits and further verification calls.
```

Load libraries and primitives through NajaEDA before verification. Existing
in-memory designs can be passed directly; no loading step is required.
The caller owns the designs, databases, and universe. Kepler does not destroy
them on success, non-equivalence, or an error.

For a high-level `najaeda.netlist.Instance`, capture its current model with
`from_najaeda()` before changing the selected top:

```python
from najaeda import netlist
from kepler_formal import from_najaeda

universe.setTopDesign(reference)
reference_handle = from_najaeda(netlist.get_top())
universe.setTopDesign(implementation)
implementation_handle = from_najaeda(netlist.get_top())
result = verify_designs(reference_handle, implementation_handle)
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

`VerificationOptions` controls the mode, solver, SEC engine/encoding/bound,
boundary handling, reports, and logging. Source formats, libraries, and
preprocessing are configured through NajaEDA when loading designs.

The call is synchronous, serialized, and holds Python's GIL. Do not mutate,
delete, or reset the shared NajaEDA universe from another native thread while
verification is running. Kepler temporarily selects and analyzes the borrowed
designs, then restores the caller's universe top selections, DNL, ordering
metadata, expression caches, configuration, and logger references.

YAML/JSON configuration and file-based verification remain available through
the standalone `kepler-formal` executable. The Python package exposes
`from_najaeda()` and `verify_designs()` for existing netlists.

## Verification options

Enum fields accept either the exported enum member or its exact string value.

| Field | Default | Meaning |
| --- | --- | --- |
| `mode` | `VerificationMode.LEC` | `lec` or `sec` |
| `solver` | `Solver.KISSAT` | `kissat`, `cadical`, or `glucose` |
| `max_k` | native default (32 for SEC) | Non-negative SEC bound |
| `sec_engine` | native default (`pdr`) | `pdr`, `k_induction`, or `imc` |
| `sec_encoding` | native default (`dual_rail_steady`) | `dual_rail_steady` or `binary` |
| `learn_internal_relations` | `True` | Learn proved internal register equalities for SEC |
| `allow_x_equality_in_internal_relations` | `True` | Allow internal X/X relations without changing the final output property |
| `allow_boundary_mismatch` | `False` | Permit supported extracted-boundary mismatches |
| `report_skipped_outputs` | `False` | Ask the native engine to write detailed skipped-output reports |
| `log_file` | `None` | Requested native log path; LEC selects a default path when omitted |
| `log_level` | `info` | Native log level (`debug` enables debug logging) |
| `set_as_boundary` | `()` | Ordered `(design1_path, design2_path)` instance-path pairs to treat as shared boundaries |

`max_k`, `sec_engine`, and `sec_encoding` are SEC-only and are rejected when
`mode` is LEC. Internal relation options may be changed only in SEC mode.
`allow_boundary_mismatch` is supported only for LEC.
`log_file` is expanded and resolved relative to the current working directory.

### Treat selected instances as shared boundaries

Use `set_as_boundary` to remove paired block implementations from the proof
and verify surrounding logic through their exposed interfaces. Each item pairs
the instance path in `design1` with its corresponding path in `design2`. Paths
are slash-separated and relative to the supplied top designs. The selected
instances must be leaves: their models must have no child instances after
loading/elaboration. Hierarchical paths to leaves are valid, but selecting a
nonleaf instance is rejected:

```python
options = VerificationOptions(
    set_as_boundary=[
        ("subsystem/memory", "u_subsystem/u_memory"),
        ("clocking/gate", "u_clocking/icg"),
    ],
    log_file="boundary-verification.log",
)
result = verify_designs(reference, implementation, options=options)
```

The outer collection and each pair may be a list or tuple. Every pair must
contain exactly two non-empty strings. Pair order is significant and duplicate
or overlapping selections are rejected by native boundary validation.

For each selected instance, Kepler treats its input pins as extra compared
outputs, so the proof checks that the two designs drive the block identically.
It treats the instance's output pins as shared inputs, so both sides see the
same unconstrained block response, without traversing the block's internals.
Scalar and bus pin names, directions, widths, and ranges must match across each
pair. `allow_boundary_mismatch` does not relax this selected-boundary interface
check. Input pins must be connected, with exactly one driver on nonconstant
input nets. Unused output pins are allowed; inout pins and aliased,
constant-connected or multiply driven output nets are rejected. Direct internal
constant-wire ties are rejected, but primitive truth-table constants are supported.

Boundary selection works for both LEC and SEC without cloning or rewiring
the netlists. The caller's designs, connectivity, selected tops, and cached
DNL remain unchanged and can be reused afterward. This option belongs to the
live-design API. The Python extension does not expose file loading or
command-line configuration parsing.

## Statuses and errors

Always use `result.status` for the semantic outcome:

| Status | Meaning |
| --- | --- |
| `EQUIVALENT` | LEC found no difference, or SEC completed a proof of equivalence under the selected model and encoding. |
| `DIFFERENT` | LEC found a difference, or SEC found a counterexample. |
| `PARTIALLY_PROVED` | SEC proved some observed outputs, but not all of them. |
| `INCONCLUSIVE` | SEC completed without either a full proof or a counterexample, commonly because a bound or engine limit was reached. |
| `UNSUPPORTED` | The selected SEC workflow cannot analyze the design pair. |
| `ERROR` | An operational step failed before a semantic verdict was produced; inspect `reason`. |

`DIFFERENT`, `PARTIALLY_PROVED`, `INCONCLUSIVE`, `UNSUPPORTED`, and ordinary
native `ERROR` outcomes are returned as values. For an `ERROR`, inspect `reason`,
the native log output, and `exit_code`; some early failures can provide only a
general reason.

Invalid arguments raise `TypeError` or `ValueError`; destroyed designs raise
`ReferenceError`. Native runtime-safety failures and unexpected binding
exceptions raise `RuntimeError`. Handles remain reusable after failed calls
while their caller-owned designs remain live.

Do not infer equivalence from `exit_code == 0`. The value preserves the native
program's historical exit convention, and LEC uses zero for both equivalent
and different designs. SEC uses zero for a proof, one for a partial proof, two
for inconclusive/unsupported, and three for a counterexample, but
`status` is the stable, mode-independent interpretation.

## Result fields

`VerificationResult` is a frozen, value-only dataclass:

| Field or property | Meaning |
| --- | --- |
| `status` | A `VerificationStatus` semantic outcome |
| `exit_code` | The native return code; do not use it alone as the verdict |
| `input_format` | `naja_design` for live-design verification |
| `verification` | Selected `lec` or `sec` mode |
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
valid after verification and after the caller later destroys the designs.
`NativeDesign` retains live NajaEDA wrappers and is valid only while the
captured design remains in the
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
- Kepler temporarily installs a process-global spdlog default/named logger and
  restores the previous loggers after the run. The mutex protects Kepler calls,
  but it cannot protect unrelated native threads. A host with C++ threads that
  concurrently use spdlog's global default logger must coordinate those
  threads or run verification in an isolated process.
