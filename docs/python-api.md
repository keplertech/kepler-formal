# Python API

Kepler Formal provides a native CPython package named `kepler_formal`. The
package calls the same C++ LEC and SEC engine as the `kepler-formal` executable
and returns an owning, structured result after each run.

This is a direct, in-process binding. It does not start the command-line
executable, use a subprocess, or communicate through MCP or another service.
The verification API is deliberately file-based. The distribution also
contains an opt-in, nested NajaEDA package for live netlist analysis and
editing, with a separate native runtime and lifecycle.

## Build and install

From a recursive source checkout:

```bash
python -m pip install .
```

The build uses `scikit-build-core`, following NajaEDA's package layout. Native
build dependencies are the same as for the CMake build. Linux and macOS are the
initially supported platforms.

## Bundled NajaEDA editor

Import the bundled editor explicitly from Kepler Formal's namespace:

```python
from kepler_formal import najaeda
from kepler_formal.najaeda import netlist

print(najaeda.__version__)

netlist.reset()
top = netlist.load_verilog("candidate.v")
top.get_net("old_name").set_name("new_name")
netlist.dump_naja_if("candidate.najaif")
```

Importing `kepler_formal` does not load this editor runtime. It is initialized
only by an explicit access such as `from kepler_formal import najaeda` or an
import below `kepler_formal.najaeda`. The usual upstream high-level modules
remain nested, for example:

- `kepler_formal.najaeda.netlist` for loading, navigation, analysis, and edits;
- `kepler_formal.najaeda.naja` for the expert raw SNL API;
- `kepler_formal.najaeda.instance_visitor`, `stats`, and `primitives` for the
  corresponding upstream helpers.

Always use these fully nested names in application code. The package does not
install or replace a `sys.modules["najaeda"]` alias, and its internal modules
resolve one another through `kepler_formal.najaeda`. A separately installed
top-level `najaeda` package is not used by the bundled editor. Prefer one API
namespace consistently in an application so object ownership remains clear.

### Two native runtimes

The editor and verifier deliberately have independent native Naja runtimes:

- `kepler_formal.najaeda` owns a persistent editing universe. Its live
  `Instance`, `Net`, `Term`, and raw SNL objects remain valid until the editor
  deletes them or `netlist.reset()` destroys that universe.
- `verify()`, `run_config()`, and `run_cli()` use Kepler's private verification
  universe. Each call creates and destroys its run state without destroying
  the nested editor's universe.

This separation allows a verification call while the nested editor has a live
design, but native objects cannot cross the boundary. `Design.files` accepts
path-like values only; it does not accept editor `Instance`, raw `SNLDesign`,
or other NajaEDA objects. Results likewise contain no live Naja objects.

Exchange a design through a file. Naja IF is the most direct snapshot format:

```python
from kepler_formal import Design, InputFormat, VerificationOptions, verify
from kepler_formal.najaeda import netlist

netlist.dump_naja_if("edited.najaif")
result = verify(
    Design("reference.najaif"),
    Design("edited.najaif"),
    options=VerificationOptions(input_format=InputFormat.NAJA_IF),
)
```

Verilog and SystemVerilog files are also valid exchange formats under their
normal input rules. Saving a path and then passing that path is intentional:
the verifier never borrows the editor's in-memory universe.

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

Results contain only Python strings, integers, tuples, and enums. They remain
valid after the call because the verification runtime's native Naja universe,
miter state, expression caches, and run-owned objects are released before
control returns to Python. The binding restores the Kepler solver/report
settings and spdlog logger references that it changes.

Within the private verification runtime, Kepler and Naja still use global
design, solver, and logging state. The binding therefore has these constraints:

- Verification calls are synchronous, serialized by a process-wide mutex, and
  not reentrant.
- The binding intentionally keeps Python's GIL for the entire native run.
  Other Python threads cannot execute Python code until verification returns.
- The API does not currently provide an in-process timeout or cancellation
  hook. A caller that needs hard cancellation or crash isolation should place
  the Python call in a separately managed process.
- A call is rejected with `RuntimeError` if Kepler's private verification
  runtime already contains an unexpected live universe. A universe owned by
  the isolated `kepler_formal.najaeda` editor is independent and may remain
  live across verification calls.
- No `Design` or result field accepts or returns live Naja/SNL/NajaEDA objects.
  Pass files or snapshots and keep only the owning result values.
- Kepler temporarily installs a process-global spdlog default/named logger and
  restores the previous loggers after the run. The mutex protects Kepler calls,
  but it cannot protect unrelated native threads. A host with C++ threads that
  concurrently use spdlog's global default logger must coordinate those
  threads or run verification in an isolated process.

## Python technology files

`py_tech_files` are not supported by the in-process verification driver. The
driver deliberately excludes Naja's CPython technology loader. This
restriction applies even though the separate nested NajaEDA editor has its own
primitive-loading facilities. Liberty libraries are supported through
`VerificationOptions.libraries`.

If a YAML/JSON configuration supplied to `run_config()` contains
`py_tech_files`, the call returns `VerificationStatus.ERROR` with a nonzero
`exit_code` and an explanatory `reason`. Use the standalone `kepler-formal`
executable for a workflow that currently requires Python technology files.
