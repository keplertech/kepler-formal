Kepler Formal Python package
============================

``kepler-formal`` provides ``kepler_formal``, a native Python interface to the
same LEC and SEC engine used by the Kepler Formal executable. Calls run directly
in the current process: this package is not a subprocess wrapper and does not
use MCP or another service.

The API verifies source files or borrows live designs from its ``najaeda``
runtime dependency without serializing or copying them. Results are owning
Python values.

Quick start
-----------

.. code-block:: python

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
       ),
   )

   if result.status is VerificationStatus.EQUIVALENT:
       print("proved equivalent")
   else:
       print(result.status.value, result.reason)

Shared NajaEDA runtime and live designs
---------------------------------------

Import the original NajaEDA package and capture live designs explicitly:

.. code-block:: python

   import najaeda
   import kepler_formal
   from najaeda import netlist
   from kepler_formal import from_najaeda, verify_designs

   assert kepler_formal.najaeda is najaeda
   print(najaeda.__version__)
   netlist.reset()
   netlist.load_verilog("reference.v")
   universe = najaeda.naja.NLUniverse.get()
   reference_design = universe.getTopDesign()
   implementation_design = reference_design.clone("implementation")

   universe.setTopDesign(reference_design)
   reference = from_najaeda(netlist.get_top())
   universe.setTopDesign(implementation_design)
   implementation = from_najaeda(netlist.get_top())
   # Edit implementation through NajaEDA here.
   result = verify_designs(reference, implementation)

``kepler_formal.najaeda`` and its submodules are compatibility aliases to the
original package, not a bundled copy. There is one Python package, one native
runtime, and one live universe. Kepler validates NajaEDA's ABI, native build,
and runtime identity before borrowing an object.

``from_najaeda()`` accepts a raw ``SNLDesign`` or a high-level ``Instance``.
It resolves an Instance's model immediately and returns a stable
``NativeDesign`` handle. Changing the selected top later does not retarget the
handle; edits to the captured design remain visible because no snapshot or
copy is made. Destroying the design or resetting the universe invalidates the
handle. ``verify_designs()`` accepts handles or raw designs, but high-level
Instances must be captured explicitly.

Live verification supports the mode, solver, SEC engine/encoding/bound,
boundary, report, and logging options. File-only input format, Liberty path,
preprocessing, and compact options are rejected. The call is synchronous,
serialized, and holds the GIL; do not concurrently mutate or reset NajaEDA.

The file APIs ``verify()``, ``run_config()``, and ``run_cli()`` continue to
accept paths/configuration. They reject an already-live Naja universe because
their loader must own its lifecycle. Call ``netlist.reset()`` first, or use
``verify_designs()`` while the editor remains live.

``Design.files`` accepts one path or a sequence. Each design needs one or more
files, a permitted SystemVerilog flist, or both:

* ``verilog`` accepts files and an optional top on both sides, but no flists.
* ``systemverilog`` accepts files and/or a flist plus an optional top on both
  sides. It requires SEC.
* ``sv2v`` accepts SystemVerilog files and/or a flist for design 1 and Verilog
  files without a flist for design 2. It requires SEC.
* ``naja_if`` requires exactly one snapshot per side and accepts neither flists
  nor top-module selections.

Options and results
-------------------

``VerificationOptions`` selects ``InputFormat``, ``VerificationMode``,
``Solver``, ``SecEngine``, and ``SecEncoding`` and provides Liberty libraries,
the SEC bound, top-level verification flags, and logging settings. Enum fields
also accept their exact string values. SEC engine, encoding, and bound options
cannot be used with LEC.

The result status is one of ``NO_RESULT``, ``EQUIVALENT``, ``DIFFERENT``,
``PARTIALLY_PROVED``, ``INCONCLUSIVE``, ``UNSUPPORTED``, or ``ERROR``. Use
``result.status`` for the verdict: the historical native ``exit_code`` is not
mode-independent, and LEC returns zero for both equivalent and different
designs. Non-equivalence, partial proof, inconclusive, unsupported, no-result,
and ordinary operational-error outcomes are returned as values. Invalid
Python arguments raise ``TypeError`` or ``ValueError``; native safety failures
and unexpected native exceptions raise ``RuntimeError``.

``VerificationResult`` also contains the parsed format/mode, actual log path,
SEC bound, reason, extraction counters, proof counters, and output-name tuples.
The result is frozen and remains valid after native run state is released.

``coverage_percent`` is ``100 * covered_outputs / total_outputs``. It measures
observed-output extraction coverage, not proof progress. ``proven_outputs`` and
``unproven_outputs`` describe proof progress when the engine reports it;
per-output names are not available on every engine path. Always check
``status`` rather than interpreting an empty ``unproven_outputs`` tuple as a
complete proof. ``skipped_observed_outputs`` identifies outputs excluded by
extraction or coverage limitations.

Existing configuration and CLI-shaped arguments
------------------------------------------------

An existing YAML or JSON configuration can be run directly:

.. code-block:: python

   from kepler_formal import run_config

   result = run_config("verify.yml")

The lower-level ``run_cli(arguments)`` accepts the native executable's argument
sequence without ``argv[0]``. Both functions still invoke the engine in process;
neither launches the executable.

Process constraints
-------------------

Kepler and Naja use global state. Verification calls are synchronous,
serialized, non-reentrant, and hold Python's GIL until the run finishes.
``verify_designs()`` works in the shared live NajaEDA universe, so callers must
not concurrently read, mutate, delete, or reset it from native threads. File
entry points reject a live universe because their loader needs to own its
lifecycle. There is no in-process timeout or cancellation hook; callers that
require hard cancellation or crash isolation should manage the Python
invocation in a separate process.

Run state and logger references are restored after each call, but Kepler
temporarily replaces spdlog's process-global default logger. Unrelated native
threads using that global logger must be coordinated, or verification should be
run in an isolated process.

Python technology files (``py_tech_files``) are not supported by the
file-oriented in-process verification driver. A configuration containing them
returns ``ERROR`` with a nonzero exit code and an explanatory reason. NajaEDA
can load primitives into the shared universe for ``verify_designs()``;
Liberty-library paths remain supported by the file API.

See ``docs/python-api.md`` in the source tree for the complete input rules,
option table, result-field semantics, and lifecycle notes.
