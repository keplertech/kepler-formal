Kepler Formal Python package
============================

``kepler-formal`` provides ``kepler_formal``, a native Python interface to the
same LEC and SEC engine used by the Kepler Formal executable. Calls run directly
in the current process: this package is not a subprocess wrapper and does not
use MCP or another service.

The initial API is file-based. It is intended for running verification and
returning owning result values, not for editing live NajaEDA netlists.

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

Kepler and Naja use process-global native state. Verification calls are
synchronous, serialized, non-reentrant, and hold Python's GIL until the run
finishes. A call raises ``RuntimeError`` if another live Naja universe exists,
including one owned by a NajaEDA session. The API does not accept or return live
Naja/SNL/NajaEDA objects. There is no in-process timeout or cancellation hook;
callers that require hard cancellation or crash isolation should manage the
Python invocation in a separate process.

Run state and logger references are restored after each call, but Kepler
temporarily replaces spdlog's process-global default logger. Unrelated native
threads using that global logger must be coordinated, or verification should be
run in an isolated process.

Python technology files (``py_tech_files``) are not supported by this
in-process package. A configuration containing them returns ``ERROR`` with a
nonzero exit code and an explanatory reason. Liberty libraries are supported.

See ``docs/python-api.md`` in the source tree for the complete input rules,
option table, result-field semantics, and lifecycle notes.
