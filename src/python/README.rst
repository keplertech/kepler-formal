Kepler Formal Python package
============================

``kepler-formal`` provides ``kepler_formal``, a native Python interface to the
same LEC and SEC engine used by the Kepler Formal executable. NajaEDA loads or
creates netlists; Kepler verifies those existing objects directly in process.
Results are owning Python values.

Quick start
-----------

Load both designs through NajaEDA in one universe, then verify them:

.. code-block:: python

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
   # Both netlists remain available for editing and repeated verification.

The caller owns the designs, databases, and universe. Kepler borrows them
without copying or serializing them, and does not destroy them when a call
succeeds or fails. Load libraries and primitives with NajaEDA as part of design
preparation. YAML/JSON configuration and file-based verification remain in the
standalone ``kepler-formal`` executable.

Shared NajaEDA runtime and live designs
---------------------------------------

``najaeda`` is a separate runtime dependency. ``kepler_formal.najaeda`` and its
submodules are compatibility aliases to that original package. Both APIs share
one native runtime and live universe. Kepler validates NajaEDA's ABI, native
build, and runtime identity before borrowing an object.

``from_najaeda()`` accepts a raw ``SNLDesign`` or a high-level
``najaeda.netlist.Instance``. It resolves an Instance's current model
immediately and returns a ``NativeDesign`` handle:

.. code-block:: python

   from najaeda import netlist
   from kepler_formal import from_najaeda

   universe.setTopDesign(reference)
   reference_handle = from_najaeda(netlist.get_top())
   universe.setTopDesign(implementation)
   implementation_handle = from_najaeda(netlist.get_top())
   result = verify_designs(reference_handle, implementation_handle)

Changing the selected top later does not retarget a handle. Edits to the
captured design remain visible because no snapshot is made. The handle retains
the Python wrappers; ownership of the native design stays with the caller.
Destroying the design or resetting the universe invalidates its handles.
``verify_designs()`` accepts handles or raw designs; high-level Instances must
be captured explicitly.

Options and results
-------------------

``VerificationOptions`` selects ``VerificationMode``, ``Solver``, ``SecEngine``,
and ``SecEncoding`` and provides ``max_k``, ``allow_boundary_mismatch``,
``report_skipped_outputs``, ``log_file``, ``log_level``, and
``set_as_boundary``. Enum fields accept their exact string values. SEC engine,
encoding, and bound options require SEC; boundary-mismatch handling requires
LEC.

``set_as_boundary`` accepts ordered pairs of slash-separated instance paths,
one path relative to each supplied top design. Only leaf instances, whose models
have no child instances, may be selected. Hierarchical paths to leaves are
supported; nonleaf selections are rejected. For every selected instance,
its inputs become additional compared outputs and its outputs become shared
unconstrained inputs. Corresponding pin interfaces must match. These logical
verification boundaries work for both LEC and SEC without cloning or rewiring
the designs; the original designs and any caller-owned DNL remain reusable
after the call.

Use ``result.status`` for the verdict: the historical native ``exit_code`` is
not mode-independent, and LEC returns zero for both equivalent and different
designs. Non-equivalence, partial proof, inconclusive, unsupported, and ordinary
operational-error outcomes are returned as values. Invalid arguments raise
``TypeError`` or ``ValueError``; destroyed designs raise ``ReferenceError``;
runtime-safety failures and unexpected binding exceptions raise ``RuntimeError``.

``VerificationResult`` contains the status, mode, actual log path, SEC bound,
reason, extraction/proof counters, and output-name tuples. Its ``input_format``
is ``naja_design``. The result is frozen and remains valid after the caller
later destroys the designs.

``coverage_percent`` is ``100 * covered_outputs / total_outputs``. It measures
observed-output extraction coverage. ``proven_outputs`` and
``unproven_outputs`` describe proof progress when the engine reports it;
per-output names are not available on every engine path. Check ``status``
for the verdict. ``skipped_observed_outputs`` identifies outputs excluded by
extraction or coverage limitations.

Process constraints
-------------------

Verification calls are synchronous, serialized, non-reentrant, and hold
Python's GIL until the run finishes. Do not concurrently read, mutate, delete,
or reset the shared NajaEDA universe from native threads. There is no
in-process timeout or cancellation hook; callers that require hard cancellation
or crash isolation should manage the Python invocation in a separate process.

Kepler restores the caller's top selections, DNL, ordering metadata, expression
caches, configuration, and logger references after each call. It temporarily
replaces spdlog's process-global default logger; unrelated native threads using
that logger must be coordinated.

See ``docs/python-api.md`` for build instructions, the option table,
result-field semantics, and lifecycle details.
