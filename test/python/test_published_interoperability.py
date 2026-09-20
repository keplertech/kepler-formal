# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from najaeda import naja, netlist

from kepler_formal import (
    VerificationOptions,
    VerificationStatus,
    _native,
    from_najaeda,
    verify_designs,
)


def _design_ids(design):
    identifier = design.getNLID()
    return (identifier.getDBID(), identifier.getLibraryID(), identifier.getDesignID())


@unittest.skipUnless(
    _native._provider_mode == "published", "requires the published NajaEDA adapter"
)
class PublishedDesignLifetimeTest(unittest.TestCase):
    def setUp(self):
        netlist.reset()
        self.temporary = tempfile.TemporaryDirectory(prefix="kepler_published_lifetime_")
        self.root = Path(self.temporary.name)
        self._create_universe()

    def _create_universe(self):
        self.universe = naja.NLUniverse.create()
        database = naja.NLDB.create(self.universe)
        self.library = naja.NLLibrary.create(database, "work")

    def _wire(self, name):
        design = naja.SNLDesign.create(self.library, name)
        input_term = naja.SNLScalarTerm.create(design, naja.SNLTerm.Direction.Input, "a")
        output_term = naja.SNLScalarTerm.create(design, naja.SNLTerm.Direction.Output, "y")
        wire = naja.SNLScalarNet.create(design, "wire")
        input_term.setNet(wire)
        output_term.setNet(wire)
        return design

    def tearDown(self):
        netlist.reset()
        self.temporary.cleanup()

    def _check_stale_and_live(self, stale, replacement):
        options = VerificationOptions(log_file=self.root / "verification.log")
        with self.assertRaises(ReferenceError):
            verify_designs(stale, replacement, options=options)
        result = verify_designs(replacement, replacement, options=options)
        self.assertEqual(VerificationStatus.EQUIVALENT, result.status)
        self.assertEqual("replacement", replacement.getName())

    def test_destroyed_handle_does_not_follow_reused_design_ids(self):
        original = self._wire("original")
        ids = _design_ids(original)
        handle = from_najaeda(original)
        original.destroy()

        replacement = self._wire("replacement")
        self.assertEqual(ids, _design_ids(replacement))
        self.assertIs(handle.najaeda_design, original)
        self.assertIs(self.universe.getSNLDesign(ids), replacement)
        self._check_stale_and_live(handle, replacement)

    def test_destroyed_handle_does_not_follow_a_recreated_universe(self):
        original = self._wire("original")
        ids = _design_ids(original)
        handle = from_najaeda(original)
        netlist.reset()
        self._create_universe()

        replacement = self._wire("replacement")
        self.assertEqual(ids, _design_ids(replacement))
        self.assertIs(handle.najaeda_design, original)
        self._check_stale_and_live(handle, replacement)


@unittest.skipUnless(
    _native._provider_mode == "published", "requires the published NajaEDA adapter"
)
class PublishedRuntimeIdentityTest(unittest.TestCase):
    def _run_isolated(self, source):
        with tempfile.TemporaryDirectory(prefix="kepler_published_identity_") as temporary:
            completed = subprocess.run(
                [sys.executable, "-X", "faulthandler", "-c", source],
                cwd=temporary,
                env=os.environ.copy(),
                capture_output=True,
                text=True,
                timeout=30,
            )
        self.assertEqual(
            0, completed.returncode,
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
        )

    def test_wrong_native_release_identity_rejects_import(self):
        for method in ("getVersion", "getGitHash"):
            with self.subTest(method=method):
                self._run_isolated(f"""
import najaeda
setattr(najaeda.naja, {method!r}, lambda: 'incorrect')
try:
    import kepler_formal
except ImportError as error:
    assert {method!r} in str(error), str(error)
else:
    raise AssertionError('incorrect native release identity was accepted')
""")

    def test_release_identity_is_checked_again_before_using_existing_handles(self):
        for method in ("getVersion", "getGitHash"):
            with self.subTest(method=method):
                self._run_isolated(f"""
from najaeda import naja
from kepler_formal import from_najaeda, verify_designs

universe = naja.NLUniverse.create()
try:
    library = naja.NLLibrary.create(naja.NLDB.create(universe), 'work')
    design = naja.SNLDesign.create(library, 'design')
    handle = from_najaeda(design)
    setattr(naja, {method!r}, lambda: 'incorrect')
    try:
        verify_designs(handle, handle)
    except RuntimeError as error:
        assert {method!r} in str(error), str(error)
    else:
        raise AssertionError('runtime identity was not checked for existing handles')
finally:
    universe.destroy()
""")

    def test_foreign_native_types_reject_import(self):
        for name in ("SNLDesign", "NLUniverse"):
            with self.subTest(type_name=name):
                self._run_isolated(f"""
import najaeda
setattr(najaeda.naja, {name!r}, object)
try:
    import kepler_formal
except ImportError as error:
    assert 'same native runtime' in str(error), str(error)
else:
    raise AssertionError('foreign native type was accepted')
""")


if __name__ == "__main__":
    unittest.main()
