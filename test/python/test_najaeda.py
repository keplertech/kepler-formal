# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

from __future__ import annotations

import importlib
import importlib.machinery
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import kepler_formal
from kepler_formal import Design, VerificationOptions, VerificationStatus, verify
from kepler_formal import _native


def _nested_artifact_available() -> bool:
    package = Path(kepler_formal.__file__).resolve().parent / "najaeda"
    has_extension = any(
        (package / f"naja{suffix}").is_file()
        for suffix in importlib.machinery.EXTENSION_SUFFIXES
    )
    return has_extension and (package / "netlist.py").is_file()


def _run_isolated_python(source: str) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    return subprocess.run(
        [sys.executable, "-c", source],
        cwd=Path.cwd(),
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


class NestedNajaedaNamespaceTest(unittest.TestCase):
    def test_importing_kepler_formal_does_not_load_editor_runtime(self):
        completed = _run_isolated_python(
            "import sys, kepler_formal; "
            "assert 'kepler_formal.najaeda' not in sys.modules; "
            "assert 'najaeda' not in vars(kepler_formal)"
        )
        self.assertEqual(
            0,
            completed.returncode,
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
        )

    @unittest.skipUnless(
        _nested_artifact_available(),
        "nested NajaEDA native artifact is not staged in this build tree",
    )
    def test_nested_import_does_not_replace_top_level_najaeda(self):
        completed = _run_isolated_python(
            "import sys, types, kepler_formal; "
            "standalone = types.ModuleType('najaeda'); "
            "sys.modules['najaeda'] = standalone; "
            "from kepler_formal import najaeda; "
            "assert najaeda is not standalone; "
            "assert sys.modules['najaeda'] is standalone; "
            "assert najaeda.naja.__name__ == 'kepler_formal.najaeda.naja'"
        )
        self.assertEqual(
            0,
            completed.returncode,
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
        )


@unittest.skipUnless(
    _nested_artifact_available(),
    "nested NajaEDA native artifact is not staged in this build tree",
)
class NestedNajaedaRuntimeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.najaeda = importlib.import_module("kepler_formal.najaeda")
        cls.netlist = importlib.import_module("kepler_formal.najaeda.netlist")

    def setUp(self):
        self.netlist.reset()
        self.temporary = tempfile.TemporaryDirectory(
            prefix="kepler_formal_nested_najaeda_test_"
        )
        self.root = Path(self.temporary.name)
        self.reference = self.root / "reference.v"
        self.equivalent = self.root / "equivalent.v"
        self.reference.write_text(
            "module top(input a, output y); assign y = a; endmodule\n",
            encoding="utf-8",
        )
        self.equivalent.write_text(
            "module top(input a, output y); wire n; "
            "assign n = a; assign y = n; endmodule\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.netlist.reset()
        self.temporary.cleanup()

    def test_nested_version_and_import_identity(self):
        self.assertEqual(self.najaeda.version(), self.najaeda.__version__)
        self.assertTrue(self.najaeda.git_hash())
        self.assertEqual(
            "kepler_formal.najaeda.naja",
            self.najaeda.naja.__name__,
        )
        self.assertIs(self.netlist.naja, self.najaeda.naja)
        self.assertNotIn("najaeda", sys.modules)

    def test_editor_universe_survives_help_and_file_verification(self):
        editor_top = self.netlist.create_top("editor_top")
        editor_top.create_input_term("editor_input")
        universe = self.najaeda.naja.NLUniverse.get()

        help_result = _native.run(["--help"])
        self.assertEqual("no_result", help_result["status"])
        self.assertEqual(universe, self.najaeda.naja.NLUniverse.get())
        self.assertEqual("editor_top", str(editor_top))

        result = verify(
            self.reference,
            self.equivalent,
            options=VerificationOptions(log_file=self.root / "verify.log"),
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, result.status)
        self.assertEqual(universe, self.najaeda.naja.NLUniverse.get())
        self.assertEqual("editor_input", editor_top.get_term("editor_input").get_name())

    def test_live_editor_objects_are_not_verification_inputs(self):
        editor_top = self.netlist.create_top("editor_top")

        with self.assertRaises(TypeError):
            verify(
                Design(editor_top),
                self.equivalent,
                options=VerificationOptions(log_file=self.root / "unused.log"),
            )

        self.assertEqual("editor_top", str(editor_top))
        self.assertIsNotNone(self.najaeda.naja.NLUniverse.get())


if __name__ == "__main__":
    unittest.main()
