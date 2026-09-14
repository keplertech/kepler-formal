# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

from __future__ import annotations

import importlib
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import kepler_formal
import najaeda
from najaeda import netlist

from kepler_formal import VerificationOptions, verify
from kepler_formal import _native


def _run_isolated_python(source: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", source],
        cwd=Path.cwd(),
        env=os.environ.copy(),
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )


class NajaedaAliasTest(unittest.TestCase):
    def test_top_level_and_nested_names_are_the_same_runtime(self):
        self.assertIs(najaeda, kepler_formal.najaeda)
        self.assertIs(
            importlib.import_module("najaeda.naja"),
            importlib.import_module("kepler_formal.najaeda.naja"),
        )
        self.assertIs(
            importlib.import_module("najaeda.netlist"),
            importlib.import_module("kepler_formal.najaeda.netlist"),
        )

    def test_alias_identity_in_both_import_orders(self):
        programs = (
            """
import importlib
import najaeda
import najaeda.netlist
import kepler_formal
assert kepler_formal.najaeda is najaeda
assert importlib.import_module('kepler_formal.najaeda.netlist') is najaeda.netlist
assert importlib.import_module('kepler_formal.najaeda.naja') is najaeda.naja
""",
            """
import importlib
import kepler_formal
nested = importlib.import_module('kepler_formal.najaeda.netlist')
import najaeda
import najaeda.netlist
assert kepler_formal.najaeda is najaeda
assert nested is najaeda.netlist
assert importlib.import_module('kepler_formal.najaeda.naja') is najaeda.naja
""",
        )
        for index, program in enumerate(programs):
            with self.subTest(import_order=index):
                completed = _run_isolated_python(program)
                self.assertEqual(
                    0,
                    completed.returncode,
                    f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
                )

    def test_dynamic_alias_preserves_canonical_import_metadata(self):
        completed = _run_isolated_python(
            """
import importlib
import importlib.resources
import najaeda
import kepler_formal

canonical = importlib.import_module('najaeda.primitives')
canonical_spec = canonical.__spec__
nested = importlib.import_module('kepler_formal.najaeda.primitives')
assert nested is canonical
assert canonical.__spec__ is canonical_spec
assert canonical.__spec__.name == 'najaeda.primitives'
source = importlib.resources.files('najaeda.primitives').joinpath('yosys.py')
assert 'def ' in source.read_text(encoding='utf-8')
assert importlib.reload(canonical) is canonical
"""
        )
        self.assertEqual(
            0,
            completed.returncode,
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
        )

    def test_incompatible_runtime_capsule_is_rejected_at_import(self):
        completed = _run_isolated_python(
            """
import ctypes
import najaeda

class Header(ctypes.Structure):
    _fields_ = [('abi_version', ctypes.c_uint32), ('struct_size', ctypes.c_size_t)]

header = Header(999, ctypes.sizeof(Header))
name = b'najaeda.naja._C_API'
new_capsule = ctypes.pythonapi.PyCapsule_New
new_capsule.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
new_capsule.restype = ctypes.py_object
najaeda.naja._C_API = new_capsule(ctypes.addressof(header), name, None)
try:
    import kepler_formal
except ImportError as error:
    assert 'Incompatible NajaEDA native runtime API' in str(error), str(error)
else:
    raise AssertionError('incompatible NajaEDA capsule was accepted')
"""
        )
        self.assertEqual(
            0,
            completed.returncode,
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
        )


class SharedNajaedaRuntimeTest(unittest.TestCase):
    def setUp(self):
        netlist.reset()
        self.temporary = tempfile.TemporaryDirectory(
            prefix="kepler_formal_shared_najaeda_test_"
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
        netlist.reset()
        self.temporary.cleanup()

    def test_file_api_rejects_an_active_editor_universe(self):
        editor_top = netlist.create_top("editor_top")

        with self.assertRaisesRegex(RuntimeError, "universe"):
            _native.run(["--help"])
        with self.assertRaisesRegex(RuntimeError, "universe"):
            verify(
                self.reference,
                self.equivalent,
                options=VerificationOptions(log_file=self.root / "unused.log"),
            )

        self.assertEqual("editor_top", editor_top.get_name())
        self.assertIsNotNone(najaeda.naja.NLUniverse.get())


if __name__ == "__main__":
    unittest.main()
