# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

import os
import subprocess
import sys
import sysconfig
import unittest

import kepler_formal
from kepler_formal import _native


class NativeApiTest(unittest.TestCase):
    def test_version(self):
        self.assertRegex(kepler_formal.version(), r"^[0-9]+\.[0-9]+\.[0-9]+$")
        self.assertEqual(kepler_formal.version(), kepler_formal.__version__)
        self.assertEqual(kepler_formal.version(), _native.get_version())
        self.assertTrue(kepler_formal.git_hash())
        self.assertEqual(kepler_formal.git_hash(), _native.get_git_hash())

    @unittest.skipUnless(
        sysconfig.get_config_var("Py_GIL_DISABLED"),
        "requires a free-threaded Python build",
    )
    def test_forced_disabled_gil_rejects_verification(self):
        code = """
import sys
from kepler_formal import _native

assert not sys._is_gil_enabled(), "PYTHON_GIL=0 was not honored"
try:
    _native.verify_designs(None, None, {})
except RuntimeError as error:
    assert "requires Python's GIL" in str(error), str(error)
else:
    raise AssertionError("verification was allowed without the GIL")
"""
        completed = subprocess.run(
            [sys.executable, "-c", code],
            env={**os.environ, "PYTHON_GIL": "0"},
            capture_output=True,
            text=True,
            timeout=30,
        )
        self.assertEqual(0, completed.returncode, completed.stdout + completed.stderr)

    def test_native_api_only_accepts_existing_designs(self):
        self.assertFalse(hasattr(_native, "run"))
        with self.assertRaisesRegex(TypeError, "design1 must be a NativeDesign"):
            _native.verify_designs(None, None, {})

    def test_native_option_conversion(self):
        with self.assertRaisesRegex(TypeError, "options must be a dict"):
            _native.verify_designs(None, None, [])
        with self.assertRaisesRegex(ValueError, "cannot contain NUL"):
            _native.verify_designs(None, None, {"mode": "lec\0"})
        with self.assertRaisesRegex(TypeError, "unknown.*option"):
            _native.verify_designs(None, None, {"verilog": "design.v"})


if __name__ == "__main__":
    unittest.main()
