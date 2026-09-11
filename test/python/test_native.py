# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

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
    _native.run(["--help"])
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

    def test_help(self):
        result = _native.run(["--help"])
        self.assertEqual(0, result["exit_code"])
        self.assertEqual("no_result", result["status"])

    def test_argument_conversion(self):
        with self.assertRaises(TypeError):
            _native.run("--help")
        with self.assertRaises(ValueError):
            _native.run(["bad\0argument"])
        with self.assertRaises(TypeError):
            _native.run([object()])


if __name__ == "__main__":
    unittest.main()
