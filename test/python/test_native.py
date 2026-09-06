# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

import unittest

import kepler_formal
from kepler_formal import _native


class NativeApiTest(unittest.TestCase):
    def test_version(self):
        self.assertEqual("1.0.0", kepler_formal.version())
        self.assertEqual(kepler_formal.version(), kepler_formal.__version__)
        self.assertTrue(kepler_formal.git_hash())

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
