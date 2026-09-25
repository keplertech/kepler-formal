# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

import unittest

from kepler_formal import VerificationOptions
from kepler_formal.api import _build_native_design_options


class OpaquePolicyOptionsTest(unittest.TestCase):
    def test_default_is_disabled(self):
        self.assertFalse(_build_native_design_options(None)["error_on_opaque"])

    def test_sec_accepts_explicit_enablement(self):
        settings = VerificationOptions(mode="sec", error_on_opaque=True)
        self.assertTrue(_build_native_design_options(settings)["error_on_opaque"])

    def test_lec_rejects_enablement(self):
        with self.assertRaisesRegex(ValueError, "only supported for SEC"):
            _build_native_design_options(VerificationOptions(error_on_opaque=True))

    def test_non_boolean_is_rejected(self):
        for value in (None, 1, "true", []):
            with self.subTest(value=value), self.assertRaises(TypeError):
                _build_native_design_options(
                    VerificationOptions(mode="sec", error_on_opaque=value))


if __name__ == "__main__":
    unittest.main()
