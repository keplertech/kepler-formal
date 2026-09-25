# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Exercise native option validation directly, bypassing the Python wrapper."""

import unittest

from kepler_formal import _native


class NativeLatchOptionsTest(unittest.TestCase):
    def contract(self, **changes):
        return dict(mode="sec", latch_input_changes="single",
                    latch_initial_inputs=0, latch_initial_storage=0) | changes

    def rejected(self, options, error, message):
        with self.assertRaisesRegex(error, message):
            _native.verify_designs(None, None, options)

    def test_valid_contract_reaches_normal_design_validation(self):
        for changes in ("single", "any"):
            self.rejected(self.contract(latch_input_changes=changes), TypeError, "design1 must be a NativeDesign")

    def test_no_event_settings_preserve_legacy_with_default_or_explicit_gate(self):
        for options in ({}, {"latch_support": False}, {"latch_support": True}, {"mode": "sec"}):
            self.rejected(options, TypeError, "design1 must be a NativeDesign")

    def test_native_tuning_requires_contract_and_cannot_override_explicit_false(self):
        for key, value in (("latch_input_changes", "any"), ("latch_initial_inputs", 0),
                           ("latch_initial_storage", 0), ("latch_workers", 0),
                           ("latch_max_waves", 1), ("latch_max_states", 1),
                           ("latch_max_transactions", 1)):
            with self.subTest(key=key):
                self.rejected({"mode": "sec", key: value}, ValueError, "requires explicit")
                self.rejected({"mode": "sec", "latch_support": False, key: value},
                              ValueError, "requires latch_support")

    def test_native_requires_complete_contract(self):
        for key in ("latch_input_changes", "latch_initial_inputs", "latch_initial_storage"):
            for remove in (False, True):
                options = self.contract(**{key: None})
                if remove:
                    options.pop(key)
                with self.subTest(key=key, remove=remove):
                    self.rejected(options, ValueError, "requires explicit")

    def test_native_requires_real_boolean_gate(self):
        for value in (None, 0, 1, "true", []):
            with self.subTest(value=value):
                self.rejected(self.contract(latch_support=value), TypeError, "latch_support must be a bool")

    def test_native_rejects_nonbinary_initialization(self):
        for key in ("latch_initial_inputs", "latch_initial_storage"):
            for value, error in ((-1, OverflowError), (2, ValueError), (True, TypeError),
                                 (False, TypeError), ("0", TypeError), (0.0, TypeError)):
                with self.subTest(key=key, value=value):
                    self.rejected(self.contract(**{key: value}), error, ".+")

    def test_native_rejects_unknown_and_embedded_nul_keys(self):
        for key in ("latch_unknown", "latch_support\0suffix", "latch_input_changes\0"):
            with self.subTest(key=key):
                self.rejected(self.contract(**{key: True}), TypeError, "unknown.*option")

    def test_native_change_modes_are_exact(self):
        for value in ("", "ANY", "single\0suffix"):
            with self.subTest(value=value):
                self.rejected(self.contract(latch_input_changes=value), ValueError, "any or single")

    def test_native_rejects_lec_and_leaf_boundary_semantics(self):
        self.rejected(self.contract(mode="lec"), ValueError, "only supported for SEC")
        self.rejected(self.contract(set_as_boundary=(("a", "b"),)), ValueError, "complete top interface")

    def test_native_resource_limits_are_checked(self):
        for key in ("latch_max_waves", "latch_max_states", "latch_max_transactions",
                    "latch_max_symbolic_nodes", "latch_max_sat_conflicts", "latch_max_sat_decisions"):
            for value, error in ((0, ValueError), (-1, OverflowError), (2**128, OverflowError), (True, TypeError)):
                with self.subTest(key=key, value=value):
                    self.rejected(self.contract(**{key: value}), error, ".+")
        self.rejected(self.contract(latch_workers=2**31), ValueError, "nonnegative int")
        for key in ("latch_max_sat_conflicts", "latch_max_sat_decisions"):
            self.rejected(self.contract(**{key: 2**32}), ValueError, "unsigned int")


if __name__ == "__main__":
    unittest.main()
