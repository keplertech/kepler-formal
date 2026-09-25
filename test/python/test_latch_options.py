# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

import ctypes
import dataclasses
import unittest

from kepler_formal import VerificationOptions
from kepler_formal.api import _build_native_design_options


class LatchOptionsTest(unittest.TestCase):
    def contract(self, **changes):
        values = dict(mode="sec", latch_support=True, latch_input_changes="single",
                      latch_initial_inputs=0, latch_initial_storage=0)
        return VerificationOptions(**(values | changes))

    def test_default_is_disabled_without_an_invented_contract(self):
        for options in (None, VerificationOptions(), VerificationOptions(mode="sec")):
            result = _build_native_design_options(options)
            self.assertIs(result["latch_support"], False)
            for key in ("latch_input_changes", "latch_initial_inputs", "latch_initial_storage"):
                self.assertIsNone(result[key])

    def test_explicit_false_preserves_legacy(self):
        self.assertFalse(_build_native_design_options(VerificationOptions(latch_support=False))["latch_support"])

    def test_any_and_single_and_all_boolean_initial_values(self):
        for mode in ("any", "single"):
            for inputs in (0, 1):
                for storage in (0, 1):
                    with self.subTest(mode=mode, inputs=inputs, storage=storage):
                        result = _build_native_design_options(self.contract(
                            latch_input_changes=mode, latch_initial_inputs=inputs,
                            latch_initial_storage=storage))
                        self.assertTrue(result["latch_support"])
                        self.assertEqual(mode, result["latch_input_changes"])
                        self.assertEqual(inputs, result["latch_initial_inputs"])
                        self.assertEqual(storage, result["latch_initial_storage"])

    def test_tuning_does_not_enable_support(self):
        for key, value in (("latch_input_changes", "any"), ("latch_initial_inputs", 0),
                           ("latch_initial_storage", 1), ("latch_workers", 0),
                           ("latch_max_waves", 1), ("latch_max_states", 1),
                           ("latch_max_transactions", 1)):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "requires latch_support"):
                _build_native_design_options(VerificationOptions(mode="sec", **{key: value}))

    def test_enabled_requires_each_contract_field(self):
        for key in ("latch_input_changes", "latch_initial_inputs", "latch_initial_storage"):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "requires explicit"):
                _build_native_design_options(self.contract(**{key: None}))

    def test_enabled_lec_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "only supported for SEC"):
            _build_native_design_options(self.contract(mode="lec"))

    def test_boundaries_rejected_only_when_enabled(self):
        boundaries = (("left/cell", "right/cell"),)
        self.assertEqual(list(boundaries), _build_native_design_options(
            VerificationOptions(mode="sec", set_as_boundary=boundaries))["set_as_boundary"])
        with self.assertRaisesRegex(ValueError, "complete top interface"):
            _build_native_design_options(self.contract(set_as_boundary=boundaries))

    def test_master_gate_requires_bool(self):
        for value in (None, 0, 1, "true", [], {}):
            with self.subTest(value=value), self.assertRaises(TypeError):
                _build_native_design_options(self.contract(latch_support=value))

    def test_input_change_names_are_exact(self):
        for value in ("", "ANY", "Single", "any\0", "any ", "unknown"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                _build_native_design_options(self.contract(latch_input_changes=value))
        for value in (0, True, [], {}):
            with self.subTest(value=value), self.assertRaises(TypeError):
                _build_native_design_options(self.contract(latch_input_changes=value))

    def test_initialization_rejects_nonbinary(self):
        for key in ("latch_initial_inputs", "latch_initial_storage"):
            for value in (-1, 2, 256):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    _build_native_design_options(self.contract(**{key: value}))

    def test_initialization_rejects_bool_float_and_string(self):
        for key in ("latch_initial_inputs", "latch_initial_storage"):
            for value in (True, False, "0", 0.0, []):
                with self.subTest(key=key, value=value), self.assertRaises(TypeError):
                    _build_native_design_options(self.contract(**{key: value}))

    def test_resource_values_forwarded_and_zero_workers_is_auto(self):
        values = dict(latch_workers=0, latch_max_waves=12, latch_max_states=34, latch_max_transactions=56)
        result = _build_native_design_options(self.contract(**values))
        for key, value in values.items():
            self.assertEqual(value, result[key])

    def test_resource_limits_positive_workers_nonnegative(self):
        for key in ("latch_max_waves", "latch_max_states", "latch_max_transactions"):
            for value in (0, -1):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    _build_native_design_options(self.contract(**{key: value}))
        with self.assertRaises(ValueError):
            _build_native_design_options(self.contract(latch_workers=-1))

    def test_resource_limits_reject_wrong_types(self):
        for key in ("latch_workers", "latch_max_waves", "latch_max_states", "latch_max_transactions"):
            for value in (True, False, "1", 1.0, []):
                with self.subTest(key=key, value=value), self.assertRaises(TypeError):
                    _build_native_design_options(self.contract(**{key: value}))

    def test_integer_overflow_rejected(self):
        for key in ("latch_max_waves", "latch_max_states", "latch_max_transactions"):
            with self.subTest(key=key), self.assertRaises(ValueError):
                _build_native_design_options(self.contract(**{key: 1 << (8 * ctypes.sizeof(ctypes.c_size_t))}))
        with self.assertRaises(ValueError):
            _build_native_design_options(self.contract(latch_workers=1 << (8 * ctypes.sizeof(ctypes.c_int) - 1)))

    def test_symbolic_budgets_are_explicit_bounded_tuning(self):
        for key in ("latch_max_symbolic_nodes", "latch_max_sat_conflicts", "latch_max_sat_decisions"):
            with self.subTest(key=key):
                self.assertEqual(123, _build_native_design_options(self.contract(**{key: 123}))[key])
                with self.assertRaisesRegex(ValueError, "requires latch_support"):
                    _build_native_design_options(VerificationOptions(mode="sec", **{key: 123}))
                for value, error in ((0, ValueError), (-1, ValueError), (True, TypeError),
                                     (1.0, TypeError), ("1", TypeError), (2**128, ValueError)):
                    with self.subTest(value=value), self.assertRaises(error):
                        _build_native_design_options(self.contract(**{key: value}))
        for key in ("latch_max_sat_conflicts", "latch_max_sat_decisions"):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "unsigned int"):
                _build_native_design_options(self.contract(**{key: 1 << (8 * ctypes.sizeof(ctypes.c_uint))}))

    def test_opaque_policy_is_independent(self):
        self.assertTrue(_build_native_design_options(
            VerificationOptions(mode="sec", error_on_opaque=True))["error_on_opaque"])
        for strict in (False, True):
            result = _build_native_design_options(self.contract(error_on_opaque=strict))
            self.assertEqual(strict, result["error_on_opaque"])
            self.assertTrue(result["latch_support"])

    def test_options_remain_immutable(self):
        options = self.contract()
        with self.assertRaises(dataclasses.FrozenInstanceError):
            options.latch_support = False


if __name__ == "__main__":
    unittest.main()
