# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

import tempfile
import unittest
from pathlib import Path

import najaeda
from najaeda import netlist

from kepler_formal import (
    SecEncoding,
    VerificationMode,
    VerificationOptions,
    VerificationStatus,
    verify_designs,
)


class PythonApiTest(unittest.TestCase):
    def setUp(self):
        netlist.reset()
        self.temporary = tempfile.TemporaryDirectory(
            prefix="kepler_formal_python_test_"
        )
        self.root = Path(self.temporary.name)
        self.universe = najaeda.naja.NLUniverse.create()
        self.reference = self._load(
            "reference", "assign y = a;"
        )
        self.equivalent = self._load(
            "equivalent", "wire n; assign n = a; assign y = n;"
        )
        self.different = self._load("different", "assign y = 1'b0;")

    def _load(self, name, body):
        source = self.root / f"{name}.v"
        source.write_text(
            f"module top(input a, output y); {body} endmodule\n",
            encoding="utf-8",
        )
        database = najaeda.naja.NLDB.create(self.universe)
        database.loadVerilog([str(source)])
        return database.getTopDesign()

    def tearDown(self):
        netlist.reset()
        self.temporary.cleanup()

    def _options(self, name):
        return VerificationOptions(log_file=self.root / name)

    def test_equivalent_and_different_lec_are_structured(self):
        passed = verify_designs(
            self.reference, self.equivalent, options=self._options("passed.log")
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, passed.status)
        self.assertTrue(passed.equivalent)
        self.assertTrue(passed.conclusive)
        self.assertEqual("naja_design", passed.input_format)
        self.assertEqual(str((self.root / "passed.log").resolve()), passed.log_file)
        self.assertTrue((self.root / "passed.log").is_file())

        failed = verify_designs(
            self.reference, self.different, options=self._options("different.log")
        )
        self.assertEqual(VerificationStatus.DIFFERENT, failed.status)
        self.assertFalse(failed.equivalent)
        # LEC retains its historical process code for both semantic verdicts.
        self.assertEqual(0, failed.exit_code)

    def test_repeated_runs_preserve_the_callers_loaded_designs(self):
        top_before = self.universe.getTopDesign()
        statuses = []
        for index, candidate in enumerate(
            (self.equivalent, self.different, self.equivalent)
        ):
            statuses.append(verify_designs(
                self.reference, candidate,
                options=self._options(f"repeat_{index}.log"),
            ).status)
            self.assertIs(najaeda.naja.NLUniverse.get(), self.universe)
            self.assertIs(self.universe.getTopDesign(), top_before)
            self.assertEqual("top", self.reference.getName())
            self.assertEqual("top", candidate.getName())
        self.assertEqual(
            [VerificationStatus.EQUIVALENT, VerificationStatus.DIFFERENT,
             VerificationStatus.EQUIVALENT],
            statuses,
        )

    def test_sec_counterexample_does_not_report_proved_outputs(self):
        result = verify_designs(
            self.reference,
            self.different,
            options=VerificationOptions(
                mode=VerificationMode.SEC,
                sec_encoding=SecEncoding.BINARY,
                max_k=2,
                log_file=self.root / "different-sec.log",
            ),
        )
        self.assertEqual(VerificationStatus.DIFFERENT, result.status)
        self.assertEqual(3, result.exit_code)
        self.assertEqual(1, result.total_outputs)
        self.assertEqual(1, result.covered_outputs)
        self.assertEqual(0, result.proven_outputs)

    def test_validation(self):
        invalid_options = (
            (ValueError, "SEC", VerificationOptions(max_k=2)),
            (TypeError, "max_k", VerificationOptions(
                mode=VerificationMode.SEC, max_k=True)),
            (ValueError, "non-negative", VerificationOptions(
                mode=VerificationMode.SEC, max_k=-1)),
            (ValueError, "solver", VerificationOptions(solver="unknown")),
            (TypeError, "allow_boundary_mismatch", VerificationOptions(
                allow_boundary_mismatch="false")),
            (ValueError, "log_file", VerificationOptions(log_file="")),
        )
        for error_type, message, options in invalid_options:
            with self.subTest(options=options):
                with self.assertRaisesRegex(error_type, message):
                    verify_designs(self.reference, self.equivalent, options=options)
        with self.assertRaisesRegex(TypeError, "NativeDesign.*SNLDesign"):
            verify_designs("reference.v", self.equivalent)


if __name__ == "__main__":
    unittest.main()
