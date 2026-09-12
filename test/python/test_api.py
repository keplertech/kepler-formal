# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

import json
import tempfile
import unittest
from pathlib import Path

from kepler_formal import (
    Design,
    InputFormat,
    SecEncoding,
    VerificationMode,
    VerificationOptions,
    VerificationStatus,
    run_cli,
    run_config,
    verify,
)


class PythonApiTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="kepler_formal_python_test_"
        )
        self.root = Path(self.temporary.name)
        self.reference = self.root / "reference.v"
        self.equivalent = self.root / "equivalent.v"
        self.different = self.root / "different.v"
        self.reference.write_text(
            "module top(input a, output y); assign y = a; endmodule\n",
            encoding="utf-8",
        )
        self.equivalent.write_text(
            "module top(input a, output y); wire n; assign n = a; assign y = n; endmodule\n",
            encoding="utf-8",
        )
        self.different.write_text(
            "module top(input a, output y); assign y = 1'b0; endmodule\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.temporary.cleanup()

    def _options(self, name):
        return VerificationOptions(log_file=self.root / name)

    def test_equivalent_and_different_lec_are_structured(self):
        passed = verify(
            Design(self.reference, top="top"),
            Design(self.equivalent, top="top"),
            options=self._options("passed.log"),
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, passed.status)
        self.assertTrue(passed.equivalent)
        self.assertTrue(passed.conclusive)
        self.assertEqual(InputFormat.VERILOG.value, passed.input_format)
        self.assertEqual(str((self.root / "passed.log").resolve()), passed.log_file)
        self.assertTrue((self.root / "passed.log").is_file())

        failed = verify(
            self.reference,
            self.different,
            options=self._options("different.log"),
        )
        self.assertEqual(VerificationStatus.DIFFERENT, failed.status)
        self.assertFalse(failed.equivalent)
        # LEC preserves its historical process code; the structured status is
        # what distinguishes this result from equivalence.
        self.assertEqual(0, failed.exit_code)

    def test_repeated_runs_release_native_state(self):
        statuses = [
            verify(
                self.reference,
                candidate,
                options=self._options(f"repeat_{index}.log"),
            ).status
            for index, candidate in enumerate(
                (self.equivalent, self.different, self.equivalent)
            )
        ]
        self.assertEqual(
            [
                VerificationStatus.EQUIVALENT,
                VerificationStatus.DIFFERENT,
                VerificationStatus.EQUIVALENT,
            ],
            statuses,
        )

    def test_existing_config(self):
        config = self.root / "verify.json"
        config.write_text(
            json.dumps(
                {
                    "format": "verilog",
                    "verification": "lec",
                    "input_paths": [
                        [str(self.reference)],
                        [str(self.equivalent)],
                    ],
                    "log_file": str(self.root / "config.log"),
                }
            ),
            encoding="utf-8",
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, run_config(config).status)

    def test_systemverilog_flist_only_sec(self):
        first_flist = self.root / "first.f"
        second_flist = self.root / "second.f"
        first_source = self.root / "first.sv"
        second_source = self.root / "second.sv"
        source = (
            "module top(input logic a, output logic y); "
            "assign y = a; endmodule\n"
        )
        first_source.write_text(source, encoding="utf-8")
        second_source.write_text(source, encoding="utf-8")
        first_flist.write_text(f"{first_source}\n", encoding="utf-8")
        second_flist.write_text(f"{second_source}\n", encoding="utf-8")

        result = verify(
            Design(flist=first_flist, top="top"),
            Design(flist=second_flist, top="top"),
            options=VerificationOptions(
                input_format=InputFormat.SYSTEMVERILOG,
                mode=VerificationMode.SEC,
                sec_encoding=SecEncoding.BINARY,
                max_k=2,
                log_file=self.root / "flist-sec.log",
            ),
        )

        self.assertEqual(VerificationStatus.EQUIVALENT, result.status)
        self.assertEqual(InputFormat.SYSTEMVERILOG.value, result.input_format)
        self.assertEqual(1, result.total_outputs)
        self.assertEqual(1, result.covered_outputs)
        self.assertEqual(1, result.proven_outputs)

    def test_sec_counterexample_does_not_report_proved_outputs(self):
        result = verify(
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

    def test_python_technology_files_are_rejected_in_process(self):
        config = self.root / "python-tech.json"
        config.write_text(
            json.dumps(
                {
                    "format": "verilog",
                    "verification": "lec",
                    "input_paths": [
                        [str(self.reference)],
                        [str(self.equivalent)],
                    ],
                    "py_tech_files": [str(self.root / "primitives.py")],
                }
            ),
            encoding="utf-8",
        )

        result = run_config(config)
        self.assertEqual(VerificationStatus.ERROR, result.status)
        self.assertNotEqual(0, result.exit_code)
        self.assertIn("not supported", result.reason or "")

    def test_repeated_config_runs_reset_the_default_solver(self):
        common = {
            "format": "verilog",
            "verification": "sec",
            "input_paths": [
                [str(self.reference)],
                [str(self.equivalent)],
            ],
            "sec_encoding": "binary",
            "max_k": 2,
        }
        first_config = self.root / "cadical.json"
        first_log = self.root / "cadical.log"
        first_config.write_text(
            json.dumps({**common, "solver": "cadical", "log_file": str(first_log)}),
            encoding="utf-8",
        )
        second_config = self.root / "default-solver.json"
        second_log = self.root / "default-solver.log"
        second_config.write_text(
            json.dumps({**common, "log_file": str(second_log)}),
            encoding="utf-8",
        )

        self.assertEqual(
            VerificationStatus.EQUIVALENT,
            run_config(first_config).status,
        )
        self.assertEqual(
            VerificationStatus.EQUIVALENT,
            run_config(second_config).status,
        )
        self.assertIn("Solver: KISSAT", second_log.read_text(encoding="utf-8"))

    def test_native_failures_return_a_reason(self):
        result = verify(self.root / "missing.v", self.equivalent)
        self.assertEqual(VerificationStatus.ERROR, result.status)
        self.assertIsNotNone(result.reason)
        self.assertTrue(result.reason)

    def test_validation(self):
        with self.assertRaisesRegex(ValueError, "design1"):
            verify(Design(), self.equivalent)
        with self.assertRaisesRegex(ValueError, "SEC"):
            verify(
                self.reference,
                self.equivalent,
                options=VerificationOptions(max_k=2),
            )
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            verify("", self.equivalent)
        with self.assertRaisesRegex(TypeError, "max_k"):
            verify(
                self.reference,
                self.equivalent,
                options=VerificationOptions(
                    mode=VerificationMode.SEC,
                    max_k=True,
                ),
            )
        with self.assertRaisesRegex(TypeError, "compact"):
            verify(
                self.reference,
                self.equivalent,
                options=VerificationOptions(compact="false"),
            )
        with self.assertRaisesRegex(TypeError, "sequence"):
            run_cli("--help")
        with self.assertRaisesRegex(ValueError, "path"):
            run_config("")


if __name__ == "__main__":
    unittest.main()
