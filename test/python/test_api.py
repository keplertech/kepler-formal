# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

import os
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

    def test_skipped_net_diagnostics_preserve_borrowed_nets(self):
        conflicting_nets = []
        for design in (self.reference, self.equivalent):
            conflicting = najaeda.naja.SNLScalarNet.create(design, "conflicting_net")
            for name in ("driver0", "driver1"):
                driver = najaeda.naja.SNLScalarTerm.create(
                    design, najaeda.naja.SNLTerm.Direction.Input, name
                )
                driver.setNet(conflicting)
            output = najaeda.naja.SNLScalarTerm.create(
                design, najaeda.naja.SNLTerm.Direction.Output, "conflicting_output"
            )
            output.setNet(conflicting)
            conflicting_nets.append(conflicting)
        top_before = self.universe.getTopDesign()
        previous_directory = Path.cwd()
        try:
            os.chdir(self.root)
            result = verify_designs(
                self.reference, self.equivalent,
                options=VerificationOptions(
                    log_file=self.root / "skipped-nets.log",
                    report_skipped_outputs=True,
                ),
            )
        finally:
            os.chdir(previous_directory)

        self.assertEqual(VerificationStatus.EQUIVALENT, result.status)
        self.assertTrue((self.root / "skipped_multi_driver_pos.txt").is_file())
        self.assertIs(self.universe.getTopDesign(), top_before)
        for design, conflicting in zip((self.reference, self.equivalent), conflicting_nets):
            self.assertEqual("conflicting_net", conflicting.getName())
            self.assertIs(design.getScalarTerm("conflicting_output").getNet(), conflicting)
        repeated = verify_designs(
            self.reference, self.equivalent, options=self._options("after-report.log")
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, repeated.status)

    def test_sec_internal_relation_switches(self):
        for learn in (False, True):
            for allow_x in (False, True):
                with self.subTest(learn=learn, allow_x=allow_x):
                    result = verify_designs(
                        self.reference, self.equivalent,
                        options=VerificationOptions(
                            mode=VerificationMode.SEC,
                            learn_internal_relations=learn,
                            allow_x_equality_in_internal_relations=allow_x,
                            log_file=self.root / f"relations-{learn}-{allow_x}.log"))
                    self.assertEqual(VerificationStatus.EQUIVALENT, result.status)

    def test_internal_relation_option_validation(self):
        invalid_options = (
            (TypeError, "learn_internal_relations", VerificationOptions(
                mode=VerificationMode.SEC, learn_internal_relations="false")),
            (TypeError, "allow_x_equality", VerificationOptions(
                mode=VerificationMode.SEC, allow_x_equality_in_internal_relations=1)),
            (ValueError, "SEC", VerificationOptions(learn_internal_relations=False)),
        )
        for error_type, message, options in invalid_options:
            with self.subTest(options=options):
                with self.assertRaisesRegex(error_type, message):
                    verify_designs(self.reference, self.equivalent, options=options)

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
            (TypeError, "set_as_boundary", VerificationOptions(
                set_as_boundary="opaque")),
            (TypeError, r"set_as_boundary\[0\]", VerificationOptions(
                set_as_boundary=["opaque"])),
            (ValueError, "exactly two", VerificationOptions(
                set_as_boundary=[("opaque",)])),
            (TypeError, r"set_as_boundary\[0\]\[1\]", VerificationOptions(
                set_as_boundary=[("opaque", 1)])),
            (ValueError, "must not be empty", VerificationOptions(
                set_as_boundary=[("", "opaque")])),
            (ValueError, "NUL", VerificationOptions(
                set_as_boundary=[("opaque\0child", "opaque")])),
            (ValueError, "log_file", VerificationOptions(log_file="")),
        )
        for error_type, message, options in invalid_options:
            with self.subTest(options=options):
                with self.assertRaisesRegex(error_type, message):
                    verify_designs(self.reference, self.equivalent, options=options)
        with self.assertRaisesRegex(TypeError, "NativeDesign.*SNLDesign"):
            verify_designs("reference.v", self.equivalent)

    def test_verification_options_preserves_positional_compatibility(self):
        options = VerificationOptions(
            VerificationMode.LEC,
            "kissat",
            None,
            None,
            None,
            True,
            True,
            "verification.log",
            "debug",
        )
        self.assertTrue(options.allow_boundary_mismatch)
        self.assertTrue(options.report_skipped_outputs)
        self.assertEqual("verification.log", options.log_file)
        self.assertEqual("debug", options.log_level)
        self.assertEqual((), options.set_as_boundary)


if __name__ == "__main__":
    unittest.main()
