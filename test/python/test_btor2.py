# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

import json
import tempfile
import unittest
from pathlib import Path

from kepler_formal import VerificationStatus, run_cli, run_config


class Btor2ExportTest(unittest.TestCase):
    def test_dump_only_returns_exported_without_claiming_equivalence(self):
        with tempfile.TemporaryDirectory(prefix="kepler_python_btor2_") as directory:
            root = Path(directory)
            reference = root / "reference.v"
            different = root / "different.v"
            reference.write_text(
                "module top(input a, input b, output y); or (y, a, b); endmodule\n",
                encoding="utf-8",
            )
            different.write_text(
                "module top(input a, input b, output y); and (y, a, b); endmodule\n",
                encoding="utf-8",
            )
            exported = root / "problem.btor2"
            result = run_cli((
                "-verilog", "-v", "sec", "--sec-encoding", "binary",
                reference, different, "--dump-btor2", exported, "--dump-only",
            ))
            self.assertEqual(VerificationStatus.EXPORTED, result.status)
            self.assertEqual(0, result.exit_code)
            self.assertFalse(result.equivalent)
            self.assertFalse(result.conclusive)
            self.assertEqual(1, result.covered_outputs)
            self.assertEqual(1, result.total_outputs)
            self.assertEqual(0, result.proven_outputs)
            self.assertIn("proof not run", result.reason)
            self.assertIn(" bad ", exported.read_text(encoding="utf-8"))

            # Existing YAML/JSON callers use the same shared driver. Repeating
            # a run must release native state and return a fresh proof result.
            config = root / "verify.json"
            config.write_text(json.dumps({
                "format": "verilog",
                "verification": "sec",
                "sec_encoding": "binary",
                "input_paths": [[str(reference)], [str(different)]],
                "btor2_export": True,
                "btor2_export_path": str(exported),
                "dump_only": False,
                "log_file": str(root / "verify.log"),
            }), encoding="utf-8")
            solved = run_config(config)
            self.assertEqual(VerificationStatus.DIFFERENT, solved.status)
            self.assertEqual(3, solved.exit_code)
            self.assertEqual(0, solved.proven_outputs)
            self.assertIn(" bad ", exported.read_text(encoding="utf-8"))
