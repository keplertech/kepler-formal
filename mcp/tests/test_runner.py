from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kepler_formal_mcp.models import BinaryResolution, CheckKind, ResultStatus
from kepler_formal_mcp.resolver import BinarySetupError
from kepler_formal_mcp.runner import CheckRunner, classify_verdict


class StubResolver:
    def __init__(self, binary: Path):
        self.binary = binary

    def resolve(self) -> BinaryResolution:
        return BinaryResolution(self.binary, "path")


class MissingResolver:
    def resolve(self) -> BinaryResolution:
        raise BinarySetupError("binary not found after PATH, cache, and release lookup")


class VerdictTests(unittest.TestCase):
    def test_lec_difference_is_a_normal_check_failure_despite_zero_exit(self):
        status, _ = classify_verdict(
            CheckKind.GATE_LEC, 0, "Difference was found. Please refer to the log", ""
        )
        self.assertEqual(status, ResultStatus.CHECK_FAILED)

    def test_sec_stable_codes_are_distinct(self):
        self.assertEqual(
            classify_verdict(CheckKind.GATE_SEC, 0, "", "")[0],
            ResultStatus.PASSED,
        )
        self.assertEqual(
            classify_verdict(CheckKind.GATE_SEC, 2, "", "")[0],
            ResultStatus.INCONCLUSIVE,
        )
        self.assertEqual(
            classify_verdict(CheckKind.GATE_SEC, 3, "", "")[0],
            ResultStatus.CHECK_FAILED,
        )

    def test_generic_exit_one_is_a_crash_not_a_partial_proof(self):
        status, _ = classify_verdict(
            CheckKind.RTL_SEC, 1, "", "Netlist loading failed: bad input"
        )
        self.assertEqual(status, ResultStatus.CRASH)

    def test_tool_mode_mismatch_is_a_wrapper_crash(self):
        status, summary = classify_verdict(
            CheckKind.GATE_SEC,
            0,
            "Verification: lec\nNo difference was found.",
            "",
        )
        self.assertEqual(status, ResultStatus.CRASH)
        self.assertIn("expected verification mode sec", summary)


class RunnerTests(unittest.IsolatedAsyncioTestCase):
    async def test_binary_resolution_failure_is_a_setup_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "case.yaml"
            config.write_text("verification: sec\n")
            result = await CheckRunner(root, MissingResolver()).run(
                CheckKind.GATE_SEC, "case.yaml", 1
            )
            self.assertEqual(result.status, "setup_error")
            self.assertEqual(result.error_kind, "setup_error")

    async def test_non_finite_timeout_stays_json_safe(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "case.yaml"
            config.write_text("verification: sec\n")
            result = await CheckRunner(root, MissingResolver()).run(
                CheckKind.GATE_SEC, "case.yaml", float("inf")
            )
            self.assertEqual(result.status, "crash")
            self.assertEqual(result.timeout_seconds, 0.0)

    async def test_relative_config_is_resolved_against_project_root(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "case.yaml"
            config.write_text("verification: lec\n")
            binary = root / "kepler-formal"
            binary.write_text("#!/bin/sh\necho 'No difference was found.'\n")
            binary.chmod(0o755)
            runner = CheckRunner(root, StubResolver(binary))
            result = await runner.run(CheckKind.GATE_LEC, "case.yaml", 2)
            self.assertEqual(result.status, "passed")
            self.assertEqual(result.config_path, str(config.resolve()))
            self.assertEqual(
                result.command, (str(binary), "--config", str(config.resolve()))
            )

    async def test_timeout_kills_the_process_and_reports_progress(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "case.yaml"
            config.write_text("verification: sec\n")
            binary = root / "kepler-formal"
            binary.write_text("#!/bin/sh\nsleep 5\n")
            binary.chmod(0o755)
            progress = []

            async def report(current, total, message):
                progress.append((current, total, message))

            runner = CheckRunner(
                root, StubResolver(binary), progress_interval_seconds=0.01
            )
            result = await runner.run(CheckKind.RTL_SEC, str(config), 0.05, report)
            self.assertEqual(result.status, "timeout")
            self.assertEqual(result.error_kind, "timeout")
            self.assertTrue(progress)


if __name__ == "__main__":
    unittest.main()
