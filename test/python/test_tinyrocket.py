# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Exercise the real NajaEDA-to-Kepler TinyRocket example in fresh processes."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


EXAMPLE = Path(__file__).resolve().parents[2] / "examples" / "tinyrocket"
SCRIPT = EXAMPLE / "verify_python.py"


def _tail(text):
    if isinstance(text, bytes):
        text = text.decode("utf-8", errors="replace")
    return "\n".join((text or "").splitlines()[-80:])


class TinyRocketPythonTest(unittest.TestCase):
    def _run_example(self, arguments, expected_status, expected_exit):
        with tempfile.TemporaryDirectory(prefix="kepler_formal_tinyrocket_") as work:
            log = Path(work) / "tinyrocket_python.log"

            def diagnostics(stdout, stderr):
                log_text = (
                    log.read_text(encoding="utf-8", errors="replace")
                    if log.is_file() else "(no log produced)"
                )
                return (
                    f"\nstdout (last 80 lines):\n{_tail(stdout)}"
                    f"\nstderr (last 80 lines):\n{_tail(stderr)}"
                    f"\nKepler log (last 80 lines):\n{_tail(log_text)}"
                )

            try:
                # Inherit PYTHONPATH so CTest can supply its locally staged
                # package, while wheel tests use their installed packages.
                completed = subprocess.run(
                    [sys.executable, "-X", "faulthandler", str(SCRIPT), *arguments],
                    cwd=work,
                    capture_output=True,
                    text=True,
                    errors="replace",
                    timeout=120,
                )
            except subprocess.TimeoutExpired as error:
                self.fail(
                    "TinyRocket Python example timed out after 120 seconds"
                    + diagnostics(error.stdout, error.stderr)
                )

            details = diagnostics(completed.stdout, completed.stderr)
            self.assertEqual(expected_exit, completed.returncode, details)
            self.assertIn(f"Result: {expected_status}\n", completed.stdout, details)
            self.assertIn(
                "Designs still available in NajaEDA:", completed.stdout, details
            )
            self.assertTrue(log.is_file(), details)
            self.assertGreater(log.stat().st_size, 0, details)

    def test_default_edited_design_is_different(self):
        # No design arguments: cover the example's default paths from a
        # working directory outside the repository.
        self._run_example([], "different", 1)

    def test_two_original_designs_are_equivalent(self):
        original = str(EXAMPLE / "tinyrocket.v")
        self._run_example([original, original], "equivalent", 0)


if __name__ == "__main__":
    unittest.main()
