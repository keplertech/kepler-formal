# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "shared_naja_wheels", Path(__file__).resolve().parents[2] / "ci/shared_naja_wheels.py")
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)


class SharedWheelPackagingTests(unittest.TestCase):
    def test_release_uses_published_provider_instead_of_rebuilding_it(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.object(helper, "PROVIDER_REQUIREMENT", "najaeda==1.2.3"), \
                 patch.object(helper, "run") as run, \
                 patch.object(helper, "repair") as repair:
                helper.build_provider(root)
            commands = [call.args for call in run.call_args_list]
            downloads = [command for command in commands if "download" in command]
            self.assertEqual(1, len(downloads))
            self.assertIn("--only-binary=:all:", downloads[0])
            self.assertIn("najaeda==1.2.3", downloads[0])
            self.assertFalse(any("--wheel-dir" in command for command in commands))
            repair.assert_not_called()

    def test_linux_excludes_provider_closure_including_allocator(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dependencies = {name: root / name for name in
                            ("libnaja_nl-123.so", "libtbb-456.so.12")}
            with patch.object(helper.platform, "system", return_value="Linux"), \
                 patch.object(helper.platform, "machine", return_value="aarch64"), \
                 patch.object(helper, "run") as run:
                helper.repair_external(root / "input.whl", root / "out", dependencies)
            command = run.call_args.args
            self.assertIn("manylinux_2_28_aarch64", command)
            for name in dependencies:
                self.assertEqual(command[command.index(name) - 1], "--exclude")

    def test_macos_excludes_provider_without_ignoring_missing_libraries(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.object(helper.platform, "system", return_value="Darwin"), \
                 patch.object(helper, "run") as run:
                helper.repair_external(root / "input.whl", root / "out",
                                       {"libnaja_nl.dylib": root / "libnaja_nl.dylib"})
            command = run.call_args.args
            self.assertIn("--exclude", command)
            self.assertIn("libnaja_nl.dylib", command)
            self.assertNotIn("--ignore-missing-dependencies", command)

    def test_windows_keeps_provider_import_library_names(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.object(helper.platform, "system", return_value="Windows"), \
                 patch.dict(os.environ, {"USERPROFILE": str(root)}), \
                 patch.object(helper, "run") as run:
                helper.repair_external(root / "provider.whl", root / "out", {})
            self.assertIn("--no-mangle", run.call_args.args)
            self.assertIn("naja_*.dll;libnaja_*.dll", run.call_args.args)

    def test_windows_consumer_excludes_provider_instead_of_copying(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with patch.object(helper.platform, "system", return_value="Windows"), \
                 patch.dict(os.environ, {"USERPROFILE": str(root)}), \
                 patch.object(helper, "run") as run:
                helper.repair_external(root / "consumer.whl", root / "out",
                                       {"naja_nl.dll": root / "naja_nl.dll"})
            self.assertIn("--exclude", run.call_args.args)
            self.assertNotIn("--no-mangle", run.call_args.args)


if __name__ == "__main__":
    unittest.main()
