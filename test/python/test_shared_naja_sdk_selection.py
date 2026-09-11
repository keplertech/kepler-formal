# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


_MODULE = Path(__file__).resolve().parents[2] / "src/python/KeplerSharedNaja.cmake"


@unittest.skipUnless(shutil.which("cmake"), "CMake is required")
class SharedNajaSDKSelectionTests(unittest.TestCase):
    def test_python_provider_overrides_cached_or_user_sdk_directory(self):
        self._check_provider_selection(windows_separators=False)

    def test_windows_sdk_path_through_package_wrapper(self):
        self._check_provider_selection(windows_separators=True)

    def _check_provider_selection(self, *, windows_separators):
        with tempfile.TemporaryDirectory(prefix="kepler sdk selection ") as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            wrong = root / "wrong sdk"
            wrong.mkdir()
            (wrong / "NajaEDAConfig.cmake").write_text(
                'message(FATAL_ERROR "Loaded a foreign cached SDK")\n', encoding="utf-8")
            (source / "FindTBB.cmake").write_text("set(TBB_FOUND TRUE)\n", encoding="utf-8")
            (source / "boost").mkdir()
            (source / "boost/version.hpp").touch()
            # vcpkg wraps find_package in a macro and expands ARGN again in
            # set(). Native Windows paths can then turn '\\U' into an invalid
            # CMake escape, even when the original call quoted its argument.
            package_wrapper = (
                "if(POLICY CMP0219)\n"
                "  cmake_policy(SET CMP0219 OLD)\n"
                "endif()\n"
                "macro(find_package)\n"
                '  set(wrapper_arguments "${ARGN}")\n'
                "  _find_package(${wrapper_arguments})\n"
                "endmacro()\n"
            ) if windows_separators else ""
            (source / "CMakeLists.txt").write_text(
                "cmake_minimum_required(VERSION 3.30)\n"
                "project(SharedSDKSelection NONE)\n"
                'set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")\n'
                + package_wrapper +
                f"include([==[{_MODULE.as_posix()}]==])\n"
                'if(NOT NajaEDA_CONFIG STREQUAL "${provider_cmake_dir}/NajaEDAConfig.cmake")\n'
                '  message(FATAL_ERROR "Selected SDK does not belong to Python provider")\n'
                "endif()\n", encoding="utf-8")
            for provider_name in ("provider one", "provider two"):
                with self.subTest(provider=provider_name):
                    provider = root / "Users" / provider_name
                    package = provider / "najaeda"
                    config = package / "sdk/cmake"
                    config.mkdir(parents=True)
                    (package / "__init__.py").touch()
                    (package / "sdk.py").write_text(
                        "from pathlib import Path\n"
                        "def get_cmake_dir():\n"
                        "    path = str(Path(__file__).resolve().parent / 'sdk/cmake')\n"
                        + ("    return path.replace('/', chr(92))\n"
                           if windows_separators else "    return path\n"),
                        encoding="utf-8")
                    (config / "NajaEDAConfig.cmake").write_text(
                        "set(NajaEDA_FOUND TRUE)\n", encoding="utf-8")
                    environment = os.environ.copy()
                    environment.update(PYTHONPATH=str(provider), PYTHONDONTWRITEBYTECODE="1")
                    result = subprocess.run([
                        shutil.which("cmake"), "-S", str(source), "-B", str(root / "build"),
                        f"-DPython3_EXECUTABLE={sys.executable}",
                        f"-DBoost_INCLUDE_DIR={source}",
                        # Keep a real stale/user-selected cache entry while
                        # changing the Python provider between configurations.
                        f"-DNajaEDA_DIR={wrong}",
                    ], env=environment, capture_output=True, text=True)
                    self.assertEqual(0, result.returncode, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
