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
_TBB_FINDER = """
set(TBB_FOUND TRUE)
foreach(dependency IN ITEMS tbb tbbmalloc)
  if(NOT TARGET TBB::${dependency})
    add_library(TBB::${dependency} SHARED IMPORTED)
    foreach(suffix IN ITEMS "" _RELEASE _DEBUG _RELWITHDEBINFO _MINSIZEREL)
      set_property(TARGET TBB::${dependency} PROPERTY IMPORTED_LOCATION${suffix}
        "${CMAKE_CURRENT_LIST_DIR}/${dependency}-provider-hash.dll")
      set_property(TARGET TBB::${dependency} PROPERTY IMPORTED_IMPLIB${suffix}
        "${CMAKE_CURRENT_LIST_DIR}/${dependency}-provider-hash.lib")
    endforeach()
  endif()
endforeach()
"""


@unittest.skipUnless(shutil.which("cmake"), "CMake is required")
class SharedNajaSDKSelectionTests(unittest.TestCase):
    def test_python_provider_overrides_cached_or_user_sdk_directory(self):
        self._check_provider_selection(windows_separators=False)

    def test_windows_sdk_path_through_package_wrapper(self):
        self._check_provider_selection(windows_separators=True)

    @unittest.skipUnless(
        any(shutil.which(compiler) for compiler in
            (os.environ.get("CXX", ""), "c++", "g++", "clang++", "clang-cl", "cl")
            if compiler), "A C++ compiler is required")
    def test_tbb_auto_link_policy_propagates_only_to_windows_consumers(self):
        with tempfile.TemporaryDirectory(prefix="kepler tbb policy ") as temporary:
            root = Path(temporary)
            package = root / "provider/najaeda"
            config = package / "sdk/cmake"
            config.mkdir(parents=True)
            (package / "__init__.py").touch()
            (package / "sdk.py").write_text(
                "from pathlib import Path\n"
                "def get_cmake_dir():\n"
                "    return str(Path(__file__).resolve().parent / 'sdk/cmake')\n",
                encoding="utf-8")
            (config / "NajaEDAConfig.cmake").write_text(
                "set(NajaEDA_FOUND TRUE)\n", encoding="utf-8")
            source = root / "source"
            (source / "boost").mkdir(parents=True)
            (source / "boost/version.hpp").touch()
            (source / "FindTBB.cmake").write_text(_TBB_FINDER, encoding="utf-8")
            (source / "probe.cpp").write_text("""
#if EXPECT_TBB_NO_IMPLICIT_LINKAGE
#if !defined(__TBB_NO_IMPLICIT_LINKAGE) || __TBB_NO_IMPLICIT_LINKAGE != 1
#error "Windows TBB consumers must disable implicit linkage"
#endif
#elif defined(__TBB_NO_IMPLICIT_LINKAGE)
#error "TBB linkage policy leaked to an unrelated or non-Windows target"
#endif
int tbb_linkage_probe() { return 0; }
""", encoding="utf-8")
            (source / "CMakeLists.txt").write_text(
                "cmake_minimum_required(VERSION 3.30)\n"
                "project(TBBAutoLinkPolicy LANGUAGES CXX)\n"
                'set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")\n'
                # Emulate only the policy branch, retaining the host compiler
                # and platform for the actual object/static-library builds.
                'set(host_is_windows "${WIN32}")\n'
                'set(WIN32 "${TEST_WINDOWS_POLICY}")\n'
                f"include([==[{_MODULE.as_posix()}]==])\n"
                'set(WIN32 "${host_is_windows}")\n'
                """
foreach(dependency IN ITEMS tbb tbbmalloc)
  foreach(suffix IN ITEMS "" _RELEASE _DEBUG _RELWITHDEBINFO _MINSIZEREL)
    get_target_property(implib TBB::${dependency} IMPORTED_IMPLIB${suffix})
    if(NOT implib STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}/${dependency}-provider-hash.lib")
      message(FATAL_ERROR "Provider import library changed: ${implib}")
    endif()
  endforeach()
  add_library(${dependency}_dependency STATIC probe.cpp)
  target_link_libraries(${dependency}_dependency PUBLIC TBB::${dependency})
  add_library(${dependency}_consumer OBJECT probe.cpp)
  target_link_libraries(${dependency}_consumer PRIVATE ${dependency}_dependency)
  foreach(target IN ITEMS ${dependency}_dependency ${dependency}_consumer)
    target_compile_definitions(${target} PRIVATE
      EXPECT_TBB_NO_IMPLICIT_LINKAGE=$<BOOL:${TEST_WINDOWS_POLICY}>)
  endforeach()
endforeach()
add_library(unrelated OBJECT probe.cpp)
target_compile_definitions(unrelated PRIVATE EXPECT_TBB_NO_IMPLICIT_LINKAGE=0)
""", encoding="utf-8")
            environment = os.environ.copy()
            environment.update(PYTHONPATH=str(package.parent), PYTHONDONTWRITEBYTECODE="1")
            for windows_policy in (False, True):
                with self.subTest(windows_policy=windows_policy):
                    build = root / ("windows" if windows_policy else "non-windows")
                    commands = ([
                        shutil.which("cmake"), "-S", str(source), "-B", str(build),
                        f"-DPython3_EXECUTABLE={sys.executable}",
                        f"-DBoost_INCLUDE_DIR={source}",
                        f"-DTEST_WINDOWS_POLICY={'ON' if windows_policy else 'OFF'}",
                    ], [shutil.which("cmake"), "--build", str(build),
                        "--config", "Release", "--parallel", "2"])
                    for command in commands:
                        result = subprocess.run(command, env=environment,
                                                capture_output=True, text=True, timeout=120)
                        self.assertEqual(0, result.returncode, result.stdout + result.stderr)

    def _check_provider_selection(self, *, windows_separators):
        with tempfile.TemporaryDirectory(prefix="kepler sdk selection ") as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            wrong = root / "wrong sdk"
            wrong.mkdir()
            (wrong / "NajaEDAConfig.cmake").write_text(
                'message(FATAL_ERROR "Loaded a foreign cached SDK")\n', encoding="utf-8")
            (source / "FindTBB.cmake").write_text(_TBB_FINDER, encoding="utf-8")
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
