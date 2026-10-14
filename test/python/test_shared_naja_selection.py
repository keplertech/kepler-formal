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
_RUNTIMES = ("naja_nl", "naja_dnl", "naja_bne", "naja_opt", "naja_metrics", "naja_python")
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


def _make_provider(site: Path, version: str = "0.7.28", git_hash: str = "abc1234") -> Path:
    """Create an importable fake NajaEDA development wheel installed in ``site``."""
    package = site / "najaeda"
    (package / "sdk/include/naja/core").mkdir(parents=True)
    (package / "__init__.py").touch()
    (package / "naja.py").write_text(
        f"def getVersion():\n    return {version!r}\n"
        f"def getGitHash():\n    return {git_hash!r}\n", encoding="utf-8")
    for name in ("NajaVersion.h", "NLUniverse.h"):
        (package / "sdk/include/naja/core" / name).touch()
    windows = sys.platform == "win32"
    for name in _RUNTIMES:
        (package / (f"{name}.dll" if windows else f"lib{name}.so")).write_bytes(b"provider " + name.encode())
        if windows:
            (package / "sdk/lib").mkdir(exist_ok=True)
            (package / "sdk/lib" / f"{name}.lib").touch()
    info = site / f"najaeda-{version}.dist-info"
    info.mkdir()
    (info / "METADATA").write_text(
        f"Metadata-Version: 2.1\nName: najaeda\nVersion: {version}\n", encoding="utf-8")
    return package


def _write_fixture(source: Path, body: str, *, languages: str = "NONE") -> None:
    (source / "FindTBB.cmake").write_text(_TBB_FINDER, encoding="utf-8")
    (source / "boost").mkdir(exist_ok=True)
    (source / "boost/version.hpp").touch()
    (source / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.30)\n"
        f"project(SharedProviderSelection LANGUAGES {languages})\n"
        'set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")\n'
        "add_library(Python3::Module INTERFACE IMPORTED)\n"
        f"include([==[{_MODULE.as_posix()}]==])\n" + body, encoding="utf-8")


@unittest.skipUnless(shutil.which("cmake"), "CMake is required")
class SharedNajaSelectionTests(unittest.TestCase):
    def test_build_links_the_provider_installed_in_the_build_interpreter(self):
        with tempfile.TemporaryDirectory(prefix="kepler provider selection ") as temporary:
            root = Path(temporary).resolve()
            source = root / "source"
            source.mkdir()
            _write_fixture(source, """
string(FIND "${NajaEDA_naja_nl_LIBRARY}" "${EXPECTED_PROVIDER}/" position)
if(NOT position EQUAL 0)
  message(FATAL_ERROR "Linked ${NajaEDA_naja_nl_LIBRARY} instead of ${EXPECTED_PROVIDER}")
endif()
get_target_property(includes naja_nl INTERFACE_INCLUDE_DIRECTORIES)
if(NOT "${EXPECTED_PROVIDER}/sdk/include/naja/core" IN_LIST includes)
  message(FATAL_ERROR "Provider headers are not used: ${includes}")
endif()
""")
            # Reconfigure the same build directory with another provider: the
            # interpreter's current package wins, not a cached selection.
            for site_name in ("provider one", "provider two"):
                with self.subTest(provider=site_name):
                    site = root / "Users" / site_name
                    package = _make_provider(site)
                    environment = os.environ.copy()
                    environment.update(PYTHONPATH=str(site), PYTHONDONTWRITEBYTECODE="1")
                    result = subprocess.run([
                        shutil.which("cmake"), "-S", str(source), "-B", str(root / "build"),
                        f"-DPython3_EXECUTABLE={sys.executable}",
                        f"-DBoost_INCLUDE_DIR={source}",
                        f"-DEXPECTED_PROVIDER={package.as_posix()}",
                    ], env=environment, capture_output=True, text=True)
                    self.assertEqual(0, result.returncode, result.stdout + result.stderr)

    def test_published_release_is_rejected_without_the_option(self):
        with tempfile.TemporaryDirectory(prefix="kepler provider selection ") as temporary:
            root = Path(temporary).resolve()
            source = root / "source"
            source.mkdir()
            _write_fixture(source, "")
            site = root / "site"
            package = _make_provider(site, version="0.7.24", git_hash="2263958")
            shutil.rmtree(package / "sdk")
            environment = os.environ.copy()
            environment.update(PYTHONPATH=str(site), PYTHONDONTWRITEBYTECODE="1")
            result = subprocess.run([
                shutil.which("cmake"), "-S", str(source), "-B", str(root / "build"),
                f"-DPython3_EXECUTABLE={sys.executable}",
                f"-DBoost_INCLUDE_DIR={source}",
            ], env=environment, capture_output=True, text=True)
            self.assertNotEqual(0, result.returncode)
            self.assertIn("ships no headers", result.stderr)
            self.assertIn("KEPLER_USE_PUBLISHED_NAJAEDA", result.stderr)

    @unittest.skipUnless(
        any(shutil.which(compiler) for compiler in
            (os.environ.get("CXX", ""), "c++", "g++", "clang++", "clang-cl", "cl")
            if compiler), "A C++ compiler is required")
    def test_tbb_auto_link_policy_propagates_only_to_windows_consumers(self):
        with tempfile.TemporaryDirectory(prefix="kepler tbb policy ") as temporary:
            root = Path(temporary).resolve()
            site = root / "provider"
            _make_provider(site)
            source = root / "source"
            source.mkdir()
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
            (source / "FindTBB.cmake").write_text(_TBB_FINDER, encoding="utf-8")
            (source / "boost").mkdir()
            (source / "boost/version.hpp").touch()
            (source / "CMakeLists.txt").write_text(
                "cmake_minimum_required(VERSION 3.30)\n"
                "project(TBBAutoLinkPolicy LANGUAGES CXX)\n"
                'set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")\n'
                "add_library(Python3::Module INTERFACE IMPORTED)\n"
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
            environment.update(PYTHONPATH=str(site), PYTHONDONTWRITEBYTECODE="1")
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


if __name__ == "__main__":
    unittest.main()
