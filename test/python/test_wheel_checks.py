# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

from __future__ import annotations

import importlib.machinery
import importlib.util
import struct
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


CHECKER_PATH = Path(__file__).resolve().parents[2] / "ci" / "check_kepler_wheel.py"
CHECKER_SPEC = importlib.util.spec_from_file_location("check_kepler_wheel", CHECKER_PATH)
if CHECKER_SPEC is None or CHECKER_SPEC.loader is None:
    raise RuntimeError(f"cannot load wheel checker from {CHECKER_PATH}")
wheel_check = importlib.util.module_from_spec(CHECKER_SPEC)
CHECKER_SPEC.loader.exec_module(wheel_check)


def _write_pe(path: Path, imports: tuple[str, ...]) -> None:
    data = bytearray(0x600)
    pe_offset = 0x80
    optional_offset = pe_offset + 24
    optional_size = 0xF0
    section_offset = optional_offset + optional_size
    raw_offset = 0x200
    section_rva = 0x1000

    data[:2] = b"MZ"
    struct.pack_into("<I", data, 0x3C, pe_offset)
    data[pe_offset : pe_offset + 4] = b"PE\0\0"
    struct.pack_into(
        "<HHIIIHH",
        data,
        pe_offset + 4,
        0x8664,
        1,
        0,
        0,
        0,
        optional_size,
        0x2022,
    )
    struct.pack_into("<H", data, optional_offset, 0x20B)
    struct.pack_into("<I", data, optional_offset + 108, 16)
    struct.pack_into(
        "<II",
        data,
        optional_offset + 120,
        section_rva,
        (len(imports) + 1) * 20,
    )
    data[section_offset : section_offset + 8] = b".idata\0\0"
    struct.pack_into(
        "<IIII",
        data,
        section_offset + 8,
        0x400,
        section_rva,
        0x400,
        raw_offset,
    )

    name_offset = raw_offset + (len(imports) + 1) * 20
    for index, dependency in enumerate(imports):
        encoded = dependency.encode("ascii") + b"\0"
        name_rva = section_rva + name_offset - raw_offset
        struct.pack_into(
            "<IIIII", data, raw_offset + index * 20, 0, 0, 0, name_rva, 0
        )
        data[name_offset : name_offset + len(encoded)] = encoded
        name_offset += len(encoded)
    path.write_bytes(data)


class FakeDistribution:
    def __init__(self, root: Path, files: tuple[str, ...]):
        self.root = root
        self.files = tuple(Path(item) for item in files)

    def locate_file(self, item: Path) -> Path:
        return self.root / item


class WheelCheckerTest(unittest.TestCase):
    def test_native_files_include_unix_and_windows_formats(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            distribution = FakeDistribution(
                root,
                (
                    "kepler_formal/_native.so",
                    "kepler_formal/libnaja_nl.so.1",
                    "kepler_formal/libnaja_nl.dylib",
                    "kepler_formal/_native.PYD",
                    "kepler_formal/naja_nl.DLL",
                    "kepler_formal/naja_nl.lib",
                ),
            )
            native_names = {
                path.name for path in wheel_check._native_files(distribution)
            }
        self.assertEqual(
            {
                "_native.so",
                "libnaja_nl.so.1",
                "libnaja_nl.dylib",
                "_native.PYD",
                "naja_nl.DLL",
            },
            native_names,
        )

    def test_pe_import_reader_and_python_runtime_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            native_file = Path(temporary) / "naja.pyd"
            expected = (
                "python314t.dll",
                "libkepler_najaeda_naja_nl.dll",
                "KERNEL32.dll",
            )
            _write_pe(native_file, expected)
            self.assertEqual(expected, wheel_check._read_pe_imports(native_file))
            wheel_check._check_windows_python_imports(
                native_file,
                expected,
                version_info=(3, 14),
                gil_disabled=True,
            )

            runtime_version = (
                wheel_check.sys.version_info.major,
                wheel_check.sys.version_info.minor,
            )
            runtime_dll = wheel_check._expected_windows_python_dll(
                runtime_version,
                bool(wheel_check.sysconfig.get_config_var("Py_GIL_DISABLED")),
            )
            _write_pe(native_file, (runtime_dll, "KERNEL32.dll"))
            with patch.object(wheel_check.platform, "system", return_value="Windows"):
                self.assertEqual(
                    (runtime_dll, "KERNEL32.dll"),
                    wheel_check._check_linkage((native_file,))[native_file],
                )

            with self.assertRaisesRegex(RuntimeError, "expected python314.dll"):
                wheel_check._check_windows_python_imports(
                    native_file,
                    expected,
                    version_info=(3, 14),
                    gil_disabled=False,
                )

    def test_pe_import_reader_rejects_non_pe_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            native_file = Path(temporary) / "broken.dll"
            native_file.write_bytes(b"not a portable executable")
            with self.assertRaisesRegex(RuntimeError, "not a PE image"):
                wheel_check._read_pe_imports(native_file)

    def test_supported_platform_tags(self):
        cases = (
            ("Linux", "x86_64", "manylinux_2_28_x86_64"),
            ("Linux", "aarch64", "manylinux_2_28_aarch64"),
            ("Darwin", "arm64", "macosx_11_0_arm64"),
            ("Windows", "AMD64", "win_amd64"),
        )
        for system, machine, expected in cases:
            with self.subTest(system=system, machine=machine):
                self.assertEqual(
                    expected,
                    wheel_check._expected_platform_tag(system, machine),
                )
                wheel_check._check_platform_tag(
                    (f"cp315-cp315-{expected}",),
                    system=system,
                    machine=machine,
                )

    def test_library_names_normalize_across_platforms(self):
        cases = {
            "libnaja_nl.so.1": "naja_nl",
            "@rpath/libnaja_nl.2.dylib": "naja_nl",
            "naja_nl.dll": "naja_nl",
            "LIBNAJA_NL.DLL": "naja_nl",
            "libkepler_najaeda_naja_nl.dll": "kepler_najaeda_naja_nl",
        }
        for library, expected in cases.items():
            with self.subTest(library=library):
                self.assertEqual(
                    expected,
                    wheel_check._canonical_library_name(library),
                )

    def test_supported_cpython_tags_include_free_threaded_314(self):
        cases = (
            ((3, 10), False, ("cp310", "cp310")),
            ((3, 15), False, ("cp315", "cp315")),
            ((3, 14), True, ("cp314", "cp314t")),
        )
        for version_info, gil_disabled, expected in cases:
            with self.subTest(version_info=version_info, gil_disabled=gil_disabled):
                self.assertEqual(
                    expected,
                    wheel_check._expected_python_tag(version_info, gil_disabled),
                )
                wheel_check._check_python_tag(
                    (f"{expected[0]}-{expected[1]}-manylinux_2_28_x86_64",),
                    version_info=version_info,
                    gil_disabled=gil_disabled,
                )

    def test_shared_layout_requires_provider_and_rejects_duplicate_runtime(self):
        for suffix in (".so", ".dylib", ".dll"):
            with self.subTest(suffix=suffix):
                extension = Path("/installed/kepler_formal/_native.so")
                provider = Path("/installed/najaeda") / ("libnaja_nl" + suffix)
                links = {extension: (provider.name,)}
                wheel_check._check_shared_native_layout(
                    (extension,), links, (provider,))
                with self.assertRaisesRegex(RuntimeError, "second Naja runtime"):
                    wheel_check._check_shared_native_layout(
                        (extension, extension.parent / provider.name),
                        links, (provider,))
                with self.assertRaisesRegex(RuntimeError, "not owned"):
                    wheel_check._check_shared_native_layout(
                        (extension,), links, ())
                with self.assertRaisesRegex(RuntimeError, "does not link"):
                    wheel_check._check_shared_native_layout(
                        (extension,), {extension: ()}, (provider,))


if __name__ == "__main__":
    unittest.main()
