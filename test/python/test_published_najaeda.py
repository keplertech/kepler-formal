# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch


_SOURCE = Path(__file__).resolve().parents[2] / "src/python"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = _load("published_najaeda_builder", _SOURCE / "published_najaeda.py")
guard = _load("published_najaeda_guard", _SOURCE / "kepler_formal/_published_runtime.py")


class PublishedProviderGuardTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="kepler published provider ")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        package = self.root / "najaeda"
        package.mkdir()
        (package / "__init__.py").touch()
        self.extension = package / "naja.so"
        self.extension.write_bytes(b"provider native module")
        self.library = package / "libnaja_nl.so"
        self.library.write_bytes(b"provider runtime")
        self.manifest = {
            path.relative_to(self.root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (self.extension, self.library)
        }
        self.naja = ModuleType("najaeda.naja")
        self.naja.__file__ = str(self.extension)
        self.naja.getGitHash = lambda: "2263958"
        najaeda = ModuleType("najaeda")
        najaeda.__file__ = str(package / "__init__.py")
        najaeda.naja = self.naja
        self._patch(patch.dict(sys.modules, {"najaeda": najaeda, "najaeda.naja": self.naja}))
        self.version = self._patch(patch.object(guard.importlib.metadata, "version", return_value="0.7.24"))
        self.distribution = self._patch(patch.object(
            guard.importlib.metadata, "distribution",
            return_value=SimpleNamespace(locate_file=lambda path: self.root / path)))
        guard.validate_provider.cache_clear()
        self.addCleanup(guard.validate_provider.cache_clear)

    def _patch(self, patcher):
        result = patcher.start()
        self.addCleanup(patcher.stop)
        return result

    def validate(self, manifest=None):
        guard.validate_provider(json.dumps(self.manifest if manifest is None else manifest))

    def test_accepts_matching_native_files(self):
        self.validate()
        self.assertEqual((self.root / "najaeda", self.extension), builder._provider())

    def test_rejects_replaced_native_library(self):
        self.library.write_bytes(b"different binary with same release number")
        with self.assertRaisesRegex(ImportError, "runtime differs"):
            self.validate()

    def test_rejects_version_and_revision_mismatch(self):
        for version, revision in (("0.7.25", "2263958"), ("0.7.24", "unknown")):
            with self.subTest(version=version, revision=revision):
                self.version.return_value = version
                self.naja.getGitHash = lambda: revision
                with self.assertRaisesRegex(ImportError, "requires NajaEDA"):
                    self.validate()
                with self.assertRaisesRegex(RuntimeError, "requires 0.7.24"):
                    builder._provider()

    def test_rejects_wrong_distribution_and_foreign_extension(self):
        self.distribution.return_value = SimpleNamespace(locate_file=lambda path: self.root / "foreign" / path)
        with self.assertRaisesRegex(ImportError, "installed distribution"):
            self.validate()
        self.distribution.return_value = SimpleNamespace(locate_file=lambda path: self.root / path)
        self.naja.__file__ = str(self.root / "foreign/naja.so")
        with self.assertRaisesRegex(ImportError, "outside its package"):
            self.validate()

    def test_rejects_missing_file_and_manifest_without_extension(self):
        self.library.unlink()
        with self.assertRaisesRegex(ImportError, "Cannot validate"):
            self.validate()
        with self.assertRaisesRegex(ImportError, "does not identify"):
            self.validate({"najaeda/libnaja_nl.so": "0" * 64})

    def test_rejects_empty_manifest_and_path_escapes(self):
        with self.assertRaisesRegex(ImportError, "invalid NajaEDA provider manifest"):
            self.validate({})
        for path in ("../outside", "/outside", "najaeda/../outside", "najaeda\\outside", "other/package"):
            with self.subTest(path=path):
                manifest = {**self.manifest, path: "0" * 64}
                with self.assertRaisesRegex(ImportError, "invalid NajaEDA provider manifest entry"):
                    self.validate(manifest)


class PublishedProviderBuildTest(unittest.TestCase):
    def test_library_discovery_uses_wheel_repaired_names_and_bundled_tbb(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve() / "najaeda"
            root.mkdir()
            siblings = root.parent / "najaeda.libs"
            siblings.mkdir()
            expected = {}
            for name in (*builder.RUNTIMES, "tbb", "tbbmalloc"):
                path = siblings / f"lib{name}-123abc.so.1"
                path.touch()
                expected[name] = path
            self.assertEqual(expected, builder._libraries(root))
            duplicate = root / "libnaja_nl.so"
            duplicate.touch()
            with self.assertRaisesRegex(RuntimeError, "Expected one naja_nl"):
                builder._libraries(root)
            duplicate.unlink()
            expected["naja_nl"].unlink()
            with self.assertRaisesRegex(RuntimeError, "Expected one naja_nl"):
                builder._libraries(root)

    def test_rejects_source_archive_with_wrong_checksum(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "wrong.tar.gz"
            archive.write_bytes(b"not the pinned release archive")
            with self.assertRaisesRegex(RuntimeError, "checksum mismatch"):
                builder._headers(root, archive)
            self.assertFalse((root / "release-headers").exists())

    def test_cmake_values_preserve_spaces_and_reject_control_delimiters(self):
        self.assertEqual("[==[C:/Program Files/provider]==]", builder._cmake_value(r"C:\Program Files\provider"))
        for value in ("path;other", "path\nother", "path\rother", "path]==]other"):
            with self.subTest(value=value), self.assertRaises(RuntimeError):
                builder._cmake_value(value)

    def test_pe_export_reader_distinguishes_function_and_uninitialized_data(self):
        # The cache singleton is an exported data symbol. Misclassifying it as
        # a function creates a bad Windows import library even when linking succeeds.
        data = bytearray(0x510)
        pe, optional, optional_size = 0x80, 0x98, 0xF0
        data[:2] = b"MZ"
        struct.pack_into("<I", data, 0x3C, pe)
        data[pe:pe + 4] = b"PE\0\0"
        struct.pack_into("<HHIIIHH", data, pe + 4, 0x8664, 3, 0, 0, 0, optional_size, 0x2022)
        struct.pack_into("<H", data, optional, 0x20B)
        struct.pack_into("<II", data, optional + 112, 0x1000, 0x100)
        for index, (rva, raw, raw_size, virtual_size, flags) in enumerate((
            (0x1000, 0x200, 0x200, 0x200, 0x40000040),
            (0x2000, 0x400, 0x100, 0x100, 0x60000020),
            (0x3000, 0x500, 0x10, 0x100, 0xC0000040),
        )):
            section = optional + optional_size + 40 * index
            struct.pack_into("<IIII", data, section + 8, virtual_size, rva, raw_size, raw)
            struct.pack_into("<I", data, section + 36, flags)
        struct.pack_into("<IIHHIIIIIII", data, 0x200, 0, 0, 0, 0, 0, 0, 2, 2, 0x1050, 0x1058, 0x1060)
        struct.pack_into("<II", data, 0x250, 0x2000, 0x3080)
        struct.pack_into("<II", data, 0x258, 0x1070, 0x1080)
        struct.pack_into("<HH", data, 0x260, 0, 1)
        data[0x270:0x279] = b"function\0"
        data[0x280:0x28E] = b"cache_pointer\0"
        with tempfile.TemporaryDirectory() as temporary:
            library = Path(temporary) / "provider.dll"
            library.write_bytes(data)
            self.assertEqual([("function", False), ("cache_pointer", True)], builder._pe_exports(library))
            library.write_bytes(data[:0x90])
            with self.assertRaisesRegex(RuntimeError, "Truncated PE"):
                builder._pe_exports(library)


if __name__ == "__main__":
    unittest.main()
