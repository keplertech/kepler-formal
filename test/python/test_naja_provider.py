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
_VERSION = "0.7.27"
_GIT_HASH = "abc1234"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = _load("naja_provider_builder", _SOURCE / "naja_provider.py")
guard = _load("naja_provider_guard", _SOURCE / "kepler_formal/_provider_check.py")


class _FakeProviderCase(unittest.TestCase):
    """An importable fake NajaEDA package with native files and metadata."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="kepler naja provider ")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.package = self.root / "najaeda"
        self.package.mkdir()
        (self.package / "__init__.py").touch()
        self.extension = self.package / "naja.so"
        self.extension.write_bytes(b"provider native module")
        windows = sys.platform == "win32"
        self.libraries = {}
        for name in builder.RUNTIMES:
            library = self.package / (f"{name}.dll" if windows else f"lib{name}.so")
            library.write_bytes(b"provider runtime " + name.encode())
            self.libraries[name] = library
            if windows:
                (self.package / "sdk/lib").mkdir(parents=True, exist_ok=True)
                (self.package / "sdk/lib" / f"{name}.lib").touch()
        self.manifest = {
            path.relative_to(self.root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (self.extension, *self.libraries.values())
        }
        self.naja = ModuleType("najaeda.naja")
        self.naja.__file__ = str(self.extension)
        self.naja.getGitHash = lambda: _GIT_HASH
        self.naja.getVersion = lambda: _VERSION
        najaeda = ModuleType("najaeda")
        najaeda.__file__ = str(self.package / "__init__.py")
        najaeda.naja = self.naja
        self._patch(patch.dict(sys.modules, {"najaeda": najaeda, "najaeda.naja": self.naja}))
        self.distribution = self._patch(patch.object(
            guard.importlib.metadata, "distribution", return_value=self.installed(_VERSION)))
        guard.validate_provider.cache_clear()
        self.addCleanup(guard.validate_provider.cache_clear)

    def _patch(self, patcher):
        result = patcher.start()
        self.addCleanup(patcher.stop)
        return result

    def installed(self, version, prefix=None):
        """Distribution metadata for a pip-installed package under ``prefix``."""
        prefix = self.root if prefix is None else prefix
        return SimpleNamespace(version=version, locate_file=lambda path: prefix / path)

    def without_metadata(self):
        """A plain CMake install, as in the source regression, has no metadata."""
        self.distribution.side_effect = guard.importlib.metadata.PackageNotFoundError("najaeda")

    def add_wheel_headers(self):
        include = self.package / "sdk/include"
        for relative in ("naja/core/NajaVersion.h", "naja/core/NLUniverse.h",
                         "naja/dnl/DNL.h", "spdlog/spdlog.h"):
            (include / relative).parent.mkdir(parents=True, exist_ok=True)
            (include / relative).touch()
        return include


class ProviderGuardTest(_FakeProviderCase):
    def validate(self, manifest=None, version=_VERSION, git_hash=_GIT_HASH):
        guard.validate_provider(
            json.dumps(self.manifest if manifest is None else manifest), version, git_hash)

    def test_accepts_matching_native_files(self):
        self.validate()

    def test_accepts_a_cmake_install_without_distribution_metadata(self):
        self.without_metadata()
        self.validate()
        self.naja.getVersion = lambda: "0.7.25"
        guard.validate_provider.cache_clear()
        with self.assertRaisesRegex(ImportError, "found 0.7.25"):
            self.validate()

    def test_rejects_replaced_native_library(self):
        self.libraries["naja_nl"].write_bytes(b"different binary with same release number")
        with self.assertRaisesRegex(ImportError, "runtime differs"):
            self.validate()

    def test_rejects_version_and_revision_mismatch(self):
        for version, revision in (("0.7.25", _GIT_HASH), (_VERSION, "unknown")):
            with self.subTest(version=version, revision=revision):
                with self.assertRaisesRegex(ImportError, "built against NajaEDA"):
                    self.validate(version=version, git_hash=revision)

    def test_rejects_wrong_distribution_and_foreign_extension(self):
        self.distribution.return_value = self.installed(_VERSION, self.root / "foreign")
        with self.assertRaisesRegex(ImportError, "installed distribution"):
            self.validate()
        self.distribution.return_value = self.installed(_VERSION)
        self.naja.__file__ = str(self.root / "foreign/naja.so")
        with self.assertRaisesRegex(ImportError, "outside its package"):
            self.validate()

    def test_rejects_missing_file_and_manifest_without_extension(self):
        self.libraries["naja_nl"].unlink()
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


class ProviderDescriptionTest(_FakeProviderCase):
    def test_development_provider_requires_wheel_headers(self):
        with self.assertRaisesRegex(RuntimeError, "ships no headers"):
            builder._provider()
        self.add_wheel_headers()
        expected = builder.Provider(self.package, self.extension, _VERSION, _GIT_HASH)
        self.assertEqual(expected, builder._provider())
        self.without_metadata()
        self.assertEqual(expected, builder._provider())

    def test_published_provider_requires_the_pinned_release(self):
        with self.assertRaisesRegex(RuntimeError, f"requires {builder.PUBLISHED_VERSION}"):
            builder._provider(published=True)
        self.distribution.return_value = self.installed(builder.PUBLISHED_VERSION)
        self.naja.getGitHash = lambda: builder.PUBLISHED_GIT_COMMIT
        self.naja.getVersion = lambda: builder.PUBLISHED_VERSION
        # The published release ships no headers; they come from its archive.
        self.assertEqual(builder.PUBLISHED_VERSION, builder._provider(published=True).version)

    def test_native_version_must_match_the_distribution(self):
        self.add_wheel_headers()
        self.naja.getVersion = lambda: "0.0.0"
        with self.assertRaisesRegex(RuntimeError, "reports native version"):
            builder._provider()

    def test_wheel_headers_list_every_naja_header_directory(self):
        include = self.add_wheel_headers()
        self.assertEqual([include, include / "naja/core", include / "naja/dnl"],
                         builder._wheel_headers(self.package))
        (include / "naja/core/NajaVersion.h").unlink()
        with self.assertRaisesRegex(RuntimeError, "no version header"):
            builder._wheel_headers(self.package)

    def test_cmake_config_describes_the_development_provider(self):
        self.add_wheel_headers()
        output = self.root / "consumer build"
        text = builder.cmake_config(output)
        self.assertIn(f"set(NajaEDA_VERSION [==[{_VERSION}]==])", text)
        self.assertIn(f"set(NajaEDA_GIT_COMMIT [==[{_GIT_HASH}]==])", text)
        self.assertIn(f"set(NajaEDA_naja_nl_LIBRARY [==[{self.libraries['naja_nl'].as_posix()}]==])", text)
        self.assertIn((self.package / "sdk/include/naja/core").as_posix(), text)
        header = (output / "include/KeplerNajaProviderBuild.h").read_text(encoding="utf-8")
        self.assertIn(f'#define KEPLER_NAJA_PROVIDER_VERSION "{_VERSION}"', header)
        self.assertIn(f'#define KEPLER_NAJA_PROVIDER_GIT_HASH "{_GIT_HASH}"', header)
        manifest_line = next(line for line in header.splitlines() if "PROVIDER_MANIFEST" in line)
        manifest = json.loads(json.loads(manifest_line.split(" ", 2)[2]))
        self.assertEqual(self.manifest, manifest)


class ProviderBuildTest(unittest.TestCase):
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
                builder._published_headers(root, archive)
            self.assertFalse((root / "release-headers").exists())

    def test_import_library_prefers_the_wheel_import_library(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "sdk/lib").mkdir(parents=True)
            original = root / "sdk/lib/naja_nl.lib"
            original.touch()
            library = root / "naja_nl.dll"
            library.write_bytes(b"MZ")
            self.assertEqual(original, builder._import_library("naja_nl", library, root / "out", root))
            # A repair-renamed DLL needs a regenerated import library instead.
            renamed = root / "naja_nl-1a2b3c.dll"
            renamed.write_bytes(b"MZ")
            with self.assertRaises(RuntimeError):
                builder._import_library("naja_nl", renamed, root / "out", root)

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
