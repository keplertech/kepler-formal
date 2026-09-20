# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Provider selection and publication guards, without building or uploading."""

import importlib.util
from pathlib import Path
import re
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "prepare_python_release", ROOT / "ci/prepare_python_release.py")
release = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release)

METADATA = '''# Keep unrelated metadata and commands intact.
[build-system]
requires = ["scikit-build-core>=0.11.3,<0.12", "najaeda==0.7.24.dev0"]

[project]
name = "kepler-formal"
dependencies = ["najaeda==0.7.24.dev0"]
dynamic = ["version"]

[tool.scikit-build.cmake.define]
BUILD_KEPLER_PYTHON = "ON"
KEPLER_USE_PUBLISHED_NAJAEDA = "OFF"

[tool.cibuildwheel]
before-build = "python {project}/ci/shared_naja_wheels.py build-provider --project {project}"
'''


class PythonReleaseTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="kepler_python_release_")
        self.addCleanup(self.temporary.cleanup)
        self.project = Path(self.temporary.name)
        self.metadata = self.project / "pyproject.toml"
        self.metadata.write_text(METADATA, encoding="utf-8")
        self.version_file = self.project / "src/bin/KeplerVersion.h.in"
        self.version_file.parent.mkdir(parents=True)
        self.version_file.write_text(
            'inline constexpr std::string_view KEPLER_VERSION { "0.5.0" };\n',
            encoding="utf-8")

    def _validate(self, **changes):
        request = dict(repository="keplertech/kepler-formal", ref="refs/heads/main",
                       confirmed_version="0.5.0")
        request.update(changes)
        return release.validate_release(self.project, **request)

    def test_version_reader_matches_package_metadata(self):
        metadata = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        block = release._section(metadata, "tool.scikit-build.metadata.version").group("body")
        regex = re.search(r"(?m)^regex\s*=\s*'(.*)'$", block).group(1)
        self.assertEqual(regex, release.VERSION_PATTERN)
        self.assertIn('input = "src/bin/KeplerVersion.h.in"', block)
        self.version_file.write_text('KEPLER_VERSION\n{\n "12.34.56" };', encoding="utf-8")
        self.assertEqual("12.34.56", release.read_version(self.project))

    def test_missing_or_ambiguous_version_is_rejected(self):
        for source in ('KEPLER_VERSION { "0.5" };',
                       'KEPLER_VERSION { "0.5.0" }; KEPLER_VERSION { "0.6.0" };'):
            with self.subTest(source=source):
                self.version_file.write_text(source, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "one package version"):
                    release.read_version(self.project)

    def test_default_preparation_does_not_modify_checkout(self):
        before = self.metadata.stat().st_mtime_ns
        release.prepare_project(self.project)
        self.assertEqual(METADATA, self.metadata.read_text(encoding="utf-8"))
        self.assertEqual(before, self.metadata.stat().st_mtime_ns)

    def test_published_preparation_changes_only_pins_and_cmake_option(self):
        release.prepare_project(self.project, published_najaeda=True)
        expected = METADATA.replace("najaeda==0.7.24.dev0", "najaeda==0.7.24").replace(
            'KEPLER_USE_PUBLISHED_NAJAEDA = "OFF"',
            'KEPLER_USE_PUBLISHED_NAJAEDA = "ON"')
        self.assertEqual(expected, self.metadata.read_text(encoding="utf-8"))
        release.prepare_project(self.project, published_najaeda=True)
        self.assertEqual(expected, self.metadata.read_text(encoding="utf-8"))

    def test_published_preparation_adds_missing_cmake_option(self):
        self.metadata.write_text(METADATA.replace(
            'KEPLER_USE_PUBLISHED_NAJAEDA = "OFF"\n', ""), encoding="utf-8")
        release.prepare_project(self.project, published_najaeda=True)
        self.assertEqual("0.5.0", self._validate())

    def test_mismatched_or_unsupported_provider_pins_do_not_modify_checkout(self):
        for text in (METADATA.replace("najaeda==0.7.24.dev0", "najaeda==0.7.24", 1),
                     METADATA.replace("najaeda==0.7.24.dev0", "najaeda>=0.7.24")):
            with self.subTest(text=text):
                self.metadata.write_text(text, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "supported provider"):
                    release.prepare_project(self.project, published_najaeda=True)
                self.assertEqual(text, self.metadata.read_text(encoding="utf-8"))

    def test_default_mode_rejects_a_prepared_release_checkout(self):
        release.prepare_project(self.project, published_najaeda=True)
        with self.assertRaisesRegex(ValueError, "Development builds require"):
            release.prepare_project(self.project)

    def test_release_requires_matching_published_metadata(self):
        with self.assertRaisesRegex(ValueError, "published NajaEDA"):
            self._validate()
        release.prepare_project(self.project, published_najaeda=True)
        self.assertEqual("0.5.0", self._validate())
        prepared = self.metadata.read_text(encoding="utf-8")
        self.metadata.write_text(prepared.replace(
            'KEPLER_USE_PUBLISHED_NAJAEDA = "ON"',
            'KEPLER_USE_PUBLISHED_NAJAEDA = "OFF"'), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "KEPLER_USE_PUBLISHED_NAJAEDA=ON"):
            self._validate()

    def test_release_rejects_fork_branch_tag_and_wrong_version(self):
        release.prepare_project(self.project, published_najaeda=True)
        for changes, message in (
            ({"repository": "nanocoh/kepler-formal"}, "only allowed"),
            ({"ref": "refs/heads/fix7"}, "only allowed"),
            ({"ref": "refs/tags/v0.5.0"}, "only allowed"),
            ({"confirmed_version": ""}, "confirm version 0.5.0"),
            ({"confirmed_version": "1.0.0"}, "confirm version 0.5.0"),
        ):
            with self.subTest(changes=changes):
                with self.assertRaisesRegex(ValueError, message):
                    self._validate(**changes)


if __name__ == "__main__":
    unittest.main()
