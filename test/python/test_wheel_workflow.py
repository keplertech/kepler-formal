# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Release dependency and artifact guards without building or publishing wheels."""

from copy import deepcopy
import importlib.util
from pathlib import Path
import unittest


SPEC = importlib.util.spec_from_file_location(
    "check_wheel_matrix", Path(__file__).resolve().parents[2] / "ci/check_wheel_matrix.py"
)
checker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checker)


def workflows():
    gate = "github.event_name == 'workflow_dispatch' && inputs.publish"
    rows = [
        {"os": "ubuntu-24.04", "platform": "manylinux_2_28", "arch": "x86_64"},
        {"os": "windows-latest", "platform": "Windows", "arch": "AMD64"},
    ]
    outer = {
        "on": {
            "pull_request": {}, "push": {"branches": ["main"]},
            "workflow_dispatch": {"inputs": {"publish": {"type": "boolean", "default": False}}},
        },
        "jobs": {
            provider: {
                "uses": "./.github/workflows/python-wheel-build.yml",
                "needs": "generate-parser", "with": {"provider": provider},
            } for provider in ("development", "published")
        },
    }
    outer["jobs"]["published"].update({"if": "${{ " + gate + " }}"})
    outer["jobs"]["published"]["with"]["version"] = "${{ inputs.version }}"
    outer["jobs"]["publish"] = {
        "needs": "published",
        "if": gate + " && github.repository == 'keplertech/kepler-formal' && github.ref == 'refs/heads/main'",
        "steps": [
            {"uses": "actions/download-artifact@v8", "with": {
                "name": f"kepler-formal-published-{row['platform']}-{row['arch']}"
            }} for row in rows
        ],
    }
    reusable = {
        "on": {"workflow_call": {"inputs": {"provider": {"type": "string", "required": True}}}},
        "jobs": {"build": {
            "strategy": {"matrix": {"os": [row["os"] for row in rows], "include": rows}},
            "env": {"KEPLER_USE_PUBLISHED_NAJAEDA": "${{ inputs.provider == 'published' && '1' || '0' }}"},
            "steps": [{"uses": "actions/upload-artifact@v4", "with": {
                "name": "kepler-formal-${{ inputs.provider }}-${{ matrix.platform }}-${{ matrix.arch }}"
            }}],
        }},
    }
    return outer, reusable


class WheelWorkflowTest(unittest.TestCase):
    def setUp(self):
        self.outer, self.reusable = workflows()

    def check(self):
        return checker.check_release_jobs(self.outer, self.reusable)

    def test_independent_release_gate_accepts_both_yaml_event_key_forms(self):
        expected = {
            "kepler-formal-published-manylinux_2_28-x86_64",
            "kepler-formal-published-Windows-AMD64",
        }
        self.assertEqual(expected, self.check())
        for workflow in (self.outer, self.reusable):
            workflow[True] = workflow.pop("on")
        self.assertEqual(expected, self.check())

    def test_development_cannot_gate_publication_directly_or_indirectly(self):
        for job, needs in (("publish", ["development", "published"]),
                           ("published", ["generate-parser", "development"]),
                           ("development", ["generate-parser", "published"])):
            with self.subTest(job=job):
                outer = deepcopy(self.outer)
                outer["jobs"][job]["needs"] = needs
                with self.assertRaises(RuntimeError):
                    checker.check_release_jobs(outer, self.reusable)

    def test_default_regression_requires_pull_requests_main_and_unconditional_job(self):
        for missing in ("pull_request", "main", "unconditional"):
            with self.subTest(missing=missing):
                outer = deepcopy(self.outer)
                if missing == "pull_request":
                    del outer["on"]["pull_request"]
                elif missing == "main":
                    outer["on"]["push"] = {"tags": ["v*"]}
                else:
                    outer["jobs"]["development"]["if"] = "inputs.publish"
                with self.assertRaisesRegex(RuntimeError, "Development wheel regression"):
                    checker.check_release_jobs(outer, self.reusable)

    def test_published_provider_requires_explicit_release_dispatch(self):
        for gate in ("", "inputs.publish", "github.event_name == 'workflow_dispatch'", "always()"):
            with self.subTest(gate=gate):
                self.outer["jobs"]["published"]["if"] = gate
                with self.assertRaisesRegex(RuntimeError, "Only publish=true"):
                    self.check()

    def test_failures_of_published_tests_cannot_be_ignored(self):
        self.outer["jobs"]["publish"]["if"] = "always() && " + self.outer["jobs"]["publish"]["if"]
        with self.assertRaisesRegex(RuntimeError, "successful release dispatch"):
            self.check()
        self.outer, self.reusable = workflows()
        self.reusable["jobs"]["build"]["continue-on-error"] = True
        with self.assertRaisesRegex(RuntimeError, "matrix must succeed"):
            self.check()

    def test_only_complete_published_artifacts_are_uploaded(self):
        downloads = self.outer["jobs"]["publish"]["steps"]
        downloads[0]["with"]["name"] = downloads[0]["with"]["name"].replace("published", "development")
        with self.assertRaisesRegex(RuntimeError, "only published wheels"):
            self.check()
        self.outer, self.reusable = workflows()
        self.outer["jobs"]["publish"]["steps"].pop()
        with self.assertRaisesRegex(RuntimeError, "only published wheels"):
            self.check()
        self.outer, self.reusable = workflows()
        self.reusable["jobs"]["build"]["steps"][0]["with"]["name"] = "kepler-formal-${{ matrix.platform }}"
        with self.assertRaisesRegex(RuntimeError, "distinguish providers"):
            self.check()


if __name__ == "__main__":
    unittest.main()
