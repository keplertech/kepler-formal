#!/usr/bin/env python3
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Check configured wheel IDs against vendored NajaEDA, without building wheels.

CI-only dependencies: cibuildwheel (the workflow's version) and PyYAML.
"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


def _conditions(expression: str) -> set[str]:
    """Accept only the simple conjunctions used by this release workflow."""
    expression = expression.strip()
    if expression.startswith("${{") and expression.endswith("}}"):
        expression = expression[3:-2].strip()
    return {part.strip() for part in expression.split("&&")}


def check_release_jobs(workflow: dict, build_workflow: dict) -> set[str]:
    """Keep default regression independent of the published-provider release gate."""
    events = workflow.get("on", workflow.get(True, {}))  # PyYAML's YAML 1.1 'on'.
    inputs = events["workflow_dispatch"]["inputs"]
    if ("published_najaeda" in inputs or inputs["publish"].get("type") != "boolean"
            or inputs["publish"].get("default", False) is not False):
        raise RuntimeError("The boolean publish input must control provider selection")
    if "pull_request" not in events or "main" not in events.get("push", {}).get("branches", []):
        raise RuntimeError("Development wheel regression must run on pull requests and main pushes")
    jobs = workflow["jobs"]
    for provider in ("development", "published"):
        caller = jobs[provider]
        if (caller.get("uses") != "./.github/workflows/python-wheel-build.yml"
                or caller.get("with", {}).get("provider") != provider
                or "strategy" in caller
                or caller.get("needs") not in ("generate-parser", ["generate-parser"])):
            raise RuntimeError("Both providers must independently call the shared wheel build")
    if jobs["development"].get("if"):
        raise RuntimeError("Development wheel regression must run by default")
    release_gate = {"github.event_name == 'workflow_dispatch'", "inputs.publish"}
    if _conditions(jobs["published"].get("if", "")) != release_gate:
        raise RuntimeError("Only publish=true dispatches may run published-provider tests")
    if jobs["published"]["with"].get("version") != "${{ inputs.version }}":
        raise RuntimeError("Published builds must validate the requested release version")
    reusable_events = build_workflow.get("on", build_workflow.get(True, {}))
    provider_input = reusable_events["workflow_call"]["inputs"]["provider"]
    if provider_input.get("type") != "string" or provider_input.get("required") is not True:
        raise RuntimeError("The shared build must require an explicit provider")
    build = build_workflow["jobs"]["build"]
    matrix = build["strategy"]["matrix"]
    rows = matrix["include"]
    if (set(matrix) != {"os", "include"}
            or sorted(matrix["os"]) != sorted(row["os"] for row in rows)
            or len(set(matrix["os"])) != len(rows)
            or any("provider" in row for row in rows)):
        raise RuntimeError("Each OS must have exactly one provider-independent platform row")
    if build.get("continue-on-error", False) or build.get("if"):
        raise RuntimeError("The complete published-provider matrix must succeed before publication")
    if build["env"]["KEPLER_USE_PUBLISHED_NAJAEDA"] != (
        "${{ inputs.provider == 'published' && '1' || '0' }}"
    ):
        raise RuntimeError("Each provider must select its own build mode")
    uploaded = [
        step["with"]["name"] for step in build["steps"]
        if step.get("uses", "").startswith("actions/upload-artifact@")
    ]
    if uploaded != [
        "kepler-formal-${{ inputs.provider }}-${{ matrix.platform }}-${{ matrix.arch }}"
    ]:
        raise RuntimeError("Wheel artifacts must distinguish providers and platforms")
    artifacts = {f"kepler-formal-published-{row['platform']}-{row['arch']}" for row in rows}
    publisher = workflow["jobs"]["publish"]
    downloaded = [
        step["with"]["name"] for step in publisher["steps"]
        if step.get("uses", "").startswith("actions/download-artifact@")
    ]
    if (set(downloaded) != artifacts or len(downloaded) != len(rows)
            or publisher["needs"] not in ("published", ["published"])):
        raise RuntimeError("Publisher must await only published tests and download only published wheels")
    # No status override may bypass a failure of the published-provider tests.
    if _conditions(publisher["if"]) != release_gate | {
        "github.repository == 'keplertech/kepler-formal'", "github.ref == 'refs/heads/main'",
    }:
        raise RuntimeError("Publication requires a successful release dispatch on upstream main")
    return artifacts


def main() -> None:
    import yaml

    root = Path(__file__).resolve().parents[1]
    workflow = yaml.safe_load(
        (root / ".github/workflows/python-wheels.yml").read_text(encoding="utf-8")
    )
    build_workflow = yaml.safe_load(
        (root / ".github/workflows/python-wheel-build.yml").read_text(encoding="utf-8")
    )
    naja_workflow = yaml.safe_load(
        (root / "thirdparty/naja/.github/workflows/wheels.yml").read_text(
            encoding="utf-8"
        )
    )
    expected = {
        f"cp{row['python']}-{row['platform_id']}"
        for row in naja_workflow["jobs"]["build_wheels"]["strategy"]["matrix"]["include"]
    }
    platforms = {"manylinux_2_28": "linux", "macOS": "macos", "Windows": "windows"}
    # Compare the checked-in configuration, not a caller's local CIBW overrides.
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith("CIBW_")}
    actual: set[str] = set()
    artifacts = check_release_jobs(workflow, build_workflow)
    for row in build_workflow["jobs"]["build"]["strategy"]["matrix"]["include"]:
        completed = subprocess.run(
            [sys.executable, "-m", "cibuildwheel", "--print-build-identifiers",
             "--platform", platforms[row["platform"]], "--archs", row["arch"]],
            cwd=root, env=environment, check=True, capture_output=True, text=True,
        )
        identifiers = set(completed.stdout.split())
        if not identifiers or actual.intersection(identifiers):
            raise RuntimeError(f"Empty or duplicate wheel selection: {row}")
        actual.update(identifiers)
    if actual != expected:
        raise RuntimeError(
            f"Wheel matrix differs from NajaEDA. Missing: {sorted(expected - actual)}; "
            f"extra: {sorted(actual - expected)}"
        )
    print(f"Wheel matrix matches NajaEDA: {len(actual)} identifiers, {len(artifacts)} platforms.")
    for identifier in sorted(actual):
        print(f"  {identifier}")


if __name__ == "__main__":
    main()
