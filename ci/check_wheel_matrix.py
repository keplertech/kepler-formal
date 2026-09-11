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

import yaml


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    workflow = yaml.safe_load(
        (root / ".github/workflows/python-wheels.yml").read_text(encoding="utf-8")
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
    artifacts: set[str] = set()
    for row in workflow["jobs"]["build"]["strategy"]["matrix"]["include"]:
        completed = subprocess.run(
            [sys.executable, "-m", "cibuildwheel", "--print-build-identifiers",
             "--platform", platforms[row["platform"]], "--archs", row["arch"]],
            cwd=root, env=environment, check=True, capture_output=True, text=True,
        )
        identifiers = set(completed.stdout.split())
        if not identifiers or actual.intersection(identifiers):
            raise RuntimeError(f"Empty or duplicate wheel selection: {row}")
        actual.update(identifiers)
        artifacts.add(f"kepler-formal-{row['platform']}-{row['arch']}")
    if actual != expected:
        raise RuntimeError(
            f"Wheel matrix differs from NajaEDA. Missing: {sorted(expected - actual)}; "
            f"extra: {sorted(actual - expected)}"
        )
    publisher = workflow["jobs"]["publish"]
    downloaded = {
        step["with"]["name"] for step in publisher["steps"]
        if step.get("uses", "").startswith("actions/download-artifact@")
    }
    if downloaded != artifacts or publisher["needs"] != "build":
        raise RuntimeError("Publisher must wait for and download the complete wheel matrix")
    print(f"Wheel matrix matches NajaEDA: {len(actual)} identifiers, {len(artifacts)} platforms.")
    for identifier in sorted(actual):
        print(f"  {identifier}")


if __name__ == "__main__":
    main()
