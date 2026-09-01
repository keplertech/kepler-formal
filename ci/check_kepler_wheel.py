#!/usr/bin/env python3

# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Validate an installed, repaired Kepler Formal wheel.

cibuildwheel runs this file from outside the source checkout after installing
the repaired wheel.  Keep the check dependency-free so it exercises the same
environment an end user gets from ``pip install kepler-formal``.
"""

from __future__ import annotations

import importlib.machinery
import importlib.metadata
import importlib.util
import os
import platform
import subprocess
import tempfile
from pathlib import Path


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _wheel_tags(distribution: importlib.metadata.Distribution) -> tuple[str, ...]:
    metadata = distribution.read_text("WHEEL")
    _require(metadata is not None, "installed distribution has no WHEEL metadata")
    return tuple(
        line.removeprefix("Tag:").strip()
        for line in metadata.splitlines()
        if line.startswith("Tag:")
    )


def _native_files(distribution: importlib.metadata.Distribution) -> tuple[Path, ...]:
    files = distribution.files
    _require(files is not None, "installed distribution has no file manifest")
    paths = tuple(
        Path(distribution.locate_file(item)).resolve()
        for item in files
        if ".so" in Path(str(item)).name or str(item).endswith(".dylib")
    )
    _require(paths, "wheel contains no native libraries")
    return paths


def _check_linkage(native_files: tuple[Path, ...]) -> None:
    system = platform.system()
    if system == "Darwin":
        command = ("otool", "-L")
    elif system == "Linux":
        command = ("ldd",)
    else:
        raise RuntimeError(f"unsupported wheel test platform: {system}")

    for native_file in native_files:
        completed = subprocess.run(
            (*command, str(native_file)),
            check=True,
            capture_output=True,
            text=True,
        )
        linkage = completed.stdout + completed.stderr
        _require("not found" not in linkage, f"unresolved dependency for {native_file}")
        _require("libpython" not in linkage.lower(), f"{native_file} links libpython")
        if system == "Darwin":
            install_name = subprocess.run(
                ("otool", "-D", str(native_file)),
                check=True,
                capture_output=True,
                text=True,
            ).stdout.splitlines()[1:]
            for line in linkage.splitlines()[1:]:
                dependency = line.strip().split(" (compatibility version", 1)[0]
                if dependency in install_name:
                    continue
                if dependency.startswith("/"):
                    _require(
                        dependency.startswith(("/System/Library/", "/usr/lib/")),
                        f"{native_file} retains an absolute dependency: {dependency}",
                    )


def _check_platform_tag(tags: tuple[str, ...]) -> None:
    system = platform.system()
    machine = platform.machine().lower()
    if system == "Darwin":
        _require(machine == "arm64", f"expected macOS arm64, got {machine}")
        expected = "macosx_11_0_arm64"
    elif system == "Linux":
        _require(machine in {"x86_64", "amd64"}, f"expected Linux x86_64, got {machine}")
        expected = "manylinux_2_28_x86_64"
    else:
        raise RuntimeError(f"unsupported wheel test platform: {system}")
    expected = os.environ.get("KEPLER_WHEEL_EXPECTED_PLATFORM", expected)
    _require(
        any(expected in tag for tag in tags),
        f"wheel tags {tags!r} do not contain the required platform {expected}",
    )


def _check_installed_api() -> None:
    import kepler_formal
    from kepler_formal import VerificationOptions, VerificationStatus, run_cli, verify

    distribution = importlib.metadata.distribution("kepler-formal")
    package_root = Path(kepler_formal.__file__).resolve().parent
    source_root = Path(__file__).resolve().parents[1] / "src" / "python" / "kepler_formal"
    _require(package_root != source_root.resolve(), "import resolved to the source checkout")
    _require(
        kepler_formal.__version__ == distribution.version,
        "module and distribution versions differ",
    )
    _require(bool(kepler_formal.git_hash()), "native module has no git revision")

    native_spec = importlib.util.find_spec("kepler_formal._native")
    _require(native_spec is not None and native_spec.origin is not None, "_native is missing")
    _require(
        any(
            native_spec.origin.endswith(suffix)
            for suffix in importlib.machinery.EXTENSION_SUFFIXES
        ),
        f"_native is not a CPython extension: {native_spec.origin}",
    )

    files = distribution.files
    _require(files is not None, "installed distribution has no file manifest")
    installed_names = tuple(str(item).lower() for item in files)
    for forbidden in ("naja_python", "naja_snl_pyloader"):
        _require(
            not any(forbidden in name for name in installed_names),
            f"wheel unexpectedly bundles {forbidden}",
        )

    tags = _wheel_tags(distribution)
    _check_platform_tag(tags)
    native_files = _native_files(distribution)
    _check_linkage(native_files)

    help_result = run_cli(("--help",))
    _require(help_result.exit_code == 0, "native --help failed")
    _require(help_result.status is VerificationStatus.NO_RESULT, "unexpected --help status")

    with tempfile.TemporaryDirectory(prefix="kepler_formal_wheel_") as temporary:
        root = Path(temporary)
        reference = root / "reference.v"
        equivalent = root / "equivalent.v"
        different = root / "different.v"
        reference.write_text(
            "module top(input a, output y); assign y = a; endmodule\n",
            encoding="utf-8",
        )
        equivalent.write_text(
            "module top(input a, output y); wire n; assign n = a; assign y = n; endmodule\n",
            encoding="utf-8",
        )
        different.write_text(
            "module top(input a, output y); assign y = 1'b0; endmodule\n",
            encoding="utf-8",
        )

        expected = (
            (equivalent, VerificationStatus.EQUIVALENT),
            (different, VerificationStatus.DIFFERENT),
            (equivalent, VerificationStatus.EQUIVALENT),
        )
        for index, (candidate, status) in enumerate(expected):
            result = verify(
                reference,
                candidate,
                options=VerificationOptions(log_file=root / f"run-{index}.log"),
            )
            _require(
                result.status is status,
                f"expected {status.value}, got {result.status.value}",
            )
            _require(
                Path(result.log_file or "").is_file(),
                "verification log was not created",
            )

    print(
        "validated installed kepler-formal "
        f"{distribution.version} ({', '.join(tags)}) with {len(native_files)} native files"
    )


if __name__ == "__main__":
    _check_installed_api()
