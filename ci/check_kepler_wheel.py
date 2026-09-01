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
import sys
import tempfile
from pathlib import Path
from typing import Mapping


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


def _check_linkage(native_files: tuple[Path, ...]) -> Mapping[Path, str]:
    system = platform.system()
    if system == "Darwin":
        command = ("otool", "-L")
    elif system == "Linux":
        command = ("ldd",)
    else:
        raise RuntimeError(f"unsupported wheel test platform: {system}")

    linkages: dict[Path, str] = {}
    for native_file in native_files:
        completed = subprocess.run(
            (*command, str(native_file)),
            check=True,
            capture_output=True,
            text=True,
        )
        linkage = completed.stdout + completed.stderr
        linkages[native_file] = linkage
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
    return linkages


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def _check_nested_package_files(
    distribution: importlib.metadata.Distribution,
    package_root: Path,
) -> None:
    files = distribution.files
    _require(files is not None, "installed distribution has no file manifest")
    installed = {Path(str(item)) for item in files}
    nested_prefix = Path("kepler_formal") / "najaeda"

    _require(
        nested_prefix / "__init__.py" in installed,
        "wheel is missing the nested kepler_formal.najaeda package",
    )
    _require(
        not any(item.parts and item.parts[0] == "najaeda" for item in installed),
        "wheel installs an unwanted top-level najaeda package",
    )
    _require(
        not any(
            any(
                part.startswith("najaeda-") and part.endswith(".dist-info")
                for part in item.parts
            )
            for item in installed
        ),
        "wheel contains duplicate najaeda distribution metadata",
    )

    vendored_root = (
        Path(__file__).resolve().parents[1]
        / "thirdparty"
        / "naja"
        / "src"
        / "najaeda"
        / "najaeda"
    )
    _require(vendored_root.is_dir(), "vendored NajaEDA Python sources are missing")
    expected_sources = {
        nested_prefix / source.relative_to(vendored_root)
        for source in vendored_root.rglob("*.py")
        if "docs" not in source.relative_to(vendored_root).parts
    }
    missing_sources = sorted(expected_sources - installed)
    _require(
        not missing_sources,
        "wheel is missing vendored NajaEDA modules: "
        + ", ".join(str(item) for item in missing_sources),
    )

    nested_root = package_root / "najaeda"
    pure_modules = (
        "kepler_formal.najaeda.netlist",
        "kepler_formal.najaeda.stats",
        "kepler_formal.najaeda.instance_visitor",
        "kepler_formal.najaeda.native.stats",
        "kepler_formal.najaeda.primitives.utils",
        "kepler_formal.najaeda.remote.serialization",
    )
    for name in pure_modules:
        module = importlib.import_module(name)
        origin = Path(module.__file__).resolve()
        _require(
            _is_relative_to(origin, nested_root),
            f"{name} resolved outside the nested package: {origin}",
        )

    _require("najaeda" not in sys.modules, "nested imports leaked top-level najaeda")
    _require("naja" not in sys.modules, "nested imports leaked top-level naja")


def _check_isolated_native_layout(
    package_root: Path,
    native_spec: importlib.machinery.ModuleSpec,
    native_files: tuple[Path, ...],
    linkages: Mapping[Path, str],
) -> None:
    nested_spec = importlib.util.find_spec("kepler_formal.najaeda.naja")
    _require(
        nested_spec is not None and nested_spec.origin is not None,
        "nested native Naja module is missing",
    )
    _require(
        any(
            nested_spec.origin.endswith(suffix)
            for suffix in importlib.machinery.EXTENSION_SUFFIXES
        ),
        f"nested naja is not a CPython extension: {nested_spec.origin}",
    )

    nested_root = package_root / "najaeda"
    nested_extension = Path(nested_spec.origin).resolve()
    _require(
        _is_relative_to(nested_extension, nested_root),
        f"nested naja extension is outside kepler_formal.najaeda: {nested_extension}",
    )

    isolated_dsos = tuple(
        path
        for path in native_files
        if path.name.lower().startswith("libkepler_najaeda_")
    )
    _require(
        len(isolated_dsos) >= 2,
        "wheel does not contain the isolated libkepler_najaeda_* DSO set",
    )
    isolated_names = tuple(path.name.lower() for path in isolated_dsos)
    _require(
        len(isolated_names) == len(set(isolated_names)),
        f"isolated NajaEDA DSO names are not unique: {isolated_names!r}",
    )
    for dso in isolated_dsos:
        _require(
            _is_relative_to(dso, nested_root),
            f"isolated NajaEDA DSO is outside the nested package: {dso}",
        )
    unisolated_nested_dsos = tuple(
        path.name
        for path in native_files
        if _is_relative_to(path, nested_root)
        and path.name.lower().startswith("libnaja_")
    )
    _require(
        not unisolated_nested_dsos,
        f"nested package contains unisolated Naja DSOs: {unisolated_nested_dsos!r}",
    )

    nested_linkage = "\n".join(
        linkages[path] for path in (nested_extension, *isolated_dsos)
    ).lower()
    _require(
        "libkepler_najaeda_" in nested_linkage,
        "nested native module does not link the isolated NajaEDA DSOs",
    )
    for unisolated in (
        "libnaja_nl",
        "libnaja_dnl",
        "libnaja_bne",
        "libnaja_opt",
        "libnaja_metrics",
        "libnaja_python",
    ):
        _require(
            unisolated not in nested_linkage,
            f"nested NajaEDA linkage collides with unisolated {unisolated}",
        )

    kepler_extension = Path(native_spec.origin).resolve()
    kepler_linkage = linkages[kepler_extension].lower()
    _require(
        "libnaja_nl" in kepler_linkage,
        "Kepler native module is not linked to its own Naja core",
    )
    _require(
        "libkepler_najaeda_" not in kepler_linkage,
        "Kepler native module unexpectedly links the nested NajaEDA copy",
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

    _require(
        importlib.util.find_spec("najaeda") is None,
        "a top-level najaeda package is importable",
    )
    _require(
        importlib.util.find_spec("naja") is None,
        "a top-level naja extension is importable",
    )
    try:
        importlib.metadata.distribution("najaeda")
    except importlib.metadata.PackageNotFoundError:
        pass
    else:
        raise RuntimeError("a separate top-level najaeda distribution is installed")

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
    for item in files:
        basename = Path(str(item)).name.lower()
        forbidden = basename.startswith(("libnaja_python", "libnaja_snl_pyloader"))
        _require(
            not forbidden,
            f"wheel unexpectedly bundles unisolated {basename}",
        )

    _check_nested_package_files(distribution, package_root)

    tags = _wheel_tags(distribution)
    _check_platform_tag(tags)
    native_files = _native_files(distribution)
    linkages = _check_linkage(native_files)
    _check_isolated_native_layout(
        package_root,
        native_spec,
        native_files,
        linkages,
    )

    nested_najaeda = importlib.import_module("kepler_formal.najaeda")
    nested_naja = nested_najaeda.naja
    _require(
        nested_naja.NLUniverse.get() is None,
        "nested NajaEDA universe was already active",
    )
    nested_universe = nested_naja.NLUniverse.create()
    nested_db = nested_naja.NLDB.create(nested_universe)
    nested_universe.setTopDB(nested_db)
    nested_library = nested_naja.NLLibrary.create(
        nested_db,
        "kepler_wheel_isolation",
    )

    def require_nested_universe() -> None:
        _require(
            nested_naja.NLUniverse.get() is nested_universe,
            "Kepler invocation destroyed or replaced the nested NajaEDA universe",
        )
        _require(
            nested_universe.getTopDB() is nested_db,
            "Kepler invocation changed the nested NajaEDA top database",
        )
        _require(
            nested_db.getLibrary("kepler_wheel_isolation") is nested_library,
            "Kepler invocation changed the nested NajaEDA database",
        )

    try:
        help_result = run_cli(("--help",))
        _require(help_result.exit_code == 0, "native --help failed")
        _require(
            help_result.status is VerificationStatus.NO_RESULT,
            "unexpected --help status",
        )
        require_nested_universe()

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
                "module top(input a, output y); wire n; "
                "assign n = a; assign y = n; endmodule\n",
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
                require_nested_universe()
    finally:
        surviving_universe = nested_naja.NLUniverse.get()
        if surviving_universe is not None:
            surviving_universe.destroy()

    print(
        "validated installed kepler-formal "
        f"{distribution.version} ({', '.join(tags)}) with {len(native_files)} native files"
    )


if __name__ == "__main__":
    _check_installed_api()
