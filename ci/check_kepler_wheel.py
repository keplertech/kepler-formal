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
import re
import struct
import subprocess
import sys
import sysconfig
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


def _is_native_filename(name: str) -> bool:
    lower_name = name.lower()
    return (
        lower_name.endswith((".so", ".dylib", ".pyd", ".dll"))
        or ".so." in lower_name
    )


def _native_files(distribution: importlib.metadata.Distribution) -> tuple[Path, ...]:
    files = distribution.files
    _require(files is not None, "installed distribution has no file manifest")
    paths = tuple(
        Path(distribution.locate_file(item)).resolve()
        for item in files
        if _is_native_filename(Path(str(item)).name)
    )
    _require(paths, "wheel contains no native libraries")
    return paths


def _read_pe_imports(native_file: Path) -> tuple[str, ...]:
    """Return the direct DLL imports from a PE image without third-party tools."""

    data = native_file.read_bytes()

    def unpack_from(fmt: str, offset: int, description: str) -> tuple[int, ...]:
        size = struct.calcsize(fmt)
        _require(
            0 <= offset and offset + size <= len(data),
            f"{native_file} has a truncated {description}",
        )
        return struct.unpack_from(fmt, data, offset)

    _require(len(data) >= 64 and data[:2] == b"MZ", f"{native_file} is not a PE image")
    (pe_offset,) = unpack_from("<I", 0x3C, "DOS header")
    _require(
        pe_offset + 24 <= len(data) and data[pe_offset : pe_offset + 4] == b"PE\0\0",
        f"{native_file} has no PE signature",
    )

    (
        machine,
        section_count,
        _timestamp,
        _symbol_table,
        _symbol_count,
        optional_header_size,
        _characteristics,
    ) = unpack_from("<HHIIIHH", pe_offset + 4, "COFF header")
    _require(
        machine == 0x8664,
        f"{native_file} is not an x86_64 PE image (machine 0x{machine:04x})",
    )
    optional_offset = pe_offset + 24
    optional_end = optional_offset + optional_header_size
    _require(optional_end <= len(data), f"{native_file} has a truncated optional header")
    (optional_magic,) = unpack_from("<H", optional_offset, "optional header")
    if optional_magic == 0x10B:  # PE32
        directory_count_offset = optional_offset + 92
        directory_offset = optional_offset + 96
    elif optional_magic == 0x20B:  # PE32+
        directory_count_offset = optional_offset + 108
        directory_offset = optional_offset + 112
    else:
        raise RuntimeError(
            f"{native_file} has unsupported PE optional-header magic "
            f"0x{optional_magic:04x}"
        )

    _require(
        directory_count_offset + 4 <= optional_end,
        f"{native_file} has no PE data-directory count",
    )
    (directory_count,) = unpack_from(
        "<I", directory_count_offset, "data-directory count"
    )
    if directory_count <= 1:
        return ()
    import_directory_offset = directory_offset + 8
    _require(
        import_directory_offset + 8 <= optional_end,
        f"{native_file} has no complete PE import directory",
    )
    import_rva, import_size = unpack_from(
        "<II", import_directory_offset, "import directory"
    )
    if import_rva == 0 or import_size == 0:
        return ()

    sections: list[tuple[int, int, int]] = []
    section_offset = optional_end
    for index in range(section_count):
        current = section_offset + index * 40
        _require(
            current + 40 <= len(data),
            f"{native_file} has a truncated section table",
        )
        virtual_size, virtual_address, raw_size, raw_offset = unpack_from(
            "<IIII", current + 8, "section header"
        )
        sections.append((virtual_address, max(virtual_size, raw_size), raw_offset))

    def rva_to_offset(rva: int) -> int:
        for virtual_address, mapped_size, raw_offset in sections:
            if virtual_address <= rva < virtual_address + mapped_size:
                result = raw_offset + rva - virtual_address
                _require(result < len(data), f"{native_file} has an invalid PE RVA")
                return result
        raise RuntimeError(f"{native_file} has an unmapped PE RVA 0x{rva:x}")

    descriptor_offset = rva_to_offset(import_rva)
    descriptor_end = min(len(data), descriptor_offset + import_size)
    imports: list[str] = []
    while descriptor_offset + 20 <= descriptor_end:
        descriptor = unpack_from("<IIIII", descriptor_offset, "import descriptor")
        if not any(descriptor):
            break
        name_rva = descriptor[3]
        _require(name_rva != 0, f"{native_file} has an import without a DLL name")
        name_offset = rva_to_offset(name_rva)
        name_end = data.find(b"\0", name_offset)
        _require(name_end >= 0, f"{native_file} has an unterminated import name")
        try:
            imports.append(data[name_offset:name_end].decode("ascii"))
        except UnicodeDecodeError as error:
            raise RuntimeError(f"{native_file} has a non-ASCII import name") from error
        descriptor_offset += 20
    return tuple(imports)


def _dependency_basename(dependency: str) -> str:
    return dependency.replace("\\", "/").rsplit("/", 1)[-1]


def _expected_windows_python_dll(
    version_info: tuple[int, int],
    gil_disabled: bool,
) -> str:
    major, minor = version_info
    return f"python{major}{minor}{'t' if gil_disabled else ''}.dll"


def _check_windows_python_imports(
    native_file: Path,
    dependencies: tuple[str, ...],
    *,
    version_info: tuple[int, int] | None = None,
    gil_disabled: bool | None = None,
) -> None:
    if version_info is None:
        version_info = (sys.version_info.major, sys.version_info.minor)
    if gil_disabled is None:
        gil_disabled = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))
    expected_python = _expected_windows_python_dll(
        version_info, gil_disabled
    ).casefold()
    python_dependencies = tuple(
        dependency
        for dependency in dependencies
        if _dependency_basename(dependency).casefold().startswith(
            ("python", "libpython")
        )
    )
    for dependency in python_dependencies:
        _require(
            _dependency_basename(dependency).casefold() == expected_python,
            f"{native_file} imports {dependency}; expected {expected_python}",
        )
    if native_file.suffix.lower() == ".pyd":
        _require(
            bool(python_dependencies),
            f"{native_file} does not import a CPython runtime DLL",
        )


def _unix_dependency_names(linkage: str, system: str) -> tuple[str, ...]:
    dependencies: list[str] = []
    lines = linkage.splitlines()[1:] if system == "Darwin" else linkage.splitlines()
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if system == "Darwin":
            dependency = stripped.split(" (compatibility version", 1)[0]
        elif "=>" in stripped:
            dependency = stripped.split("=>", 1)[0].strip()
        else:
            dependency = stripped.split(maxsplit=1)[0]
        dependencies.append(_dependency_basename(dependency))
    return tuple(dependencies)


def _check_linkage(
    native_files: tuple[Path, ...],
) -> Mapping[Path, tuple[str, ...]]:
    system = platform.system()
    if system == "Darwin":
        command = ("otool", "-L")
    elif system == "Linux":
        command = ("ldd",)
    elif system == "Windows":
        linkages: dict[Path, tuple[str, ...]] = {}
        for native_file in native_files:
            dependencies = _read_pe_imports(native_file)
            _check_windows_python_imports(native_file, dependencies)
            linkages[native_file] = dependencies
        return linkages
    else:
        raise RuntimeError(f"unsupported wheel test platform: {system}")

    linkages = {}
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
        linkages[native_file] = _unix_dependency_names(linkage, system)
    return linkages


def _check_shared_native_layout(
    native_files: tuple[Path, ...],
    linkages: Mapping[Path, tuple[str, ...]],
    provider_files: tuple[Path, ...],
) -> None:
    """Naja DSOs must be supplied only by the separate NajaEDA distribution."""
    for path in native_files:
        name = _canonical_library_name(path.name)
        _require(
            not name.startswith(("naja_", "kepler_najaeda_")) and name != "naja",
            f"Kepler wheel bundles a second Naja runtime: {path}",
        )
    provider_names = {path.name.casefold() for path in provider_files}
    naja_links = {
        _dependency_basename(dependency)
        for dependencies in linkages.values()
        for dependency in dependencies
        if _canonical_library_name(dependency).startswith("naja_")
    }
    _require(naja_links, "Kepler extension does not link the shared Naja provider")
    for dependency in naja_links:
        _require(
            dependency.casefold() in provider_names,
            f"Naja dependency is not owned by the installed provider: {dependency}",
        )


def _is_shared_library_filename(name: str) -> bool:
    lower_name = _dependency_basename(name).lower()
    return (
        lower_name.endswith((".so", ".dylib", ".dll"))
        or ".so." in lower_name
    )


def _canonical_library_name(name: str) -> str:
    basename = _dependency_basename(name).lower()
    if ".so" in basename:
        basename = basename.split(".so", 1)[0]
    elif basename.endswith(".dylib"):
        basename = basename[: -len(".dylib")]
        basename = re.sub(r"\.\d+(?:\.\d+)*$", "", basename)
    elif basename.endswith((".dll", ".pyd")):
        basename = basename.rsplit(".", 1)[0]
    if basename.startswith("lib"):
        basename = basename[3:]
    return basename


def _expected_platform_tag(system: str, machine: str) -> str:
    machine = machine.lower()
    if system == "Darwin":
        _require(machine in {"arm64", "aarch64"}, f"expected macOS arm64, got {machine}")
        return "macosx_11_0_arm64"
    if system == "Linux":
        if machine in {"x86_64", "amd64"}:
            return "manylinux_2_28_x86_64"
        if machine in {"aarch64", "arm64"}:
            return "manylinux_2_28_aarch64"
        raise RuntimeError(f"unsupported Linux wheel architecture: {machine}")
    if system == "Windows":
        _require(
            machine in {"amd64", "x86_64"},
            f"expected Windows x86_64, got {machine}",
        )
        return "win_amd64"
    raise RuntimeError(f"unsupported wheel test platform: {system}")


def _check_platform_tag(
    tags: tuple[str, ...],
    *,
    system: str | None = None,
    machine: str | None = None,
) -> None:
    expected = _expected_platform_tag(
        system or platform.system(), machine or platform.machine()
    )
    expected = os.environ.get("KEPLER_WHEEL_EXPECTED_PLATFORM", expected)
    _require(
        any(expected in tag for tag in tags),
        f"wheel tags {tags!r} do not contain the required platform {expected}",
    )


def _expected_python_tag(
    version_info: tuple[int, int],
    gil_disabled: bool,
) -> tuple[str, str]:
    major, minor = version_info
    _require(
        major == 3 and 10 <= minor <= 15,
        f"unsupported CPython wheel interpreter: {major}.{minor}",
    )
    interpreter = f"cp{major}{minor}"
    abi = f"{interpreter}{'t' if gil_disabled else ''}"
    return interpreter, abi


def _check_python_tag(
    tags: tuple[str, ...],
    *,
    version_info: tuple[int, int] | None = None,
    gil_disabled: bool | None = None,
) -> None:
    if version_info is None:
        version_info = (sys.version_info.major, sys.version_info.minor)
    if gil_disabled is None:
        gil_disabled = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))
    interpreter, abi = _expected_python_tag(version_info, gil_disabled)
    _require(
        any(tag.startswith(f"{interpreter}-{abi}-") for tag in tags),
        f"wheel tags {tags!r} do not contain required tag {interpreter}-{abi}",
    )


def _check_installed_api() -> None:
    import kepler_formal
    import najaeda
    from kepler_formal import (
        Solver, VerificationOptions, VerificationStatus, run_cli, verify,
        verify_designs,
    )

    distribution = importlib.metadata.distribution("kepler-formal")
    provider = importlib.metadata.distribution("najaeda")
    _require(kepler_formal.__version__ == distribution.version,
             "module and distribution versions differ")
    _require(bool(kepler_formal.git_hash()), "native module has no git revision")
    _require(importlib.import_module("kepler_formal.najaeda") is najaeda,
             "nested import is not an alias of the original NajaEDA package")
    native_spec = importlib.util.find_spec("kepler_formal._native")
    _require(native_spec is not None and native_spec.origin is not None,
             "_native is missing")
    _require(any(native_spec.origin.endswith(suffix)
                 for suffix in importlib.machinery.EXTENSION_SUFFIXES),
             "_native is not a CPython extension")
    tags = _wheel_tags(distribution)
    _check_python_tag(tags)
    _check_platform_tag(tags)
    native_files = _native_files(distribution)
    _check_shared_native_layout(
        native_files, _check_linkage(native_files), _native_files(provider))

    # First exercise the existing owning/file API without an editor universe.
    _require(najaeda.naja.NLUniverse.get() is None, "unexpected active universe")
    _require(run_cli(("--help",)).exit_code == 0, "native --help failed")
    with tempfile.TemporaryDirectory(prefix="kepler_formal_wheel_") as temporary:
        root = Path(temporary)
        reference = root / "reference.v"
        equivalent = root / "equivalent.v"
        different = root / "different.v"
        reference.write_text(
            "module top(input a, output y); assign y = a; endmodule\n",
            encoding="utf-8")
        equivalent.write_text(
            "module top(input a, output y); wire n; assign n = a; "
            "assign y = n; endmodule\n", encoding="utf-8")
        different.write_text(
            "module top(input a, output y); assign y = 1'b0; endmodule\n",
            encoding="utf-8")
        for index, (candidate, expected) in enumerate((
            (equivalent, VerificationStatus.EQUIVALENT),
            (different, VerificationStatus.DIFFERENT),
            (equivalent, VerificationStatus.EQUIVALENT),
        )):
            result = verify(reference, candidate, options=VerificationOptions(
                log_file=root / f"file-{index}.log"))
            _require(result.status is expected,
                     f"expected {expected.value}, got {result.status.value}")
            _require(najaeda.naja.NLUniverse.get() is None,
                     "file API leaked its owned universe")

        # Construct original Naja objects and pass their wrappers directly.
        naja = najaeda.naja
        universe = naja.NLUniverse.create()
        db = naja.NLDB.create(universe)
        library = naja.NLLibrary.create(db, "shared")
        def make_design(name: str):
            design = naja.SNLDesign.create(library, name)
            a = naja.SNLScalarTerm.create(design, naja.SNLTerm.Direction.Input, "a")
            y = naja.SNLScalarTerm.create(design, naja.SNLTerm.Direction.Output, "y")
            net = naja.SNLScalarNet.create(design, "a")
            a.setNet(net)
            y.setNet(net)
            return design
        try:
            first, second = make_design("first"), make_design("second")
            universe.setTopDesign(first)
            for solver in (Solver.KISSAT, Solver.CADICAL, Solver.GLUCOSE):
                result = verify_designs(first, second, options=VerificationOptions(
                    solver=solver, log_file=root / f"direct-{solver.value}.log"))
                _require(result.status is VerificationStatus.EQUIVALENT,
                         f"direct {solver.value} returned {result.status.value}")
                _require(naja.NLUniverse.get() is universe,
                         "borrowed call replaced the universe")
                _require(universe.getTopDesign() is first,
                         "borrowed call changed the current top")
                _require(library.getSNLDesign("second") is second,
                         "borrowed call replaced/destroyed its input")
        finally:
            universe.destroy()

    print(f"validated shared-runtime kepler-formal {distribution.version} "
          f"with najaeda {provider.version} ({', '.join(tags)})")


if __name__ == "__main__":
    _check_installed_api()
