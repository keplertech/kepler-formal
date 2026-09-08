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
    linkages: Mapping[Path, tuple[str, ...]],
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
        if _is_shared_library_filename(path.name)
        and path.name.lower().startswith("libkepler_najaeda_")
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
        and _is_shared_library_filename(path.name)
        and _canonical_library_name(path.name) in _UNISOLATED_NAJA_LIBRARIES
    )
    _require(
        not unisolated_nested_dsos,
        f"nested package contains unisolated Naja DSOs: {unisolated_nested_dsos!r}",
    )

    nested_dependencies = {
        dependency
        for path in (nested_extension, *isolated_dsos)
        for dependency in linkages[path]
    }
    _require(
        any(
            _dependency_basename(dependency).lower().startswith(
                "libkepler_najaeda_"
            )
            for dependency in nested_dependencies
        ),
        "nested native module does not link the isolated NajaEDA DSOs",
    )
    nested_unisolated = sorted(
        {
            _dependency_basename(dependency)
            for dependency in nested_dependencies
            if _canonical_library_name(dependency) in _UNISOLATED_NAJA_LIBRARIES
        }
    )
    _require(
        not nested_unisolated,
        "nested NajaEDA linkage collides with unisolated libraries: "
        + ", ".join(nested_unisolated),
    )

    kepler_extension = Path(native_spec.origin).resolve()
    kepler_dependencies = linkages[kepler_extension]
    kepler_library_names = {
        _canonical_library_name(dependency) for dependency in kepler_dependencies
    }
    _require(
        "naja_nl" in kepler_library_names,
        "Kepler native module is not linked to its own Naja core",
    )
    _require(
        not any(
            _dependency_basename(dependency).lower().startswith(
                "libkepler_najaeda_"
            )
            for dependency in kepler_dependencies
        ),
        "Kepler native module unexpectedly links the nested NajaEDA copy",
    )


_UNISOLATED_NAJA_LIBRARIES = {
    "naja_nl",
    "naja_dnl",
    "naja_bne",
    "naja_opt",
    "naja_metrics",
    "naja_python",
    "naja_snl_pyloader",
}


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
    from kepler_formal import (
        Solver,
        VerificationOptions,
        VerificationStatus,
        run_cli,
        verify,
    )

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
        forbidden = (
            _is_shared_library_filename(basename)
            and _canonical_library_name(basename)
            in {"naja_python", "naja_snl_pyloader"}
        )
        _require(
            not forbidden,
            f"wheel unexpectedly bundles unisolated {basename}",
        )

    _check_nested_package_files(distribution, package_root)

    tags = _wheel_tags(distribution)
    _check_python_tag(tags)
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

            for solver in (Solver.KISSAT, Solver.CADICAL, Solver.GLUCOSE):
                result = verify(
                    reference,
                    equivalent,
                    options=VerificationOptions(
                        solver=solver,
                        log_file=root / f"solver-{solver.value}.log",
                    ),
                )
                _require(
                    result.status is VerificationStatus.EQUIVALENT,
                    f"{solver.value} wheel smoke test returned {result.status.value}",
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
