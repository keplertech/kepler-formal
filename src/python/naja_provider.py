# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0
"""Describe the installed NajaEDA provider to Kepler's CMake build.

Kepler's Python extension links the native libraries of the NajaEDA package
installed in the build interpreter. The default provider is a wheel built from
``thirdparty/naja``, which ships its headers under ``najaeda/sdk/include``.
``--published`` selects the unmodified NajaEDA 0.7.24 release instead and takes
its headers from the verified source archive. This helper never builds or
edits Naja. PE export/import-library handling and Mach-O fixup are adapted from
the former Apache-2.0 licensed najaeda/sdk.py (The Naja authors).
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import struct
import subprocess
import sys
import tarfile
import tempfile
from typing import NamedTuple
import urllib.request

PUBLISHED_VERSION = "0.7.24"
PUBLISHED_GIT_COMMIT = "2263958"
PUBLISHED_SOURCE_SHA256 = "00145fd267e50031e92637160df3b164e9da5df2f75591ce5bd3f2f786a7abc0"
PUBLISHED_SOURCE_URL = "https://files.pythonhosted.org/packages/source/n/najaeda/najaeda-0.7.24.tar.gz"
RUNTIMES = ("naja_nl", "naja_dnl", "naja_bne", "naja_opt", "naja_metrics", "naja_python")


class Provider(NamedTuple):
    root: Path
    extension: Path
    version: str
    git_commit: str


def _provider_root() -> tuple[Path, Path, str | None]:
    """Return the imported package root, its extension and distribution version."""
    import najaeda
    from najaeda import naja

    root = Path(najaeda.__file__).resolve().parent
    extension = Path(naja.__file__).resolve()
    try:
        distribution = importlib.metadata.distribution("najaeda")
    except importlib.metadata.PackageNotFoundError:
        # A plain CMake install, as in the source regression, has no metadata.
        distribution = None
    if (distribution is not None
            and distribution.locate_file("najaeda/__init__.py").resolve() != root / "__init__.py"):
        raise RuntimeError("Imported NajaEDA does not match the installed distribution")
    if extension.parent != root:
        raise RuntimeError("NajaEDA imported a native extension outside its package")
    return root, extension, None if distribution is None else distribution.version


def _provider(published: bool = False) -> Provider:
    from najaeda import naja

    root, extension, version = _provider_root()
    if version is None:
        version = naja.getVersion()
    git_commit = naja.getGitHash()
    if published:
        if version != PUBLISHED_VERSION or git_commit != PUBLISHED_GIT_COMMIT:
            raise RuntimeError(
                f"Published NajaEDA adapter requires {PUBLISHED_VERSION} ({PUBLISHED_GIT_COMMIT}); "
                f"found {version} ({git_commit})")
    elif not (root / "sdk/include/naja/core/NajaVersion.h").is_file():
        raise RuntimeError(
            f"NajaEDA {version} ({git_commit}) in {root} ships no headers. Install the "
            "development wheel built from thirdparty/naja, or configure with "
            f"KEPLER_USE_PUBLISHED_NAJAEDA=ON for the published {PUBLISHED_VERSION} release")
    if naja.getVersion() != version:
        raise RuntimeError(
            f"NajaEDA distribution {version} reports native version {naja.getVersion()}")
    return Provider(root, extension, version, git_commit)


def _provider_files(root: Path) -> tuple[Path, ...]:
    roots = (root, root.parent / "najaeda.libs")
    return tuple(sorted({path.resolve() for directory in roots
                         for path in directory.rglob("*") if path.is_file()}))


def _libraries(root: Path) -> dict[str, Path]:
    files = _provider_files(root)
    libraries = {}
    for name in (*RUNTIMES, "tbb", "tbbmalloc"):
        suffix = r"(?:[0-9]+)?(?:[-.][a-zA-Z0-9_]+)*" if name.startswith("tbb") else r"(?:-[a-zA-Z0-9_]+)?"
        pattern = re.compile(
            rf"(?:lib)?{name}{suffix}(?:\.so(?:\.[0-9.]+)?|\.dylib|\.dll)$",
            re.IGNORECASE)
        matches = [path for path in files if pattern.fullmatch(path.name)]
        if len(matches) > 1 or (not matches and name in RUNTIMES):
            raise RuntimeError(f"Expected one {name} provider library, found {matches}")
        if matches:
            libraries[name] = matches[0]
    return libraries


def _wheel_headers(root: Path) -> list[Path]:
    """Use the headers the development wheel installed next to its libraries."""
    include_root = root / "sdk/include"
    include_dirs = {include_root}
    include_dirs.update(path.parent for path in (include_root / "naja").rglob("*.h"))
    if not (include_root / "naja/core/NajaVersion.h").is_file():
        raise RuntimeError("NajaEDA wheel headers contain no version header")
    return sorted(include_dirs)


def _published_headers(output_dir: Path, source_archive: Path | None) -> list[Path]:
    """Extract the published release's headers from its verified source archive."""
    if source_archive is None:
        source_archive = output_dir / f"najaeda-{PUBLISHED_VERSION}.tar.gz"
        if not source_archive.is_file():
            with urllib.request.urlopen(PUBLISHED_SOURCE_URL, timeout=60) as response:
                content = response.read()
            if hashlib.sha256(content).hexdigest() != PUBLISHED_SOURCE_SHA256:
                raise RuntimeError("Published NajaEDA source archive checksum mismatch")
            source_archive.write_bytes(content)
    if hashlib.sha256(source_archive.read_bytes()).hexdigest() != PUBLISHED_SOURCE_SHA256:
        raise RuntimeError("Published NajaEDA source archive checksum mismatch")
    destination = output_dir / "release-headers"
    include_dirs = {output_dir / "include"}
    with tarfile.open(source_archive, "r:gz") as archive:
        for member in archive:
            path = PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts or "\\" in member.name:
                raise RuntimeError(f"Unsafe path in NajaEDA source archive: {member.name}")
            if not path.parts or path.parts[0] != f"najaeda-{PUBLISHED_VERSION}":
                raise RuntimeError("Unexpected NajaEDA source archive root")
            relative = PurePosixPath(*path.parts[1:])
            if not relative.parts:
                continue
            is_source = relative.parts[0] == "src"
            is_dependency = (len(relative.parts) > 2 and relative.parts[0] == "thirdparty"
                             and relative.parts[1] in ("spdlog-1.17.0", "argparse-3.1")
                             and relative.parts[2] == "include")
            if not member.isfile() or not (is_source or is_dependency):
                continue
            if relative.suffix not in (".h", ".hpp") and relative.name != "NajaVersion.h.in":
                continue
            content = archive.extractfile(member).read()
            target = destination.joinpath(*relative.parts)
            if relative.name == "NajaVersion.h.in":
                target = output_dir / "include" / "NajaVersion.h"
                content = content.replace(b"@NAJA_GIT_HASH@", PUBLISHED_GIT_COMMIT.encode("ascii"))
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
            if is_source:
                include_dirs.add(target.parent)
            else:
                include_dirs.add(destination / "thirdparty" / relative.parts[1] / "include")
    if not (output_dir / "include" / "NajaVersion.h").is_file():
        raise RuntimeError("NajaEDA release archive contains no version header")
    return sorted(include_dirs)


def fixup_consumer(consumer: Path) -> None:
    """Resolve repaired Mach-O install IDs without modifying provider files."""
    if sys.platform != "darwin":
        return
    root, _, _ = _provider_root()
    consumer = consumer.resolve(strict=True)
    files = _provider_files(root)
    if consumer in files:
        raise RuntimeError("Consumer fixup must not modify the NajaEDA provider")
    libraries = {}
    for library in files:
        if library.suffix != ".dylib":
            continue
        if library.name in libraries:
            raise RuntimeError(f"Ambiguous provider library: {library.name}")
        libraries[library.name] = library
    linkage = subprocess.check_output(["otool", "-L", str(consumer)], text=True)
    changed = False
    for line in linkage.splitlines()[1:]:
        old = line.strip().split(" (compatibility version", 1)[0]
        provider = libraries.get(Path(old).name)
        if provider is not None and old != str(provider):
            subprocess.run(["install_name_tool", "-change", old, str(provider),
                            str(consumer)], check=True)
            changed = True
    if changed:
        subprocess.run(["codesign", "--force", "--sign", "-", str(consumer)], check=True)


def _pe_exports(library: Path) -> list[tuple[str, bool]]:
    """Read named AMD64 PE exports, including the data/function distinction."""
    data = library.read_bytes()

    def unpack(fmt: str, offset: int) -> tuple[int, ...]:
        if offset < 0 or offset + struct.calcsize(fmt) > len(data):
            raise RuntimeError(f"Truncated PE export table: {library}")
        return struct.unpack_from(fmt, data, offset)

    if data[:2] != b"MZ":
        raise RuntimeError(f"Not a PE library: {library}")
    (pe,) = unpack("<I", 0x3C)
    if data[pe:pe + 4] != b"PE\0\0":
        raise RuntimeError(f"Invalid PE signature: {library}")
    machine, count, _, _, _, optional_size, _ = unpack("<HHIIIHH", pe + 4)
    optional = pe + 24
    if machine != 0x8664 or unpack("<H", optional)[0] != 0x20B:
        raise RuntimeError(f"Kepler requires an AMD64 PE32+ NajaEDA library: {library}")
    export_rva, export_size = unpack("<II", optional + 112)
    sections = []
    for index in range(count):
        section = optional + optional_size + index * 40
        virtual_size, rva, raw_size, raw = unpack("<IIII", section + 8)
        (flags,) = unpack("<I", section + 36)
        sections.append((rva, max(virtual_size, raw_size), raw, raw_size, flags))

    def locate(rva: int, *, require_data: bool = True) -> tuple[int, int]:
        for start, size, raw, raw_size, flags in sections:
            if start <= rva < start + size and (not require_data or rva - start < raw_size):
                return raw + rva - start, flags
        raise RuntimeError(f"Unmapped PE export RVA in {library}")

    offset, _ = locate(export_rva)
    _, _, _, _, _, _, function_count, name_count, functions, names, ordinals = unpack("<IIHHIIIIIII", offset)
    function_table, _ = locate(functions)
    name_table, _ = locate(names)
    ordinal_table, _ = locate(ordinals)
    exports = []
    for index in range(name_count):
        (name_rva,) = unpack("<I", name_table + 4 * index)
        name_offset, _ = locate(name_rva)
        end = data.find(b"\0", name_offset)
        if end < 0:
            raise RuntimeError(f"Unterminated PE export name in {library}")
        name = data[name_offset:end].decode("ascii")
        if not name or any(character in name for character in ('"', "\n", "\r")):
            raise RuntimeError(f"Invalid PE export name in {library}")
        (ordinal,) = unpack("<H", ordinal_table + 2 * index)
        if ordinal >= function_count:
            raise RuntimeError(f"Invalid PE export ordinal in {library}")
        (address,) = unpack("<I", function_table + 4 * ordinal)
        # An exported zero-initialized variable can live in a section whose
        # virtual size exceeds its file-backed bytes. Only its flags are
        # needed here; the export/name tables above still require real bytes.
        _, flags = locate(address, require_data=False)
        forwarded = export_rva <= address < export_rva + export_size
        exports.append((name, not forwarded and not flags & 0x20000020))
    if not exports:
        raise RuntimeError(f"NajaEDA provider exports no named symbols: {library}")
    return exports


def _import_library(name: str, library: Path, output_dir: Path, root: Path | None = None) -> Path:
    # The development wheel ships the original import libraries. They only
    # match DLLs whose names wheel repair left unchanged.
    if root is not None and library.name.lower() == f"{name}.dll":
        original = root / "sdk/lib" / f"{name}.lib"
        if original.is_file():
            return original
    digest = hashlib.sha256(library.read_bytes()).hexdigest()[:20]
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / f"{name}-{digest}.lib"
    if destination.is_file():
        return destination
    dlltool = shutil.which("llvm-dlltool")
    if dlltool is None:
        candidate = Path(os.environ.get("ProgramFiles", "C:/Program Files")) / "LLVM/bin/llvm-dlltool.exe"
        if candidate.is_file():
            dlltool = str(candidate)
    librarian = shutil.which("lib.exe") if dlltool is None else None
    if dlltool is None and librarian is None:
        raise RuntimeError("Repaired NajaEDA DLLs require llvm-dlltool or MSVC lib.exe to build import libraries")
    exports = _pe_exports(library)
    with tempfile.TemporaryDirectory(prefix=f"{name}-", dir=output_dir) as temporary:
        definition = Path(temporary) / "provider.def"
        temporary_library = Path(temporary) / "provider.lib"
        definition.write_text(
            f'LIBRARY "{library.name}"\nEXPORTS\n' + "".join(
                f'  "{symbol}"' + (" DATA" if is_data else "") + "\n"
                for symbol, is_data in exports), encoding="ascii")
        command = ([dlltool, "-m", "i386:x86-64", "-d", str(definition), "-l", str(temporary_library)]
                   if dlltool else [librarian, "/nologo", "/machine:x64", f"/def:{definition}", f"/out:{temporary_library}"])
        subprocess.run(command, check=True, capture_output=True, text=True)
        temporary_library.replace(destination)
    return destination


def _cmake_value(value: object) -> str:
    value = str(value).replace("\\", "/")
    if any(character in value for character in (";", "\n", "\r", "]==]")):
        raise RuntimeError(f"Unsupported character in CMake value: {value!r}")
    return f"[==[{value}]==]"


def cmake_config(output_dir: Path, *, published: bool = False,
                 source_archive: Path | None = None) -> str:
    provider = _provider(published)
    libraries = _libraries(provider.root)
    output_dir = output_dir.resolve()
    (output_dir / "include").mkdir(parents=True, exist_ok=True)
    if published:
        include_dirs = _published_headers(output_dir, source_archive)
    else:
        include_dirs = [output_dir / "include", *_wheel_headers(provider.root)]
    manifest = {}
    for path in (*libraries.values(), provider.extension):
        try:
            relative = path.relative_to(provider.root.parent).as_posix()
        except ValueError as error:
            raise RuntimeError(f"Provider library is outside its installation: {path}") from error
        manifest[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_json = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    # The extension validates these at import: the same release and the same
    # native files it was compiled against. Package versions alone do not
    # identify a C++ ABI.
    (output_dir / "include" / "KeplerNajaProviderBuild.h").write_text(
        "// Generated from the installed NajaEDA provider.\n#pragma once\n"
        f"#define KEPLER_NAJA_PROVIDER_VERSION {json.dumps(provider.version)}\n"
        f"#define KEPLER_NAJA_PROVIDER_GIT_HASH {json.dumps(provider.git_commit)}\n"
        f"#define KEPLER_NAJA_PROVIDER_MANIFEST {json.dumps(manifest_json)}\n",
        encoding="utf-8")
    values = {"NajaEDA_VERSION": provider.version, "NajaEDA_GIT_COMMIT": provider.git_commit,
              "NajaEDA_BUILD_ID": hashlib.sha256(manifest_json.encode()).hexdigest()}
    lines = [f"set({name} {_cmake_value(value)})" for name, value in values.items()]
    lines.append("set(NajaEDA_INCLUDE_DIRS " + " ".join(map(_cmake_value, include_dirs)) + ")")
    lines.append("set(NajaEDA_RUNTIME_LIBRARIES " + " ".join(map(_cmake_value, libraries.values())) + ")")
    for name, library in libraries.items():
        lines.append(f"set(NajaEDA_{name}_LIBRARY {_cmake_value(library)})")
        if sys.platform == "win32":
            implib = _import_library(name, library, output_dir / "import-libs",
                                     None if published else provider.root)
            lines.append(f"set(NajaEDA_{name}_IMPLIB {_cmake_value(implib)})")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--cmake", action="store_true")
    action.add_argument("--fixup-consumer", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--published", action="store_true",
                        help=f"use the published NajaEDA {PUBLISHED_VERSION} release")
    parser.add_argument("--source-archive", type=Path,
                        help="cached source archive of the published release")
    args = parser.parse_args()
    if args.cmake and args.output_dir is None:
        parser.error("--cmake requires --output-dir")
    if args.source_archive is not None and not args.published:
        parser.error("--source-archive requires --published")
    if args.cmake:
        print(cmake_config(args.output_dir, published=args.published,
                           source_archive=args.source_archive), end="")
    else:
        fixup_consumer(args.fixup_consumer)


if __name__ == "__main__":
    main()
