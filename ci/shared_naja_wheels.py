#!/usr/bin/env python3
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Build a local SDK provider for CI, and repair consumers without copying it.

The provider is a separate wheel. This helper never publishes either package.
After the matching NajaEDA release exists, ordinary isolated pip builds can
resolve it from the package index instead of using these CI-local artifacts.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

PROVIDER_REQUIREMENT = "najaeda==0.7.24.dev0"


def run(*arguments: str, env: dict[str, str] | None = None) -> None:
    subprocess.run(arguments, check=True, env=env)


def provider_files() -> dict[str, Path]:
    distribution = importlib.metadata.distribution("najaeda")
    files: dict[str, Path] = {}
    for item in distribution.files or ():
        path = Path(distribution.locate_file(item)).resolve()
        if path.suffix.lower() in (".dll", ".dylib", ".so") or ".so." in path.name:
            if path.name in files and files[path.name] != path:
                raise RuntimeError(f"Ambiguous provider library: {path.name}")
            files[path.name] = path
    if not files:
        raise RuntimeError("NajaEDA provider contains no native libraries")
    return files


def rewrite_consumer(wheel: Path, destination: Path, dependencies: dict[str, Path],
                     *, absolute: bool) -> Path:
    """Relocate provider references; wheel pack regenerates RECORD hashes."""
    provider_root = Path(importlib.metadata.distribution("najaeda").locate_file("")).resolve()
    with tempfile.TemporaryDirectory(prefix="kepler-wheel-links-") as temporary:
        root = Path(temporary)
        run(sys.executable, "-m", "wheel", "unpack", str(wheel), "--dest", str(root))
        unpacked, = [path for path in root.iterdir() if path.is_dir()]
        for native in unpacked.rglob("*"):
            if not native.is_file() or native.suffix not in (".so", ".dylib"):
                continue
            if platform.system() == "Darwin":
                linkage = subprocess.check_output(["otool", "-L", str(native)], text=True)
                changed = False
                for line in linkage.splitlines()[1:]:
                    old = line.strip().split(" (compatibility version", 1)[0]
                    provider = dependencies.get(Path(old).name)
                    if provider is None:
                        continue
                    relative = os.path.relpath(provider.relative_to(provider_root),
                                               native.relative_to(unpacked).parent)
                    new = str(provider) if absolute else "@loader_path/" + relative
                    if old != new:
                        run("install_name_tool", "-change", old, new, str(native))
                        changed = True
                if changed:
                    run("codesign", "--force", "--sign", "-", str(native))
            elif platform.system() == "Linux" and not absolute:
                old = subprocess.check_output(["patchelf", "--print-rpath", str(native)],
                                              text=True).strip()
                paths = {path for path in old.split(":") if path.startswith("$ORIGIN")}
                for provider in dependencies.values():
                    relative = os.path.relpath(provider.relative_to(provider_root).parent,
                                               native.relative_to(unpacked).parent)
                    paths.add("$ORIGIN/" + relative)
                run("patchelf", "--set-rpath", ":".join(sorted(paths)), str(native))
        destination.mkdir(parents=True, exist_ok=True)
        run(sys.executable, "-m", "wheel", "pack", str(unpacked), "--dest-dir", str(destination))
    return destination / wheel.name


def repair(wheel: Path, destination: Path, *, provider: bool = False) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    dependencies = {} if provider else provider_files()
    system = platform.system()
    if dependencies and system in ("Darwin", "Linux"):
        # Resolve excluded provider libraries while repairing other dependencies,
        # then point the finished wheel at its installed sibling distribution.
        with tempfile.TemporaryDirectory(prefix="kepler-repair-") as temporary:
            root = Path(temporary)
            prepared = rewrite_consumer(wheel, root / "input", dependencies, absolute=True)
            repair_external(prepared, root / "repaired", dependencies)
            repaired, = (root / "repaired").glob("*.whl")
            rewrite_consumer(repaired, destination, dependencies, absolute=False)
    else:
        repair_external(wheel, destination, dependencies)


def repair_external(wheel: Path, destination: Path, dependencies: dict[str, Path]) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    system = platform.system()
    if system == "Linux":
        command = [sys.executable, "-m", "auditwheel", "repair", "--plat",
                   "manylinux_2_28_" + platform.machine(), "-w", str(destination)]
        for name in dependencies:
            command += ["--exclude", name]
        run(*command, str(wheel))
    elif system == "Darwin":
        # delocate resolves the provider's existing install names before it
        # filters dependencies. Give it the installed provider directories;
        # never suppress missing unrelated dependencies.
        environment = os.environ.copy()
        environment["DYLD_LIBRARY_PATH"] = os.pathsep.join(sorted({
            str(path.parent) for path in dependencies.values()
        }))
        command = [sys.executable, "-m", "delocate.cmd.delocate_wheel",
                   "--require-archs", "arm64", "-w", str(destination)]
        for name in dependencies:
            command += ["--exclude", name]
        run(*command, str(wheel), env=environment)
    elif system == "Windows":
        command = [sys.executable, "-m", "delvewheel", "repair",
                   "--ignore-existing", "--analyze-existing", "--include-imports",
                   "-w", str(destination)]
        extra = Path(os.environ["USERPROFILE"]) / "vcpkg/installed/x64-windows/bin"
        command += ["--add-path", os.pathsep.join([
            str(extra), *(str(path.parent) for path in dependencies.values())])]
        if dependencies:
            command += ["--exclude", os.pathsep.join(dependencies)]
        else:
            # Keep the provider's SDK import-library DLL names stable.
            command += ["--no-mangle", "naja_*.dll;libnaja_*.dll"]
        run(*command, str(wheel))
    else:
        raise RuntimeError(f"Unsupported wheel platform: {system}")


def build_provider(project: Path) -> None:
    run(sys.executable, "-m", "pip", "install", "scikit-build-core>=0.11.3,<0.12",
        "build", "wheel")
    repair_package = {"Linux": "auditwheel", "Darwin": "delocate",
                      "Windows": "delvewheel"}[platform.system()]
    run(sys.executable, "-m", "pip", "install", repair_package)
    destination = project / ".kepler-provider-wheels"
    if ".dev" not in PROVIDER_REQUIREMENT:
        # Release consumers must link the exact distributed provider, not a
        # locally rebuilt same-version runtime with a different native build
        # identity. Download it for the isolated cibuildwheel test environment.
        destination.mkdir(parents=True, exist_ok=True)
        run(sys.executable, "-m", "pip", "download", "--only-binary=:all:",
            "--no-deps", "--dest", str(destination), PROVIDER_REQUIREMENT)
        install_provider(project)
        return
    with tempfile.TemporaryDirectory(prefix="kepler-naja-wheel-") as temporary:
        environment = os.environ.copy()
        if sys.platform == "win32":
            toolchain = Path(os.environ["USERPROFILE"]) / "vcpkg/scripts/buildsystems/vcpkg.cmake"
            environment["CMAKE_ARGS"] = (
                f'-DCMAKE_TOOLCHAIN_FILE="{toolchain.as_posix()}" '
                '-DPREGENERATED_PARSER_SOURCES=ON')
        else:
            environment["CMAKE_ARGS"] = (
                environment.get("CMAKE_ARGS", "") + " -DPREGENERATED_PARSER_SOURCES=OFF")
        run(sys.executable, "-m", "pip", "wheel", "--no-build-isolation", "--no-deps",
            "--wheel-dir", temporary, str(project / "thirdparty/naja"), env=environment)
        wheels = list(Path(temporary).glob("najaeda-*.whl"))
        if len(wheels) != 1:
            raise RuntimeError(f"Expected one provider wheel, found {wheels}")
        repair(wheels[0], destination, provider=True)
    install_provider(project)


def install_provider(project: Path) -> None:
    run(sys.executable, "-m", "pip", "install", "--no-index", "--force-reinstall",
        "--find-links", str(project / ".kepler-provider-wheels"), PROVIDER_REQUIREMENT)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build-provider", "install-provider", "repair"))
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--wheel", type=Path)
    parser.add_argument("--dest", type=Path)
    args = parser.parse_args()
    if args.action == "build-provider":
        build_provider(args.project)
    elif args.action == "install-provider":
        install_provider(args.project)
    else:
        if args.wheel is None or args.dest is None:
            parser.error("repair requires --wheel and --dest")
        repair(args.wheel, args.dest)


if __name__ == "__main__":
    main()
