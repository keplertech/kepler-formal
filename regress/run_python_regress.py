#!/usr/bin/env python3
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Build NajaEDA and Kepler Python from this checkout and run their regression.

Uses direct CMake installs in a private directory: no pip, wheels, or publishing.
Run with the Python interpreter you want to build and test (Python >= 3.10).
"""

import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import venv


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build/python-regress")
    parser.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument(
        "--cmake-arg", action="append", default=[], metavar="ARG",
        help="extra option for both CMake configures; use --cmake-arg=-DNAME=VALUE",
    )
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    for tool in ("cmake", "ctest", "ninja"):
        if shutil.which(tool) is None:
            parser.error(f"{tool} is required on PATH")
    if not (ROOT / "thirdparty/naja/CMakeLists.txt").is_file():
        parser.error("initialize the pinned sources with git submodule update --init --recursive")

    build = args.build_dir.resolve()
    build.mkdir(parents=True, exist_ok=True)
    stage = build / "packages"
    environment = build / "venv"
    # Keep system and user Python packages out of the test environment. Both
    # projects use this exact interpreter, compiler configuration, and SDK.
    venv.EnvBuilder(with_pip=False).create(environment)
    python = environment / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    env = dict(os.environ, PYTHONPATH=str(stage), PYTHONNOUSERSITE="1",
               PYTHONFAULTHANDLER="1", PYTHONUNBUFFERED="1")
    env.pop("PYTHONHOME", None)

    log_path = build / "regression.log"
    with log_path.open("w", encoding="utf-8") as log:
        def run(*command):
            command = [str(item) for item in command]
            line = "+ " + shlex.join(command)
            print(line, flush=True)
            log.write(line + "\n")
            with subprocess.Popen(
                command, cwd=build, env=env, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, errors="replace",
            ) as process:
                for line in process.stdout:
                    print(line, end="", flush=True)
                    log.write(line)
                log.flush()
                if process.wait():
                    raise subprocess.CalledProcessError(process.returncode, command)

        common = [
            "-G", "Ninja", *args.cmake_arg,
            "-DCMAKE_BUILD_TYPE=Release", f"-DPython3_EXECUTABLE={python}",
            f"-DCMAKE_INSTALL_PREFIX={stage}",
        ]
        try:
            # Fail early rather than reaching the SDK-selection regression
            # after compiling with an older CMake.
            version_check = build / "check_cmake.cmake"
            version_check.write_text("cmake_minimum_required(VERSION 3.30)\n")
            run("cmake", "-P", version_check)
            provider = build / "naja"
            run("cmake", "-S", ROOT / "thirdparty/naja", "-B", provider,
                *common, "-DBUILD_NAJA_PYTHON=ON", "-DPREGENERATED_PARSER_SOURCES=OFF")
            run("cmake", "--build", provider, "--target", "naja", "--parallel", args.jobs)
            run("cmake", "--install", provider)

            # Import before configuring Kepler so a missing provider cannot
            # silently fall back to an unrelated installation.
            run(python, "-c", """
import sys
from pathlib import Path
import najaeda
from najaeda import sdk
stage = Path(sys.argv[1]).resolve()
assert Path(najaeda.__file__).resolve().is_relative_to(stage), najaeda.__file__
assert Path(sdk.get_cmake_dir()).resolve().is_relative_to(stage), sdk.get_cmake_dir()
print('Source-built NajaEDA:', najaeda.__file__)
""", stage)

            consumer = build / "kepler"
            run("cmake", "-S", ROOT, "-B", consumer, *common,
                "-DBUILD_KEPLER_PYTHON=ON", "-DENABLE_UNIT_TESTS=ON")
            run("cmake", "--build", consumer, "--target", "kepler_formal_native",
                "kepler-borrowed-native-tests", "--parallel", args.jobs)
            run("cmake", "--install", consumer, "--component", "python")
            if sys.platform == "darwin":
                # The SDK signs the build artifact. CMake then adjusts its
                # install RPATH, invalidating that signature on Apple Silicon.
                for extension in (stage / "kepler_formal").glob("_native*.so"):
                    run("codesign", "--force", "--sign", "-", extension)
            run(python, "-c", """
import sys
from pathlib import Path
import najaeda
from najaeda import naja
import kepler_formal
from kepler_formal import _native, from_najaeda
stage = Path(sys.argv[1]).resolve()
for module in (najaeda, naja, kepler_formal, _native):
    assert Path(module.__file__).resolve().is_relative_to(stage), module.__file__
    print(module.__name__ + ':', module.__file__)
assert kepler_formal.najaeda is najaeda
universe = naja.NLUniverse.create()
database = naja.NLDB.create(universe)
library = naja.NLLibrary.create(database, 'runtime_check')
design = naja.SNLDesign.create(library, 'borrowed')
assert from_najaeda(design).najaeda_design is design
print('Shared runtime and borrowed design check passed.')
""", stage)

            run("ctest", "--test-dir", consumer, "--output-on-failure",
                "-R", "^kepler-formal-borrowed-native-tests$")
            # CTest's Python fixture replaces PYTHONPATH. Test the installed
            # sibling packages directly instead, including the real example.
            run(python, "-m", "unittest", "discover", "-v", "-s", ROOT / "test/python")
        except subprocess.CalledProcessError as error:
            print(f"Regression failed (exit {error.returncode}); see {log_path}", file=sys.stderr)
            return 1

    print(f"Source-built Python regression passed. Packages: {stage}")
    print(f"Log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
