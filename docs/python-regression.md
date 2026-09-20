# Python regression from source

The regression builds NajaEDA from the pinned `thirdparty/naja` submodule,
then builds Kepler's Python extension against that NajaEDA SDK. Both packages
are installed with CMake into a private directory. No wheels, package index,
or publishing are involved.

From the repository root on Linux or macOS:

```bash
git submodule update --init --recursive
python3 regress/run_python_regress.py --jobs 4
```

Use Python 3.10 or newer, CMake 3.30 or newer, Ninja, a C/C++ compiler, Boost,
TBB, Cap'n Proto, zlib, Flex, and Bison. On Ubuntu 24.04, the native packages are:

```bash
sudo apt-get install build-essential ninja-build clang pkg-config \
  libboost-dev libboost-iostreams-dev libtbb-dev capnproto libcapnp-dev \
  zlib1g-dev flex bison libfl-dev python3-dev python3-venv
```

Install CMake 3.30 or newer separately if the system version is older.
On macOS, use the same native dependencies as the normal CMake build;
Homebrew Flex and Bison must be on `PATH`. Extra CMake options are passed
to **both** builds, for example:

```bash
python3 regress/run_python_regress.py --jobs 4 \
  --cmake-arg=-DCMAKE_PREFIX_PATH=/opt/homebrew
```

`CC` and `CXX` can select the compiler for a fresh build. Use a fresh
`--build-dir PATH` when changing Python, compiler, or native dependency versions.
CMake may download its third-party build dependencies on the first run.

The runner creates `build/python-regress/venv` without pip and installs both
packages under `build/python-regress/packages`. It checks the actual imported
package paths and passes a live NajaEDA design to Kepler before testing. Kepler's
native API validates that the packages share the same runtime and universe.

It runs the native borrowed-design tests and the full Python test suite,
including the actual [TinyRocket example](../examples/tinyrocket/verify_python.py):

- Original versus edited: must report `different` and exit 1.
- Original loaded into two databases: must report `equivalent` and exit 0.
- Both designs must remain accessible in NajaEDA after verification.

The overall regression exits 0 only when all tests pass. Build and test output
is saved in `build/python-regress/regression.log`.

To run the example yourself with these compiled packages, from the repository
root:

```bash
export PYTHONPATH="$PWD/build/python-regress/packages"
build/python-regress/venv/bin/python examples/tinyrocket/verify_python.py
# Exit 1 is expected: this pair is intentionally different.

build/python-regress/venv/bin/python examples/tinyrocket/verify_python.py \
  examples/tinyrocket/tinyrocket.v examples/tinyrocket/tinyrocket.v
# Exit 0 is expected: both loaded designs are equivalent.
```

The [Python source regression workflow](../.github/workflows/regress-python.yml)
runs this same runner on pushes and pull requests. Wheel builds and release
workflows remain separate.
