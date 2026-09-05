# Bazel Central Registry (BCR) Roadmap

This document tracks the plan to publish kepler-formal to the
[Bazel Central Registry](https://registry.bazel.build/) so that
downstream projects like [bazel-orfs](https://github.com/The-OpenROAD-Project/bazel-orfs)
can depend on it with a simple `bazel_dep()`.

## Build systems

kepler-formal supports two build systems:

- **CMake** (primary): Uses git submodules in `thirdparty/` directly.
  This is the established flow and remains a first-class citizen.
- **Bazel** (new): Uses `http_archive` to fetch the same dependencies
  as source archives from GitHub. This avoids requiring git submodules
  and is compatible with BCR distribution.

Both build systems use the same source code and overlay BUILD files
(in `bazel/`). The cmake flow is unaffected by Bazel changes.

### Current Bazel usage

Bazel build support is available alongside CMake. It uses [bzlmod](https://bazel.build/external/overview#bzlmod) (`MODULE.bazel`) and pulls most dependencies from the [Bazel Central Registry](https://registry.bazel.build/).

Build with Bazel:

```bash
git clone --recurse-submodules https://github.com/keplertech/kepler-formal.git
cd kepler-formal
bazel build //src/bin:kepler-formal
```

Run tests:

```bash
bazel test //test/...
```

### Dependency strategy

- **BCR**: `yaml-cpp`, `googletest`, `zlib`, `spdlog`
- **Native BUILD files**: `kissat`, `glucose`, `boost` (headers),
  `FlexLexer.h`
- **`rules_foreign_cc`**: `naja` and its nested sources, `capnproto`,
  `onetbb` (cmake package trees consumed by naja's `find_package`)
- **System packages still required**: build tools only (cmake, make,
  bison, flex) and Python 3 headers — see the hermeticity inventory below

### Future work: bazel-orfs integration

Once kepler-formal is fully buildable with Bazel, it can be consumed directly by [bazel-orfs](https://github.com/The-OpenROAD-Project/bazel-orfs) as a proper Bazel dependency instead of the current `$PATH`-based wrapper
([bazel-orfs#523](https://github.com/The-OpenROAD-Project/bazel-orfs/pull/523)).
This enables:

- **`bazel_dep` or `git_override`** in bazel-orfs `MODULE.bazel` to pin kepler-formal to a specific version or commit
- **Hermetic LEC tests** in bazel-orfs CI without requiring a pre-installed kepler-formal binary
- **Remote caching** of kepler-formal build artifacts shared across bazel-orfs users
- **Eventual BCR publication** of kepler-formal as a first-class Bazel module

## Phases

### Phase 1: http_archive deps (done)

The Bazel build fetches glucose, kissat, and naja via `http_archive`
instead of `new_local_repository`. This means `bazel build` works
from a source tarball without `git submodule init`.

Dependency versions are pinned in `bazel/deps.bzl` by commit SHA.

### Phase 2: git_override (next)

Downstream projects can depend on kepler-formal via `git_override`
in their `MODULE.bazel`:

```python
bazel_dep(name = "kepler-formal", version = "1.0.0")

git_override(
    module_name = "kepler-formal",
    remote = "https://github.com/keplertech/kepler-formal.git",
    commit = "<sha>",
)
```

This works today — no BCR submission required.

### Phase 2.5: Binary releases (done)

Pre-built binaries are published as GitHub Releases.  Run
`bazelisk run //:release` to tag and push; CI builds an optimized binary
and uploads it.  See `docs/releasing.md` for full instructions.

Cap'n Proto is statically linked.  Naja shared libraries and TBB
(which has no static libs) are bundled in the tarball with a wrapper
script.

Downstream projects can fetch the tarball via `http_archive` for fast
CI without building from source.

### Phase 3: BCR publication

To publish to BCR:

1. Ensure the version in `MODULE.bazel` follows semver (currently `1.0.0`).
2. Create a GitHub release/tag matching the version (e.g. `v1.0.0`) —
   this now happens automatically via `bazelisk run //:release`.
3. Submit to BCR via the
   [publish-to-bcr](https://github.com/bazel-contrib/publish-to-bcr) GitHub App
   or manual PR to [bazel-central-registry](https://github.com/bazelbuild/bazel-central-registry).
4. Once published, downstream projects use plain `bazel_dep` with no override.

## Updating pinned dependencies

Dependency commit SHAs are pinned in `bazel/deps.bzl`. To update:

1. Update the `_*_COMMIT` constants to the new commit SHA.
2. Run `bazel fetch @glucose @kissat @naja` — it will fail with a SHA
   mismatch and print the correct `sha256`. Update the hash.
3. Verify: `bazel build //src/bin:kepler-formal && bazel test //test/...`

The cmake flow uses whatever submodule commit is checked out in
`thirdparty/`. Keep the Bazel pins in sync with the submodule SHAs
when updating submodules.

## Hermeticity: local-system dependency inventory

All C/C++ *libraries* are now provided through Bazel; what remains from the
host is build *tools*. Full inventory with per-item status:

| Dependency | Where used | Status / hermetic path |
|---|---|---|
| libcapnp/libkj + capnp compiler | src/bin, tests, naja cmake | **hermetic**: cmake-built `@capnproto` 1.4.0 (`bazel/capnproto.BUILD.bazel`) |
| libtbb/libtbbmalloc | src/bin, tests, naja cmake, release bundle | **hermetic**: cmake-built `@onetbb` 2022.3.0 — naja's `.so`s and the binary share one `libtbb.so.12` |
| zlib | src/bin, tests, naja cmake, glucose | **hermetic**: BCR `@zlib` everywhere |
| Boost headers | naja cmake + kepler compiles (via naja headers) | **hermetic**: `@boost_headers` 1.89.0 |
| FlexLexer.h (was libfl-dev) | naja-verilog scanner (previously leaked `-I/usr/include` via FindFLEX) | **hermetic**: `@flex_src//:flexlexer` + `FLEX_INCLUDE_DIR` cache entry |
| Host C++ toolchain | whole Bazel build | **hermetic** (Linux): hermetic-llvm (BCR `llvm`, zero-sysroot clang 22/libc++, glibc 2.28 floor) — see OpenROAD PR 10812. macOS still uses the host Xcode toolchain |
| Python3 interpreter + dev headers | naja cmake `find_package(Python3 ... Development.Embed REQUIRED)` runs even with `BUILD_NAJA_PYTHON=OFF`; `libnaja_python.so` links system libpython | future: rules_python / python-build-standalone + `Python3_ROOT_DIR`, or an upstream naja option to skip the probe |
| bison/flex executables | naja-verilog codegen (Homebrew paths hardcoded on macOS) | future: BCR `bison`/`flex`/`m4` modules fed via `BISON_EXECUTABLE`/`FLEX_EXECUTABLE` |
| cmake/make/pkg-config binaries | rules_foreign_cc preinstalled toolchains | future: rules_foreign_cc built/prebuilt toolchains (built pkg-config currently trips over a glib C23 issue with GCC 15) |
| Network during naja build | slang FetchContent (fmt) — the `requires-network` tag on `@naja` | future: pre-fetch in `naja_repo`, set `FETCHCONTENT_SOURCE_DIR_FMT` |
| macOS host toolchain | macOS Bazel build | future: hermetic-llvm macOS (SDK subset; CLI defaults suffice) |
| clang-tidy | not wired | future: `@llvm//tools:clang-tidy` lint aspect |

The BCR `onetbb` module is Bazel-8-compatible nowadays; the cmake-built
`@onetbb` is used instead because naja's `find_package(TBB)` needs
`TBBConfig.cmake` and both naja and the binary must share a single TBB
runtime.

### Hermetic dependencies for the CMake flow: `bazelisk run //:deps`

The same install trees serve the plain CMake flow, opt-in:

```bash
bazelisk run //:deps
cmake -B build -DCMAKE_TOOLCHAIN_FILE=$PWD/deps/toolchain.cmake
```

This exports capnproto, oneTBB, boost, zlib and FlexLexer.h into `deps/`
(git- and bazel-ignored) with a generated `toolchain.cmake`, replacing the
platform-specific apt/homebrew install step. On Linux the hermetic
compiler is exported too (`deps/llvm/`, with `toolchain.cmake` generated
from the resolved Bazel cc_toolchain by `bazel/toolchain_cmake.bzl`), so
the CMake flow compiles with the identical clang/libc++ as the Bazel flow
and needs no host compiler at all. The system-package CMake flow keeps
working unchanged. Host tools still required: cmake, make, bison, flex,
python3 (dev headers).
