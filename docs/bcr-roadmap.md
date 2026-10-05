# Bazel and the Bazel Central Registry (BCR)

kepler-formal's Bazel build is a plain bzlmod module, ready to be
published to the [Bazel Central Registry](https://registry.bazel.build/)
so that downstream projects like
[bazel-orfs](https://github.com/The-OpenROAD-Project/bazel-orfs) can
depend on it with a single `bazel_dep()`.

The Bazel build is independent of the CMake flow: it does not use the
`thirdparty/` git submodules, and it never runs cmake.

## Build

```bash
bazelisk build //src/bin:kepler-formal
bazelisk test //...
```

No host packages are needed on Linux: the C++ toolchain (hermetic-llvm),
every library, and every code generator (bison, flex, the capnp compiler,
slang's Python generators) come from Bazel modules.

## Dependencies

`MODULE.bazel` contains only `bazel_dep`s. Everything on BCR comes from
BCR: `onetbb`, `spdlog`, `yaml-cpp`, `googletest`, and, through naja,
`capnp-cpp`, `boost.*`, `fmt`, `bison`, `flex`, `rules_python`, and
others.

What is not on BCR yet comes from its open bazel-central-registry pull
request, which `.bazelrc` lists by commit ahead of BCR:

| Module | BCR pull request |
|---|---|
| `bison` 3.8.2.bcr.10 | bazelbuild/bazel-central-registry#10879 |
| `sv-lang` 11.0.0-20260701-b60d729d.bcr.1 | bazelbuild/bazel-central-registry#10882 |
| `naja-if` | bazelbuild/bazel-central-registry#10881 |
| `naja-verilog` | bazelbuild/bazel-central-registry#10885 |
| `kissat` 4.0.4 | bazelbuild/bazel-central-registry#10880 |
| `cadical` 3.0.0 | bazelbuild/bazel-central-registry#10883 |
| `glucose` 4.2.1-20251230-674dbba | bazelbuild/bazel-central-registry#10884 |

`naja` itself is served from the in-tree
[`bazel/registry/`](../bazel/registry/README.md) until it is released and
on BCR.

## Depending on kepler-formal

Until kepler-formal and the modules above are on BCR, a downstream
module lists the same registries as kepler-formal's `.bazelrc`, with
kepler-formal's `bazel/registry/` by pinned URL:

```
# .bazelrc
common --registry=https://raw.githubusercontent.com/keplertech/kepler-formal/<commit>/bazel/registry/
# ...the BCR pull request registries from kepler-formal's .bazelrc...
common --registry=https://bcr.bazel.build/
build --cxxopt=-std=c++20
```

```python
# MODULE.bazel
bazel_dep(name = "kepler-formal", version = "0.5.0")
git_override(
    module_name = "kepler-formal",
    remote = "https://github.com/keplertech/kepler-formal.git",
    commit = "<commit>",
)
```

`.github/consumer-test` is such a consumer and runs in CI.

## Publishing to BCR

1. Get the BCR pull requests above merged, dropping each `.bazelrc`
   line as it lands, then naja (released, via publish-to-bcr), deleting
   it from `bazel/registry/`.
2. Ensure the version in `MODULE.bazel` matches
   `src/bin/KeplerVersion.h.in`, and tag a release with
   `bazelisk run //:release` (see `docs/releasing.md`).
3. Submit kepler-formal with the
   [publish-to-bcr](https://github.com/bazel-contrib/publish-to-bcr)
   GitHub App, which uses the templates in `.bcr/`.

## Binary releases

`bazelisk run //:release` tags a version; CI builds an optimized binary
and publishes it as a GitHub Release. See `docs/releasing.md`.
