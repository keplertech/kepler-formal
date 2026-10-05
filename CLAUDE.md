# kepler-formal — Agent Guide

Guidance for AI coding agents (Claude Code reads this as `CLAUDE.md`;
Codex reads it as `AGENTS.md` via symlink) working in this repository.

## Build systems

CMake (`CMakeLists.txt`, `thirdparty/` git submodules) and Bazel
(`MODULE.bazel`, `BUILD.bazel` files) are independent. The Bazel build
does not use the submodules, and Bazel changes must not touch the CMake
flow.

## Bazel: a BCR-ready module

kepler-formal's Bazel build is meant to be published to the Bazel Central
Registry (BCR) and consumed by other modules (e.g. bazel-orfs) with a
plain `bazel_dep`. Keep it that way. Invariants:

- `MODULE.bazel` contains only `bazel_dep`s, plus the hermetic-llvm
  toolchain setup, whose registration is a `dev_dependency`. No
  `http_archive`, `git_override`, `archive_override`, module extensions
  or repository rules of our own; overrides are ignored for non-root
  modules, so they would only hide breakage that consumers then hit.
- Nothing runs cmake/make inside Bazel (no `rules_foreign_cc`), and
  nothing is found on the host (`PATH`, pkg-config, system headers or
  libraries). CI installs no packages for Bazel.
- Depend on naja's own targets (`@naja//src/dnl:naja_dnl`, …). If naja
  (or any dependency) needs a fix, fix it there — in naja itself or in
  its registry entry — never by patching its BUILD files from here.
- Load every rule (`@rules_cc//cc:cc_library.bzl`,
  `@rules_shell//shell:sh_test.bzl`, …); Bazel 9 has no native ones.

**Where dependencies come from.** Everything comes from BCR. A module
version not on BCR yet comes from its open bazel-central-registry pull
request (one module per PR): `.bazelrc` lists the PR's commit
(`https://raw.githubusercontent.com/<fork>/bazel-central-registry/<sha>/`)
ahead of BCR, and Bazel takes each `name@version` from the first
registry that has it. Drop the line once the PR merges; when a PR needs
a fix, push a new commit (no force-push, pinned commits must stay
reachable) and move the pin. The one exception is naja's development
pin (#457 + #455 merged, `oharboe/naja@kepler-pin`), served from the
in-tree `bazel/registry/` until naja is released and on BCR; bump it with
`bazel/registry/update_source.py naja <version> <url> [<strip_prefix>]`.
Unreleased commits use `<release>-<YYYYMMDD>-<commit>` versions. Never
change a version's contents once something depends on it; add a new
version.

To try a local naja checkout without touching the registry:
`bazel test //... --override_module=naja=<path>`. Don't commit such
overrides.

**Publishing to BCR**: one module per bazel-central-registry pull
request, each based on BCR main, so each goes in on its own. Then naja
(publish-to-bcr from a naja release), then kepler-formal itself through
the publish-to-bcr app (`.bcr/` templates; maintainers: xtofalex,
nanocoh). BCR rejects symlinks in entries, so `overlay/MODULE.bazel` is a
copy. See `docs/bcr-roadmap.md`.

**Known traps** (each already fixed; don't reintroduce):

- The hermetic-llvm toolchain passes libc headers as early `-isystem`
  flags; libraries whose headers must shadow libc's need `-I` (see the
  registry's bison). It also exposes headers that a host-toolchain build
  silently takes from `/usr/include`.
- `kepler-formal` links naja's `naja_runtime` as a `cc_shared_library`
  (`dynamic_deps`), shared with the `naja.so` Python module. Anything it
  links that the binary also uses must be exported by `naja_runtime`
  (TBB, zlib), and linker inputs must be owned by targets the
  `cc_shared_library` graph aspect can see — naja's `python_libs` exists
  for libpython.
- The embedded Python is rules_python's hermetic runtime. Its compiled-in
  prefix does not exist at run time: `src/bin/CppDriver.cpp` points
  `PYTHONHOME` at the runtime in runfiles (Bazel builds only), and the
  release tarball bundles it and sets `PYTHONHOME` in its wrapper.
- `.github/consumer-test` is a separate workspace (in `.bazelignore`)
  that builds kepler-formal as a dependency, on the newest Bazel. Run it
  for any change to `MODULE.bazel`, the registry, or labels that cross
  repositories.

## Build & test

```bash
bazelisk build //src/bin:kepler-formal
bazelisk test //...
(cd .github/consumer-test && bazelisk build @kepler-formal//src/bin:kepler-formal)
```
