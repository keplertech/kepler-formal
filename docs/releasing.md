# Releasing kepler-formal

This document explains how to set up and perform binary releases.

## How it works

1. A maintainer runs `bazelisk run //:release` locally.
2. The script validates that MODULE.bazel and CMakeLists.txt have the same
   version, creates an annotated git tag (`v1.0.0`), and pushes it.
3. GitHub Actions (`.github/workflows/release.yml`) triggers on the tag,
   builds an optimized binary with `bazelisk build -c opt`, packages it
   into a tarball with bundled naja shared libraries, and creates a GitHub
   Release with the tarball attached.

The binary statically links libstdc++, libgcc, TBB, and Cap'n Proto.
Only the naja shared libraries (which naja builds as explicitly `SHARED`)
are bundled alongside the binary.  A small wrapper script sets
`LD_LIBRARY_PATH` so the tarball is relocatable.

## One-time GitHub setup

No special setup is required beyond the defaults:

- The release workflow uses `${{ github.token }}` (the built-in
  `GITHUB_TOKEN`), which has write access to repository contents by
  default.  No secrets to configure.
- The workflow triggers on tag push (`v*`), which is enabled by default.

If your repository has branch protection rules that restrict who can push
tags, ensure that the maintainers performing releases have the necessary
permissions.

## Performing a release

### 1. Bump the version

Edit both files to the new version:

```
MODULE.bazel:  version = "1.1.0"
CMakeLists.txt: VERSION 1.1.0
```

Commit the version bump:

```bash
git add MODULE.bazel CMakeLists.txt
git commit -m "Bump version to 1.1.0"
git push origin main
```

### 2. Create the release

```bash
bazelisk run //:release
```

The script will:
- Verify MODULE.bazel and CMakeLists.txt versions match
- Check that the working tree is clean
- Check that the tag doesn't already exist
- Create and push `v1.1.0`

### 3. Wait for CI

The GitHub Actions workflow will:
- Build `//src/bin:kepler-formal` with `-c opt`
- Package the binary with bundled naja shared libraries
- Create a GitHub Release at
  `https://github.com/keplertech/kepler-formal/releases/tag/v1.1.0`
  with auto-generated release notes

### 4. Verify

Download the tarball and test on a clean system:

```bash
tar xzf kepler-formal-1.1.0-linux-x86_64.tar.gz
./kepler-formal-1.1.0-linux-x86_64/kepler-formal --help
```

## Versioning

The project uses [semantic versioning](https://semver.org/).  The version
must match in both `MODULE.bazel` and `CMakeLists.txt`.  The release script
validates this.

Tags use the `v` prefix (`v1.0.0`).  This is the convention expected by
the [publish-to-bcr](https://github.com/bazel-contrib/publish-to-bcr)
GitHub App for future Bazel Central Registry publication.

## Downstream usage

### Bazel (pre-built binary)

```python
http_archive(
    name = "kepler-formal-bin",
    url = "https://github.com/keplertech/kepler-formal/releases/download/v1.0.0/kepler-formal-1.0.0-linux-x86_64.tar.gz",
    sha256 = "...",
    build_file_content = 'exports_files(["kepler-formal"])',
)
```

### Bazel (build from source)

```python
bazel_dep(name = "kepler-formal", version = "1.0.0")

git_override(
    module_name = "kepler-formal",
    remote = "https://github.com/keplertech/kepler-formal.git",
    commit = "<sha>",
)
```

## Troubleshooting

**"Version mismatch" error**: Edit both MODULE.bazel and CMakeLists.txt
to have the same version string, commit, and retry.

**"Tag already exists" error**: The version has already been released.
Bump the version to a new number and retry.

**"Working tree is dirty" error**: Commit or stash all changes first.

**Release workflow fails**: Check the Actions tab for logs.  The most
common issue is a stale dependency hash in `bazel/deps.bzl` — update
per the instructions in `docs/bcr-roadmap.md`.
