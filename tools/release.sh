#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only
#
# Create a release tag and push it to trigger the CI release workflow.
# Usage: bazelisk run //:release
set -euo pipefail

# Parse version from MODULE.bazel
BAZEL_VERSION=$(grep -oP 'version\s*=\s*"\K[^"]+' MODULE.bazel)
# Parse version from CMakeLists.txt
CMAKE_VERSION=$(grep -oP 'VERSION\s+\K[0-9]+\.[0-9]+\.[0-9]+' CMakeLists.txt | head -1)

if [ -z "$BAZEL_VERSION" ]; then
  echo "ERROR: Could not parse version from MODULE.bazel"
  exit 1
fi

if [ -z "$CMAKE_VERSION" ]; then
  echo "ERROR: Could not parse version from CMakeLists.txt"
  exit 1
fi

if [ "$BAZEL_VERSION" != "$CMAKE_VERSION" ]; then
  echo "ERROR: Version mismatch: MODULE.bazel=${BAZEL_VERSION} CMakeLists.txt=${CMAKE_VERSION}"
  echo "Update both files to the same version before releasing."
  exit 1
fi

TAG="v${BAZEL_VERSION}"

if [ -n "$(git status --porcelain)" ]; then
  echo "ERROR: Working tree is dirty. Commit or stash changes first."
  exit 1
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
  echo "ERROR: Tag $TAG already exists."
  echo "Bump the version in MODULE.bazel and CMakeLists.txt, commit, then retry."
  exit 1
fi

echo "Creating tag $TAG ..."
git tag -a "$TAG" -m "Release $TAG"

echo "Pushing tag $TAG ..."
git push origin "$TAG"

echo ""
echo "Done. GitHub Actions will build and publish the release at:"
echo "  https://github.com/keplertech/kepler-formal/releases/tag/$TAG"
