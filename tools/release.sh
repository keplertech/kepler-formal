#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only
#
# Create a release tag and push it to trigger the CI release workflow.
# Usage: bazelisk run //:release
set -euo pipefail

# Parse version from MODULE.bazel
BAZEL_VERSION=$(
  sed -nE 's/.*version[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/p' MODULE.bazel |
    head -1
)
# Parse version from CMakeLists.txt
CMAKE_VERSION=$(
  sed -nE 's/.*VERSION[[:space:]]+([0-9]+\.[0-9]+\.[0-9]+).*/\1/p' CMakeLists.txt |
    head -1
)
# Parse version from the optional MCP add-on.
MCP_VERSION=$(
  sed -nE 's/^version[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/p' mcp/pyproject.toml |
    head -1
)

if [ -z "$BAZEL_VERSION" ]; then
  echo "ERROR: Could not parse version from MODULE.bazel"
  exit 1
fi

if [ -z "$CMAKE_VERSION" ]; then
  echo "ERROR: Could not parse version from CMakeLists.txt"
  exit 1
fi

if [ -z "$MCP_VERSION" ]; then
  echo "ERROR: Could not parse version from mcp/pyproject.toml"
  exit 1
fi

if [ "$BAZEL_VERSION" != "$CMAKE_VERSION" ] || [ "$BAZEL_VERSION" != "$MCP_VERSION" ]; then
  echo "ERROR: Version mismatch:"
  echo "  MODULE.bazel=${BAZEL_VERSION}"
  echo "  CMakeLists.txt=${CMAKE_VERSION}"
  echo "  mcp/pyproject.toml=${MCP_VERSION}"
  echo "Update all three files to the same version before releasing."
  exit 1
fi

TAG="v${BAZEL_VERSION}"

if [ -n "$(git status --porcelain)" ]; then
  echo "ERROR: Working tree is dirty. Commit or stash changes first."
  exit 1
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
  echo "ERROR: Tag $TAG already exists."
  echo "Bump the version in MODULE.bazel, CMakeLists.txt, and"
  echo "mcp/pyproject.toml, commit, then retry."
  exit 1
fi

echo "Creating tag $TAG ..."
git tag -a "$TAG" -m "Release $TAG"

echo "Pushing tag $TAG ..."
git push origin "$TAG"

echo ""
echo "Done. GitHub Actions will build and publish the release at:"
echo "  https://github.com/keplertech/kepler-formal/releases/tag/$TAG"
