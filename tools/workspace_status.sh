#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only
set -euo pipefail

git_hash=unknown
if [ -e .git ]; then
  git_hash=$(git rev-parse HEAD 2>/dev/null) || git_hash=unknown
fi
printf 'STABLE_KEPLER_GIT_HASH %s\n' "$git_hash"
