#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only
set -euo pipefail

binary="$1"
# The command must work with no designs, config, or Python module search path.
output=$(PYTHONPATH=/nonexistent "$binary" --version)
short_output=$(PYTHONPATH=/nonexistent "$binary" -V)
[ "$output" = "$short_output" ]
[ "$(printf '%s\n' "$output" | wc -l | tr -d ' ')" = 4 ]
for project in kepler-formal naja; do
  printf '%s\n' "$output" | grep -Eq "^${project} version: [0-9]+\.[0-9]+\.[0-9]+$"
  printf '%s\n' "$output" | grep -Eq "^${project} git hash: ([0-9a-f]{7,40}|unknown)$"
done

# Informational options must not accidentally launch a verification workflow.
if "$binary" --version --config /nonexistent >/dev/null 2>&1; then
  echo 'Version mode unexpectedly accepted a config' >&2
  exit 1
fi
