#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0
#
# Install kepler-formal into <prefix>/bin, its licenses into
# <prefix>/share/kepler-formal.
# Usage: bazelisk run -c opt //:install -- <prefix>

# --- begin runfiles.bash initialization v3 ---
set -uo pipefail; set +e; f=bazel_tools/tools/bash/runfiles/runfiles.bash
# shellcheck disable=SC1090
source "${RUNFILES_DIR:-/dev/null}/$f" 2>/dev/null || \
  source "$(grep -sm1 "^$f " "${RUNFILES_MANIFEST_FILE:-/dev/null}" | cut -f2- -d' ')" 2>/dev/null || \
  source "$0.runfiles/$f" 2>/dev/null || \
  source "$(grep -sm1 "^$f " "$0.runfiles_manifest" | cut -f2- -d' ')" 2>/dev/null || \
  source "$(grep -sm1 "^$f " "$0.exe.runfiles_manifest" | cut -f2- -d' ')" 2>/dev/null || \
  { echo>&2 "ERROR: cannot find $f"; exit 1; }; f=; set -e
# --- end runfiles.bash initialization v3 ---

if [[ $# -ne 1 ]]; then
  echo "Usage: bazelisk run -c opt //:install -- <prefix>" >&2
  exit 2
fi

# A relative prefix is relative to where bazelisk was run, not to the
# runfiles directory this script runs in.
prefix="$1"
if [[ "${prefix}" != /* ]]; then
  prefix="${BUILD_WORKING_DIRECTORY}/${prefix}"
fi

tarfile="$(rlocation kepler-formal/packaging/kepler-formal.tar)"

# Bazel outputs are read-only; make a previous install removable.
for path in bin/kepler-formal bin/kepler-formal.runfiles \
  bin/kepler-formal.repo_mapping bin/naja.so share/kepler-formal; do
  if [[ -e "${prefix}/${path}" ]]; then
    chmod -R u+w "${prefix}/${path}"
    rm -rf "${prefix:?}/${path}"
  fi
done
mkdir -p "${prefix}"
tar -xf "${tarfile}" -C "${prefix}"
chmod -R u+w "${prefix}/bin" "${prefix}/share/kepler-formal"

echo "kepler-formal installed to $(cd "${prefix}" && pwd)/bin/kepler-formal"
