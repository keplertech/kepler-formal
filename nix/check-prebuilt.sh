#!/usr/bin/env bash
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <published-flake-reference>" >&2
  exit 1
fi
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cache=https://keplertech.cachix.org
no_build=(--max-jobs 0 --option builders '' --option allow-import-from-derivation false)

# Resolve once, so the cache check and installation cannot use different
# revisions of a moving branch. Evaluation does not build the package.
locked_flake=$(nix flake metadata "${no_build[@]}" --json --no-update-lock-file "$1" | jq -er '.url')
system=$(nix eval "${no_build[@]}" --impure --raw --expr builtins.currentSystem)
package=$(nix eval "${no_build[@]}" --raw --no-update-lock-file \
  "$locked_flake#packages.$system.default.outPath")
echo "Testing prebuilt package: $package ($system)"
echo "Locked flake: $locked_flake"

# Fail explicitly on an unpublished revision or a reused store. Otherwise a
# pre-existing local build could make this test pass without any cache download.
if ! nix path-info --store "$cache" --option narinfo-cache-negative-ttl 0 "$package"; then
  echo "Package is not available in the public keplertech cache; publish this revision first." >&2
  exit 1
fi
if nix path-info "$package" >/dev/null 2>&1; then
  echo "Package is already in the local store; this test requires a fresh runner." >&2
  exit 1
fi

install_root=$(mktemp -d)
trap 'rm -rf "$install_root"' EXIT
profile="$install_root/profile"

# Use the user-facing profile installation, not nix build/nix copy. Disable
# local AND remote builders, and allow only the public package/dependency caches.
nix profile install --profile "$profile" --no-update-lock-file \
  "${no_build[@]}" \
  --option substituters "$cache https://cache.nixos.org" \
  --option extra-substituters '' --option fallback false \
  --option narinfo-cache-negative-ttl 0 "$locked_flake"

nix profile list --profile "$profile" --json | \
  jq -e --arg package "$package" '[.elements[].storePaths[]] == [$package]'
test -x "$profile/bin/kepler-formal"
cd "$install_root"
"$profile/bin/kepler-formal" --help >/dev/null
bash "$script_dir/check.sh" "$package"
echo "Prebuilt profile installation passed without building."
