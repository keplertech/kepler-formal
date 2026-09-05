#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <installed-runtime-directory>" >&2
  exit 2
fi

stage="$(readlink -f "$1")"
lib_dir="${stage}/lib"

if [[ ! -x "${stage}/bin/kepler-formal" ]]; then
  echo "Missing runtime executable: ${stage}/bin/kepler-formal" >&2
  exit 2
fi

mkdir -p "${lib_dir}"

is_host_runtime_library() {
  case "$1" in
    linux-vdso.so.*|ld-linux-*.so.*|libanl.so.*|libBrokenLocale.so.*|\
    libc.so.*|libdl.so.*|libgcc_s.so.*|libm.so.*|libmvec.so.*|\
    libnss_*.so.*|libpthread.so.*|libresolv.so.*|librt.so.*|\
    libstdc++.so.*|libthread_db.so.*|libutil.so.*)
      return 0
      ;;
  esac
  return 1
}

declare -a queue=()
declare -a roots=()
declare -A scanned=()
dependency_pattern='^[[:space:]]*([^[:space:]]+)[[:space:]]+=>[[:space:]]+([^[:space:]]+)'

while IFS= read -r -d '' candidate; do
  roots+=("${candidate}")
  queue+=("${candidate}")
done < <(
  find -L "${stage}/bin" "${stage}/lib" \
    -type f \( -name kepler-formal -o -name '*.so' -o -name '*.so.*' \) \
    -print0
)

copied=0
for ((index = 0; index < ${#queue[@]}; ++index)); do
  target="${queue[index]}"
  canonical="$(readlink -f "${target}")"
  if [[ -n "${scanned[${canonical}]:-}" ]]; then
    continue
  fi
  scanned["${canonical}"]=1

  if ! dependencies="$(LD_LIBRARY_PATH="${lib_dir}:${LD_LIBRARY_PATH:-}" ldd "${target}" 2>&1)"; then
    if [[ "${dependencies}" == *"not a dynamic executable"* ]]; then
      continue
    fi
    echo "Failed to inspect runtime dependencies for ${target}:" >&2
    echo "${dependencies}" >&2
    exit 1
  fi

  while IFS= read -r line; do
    if [[ "${line}" != *"=>"* ]]; then
      continue
    fi
    if [[ "${line}" == *"=> not found"* ]]; then
      echo "Unresolved runtime dependency for ${target}: ${line}" >&2
      exit 1
    fi
    if [[ ! "${line}" =~ ${dependency_pattern} ]]; then
      continue
    fi

    soname="${BASH_REMATCH[1]}"
    resolved="${BASH_REMATCH[2]}"
    if is_host_runtime_library "${soname}"; then
      continue
    fi

    resolved="$(readlink -f "${resolved}")"
    if [[ "${resolved}" == "${stage}/"* ]]; then
      queue+=("${resolved}")
      continue
    fi

    destination="${lib_dir}/${soname}"
    if [[ ! -e "${destination}" && ! -L "${destination}" ]]; then
      cp -L -- "${resolved}" "${destination}"
      chmod u+w "${destination}"
      ((copied += 1))
    fi
    queue+=("${destination}")
  done <<< "${dependencies}"
done

for root in "${roots[@]}"; do
  if unresolved="$(LD_LIBRARY_PATH="${lib_dir}:${LD_LIBRARY_PATH:-}" ldd "${root}" 2>&1 | grep '=> not found')"; then
    echo "Runtime artifact still has unresolved dependencies for ${root}:" >&2
    echo "${unresolved}" >&2
    exit 1
  fi
done

echo "Bundled ${copied} non-system runtime libraries in ${lib_dir}."
