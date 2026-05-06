#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 <test-name> <case-dir> <kepler-formal-bin> <config-path> [expect-identical]" >&2
  exit 2
fi

test_name="$1"
case_dir="$2"
kepler_formal_bin="$3"
config_path="$4"
expect_identical="${5:-}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

output_dir="${repo_root}/regress-output/${test_name}"
output_cnf="${output_dir}/miter.cnf"
output_log="${output_dir}/miter.log"
tmp_config="${output_dir}/config.with_cnf.yaml"
golden_cnf="${repo_root}/regress/${test_name}/miter.cnf"

mkdir -p "${output_dir}"

(
  cd "${case_dir}"
  awk '
    /^[[:space:]]*cnf_export:/ { next }
    /^[[:space:]]*cnf_export_path:/ { next }
    /^[[:space:]]*log_file:/ { next }
    { print }
  ' "${config_path}" > "${tmp_config}"
  {
    echo
    echo "log_file: ${output_log}"
    echo "cnf_export: true"
    echo "cnf_export_path: ${output_cnf}"
  } >> "${tmp_config}"

  "${kepler_formal_bin}" --config "${tmp_config}"

  if [[ -n "${expect_identical}" ]]; then
    grep "Circuits are IDENTICAL" "${output_log}"
  fi
)

if [[ ! -s "${output_cnf}" ]]; then
  echo "CNF dump was not produced: ${output_cnf}" >&2
  exit 1
fi

diff -u "${golden_cnf}" "${output_cnf}"
