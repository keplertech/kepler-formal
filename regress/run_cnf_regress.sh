#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <test-name> <case-dir> <kepler-formal-bin> <config-path> [expect-identical] [po-cnfs]" >&2
  exit 2
fi

test_name="$1"
case_dir="$2"
kepler_formal_bin="$3"
config_path="$4"
expect_identical=""
dump_po_cnfs=""

for option in "${@:5}"; do
  case "${option}" in
    expect-identical)
      expect_identical="1"
      ;;
    po-cnfs)
      dump_po_cnfs="1"
      ;;
    *)
      echo "Unknown option: ${option}" >&2
      exit 2
      ;;
  esac
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

output_dir="${repo_root}/regress-output/${test_name}"
output_cnf="${output_dir}/miter.cnf"
output_cnf_archive="${output_dir}/miter.cnf.tar.gz"
output_po_cnf_dir="${output_dir}/po_cnfs"
output_po_cnf_archive="${output_dir}/po_cnfs.tar.gz"
output_log="${output_dir}/miter.log"
tmp_config="${output_dir}/config.with_cnf.yaml"
golden_cnf_archive="${repo_root}/regress/${test_name}/miter.cnf.tar.gz"
golden_po_cnf_archive="${repo_root}/regress/${test_name}/po_cnfs.tar.gz"

create_deterministic_tar_gz() {
  local source="$1"
  local archive_root="$2"
  local destination="$3"

  python3 - "${source}" "${archive_root}" "${destination}" <<'PY'
import gzip
import os
import sys
import tarfile

source, archive_root, destination = sys.argv[1:]

def reset_metadata(info):
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = 0
    return info

def add_file(archive, path, arcname):
    info = reset_metadata(archive.gettarinfo(path, arcname=arcname))
    with open(path, "rb") as handle:
        archive.addfile(info, handle)

with open(destination, "wb") as raw:
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode="w") as archive:
            if os.path.isdir(source):
                for root, dirs, files in os.walk(source):
                    dirs.sort()
                    for name in sorted(files):
                        path = os.path.join(root, name)
                        arcname = os.path.join(archive_root, os.path.relpath(path, source))
                        add_file(archive, path, arcname)
            else:
                add_file(archive, source, archive_root)
PY
}

mkdir -p "${output_dir}"
rm -rf "${output_po_cnf_dir}"
rm -f "${output_cnf}" "${output_cnf_archive}" "${output_po_cnf_archive}"

(
  cd "${case_dir}"
  awk '
    /^[[:space:]]*cnf_export:/ { next }
    /^[[:space:]]*cnf_export_path:/ { next }
    /^[[:space:]]*po_cnf_export:/ { next }
    /^[[:space:]]*po_cnf_export_path:/ { next }
    /^[[:space:]]*log_file:/ { next }
    { print }
  ' "${config_path}" > "${tmp_config}"
  {
    echo
    echo "log_file: ${output_log}"
    echo "cnf_export: true"
    echo "cnf_export_path: ${output_cnf}"
    if [[ -n "${dump_po_cnfs}" ]]; then
      echo "po_cnf_export: true"
      echo "po_cnf_export_path: ${output_po_cnf_dir}"
    fi
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

create_deterministic_tar_gz "${output_cnf}" "miter.cnf" "${output_cnf_archive}"
cmp "${golden_cnf_archive}" "${output_cnf_archive}"

if [[ -n "${dump_po_cnfs}" ]]; then
  if [[ ! -d "${output_po_cnf_dir}" ]]; then
    echo "PO CNF dump directory was not produced: ${output_po_cnf_dir}" >&2
    exit 1
  fi

  create_deterministic_tar_gz "${output_po_cnf_dir}" "po_cnfs" "${output_po_cnf_archive}"
  cmp "${golden_po_cnf_archive}" "${output_po_cnf_archive}"
fi
