#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <test-name> <case-dir> <kepler-formal-bin> <config-path> [expect-equivalent|expect-equivalent-or-partial|expect-different|expect-unsupported|expect-full-coverage|allow-inconclusive|allow-unset-state-inconclusive] [max-k=<n>] [compact] [engine=<name>] [sec-encoding=<name>]" >&2
  exit 2
fi

test_name="$1"
case_dir="$2"
kepler_formal_bin="$3"
config_path="$4"
expectation=""
max_k_override=""
compact_mode=""
# The CLI default is dual-rail SEC.  This regression helper keeps the
# historical binary workflow behavior unless a caller explicitly opts into
# dual_rail_steady, which prevents old regressions from silently changing mode.
sec_encoding="${SEC_ENCODING:-binary}"
# By default the helper is useful for local all-engine smoke checks.  CI passes
# engine=<name> from the split regress workflows so each job owns one strategy.
engines=(k_induction imc pdr)

for option in "${@:5}"; do
  case "${option}" in
    expect-equivalent|expect-equivalent-or-partial|expect-different|expect-unsupported|expect-full-coverage|allow-inconclusive|allow-unset-state-inconclusive)
      expectation="${option}"
      ;;
    compact)
      compact_mode="1"
      ;;
    engine=*)
      engine="${option#engine=}"
      case "${engine}" in
        k_induction|imc|pdr)
          engines=("${engine}")
          ;;
        *)
          echo "Invalid SEC engine override: ${engine}" >&2
          exit 2
          ;;
      esac
      ;;
    sec-encoding=*|encoding=*)
      sec_encoding="${option#*=}"
      ;;
    max-k=*)
      max_k_override="${option#max-k=}"
      if [[ ! "${max_k_override}" =~ ^[0-9]+$ ]]; then
        echo "Invalid max-k override: ${max_k_override}" >&2
        exit 2
      fi
      ;;
    *)
      echo "Unknown option: ${option}" >&2
      exit 2
      ;;
  esac
done

case "${sec_encoding}" in
  ""|binary|dual_rail_steady)
    ;;
  *)
    echo "Invalid SEC encoding override: ${sec_encoding}" >&2
    exit 2
    ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
output_dir="${repo_root}/regress-output/${test_name}/sec"

mkdir -p "${output_dir}"

format_bytes() {
  local bytes="$1"
  if [[ -z "${bytes}" || "${bytes}" == "max" ]]; then
    printf '%s' "${bytes}"
    return
  fi
  awk -v bytes="${bytes}" '
    BEGIN {
      value = bytes + 0
      split("B KiB MiB GiB TiB", units, " ")
      unit = 1
      while (value >= 1024 && unit < 5) {
        value /= 1024
        unit += 1
      }
      if (unit == 1) {
        printf "%d%s", value, units[unit]
      } else {
        printf "%.1f%s", value, units[unit]
      }
    }'
}

print_regress_memory_snapshot() {
  local stage="$1"
  local engine="$2"
  local status="${3:-running}"
  local prefix="[sec-regress-memory] stage=${stage} test=${test_name} engine=${engine} status=${status}"

  if [[ -r /proc/meminfo ]]; then
    awk -v prefix="${prefix}" '
      function fmt(kib, value, unit, units) {
        split("KiB MiB GiB TiB", units, " ")
        value = kib
        unit = 1
        while (value >= 1024 && unit < 4) {
          value /= 1024
          unit += 1
        }
        return sprintf("%.1f%s", value, units[unit])
      }
      /^(MemTotal|MemAvailable|MemFree|Buffers|Cached|SwapTotal|SwapFree):/ {
        key = $1
        sub(":", "", key)
        value[key] = $2 + 0
      }
      END {
        if ("MemAvailable" in value) {
          printf "%s meminfo MemAvailable=%s MemFree=%s Buffers=%s Cached=%s SwapFree=%s SwapTotal=%s MemTotal=%s\n",
              prefix,
              fmt(value["MemAvailable"]),
              fmt(value["MemFree"]),
              fmt(value["Buffers"]),
              fmt(value["Cached"]),
              fmt(value["SwapFree"]),
              fmt(value["SwapTotal"]),
              fmt(value["MemTotal"])
        }
      }' /proc/meminfo || true
  elif command -v vm_stat >/dev/null 2>&1; then
    vm_stat | awk -v prefix="${prefix}" '
      /page size of/ {
        page_size = $8
        gsub("\\.", "", page_size)
      }
      /Pages free:/ {
        free = $3
        gsub("\\.", "", free)
      }
      /Pages inactive:/ {
        inactive = $3
        gsub("\\.", "", inactive)
      }
      /Pages speculative:/ {
        speculative = $3
        gsub("\\.", "", speculative)
      }
      END {
        if (page_size == "") {
          page_size = 4096
        }
        available = (free + inactive + speculative) * page_size
        printf "%s vm_stat estimated_available=%s pages_free=%d pages_inactive=%d pages_speculative=%d page_size=%d\n",
            prefix,
            fmt(available),
            free,
            inactive,
            speculative,
            page_size
      }
      function fmt(bytes, value, unit, units) {
        split("B KiB MiB GiB TiB", units, " ")
        value = bytes
        unit = 1
        while (value >= 1024 && unit < 5) {
          value /= 1024
          unit += 1
        }
        return sprintf("%.1f%s", value, units[unit])
      }' || true
  fi

  if [[ -r /sys/fs/cgroup/memory.current ]]; then
    local current peak max swap_current swap_peak
    current="$(< /sys/fs/cgroup/memory.current)"
    peak="$(< /sys/fs/cgroup/memory.peak)"
    max="$(< /sys/fs/cgroup/memory.max)"
    swap_current=""
    swap_peak=""
    if [[ -r /sys/fs/cgroup/memory.swap.current ]]; then
      swap_current="$(< /sys/fs/cgroup/memory.swap.current)"
    fi
    if [[ -r /sys/fs/cgroup/memory.swap.peak ]]; then
      swap_peak="$(< /sys/fs/cgroup/memory.swap.peak)"
    fi
    printf '%s cgroup_v2 memory.current=%s memory.peak=%s memory.max=%s swap.current=%s swap.peak=%s\n' \
      "${prefix}" \
      "$(format_bytes "${current}")" \
      "$(format_bytes "${peak}")" \
      "$(format_bytes "${max}")" \
      "$(format_bytes "${swap_current}")" \
      "$(format_bytes "${swap_peak}")"
  elif [[ -r /sys/fs/cgroup/memory/memory.usage_in_bytes ]]; then
    local current peak limit
    current="$(< /sys/fs/cgroup/memory/memory.usage_in_bytes)"
    peak="$(< /sys/fs/cgroup/memory/memory.max_usage_in_bytes)"
    limit="$(< /sys/fs/cgroup/memory/memory.limit_in_bytes)"
    printf '%s cgroup_v1 memory.usage=%s memory.max_usage=%s memory.limit=%s\n' \
      "${prefix}" \
      "$(format_bytes "${current}")" \
      "$(format_bytes "${peak}")" \
      "$(format_bytes "${limit}")"
  fi

  if [[ -r /proc/buddyinfo ]]; then
    local page_kib=4
    if command -v getconf >/dev/null 2>&1; then
      local page_size
      page_size="$(getconf PAGESIZE 2>/dev/null || true)"
      if [[ "${page_size}" =~ ^[0-9]+$ && "${page_size}" -gt 0 ]]; then
        page_kib=$((page_size / 1024))
      fi
    fi
    awk -v prefix="${prefix}" -v page_kib="${page_kib}" '
      {
        for (i = 5; i <= NF; ++i) {
          count = $i + 0
          order = i - 5
          if (count > 0 && order > largest) {
            largest = order
          }
          free_pages += count * (2 ^ order)
        }
      }
      END {
        if (largest >= 0) {
          printf "%s buddyinfo largest_free_order=%d largest_free_block=%s estimated_free=%s\n",
              prefix,
              largest,
              fmt((2 ^ largest) * page_kib),
              fmt(free_pages * page_kib)
        }
      }
      function fmt(kib, value, unit, units) {
        split("KiB MiB GiB TiB", units, " ")
        value = kib
        unit = 1
        while (value >= 1024 && unit < 4) {
          value /= 1024
          unit += 1
        }
        return sprintf("%.1f%s", value, units[unit])
      }' /proc/buddyinfo || true
  fi
}

run_engine() {
  local engine="$1"
  local tmp_config="${output_dir}/config.${engine}.yaml"
  local stdout_log="${output_dir}/${engine}.stdout"
  local memory_snapshot_recorded=0

  (
    report_memory_on_exit() {
      local status=$?
      if [[ "${memory_snapshot_recorded}" -eq 0 ]]; then
        print_regress_memory_snapshot "exit" "${engine}" "${status}"
      fi
    }
    trap report_memory_on_exit EXIT

    cd "${case_dir}"
    # SEC currently rejects CNF-export options. Keep the design, library, and
    # solver settings from the original regression config, then override only
    # the verification mode and selected SEC strategy. Drop LEC max_k by
    # default: PDR proves by frame convergence and the split SEC regressions
    # should not accidentally inherit toy bounded-check depths such as k=4.
    # Do not force a log_file here: repeated local regress runs are easier to
    # inspect when kepler-formal keeps its own per-run log naming.
    awk -v max_k_override="${max_k_override}" -v compact_mode="${compact_mode}" '
      /^[[:space:]]*verification:/ { next }
      /^[[:space:]]*sec_engine:/ { next }
      /^[[:space:]]*sec_encoding:/ { next }
      /^[[:space:]]*max_k:/ { next }
      /^[[:space:]]*compact_mode:/ {
        if (compact_mode != "") {
          next
        }
      }
      /^[[:space:]]*cnf_export:/ { next }
      /^[[:space:]]*cnf_export_path:/ { next }
      /^[[:space:]]*po_cnf_export:/ { next }
      /^[[:space:]]*po_cnf_export_path:/ { next }
      /^[[:space:]]*log_file:/ { next }
      { print }
    ' "${config_path}" > "${tmp_config}"
    {
      echo
      echo "verification: sec"
      echo "sec_engine: ${engine}"
      if [[ -n "${sec_encoding}" ]]; then
        echo "sec_encoding: ${sec_encoding}"
      fi
      if [[ -n "${max_k_override}" ]]; then
        echo "max_k: ${max_k_override}"
      fi
      if [[ -n "${compact_mode}" ]]; then
        echo "compact_mode: true"
      fi
    } >> "${tmp_config}"

    echo "Running SEC ${engine} for ${test_name}"
    print_regress_memory_snapshot "before" "${engine}" "starting"
    : > "${stdout_log}"
    # Stream the tool output as it is produced instead of only dumping it after
    # completion.  Large SEC/PDR cases can run for minutes between solver
    # decisions, so emit a lightweight heartbeat to keep GitHub logs obviously
    # alive and to make a true hang easier to distinguish from solver work.
    "${kepler_formal_bin}" --config "${tmp_config}" > "${stdout_log}" 2>&1 &
    local kepler_pid=$!
    tail -n +1 -f "${stdout_log}" &
    local tail_pid=$!
    local kepler_status=0
    local heartbeat_elapsed=0
    while kill -0 "${kepler_pid}" 2>/dev/null; do
      sleep 5
      heartbeat_elapsed=$((heartbeat_elapsed + 5))
      if [[ "${heartbeat_elapsed}" -ge 60 ]] &&
          kill -0 "${kepler_pid}" 2>/dev/null; then
        printf '[%s] SEC %s for %s is still running...\n' \
          "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "${engine}" "${test_name}"
        heartbeat_elapsed=0
      fi
    done
    wait "${kepler_pid}" || kepler_status=$?
    sleep 1
    kill "${tail_pid}" 2>/dev/null || true
    wait "${tail_pid}" 2>/dev/null || true
    print_regress_memory_snapshot "after" "${engine}" "${kepler_status}"
    memory_snapshot_recorded=1
    if [[ "${expectation}" == "expect-different" ]]; then
      if [[ "${kepler_status}" -eq 3 ]] &&
          grep -q "SEC found a counterexample" "${stdout_log}"; then
        grep "SEC found a counterexample" "${stdout_log}"
        return 0
      fi
      if [[ "${kepler_status}" -ne 0 ]]; then
        return "${kepler_status}"
      fi
      # Negative SEC regressions must be witnessed by the selected SEC engine.
      # A partial-coverage proof or an LEC-only structural mismatch is not
      # enough for these workflows.
      echo "Expected SEC counterexample for ${test_name} (${engine})" >&2
      return 1
    fi

    # Treat any discovered counterexample as a hard failure unless the caller
    # explicitly asked for one.  This keeps coverage-only checks from masking a
    # real SEC mismatch.
    if grep -q "SEC found a counterexample" "${stdout_log}"; then
      echo "Unexpected SEC counterexample for ${test_name} (${engine})" >&2
      return 1
    fi

    if [[ "${expectation}" == "expect-unsupported" ]]; then
      if [[ "${kepler_status}" -eq 2 ]] &&
          grep -q "SEC cannot run on this design pair" "${stdout_log}"; then
        grep "SEC cannot run on this design pair" "${stdout_log}"
        return 0
      fi
      if [[ "${kepler_status}" -ne 0 ]]; then
        return "${kepler_status}"
      fi
      echo "Expected unsupported SEC result for ${test_name} (${engine})" >&2
      return 1
    fi

    # A partial proof is inconclusive for its remaining outputs and deliberately
    # exits with status 1. Positive regressions may explicitly accept that
    # distinct verdict without accepting a fully inconclusive result.
    if [[ "${kepler_status}" -eq 1 ]] &&
       [[ "${expectation}" == "expect-equivalent-or-partial" ||
          "${expectation}" == "allow-inconclusive" ||
          "${expectation}" == "allow-unset-state-inconclusive" ]] &&
        grep -q "SEC partially proved equivalence" "${stdout_log}"; then
      grep "SEC partially proved equivalence" "${stdout_log}"
      return 0
    fi

    # Measurement-only SEC runs still fail on real counterexamples above, but
    # allow inconclusive positive proofs so one hard design does not stop the
    # rest of the regression from reporting its current behavior.
    if [[ "${expectation}" == "allow-inconclusive" ]]; then
      if [[ "${kepler_status}" -eq 0 ]] &&
          grep -q "SEC proved equivalence" "${stdout_log}"; then
        grep "SEC proved equivalence" "${stdout_log}"
        return 0
      fi
      if [[ "${kepler_status}" -eq 2 ]] &&
          grep -q "SEC was inconclusive" "${stdout_log}"; then
        grep "SEC was inconclusive" "${stdout_log}"
        return 0
      fi
      if [[ "${kepler_status}" -ne 0 ]]; then
        return "${kepler_status}"
      fi
      echo "Expected SEC equivalence or inconclusive result for ${test_name} (${engine})" >&2
      return 1
    fi

    # Non-dual positive SEC regressions may have all observed outputs skipped
    # because both sides depend on reset-unanchored internal state. Treat that
    # as measurement-only inconclusive when the workflow explicitly asks for it.
    if [[ "${expectation}" == "allow-unset-state-inconclusive" ]]; then
      if [[ "${kepler_status}" -eq 0 ]] &&
          grep -q "SEC proved equivalence" "${stdout_log}"; then
        grep "SEC proved equivalence" "${stdout_log}"
        return 0
      fi
      if [[ "${kepler_status}" -eq 2 ]] &&
          grep -q "SEC was inconclusive" "${stdout_log}"; then
        grep "SEC was inconclusive" "${stdout_log}"
        return 0
      fi
      if [[ "${kepler_status}" -eq 2 ]] &&
          grep -q "No aligned observed outputs remain after skipping cones that depend on reset-unanchored internal state" "${stdout_log}"; then
        grep "SEC cannot run on this design pair" "${stdout_log}"
        return 0
      fi
      if [[ "${kepler_status}" -ne 0 ]]; then
        return "${kepler_status}"
      fi
      echo "Expected SEC equivalence, inconclusive, or reset-unanchored no-output result for ${test_name} (${engine})" >&2
      return 1
    fi

    if [[ "${expectation}" == "expect-full-coverage" ]]; then
      if [[ "${kepler_status}" -ne 0 ]]; then
        return "${kepler_status}"
      fi
      grep -E "SEC (checked-output|output) coverage: 100.00%" "${stdout_log}"
      grep "SEC proved equivalence" "${stdout_log}"
      return 0
    fi

    if [[ "${kepler_status}" -ne 0 ]]; then
      return "${kepler_status}"
    fi

    if [[ "${expectation}" == "expect-equivalent" ||
          "${expectation}" == "expect-equivalent-or-partial" ]]; then
      grep "SEC proved equivalence" "${stdout_log}"
    else
      grep -E "SEC proved equivalence|SEC found a counterexample" "${stdout_log}"
    fi
  )
}

for engine in "${engines[@]}"; do
  run_engine "${engine}"
done
