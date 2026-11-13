#!/usr/bin/env bash
# Sequential equivalence check of an upstream Ibex instruction-cache change
# (commit cb361b2b against its parent, cache configuration) that is a known
# inequivalence: a line filled into two ways makes the cache return different
# data. The difference only appears after the cache's start-up sequence, so it
# exercises scoped SEC with the simulation-seeded search. The check passes when
# kepler-formal reports a counterexample (exit code 3) with an output mismatch.
#
# Usage: ci/ibex_icache_sec.sh <kepler-formal> <work-dir> [kissat|cadical]
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "usage: $0 <kepler-formal> <work-dir> [kissat|cadical]" >&2
  exit 2
fi
kf=$1
solver=${3:-kissat}
mkdir -p "$2"
work=$(cd "$2" && pwd)

ibex_url=https://github.com/keplertech/ibex
flist_rev=58a1333caf0eeed475c767b4261c19a6d4f2106c  # keplertech/ file lists
before_rev=65368a4c41d8e6c2e56ce816e418e21fc4e7b328
after_rev=cb361b2b6cc9ace09eab93efaf9fd7ea3931ab02

repo=$work/ibex
if [[ ! -d $repo/.git ]]; then
  git init -q "$repo"
  git -C "$repo" remote add origin "$ibex_url"
fi
git -C "$repo" fetch -q --depth 1 origin "$flist_rev" "$before_rev" "$after_rev"

# Design sources from each revision, file lists from the keplertech/ tree.
for side in "before:$before_rev" "after:$after_rev"; do
  name=${side%%:*}
  rev=${side#*:}
  rm -rf "${work:?}/$name"
  mkdir -p "$work/$name/keplertech"
  git -C "$repo" archive "$rev" rtl vendor | tar -x -C "$work/$name"
  for flist in ibex_top.f ibex_top_bisect_icache.f; do
    git -C "$repo" show "$flist_rev:keplertech/$flist" > "$work/$name/keplertech/$flist"
  done
  echo "-F $work/$name/keplertech/ibex_top_bisect_icache.f" > "$work/$name.f"
done

cat > "$work/config.yaml" <<EOF
format: systemverilog
verification: sec
sv_design1_flist: $work/before.f
sv_design2_flist: $work/after.f
sv_design1_top: ibex_top
sv_design2_top: ibex_top
use_scopes: true
scope_refinement: false
report_skipped_pos: true
solver: $solver
sec_reset:
  cycles: 2
  ports:
    - name: rst_ni
      active_value: 0
sec_seed_cycles: 300
sec_seed_trials: 32
sec_seed_depth: 40
sec_seed_probe_depth: 8
sec_seed_probe_conflicts: 200000
log_file: $work/sec.log
EOF

cd "$work"
set +e
KEPLER_SEC_DIAG=1 "$kf" --config "$work/config.yaml" > "$work/run.log" 2>&1
status=$?
set -e

grep -sE "finished with status|Observed output mismatches at cycle" "$work/sec.log" || true
grep -sE "SEC seeded search: confirmed difference" "$work/run.log" || true
if [[ $status -ne 3 ]] || ! grep -qs "Observed output mismatches at cycle" "$work/sec.log"; then
  echo "expected a counterexample (exit code 3), got exit code $status" >&2
  tail -n 40 "$work/sec.log" >&2 2>/dev/null || true
  exit 1
fi
echo "counterexample found (exit code $status)"
