#!/usr/bin/env bash
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set -euo pipefail

helper="$1"
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/kepler-sec-regress-test.XXXXXX")"
trap 'rm -rf "${tmp_dir}"' EXIT

# Copy the helper so its output directory also stays inside the temporary tree.
mkdir -p "${tmp_dir}/repo/regress" "${tmp_dir}/case"
cp "${helper}" "${tmp_dir}/repo/regress/run_sec_strategies_regress.sh"
printf 'verification: sec\n' > "${tmp_dir}/case/config.yaml"

cat > "${tmp_dir}/fake-kepler-formal" <<'EOF'
#!/usr/bin/env bash
echo "SEC partially proved equivalence at k = 1: 1/2 outputs proved; remaining outputs are inconclusive."
exit 1
EOF
chmod +x "${tmp_dir}/fake-kepler-formal"

cat > "${tmp_dir}/fake-equivalent-kepler-formal" <<'EOF'
#!/usr/bin/env bash
echo "No difference was found. SEC proved equivalence at k = 1."
exit 0
EOF
chmod +x "${tmp_dir}/fake-equivalent-kepler-formal"

cat > "${tmp_dir}/fake-inconclusive-kepler-formal" <<'EOF'
#!/usr/bin/env bash
echo "SEC was inconclusive up to max_k = 1: no proof or counterexample"
exit 2
EOF
chmod +x "${tmp_dir}/fake-inconclusive-kepler-formal"

cat > "${tmp_dir}/fake-different-kepler-formal" <<'EOF'
#!/usr/bin/env bash
echo "Difference was found. SEC found a counterexample at k = 1."
exit 3
EOF
chmod +x "${tmp_dir}/fake-different-kepler-formal"

for expectation in expect-equivalent-or-partial allow-inconclusive allow-unset-state-inconclusive; do
  bash "${tmp_dir}/repo/regress/run_sec_strategies_regress.sh" \
    "partial-${expectation}" \
    "${tmp_dir}/case" \
    "${tmp_dir}/fake-kepler-formal" \
    config.yaml \
    "${expectation}" \
    engine=pdr
done

bash "${tmp_dir}/repo/regress/run_sec_strategies_regress.sh" \
  equivalent-positive \
  "${tmp_dir}/case" \
  "${tmp_dir}/fake-equivalent-kepler-formal" \
  config.yaml \
  expect-equivalent-or-partial \
  engine=pdr

bash "${tmp_dir}/repo/regress/run_sec_strategies_regress.sh" \
  inconclusive-measurement \
  "${tmp_dir}/case" \
  "${tmp_dir}/fake-inconclusive-kepler-formal" \
  config.yaml \
  allow-inconclusive \
  engine=pdr

bash "${tmp_dir}/repo/regress/run_sec_strategies_regress.sh" \
  different-negative \
  "${tmp_dir}/case" \
  "${tmp_dir}/fake-different-kepler-formal" \
  config.yaml \
  expect-different \
  engine=pdr

if bash "${tmp_dir}/repo/regress/run_sec_strategies_regress.sh" \
    inconclusive-positive \
    "${tmp_dir}/case" \
    "${tmp_dir}/fake-inconclusive-kepler-formal" \
    config.yaml \
    expect-equivalent-or-partial \
    engine=pdr; then
  echo "Positive equivalence unexpectedly accepted an inconclusive result" >&2
  exit 1
fi

# Strict equivalence remains available for regressions that require a full
# proof; it must continue to reject the distinct partial-proof exit code.
if bash "${tmp_dir}/repo/regress/run_sec_strategies_regress.sh" \
    partial-strict \
    "${tmp_dir}/case" \
    "${tmp_dir}/fake-kepler-formal" \
    config.yaml \
    expect-equivalent \
    engine=pdr; then
  echo "Strict equivalence unexpectedly accepted a partial proof" >&2
  exit 1
fi
