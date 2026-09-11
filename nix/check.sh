#!/usr/bin/env bash
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <installed-package>" >&2
  exit 1
fi
package=$(cd "$1" && pwd -P)
test -x "$package/bin/kepler-formal"
test -f "$package/bin/naja.so"
workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT
cd "$workdir"

# Exercise the installed executable from outside the source/build directories.
run_case() {
  local name=$1 expected_status=$2 expected_message=$3 status=0
  env -u PYTHONPATH -u PYTHONHOME -u LD_LIBRARY_PATH -u DYLD_LIBRARY_PATH \
    PYTHONNOUSERSITE=1 KEPLER_SMOKE_NAJA="$package/bin/naja.so" \
    "$package/bin/kepler-formal" --config "$name.yaml" > "$name.output" 2>&1 \
    || status=$?
  if [[ $status -ne $expected_status ]] ||
      ! grep -Fq "$expected_message" "$name.output"; then
    cat "$name.output" >&2
    echo "$name: expected exit $expected_status, got $status" >&2
    exit 1
  fi
}

write_register() {
  cat > "$1.sv" <<EOF
module top(input logic clk, reset, d, output logic q);
  always_ff @(posedge clk) q <= reset ? 1'b0 : $2;
endmodule
EOF
}
write_register reference d
write_register equivalent "(d ^ 1'b0)"
write_register different "~d"
for name in equivalent different; do
  cat > "$name.yaml" <<EOF
format: systemverilog
verification: sec
sec_engine: pdr
sec_encoding: dual_rail_steady
max_k: 4
sv_design1_top: top
sv_design2_top: top
input_paths:
  - reference.sv
  - $name.sv
log_file: $name.log
EOF
done
run_case equivalent 0 "SEC proved equivalence"
run_case different 3 "Difference was found. SEC found a counterexample"
for name in equivalent different; do
  grep -Fq "SEC checked-output coverage: 100.00% (1/1 covered/existing outputs)." \
    "$name.output" || { cat "$name.output" >&2; exit 1; }
done

# A real Python-defined cell checks the adjacent extension and embedded Python.
cat > primitives.py <<'EOF'
import os
import naja

assert os.path.realpath(naja.__file__) == os.path.realpath(os.environ['KEPLER_SMOKE_NAJA'])

def constructPrimitives(lib):
    cell = naja.SNLDesign.createPrimitive(lib, 'BUF')
    naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Input, 'A')
    naja.SNLScalarTerm.create(cell, naja.SNLTerm.Direction.Output, 'Z')
    cell.setTruthTable(0b10)
EOF
cat > mapped.v <<'EOF'
module top(input a, output y);
  BUF u_buf(.A(a), .Z(y));
endmodule
EOF
cat > direct.v <<'EOF'
module top(input a, output y);
  assign y = a;
endmodule
EOF
cat > python-tech.yaml <<'EOF'
format: verilog
verification: lec
input_paths:
  - mapped.v
  - direct.v
py_tech_files:
  - primitives.py
log_file: python-tech.log
EOF
run_case python-tech 0 "No difference was found."
echo "Installed package smoke checks passed (SEC equivalence, SEC mismatch, Python technology)."
