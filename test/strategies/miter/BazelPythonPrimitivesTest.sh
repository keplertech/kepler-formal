#!/usr/bin/env bash

set -eu

workspace="${TEST_SRCDIR}/${TEST_WORKSPACE}"
workdir="${TEST_TMPDIR}/xilinx-python-primitives"
mkdir -p "${workdir}"
cp -R "${workspace}/examples/xilinx" "${workdir}/"

cd "${workdir}/xilinx/register_slice"
"${workspace}/src/bin/kepler-formal" \
  --config test_config_verilog_xilinx_sec.yaml

cd "${workdir}/xilinx/vexriscv"
set +e
difference_output="$("${workspace}/src/bin/kepler-formal" \
  --config test_config_verilog_vexriscv_genfull_xilinx_lec_different.yaml 2>&1)"
difference_status=$?
set -e
if [[ ${difference_status} -ne 0 || "${difference_output}" != *"Difference was found."* ]]; then
  printf '%s\n' "${difference_output}" >&2
  exit 1
fi
