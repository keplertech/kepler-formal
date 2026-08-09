#!/usr/bin/env bash

set -eu

workspace="${TEST_SRCDIR}/${TEST_WORKSPACE}"
workdir="${TEST_TMPDIR}/xilinx-python-primitives"
mkdir -p "${workdir}"
cp "${workspace}/example/xilinx.py" "${workdir}/"
cp "${workspace}/example/xilinx_register_slice_compact.v" "${workdir}/"
cp "${workspace}/example/xilinx_register_slice_mapped.v" "${workdir}/"
cp "${workspace}/example/test_config_verilog_xilinx_sec.yaml" "${workdir}/"
cp "${workspace}/example/vexriscv_genfull_xilinx.v" "${workdir}/"
cp "${workspace}/example/vexriscv_genfull_xilinx_different.v" "${workdir}/"
cp "${workspace}/example/test_config_verilog_vexriscv_genfull_xilinx_lec_different.yaml" "${workdir}/"

cd "${workdir}"
"${workspace}/src/bin/kepler-formal" \
  --config test_config_verilog_xilinx_sec.yaml

set +e
difference_output="$("${workspace}/src/bin/kepler-formal" \
  --config test_config_verilog_vexriscv_genfull_xilinx_lec_different.yaml 2>&1)"
difference_status=$?
set -e
if [[ ${difference_status} -ne 0 || "${difference_output}" != *"Difference was found."* ]]; then
  printf '%s\n' "${difference_output}" >&2
  exit 1
fi
