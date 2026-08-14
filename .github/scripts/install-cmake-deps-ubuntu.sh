#!/usr/bin/env bash
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set -euo pipefail

llvm_version="${KEPLER_LLVM_VERSION:-22}"
# shellcheck disable=SC1091
source /etc/os-release
ubuntu_codename="${VERSION_CODENAME:-jammy}"
llvm_keyring="/usr/share/keyrings/llvm-snapshot.gpg"
llvm_list="/etc/apt/sources.list.d/llvm-toolchain-${ubuntu_codename}-${llvm_version}.list"

sudo apt-get update
sudo apt-get install -yq ca-certificates curl gnupg lsb-release

if [[ ! -f "${llvm_keyring}" ]]; then
  curl -fsSL https://apt.llvm.org/llvm-snapshot.gpg.key \
    | sudo gpg --dearmor -o "${llvm_keyring}"
fi
printf 'deb [signed-by=%s] https://apt.llvm.org/%s/ llvm-toolchain-%s-%s main\n' \
  "${llvm_keyring}" "${ubuntu_codename}" "${ubuntu_codename}" "${llvm_version}" \
  | sudo tee "${llvm_list}" >/dev/null

sudo apt-get update
sudo apt-get install -yq \
  build-essential cmake ninja-build pkg-config curl ca-certificates \
  bison flex doxygen python3-dev \
  "clang-${llvm_version}" "clang-tools-${llvm_version}" \
  "llvm-${llvm_version}-dev" \
  "libclang-${llvm_version}-dev" "libclang-cpp${llvm_version}-dev" \
  libboost-dev libboost-iostreams-dev libfl-dev \
  capnproto libcapnp-dev libtbb-dev libspdlog-dev libfmt-dev \
  libgtest-dev libssl-dev \
  libz3-dev zlib1g-dev \
  "$@"

absl_version="20260107.0"
protobuf_version="${KEPLER_PROTOBUF_VERSION:-35.1}"
re2_version="${KEPLER_RE2_VERSION:-2025-11-05}"
ortools_version="${KEPLER_ORTOOLS_VERSION:-9.15}"
absl_prefix="/usr/local"
absl_config="${absl_prefix}/lib/cmake/absl/abslConfig.cmake"

if [[ ! -f "${absl_config}" ]]; then
  work_dir="${RUNNER_TEMP:-/tmp}/abseil-cpp-${absl_version}"
  src_dir="${work_dir}/src"
  build_dir="${work_dir}/build"

  rm -rf "${work_dir}"
  mkdir -p "${work_dir}"
  curl -L "https://github.com/abseil/abseil-cpp/archive/refs/tags/${absl_version}.tar.gz" \
    -o "${work_dir}/abseil-cpp.tar.gz"
  tar -xzf "${work_dir}/abseil-cpp.tar.gz" -C "${work_dir}"
  mv "${work_dir}/abseil-cpp-${absl_version}" "${src_dir}"

  cmake -S "${src_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=20 \
    -DCMAKE_C_COMPILER="clang-${llvm_version}" \
    -DCMAKE_CXX_COMPILER="clang++-${llvm_version}" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DABSL_BUILD_TESTING=OFF \
    -DABSL_ENABLE_INSTALL=ON \
    -DABSL_PROPAGATE_CXX_STD=ON \
    -DCMAKE_INSTALL_PREFIX="${absl_prefix}"
  cmake --build "${build_dir}" -j "$(nproc)"
  sudo cmake --install "${build_dir}"
  sudo ldconfig
fi

re2_prefix="/usr/local"
re2_header="${re2_prefix}/include/re2/re2.h"

if [[ ! -f "${re2_header}" ]]; then
  work_dir="${RUNNER_TEMP:-/tmp}/re2-${re2_version}"
  src_dir="${work_dir}/src"
  build_dir="${work_dir}/build"

  rm -rf "${work_dir}"
  mkdir -p "${work_dir}"
  curl -L "https://github.com/google/re2/archive/refs/tags/${re2_version}.tar.gz" \
    -o "${work_dir}/re2.tar.gz"
  tar -xzf "${work_dir}/re2.tar.gz" -C "${work_dir}"
  mv "${work_dir}/re2-${re2_version}" "${src_dir}"

  cmake -S "${src_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=20 \
    -DCMAKE_C_COMPILER="clang-${llvm_version}" \
    -DCMAKE_CXX_COMPILER="clang++-${llvm_version}" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DRE2_BUILD_TESTING=OFF \
    -DCMAKE_PREFIX_PATH="${absl_prefix}" \
    -DCMAKE_INSTALL_PREFIX="${re2_prefix}"
  cmake --build "${build_dir}" -j "$(nproc)"
  sudo cmake --install "${build_dir}"
  sudo ldconfig
fi

protobuf_prefix="/usr/local"
protobuf_config="${protobuf_prefix}/lib/cmake/protobuf/protobuf-config.cmake"

if [[ ! -f "${protobuf_config}" ]]; then
  work_dir="${RUNNER_TEMP:-/tmp}/protobuf-${protobuf_version}"
  src_dir="${work_dir}/src"
  build_dir="${work_dir}/build"

  rm -rf "${work_dir}"
  mkdir -p "${work_dir}"
  curl -L "https://github.com/protocolbuffers/protobuf/releases/download/v${protobuf_version}/protobuf-${protobuf_version}.tar.gz" \
    -o "${work_dir}/protobuf.tar.gz"
  tar -xzf "${work_dir}/protobuf.tar.gz" -C "${work_dir}"
  mv "${work_dir}/protobuf-${protobuf_version}" "${src_dir}"

  cmake -S "${src_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=20 \
    -DCMAKE_C_COMPILER="clang-${llvm_version}" \
    -DCMAKE_CXX_COMPILER="clang++-${llvm_version}" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -Dprotobuf_BUILD_TESTS=OFF \
    -DBUILD_SHARED_LIBS=ON \
    -Dprotobuf_BUILD_SHARED_LIBS=ON \
    -Dprotobuf_ABSL_PROVIDER=package \
    -Dprotobuf_BUILD_PROTOC_BINARIES=ON \
    -Dprotobuf_INSTALL=ON \
    -DCMAKE_PREFIX_PATH="${absl_prefix}" \
    -DCMAKE_INSTALL_PREFIX="${protobuf_prefix}"
  cmake --build "${build_dir}" -j "$(nproc)"
  sudo cmake --install "${build_dir}"
  sudo ldconfig
fi

ortools_prefix="/usr/local"
ortools_header="${ortools_prefix}/include/ortools/graph/graph.h"

if [[ ! -f "${ortools_header}" ]]; then
  work_dir="${RUNNER_TEMP:-/tmp}/or-tools-${ortools_version}"
  src_dir="${work_dir}/src"

  rm -rf "${work_dir}"
  mkdir -p "${src_dir}"
  curl -L \
    "https://github.com/google/or-tools/releases/download/v${ortools_version}/or-tools-${ortools_version}.tar.gz" \
    -o "${work_dir}/or-tools.tar.gz"
  tar -xzf "${work_dir}/or-tools.tar.gz" \
    -C "${src_dir}" \
    --strip-components=1

  sudo mkdir -p "${ortools_prefix}/include"
  sudo cp -R "${src_dir}/ortools" "${ortools_prefix}/include/"
fi

if [[ -n "${GITHUB_ENV:-}" ]]; then
  echo "CPATH=${ortools_prefix}/include${CPATH:+:${CPATH}}" >> "${GITHUB_ENV}"
  echo "CPLUS_INCLUDE_PATH=${ortools_prefix}/include${CPLUS_INCLUDE_PATH:+:${CPLUS_INCLUDE_PATH}}" >> "${GITHUB_ENV}"
  echo "CXXFLAGS=-I${ortools_prefix}/include${CXXFLAGS:+ ${CXXFLAGS}}" >> "${GITHUB_ENV}"
fi
