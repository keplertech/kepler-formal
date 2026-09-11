# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

{
  lib,
  stdenv,
  cmake,
  ninja,
  pkg-config,
  bison,
  flex,
  python312,
  boost,
  tbb,
  capnproto,
  zlib,
  fmt,
  tomlplusplus,
  autoPatchelfHook,
  makeWrapper,
  src,
}:

assert lib.assertMsg (lib.versionAtLeast boost.version "1.87") "Kepler requires Boost >= 1.87 for an offline Slang build";
assert lib.assertMsg (lib.versionAtLeast fmt.version "12.2") "Kepler requires fmt >= 12.2 for an offline Slang build";
assert lib.assertMsg (lib.versionAtLeast tomlplusplus.version "3.4") "Kepler requires tomlplusplus >= 3.4 for an offline Slang build";

stdenv.mkDerivation {
  pname = "kepler-formal";
  version = builtins.head (
    builtins.match
      ".*project\\(kepler-formal[[:space:]]+VERSION[[:space:]]+([0-9]+\\.[0-9]+\\.[0-9]+).*"
      (builtins.readFile "${src}/CMakeLists.txt")
  );
  inherit src;

  nativeBuildInputs = [
    cmake
    ninja
    pkg-config
    bison
    flex
    python312
    autoPatchelfHook
    makeWrapper
  ];
  buildInputs = [
    boost
    tbb
    capnproto
    zlib
    fmt
    tomlplusplus
    python312
    (lib.getLib stdenv.cc.cc)
  ];

  postPatch = ''
    # Both vendored SAT solvers configure and generate headers in their source
    # directories. The unpacked source is writable; all tools stay in the store.
    patchShebangs --build \
      thirdparty/kissat/configure thirdparty/kissat/scripts \
      thirdparty/cadical/configure thirdparty/cadical/scripts

    # Keep solver build descriptions independent of wall time and builder host.
    substituteInPlace \
      thirdparty/kissat/scripts/generate-build-header.sh \
      thirdparty/cadical/scripts/make-build-header.sh \
      --replace-fail 'date 2>/dev/null' 'date -u -d "@$SOURCE_DATE_EPOCH" 2>/dev/null' \
      --replace-fail 'uname -srmn 2>/dev/null' 'printf "%s" "${stdenv.hostPlatform.system}"'
    substituteInPlace thirdparty/kissat/scripts/generate-build-header.sh \
      --replace-fail 'DIR="`pwd`"' 'DIR="."'
  '';

  cmakeFlags = [
    "-DENABLE_UNIT_TESTS=OFF"
    "-DBUILD_KEPLER_PYTHON=OFF"
    "-DBUILD_NAJA_PYTHON=OFF"
    "-DBUILD_BENCHMARKS=OFF"
    "-DSLANG_INCLUDE_TESTS=OFF"
    "-DSLANG_INCLUDE_TOOLS=OFF"
    "-DSLANG_INCLUDE_INSTALL=OFF"
    "-DFETCHCONTENT_FULLY_DISCONNECTED=ON"
    "-DPython3_EXECUTABLE=${lib.getExe python312}"
    "-DPython_EXECUTABLE=${lib.getExe python312}"
    "-DCMAKE_INSTALL_LIBDIR=lib"
  ];

  enableParallelBuilding = true;

  # CMake copies naja.so beside the executable without rewriting its build
  # RPATH. Repair that module and the installed shared libraries together.
  preFixup = ''
    # Remove build paths before stdenv checks them; autoPatchelfHook then finds
    # the installed libraries and writes the module's runtime search path.
    patchelf --remove-rpath "$out/bin/naja.so"
    addAutoPatchelfSearchPath "$out/lib"
  '';
  postFixup = ''
    # The CLI embeds Python for technology primitives; its standard library
    # must be available even when no Python is installed in the user's profile.
    wrapProgram "$out/bin/kepler-formal" \
      --set PYTHONHOME "${python312}" \
      --prefix PYTHONPATH : "$out/lib/python"
  '';

  meta = {
    description = "Formal equivalence checking for digital circuits";
    homepage = "https://github.com/keplertech/kepler-formal";
    license = lib.licenses.gpl3Only;
    mainProgram = "kepler-formal";
    platforms = [ "x86_64-linux" ];
  };
}
