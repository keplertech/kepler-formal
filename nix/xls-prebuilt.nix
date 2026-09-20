# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

{
  lib,
  stdenv,
  stdenvNoCC,
  fetchurl,
  python3,
  zstd,
  autoPatchelfHook,
  darwin,
  libcxx,
  zlib,
}:

let
  # These conda-forge builds share Abseil's lts_20260526 ABI. Pin the
  # complete set: substituting only Abseil breaks Protobuf and RE2 linkage.
  packages = {
    x86_64-linux = {
      subdir = "linux-64";
      files = {
        "libabseil-20260526.0-cxx17_h0dc7533_2.conda" = "7bb3d495411b7059e85225aae500d84e7846a4bb02ebd77e4450200b20c180f0";
        "libprotobuf-7.35.1-h2840a7c_2.conda" = "b5ac3938186516c1091b96414c34f90255da2a3bbd736a0712b22bbf4b5ebf02";
        "libre2-11-2025.11.05-h60473fc_2.conda" = "43132eb45a569b1f63199b0e7d5874d94227e9be950b327a323e7dcc1f8cd58c";
        "re2-2025.11.05-h94463f1_2.conda" = "7279908454f0c87b84ebc4aec6924f02f41a105d88420fb35dfaf5ecd976e0f5";
      };
    };
    aarch64-darwin = {
      subdir = "osx-arm64";
      files = {
        "libabseil-20260526.0-cxx17_h4c27e2a_2.conda" = "d9c62dae9d8f5bebadf72fcf369c00cc353183cb2521f9fffbf4c2f70b58bef9";
        "libprotobuf-7.35.1-h8daa630_2.conda" = "782d30325e61249e3240979e0f76d9dbf67867b1f38424cfdb1abb3a0f2cfa95";
        "libre2-11-2025.11.05-h5d15848_2.conda" = "06f49291605c379467987989ff08d44b470a23afeb690d7e078517871969bff7";
        "re2-2025.11.05-h5d60da5_2.conda" = "ba1154085d70e8e8f94cfc38a018a13d1a734ff1bcc2ab384e4fc3e532b23066";
      };
    };
  }.${stdenv.hostPlatform.system};
  archives = lib.mapAttrsToList (name: sha256: fetchurl {
    url = "https://conda.anaconda.org/conda-forge/${packages.subdir}/${name}";
    inherit sha256;
  }) packages.files;
in

assert lib.assertMsg (lib.versionAtLeast zlib.version "1.3.2") "Prebuilt XLS dependencies require zlib >= 1.3.2";
assert lib.assertMsg (!stdenv.hostPlatform.isDarwin || lib.versionAtLeast libcxx.version "19") "Prebuilt XLS dependencies require libc++ >= 19";

stdenvNoCC.mkDerivation {
  pname = "kepler-xls-prebuilt-dependencies";
  version = "20260526-protobuf-35.1-re2-2025.11.05";

  dontUnpack = true;
  dontConfigure = true;
  dontBuild = true;
  dontStrip = true;

  nativeBuildInputs = [ python3 zstd ]
    ++ lib.optional stdenv.hostPlatform.isLinux autoPatchelfHook
    ++ lib.optionals stdenv.hostPlatform.isDarwin [
      darwin.cctools
      darwin.autoSignDarwinBinariesHook
    ];
  buildInputs = [ zlib ]
    ++ lib.optional stdenv.hostPlatform.isLinux (lib.getLib stdenv.cc.cc)
    ++ lib.optional stdenv.hostPlatform.isDarwin libcxx;

  installPhase = ''
    runHook preInstall
    python3 ${./unpack-conda.py} "$out" ${lib.escapeShellArgs archives}
    runHook postInstall
  '';

  preFixup = lib.optionalString stdenv.hostPlatform.isLinux ''
    addAutoPatchelfSearchPath "$out/lib"
  '' + lib.optionalString stdenv.hostPlatform.isDarwin ''
    # Replace conda's @rpath references with store paths, including dylib IDs
    # used when Kepler links these libraries. The signing hook runs afterward.
    changes=(
      -change @rpath/libc++.1.dylib ${lib.getLib libcxx}/lib/libc++.1.dylib
      -change @rpath/libz.1.dylib ${lib.getLib zlib}/lib/libz.1.dylib
    )
    for dylib in "$out"/lib/*.dylib; do
      changes+=(-change "@rpath/$(basename "$dylib")" "$dylib")
    done
    for dylib in "$out"/lib/*.dylib; do
      if [ ! -L "$dylib" ]; then
        install_name_tool -id "$dylib" "''${changes[@]}" "$dylib"
      fi
    done
    for executable in "$out"/bin/*; do
      if [ ! -L "$executable" ]; then
        install_name_tool "''${changes[@]}" "$executable"
      fi
    done
  '';

  meta = {
    description = "Prebuilt, ABI-matched Abseil, Protobuf and RE2 libraries for XLS";
    homepage = "https://conda-forge.org/";
    license = [ lib.licenses.asl20 lib.licenses.bsd3 ];
    platforms = [ "x86_64-linux" "aarch64-darwin" ];
  };
}
