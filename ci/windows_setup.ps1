# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Invoke-Checked {
    param([string]$Program, [string[]]$Arguments)
    & $Program @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Program failed with exit code $LASTEXITCODE"
    }
}

# Match NajaEDA's x64-windows (MSVC ABI, dynamic CRT) dependency layout.
$keplerVcpkgRoot = Join-Path $env:USERPROFILE 'vcpkg'
if (-not (Test-Path (Join-Path $keplerVcpkgRoot '.git'))) {
    if (Test-Path $keplerVcpkgRoot) {
        throw "Refusing to overwrite non-repository directory $keplerVcpkgRoot"
    }
    Invoke-Checked 'git' @('clone', 'https://github.com/microsoft/vcpkg.git', $keplerVcpkgRoot)
}
Invoke-Checked (Join-Path $keplerVcpkgRoot 'bootstrap-vcpkg.bat') @('-disableMetrics')
Invoke-Checked (Join-Path $keplerVcpkgRoot 'vcpkg.exe') @(
    'install', '--triplet=x64-windows',
    'capnproto', 'tbb', 'zlib',
    'boost-intrusive', 'boost-dynamic-bitset', 'boost-unordered', 'boost-regex'
)

# GitHub's Windows image includes LLVM and an MSVC SDK. The LLVM frontend is
# needed for GNU language extensions in the solvers; its ABI remains MSVC.
$keplerClangCl = Join-Path $env:ProgramFiles 'LLVM/bin/clang-cl.exe'
if (-not (Test-Path $keplerClangCl)) {
    throw "clang-cl is required at $keplerClangCl; install LLVM and the MSVC C++ SDK"
}
Invoke-Checked $keplerClangCl @('--version')

$keplerEnvironment = @{
    CMAKE_TOOLCHAIN_FILE = (Join-Path $keplerVcpkgRoot 'scripts/buildsystems/vcpkg.cmake')
    VCPKG_ROOT = $keplerVcpkgRoot
    VCPKG_DEFAULT_TRIPLET = 'x64-windows'
    CMAKE_GENERATOR = 'Ninja'
    CC = $keplerClangCl
    CXX = $keplerClangCl
}
foreach ($entry in $keplerEnvironment.GetEnumerator()) {
    [Environment]::SetEnvironmentVariable($entry.Key, $entry.Value, 'Process')
    if ($env:GITHUB_ENV) {
        Add-Content -Path $env:GITHUB_ENV -Value "$($entry.Key)=$($entry.Value)"
    }
}
Write-Host "Windows dependencies ready in $keplerVcpkgRoot (x64-windows)"
