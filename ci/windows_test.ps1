# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

$ErrorActionPreference = 'Stop'
$env:PYTHONPATH = ''
$env:PYTHONFAULTHANDLER = '1'

python "$PSScriptRoot/check_kepler_wheel.py"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

# Stop on the first ordinary error so its traceback is printed before another
# test can crash the process. Faulthandler also reports native access violations.
python -m unittest discover -v -f "$PSScriptRoot/../test/python"
if ($LASTEXITCODE -ne 0) {
    $keplerTestExitCode = $LASTEXITCODE
    # A native crash bypasses TemporaryDirectory cleanup. Keep the last native
    # messages in CI output even when unittest cannot finish its error report.
    Get-ChildItem -LiteralPath $env:TEMP -Directory -Filter 'kepler_formal_*' -ErrorAction SilentlyContinue |
        ForEach-Object {
            Get-ChildItem -LiteralPath $_.FullName -File -Recurse -Filter '*.log' -ErrorAction SilentlyContinue
        } |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 5 |
        ForEach-Object {
            Write-Host "Native test log: $($_.FullName)"
            Get-Content -LiteralPath $_.FullName -Tail 80 -ErrorAction Continue
        }
    exit $keplerTestExitCode
}
