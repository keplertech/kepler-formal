# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

$ErrorActionPreference = 'Stop'
$env:PYTHONPATH = ''

python "$PSScriptRoot/check_kepler_wheel.py"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

python -m unittest discover -v "$PSScriptRoot/../test/python"
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
