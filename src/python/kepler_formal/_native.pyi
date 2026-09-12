# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

from collections.abc import Sequence
from os import PathLike
from typing import TypedDict


class NativeResult(TypedDict):
    status: str
    exit_code: int
    input_format: str | None
    verification: str | None
    log_file: str | None
    bound: int
    reason: str | None
    covered_outputs: int
    total_outputs: int
    proven_outputs: int
    unproven_outputs: list[str]
    skipped_observed_outputs: list[str]


def run(arguments: Sequence[str | PathLike[str]]) -> NativeResult: ...
def get_version() -> str: ...
def get_git_hash() -> str: ...
