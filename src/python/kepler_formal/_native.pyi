# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

from typing import Mapping, TypedDict


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


class NativeDesign:
    @property
    def source(self) -> object: ...

    @property
    def najaeda_design(self) -> object: ...


def from_najaeda(design: object, source: object = ...) -> NativeDesign: ...
def verify_designs(
    design1: NativeDesign,
    design2: NativeDesign,
    options: Mapping[str, object],
) -> NativeResult: ...
def get_version() -> str: ...
def get_git_hash() -> str: ...
