# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Owning result values returned by the native verifier."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any, Mapping


class VerificationStatus(str, Enum):
    """Semantic outcome of one completed invocation."""

    NO_RESULT = "no_result"
    EQUIVALENT = "equivalent"
    DIFFERENT = "different"
    PARTIALLY_PROVED = "partially_proved"
    INCONCLUSIVE = "inconclusive"
    UNSUPPORTED = "unsupported"
    ERROR = "error"


@dataclass(frozen=True, slots=True)
class VerificationResult:
    """A value-only result that remains valid after native state is released."""

    status: VerificationStatus
    exit_code: int
    input_format: str | None = None
    verification: str | None = None
    log_file: str | None = None
    bound: int = 0
    reason: str | None = None
    covered_outputs: int = 0
    total_outputs: int = 0
    proven_outputs: int = 0
    unproven_outputs: tuple[str, ...] = ()
    skipped_observed_outputs: tuple[str, ...] = ()

    @property
    def equivalent(self) -> bool:
        return self.status is VerificationStatus.EQUIVALENT

    @property
    def conclusive(self) -> bool:
        return self.status in {
            VerificationStatus.EQUIVALENT,
            VerificationStatus.DIFFERENT,
        }

    @property
    def coverage_percent(self) -> float | None:
        if self.total_outputs == 0:
            return None
        return 100.0 * self.covered_outputs / self.total_outputs

    @classmethod
    def _from_native(cls, value: Mapping[str, Any]) -> "VerificationResult":
        return cls(
            status=VerificationStatus(value["status"]),
            exit_code=int(value["exit_code"]),
            input_format=value.get("input_format"),
            verification=value.get("verification"),
            log_file=value.get("log_file"),
            bound=int(value.get("bound", 0)),
            reason=value.get("reason"),
            covered_outputs=int(value.get("covered_outputs", 0)),
            total_outputs=int(value.get("total_outputs", 0)),
            proven_outputs=int(value.get("proven_outputs", 0)),
            unproven_outputs=tuple(value.get("unproven_outputs", ())),
            skipped_observed_outputs=tuple(
                value.get("skipped_observed_outputs", ())
            ),
        )
