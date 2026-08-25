"""Data shared by binary resolution, execution, and the MCP facade."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Literal, cast

from typing_extensions import TypedDict


class CheckKind(str, Enum):
    GATE_LEC = "gate_lec"
    GATE_SEC = "gate_sec"
    RTL_SEC = "rtl_sec"


class ResultStatus(str, Enum):
    PASSED = "passed"
    CHECK_FAILED = "check_failed"
    INCONCLUSIVE = "inconclusive"
    SETUP_ERROR = "setup_error"
    CRASH = "crash"
    TIMEOUT = "timeout"


class ErrorKind(str, Enum):
    SETUP_ERROR = "setup_error"
    CHECK_FAILED = "check_failed"
    CRASH = "crash"
    TIMEOUT = "timeout"


class CheckResultPayload(TypedDict):
    comparison: Literal["gate_lec", "gate_sec", "rtl_sec"]
    status: Literal[
        "passed",
        "check_failed",
        "inconclusive",
        "setup_error",
        "crash",
        "timeout",
    ]
    passed: bool
    summary: str
    error_kind: Literal["setup_error", "check_failed", "crash", "timeout"] | None
    config_path: str
    timeout_seconds: float
    duration_seconds: float
    exit_code: int | None
    binary_path: str | None
    binary_source: Literal["path", "cache", "download"] | None
    command: list[str]
    stdout_tail: str
    stderr_tail: str
    schema_version: Literal["1"]


@dataclass(frozen=True, slots=True)
class BinaryResolution:
    path: Path
    source: str
    release_tag: str | None = None


@dataclass(frozen=True, slots=True)
class CheckResult:
    comparison: str
    status: str
    passed: bool
    summary: str
    error_kind: str | None
    config_path: str
    timeout_seconds: float
    duration_seconds: float
    exit_code: int | None = None
    binary_path: str | None = None
    binary_source: str | None = None
    command: tuple[str, ...] = ()
    stdout_tail: str = ""
    stderr_tail: str = ""
    schema_version: str = "1"

    def to_dict(self) -> CheckResultPayload:
        result = asdict(self)
        result["command"] = list(self.command)
        return cast(CheckResultPayload, result)


def error_kind_for(status: ResultStatus) -> str | None:
    if status is ResultStatus.CHECK_FAILED:
        return ErrorKind.CHECK_FAILED.value
    if status is ResultStatus.SETUP_ERROR:
        return ErrorKind.SETUP_ERROR.value
    if status is ResultStatus.CRASH:
        return ErrorKind.CRASH.value
    if status is ResultStatus.TIMEOUT:
        return ErrorKind.TIMEOUT.value
    return None
