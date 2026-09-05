"""One-shot execution and verdict parsing for kepler-formal."""

from __future__ import annotations

import asyncio
import contextlib
import math
import os
import re
import signal
import time
from collections.abc import Awaitable, Callable
from pathlib import Path

from .models import (
    BinaryResolution,
    CheckKind,
    CheckResult,
    ResultStatus,
    error_kind_for,
)
from .resolver import BinaryResolver, BinarySetupError


MAX_CAPTURE_BYTES = 64 * 1024
DEFAULT_PROGRESS_INTERVAL_SECONDS = 15.0
ProgressCallback = Callable[[float, float, str], Awaitable[None]]


class CheckRunner:
    def __init__(
        self,
        project_root: Path,
        resolver: BinaryResolver | None = None,
        *,
        progress_interval_seconds: float = DEFAULT_PROGRESS_INTERVAL_SECONDS,
    ) -> None:
        self.project_root = project_root.expanduser().resolve()
        self.resolver = resolver or BinaryResolver()
        self.progress_interval_seconds = progress_interval_seconds

    async def run(
        self,
        kind: CheckKind,
        config_path: str,
        timeout_seconds: float,
        progress: ProgressCallback | None = None,
    ) -> CheckResult:
        started = time.monotonic()
        try:
            resolved_config = self._resolve_path(config_path)
        except (OSError, ValueError) as exc:
            return self._result(
                kind,
                ResultStatus.CRASH,
                f"invalid configuration path {config_path!r}: {exc}",
                Path(config_path),
                timeout_seconds if math.isfinite(timeout_seconds) else 0.0,
                started,
            )
        if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
            return self._result(
                kind,
                ResultStatus.CRASH,
                f"timeout_seconds must be a finite positive number, got {timeout_seconds!r}",
                resolved_config,
                timeout_seconds if math.isfinite(timeout_seconds) else 0.0,
                started,
            )
        if not resolved_config.is_file():
            return self._result(
                kind,
                ResultStatus.CRASH,
                f"configuration file does not exist: {resolved_config}",
                resolved_config,
                timeout_seconds,
                started,
            )

        try:
            binary = await asyncio.to_thread(self.resolver.resolve)
        except BinarySetupError as exc:
            return self._result(
                kind,
                ResultStatus.SETUP_ERROR,
                str(exc),
                resolved_config,
                timeout_seconds,
                started,
            )
        except Exception as exc:  # defensive wrapper boundary
            return self._result(
                kind,
                ResultStatus.SETUP_ERROR,
                f"unexpected binary resolution failure: {exc}",
                resolved_config,
                timeout_seconds,
                started,
            )

        command = (str(binary.path), "--config", str(resolved_config))
        try:
            process = await asyncio.create_subprocess_exec(
                *command,
                cwd=self.project_root,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                start_new_session=(os.name != "nt"),
            )
        except (OSError, ValueError) as exc:
            return self._result(
                kind,
                ResultStatus.CRASH,
                f"failed to start kepler-formal: {exc}",
                resolved_config,
                timeout_seconds,
                started,
                binary=binary,
                command=command,
            )

        stdout_task = asyncio.create_task(self._read_tail(process.stdout))
        stderr_task = asyncio.create_task(self._read_tail(process.stderr))
        progress_task = asyncio.create_task(
            self._report_progress(progress, timeout_seconds, started)
        )
        timed_out = False
        try:
            await asyncio.wait_for(process.wait(), timeout=timeout_seconds)
        except asyncio.TimeoutError:
            timed_out = True
            self._terminate(process)
            await process.wait()
        except asyncio.CancelledError:
            self._terminate(process)
            await process.wait()
            await asyncio.gather(stdout_task, stderr_task, return_exceptions=True)
            raise
        finally:
            progress_task.cancel()
            with contextlib.suppress(asyncio.CancelledError):
                await progress_task

        stdout = (await stdout_task).decode("utf-8", errors="replace")
        stderr = (await stderr_task).decode("utf-8", errors="replace")
        if timed_out:
            return self._result(
                kind,
                ResultStatus.TIMEOUT,
                f"{kind.value} exceeded its {timeout_seconds:g}-second timeout and was terminated",
                resolved_config,
                timeout_seconds,
                started,
                exit_code=process.returncode,
                binary=binary,
                command=command,
                stdout=stdout,
                stderr=stderr,
            )

        status, summary = classify_verdict(kind, process.returncode, stdout, stderr)
        return self._result(
            kind,
            status,
            summary,
            resolved_config,
            timeout_seconds,
            started,
            exit_code=process.returncode,
            binary=binary,
            command=command,
            stdout=stdout,
            stderr=stderr,
        )

    def _resolve_path(self, value: str) -> Path:
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = self.project_root / path
        return path.resolve()

    async def _report_progress(
        self,
        callback: ProgressCallback | None,
        timeout_seconds: float,
        started: float,
    ) -> None:
        if callback is None:
            return
        while True:
            elapsed = time.monotonic() - started
            try:
                await callback(
                    min(elapsed, timeout_seconds),
                    timeout_seconds,
                    f"kepler-formal is still running ({elapsed:.0f}s elapsed)",
                )
            except Exception:
                return
            await asyncio.sleep(self.progress_interval_seconds)

    @staticmethod
    async def _read_tail(reader: asyncio.StreamReader | None) -> bytes:
        if reader is None:
            return b""
        tail = b""
        while True:
            chunk = await reader.read(8192)
            if not chunk:
                return tail
            tail = (tail + chunk)[-MAX_CAPTURE_BYTES:]

    @staticmethod
    def _terminate(process: asyncio.subprocess.Process) -> None:
        if os.name != "nt" and process.pid is not None:
            with contextlib.suppress(ProcessLookupError):
                os.killpg(process.pid, signal.SIGKILL)
                return
        with contextlib.suppress(ProcessLookupError):
            process.kill()

    @staticmethod
    def _result(
        kind: CheckKind,
        status: ResultStatus,
        summary: str,
        config_path: Path,
        timeout_seconds: float,
        started: float,
        *,
        exit_code: int | None = None,
        binary: BinaryResolution | None = None,
        command: tuple[str, ...] = (),
        stdout: str = "",
        stderr: str = "",
    ) -> CheckResult:
        return CheckResult(
            comparison=kind.value,
            status=status.value,
            passed=status is ResultStatus.PASSED,
            summary=summary,
            error_kind=error_kind_for(status),
            config_path=str(config_path),
            timeout_seconds=timeout_seconds,
            duration_seconds=round(time.monotonic() - started, 3),
            exit_code=exit_code,
            binary_path=str(binary.path) if binary else None,
            binary_source=binary.source if binary else None,
            command=command,
            stdout_tail=stdout,
            stderr_tail=stderr,
        )


_ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def classify_verdict(
    kind: CheckKind, exit_code: int | None, stdout: str, stderr: str
) -> tuple[ResultStatus, str]:
    output = _ANSI_ESCAPE.sub("", "\n".join((stdout, stderr)))
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    mode_line = _last_matching(lines, r"\bVerification:\s*(lec|sec)\b")
    expected_mode = "lec" if kind is CheckKind.GATE_LEC else "sec"
    if mode_line:
        mode_match = re.search(r"\bVerification:\s*(lec|sec)\b", mode_line, re.I)
        actual_mode = mode_match.group(1).lower() if mode_match else "unknown"
        if actual_mode != expected_mode:
            return (
                ResultStatus.CRASH,
                f"{kind.value} expected verification mode {expected_mode}, but the "
                f"checker ran in {actual_mode} mode",
            )

    if kind is CheckKind.GATE_LEC:
        difference = _last_matching(lines, r"(?<!No )Difference was found")
        equivalent = _last_matching(lines, r"No difference was found")
        if difference:
            return ResultStatus.CHECK_FAILED, difference
        if exit_code == 0 and equivalent:
            return ResultStatus.PASSED, equivalent
        return ResultStatus.CRASH, _crash_summary(exit_code, lines)

    if exit_code == 0:
        return ResultStatus.PASSED, _last_matching(
            lines, r"SEC proved equivalence|No difference was found"
        ) or "SEC proved the designs equivalent."
    if exit_code == 3:
        return ResultStatus.CHECK_FAILED, _last_matching(
            lines, r"SEC found a counterexample|Difference was found"
        ) or "SEC found a counterexample."
    if exit_code == 2:
        return ResultStatus.INCONCLUSIVE, _last_matching(
            lines, r"SEC was inconclusive|SEC cannot run"
        ) or "SEC completed without a proof or counterexample."
    if exit_code == 1 and _last_matching(lines, r"SEC partially proved equivalence"):
        return ResultStatus.INCONCLUSIVE, _last_matching(
            lines, r"SEC partially proved equivalence"
        ) or "SEC partially proved equivalence; remaining outputs are inconclusive."
    return ResultStatus.CRASH, _crash_summary(exit_code, lines)


def _last_matching(lines: list[str], pattern: str) -> str | None:
    matcher = re.compile(pattern, re.IGNORECASE)
    return next((line for line in reversed(lines) if matcher.search(line)), None)


def _crash_summary(exit_code: int | None, lines: list[str]) -> str:
    detail = lines[-1] if lines else "no diagnostic output"
    return f"kepler-formal exited abnormally with code {exit_code}: {detail}"
