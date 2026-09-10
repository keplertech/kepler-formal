# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Pythonic, file-based interface to the native Kepler Formal engine."""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Sequence, TypeVar

from . import _native
from .result import VerificationResult

PathLike = str | os.PathLike[str]
_EnumType = TypeVar("_EnumType", bound=Enum)


class InputFormat(str, Enum):
    VERILOG = "verilog"
    SYSTEMVERILOG = "systemverilog"
    SV2V = "sv2v"
    NAJA_IF = "naja_if"


class VerificationMode(str, Enum):
    LEC = "lec"
    SEC = "sec"


class Solver(str, Enum):
    KISSAT = "kissat"
    CADICAL = "cadical"
    GLUCOSE = "glucose"


class SecEngine(str, Enum):
    PDR = "pdr"
    K_INDUCTION = "k_induction"
    IMC = "imc"


class SecEncoding(str, Enum):
    DUAL_RAIL_STEADY = "dual_rail_steady"
    BINARY = "binary"


@dataclass(frozen=True, slots=True)
class Design:
    """One side of a comparison."""

    files: PathLike | Sequence[PathLike] = field(default_factory=tuple)
    top: str | None = None
    flist: PathLike | None = None


@dataclass(frozen=True, slots=True)
class VerificationOptions:
    """Verification settings shared by the two designs."""

    input_format: InputFormat | str = InputFormat.VERILOG
    mode: VerificationMode | str = VerificationMode.LEC
    libraries: PathLike | Sequence[PathLike] = field(default_factory=tuple)
    solver: Solver | str = Solver.KISSAT
    max_k: int | None = None
    sec_engine: SecEngine | str | None = None
    sec_encoding: SecEncoding | str | None = None
    verilog_preprocessing: bool = False
    compact: bool = False
    allow_boundary_mismatch: bool = False
    report_skipped_outputs: bool = False
    log_file: PathLike | None = None
    log_level: str | None = None


def run_cli(arguments: Sequence[PathLike]) -> VerificationResult:
    """Run the native engine with CLI-style arguments (without ``argv[0]``).

    Calls are synchronous, serialized, and hold the GIL for the duration of
    the native run. Normal non-equivalence, partial, inconclusive, and
    unsupported verdicts are returned rather than raised.
    """

    if isinstance(arguments, (str, bytes, os.PathLike)):
        raise TypeError("arguments must be a sequence, not a single path")
    native_result = _native.run(
        [_raw_path(argument, f"arguments[{index}]", resolve=False)
         for index, argument in enumerate(arguments)]
    )
    return VerificationResult._from_native(native_result)


def run_config(path: PathLike) -> VerificationResult:
    """Run an existing Kepler Formal YAML or JSON configuration file."""

    return run_cli(("--config", _raw_path(path, "path")))


def verify(
    design1: Design | PathLike | Sequence[PathLike],
    design2: Design | PathLike | Sequence[PathLike],
    *,
    options: VerificationOptions | None = None,
) -> VerificationResult:
    """Compare two file-based designs and return a structured result."""

    first = design1 if isinstance(design1, Design) else Design(design1)
    second = design2 if isinstance(design2, Design) else Design(design2)
    settings = options or VerificationOptions()
    config = _build_config(first, second, settings)

    with tempfile.TemporaryDirectory(prefix="kepler_formal_python_") as directory:
        config_path = Path(directory) / "verify.json"
        config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
        return run_config(config_path)


def _enum_value(value: _EnumType | str, enum_type: type[_EnumType], label: str) -> str:
    try:
        return enum_type(value).value
    except (TypeError, ValueError) as error:
        choices = ", ".join(item.value for item in enum_type)
        raise ValueError(f"{label} must be one of: {choices}") from error


def _raw_path(value: PathLike, label: str, *, resolve: bool = True) -> str:
    try:
        raw = os.fspath(value)
    except TypeError as error:
        raise TypeError(f"{label} must be a string or path-like object") from error
    if not isinstance(raw, str):
        raise TypeError(f"{label} must resolve to a string path")
    if not raw:
        raise ValueError(f"{label} must not be empty")
    path = Path(raw).expanduser()
    return str(path.resolve()) if resolve else str(path)


def _paths(value: PathLike | Sequence[PathLike], label: str) -> list[str]:
    values: Sequence[PathLike]
    if isinstance(value, (str, os.PathLike)):
        values = (value,)
    else:
        values = value
    return [_raw_path(item, f"{label}[{index}]") for index, item in enumerate(values)]


def _optional_path(value: PathLike | None, label: str) -> str | None:
    if value is None:
        return None
    return _raw_path(value, label)


def _optional_text(value: str | None, label: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string")
    if not value:
        raise ValueError(f"{label} must not be empty")
    return value


def _boolean(value: bool, label: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{label} must be a bool")
    return value


def _build_config(
    design1: Design,
    design2: Design,
    options: VerificationOptions,
) -> Mapping[str, Any]:
    input_format = _enum_value(options.input_format, InputFormat, "input_format")
    mode = _enum_value(options.mode, VerificationMode, "mode")
    solver = _enum_value(options.solver, Solver, "solver")
    design1_files = _paths(design1.files, "design1.files")
    design2_files = _paths(design2.files, "design2.files")

    if not design1_files and design1.flist is None:
        raise ValueError("design1 needs at least one file or a flist")
    if not design2_files and design2.flist is None:
        raise ValueError("design2 needs at least one file or a flist")
    if input_format in {InputFormat.VERILOG.value, InputFormat.NAJA_IF.value}:
        if design1.flist is not None or design2.flist is not None:
            raise ValueError(f"{input_format} does not accept SystemVerilog flists")
    if input_format == InputFormat.NAJA_IF.value and (
        len(design1_files) != 1 or len(design2_files) != 1
    ):
        raise ValueError("naja_if requires exactly one snapshot per design")
    if mode == VerificationMode.LEC.value and input_format in {
        InputFormat.SYSTEMVERILOG.value,
        InputFormat.SV2V.value,
    }:
        raise ValueError("SystemVerilog input formats require SEC verification")
    if options.max_k is not None:
        if isinstance(options.max_k, bool) or not isinstance(options.max_k, int):
            raise TypeError("max_k must be an integer")
        if options.max_k < 0:
            raise ValueError("max_k must be non-negative")
    if mode == VerificationMode.LEC.value and any(
        value is not None
        for value in (options.max_k, options.sec_engine, options.sec_encoding)
    ):
        raise ValueError("SEC engine, encoding, and max_k cannot be used with LEC")

    config: dict[str, Any] = {
        "format": input_format,
        "verification": mode,
        "solver": solver,
        "input_paths": [design1_files, design2_files],
        "liberty_files": _paths(options.libraries, "libraries"),
        "verilog_preprocessing": _boolean(
            options.verilog_preprocessing, "verilog_preprocessing"
        ),
        "compact_mode": _boolean(options.compact, "compact"),
        "allow-boundary-mismatch": _boolean(
            options.allow_boundary_mismatch, "allow_boundary_mismatch"
        ),
        "report_skipped_pos": _boolean(
            options.report_skipped_outputs, "report_skipped_outputs"
        ),
    }

    if options.max_k is not None:
        config["max_k"] = options.max_k
    if options.sec_engine is not None:
        config["sec_engine"] = _enum_value(
            options.sec_engine, SecEngine, "sec_engine"
        )
    if options.sec_encoding is not None:
        config["sec_encoding"] = _enum_value(
            options.sec_encoding, SecEncoding, "sec_encoding"
        )
    if options.log_file is not None:
        config["log_file"] = _optional_path(options.log_file, "log_file")
    if options.log_level is not None:
        config["log_level"] = _optional_text(options.log_level, "log_level")

    if input_format == InputFormat.SYSTEMVERILOG.value:
        _set_optional(
            config, "sv_design1_flist", _optional_path(design1.flist, "design1.flist")
        )
        _set_optional(
            config, "sv_design2_flist", _optional_path(design2.flist, "design2.flist")
        )
        _set_optional(config, "sv_design1_top", _optional_text(design1.top, "design1.top"))
        _set_optional(config, "sv_design2_top", _optional_text(design2.top, "design2.top"))
    elif input_format == InputFormat.SV2V.value:
        _set_optional(
            config, "sv_design1_flist", _optional_path(design1.flist, "design1.flist")
        )
        if design2.flist is not None:
            raise ValueError("sv2v design2 is Verilog and cannot use an SV flist")
        _set_optional(config, "sv_design1_top", _optional_text(design1.top, "design1.top"))
        _set_optional(
            config,
            "verilog_design2_top",
            _optional_text(design2.top, "design2.top"),
        )
    elif input_format == InputFormat.VERILOG.value:
        _set_optional(
            config,
            "verilog_design1_top",
            _optional_text(design1.top, "design1.top"),
        )
        _set_optional(
            config,
            "verilog_design2_top",
            _optional_text(design2.top, "design2.top"),
        )
    elif design1.top is not None or design2.top is not None:
        raise ValueError("naja_if does not accept top-module selections")

    return config


def _set_optional(config: dict[str, Any], key: str, value: Any | None) -> None:
    if value is not None:
        config[key] = value
