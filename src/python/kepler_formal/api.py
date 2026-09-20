# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Verification of live NajaEDA designs through the shared native engine."""

from __future__ import annotations

import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, TypeVar

import najaeda as _najaeda

from . import _native
from .result import VerificationResult

PathLike = str | os.PathLike[str]
_EnumType = TypeVar("_EnumType", bound=Enum)


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
class VerificationOptions:
    """Verification settings shared by the two designs.

    ``set_as_boundary`` pairs top-relative paths to leaf instances, whose
    models have no child instances. Hierarchical paths to leaves are valid;
    selecting a nonleaf instance is rejected.
    """

    mode: VerificationMode | str = VerificationMode.LEC
    solver: Solver | str = Solver.KISSAT
    max_k: int | None = None
    sec_engine: SecEngine | str | None = None
    sec_encoding: SecEncoding | str | None = None
    allow_boundary_mismatch: bool = False
    report_skipped_outputs: bool = False
    log_file: PathLike | None = None
    log_level: str | None = None
    set_as_boundary: (
        list[tuple[str, str]] | tuple[tuple[str, str], ...]
    ) = ()


NativeDesign = _native.NativeDesign


def from_najaeda(design: object) -> NativeDesign:
    """Capture a live NajaEDA design as a stable, zero-copy native handle.

    ``design`` may be a raw ``najaeda.naja.SNLDesign`` or a high-level
    ``najaeda.netlist.Instance``.  An ``Instance`` is resolved to its current
    model when this function is called, so later changes to NajaEDA's selected
    top design do not retarget the returned handle.  The native netlist is not
    cloned or serialized, and later edits to that same design remain visible.
    """

    if isinstance(design, _najaeda.naja.SNLDesign):
        raw_design = design
    else:
        from najaeda import netlist

        if not isinstance(design, netlist.Instance):
            raise TypeError(
                "design must be a najaeda.naja.SNLDesign or "
                "najaeda.netlist.Instance"
            )
        universe = _najaeda.naja.NLUniverse.get()
        if universe is None:
            raise ReferenceError("the NajaEDA instance has no live universe")
        try:
            identity = design.get_model_id()
            raw_design = universe.getSNLDesign(identity)
        except (AttributeError, RuntimeError) as error:
            raise ReferenceError(
                "the NajaEDA instance no longer resolves to a live design"
            ) from error
        if raw_design is None:
            raise ReferenceError(
                "the NajaEDA instance no longer resolves to a live design"
            )
    return _native.from_najaeda(raw_design, design)


def verify_designs(
    design1: NativeDesign | object,
    design2: NativeDesign | object,
    *,
    options: VerificationOptions | None = None,
) -> VerificationResult:
    """Compare two live NajaEDA designs without files, copying, or rebuilding.

    Each argument must be a :class:`NativeDesign` or a raw
    ``najaeda.naja.SNLDesign``.  Capture high-level ``Instance`` objects first
    with :func:`from_najaeda`; this freezes which model the instance denotes,
    while retaining the original object for the synchronous native call.
    Selected instance pins act as logical verification boundaries; the native
    designs are analyzed directly without modification.
    The caller owns both netlists; verification leaves them available for
    further edits and calls, including when verification reports an error.
    """

    native_options = _build_native_design_options(options)
    first = _as_native_design(design1, "design1")
    second = _as_native_design(design2, "design2")
    return VerificationResult._from_native(
        _native.verify_designs(first, second, native_options)
    )


def _as_native_design(value: object, label: str) -> NativeDesign:
    if isinstance(value, NativeDesign):
        return value
    if isinstance(value, _najaeda.naja.SNLDesign):
        return from_najaeda(value)

    from najaeda import netlist

    if isinstance(value, netlist.Instance):
        raise TypeError(
            f"{label} is a najaeda.netlist.Instance; capture it with "
            "from_najaeda() before changing the current top design"
        )
    raise TypeError(f"{label} must be a NativeDesign or najaeda.naja.SNLDesign")


def _build_native_design_options(
    options: VerificationOptions | None,
) -> Mapping[str, Any]:
    if options is None:
        settings = VerificationOptions()
    elif isinstance(options, VerificationOptions):
        settings = options
    else:
        raise TypeError("options must be a VerificationOptions instance or None")

    mode = _enum_value(settings.mode, VerificationMode, "mode")
    solver = _enum_value(settings.solver, Solver, "solver")
    allow_boundary_mismatch = _boolean(
        settings.allow_boundary_mismatch, "allow_boundary_mismatch"
    )
    report_skipped_outputs = _boolean(
        settings.report_skipped_outputs, "report_skipped_outputs"
    )
    set_as_boundary = _boundary_pairs(settings.set_as_boundary)
    if settings.max_k is not None:
        if isinstance(settings.max_k, bool) or not isinstance(settings.max_k, int):
            raise TypeError("max_k must be an integer")
        if settings.max_k < 0:
            raise ValueError("max_k must be non-negative")
    if mode == VerificationMode.LEC.value and any(
        value is not None
        for value in (settings.max_k, settings.sec_engine, settings.sec_encoding)
    ):
        raise ValueError("SEC engine, encoding, and max_k cannot be used with LEC")
    if mode == VerificationMode.SEC.value and allow_boundary_mismatch:
        raise ValueError("allow_boundary_mismatch is only supported for LEC")

    sec_engine = (
        SecEngine.PDR.value
        if settings.sec_engine is None
        else _enum_value(settings.sec_engine, SecEngine, "sec_engine")
    )
    sec_encoding = (
        SecEncoding.DUAL_RAIL_STEADY.value
        if settings.sec_encoding is None
        else _enum_value(settings.sec_encoding, SecEncoding, "sec_encoding")
    )
    log_file = _optional_path(settings.log_file, "log_file")
    log_level = _optional_text(settings.log_level, "log_level")
    if log_level not in {None, "debug", "info"}:
        raise ValueError("log_level must be 'debug', 'info', or None")

    return {
        "mode": mode,
        "solver": solver,
        "max_k": 32 if settings.max_k is None else settings.max_k,
        "sec_engine": sec_engine,
        "sec_encoding": sec_encoding,
        "allow_boundary_mismatch": allow_boundary_mismatch,
        "set_as_boundary": set_as_boundary,
        "report_skipped_outputs": report_skipped_outputs,
        "log_file": log_file,
        "log_level": log_level,
    }


def _enum_value(value: _EnumType | str, enum_type: type[_EnumType], label: str) -> str:
    try:
        return enum_type(value).value
    except (TypeError, ValueError) as error:
        choices = ", ".join(item.value for item in enum_type)
        raise ValueError(f"{label} must be one of: {choices}") from error


def _argument_string(value: PathLike, label: str) -> str:
    try:
        raw = os.fspath(value)
    except TypeError as error:
        raise TypeError(f"{label} must be a string or path-like object") from error
    if not isinstance(raw, str):
        raise TypeError(f"{label} must resolve to a string")
    return raw


def _raw_path(value: PathLike, label: str) -> str:
    raw = _argument_string(value, label)
    if not raw:
        raise ValueError(f"{label} must not be empty")
    return str(Path(raw).expanduser().resolve())


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


def _boundary_pairs(value: object) -> list[tuple[str, str]]:
    if not isinstance(value, (list, tuple)):
        raise TypeError("set_as_boundary must be a list or tuple of path pairs")
    result: list[tuple[str, str]] = []
    for index, pair in enumerate(value):
        if not isinstance(pair, (list, tuple)):
            raise TypeError(
                f"set_as_boundary[{index}] must be a list or tuple of two paths"
            )
        if len(pair) != 2:
            raise ValueError(
                f"set_as_boundary[{index}] must contain exactly two paths"
            )
        paths: list[str] = []
        for side, path in enumerate(pair):
            label = f"set_as_boundary[{index}][{side}]"
            if not isinstance(path, str):
                raise TypeError(f"{label} must be a string")
            if not path:
                raise ValueError(f"{label} must not be empty")
            if "\0" in path:
                raise ValueError(f"{label} cannot contain NUL bytes")
            paths.append(path)
        result.append((paths[0], paths[1]))
    return result


def _boolean(value: bool, label: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{label} must be a bool")
    return value
