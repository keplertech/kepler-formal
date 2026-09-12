# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Native Python interface for Kepler Formal."""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING, Any

from ._version import git_hash, version
from .api import (
    Design,
    InputFormat,
    SecEncoding,
    SecEngine,
    Solver,
    VerificationMode,
    VerificationOptions,
    run_cli,
    run_config,
    verify,
)
from .result import VerificationResult, VerificationStatus

__version__ = version()

if TYPE_CHECKING:
    from . import najaeda as najaeda


def __getattr__(name: str) -> Any:
    """Lazily expose the bundled, isolated NajaEDA editor package."""

    if name == "najaeda":
        module = importlib.import_module(".najaeda", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted({*globals(), "najaeda"})

__all__ = [
    "Design",
    "InputFormat",
    "SecEncoding",
    "SecEngine",
    "Solver",
    "VerificationMode",
    "VerificationOptions",
    "VerificationResult",
    "VerificationStatus",
    "git_hash",
    "najaeda",
    "run_cli",
    "run_config",
    "verify",
    "version",
    "__version__",
]
