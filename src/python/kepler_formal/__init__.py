# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Native Python interface for Kepler Formal."""

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
    "run_cli",
    "run_config",
    "verify",
    "version",
    "__version__",
]
