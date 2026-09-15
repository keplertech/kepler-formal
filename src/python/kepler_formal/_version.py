# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

from . import _native


def version() -> str:
    """Return the Kepler Formal version compiled into the extension."""

    return _native.get_version()


def git_hash() -> str:
    """Return the source revision recorded at build time."""

    return _native.get_git_hash()
