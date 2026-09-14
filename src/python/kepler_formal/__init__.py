# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Native Python interface for Kepler Formal."""

from __future__ import annotations

import importlib
import importlib.abc
import importlib.util
import sys

import najaeda as najaeda


class _NajaedaAliasLoader(importlib.abc.Loader):
    def __init__(self, canonical_name: str):
        self._canonical_name = canonical_name
        self._canonical_metadata = None

    def create_module(self, spec):
        module = importlib.import_module(self._canonical_name)
        self._canonical_metadata = (
            module.__spec__,
            module.__loader__,
            module.__package__,
        )
        return module

    def exec_module(self, module) -> None:
        # Import machinery temporarily initializes attributes for the alias
        # spec even though create_module() returned the canonical module.
        # Restore them so reload and importlib.resources keep using NajaEDA's
        # real loader and package search locations.
        module.__spec__, module.__loader__, module.__package__ = (
            self._canonical_metadata
        )


class _NajaedaAliasFinder(importlib.abc.MetaPathFinder):
    _kepler_najaeda_alias = True
    _prefix = f"{__name__}.najaeda."

    def find_spec(self, fullname, path=None, target=None):
        if not fullname.startswith(self._prefix):
            return None
        canonical_name = "najaeda." + fullname.removeprefix(self._prefix)
        canonical_spec = importlib.util.find_spec(canonical_name)
        if canonical_spec is None:
            return None
        return importlib.util.spec_from_loader(
            fullname,
            _NajaedaAliasLoader(canonical_name),
            is_package=canonical_spec.submodule_search_locations is not None,
        )


# NajaEDA must initialize first so its package can arrange platform DLL search
# paths and publish the native runtime capsule before Kepler's extension loads.
# Both package names intentionally resolve to the same module objects.
sys.modules[f"{__name__}.najaeda"] = najaeda
for _module_name, _module in tuple(sys.modules.items()):
    if _module_name.startswith("najaeda."):
        sys.modules[f"{__name__}.{_module_name}"] = _module
if not any(
    getattr(finder, "_kepler_najaeda_alias", False) for finder in sys.meta_path
):
    sys.meta_path.insert(0, _NajaedaAliasFinder())

from ._version import git_hash, version
from .api import (
    Design,
    InputFormat,
    NativeDesign,
    SecEncoding,
    SecEngine,
    Solver,
    VerificationMode,
    VerificationOptions,
    from_najaeda,
    run_cli,
    run_config,
    verify,
    verify_designs,
)
from .result import VerificationResult, VerificationStatus

__version__ = version()

__all__ = [
    "Design",
    "InputFormat",
    "NativeDesign",
    "SecEncoding",
    "SecEngine",
    "Solver",
    "VerificationMode",
    "VerificationOptions",
    "VerificationResult",
    "VerificationStatus",
    "from_najaeda",
    "git_hash",
    "najaeda",
    "run_cli",
    "run_config",
    "verify",
    "verify_designs",
    "version",
    "__version__",
]
