# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0
"""Fail closed when Kepler's native extension loads another NajaEDA runtime."""
from __future__ import annotations

from functools import lru_cache
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import re


@lru_cache(maxsize=4)
def validate_provider(manifest_json: str, version: str, git_hash: str) -> None:
    """Validate the exact provider binaries used to compile this extension.

    Python package versions alone do not identify a C++ ABI. The extension is
    compiled for a particular wheel's native files; their hashes travel with
    Kepler's extension. Successful validation is cached for this process.
    """
    import najaeda
    from najaeda import naja

    try:
        try:
            distribution = importlib.metadata.distribution("najaeda")
        except importlib.metadata.PackageNotFoundError:
            # A plain CMake install, as in the source regression, has no metadata.
            distribution = None
        installed_version = naja.getVersion() if distribution is None else distribution.version
        if installed_version != version or naja.getGitHash() != git_hash:
            raise ImportError(
                f"Kepler Formal was built against NajaEDA {version} ({git_hash}); "
                f"found {installed_version} ({naja.getGitHash()})")
        root = Path(najaeda.__file__).resolve().parent
        if (distribution is not None
                and distribution.locate_file("najaeda/__init__.py").resolve() != root / "__init__.py"):
            raise ImportError("Imported NajaEDA does not match its installed distribution")
        extension = Path(naja.__file__).resolve()
        if extension.parent != root:
            raise ImportError("NajaEDA imported a native extension outside its package")
        manifest = json.loads(manifest_json)
        if not isinstance(manifest, dict) or not manifest:
            raise ImportError("Kepler contains an invalid NajaEDA provider manifest")
        extension_key = extension.relative_to(root.parent).as_posix()
        if extension_key not in manifest:
            raise ImportError("Kepler's provider manifest does not identify the NajaEDA extension")
        for relative, expected in manifest.items():
            path = PurePosixPath(relative)
            if (path.is_absolute() or ".." in path.parts or "\\" in relative
                    or not path.parts or path.parts[0] not in ("najaeda", "najaeda.libs")
                    or not isinstance(expected, str) or not re.fullmatch(r"[a-f0-9]{64}", expected)):
                raise ImportError("Kepler contains an invalid NajaEDA provider manifest entry")
            installed = root.parent.joinpath(*path.parts).resolve(strict=True)
            installed.relative_to(root.parent)
            actual = hashlib.sha256(installed.read_bytes()).hexdigest()
            if actual != expected:
                raise ImportError(
                    f"NajaEDA runtime differs from the wheel Kepler was built against: {relative}. "
                    "Install matching Kepler Formal and NajaEDA wheels.")
    except (AttributeError, OSError, TypeError, ValueError) as error:
        raise ImportError(f"Cannot validate Kepler's NajaEDA provider: {error}") from error
