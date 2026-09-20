#!/usr/bin/env python3
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Select published-provider metadata in a build checkout, or validate a release.

No build, upload, or version bump is performed. The default development checkout
is left unchanged. Uses only the standard library, including on Python 3.10.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re


DEVELOPMENT_REQUIREMENT = "najaeda==0.7.24.dev0"
PUBLISHED_REQUIREMENT = "najaeda==0.7.24"
CMAKE_OPTION = "KEPLER_USE_PUBLISHED_NAJAEDA"
VERSION_PATTERN = r'KEPLER_VERSION\s*\{\s*"(?P<value>[0-9]+\.[0-9]+\.[0-9]+)"'


def read_version(project: Path) -> str:
    """Read the same source and regex used by pyproject's version provider."""
    source = project / "src/bin/KeplerVersion.h.in"
    matches = list(re.finditer(VERSION_PATTERN, source.read_text(encoding="utf-8")))
    if len(matches) != 1:
        raise ValueError(f"Expected one package version in {source}")
    return matches[0].group("value")


def _section(text: str, name: str) -> re.Match[str]:
    # Match only the checked-in TOML sections we edit. Other metadata,
    # comments, and multiline platform commands are preserved verbatim.
    pattern = rf"(?ms)^\[{re.escape(name)}\][ \t]*\n(?P<body>.*?)(?=^\[|\Z)"
    matches = list(re.finditer(pattern, text))
    if len(matches) != 1:
        raise ValueError(f"Expected one [{name}] section in pyproject.toml")
    return matches[0]


def _requirements(text: str) -> tuple[str, str]:
    result = []
    for section in ("build-system", "project"):
        requirements = re.findall(r'["\'](najaeda[^"\']*)["\']',
                                  _section(text, section).group("body"))
        if len(requirements) != 1:
            raise ValueError(f"Expected one exact NajaEDA pin in [{section}]")
        result.append(requirements[0])
    return tuple(result)


def _cmake_option(text: str) -> str | None:
    body = _section(text, "tool.scikit-build.cmake.define").group("body")
    values = re.findall(rf'(?m)^{CMAKE_OPTION}\s*=\s*["\'](ON|OFF)["\']\s*(?:#.*)?$', body)
    if not values and CMAKE_OPTION not in body:
        return None
    if len(values) != 1:
        raise ValueError(f"Expected one ON/OFF {CMAKE_OPTION} CMake option")
    return values[0]


def prepare_project(project: Path, *, published_najaeda: bool = False) -> None:
    """Change only the two provider pins and CMake option when opted in."""
    metadata = project / "pyproject.toml"
    original = metadata.read_text(encoding="utf-8")
    requirements = _requirements(original)
    option = _cmake_option(original)
    if not published_najaeda:
        if requirements != (DEVELOPMENT_REQUIREMENT,) * 2 or option not in (None, "OFF"):
            raise ValueError("Development builds require development NajaEDA pins and the published option OFF")
        return
    if requirements not in ((DEVELOPMENT_REQUIREMENT,) * 2, (PUBLISHED_REQUIREMENT,) * 2):
        raise ValueError("Build and runtime NajaEDA pins must both select the supported provider")

    prepared = original
    for section in ("build-system", "project"):
        match = _section(prepared, section)
        body = match.group("body").replace(DEVELOPMENT_REQUIREMENT, PUBLISHED_REQUIREMENT)
        prepared = prepared[:match.start("body")] + body + prepared[match.end("body"):]
    match = _section(prepared, "tool.scikit-build.cmake.define")
    body = match.group("body")
    if option is None:
        body = body.rstrip() + f'\n{CMAKE_OPTION} = "ON"\n\n'
    else:
        body = re.sub(rf'(?m)^({CMAKE_OPTION}\s*=\s*["\'])(ON|OFF)(["\'])',
                      r'\g<1>ON\3', body)
    prepared = prepared[:match.start("body")] + body + prepared[match.end("body"):]
    if prepared != original:
        metadata.write_text(prepared, encoding="utf-8")


def validate_release(project: Path, *, repository: str, ref: str,
                     confirmed_version: str) -> str:
    if repository != "keplertech/kepler-formal" or ref != "refs/heads/main":
        raise ValueError("PyPI publishing is only allowed from keplertech/kepler-formal main")
    metadata = (project / "pyproject.toml").read_text(encoding="utf-8")
    if _requirements(metadata) != (PUBLISHED_REQUIREMENT,) * 2:
        raise ValueError("Publishing requires matching published NajaEDA build and runtime pins")
    if _cmake_option(metadata) != "ON":
        raise ValueError(f"Publishing requires {CMAKE_OPTION}=ON")
    version = read_version(project)
    if confirmed_version != version:
        raise ValueError(f"To publish, confirm version {version} in the workflow's version input")
    return version


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--published-najaeda", action="store_true",
                        default=os.environ.get(CMAKE_OPTION) == "1")
    parser.add_argument("--validate-publish", action="store_true")
    parser.add_argument("--version", default=os.environ.get("RELEASE_VERSION", ""))
    args = parser.parse_args()
    try:
        if args.validate_publish:
            version = validate_release(
                args.project, repository=os.environ.get("GITHUB_REPOSITORY", ""),
                ref=os.environ.get("GITHUB_REF", ""), confirmed_version=args.version,
            )
            summary = f"Release request: kepler-formal {version} from {os.environ.get('GITHUB_SHA', 'unknown')}\n"
            print(summary, end="")
            if path := os.environ.get("GITHUB_STEP_SUMMARY"):
                with open(path, "a", encoding="utf-8") as output:
                    output.write(summary)
        else:
            prepare_project(args.project, published_najaeda=args.published_najaeda)
            requirement = PUBLISHED_REQUIREMENT if args.published_najaeda else DEVELOPMENT_REQUIREMENT
            print(f"Python build provider: {requirement}")
    except (ValueError, OSError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
