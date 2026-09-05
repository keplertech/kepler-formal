"""Resolve or install the native kepler-formal release binary."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import stat
import tarfile
import tempfile
import urllib.request
import zipfile
from collections.abc import Callable
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO

from . import __version__
from .models import BinaryResolution


DEFAULT_RELEASE_API = (
    "https://api.github.com/repos/keplertech/kepler-formal/releases/latest"
)
DEFAULT_CACHE_DIR = Path.home() / ".cache" / "kepler-formal"
DOWNLOAD_TIMEOUT_SECONDS = 30


class BinarySetupError(RuntimeError):
    """The checker could not be found or installed."""


class BinaryResolver:
    """Resolve PATH, then the local cache, then a checksummed GitHub release."""

    def __init__(
        self,
        cache_dir: Path | None = None,
        *,
        release_api: str = DEFAULT_RELEASE_API,
        which: Callable[[str], str | None] = shutil.which,
        urlopen: Callable[..., BinaryIO] = urllib.request.urlopen,
        system_name: str | None = None,
        machine: str | None = None,
    ) -> None:
        self.cache_dir = (cache_dir or DEFAULT_CACHE_DIR).expanduser().resolve()
        self.release_api = release_api
        self._which = which
        self._urlopen = urlopen
        self._system_name = system_name or platform.system()
        self._machine = machine or platform.machine()

    def resolve(self) -> BinaryResolution:
        path_candidate = self._which("kepler-formal")
        if path_candidate:
            path = Path(path_candidate).expanduser().resolve()
            if self._is_executable(path):
                return BinaryResolution(path=path, source="path")

        cached = self._find_cached()
        if cached is not None:
            return BinaryResolution(path=cached, source="cache")

        try:
            return self._download_release()
        except Exception as exc:
            detail = str(exc) or exc.__class__.__name__
            raise BinarySetupError(
                "kepler-formal binary resolution failed: no executable named "
                "`kepler-formal` was found on PATH; no cached executable was "
                f"found under {self.cache_dir / 'bin'}; GitHub Releases fallback "
                f"failed: {detail}. Install/build kepler-formal on PATH or place "
                "a compatible release in the cache."
            ) from exc

    def _find_cached(self) -> Path | None:
        binary_name = "kepler-formal.exe" if os.name == "nt" else "kepler-formal"
        bin_root = self.cache_dir / "bin"
        direct = bin_root / binary_name
        if self._is_executable(direct):
            return direct.resolve()
        if not bin_root.is_dir():
            return None

        candidates = [
            path
            for path in bin_root.glob(f"*/*/kepler-formal-*/{binary_name}")
            if self._is_executable(path)
        ]
        candidates.sort(key=lambda path: path.stat().st_mtime, reverse=True)
        return candidates[0].resolve() if candidates else None

    def _download_release(self) -> BinaryResolution:
        platform_key = self._platform_key()
        release = self._read_json(self.release_api)
        tag = self._safe_component(str(release.get("tag_name", "")), "release tag")
        assets = release.get("assets")
        if not isinstance(assets, list):
            raise BinarySetupError("latest GitHub release has no asset list")

        archive_asset = self._select_archive(assets, platform_key)
        archive_name = archive_asset["name"]
        checksum_name = f"{archive_name}.sha256"
        checksum_asset = next(
            (asset for asset in assets if asset.get("name") == checksum_name), None
        )
        if checksum_asset is None:
            raise BinarySetupError(
                f"release asset {archive_name!r} has no checksum asset "
                f"{checksum_name!r}"
            )

        checksum_text = self._read_bytes(checksum_asset["browser_download_url"]).decode(
            "utf-8", errors="replace"
        )
        expected = checksum_text.split(maxsplit=1)[0].lower() if checksum_text else ""
        if not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise BinarySetupError(f"invalid SHA-256 checksum file {checksum_name!r}")

        bin_root = self.cache_dir / "bin"
        staging_root = self.cache_dir / ".staging"
        bin_root.mkdir(parents=True, exist_ok=True)
        staging_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=staging_root) as temporary:
            temp_root = Path(temporary)
            archive_path = temp_root / archive_name
            actual = self._download_to_file(
                archive_asset["browser_download_url"], archive_path
            )
            if actual != expected:
                raise BinarySetupError(
                    f"SHA-256 mismatch for {archive_name}: expected {expected}, "
                    f"downloaded {actual}"
                )

            extract_root = temp_root / "extract"
            extract_root.mkdir()
            self._extract_archive(archive_path, extract_root)
            distribution_root, executable = self._find_extracted_binary(extract_root)
            destination = bin_root / tag / platform_key / distribution_root.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                os.replace(distribution_root, destination)
            cached_executable = destination / executable.relative_to(distribution_root)
            cached_executable.chmod(cached_executable.stat().st_mode | stat.S_IXUSR)
            if not self._is_executable(cached_executable):
                raise BinarySetupError(
                    f"downloaded release did not yield an executable at {cached_executable}"
                )
            return BinaryResolution(
                path=cached_executable.resolve(), source="download", release_tag=tag
            )

    def _select_archive(
        self, assets: list[dict[str, Any]], platform_key: str
    ) -> dict[str, Any]:
        suffixes = (f"-{platform_key}.tar.gz", f"-{platform_key}.zip")
        matches = [
            asset
            for asset in assets
            if isinstance(asset, dict)
            and isinstance(asset.get("name"), str)
            and asset["name"].startswith("kepler-formal-")
            and asset["name"].endswith(suffixes)
            and isinstance(asset.get("browser_download_url"), str)
        ]
        if len(matches) != 1:
            available = ", ".join(
                str(asset.get("name")) for asset in assets if asset.get("name")
            )
            raise BinarySetupError(
                f"no unique release build exists for platform {platform_key!r}; "
                f"available assets: {available or 'none'}"
            )
        return matches[0]

    def _platform_key(self) -> str:
        systems = {"linux": "linux", "darwin": "macos", "windows": "windows"}
        machines = {
            "x86_64": "x86_64",
            "amd64": "x86_64",
            "aarch64": "aarch64",
            "arm64": "arm64" if self._system_name.lower() == "darwin" else "aarch64",
        }
        system = systems.get(self._system_name.lower())
        machine = machines.get(self._machine.lower())
        if system is None or machine is None:
            raise BinarySetupError(
                f"unsupported platform {self._system_name}/{self._machine}"
            )
        return f"{system}-{machine}"

    def _read_json(self, url: str) -> dict[str, Any]:
        value = json.loads(self._read_bytes(url).decode("utf-8"))
        if not isinstance(value, dict):
            raise BinarySetupError("GitHub release response was not a JSON object")
        return value

    def _read_bytes(self, url: str) -> bytes:
        request = urllib.request.Request(
            url,
            headers={
                "Accept": "application/vnd.github+json",
                "User-Agent": f"kepler-formal-mcp/{__version__}",
                "X-GitHub-Api-Version": "2022-11-28",
            },
        )
        with self._urlopen(request, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response:
            return response.read()

    def _download_to_file(self, url: str, destination: Path) -> str:
        request = urllib.request.Request(
            url, headers={"User-Agent": f"kepler-formal-mcp/{__version__}"}
        )
        digest = hashlib.sha256()
        with self._urlopen(request, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response:
            with destination.open("wb") as output:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)
                    digest.update(chunk)
        return digest.hexdigest()

    def _extract_archive(self, archive: Path, destination: Path) -> None:
        if archive.name.endswith(".tar.gz"):
            with tarfile.open(archive, "r:gz") as tar:
                for member in tar.getmembers():
                    self._validate_archive_name(member.name)
                    if member.issym() or member.islnk() or not (
                        member.isfile() or member.isdir()
                    ):
                        raise BinarySetupError(
                            f"unsafe member {member.name!r} in release archive"
                        )
                # Names and member types are validated above. Avoid the newer
                # filter= API so the add-on remains compatible with Python 3.10.
                tar.extractall(destination)
            return
        if archive.name.endswith(".zip"):
            with zipfile.ZipFile(archive) as zipped:
                for member in zipped.infolist():
                    self._validate_archive_name(member.filename)
                    unix_mode = member.external_attr >> 16
                    if stat.S_ISLNK(unix_mode):
                        raise BinarySetupError(
                            f"unsafe member {member.filename!r} in release archive"
                        )
                zipped.extractall(destination)
            return
        raise BinarySetupError(f"unsupported release archive {archive.name!r}")

    @staticmethod
    def _validate_archive_name(name: str) -> None:
        path = PurePosixPath(name)
        if (
            path.is_absolute()
            or ".." in path.parts
            or "\\" in name
            or (path.parts and ":" in path.parts[0])
        ):
            raise BinarySetupError(f"unsafe path {name!r} in release archive")

    @staticmethod
    def _find_extracted_binary(extract_root: Path) -> tuple[Path, Path]:
        names = ("kepler-formal", "kepler-formal.exe")
        for distribution_root in extract_root.iterdir():
            if not distribution_root.is_dir():
                continue
            for name in names:
                executable = distribution_root / name
                if executable.is_file():
                    return distribution_root, executable
        raise BinarySetupError(
            "release archive contains no top-level kepler-formal launcher"
        )

    @staticmethod
    def _safe_component(value: str, label: str) -> str:
        if not value or not re.fullmatch(r"[A-Za-z0-9._-]+", value):
            raise BinarySetupError(f"invalid {label}: {value!r}")
        return value

    @staticmethod
    def _is_executable(path: Path) -> bool:
        return path.is_file() and (os.name == "nt" or os.access(path, os.X_OK))
