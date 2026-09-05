from __future__ import annotations

import hashlib
import io
import json
import os
import tarfile
import tempfile
import unittest
from pathlib import Path

from kepler_formal_mcp.resolver import BinaryResolver, BinarySetupError


class FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()


def executable(path: Path, text: str = "#!/bin/sh\nexit 0\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(0o755)
    return path


def release_archive() -> bytes:
    payload = b"#!/bin/sh\nexit 0\n"
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as tar:
        info = tarfile.TarInfo("kepler-formal-1.0.0-linux-x86_64/kepler-formal")
        info.mode = 0o755
        info.size = len(payload)
        tar.addfile(info, io.BytesIO(payload))
    return output.getvalue()


class BinaryResolverTests(unittest.TestCase):
    def test_path_precedes_cache_and_network(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path_binary = executable(root / "path" / "kepler-formal")
            executable(root / "cache" / "bin" / "kepler-formal")
            resolver = BinaryResolver(
                root / "cache",
                which=lambda _: str(path_binary),
                urlopen=lambda *_args, **_kwargs: self.fail("network used"),
            )
            result = resolver.resolve()
            self.assertEqual(result.source, "path")
            self.assertEqual(result.path, path_binary.resolve())

    def test_cache_precedes_network(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cached = executable(root / "cache" / "bin" / "kepler-formal")
            resolver = BinaryResolver(
                root / "cache",
                which=lambda _: None,
                urlopen=lambda *_args, **_kwargs: self.fail("network used"),
            )
            result = resolver.resolve()
            self.assertEqual(result.source, "cache")
            self.assertEqual(result.path, cached.resolve())

    def test_download_verifies_checksum_and_is_reused_from_cache(self):
        archive = release_archive()
        digest = hashlib.sha256(archive).hexdigest()
        archive_name = "kepler-formal-1.0.0-linux-x86_64.tar.gz"
        release = {
            "tag_name": "v1.0.0",
            "assets": [
                {"name": archive_name, "browser_download_url": "memory:archive"},
                {
                    "name": f"{archive_name}.sha256",
                    "browser_download_url": "memory:checksum",
                },
            ],
        }
        responses = {
            "memory:api": json.dumps(release).encode(),
            "memory:archive": archive,
            "memory:checksum": f"{digest}  {archive_name}\n".encode(),
        }

        def urlopen(request, timeout):
            return FakeResponse(responses[request.full_url])

        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary) / "cache"
            resolver = BinaryResolver(
                cache,
                release_api="memory:api",
                which=lambda _: None,
                urlopen=urlopen,
                system_name="Linux",
                machine="x86_64",
            )
            result = resolver.resolve()
            self.assertEqual(result.source, "download")
            self.assertEqual(result.release_tag, "v1.0.0")
            self.assertTrue(os.access(result.path, os.X_OK))

            cached = BinaryResolver(
                cache,
                which=lambda _: None,
                urlopen=lambda *_args, **_kwargs: self.fail("network used"),
                system_name="Linux",
                machine="x86_64",
            ).resolve()
            self.assertEqual(cached.source, "cache")
            self.assertEqual(cached.path, result.path)

    def test_checksum_mismatch_is_a_specific_setup_error(self):
        archive = release_archive()
        archive_name = "kepler-formal-1.0.0-linux-x86_64.tar.gz"
        release = {
            "tag_name": "v1.0.0",
            "assets": [
                {"name": archive_name, "browser_download_url": "memory:archive"},
                {
                    "name": f"{archive_name}.sha256",
                    "browser_download_url": "memory:checksum",
                },
            ],
        }
        responses = {
            "memory:api": json.dumps(release).encode(),
            "memory:archive": archive,
            "memory:checksum": ("0" * 64).encode(),
        }

        def urlopen(request, timeout):
            return FakeResponse(responses[request.full_url])

        with tempfile.TemporaryDirectory() as temporary:
            resolver = BinaryResolver(
                Path(temporary) / "cache",
                release_api="memory:api",
                which=lambda _: None,
                urlopen=urlopen,
                system_name="Linux",
                machine="x86_64",
            )
            with self.assertRaisesRegex(BinarySetupError, "SHA-256 mismatch"):
                resolver.resolve()


if __name__ == "__main__":
    unittest.main()
