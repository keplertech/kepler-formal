# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Unpack pinned conda libraries into one prefix without installing Conda."""

import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile


def unpack(archive, output):
    with tempfile.TemporaryDirectory() as temporary, zipfile.ZipFile(archive) as package:
        root = Path(temporary)
        for name in package.namelist():
            if name.endswith(".tar.zst"):
                data = subprocess.run(
                    ["zstd", "--decompress", "--stdout"],
                    input=package.read(name), capture_output=True, check=True,
                ).stdout
                with tarfile.open(fileobj=io.BytesIO(data)) as contents:
                    contents.extractall(root, filter="data")

        metadata = json.loads((root / "info/index.json").read_text())
        paths = json.loads((root / "info/paths.json").read_text())["paths"]
        for entry in paths:
            if "prefix_placeholder" not in entry:
                continue
            path = root / entry["_path"]
            old, new = entry["prefix_placeholder"].encode(), os.fsencode(output)
            data = path.read_bytes()
            if entry["file_mode"] == "text":
                data = data.replace(old, new)
            else:
                # Preserve binary offsets and the suffix after the prefix;
                # padding belongs after the complete NUL-terminated string.
                def replace(match):
                    changed = match[0].replace(old, new)
                    if len(changed) > len(match[0]):
                        raise ValueError(f"Installation prefix too long for {path}")
                    return changed.ljust(len(match[0]), b"\0")

                data = re.sub(re.escape(old) + rb"[^\0]*\0", replace, data)
            path.write_bytes(data)

        licenses = root / "info/licenses"
        if licenses.exists():
            shutil.copytree(licenses, output / "share/licenses" / metadata["name"], dirs_exist_ok=True)
        shutil.rmtree(root / "info")
        shutil.copytree(root, output, symlinks=True, dirs_exist_ok=True)


if __name__ == "__main__":
    destination = Path(sys.argv[1])
    for source in sys.argv[2:]:
        unpack(source, destination)
