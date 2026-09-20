# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Check the built extension's imports, including stripped wheel binaries."""

from pathlib import Path
import shutil
import subprocess
import sys
import unittest

from kepler_formal import _native


def _import_command(module_path):
    if sys.platform == "win32":
        # llvm-readobj accompanies the clang-cl wheel toolchain, but some
        # environments expose only clang-cl itself on PATH.
        reader = shutil.which("llvm-readobj")
        compiler = shutil.which("clang-cl")
        if reader is None and compiler is not None:
            sibling = Path(compiler).with_name("llvm-readobj.exe")
            if sibling.is_file():
                reader = str(sibling)
        if reader is not None:
            return [reader, "--coff-imports", module_path]
        dumpbin = shutil.which("dumpbin")
        if dumpbin is not None:
            return [dumpbin, "/nologo", "/imports", module_path]
        raise unittest.SkipTest("requires llvm-readobj or dumpbin to inspect PE imports")

    if sys.platform == "darwin" or sys.platform.startswith("linux"):
        nm = shutil.which("nm") or shutil.which("llvm-nm")
        if nm is None:
            raise unittest.SkipTest("requires nm to inspect native imports")
        # Linux wheels may have their regular symbol table stripped; imported
        # functions are still present in the dynamic symbol table.
        flags = ["-u"] if sys.platform == "darwin" else ["-D", "-u"]
        return [nm, *flags, module_path]
    raise unittest.SkipTest("native import inspection is unavailable on this platform")


class NativeDependenciesTest(unittest.TestCase):
    def test_extension_does_not_import_the_verilog_parser(self):
        completed = subprocess.run(
            _import_command(_native.__file__),
            capture_output=True,
            text=True,
            timeout=30,
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        imports = completed.stdout
        # Confirm we inspected actual C++ imports rather than accepting an
        # empty or unsupported tool response as evidence of independence.
        self.assertIn("NLUniverse", imports)
        parser_imports = [
            line.strip()
            for line in imports.splitlines()
            if "VerilogConstructor" in line or "SNLVRLConstructor" in line
        ]
        self.assertEqual([], parser_imports, "KF Python must leave parsing to NajaEDA")


if __name__ == "__main__":
    unittest.main()
