# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

import sys

from .api import run_cli


def main() -> int:
    return run_cli(sys.argv[1:]).exit_code


if __name__ == "__main__":
    raise SystemExit(main())
