#!/usr/bin/env python3
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

SKIP_RETURN_CODE = 77


def load_manifest(benchmarks_root):
    manifest = benchmarks_root / "manifest.tsv"
    cases = {}
    with manifest.open() as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts))
            cases[row["case"]] = row
    return cases


def yaml_quote(path):
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def uses_supported_c2rtl_subset(path):
    return "metron/metron_tools.h" in path.read_text(errors="ignore")


def run_case(args):
    env = os.environ.copy()
    benchmarks_root = Path(args.benchmarks_root).resolve()

    cases = load_manifest(benchmarks_root)
    if args.case not in cases:
        print(f"Unknown benchmark case {args.case!r}", file=sys.stderr)
        return 2

    row = cases[args.case]
    if row["status"] != "complete":
        print(f"Benchmark case {args.case} is not runnable: {row['note']}", file=sys.stderr)
        return SKIP_RETURN_CODE

    case_dir = benchmarks_root / "cases" / args.case
    c_path = case_dir / row["c_source"]
    rtl_path = case_dir / row["rtl_source"]
    if not c_path.is_file() or not rtl_path.is_file():
        print(f"Benchmark fixture is incomplete for {args.case}", file=sys.stderr)
        print(f"  C:   {c_path}", file=sys.stderr)
        print(f"  RTL: {rtl_path}", file=sys.stderr)
        return 2
    if not uses_supported_c2rtl_subset(c_path):
        print(
            f"Skipping {args.case}: upstream verilog-c source is ANSI C verifier "
            "code outside the currently supported synthesizable C2RTL subset",
            file=sys.stderr,
        )
        return SKIP_RETURN_CODE

    with tempfile.TemporaryDirectory(prefix="kepler_verilog_c_") as tmp:
        tmp_dir = Path(tmp)
        c_work = tmp_dir / "c2rtl"
        cfg_path = tmp_dir / "config.yaml"
        cfg_path.write_text(
            "verification: lec\n"
            "design1:\n"
            "  format: c\n"
            f"  input_paths: [\"{yaml_quote(c_path)}\"]\n"
            f"  top: {row['c_top']}\n"
            f"  work_dir: \"{yaml_quote(c_work)}\"\n"
            "design2:\n"
            "  format: systemverilog\n"
            f"  input_paths: [\"{yaml_quote(rtl_path)}\"]\n"
            f"  top: {row['rtl_top']}\n"
        )

        cmd = [str(Path(args.kepler_formal).resolve()), "--config", str(cfg_path)]
        print(f"Running {args.case}: {' '.join(cmd)}")
        completed = subprocess.run(
            cmd,
            cwd=tmp_dir,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=args.timeout,
        )
        if completed.stdout:
            print(completed.stdout, end="")
        if completed.stderr:
            print(completed.stderr, end="", file=sys.stderr)
        if completed.returncode != 0:
            print(
                f"Benchmark {args.case} failed: kepler-formal exited {completed.returncode}",
                file=sys.stderr,
            )
            print(f"Upstream safety expectation: {row['upstream_expected_safety']}", file=sys.stderr)
        return completed.returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kepler-formal", required=True)
    parser.add_argument("--benchmarks-root", required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    try:
        return run_case(args)
    except subprocess.TimeoutExpired as exc:
        print(f"Benchmark {args.case} timed out after {exc.timeout}s", file=sys.stderr)
        return 124


if __name__ == "__main__":
    sys.exit(main())
