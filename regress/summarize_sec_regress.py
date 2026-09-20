#!/usr/bin/env python3
# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0
"""Summarize SEC matrix jobs and their text artifacts without executing them."""

import argparse
from datetime import datetime
import html
import json
import math
from pathlib import Path
import re
from urllib.parse import quote, urlsplit


MATRIX_JOB = re.compile(r"^(?:sv2v (?:negative )?)?(imc|ki|pdr)-[^/]+ / .+$", re.DOTALL)
SEC_STEP = re.compile(r"^Run (?:SV2V )?SEC (?:negative )?regression$")
ANSI = re.compile(r"\x1b\[[0-9;]*m")
COVERAGE = re.compile(
    r"SEC (?:checked-output|output) coverage:\s*(\d+(?:\.\d+)?%)"
    r"(?:\s*\((\d+)\s*/\s*(\d+)[^)]*\))?"
)
VERDICTS = (
    ("SEC partially proved equivalence", "Partially proved"),
    ("SEC proved equivalence", "Proved"),
    ("SEC found a counterexample", "Counterexample"),
    ("SEC was inconclusive", "Inconclusive"),
    ("SEC cannot run on this design pair", "Unsupported"),
    ("SEC BTOR2 exported", "Exported; proof not run"),
)
ERROR = re.compile(
    r"\[(?:error|critical)\]|(?:SEC (?:compact )?)?[Ww]orkflow failed:|\bError:"
)


def escape(value):
    """Keep artifact/job text literal inside Markdown table cells."""
    value = re.sub(r"[\x00-\x1f\x7f]", " ", str(value))
    if len(value) > 240:
        value = value[:237] + "..."
    value = html.escape(value, quote=True)
    for char in "\\`*_[]|":
        value = value.replace(char, "\\" + char)
    return value


def latest_jobs(jobs, run_attempt=None):
    latest = {}
    for job in jobs:
        attempt = job.get("run_attempt", 1)
        if run_attempt is not None and attempt > run_attempt:
            continue
        previous = latest.get(job["name"])
        if previous is None or (attempt, job["id"]) > (
            previous.get("run_attempt", 1), previous["id"]
        ):
            latest[job["name"]] = job
    return latest.values()


def artifact_file(root, job, engine, extension):
    name = job["name"].replace(" / ", "-").replace(" ", "-")
    directory = root / f"sec-regress-results-{job.get('run_attempt', 1)}-{name}"
    path = directory / f"{engine}.{extension}"
    # Artifact contents and job names are untrusted. Do not follow paths or
    # symlinks outside the expected artifact directory.
    try:
        if directory.resolve().parent != root.resolve():
            return None
        if path.resolve().parent != directory.resolve() or not path.is_file():
            return None
    except (OSError, ValueError):
        return None
    return path


def read_result(path):
    if path is None:
        return "No result (log unavailable)", "—"
    result, coverage = "No final SEC result", "—"
    try:
        with path.open(encoding="utf-8", errors="replace") as stream:
            for line in stream:
                line = ANSI.sub("", line)
                match = COVERAGE.search(line)
                if match:
                    coverage = match[1]
                    if match[2] is not None:
                        coverage += f" ({match[2]}/{match[3]})"
                verdict = next((label for text, label in VERDICTS if text in line), None)
                if verdict is not None:
                    result = verdict
                elif ERROR.search(line):
                    result = "Error"
    except OSError:
        return "No result (log unreadable)", "—"
    return result, coverage


def read_runtime(path):
    if path is None:
        return "—"
    try:
        # Bash can include a signal diagnostic before the final timing line.
        for line in reversed(path.read_text(encoding="utf-8").splitlines()):
            if re.fullmatch(r"\d+(?:\.\d+)?", line.strip()):
                seconds = float(line)
                if math.isfinite(seconds):
                    return f"{seconds:.2f} s"
    except (OSError, UnicodeError, ValueError):
        pass
    return "—"


def step_runtime(step):
    if step is not None:
        try:
            start = datetime.fromisoformat(step["started_at"].replace("Z", "+00:00"))
            end = datetime.fromisoformat(step["completed_at"].replace("Z", "+00:00"))
            seconds = (end - start).total_seconds()
            if seconds >= 0:
                return f"{seconds:.2f} s (SEC step)"
        except (KeyError, TypeError, ValueError, AttributeError):
            pass
    return "—"


def summarize(jobs, artifacts, run_url=None, run_attempt=None):
    selected = list(latest_jobs(jobs, run_attempt))
    rows = []
    for job in selected:
        match = MATRIX_JOB.match(job["name"])
        if match is None:
            continue
        flow, case = job["name"].split(" / ", 1)
        engine = {"imc": "imc", "ki": "k_induction", "pdr": "pdr"}[match[1]]
        result, coverage = read_result(artifact_file(artifacts, job, engine, "stdout"))
        runtime = read_runtime(artifact_file(artifacts, job, engine, "seconds"))
        step = next((s for s in job.get("steps", []) if SEC_STEP.match(s["name"])), None)
        if runtime == "—":
            runtime = step_runtime(step)
        if job.get("conclusion") == "skipped" or (
            step is not None and step.get("conclusion") == "skipped"
        ):
            result, coverage, runtime = "Skipped (SEC not run)", "—", "—"
        status = job.get("conclusion") or job.get("status", "unknown")
        rows.append((case, flow, runtime, result, coverage, status))

    lines = ["# SEC regression results", ""]
    if run_url and urlsplit(run_url).scheme in ("http", "https"):
        lines.extend([f"[Source workflow run]({quote(run_url, safe=':/?=&%#')})", ""])
    builds = [j for j in selected if j["name"] == "build SEC runtime"]
    if builds:
        status = builds[0].get("conclusion") or builds[0].get("status", "unknown")
        lines.extend([f"Runtime build: {escape(status)}.", ""])
    lines.extend([
        "Runtime is kepler-formal wall time; labeled SEC-step fallbacks include setup/validation.",
        "Output coverage is checked outputs, not proven outputs.",
        "CI status reports whether the regression expectation passed; it is not the SEC verdict.",
        "Each row uses the latest job attempt; missing measurements are shown as —.",
        "",
    ])
    if rows:
        lines.extend([
            "| Case | Flow | Runtime | SEC result | Output coverage | CI status |",
            "| --- | --- | ---: | --- | ---: | --- |",
        ])
        lines.extend("| " + " | ".join(map(escape, row)) + " |" for row in sorted(rows))
    else:
        lines.append("No SEC matrix jobs were available; no proof results can be reported.")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", required=True, type=Path)
    parser.add_argument("--artifacts", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run-url")
    parser.add_argument("--run-attempt", type=int)
    args = parser.parse_args()
    jobs = json.loads(args.jobs.read_text(encoding="utf-8"))
    args.output.write_text(
        summarize(jobs, args.artifacts, args.run_url, args.run_attempt), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
