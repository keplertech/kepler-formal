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
FIELDS = ("case", "flow", "runtime", "result", "coverage", "status")
RUNTIME_SECONDS = re.compile(r"^(\d+(?:\.\d+)?) s")
STEP_RUNTIME_SUFFIX = " (SEC step)"
COMMIT = re.compile(r"^[0-9a-f]{7,40}$")


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
                return f"{seconds:.2f} s{STEP_RUNTIME_SUFFIX}"
        except (KeyError, TypeError, ValueError, AttributeError):
            pass
    return "—"


def collect(jobs, artifacts, run_attempt=None):
    """Return the selected jobs and one sorted row of FIELDS per matrix job."""
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
    return selected, sorted(rows)


def load_baseline(path):
    """Read another run's JSON rows as data: plain strings, keyed by case and flow."""
    try:
        rows = json.loads(path.read_text(encoding="utf-8"))["rows"]
        baseline = {}
        for row in rows:
            values = tuple(row[field] for field in FIELDS)
            if not all(isinstance(value, str) for value in values):
                return None
            baseline[values[:2]] = values
        return baseline
    except (OSError, UnicodeError, ValueError, KeyError, TypeError):
        return None


def cli_seconds(runtime):
    """Only kepler-formal wall time is comparable; SEC-step fallbacks are not."""
    match = RUNTIME_SECONDS.match(runtime)
    if match is None or runtime.endswith(STEP_RUNTIME_SUFFIX):
        return None
    return float(match[1])


def runtime_change(now, before):
    if now is None or before is None or before <= 0:
        return "—"
    return f"{(now - before) / before * 100:+.1f}%"


def compare(rows, baseline, baseline_url=None, baseline_sha=None):
    """Report rows whose result, coverage, or CI status differ from the baseline."""
    source = "latest main run"
    if baseline_url and urlsplit(baseline_url).scheme in ("http", "https"):
        source = f"[{source}]({quote(baseline_url, safe=':/?=&%#')})"
    if baseline_sha and COMMIT.match(baseline_sha):
        source += f" at `{baseline_sha[:7]}`"
    lines = ["## Changes vs main", "", f"Baseline: {source}.", ""]

    current = {row[:2]: row for row in rows}
    shared = [key for key in current if key in baseline]
    timed = [(cli_seconds(current[key][2]), cli_seconds(baseline[key][2])) for key in shared]
    timed = [(now, before) for now, before in timed if now is not None and before is not None]
    if timed:
        now, before = sum(t[0] for t in timed), sum(t[1] for t in timed)
        lines.extend([
            f"Total kepler-formal runtime over {len(timed)} rows measured in both runs: "
            f"{now:.2f} s vs {before:.2f} s on main ({runtime_change(now, before)}).", ""])

    def cell(before, now):
        return escape(now) if before == now else f"{escape(before)} → {escape(now)}"

    changes = []
    for key in sorted(set(current) | set(baseline)):
        now, before = current.get(key), baseline.get(key)
        if now is None:
            change, cells = "Not in this run", [escape(value) for value in before[3:]]
        elif before is None:
            change, cells = "Not in main", [escape(value) for value in now[3:]]
        elif now[3:] != before[3:]:
            change, cells = "Changed", [cell(b, n) for b, n in zip(before[3:], now[3:])]
        else:
            continue
        changes.append((escape(key[0]), escape(key[1]), change, *cells))
    if changes:
        lines.extend([
            "| Case | Flow | Change | SEC result | Output coverage | CI status |",
            "| --- | --- | --- | --- | ---: | --- |",
        ])
        lines.extend("| " + " | ".join(row) + " |" for row in changes)
    else:
        lines.append("No SEC result, output coverage, or CI status differences from main "
                     f"across {len(shared)} shared rows.")
    return lines + [""]


def summarize(jobs, artifacts, run_url=None, run_attempt=None,
              baseline=None, baseline_url=None, baseline_sha=None):
    selected, rows = collect(jobs, artifacts, run_attempt)

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
    if rows and baseline is not None:
        lines.extend(compare(rows, baseline, baseline_url, baseline_sha))
        lines.extend([
            "## All results", "",
            "| Case | Flow | Runtime | SEC result | Output coverage | CI status "
            "| Main runtime | Runtime vs main |",
            "| --- | --- | ---: | --- | ---: | --- | ---: | ---: |",
        ])
        for row in rows:
            before = baseline.get(row[:2])
            main_runtime = before[2] if before else "—"
            change = runtime_change(cli_seconds(row[2]), cli_seconds(main_runtime))
            lines.append("| " + " | ".join(map(escape, (*row, main_runtime, change))) + " |")
    elif rows:
        lines.extend([
            "| Case | Flow | Runtime | SEC result | Output coverage | CI status |",
            "| --- | --- | ---: | --- | ---: | --- |",
        ])
        lines.extend("| " + " | ".join(map(escape, row)) + " |" for row in rows)
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
    parser.add_argument("--json-output", type=Path,
                        help="also write the rows as JSON, for use as a later --baseline")
    parser.add_argument("--baseline", type=Path,
                        help="--json-output of the main run to compare against")
    parser.add_argument("--baseline-url")
    parser.add_argument("--baseline-sha")
    args = parser.parse_args()
    jobs = json.loads(args.jobs.read_text(encoding="utf-8"))
    baseline = load_baseline(args.baseline) if args.baseline else None
    report = summarize(jobs, args.artifacts, args.run_url, args.run_attempt,
                       baseline, args.baseline_url, args.baseline_sha)
    if args.baseline and baseline is None:
        report += "\nNo main baseline summary was available, so this run is not compared with main.\n"
    args.output.write_text(report, encoding="utf-8")
    if args.json_output:
        rows = collect(jobs, args.artifacts, args.run_attempt)[1]
        args.json_output.write_text(
            json.dumps({"version": 1, "rows": [dict(zip(FIELDS, row)) for row in rows]}),
            encoding="utf-8")


if __name__ == "__main__":
    main()
