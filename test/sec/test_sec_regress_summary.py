# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[2] / "regress/summarize_sec_regress.py"
SPEC = importlib.util.spec_from_file_location("sec_summary", SCRIPT)
SUMMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMMARY)


class SecRegressionSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def job(self, name="pdr-dual-rail / example", conclusion="success", attempt=1, job_id=1):
        return {"name": name, "conclusion": conclusion, "run_attempt": attempt, "id": job_id}

    def artifact(self, job, stdout, seconds="1.234", engine="pdr"):
        name = job["name"].replace(" / ", "-").replace(" ", "-")
        directory = self.root / f"sec-regress-results-{job['run_attempt']}-{name}"
        directory.mkdir(parents=True, exist_ok=True)
        (directory / f"{engine}.stdout").write_text(stdout, encoding="utf-8")
        if seconds is not None:
            (directory / f"{engine}.seconds").write_text(seconds, encoding="utf-8")
        return directory

    def test_all_sec_verdicts_are_independent_of_ci_success(self):
        messages = {
            "SEC proved equivalence at k = 1.": "Proved",
            "SEC partially proved equivalence at k = 1: 1/2 outputs proved.": "Partially proved",
            "SEC found a counterexample at k = 1.": "Counterexample",
            "SEC was inconclusive up to max_k = 1.": "Inconclusive",
            "[critical] SEC cannot run on this design pair: no aligned outputs": "Unsupported",
            "[error] SEC compact workflow failed: missing input": "Error",
            "SEC BTOR2 exported to model.btor; proof not run.": "Exported; proof not run",
        }
        for message, verdict in messages.items():
            with self.subTest(verdict=verdict):
                job = self.job()
                self.artifact(job, message)
                report = SUMMARY.summarize([job], self.root)
                self.assertIn(f"| 1.23 s | {verdict} | — | success |", report)

    def test_final_coverage_and_verdict_win(self):
        job = self.job()
        self.artifact(job, """SEC output coverage: 25.00% (1/4 outputs).
SEC was inconclusive up to max_k = 1.
\x1b[32mSEC checked-output coverage: 75.00% (3/4 covered/existing outputs).\x1b[0m
SEC partially proved equivalence at k = 2: 2/4 outputs proved.
""")
        report = SUMMARY.summarize([job], self.root)
        self.assertIn("| Partially proved | 75.00% (3/4) | success |", report)
        self.assertNotIn("25.00%", report)

    def test_legacy_coverage(self):
        job = self.job()
        self.artifact(job, "SEC output coverage: 100.00% (4/4 outputs).\nSEC proved equivalence")
        self.assertIn("100.00% (4/4)", SUMMARY.summarize([job], self.root))

    def test_missing_failed_timeout_and_cancelled_jobs_are_visible(self):
        for conclusion in ("success", "failure", "timed_out", "cancelled"):
            with self.subTest(conclusion=conclusion):
                report = SUMMARY.summarize([self.job(conclusion=conclusion)], self.root)
                self.assertIn(f"| — | No result (log unavailable) | — | {SUMMARY.escape(conclusion)} |", report)
                self.assertNotIn("| Proved |", report)

    def test_interrupted_log_retains_coverage_but_not_a_proof(self):
        job = self.job(conclusion="timed_out")
        self.artifact(job, "SEC checked-output coverage: 50.00% (2/4 covered/existing outputs).", None)
        self.assertIn("| — | No final SEC result | 50.00% (2/4) | timed\\_out |",
                      SUMMARY.summarize([job], self.root))

    def test_skipped_regression_step_is_not_success(self):
        for name in ("Run SEC regression", "Run SV2V SEC regression", "Run SV2V SEC negative regression"):
            with self.subTest(step=name):
                job = self.job("sv2v negative ki-dual-rail / disabled")
                job["steps"] = [{"name": name, "conclusion": "skipped"}]
                self.artifact(job, "SEC proved equivalence", engine="k_induction")
                self.assertIn("| — | Skipped (SEC not run) | — | success |",
                              SUMMARY.summarize([job], self.root))

    def test_skipped_job(self):
        report = SUMMARY.summarize([self.job(conclusion="skipped")], self.root)
        self.assertIn("Skipped (SEC not run)", report)

    def test_latest_rerun_replaces_previous_result_and_artifact(self):
        old = self.job(attempt=1)
        new = self.job(attempt=2, job_id=3, conclusion="failure")
        future = self.job(attempt=3, job_id=5)
        self.artifact(old, "SEC proved equivalence")
        self.artifact(new, "SEC found a counterexample", "3")
        self.artifact(future, "SEC proved equivalence")
        report = SUMMARY.summarize([future, new, old], self.root, run_attempt=2)
        self.assertEqual(report.count("| example |"), 1)
        self.assertIn("| 3.00 s | Counterexample | — | failure |", report)
        self.assertNotIn("| Proved |", report)
        same_attempt = self.job(attempt=2, job_id=4, conclusion="cancelled")
        selected = list(SUMMARY.latest_jobs([new, old, same_attempt], 2))
        self.assertEqual(selected, [same_attempt])

    def test_rerun_without_artifact_never_uses_older_artifact(self):
        old, new = self.job(), self.job(attempt=2, job_id=2)
        self.artifact(old, "SEC proved equivalence")
        report = SUMMARY.summarize([old, new], self.root)
        self.assertIn("No result (log unavailable)", report)

    def test_all_flow_names_and_deterministic_sorting(self):
        jobs = [self.job("sv2v negative ki-dual-rail / zebra"),
                self.job("sv2v imc-binary / alpha"), self.job("pdr-binary / beta")]
        for job, engine in zip(jobs, ("k_induction", "imc", "pdr")):
            self.artifact(job, "SEC proved equivalence", engine=engine)
        report = SUMMARY.summarize(jobs, self.root)
        self.assertEqual(report, SUMMARY.summarize(list(reversed(jobs)), self.root))
        self.assertEqual(report.count("| Proved |"), 3)
        self.assertLess(report.index("| alpha |"), report.index("| beta |"))

    def test_invalid_and_missing_runtime_are_not_zero(self):
        job = self.job()
        for seconds in (None, "garbage", "nan", "inf", "-1"):
            with self.subTest(seconds=seconds):
                directory = self.artifact(job, "SEC proved equivalence", seconds)
                report = SUMMARY.summarize([job], self.root)
                self.assertIn("| — | Proved |", report)
                (directory / "pdr.seconds").unlink(missing_ok=True)
        self.artifact(job, "SEC proved equivalence", "0")
        self.assertIn("| 0.00 s |", SUMMARY.summarize([job], self.root))

    def test_bash_signal_diagnostic_before_elapsed_time(self):
        job = self.job(conclusion="failure")
        self.artifact(job, "running", "line 321: 1234 Killed kepler-formal\n18.123\n")
        self.assertIn("| 18.12 s | No final SEC result |", SUMMARY.summarize([job], self.root))

    def test_step_runtime_fallback_is_labeled_and_cli_runtime_is_preferred(self):
        job = self.job(conclusion="timed_out")
        job["steps"] = [{"name": "Run SEC regression", "conclusion": "timed_out",
                         "started_at": "2026-09-19T12:00:00Z",
                         "completed_at": "2026-09-19T12:01:05Z"}]
        self.assertIn("| 65.00 s (SEC step) |", SUMMARY.summarize([job], self.root))
        self.artifact(job, "running", "61")
        self.assertIn("| 61.00 s |", SUMMARY.summarize([job], self.root))
        self.assertNotIn("65.00 s", SUMMARY.summarize([job], self.root))

    def test_names_are_escaped_and_artifact_paths_are_confined(self):
        name = "pdr-binary / <img>|[link](url)_*`\\\ncase"
        report = SUMMARY.summarize([self.job(name)], self.root)
        self.assertIn("&lt;img&gt;\\|\\[link\\](url)\\_\\*\\`\\\\ case", report)
        self.assertNotIn("<img>", report)
        job = self.job("pdr-binary / ../../escape")
        self.assertIsNone(SUMMARY.artifact_file(self.root, job, "pdr", "stdout"))
        report = SUMMARY.summarize([self.job("pdr-binary / " + "a" * 1000)], self.root)
        self.assertIn("a" * 237 + "...", report)

    def test_artifact_symlink_is_not_followed(self):
        job = self.job()
        directory = self.artifact(job, "SEC proved equivalence")
        outside = self.root / "outside.txt"
        outside.write_text("SEC proved equivalence")
        (directory / "pdr.stdout").unlink()
        (directory / "pdr.stdout").symlink_to(outside)
        self.assertIn("No result (log unavailable)", SUMMARY.summarize([job], self.root))

    def test_failed_build_and_no_cases(self):
        report = SUMMARY.summarize([self.job("build SEC runtime", conclusion="failure")], self.root)
        self.assertIn("Runtime build: failure.", report)
        self.assertIn("No SEC matrix jobs were available", report)

    def test_cli_writes_report_and_safe_run_link(self):
        jobs_path, output = self.root / "jobs.json", self.root / "report.md"
        jobs_path.write_text(json.dumps([self.job()]), encoding="utf-8")
        subprocess.run([sys.executable, str(SCRIPT), "--jobs", str(jobs_path),
                        "--artifacts", str(self.root), "--output", str(output),
                        "--run-url", "https://github.com/org/repo/actions/runs/1",
                        "--run-attempt", "1"], check=True)
        self.assertIn("[Source workflow run](https://github.com/org/repo/actions/runs/1)",
                      output.read_text(encoding="utf-8"))
        self.assertNotIn("javascript:", SUMMARY.summarize([], self.root, "javascript:alert(1)"))

    def baseline(self, *rows):
        return {row[:2]: row for row in rows}

    def test_baseline_reports_only_changed_new_and_removed_rows(self):
        same, changed, added = (self.job("pdr-dual-rail / same"),
                                self.job("pdr-dual-rail / changed", conclusion="failure"),
                                self.job("pdr-dual-rail / added"))
        self.artifact(same, "SEC proved equivalence", "2")
        self.artifact(changed, "SEC checked-output coverage: 87.50% (7/8 covered/existing outputs).\n"
                               "SEC partially proved equivalence", "30")
        self.artifact(added, "SEC proved equivalence")
        baseline = self.baseline(
            ("same", "pdr-dual-rail", "1.00 s", "Proved", "—", "success"),
            ("changed", "pdr-dual-rail", "10.00 s", "Proved", "100.00% (8/8)", "success"),
            ("removed", "pdr-dual-rail", "5.00 s", "Proved", "—", "success"))
        report = SUMMARY.summarize([same, changed, added], self.root, baseline=baseline,
                                   baseline_url="https://github.com/org/repo/actions/runs/7",
                                   baseline_sha="c7b8fe6db94d605acd6c4c2ae9b349dc272b1106")
        changes, results = report.split("## All results")
        self.assertIn("Baseline: [latest main run](https://github.com/org/repo/actions/runs/7)"
                      " at `c7b8fe6`.", changes)
        self.assertIn("| changed | pdr-dual-rail | Changed | Proved → Partially proved "
                      "| 100.00% (8/8) → 87.50% (7/8) | success → failure |", changes)
        self.assertIn("| added | pdr-dual-rail | Not in main | Proved | — | success |", changes)
        self.assertIn("| removed | pdr-dual-rail | Not in this run | Proved | — | success |", changes)
        self.assertNotIn("| same |", changes)
        self.assertIn("over 2 rows measured in both runs: 32.00 s vs 11.00 s on main (+190.9%)",
                      changes)
        self.assertIn("| same | pdr-dual-rail | 2.00 s | Proved | — | success | 1.00 s | +100.0% |",
                      results)
        self.assertIn("| added | pdr-dual-rail | 1.23 s | Proved | — | success | — | — |", results)

    def test_baseline_without_differences_and_incomparable_runtimes(self):
        job = self.job(conclusion="timed_out")
        job["steps"] = [{"name": "Run SEC regression", "conclusion": "timed_out",
                         "started_at": "2026-09-19T12:00:00Z",
                         "completed_at": "2026-09-19T12:01:05Z"}]
        baseline = self.baseline(("example", "pdr-dual-rail", "0.00 s",
                                  "No result (log unavailable)", "—", "timed_out"))
        report = SUMMARY.summarize([job], self.root, baseline=baseline,
                                   baseline_url="javascript:alert(1)", baseline_sha="`x`")
        self.assertIn("No SEC result, output coverage, or CI status differences from main "
                      "across 1 shared rows.", report)
        self.assertIn("Baseline: latest main run.", report)
        self.assertNotIn("javascript:", report)
        self.assertNotIn("measured in both runs", report)
        self.assertIn("| 65.00 s (SEC step) | No result (log unavailable) | — | timed\\_out "
                      "| 0.00 s | — |", report)

    def test_baseline_text_is_escaped(self):
        job = self.job()
        self.artifact(job, "SEC proved equivalence")
        baseline = self.baseline(("example", "pdr-dual-rail", "<b>", "<img>|x", "—", "success"))
        report = SUMMARY.summarize([job], self.root, baseline=baseline)
        self.assertIn("&lt;img&gt;\\|x → Proved", report)
        self.assertIn("| &lt;b&gt; | — |", report)
        self.assertNotIn("<img>", report)

    def test_cli_json_output_round_trips_as_baseline(self):
        jobs_path = self.root / "jobs.json"
        job = self.job()
        self.artifact(job, "SEC proved equivalence")
        jobs_path.write_text(json.dumps([job]), encoding="utf-8")
        base = [sys.executable, str(SCRIPT), "--jobs", str(jobs_path), "--artifacts", str(self.root)]
        subprocess.run(base + ["--output", str(self.root / "main.md"),
                               "--json-output", str(self.root / "main.json")], check=True)
        self.assertEqual(json.loads((self.root / "main.json").read_text(encoding="utf-8")), {
            "version": 1, "rows": [{"case": "example", "flow": "pdr-dual-rail",
                                    "runtime": "1.23 s", "result": "Proved",
                                    "coverage": "—", "status": "success"}]})
        self.assertNotIn("vs main", (self.root / "main.md").read_text(encoding="utf-8"))
        subprocess.run(base + ["--output", str(self.root / "pr.md"),
                               "--baseline", str(self.root / "main.json")], check=True)
        report = (self.root / "pr.md").read_text(encoding="utf-8")
        self.assertIn("across 1 shared rows", report)
        self.assertIn("| 1.23 s | Proved | — | success | 1.23 s | +0.0% |", report)

    def test_missing_or_malformed_baseline_is_reported_not_trusted(self):
        jobs_path, output = self.root / "jobs.json", self.root / "report.md"
        jobs_path.write_text(json.dumps([self.job()]), encoding="utf-8")
        malformed = self.root / "malformed.json"
        for content in (None, "not json", "[]", '{"rows": [{"case": "example"}]}',
                        json.dumps({"rows": [dict.fromkeys(SUMMARY.FIELDS, 1)]})):
            with self.subTest(content=content):
                if content is not None:
                    malformed.write_text(content, encoding="utf-8")
                self.assertIsNone(SUMMARY.load_baseline(malformed))
                subprocess.run([sys.executable, str(SCRIPT), "--jobs", str(jobs_path),
                                "--artifacts", str(self.root), "--output", str(output),
                                "--baseline", str(malformed)], check=True)
                report = output.read_text(encoding="utf-8")
                self.assertIn("No main baseline summary was available", report)
                self.assertIn("| Case | Flow | Runtime | SEC result | Output coverage | CI status |\n",
                              report)


if __name__ == "__main__":
    unittest.main()
