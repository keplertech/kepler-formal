# SEC regression summary

After each `regress-sec` workflow finishes, **SEC regression summary** generates
a Markdown table in its GitHub Actions run summary and saves it as the
`sec-regression-summary` artifact. This includes unsuccessful or cancelled runs,
and covers the regular, SV2V, and SV2V negative matrix jobs.

Each row shows the case, flow (engine and encoding), CLI runtime, SEC verdict,
checked-output coverage (percentage and covered/total counts), and CI status.
Runtime measures the `kepler-formal` process, excluding checkout, build, and the
regression helper's heartbeat wait. If that measurement is missing, an available
GitHub SEC-step duration is used and explicitly labeled; it includes step setup
and validation. CI success means the configured expectation
was met; it does not necessarily mean equivalence was proved. Checked-output
coverage is not a proven-output percentage.

Unavailable measurements are shown as missing, not zero. Skipped, failed, and
timed-out jobs remain visible even when no result artifact was uploaded. Reruns
use the latest attempt for each job, retaining successful jobs from earlier
attempts. Raw result artifacts are kept for seven days; the table for thirty.
If the runtime build fails before GitHub expands the matrices, the report shows
the build failure and explains that no SEC case results are available.

The follow-up uses GitHub's `workflow_run` trigger, so it becomes automatic once
the workflow file is on the repository's default branch. It uses read-only
repository permissions and treats downloaded results as data, never executable
code.

To test the report generator locally:

```sh
python3 -m unittest discover -s test/sec -p test_sec_regress_summary.py
```
