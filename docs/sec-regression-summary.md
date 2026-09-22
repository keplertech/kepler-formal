# SEC regression summary

The last job of every `regress-sec` run, **SEC regression summary**, generates a
Markdown table in that run's GitHub Actions summary and saves it as the
`sec-regression-summary` artifact. Because it is part of the run it reports on,
each pull request gets the table for its own regression: it is listed in the
pull request's checks and opens from there, the same way as for a push to
`main`. This includes unsuccessful or cancelled runs, and covers the regular,
SV2V, and SV2V negative matrix jobs.

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

## Comparison with main

The artifact also holds the rows as `sec-regression-summary.json`. Each run
looks up the newest completed `regress-sec` push run on `main` that still has
this artifact and compares itself with it:

- **Changes vs main** lists only the rows whose SEC result, output coverage, or
  CI status differ, shown as `main → this run`, plus cases that exist on one
  side only. It also gives the total runtime over the rows measured in both.
- **All results** keeps the usual columns and adds the main runtime and the
  relative runtime change for every row.

Only `kepler-formal` wall times are compared; a labeled SEC-step fallback on
either side is shown but not compared. Runtimes come from different runners, so
small differences are noise, and pull-request jobs that run at a bounded smoke
depth are not comparable with their full-depth run on `main`. If no baseline is
available, for example before the first `main` run with this job or after the
thirty-day retention, the report says so and shows the plain table.

The job uses read-only repository permissions and treats downloaded results
and the baseline as data, never executable code.

To test the report generator locally:

```sh
python3 -m unittest discover -s test/sec -p test_sec_regress_summary.py
```
