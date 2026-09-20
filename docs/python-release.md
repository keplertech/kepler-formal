# Publishing the Python package

The **Python wheels** workflow (`.github/workflows/python-wheels.yml`) can
publish `kepler-formal` to PyPI through a manual run on `main`. Publishing is
off by default. Relevant PRs, pushes to `main`, `v*` tags, and manual builds with
`publish` unchecked use the development NajaEDA provider and only produce
tested wheel artifacts.
A manual `publish` request adds published-provider tests in parallel with the
development jobs. Only the published-provider jobs gate publication.

This is separate from [binary releases](releasing.md). Do not create a tag or
run `tools/release.sh` for a Python-only release: `v*` tags also trigger the
standalone binary release workflow. The manual Python workflow creates no
tags or GitHub Releases and does not publish the MCP add-on.

The Python package verifies live NajaEDA netlists with `verify_designs()`.
NajaEDA owns file loading and the caller owns the netlists' lifetime; Kepler
does not delete them after verification. The Python `verify()`, `run_cli()`,
and `python -m kepler_formal` entrypoints are removed. Use NajaEDA to load
designs before calling the Python API, or the standalone binary for file/YAML
workflows.

## One-time setup

Complete these account-side settings before requesting a publication. Merely
declaring an environment in the workflow does not configure its protections.

1. In `keplertech/kepler-formal`, open **Settings → Environments** and create
   an environment named `pypi`.
2. Configure **required reviewers** for release approval. Under deployment
   branches and tags, allow only the branch `main`, with no tag rule. Restrict
   administrator bypass where available. These settings protect the publisher
   independently of the workflow's own branch checks.
3. On PyPI, configure a GitHub **Trusted Publisher**. For a new project, use
   [pending publisher registration](https://pypi.org/manage/account/publishing/).
   If the project already exists and you own it, use its **Manage → Publishing**
   page instead. Enter:

   - PyPI project: `kepler-formal`
   - GitHub owner: `keplertech`
   - Repository: `kepler-formal`
   - Workflow filename: `python-wheels.yml` (not its full path)
   - Environment: `pypi`

No stored PyPI API token is needed. Only the separate publishing job receives
`id-token: write`; build jobs do not. The publisher downloads the current
run's four published-provider wheel artifacts and uploads them without
checking out or building source in the privileged job. Development-provider
artifacts are never uploaded to PyPI.

See the official [PyPI setup instructions](https://docs.pypi.org/trusted-publishers/adding-a-publisher/),
[Trusted Publishing instructions](https://docs.pypi.org/trusted-publishers/using-a-publisher/),
and [GitHub environment protection documentation](https://docs.github.com/en/actions/reference/workflows-and-actions/deployments-and-environments).

## Build without publishing

After the workflow is present on the repository's default branch, open
**Actions → Python wheels → Run workflow**. Leave `publish` unchecked and
`version` empty. These jobs use locally built NajaEDA `0.7.24.dev0` from the
pinned submodule and its shared-runtime SDK. The workflow builds, repairs,
and tests wheels, then stores them as Actions artifacts. Use this to check a
branch before merging.

| Run | Build and test providers |
| --- | --- |
| Relevant pull request, push to `main`, `v*` tag, or manual run without `publish` | Development provider only. |
| Manual run with `publish` | Development provider and published NajaEDA `0.7.24`, in parallel. |

There is no separate provider checkbox. The published-provider jobs run only
as a prerequisite to a requested publication. The [local source
regression](python-regression.md) remains separate and continues building both
packages from the pinned source without wheels or publishing.

Both wheel paths use the reusable `python-wheel-build.yml` workflow to share
their build and test steps. The publishing job remains in `python-wheels.yml`
with the same `pypi` environment and Trusted Publisher configuration.

## Published-provider integration

The published-provider jobs set `KEPLER_USE_PUBLISHED_NAJAEDA=1`
for the wheel helpers, including Linux containers. Before building,
`ci/prepare_python_release.py` updates only those jobs' checkouts: both NajaEDA
requirements become `najaeda==0.7.24`, and the CMake option
`KEPLER_USE_PUBLISHED_NAJAEDA` is enabled. These edits are not committed:
repository defaults and development jobs stay on the development provider
with the option off. Build and runtime pins must match, so installations of
a published-provider Kepler wheel select the same NajaEDA release.

Published NajaEDA does not supply the development shared-runtime SDK. This
compatibility path is deliberately pinned to `0.7.24` and uses released native
symbols, including the private DNL cache integration. It is not a general ABI
guarantee for arbitrary NajaEDA versions; any provider upgrade requires a new
compatibility review and the complete wheel tests. It does not change NajaEDA
or bundle a second Naja runtime.

## Publish a release

`publish` is the only upload switch. Its additional jobs download and link
against NajaEDA's actual distributed wheels instead of rebuilding a
same-version provider locally. Publication waits only for the
published-provider matrix and uploads only its wheels. Development-provider
jobs continue as regression checks; their status does not block publication.

1. Merge the intended code and workflow changes into `main` and confirm its
   checks pass. Choose an unused PyPI release version. The Python package
   version is read from `src/bin/KeplerVersion.h.in`, the same source used by
   CMake and Python package metadata;
   merge any necessary version change before starting the workflow. The
   publishing workflow does not bump or override any version declarations.
2. Open **Actions → Python wheels → Run workflow**, select `main`, check
   `publish`, and enter that exact version in `version`.
3. Four development jobs and four published-provider jobs build, repair, and
   test wheels in parallel. Published-provider jobs validate the version and
   print the selected commit SHA in their summaries. A mismatched or missing
   version, another branch, or another repository fails validation.
   Publication waits for the four published-provider jobs to succeed,
   independently of the development-provider results.
4. Review the commit, version, test results, and artifacts. Approve the waiting
   `pypi` environment deployment to allow the publishing job to run.
5. Confirm the release on [PyPI](https://pypi.org/project/kepler-formal/) and
   test installation in a fresh virtual environment on a supported platform:

   ```bash
   python -m pip install "kepler-formal==<released-version>"
   python -c "import kepler_formal; from kepler_formal import najaeda; print(kepler_formal.__version__)"
   ```

## Wheel coverage

Each provider matrix matches the vendored NajaEDA workflow: **26 wheels**.
A publish run tests both matrices and uploads only the 26 published-provider
wheels.

| Platform | Architecture | CPython versions | Wheels |
| --- | --- | --- | --- |
| Linux (`manylinux_2_28`) | x86_64 | 3.10–3.15, plus 3.14t | 7 |
| Linux (`manylinux_2_28`) | aarch64 | 3.10–3.15 | 6 |
| macOS | arm64 | 3.10–3.15, plus 3.14t | 7 |
| Windows | AMD64 | 3.10–3.15 | 6 |

Python 3.15 prerelease builds are explicitly enabled in cibuildwheel. The
3.14t wheels target free-threaded interpreters but require the GIL for native
verification; they do not promise concurrent, GIL-free engine calls.

`ci/check_wheel_matrix.py` compares cibuildwheel's actual selected build IDs
against NajaEDA's checked-in matrix and checks that publishing includes every
published-provider platform artifact. A Naja update that changes its matrix
requires an explicit corresponding update here. Every wheel runs the package tests and native
dependency/shared-runtime checks after repair. Windows reuses Naja's pregenerated
parser and vcpkg dependency approach, with KF-owned solver compatibility code.
The wheel smoke test loads designs through NajaEDA, verifies equivalent and
different pairs with all three solvers, and checks that the caller's designs
and selected top survive repeated verification and native errors. KF's Python
extension contains no Verilog file-loading frontend.

This workflow publishes wheels only, not a source distribution. Intel macOS,
musllinux, and Windows ARM64 wheels are not in NajaEDA's referenced matrix and
are not added here.

## Failure handling

Failed published-provider builds or missing published-provider platform
artifacts prevent publication. Development-provider failures remain visible
as regression failures and do not block the publishing job. PyPI upload keeps
metadata verification and attestations enabled, and does not silently
skip existing files. If an upload fails partway through, inspect the files
already present on PyPI before deciding how to recover; do not assume a rerun
can replace them. Use a new version for a rebuilt or changed release instead
of trying to overwrite published artifacts.

The workflow's `pypi` environment reference and PyPI publisher configuration
must match exactly. If authentication or approval fails, check the one-time
setup rather than adding credentials to the build jobs. See the
[PyPA publishing action documentation](https://github.com/pypa/gh-action-pypi-publish).
