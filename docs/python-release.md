# Publishing the Python package

The **Python wheels** workflow (`.github/workflows/python-wheels.yml`) can
publish `kepler-formal` to PyPI through a manual run on `main`. Publishing is
off by default. PRs, `v*` tags, and ordinary manual builds only produce tested
wheel artifacts; they never publish to PyPI.

This is separate from [binary releases](releasing.md). Do not create a tag or
run `tools/release.sh` for a Python-only release: `v*` tags also trigger the
standalone binary release workflow. The manual Python workflow creates no
tags or GitHub Releases and does not publish the MCP add-on.

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
run's four wheel artifacts and uploads them without checking out or building
source in the privileged job.

See the official [PyPI setup instructions](https://docs.pypi.org/trusted-publishers/adding-a-publisher/),
[Trusted Publishing instructions](https://docs.pypi.org/trusted-publishers/using-a-publisher/),
and [GitHub environment protection documentation](https://docs.github.com/en/actions/reference/workflows-and-actions/deployments-and-environments).

## Build without publishing

After the workflow is present on the repository's default branch, open
**Actions → Python wheels → Run workflow**. Leave `publish` unchecked and
`version` empty. The workflow builds, repairs, and tests wheels, then stores
them as Actions artifacts. This can be used to check a branch before merging.

## Publish a release

The shared-runtime integration currently pins the **unreleased** NajaEDA SDK
`0.7.24.dev0`. Publication is deliberately blocked. First release that SDK
with the matching platform/Python wheels, then update the provider pin in
`pyproject.toml` (build and runtime requirements) and
`ci/shared_naja_wheels.py`. Development CI builds a local provider wheel;
after switching to a release pin, CI downloads the published provider instead.
This matters because a same-version rebuild is not necessarily ABI/build
identical to the provider that users install. Kepler and Naja validate native
build identity before sharing objects.

1. Merge the intended code and workflow changes into `main` and confirm its
   checks pass. Choose an unused PyPI release version. The Python package
   version is read from `project(kepler-formal VERSION ...)` in `CMakeLists.txt`;
   merge any necessary version change before starting the workflow. The
   publishing workflow does not bump or override any version declarations.
2. Open **Actions → Python wheels → Run workflow**, select `main`, check
   `publish`, and enter that exact version in `version`.
3. All four matrix jobs validate the version and print the selected commit SHA in
   their summaries, then build, repair, and test the wheels from that run's
   checkout. A mismatched or missing version, another branch, or another
   repository fails validation. Publication waits for all four jobs to succeed.
4. Review the commit, version, test results, and artifacts. Approve the waiting
   `pypi` environment deployment to allow the publishing job to run.
5. Confirm the release on [PyPI](https://pypi.org/project/kepler-formal/) and
   test installation in a fresh virtual environment on a supported platform:

   ```bash
   python -m pip install "kepler-formal==<released-version>"
   python -c "import kepler_formal; from kepler_formal import najaeda; print(kepler_formal.__version__)"
   ```

## Wheel coverage

The matrix matches the vendored NajaEDA workflow: **26 wheels** per release.

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
platform artifact. A Naja update that changes its matrix requires an explicit
corresponding update here. Every wheel runs the package tests and native
dependency/shared-runtime checks after repair. Windows reuses Naja's pregenerated
parser and vcpkg dependency approach, with KF-owned solver compatibility code.

This workflow publishes wheels only, not a source distribution. Intel macOS,
musllinux, and Windows ARM64 wheels are not in NajaEDA's referenced matrix and
are not added here.

## Failure handling

Failed builds or missing platform artifacts prevent publication. PyPI upload
keeps metadata verification and attestations enabled, and does not silently
skip existing files. If an upload fails partway through, inspect the files
already present on PyPI before deciding how to recover; do not assume a rerun
can replace them. Use a new version for a rebuilt or changed release instead
of trying to overwrite published artifacts.

The workflow's `pypi` environment reference and PyPI publisher configuration
must match exactly. If authentication or approval fails, check the one-time
setup rather than adding credentials to the build jobs. See the
[PyPA publishing action documentation](https://github.com/pypa/gh-action-pypi-publish).
