# In-tree Bazel registry

This registry holds `naja` at the commit the `thirdparty/naja` submodule
points to, and `naja-if` and `naja-verilog` at the commits naja's own
submodules point to there, until all three are released and on the
[Bazel Central Registry](https://registry.bazel.build/). Then it goes away.
naja's `MODULE.bazel` names older registered versions of naja-if and
naja-verilog than its sources need (its own builds override both with the
submodule checkouts), so naja's entry here declares the versions in this
registry instead.

Every other module that is not on BCR yet comes from its open
bazel-central-registry pull request; `.bazelrc` lists those by commit,
after this registry and ahead of BCR.

The layout is BCR's (`modules/naja/metadata.json`,
`modules/naja/<version>/{MODULE.bazel,source.json,presubmit.yml,overlay/}`).
`overlay/MODULE.bazel` is a copy of the version's `MODULE.bazel`; BCR
rejects symlinks.

## Bumping naja

1. Add `modules/naja/<version>/`, versioned `0.7.26-<YYYYMMDD>-<commit>`,
   with naja's `MODULE.bazel` at that commit (its `module(version = ...)`
   set to the directory name), copied to `overlay/MODULE.bazel`.
2. Put the version in `modules/naja/metadata.json`.
3. `bazel/registry/update_source.py naja <version> https://github.com/<owner>/naja/archive/<commit>.tar.gz naja-<commit>`
4. Point the `bazel_dep` in `MODULE.bazel` at it, and add any new BCR
   pull request registries naja's `.bazelrc` lists.
5. If naja's `thirdparty/naja-if` or `thirdparty/naja-verilog` moved, add
   the same kind of entry for it (`0.0.0-<YYYYMMDD>-<commit>`, the
   repository's `MODULE.bazel` with the version set, sources from
   `github.com/keplertech/<module>`), and name that version in naja's
   entry (`MODULE.bazel` and `overlay/MODULE.bazel`, then rerun
   `update_source.py` for naja so the overlay hash matches).

## Depending on kepler-formal from another module

List this registry by a pinned URL, then the BCR pull request registries
from kepler-formal's `.bazelrc`, then BCR:

```
common --registry=https://raw.githubusercontent.com/keplertech/kepler-formal/<commit>/bazel/registry/
common --registry=https://raw.githubusercontent.com/oharboe/bazel-central-registry/<sha>/   # one per BCR PR
common --registry=https://bcr.bazel.build/
```
