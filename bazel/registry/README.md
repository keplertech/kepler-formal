# In-tree Bazel registry

This registry holds one module: `naja`, at the development pin
kepler-formal is tested against (najaeda/naja#457 merged with #455,
`oharboe/naja@kepler-pin`), until naja is released and on the
[Bazel Central Registry](https://registry.bazel.build/). Then it goes away.

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

## Depending on kepler-formal from another module

List this registry by a pinned URL, then the BCR pull request registries
from kepler-formal's `.bazelrc`, then BCR:

```
common --registry=https://raw.githubusercontent.com/keplertech/kepler-formal/<commit>/bazel/registry/
common --registry=https://raw.githubusercontent.com/oharboe/bazel-central-registry/<sha>/   # one per BCR PR
common --registry=https://bcr.bazel.build/
```
