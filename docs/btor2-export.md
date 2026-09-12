# BTOR2 Export

KF can dump its prepared sequential equivalence obligation after reading both
designs and before running the selected proof engine. The file contains both
designs as one transition system, with a `bad` property for a covered output
mismatch.

```text
Read both designs → extract and align SEC models → write equivalence.btor2
                                              ↘ run the selected KF engine
```

## CLI

```sh
kepler-formal -verilog -v sec design0.v design1.v \
  --dump-btor2 equivalence.btor2
```

Add `--dump-only` to stop after writing the file:

```sh
kepler-formal -verilog -v sec design0.v design1.v \
  --dump-btor2 equivalence.btor2 --dump-only
```

Both export flags can appear before or after the input format selector, or
among the normal input options. `--dump-btor2` requires a non-empty file path.

## YAML

```yaml
format: verilog
verification: sec
input_paths: [design0.v, design1.v]
sec_encoding: dual_rail_steady
btor2_export: true
btor2_export_path: equivalence.btor2
dump_only: true
```

Run this with `kepler-formal --config config.yaml`. Config mode cannot be
combined with other CLI options.

| Key | Default | Meaning |
| --- | --- | --- |
| `btor2_export` | `false` | Enable BTOR2 export. |
| `btor2_export_path` | `miter.btor2` when enabled | File to write. Relative paths resolve against the run's working directory. Requires `btor2_export: true`. |
| `dump_only` | `false` | Stop after export, without running a proof engine. Requires export to be enabled. |

Export options require `verification: sec`. They work with normal and compact
SEC loading, including the identical-input model reuse path.

## Meaning of the dump

- The output is a bit-level BTOR2 model, primarily using one-bit bit-vectors.
  Arithmetic and memories have already been lowered by KF; export does not
  recover their original word-level structure.
- Shared environment inputs, design-local state, initial-state relations, and
  next-state logic are preserved. The mismatch property uses KF's selected
  output coverage and encoding.
- In `dual_rail_steady`, a mismatch requires both outputs to be binary-defined
  and opposite. Cycles with an unknown output remain outside that property.
- Startup and reset timing follow KF's concrete bounded-counterexample
  semantics. The effective `sec_reset` bootstrap is encoded; as in KF's
  counterexample checker, a complete binary initial-state assignment can
  bypass a configured reset prefix.
- The file describes the prepared transition-system obligation, rather than
  an engine's SAT clauses or learned proof state. The selected `max_k` is a KF
  search bound; it does not bound the exported transition system.
- Skipped outputs remain excluded. Export coverage is reported, so a dump
  with partial coverage must not be treated as a full-design equivalence
  obligation.

Successful `dump_only` exits with code `0` and reports **exported; proof not
run**, together with the file path and output coverage. This is an export
result, not an equivalence verdict. With `dump_only: false`, KF writes the file
and continues with normal proof results and exit codes. An export failure
returns an error rather than continuing silently.
