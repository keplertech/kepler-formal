# verilog-c CLI Benchmark Fixtures

This directory contains C/RTL benchmark pairs copied from `rajdeep87/verilog-c` at commit `0b46ae072a477b3d8c00559e0df2483be470dcc4`.

Upstream describes these as ANSI-C benchmarks generated from Verilog RTL circuits with safety assertions. `manifest.tsv` records the complete C/RTL pairs so they can become C-vs-RTL equivalence fixtures once the C2RTL frontend supports this ANSI-C verifier style.

Entries marked `complete` have exactly one C source and one RTL source and are registered as CLI integration tests. The runner skips cases whose C source is outside the currently supported synthesizable C2RTL subset instead of using a fake translator. Entries marked `incomplete` do not have a complete C/RTL pair in the upstream snapshot, so they are retained as manifest metadata but not registered as equivalence cases.

The benchmark tests run the real `kepler-formal` executable. The C-to-RTL translator is linked into that binary through Kepler's `src/c2rtl` wrapper around the vendored Metron sources. No fake translator is used by this benchmark suite.
