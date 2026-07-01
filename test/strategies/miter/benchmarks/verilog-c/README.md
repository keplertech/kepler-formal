# verilog-c CLI Benchmark Fixtures

This directory contains C/RTL benchmark pairs copied from `rajdeep87/verilog-c` at commit `0b46ae072a477b3d8c00559e0df2483be470dcc4`.

Upstream describes these as ANSI-C benchmarks generated from Verilog RTL circuits with safety assertions. In Kepler's C2RTL flow, the runnable cases are used as C-vs-RTL equivalence fixtures: `kepler-formal` reads the C design through the hidden C frontend, translates it to SystemVerilog, and compares it against the paired upstream Verilog/SystemVerilog design.

`manifest.tsv` records all upstream `test.desc` benchmark entries. Entries marked `complete` have exactly one C source and one RTL source and are registered as CLI integration tests. The manifest also records the C and RTL top names used by the generated YAML configs. Entries marked `incomplete` do not have a complete C/RTL pair in the upstream snapshot, so they are retained as manifest metadata but not registered as equivalence cases.

The benchmark tests run the real `kepler-formal` executable. They require a real C frontend translator through `KEPLER_C_FRONTEND_TRANSLATOR`, `KEPLER_C2RTL_TRANSLATOR`, or `metron` on `PATH`. No fake translator is used by this benchmark suite.
