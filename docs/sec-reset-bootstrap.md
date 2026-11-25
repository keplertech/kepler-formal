# SEC Reset Bootstrap

SEC reset bootstrap constrains user-named top-level reset inputs before the
normal SEC property is checked. Use it for designs whose state is initialized by
a reset sequence rather than by explicit initial values.

## YAML

```yaml
verification: sec
sec_reset:
  cycles: 100
  ports:
    - name: reset
      active_value: 1
```

For multiple reset ports, list each port:

```yaml
sec_reset:
  cycles: 100
  ports:
    - name: reset
      active_value: 1
    - name: scan_reset_n
      active_value: 0
```

## CLI

```sh
kepler-formal -verilog -v sec \
  --sec-reset-cycles 100 \
  --sec-reset-port reset=1 \
  design0.v design1.v
```

Repeat `--sec-reset-port` for multiple ports:

```sh
kepler-formal -verilog -v sec \
  --sec-reset-cycles 100 \
  --sec-reset-port reset=1 \
  --sec-reset-port scan_reset_n=0 \
  design0.v design1.v
```

## Semantics

`cycles` must be greater than zero. During those cycles, each listed port is
held at its `active_value`. After the bootstrap window, the same ports are held
at the opposite value while SEC checks the post-reset behavior.

Port names are resolved after SEC aligns top-level inputs between the two
designs. A scalar top input can be named as `reset` even if the aligned bit name
is reported as `reset[0]`. For a bus reset, name each reset bit explicitly, such
as `reset_bus[0]`.

A reset sequence can leave some state uninitialized. SEC finds that state by
propagating the reset values through the next-state logic for the bootstrap
cycles with every other input unknown, then releasing the reset ports and
following whatever can still read a register left free. In the `binary`
encoding, the outputs whose cones read such a register are skipped as
reset-unanchored, exactly as without a bootstrap, and listed in
`skipped_reset_unanchored_pos.txt`: the two copies of a free register may hold
different values, so a mismatch on them would say nothing about the designs.
For the remaining outputs SEC assumes that they agree on the first frame after
the bootstrap window. Before the engine runs, SEC checks that some reset trace
can satisfy that assumption. If none can, the designs already differ on that
frame, and SEC reports the difference with a counterexample instead of proving
equivalence from an assumption that cannot hold. The dual-rail encoding keeps
every output and carries the free state as X instead.

## Incomplete reset warning

In the dual-rail encoding SEC reports a reset that does not initialise every
state bit. It takes each state bit at its boot value, X where the design gives
none, holds the reset ports active for the bootstrap cycles with every other
input a 0 or a 1, and counts the bits X can still reach. Without a reset
bootstrap the count is the state without a boot value. When the count is not
zero the log carries a warning next to the verdict, with the first eight names
of each design:

```text
SEC reset is incomplete: after the reset window of 3 cycle(s), 128 of 134 state bits of design 1 and 132 of 139 state bits of design 2 can still be unknown (X). Design 1: ..., and 120 more. Design 2: ..., and 124 more.
```

A second line follows when the binary passes may run:

```text
SEC binary passes are on: they start matched unknown registers of the two designs from the same value. Set binary_relation_round to false (--binary-relation-round false) to switch them off.
```

The binary passes are the binary round and signal correspondence on the binary
functions, described in
[sec-dual-rail-binary-passes.md](sec-dual-rail-binary-passes.md). One switch covers
both: `binary_relation_round` (`--binary-relation-round <bool>`), on by
default. With a complete reset no state is unknown when the comparison starts,
so the common start value of unknown registers plays no part in what is
compared, and there is no warning. A longer window, or reset logic that
reaches the listed state, removes the warning.

With `KEPLER_SEC_DIAG=1` or `KEPLER_SEC_SUMMARY_STATS=1` the counts are also
printed before the proof starts, as
`SEC summary: reset_unknown_state reset_cycles=... design0=.../... design1=.../...`.

Reset bootstrap is SEC-only. It does not discover internal reset nets, infer
reset polarity, or relax top-level input matching; every listed reset port must
exist as an aligned top-level input in both designs.
