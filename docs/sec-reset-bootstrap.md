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

Reset bootstrap is SEC-only. It does not discover internal reset nets, infer
reset polarity, or relax top-level input matching; every listed reset port must
exist as an aligned top-level input in both designs.
