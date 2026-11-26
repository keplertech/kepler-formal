<!--
Copyright 2026 keplertech.io
SPDX-License-Identifier: Apache-2.0
-->

# Multi-cycle swap

A four-cycle step that saves a word in the upper half of a memory and reads a
word of the lower half. `mode_swapped_order.v` trades the two accesses when the
mode bit of the step is set, so against `fixed_order.v` the outputs are equal
in every cycle while the memory and one register differ inside the step. The
pair needs `signal_correspondence_step: 4`.

`mode_swapped_order_bug.v` saves a wrong value to word 0 and is not equivalent.
