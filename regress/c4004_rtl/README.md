<!--
Copyright 2026 keplertech.io
SPDX-License-Identifier: Apache-2.0
-->

# Intel 4004 clone RTL trials

SEC inputs of issue 317: four RTL trials of one 4-bit CPU core behind the same
wrapper `cpu_top`, each with its 32x4 SRAM model.

`rtl_01.v`, `rtl_02.v`, `rtl_03.v` and `rtl_04.v` are copies of `keplertech/cfg/rtl_0N.v`
from <https://github.com/xtofalex/TR-1um_C4004_Clone2Build> at commit
`e8e387e`, licensed under Apache-2.0.

| File | Against `rtl_01.v` |
| --- | --- |
| `rtl_02.v` | Equivalent. Keeps the program counter in registers beside the stack memory. |
| `rtl_03.v` | Equivalent. Reads the index register one cycle early and writes the stack one cycle late for register instructions, so the memories agree at instruction ends only: needs `signal_correspondence_step: 8`. |
| `rtl_04.v` | Equivalent, with `signal_correspondence_step: 8`. |
| `rtl_02_bug.v` | Not equivalent. `rtl_02.v` with one line removed: the register copy stays valid after INC. |
| `rtl_03_bug.v` | Not equivalent. `rtl_03.v` with one line removed: the register copy stays valid after the early read. |

The reset `rst` is held for three cycles. The memory clock is
`~clk & rst_n & sram_ce_en`, so the memory must hold during reset.
