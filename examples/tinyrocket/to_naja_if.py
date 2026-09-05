# SPDX-FileCopyrightText: 2024 The Naja authors
# <https://github.com/najaeda/naja/blob/main/AUTHORS>
#
# SPDX-License-Identifier: Apache-2.0

from os import path
import logging

from najaeda import netlist

logging.basicConfig(level=logging.INFO)

# snippet-start: load_design_liberty
benchmarks = path.join(".")
liberty_files = [
    "NangateOpenCellLibrary_typical.lib",
    "fakeram45_1024x32.lib",
    "fakeram45_64x32.lib",
]
liberty_files = list(map(lambda p: path.join(benchmarks, p), liberty_files))

netlist.load_liberty(liberty_files)
top = netlist.load_verilog("tinyrocket.v")
netlist.dump_naja_if("tinyrocket_naja.if")
