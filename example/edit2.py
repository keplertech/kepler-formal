# SPDX-FileCopyrightText: 2024 The Naja authors
# <https://github.com/najaeda/naja/blob/main/AUTHORS>
#
# SPDX-License-Identifier: Apache-2.0

from os import path
import sys
import logging
from najaeda import netlist
from najaeda import naja

logging.basicConfig(level=logging.INFO)

# snippet-start: load_design_liberty
benchmarks = path.join('.')
liberty_files = [
    'NangateOpenCellLibrary_typical.lib',
    'fakeram45_1024x32.lib',
    'fakeram45_64x32.lib'
]
liberty_files = list(map(lambda p:path.join(benchmarks, p), liberty_files))
    
netlist.load_liberty(liberty_files)
top = netlist.load_verilog('tinyrocket.v')
u = naja.NLUniverse.get()
db = u.getTopDesign().getDB()
prims = list(db.getPrimitiveLibraries())
logic_1 = prims[0].getSNLDesign('LOGIC0_X1')
print(logic_1)
logic_1_inst = naja.SNLInstance.create(u.getTopDesign(), logic_1, "logic_1_inst")

inst = top.get_child_instance("logic_1_inst")
print(inst)
for term in inst.get_output_bit_terms():
    print(term)

net = None
index = 0
for input in top.get_input_bit_terms():
    if index == 2:
        net = input.get_lower_net()
        input.disconnect_lower_net()
        print(input)
        break
    index += 1

out = None
index = 0
for output in inst.get_output_bit_terms():
    out = output
    break

out.connect_upper_net(net) 

netlist.dump_naja_if('tinyrocket_naja_edited2.if')