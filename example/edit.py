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
netlist.get_top().dump_verilog('tinyrocket_pre_edited.v')
netlist.dump_naja_if('tinyrocket.if')
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
        ## assert input has no net
        if input.get_lower_net() is not None:
            print("net", input.get_lower_net())
            raise TypeError(f"Not disconnected: {input}")
        print(input)
        break
    index += 1

out = None
index = 0
for output in inst.get_output_bit_terms():
    out = output
    break

out.connect_upper_net(net) 
# insure net has only one driver
drivers = out.get_equipotential().get_leaf_drivers()
number_of_drivers = 0
for term in net.get_inst_terms():
    print("inst term:", term)
for term in net.get_design_terms():
    print("design term:", term)
for driver in drivers:
    print("driver", driver)
    number_of_drivers += 1
top_drivers = out.get_equipotential().get_top_drivers()
for top_driver in top_drivers:
    print("top driver", top_driver)
    number_of_drivers += 1
if number_of_drivers > 1:
    raise TypeError(f"Net has multiple drivers: {net}, {drivers}")
net.set_name("edit")
netlist.dump_naja_if('tinyrocket_naja_edited.if')
netlist.get_top().dump_verilog('tinyrocket_edited.v')
top.dump_verilog('tinyrocket_edited2.v')