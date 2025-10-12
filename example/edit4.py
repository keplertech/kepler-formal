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

child_instance = None
childI = 0
for child in top.get_child_instances():
    num_children = 0
    for _ in child.get_child_instances():
        num_children += 1
        break
    if num_children > 0 and childI > 2333:
        print(childI)
        child_instance = child
        break
    childI += 1
u = naja.NLUniverse.get()
db = u.getTopDesign().getDB()
prims = list(db.getPrimitiveLibraries())
logic_0 = prims[0].getSNLDesign('LOGIC0_X1')
inv = prims[0].getSNLDesign('INV_X1')
print(logic_0)
logic_0_inst = naja.SNLInstance.create(u.getTopDesign().getInstance(child_instance.get_name()).getModel(), logic_0, "logic_0_inst")
inv_inst = naja.SNLInstance.create(u.getTopDesign().getInstance(child_instance.get_name()).getModel(), inv, "inv_inst")
inst = child_instance.get_child_instance("logic_0_inst")
invInst = child_instance.get_child_instance("inv_inst")
inv_inst2 = naja.SNLInstance.create(u.getTopDesign().getInstance(child_instance.get_name()).getModel(), inv, "inv_inst2")
invInst2 = child_instance.get_child_instance("inv_inst2")

# create a net
inv_net = naja.SNLScalarNet.create(u.getTopDesign().getInstance(child_instance.get_name()).getModel(), "inv_net")
inv_net2 = naja.SNLScalarNet.create(u.getTopDesign().getInstance(child_instance.get_name()).getModel(), "inv_net2")
print(inst)
for term in inst.get_output_bit_terms():
    print(term)

net = None
index = 0
for input in child_instance.get_input_bit_terms():
    if index == 1:
        net = input.get_lower_net()
        input.disconnect_lower_net()
        print("input: ", input)
        break
    index += 1

out = None
index = 0
for output in inst.get_output_bit_terms():
    out = output
    break
outInv = None
for output in invInst.get_output_bit_terms():
    outInv = output
    break
outInv2 = None
for output in invInst2.get_output_bit_terms():
    outInv2 = output
    break
inInv = None
for input in invInst.get_input_bit_terms():
    inInv = input
    break
inInv2 = None
for input in invInst2.get_input_bit_terms():
    inInv2 = input
    break
inv_net_najaeda = child_instance.get_net("inv_net")
inv_net2_najaeda = child_instance.get_net("inv_net2")

out.connect_upper_net(inv_net_najaeda) 
inInv.connect_upper_net(inv_net_najaeda)
outInv.connect_upper_net(inv_net2_najaeda)
inInv2.connect_upper_net(inv_net2_najaeda)
outInv2.connect_upper_net(net)

netlist.dump_naja_if('tinyrocket_naja_edited4.if')