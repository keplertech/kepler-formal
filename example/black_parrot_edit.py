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
    'fakeram45_64x32.lib',
    'fakeram45_64x15.lib',
    'fakeram45_64x96.lib',
    'fakeram45_512x64.lib',
    'fakeram45_256x95.lib',
    'fakeram45_64x7.lib'
]
liberty_files = list(map(lambda p:path.join(benchmarks, p), liberty_files))
    
netlist.load_liberty(liberty_files)
top = netlist.load_verilog('black_parrot.v')
netlist.dump_naja_if('black_parrot_pre_edit.if')
child_instance = None
childI = 0
for child in top.get_child_instances():
    num_children = 0
    for _ in child.get_child_instances():
        num_children += 1
        break
    if num_children > 0 and childI > 4:
        print(childI)
        child_instance = child
        break
    childI += 1
u = naja.NLUniverse.get()
db = u.getTopDesign().getDB()
prims = list(db.getPrimitiveLibraries())
logic_1 = prims[0].getSNLDesign('LOGIC0_X1')
print(logic_1)
logic_1_inst = naja.SNLInstance.create(u.getTopDesign().getInstance(child_instance.get_name()).getModel(), logic_1, "logic_1_inst")

inst = child_instance.get_child_instance("logic_1_inst")
print(inst)
for term in inst.get_output_bit_terms():
    print(term)

net = None
index = 0
for input in child_instance.get_input_bit_terms():
    if input.is_bus_bit(): 
        continue
    if index == 1:
        net = input.get_lower_net()
        if net.is_bus() or net.is_bus_bit():
            print("net is bus!")
            sys.exit(1)
        
        print("input: ", input)
        print("net :", net)
        for insterm in input.get_equipotential().get_inst_terms():
            print("  insterm: ", insterm)
        for term in input.get_flat_fanout():
            print(term)
        input.disconnect_lower_net()
        break
    index += 1

out = None
index = 0
for output in inst.get_output_bit_terms():
    out = output
    break

out.connect_upper_net(net) 

netlist.dump_naja_if('black_parrot_edit.if')