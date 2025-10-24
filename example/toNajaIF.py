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
# leafCount = 0
# for leaf in top.get_leaf_children():
#     leafCount += 1
# print("Leaf children count: ", leafCount)
# u = naja.NLUniverse.get()
# db = u.getTopDesign().getDB()
# prims = list(db.getPrimitiveLibraries())
# logic_1 = prims[0].getSNLDesign('LOGIC0_X1')
# print(logic_1)
# logic_1_inst = naja.SNLInstance.create(u.getTopDesign(), logic_1, "logic_1_inst")

# inst = top.get_child_instance("logic_1_inst")
# print(inst)
# for term in inst.get_output_bit_terms():
#     print(term)

# net = None
# for input in top.get_input_bit_terms():
#     net = input.get_lower_net()
#     input.disconnect_lower_net()
#     print(input)
#     break

# out = None
# for output in inst.get_output_bit_terms():
#     out = output
#     break

# out.connect_upper_net(net) 

netlist.dump_naja_if('black_parrot.v.if')