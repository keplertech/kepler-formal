# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Emit a kepler-formal python-format script that rebuilds a flat gate-level
Verilog netlist object for object, in the order the Verilog reader created
them, so that the loaded design is identical and the LEC CNF is the same.

Usage, with the built naja module on PYTHONPATH:
  python3 verilog_to_naja_py.py <cells.lib> <netlist.v> <out.py>
"""

import sys
import naja

liberty, verilog, out_path = sys.argv[1:4]
universe = naja.NLUniverse.create()
db = naja.NLDB.create(universe)
db.loadLibertyPrimitives([liberty])
db.loadVerilog([verilog])
top = db.getTopDesign()
if top is None:
    raise SystemExit("no top design")

DIR = {naja.SNLTerm.Direction.Input: "Input", naja.SNLTerm.Direction.Output: "Output", naja.SNLTerm.Direction.InOut: "InOut"}

def q(s):
    return repr(str(s))

net_index = {}

def net_var(net):
    return "n%d" % net_index[str(net.getNLID())]

def net_ref(net):
    """Python expression that evaluates to the given bit net or bus net."""
    if isinstance(net, naja.SNLBusNetBit):
        return "%s.getBit(%d)" % (net_var(net.getBus()), net.getBit())
    return net_var(net)

def bitterm_ref(design_expr, bit_term):
    if isinstance(bit_term, naja.SNLBusTermBit):
        return "%s.getBusTerm(%s).getBusTermBit(%d)" % (design_expr, q(bit_term.getBus().getName()), bit_term.getBit())
    return "%s.getScalarTerm(%s)" % (design_expr, q(bit_term.getName()))

lines = ["# Copyright 2026 keplertech.io", "# SPDX-License-" + "Identifier: Apache-2.0", "", "import naja", "", "",
         "def primitive(lib, name):",
         "    for primitives in lib.getDB().getPrimitiveLibraries():",
         "        model = primitives.getSNLDesign(name)",
         "        if model is not None:",
         "            return model",
         "    raise RuntimeError('primitive %s was not loaded' % name)", "", "",
         "def constructLibrary(lib):",
         "    top = naja.SNLDesign.create(lib, %s)" % q(top.getName())]
for term in top.getTerms():
    d = "naja.SNLTerm.Direction." + DIR[term.getDirection()]
    if isinstance(term, naja.SNLBusTerm):
        lines.append("    naja.SNLBusTerm.create(top, %s, %d, %d, %s)" % (d, term.getMSB(), term.getLSB(), q(term.getName())))
    else:
        lines.append("    naja.SNLScalarTerm.create(top, %s, %s)" % (d, q(term.getName())))
for net in top.getNets():
    net_index[str(net.getNLID())] = len(net_index)
    name = "" if net.isUnnamed() else ", " + q(net.getName())
    if net.isUnnamed():
        bits = list(net.getBits())
        print("unnamed net %s type=%s bits=%d components=%d" % (net_var(net), bits[0].getTypeAsString() if bits else "?", len(bits), sum(len(list(b.getComponents())) for b in bits)), file=sys.stderr)
    if isinstance(net, naja.SNLBusNet):
        lines.append("    %s = naja.SNLBusNet.create(top, %d, %d%s)" % (net_var(net), net.getMSB(), net.getLSB(), name))
    else:
        lines.append("    %s = naja.SNLScalarNet.create(top%s)" % (net_var(net), name))
    for bit in net.getBits():
        if bit.getType() != naja.SNLNet.Type.Standard:
            lines.append("    %s.setType(naja.SNLNet.Type.%s)" % (net_ref(bit), bit.getTypeAsString()))
for term in top.getTerms():
    for bit in term.getBits():
        net = bit.getNet()
        if net is not None:
            lines.append("    %s.setNet(%s)" % (bitterm_ref("top", bit), net_ref(net)))
for inst in top.getInstances():
    model = inst.getModel()
    if not model.isPrimitive():
        raise SystemExit("hierarchical instance %s: not supported by this generator" % inst.getName())
    lines.append("    inst = naja.SNLInstance.create(top, primitive(lib, %s), %s)" % (q(model.getName()), q(inst.getName())))
    for it in inst.getInstTerms():
        net = it.getNet()
        if net is not None:
            lines.append("    inst.getInstTerm(%s).setNet(%s)" % (bitterm_ref("inst.getModel()", it.getBitTerm()), net_ref(net)))
with open(out_path, "w") as f:
    f.write("\n".join(lines) + "\n")
print("%s: %d terms, %d nets, %d instances -> %s (%d lines)" % (top.getName(), len(list(top.getTerms())), len(list(top.getNets())), len(list(top.getInstances())), out_path, len(lines)))
