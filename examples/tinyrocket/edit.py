# SPDX-FileCopyrightText: 2024 The Naja authors
# <https://github.com/najaeda/naja/blob/main/AUTHORS>
#
# SPDX-License-Identifier: Apache-2.0

from os import path
import logging

from najaeda import naja
from najaeda import netlist

logging.basicConfig(level=logging.INFO)

LIBERTY_FILES = [
    "NangateOpenCellLibrary_typical.lib",
    "fakeram45_1024x32.lib",
    "fakeram45_64x32.lib",
]
TARGET_INPUT_NAME = "auto_intsink_in_sync_0"
SEC_PROBE_OUTPUT_NAME = "sec_edit_probe_o"


def liberty_paths():
    return [path.join(".", liberty_file) for liberty_file in LIBERTY_FILES]


def load_tinyrocket():
    netlist.load_liberty(liberty_paths())
    return netlist.load_verilog("tinyrocket.v")


def find_input_bit_term(top, name):
    for term in top.get_input_bit_terms():
        if term.get_name() == name:
            return term
    raise RuntimeError(f"Could not find top input {name}")


def first_output_bit_term(instance):
    for term in instance.get_output_bit_terms():
        return term
    raise RuntimeError(f"Instance {instance} has no output bit terms")


def create_logic0_driver(top):
    universe = naja.NLUniverse.get()
    db = universe.getTopDesign().getDB()
    primitive_libraries = list(db.getPrimitiveLibraries())
    logic0 = primitive_libraries[0].getSNLDesign("LOGIC0_X1")
    logic0_instance = naja.SNLInstance.create(
        universe.getTopDesign(), logic0, "logic_1_inst")
    return top.get_child_instance(logic0_instance.getName())


def add_sec_probe_output(source_net):
    # Expose the edited interrupt input as a real top output instead of relying
    # on LEC or any cross-design relation between internal elements.
    probe = naja.SNLScalarTerm.create(
        naja.NLUniverse.get().getTopDesign(),
        naja.SNLTerm.Direction.Output,
        SEC_PROBE_OUTPUT_NAME,
    )
    probe.setNet(source_net.net)


def tie_target_input_to_logic0(top, target_input):
    target_net = target_input.get_lower_net()
    target_input.disconnect_lower_net()
    if target_input.get_lower_net() is not None:
        raise RuntimeError(f"Not disconnected: {target_input}")

    logic0_instance = create_logic0_driver(top)
    first_output_bit_term(logic0_instance).connect_upper_net(target_net)

    drivers = first_output_bit_term(logic0_instance).get_equipotential().get_leaf_drivers()
    top_drivers = first_output_bit_term(logic0_instance).get_equipotential().get_top_drivers()
    if len(list(drivers)) + len(list(top_drivers)) > 1:
        raise RuntimeError(f"Net has multiple drivers: {target_net}")

    target_net.set_name("edit")


def dump_legacy_edit_outputs():
    top = load_tinyrocket()

    target_input = find_input_bit_term(top, TARGET_INPUT_NAME)
    tie_target_input_to_logic0(top, target_input)
    netlist.dump_naja_if("tinyrocket_naja_edited.if")
    netlist.get_top().dump_verilog("tinyrocket_edited.v")


def dump_sec_probe_outputs():
    naja.NLUniverse.get().destroy()
    top = load_tinyrocket()
    target_input = find_input_bit_term(top, TARGET_INPUT_NAME)
    add_sec_probe_output(target_input.get_lower_net())
    netlist.dump_naja_if("tinyrocket_sec_probe_naja.if")
    netlist.get_top().dump_verilog("tinyrocket_sec_probe.v")

    tie_target_input_to_logic0(top, target_input)
    netlist.dump_naja_if("tinyrocket_sec_probe_naja_edited.if")
    netlist.get_top().dump_verilog("tinyrocket_sec_probe_edited.v")


def main():
    dump_legacy_edit_outputs()
    dump_sec_probe_outputs()


if __name__ == "__main__":
    main()
