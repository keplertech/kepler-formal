#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from najaeda import naja
from najaeda import netlist


DEFAULT_SOURCE_INPUT = "auto_intsink_in_sync_0"
DEFAULT_TARGET_OUTPUT = "auto_tl_slave_xing_in_d_valid"
DEFAULT_OUTPUT_PATH = "tinyrocket_k3_edited.v"


def load_design(example_dir: Path):
    liberty_files = [
        "NangateOpenCellLibrary_typical.lib",
        "fakeram45_1024x32.lib",
        "fakeram45_64x32.lib",
        "fakeram45_64x15.lib",
    ]
    netlist.load_liberty([str(example_dir / name) for name in liberty_files])
    return netlist.load_verilog(str(example_dir / "tinyrocket.v"))


def get_primitive(model_name: str):
    universe = naja.NLUniverse.get()
    primitive_libraries = list(universe.getTopDesign().getDB().getPrimitiveLibraries())
    if not primitive_libraries:
        raise RuntimeError("No primitive library is loaded")
    model = primitive_libraries[0].getSNLDesign(model_name)
    if model is None:
        raise RuntimeError(f"Could not find primitive {model_name}")
    return model


def create_instance(model_name: str, instance_name: str):
    return naja.SNLInstance.create(
        naja.NLUniverse.get().getTopDesign(), get_primitive(model_name), instance_name
    )


def find_top_term(terms, name: str):
    for term in terms:
        if term.get_name() == name:
            return term
    raise RuntimeError(f"Could not find top-level term {name}")


def find_instance_term(instance, pin_name: str):
    model_term = instance.getModel().getScalarTerm(pin_name)
    if model_term is None:
        raise RuntimeError(f"Could not find pin {pin_name} on {instance.getName()}")
    return instance.getInstTerm(model_term)


def connect(instance, pin_name: str, net):
    term = find_instance_term(instance, pin_name)
    term.setNet(getattr(net, "net", net))


def main():
    parser = argparse.ArgumentParser(
        description="Inject a reset-initialized 3-stage pipeline into an existing TinyRocket output."
    )
    parser.add_argument(
        "--source-input",
        default=DEFAULT_SOURCE_INPUT,
        help="Top-level input that feeds the injected 3-stage pipeline.",
    )
    parser.add_argument(
        "--target-output",
        default=DEFAULT_TARGET_OUTPUT,
        help="Existing top-level output to perturb through XOR.",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_PATH,
        help="Edited Verilog output path.",
    )
    args = parser.parse_args()

    example_dir = Path(__file__).resolve().parent
    top = load_design(example_dir)

    source_input = find_top_term(top.get_input_bit_terms(), args.source_input)
    clock_input = find_top_term(top.get_input_bit_terms(), "clock")
    reset_input = find_top_term(top.get_input_bit_terms(), "reset")
    target_output = find_top_term(top.get_output_bit_terms(), args.target_output)

    source_net = source_input.get_lower_net()
    clock_net = clock_input.get_lower_net()
    reset_net = reset_input.get_lower_net()
    original_output_net = target_output.get_lower_net()
    if source_net is None or clock_net is None or reset_net is None:
        raise RuntimeError("Expected source, clock, and reset to be connected")
    if original_output_net is None:
        raise RuntimeError(f"{args.target_output} is not connected")

    # Keep the original TinyRocket output cone intact on a private internal net
    # and drive the existing public output port from the XORed result.
    target_output.disconnect_lower_net()
    original_output_net.set_name(f"{args.target_output}_orig")

    reset_n_net = top.create_net("k3_reset_n")
    stage1_d_net = top.create_net("k3_stage1_d")
    stage1_q_net = top.create_net("k3_stage1_q")
    stage1_qn_net = top.create_net("k3_stage1_qn")
    stage2_d_net = top.create_net("k3_stage2_d")
    stage2_q_net = top.create_net("k3_stage2_q")
    stage2_qn_net = top.create_net("k3_stage2_qn")
    stage3_d_net = top.create_net("k3_stage3_d")
    stage3_q_net = top.create_net("k3_stage3_q")
    stage3_qn_net = top.create_net("k3_stage3_qn")
    edited_output_net = top.create_net("k3_edited_output")

    # Keep the added state reset-initialized by forcing every stage D input to
    # zero when the TinyRocket top-level reset is asserted.
    reset_inv = create_instance("INV_X1", "k3_reset_inv")
    gate1 = create_instance("AND2_X1", "k3_stage1_gate")
    gate2 = create_instance("AND2_X1", "k3_stage2_gate")
    gate3 = create_instance("AND2_X1", "k3_stage3_gate")
    ff1 = create_instance("DFF_X1", "k3_stage1_ff")
    ff2 = create_instance("DFF_X1", "k3_stage2_ff")
    ff3 = create_instance("DFF_X1", "k3_stage3_ff")
    xor_out = create_instance("XOR2_X1", "k3_output_xor")

    connect(reset_inv, "A", reset_net)
    connect(reset_inv, "ZN", reset_n_net)

    connect(gate1, "A1", source_net)
    connect(gate1, "A2", reset_n_net)
    connect(gate1, "ZN", stage1_d_net)

    connect(ff1, "D", stage1_d_net)
    connect(ff1, "CK", clock_net)
    connect(ff1, "Q", stage1_q_net)
    connect(ff1, "QN", stage1_qn_net)

    connect(gate2, "A1", stage1_q_net)
    connect(gate2, "A2", reset_n_net)
    connect(gate2, "ZN", stage2_d_net)

    connect(ff2, "D", stage2_d_net)
    connect(ff2, "CK", clock_net)
    connect(ff2, "Q", stage2_q_net)
    connect(ff2, "QN", stage2_qn_net)

    connect(gate3, "A1", stage2_q_net)
    connect(gate3, "A2", reset_n_net)
    connect(gate3, "ZN", stage3_d_net)

    connect(ff3, "D", stage3_d_net)
    connect(ff3, "CK", clock_net)
    connect(ff3, "Q", stage3_q_net)
    connect(ff3, "QN", stage3_qn_net)

    connect(xor_out, "A", original_output_net)
    connect(xor_out, "B", stage3_q_net)
    connect(xor_out, "Z", edited_output_net)

    target_output.connect_lower_net(edited_output_net)

    top.dump_verilog(str((example_dir / args.output).resolve()))


if __name__ == "__main__":
    main()
