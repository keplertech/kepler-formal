#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from najaeda import naja
from najaeda import netlist


def load_design(example_dir: Path):
    liberty_files = [
        "NangateOpenCellLibrary_typical.lib",
        "fakeram45_1024x32.lib",
        "fakeram45_64x32.lib",
        "fakeram45_64x15.lib",
    ]
    netlist.load_liberty([str(example_dir / name) for name in liberty_files])
    return netlist.load_verilog(str(example_dir / "tinyrocket.v"))


def create_const_instance(top, logic_model_name: str, instance_name: str):
    universe = naja.NLUniverse.get()
    primitive_libraries = list(universe.getTopDesign().getDB().getPrimitiveLibraries())
    if not primitive_libraries:
        raise RuntimeError("No primitive library is loaded")
    logic_model = primitive_libraries[0].getSNLDesign(logic_model_name)
    if logic_model is None:
        raise RuntimeError(f"Could not find primitive {logic_model_name}")
    return naja.SNLInstance.create(universe.getTopDesign(), logic_model, instance_name)


def find_instance_output_term(instance, pin_name: str):
    for term in instance.get_output_bit_terms():
        if term.get_name() == pin_name:
            return term
    available = ", ".join(term.get_name() for term in instance.get_output_bit_terms())
    raise RuntimeError(
        f"Instance {instance} has no output pin {pin_name}. Available outputs: {available}"
    )


def find_instance_input_term(instance, pin_name: str):
    for term in instance.get_input_bit_terms():
        if term.get_name() == pin_name:
            return term
    available = ", ".join(term.get_name() for term in instance.get_input_bit_terms())
    raise RuntimeError(
        f"Instance {instance} has no input pin {pin_name}. Available inputs: {available}"
    )


def rewire_input_to_constant(top, target_instance, pin_name: str, const_value: str):
    target_input = find_instance_input_term(target_instance, pin_name)
    original_net = target_input.get_upper_net()
    if original_net is None:
        raise RuntimeError(f"{target_instance.get_name()}.{pin_name} is not connected")

    safe_prefix = "".join(
        ch if ch.isalnum() or ch == "_" else "_" for ch in target_instance.get_name()
    )
    logic_model_name = "LOGIC1_X1" if const_value == "1" else "LOGIC0_X1"
    const_instance_name = f"{safe_prefix}_{pin_name}_const_{const_value}"

    target_input.disconnect_upper_net()
    const_net = top.create_net(f"{safe_prefix}_{pin_name}_forced")
    target_input.connect_upper_net(const_net)

    create_const_instance(top, logic_model_name, const_instance_name)
    const_instance = top.get_child_instance(const_instance_name)
    const_output = next(iter(const_instance.get_output_bit_terms()))
    const_output.connect_upper_net(const_net)


def rewire_output_to_constant(top, target_instance, pin_name: str, const_value: str):
    target_output = find_instance_output_term(target_instance, pin_name)
    original_net = target_output.get_upper_net()
    if original_net is None:
        raise RuntimeError(f"{target_instance.get_name()}.{pin_name} is not connected")

    # Preserve a valid single-driver topology: detach the register output from
    # the old fanout net, park it on an unused private net, and drive the
    # original net from a fresh constant cell.
    safe_prefix = "".join(
        ch if ch.isalnum() or ch == "_" else "_" for ch in target_instance.get_name()
    )

    target_output.disconnect_upper_net()
    parked_net = top.create_net(f"{safe_prefix}_{pin_name}_parked")
    target_output.connect_upper_net(parked_net)

    logic_model_name = "LOGIC1_X1" if const_value == "1" else "LOGIC0_X1"
    const_instance_name = f"{safe_prefix}_{pin_name}_const_{const_value}"
    create_const_instance(
        top,
        logic_model_name,
        const_instance_name,
    )
    const_instance = top.get_child_instance(const_instance_name)
    const_output = next(iter(const_instance.get_output_bit_terms()))
    const_output.connect_upper_net(original_net)


def main():
    parser = argparse.ArgumentParser(
        description="Replace a TinyRocket instance pin with a constant driver."
    )
    parser.add_argument("--instance", required=True, help="Hierarchical instance name")
    parser.add_argument(
        "--role",
        choices=("input", "output"),
        default="input",
        help="Whether to override an instance input pin or output pin",
    )
    parser.add_argument(
        "--pin",
        help="Pin to override. Defaults to D for input edits and Q for output edits",
    )
    parser.add_argument(
        "--const",
        choices=("0", "1"),
        default="0",
        help="Constant value that should drive the original output net",
    )
    parser.add_argument("--output", required=True, help="Edited Verilog output path")
    args = parser.parse_args()

    example_dir = Path(__file__).resolve().parent
    top = load_design(example_dir)

    target_instance = top.get_child_instance(args.instance)
    if target_instance is None:
        raise RuntimeError(f"Could not find instance {args.instance}")

    pin_name = args.pin or ("D" if args.role == "input" else "Q")
    if args.role == "input":
        rewire_input_to_constant(top, target_instance, pin_name, args.const)
    else:
        rewire_output_to_constant(top, target_instance, pin_name, args.const)

    top.dump_verilog(str(Path(args.output).resolve()))


if __name__ == "__main__":
    main()
