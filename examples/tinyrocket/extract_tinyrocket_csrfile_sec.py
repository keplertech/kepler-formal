#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

from najaeda import naja
from najaeda import netlist


TARGET_DESIGN = "CSRFile"
OUTPUT_FILE = "tinyrocket_csrfile.v"


def load_tinyrocket(example_dir: Path) -> None:
    liberty_files = [
        "NangateOpenCellLibrary_typical.lib",
        "fakeram45_1024x32.lib",
        "fakeram45_64x32.lib",
        "fakeram45_64x15.lib",
    ]
    netlist.load_liberty([str(example_dir / name) for name in liberty_files])
    netlist.load_verilog(str(example_dir / "tinyrocket.v"))


def find_design(name: str):
    universe = naja.NLUniverse.get()
    top = universe.getTopDesign()
    if top is None:
        raise RuntimeError("No top design was loaded")

    for library in top.getDB().getLibraries():
        for design in library.getSNLDesigns():
            if design.getName() == name:
                return design
    raise RuntimeError(f"Could not find design {name}")


def main() -> None:
    example_dir = Path(__file__).resolve().parent
    load_tinyrocket(example_dir)

    universe = naja.NLUniverse.get()
    target = find_design(TARGET_DESIGN)
    universe.setTopDesign(target)

    output_path = example_dir / OUTPUT_FILE
    netlist.get_top().dump_verilog(str(output_path))
    print(output_path)


if __name__ == "__main__":
    main()
