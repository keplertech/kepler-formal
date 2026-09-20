#!/usr/bin/env python3
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Load two TinyRocket designs in NajaEDA, then verify them with Kepler Python."""

import argparse
from pathlib import Path

import najaeda
from najaeda import naja
import kepler_formal
from kepler_formal import VerificationOptions, VerificationStatus, verify_designs


EXAMPLE = Path(__file__).resolve().parent
LIBERTY_FILES = [
    EXAMPLE / "NangateOpenCellLibrary_typical.lib",
    EXAMPLE / "fakeram45_1024x32.lib",
    EXAMPLE / "fakeram45_64x32.lib",
    EXAMPLE / "fakeram45_64x15.lib",
]


def load_design(universe, path):
    # Separate databases allow both netlists to use the same module names.
    database = naja.NLDB.create(universe)
    database.loadLibertyPrimitives([str(path) for path in LIBERTY_FILES])
    database.loadVerilog([str(path)])
    design = database.getTopDesign()
    if design is None:
        raise RuntimeError(f"NajaEDA did not find a top design in {path}")
    return design


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", nargs="?", type=Path,
                        default=EXAMPLE / "tinyrocket.v")
    parser.add_argument("candidate", nargs="?", type=Path,
                        default=EXAMPLE / "tinyrocket_edited.v")
    parser.add_argument("--log-file", type=Path, default=Path("tinyrocket_python.log"))
    args = parser.parse_args()
    for path in (args.reference, args.candidate):
        if not path.is_file():
            parser.error(f"Design does not exist: {path}")

    print(f"NajaEDA:       {najaeda.__file__}", flush=True)
    print(f"Kepler Formal: {kepler_formal.__file__}", flush=True)

    # NajaEDA owns both designs in one universe. No Kepler file loader is used.
    universe = naja.NLUniverse.create()
    print(f"Loading reference: {args.reference}", flush=True)
    reference = load_design(universe, args.reference)
    print(f"Loading candidate: {args.candidate}", flush=True)
    candidate = load_design(universe, args.candidate)

    print("Checking equivalence through the Kepler Python library...", flush=True)
    result = verify_designs(
        reference,
        candidate,
        options=VerificationOptions(log_file=args.log_file),
    )
    print(f"Result: {result.status.value}")
    if result.reason:
        print(f"Reason: {result.reason}")
    print(f"Log: {result.log_file}")
    print(f"Designs still available in NajaEDA: {reference.getName()}, {candidate.getName()}")

    # LEC's native exit code is zero even for a difference; use the verdict.
    if result.status is VerificationStatus.EQUIVALENT:
        return 0
    return 1 if result.status is VerificationStatus.DIFFERENT else 2


if __name__ == "__main__":
    raise SystemExit(main())
