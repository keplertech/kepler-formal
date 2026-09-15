# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

from __future__ import annotations

import gc
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import najaeda
from najaeda import netlist

from kepler_formal import (
    NativeDesign,
    SecEncoding,
    SecEngine,
    VerificationMode,
    VerificationOptions,
    VerificationStatus,
    from_najaeda,
    verify_designs,
)
from kepler_formal import _native


def _create_design(library, primitives, name: str, *, invert: bool = False):
    model = najaeda.naja.SNLDesign.createPrimitive(primitives, f"{name}_lut")
    model_input = najaeda.naja.SNLScalarTerm.create(
        model, najaeda.naja.SNLTerm.Direction.Input, "A"
    )
    model_output = najaeda.naja.SNLScalarTerm.create(
        model, najaeda.naja.SNLTerm.Direction.Output, "Y"
    )
    init = najaeda.naja.SNLParameter.create_binary(model, "INIT", 2, 0b10)
    model.setTruthTableFromParameter(model_output, [model_input], init)

    design = najaeda.naja.SNLDesign.create(library, name)
    input_term = najaeda.naja.SNLScalarTerm.create(
        design, najaeda.naja.SNLTerm.Direction.Input, "a"
    )
    output_term = najaeda.naja.SNLScalarTerm.create(
        design, najaeda.naja.SNLTerm.Direction.Output, "y"
    )
    input_net = najaeda.naja.SNLScalarNet.create(design, "a")
    input_term.setNet(input_net)
    output_net = najaeda.naja.SNLScalarNet.create(design, "y")
    output_term.setNet(output_net)
    instance = najaeda.naja.SNLInstance.create(design, model, "lut")
    instance.getInstTerm(model_input).setNet(input_net)
    instance.getInstTerm(model_output).setNet(output_net)
    parameter = najaeda.naja.SNLInstParameter.create(
        instance, init, "2'h1" if invert else "2'h2"
    )
    return design, parameter


class LiveNajaedaInteroperabilityTest(unittest.TestCase):
    def setUp(self):
        netlist.reset()
        self.universe = najaeda.naja.NLUniverse.create()
        self.database = najaeda.naja.NLDB.create(self.universe)
        self.library = najaeda.naja.NLLibrary.create(self.database, "designs")
        self.primitives = najaeda.naja.NLLibrary.createPrimitives(
            self.database, "primitives"
        )
        self.reference, self.reference_parameter = _create_design(
            self.library, self.primitives, "reference"
        )
        self.candidate, self.candidate_parameter = _create_design(
            self.library, self.primitives, "candidate", invert=True
        )
        self.temporary = tempfile.TemporaryDirectory(
            prefix="kepler_formal_live_najaeda_test_"
        )
        self.root = Path(self.temporary.name)

    def tearDown(self):
        netlist.reset()
        self.temporary.cleanup()

    def _lec_options(self, name: str = "borrowed.log") -> VerificationOptions:
        return VerificationOptions(log_file=self.root / name)

    def test_instance_capture_is_stable_and_verifies_without_copying(self):
        self.universe.setTopDesign(self.reference)
        reference_instance = netlist.get_top()
        captured_reference = from_najaeda(reference_instance)

        self.universe.setTopDesign(self.candidate)
        candidate_instance = netlist.get_top()
        captured_candidate = from_najaeda(candidate_instance)

        self.assertIs(captured_reference.source, reference_instance)
        self.assertIs(captured_reference.najaeda_design, self.reference)
        self.assertIs(captured_candidate.source, candidate_instance)
        self.assertIs(captured_candidate.najaeda_design, self.candidate)
        self.assertIsInstance(captured_reference, NativeDesign)

        different = verify_designs(
            captured_reference,
            captured_candidate,
            options=self._lec_options("different.log"),
        )
        self.assertEqual(VerificationStatus.DIFFERENT, different.status)

        # Mutating the original NajaEDA object changes the next result.  The
        # NativeDesign handle therefore points at that object; it is not a
        # serialized or reconstructed snapshot.
        self.candidate_parameter.setValue("2'h2")
        equivalent = verify_designs(
            captured_reference,
            captured_candidate,
            options=self._lec_options("equivalent.log"),
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, equivalent.status)
        self.assertIs(captured_candidate.najaeda_design, self.candidate)

    def test_raw_designs_are_accepted_and_handles_retain_their_wrappers(self):
        candidate = self.candidate
        candidate_wrapper_id = id(candidate)
        handle = from_najaeda(candidate)
        self.candidate = None
        del candidate
        gc.collect()

        self.assertEqual(candidate_wrapper_id, id(handle.source))
        self.candidate_parameter.setValue("2'h2")
        result = verify_designs(
            self.reference,
            handle,
            options=self._lec_options("raw.log"),
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, result.status)

    def test_handles_survive_validation_and_native_exceptions(self):
        reference = from_najaeda(self.reference)
        candidate = from_najaeda(self.candidate)

        with self.assertRaisesRegex(ValueError, "SEC"):
            verify_designs(
                reference,
                candidate,
                options=VerificationOptions(max_k=1),
            )
        with self.assertRaisesRegex(TypeError, "unknown"):
            _native.verify_designs(reference, candidate, {"unknown": True})

        # A native logging failure occurs after the caller's state has been
        # borrowed. It must return an error and leave both handles reusable.
        top_before = self.universe.getTopDesign()
        failed = verify_designs(
            reference,
            candidate,
            options=VerificationOptions(mode=VerificationMode.SEC, log_file=self.root),
        )
        self.assertEqual(VerificationStatus.ERROR, failed.status)
        self.assertTrue(failed.reason)
        self.assertIs(najaeda.naja.NLUniverse.get(), self.universe)
        self.assertIs(self.universe.getTopDesign(), top_before)

        self.candidate_parameter.setValue("2'h2")
        result = verify_designs(
            reference,
            candidate,
            options=self._lec_options("after-exception.log"),
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, result.status)

    def test_invalid_objects_and_destroyed_designs_are_rejected(self):
        with self.assertRaisesRegex(TypeError, "SNLDesign.*Instance"):
            from_najaeda(object())
        with self.assertRaisesRegex(TypeError, "created by from_najaeda"):
            NativeDesign()

        self.universe.setTopDesign(self.reference)
        instance = netlist.get_top()
        with self.assertRaisesRegex(TypeError, "capture it with from_najaeda"):
            verify_designs(instance, self.candidate)

        handle = from_najaeda(self.candidate)
        self.candidate.destroy()
        with self.assertRaises(ReferenceError):
            verify_designs(handle, self.reference)

    def test_options_are_validated_before_native_object_conversion(self):
        with patch("kepler_formal.api._native.verify_designs") as native_verify:
            with self.assertRaisesRegex(TypeError, "VerificationOptions"):
                verify_designs(object(), object(), options={})
            native_verify.assert_not_called()

    def test_small_same_design_sec_run(self):
        result = verify_designs(
            self.reference,
            self.reference,
            options=VerificationOptions(
                mode=VerificationMode.SEC,
                sec_engine=SecEngine.K_INDUCTION,
                sec_encoding=SecEncoding.BINARY,
                max_k=1,
                log_file=self.root / "same-sec.log",
            ),
        )
        self.assertEqual(VerificationStatus.EQUIVALENT, result.status)
        self.assertEqual("naja_design", result.input_format)
        self.assertEqual(1, result.total_outputs)


if __name__ == "__main__":
    unittest.main()
