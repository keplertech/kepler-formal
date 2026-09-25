# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""End-to-end borrowed latch verification using explicitly modeled Naja cells."""

import tempfile
import unittest
from pathlib import Path

from najaeda import naja, netlist
from kepler_formal import VerificationOptions, VerificationStatus, from_najaeda, verify_designs


class BorrowedLatchApiTest(unittest.TestCase):
    def setUp(self):
        netlist.reset()
        self.temporary = tempfile.TemporaryDirectory(prefix="kepler_latch_api_")
        self.universe = naja.NLUniverse.create()
        self.database = naja.NLDB.create(self.universe)
        self.designs = naja.NLLibrary.create(self.database, "designs")
        self.primitives = naja.NLLibrary.createPrimitives(self.database, "primitives")
        self.model = self.make_model("declared_latch")
        self.first = self.make_top("first", self.model)
        self.second = self.make_top("second", self.model)
        self.universe.setTopDesign(self.first)

    def tearDown(self):
        netlist.reset()
        self.temporary.cleanup()

    def make_model(self, name, invert=False):
        model = naja.SNLDesign.createPrimitive(self.primitives, name)
        for pin in ("D", "E"):
            naja.SNLScalarTerm.create(model, naja.SNLTerm.Direction.Input, pin)
        output = naja.SNLScalarTerm.create(model, naja.SNLTerm.Direction.Output, "Q")
        model.setSequentialModel(clocked_on="E", kind="latch",
                                 states=[{"name": "H", "next_state": "D"}],
                                 outputs=[(output, "!H" if invert else "H")])
        return model

    def make_top(self, name, model):
        top = naja.SNLDesign.create(self.designs, name)
        cell = naja.SNLInstance.create(top, model, "latch")
        for pin in model.getScalarTerms():
            net = naja.SNLScalarNet.create(top, pin.getName())
            naja.SNLScalarTerm.create(top, pin.getDirection(), pin.getName()).setNet(net)
            cell.getInstTerm(pin).setNet(net)
        return top

    def options(self, enabled=True, **changes):
        values = dict(mode="sec", sec_engine="k_induction", sec_encoding="binary",
                      log_file=Path(self.temporary.name) / "verification.log")
        if enabled:
            values.update(latch_support=True, latch_input_changes="single",
                          latch_initial_inputs=0, latch_initial_storage=0)
        return VerificationOptions(**(values | changes))

    def unchanged(self):
        self.assertIs(naja.NLUniverse.get(), self.universe)
        self.assertIs(self.universe.getTopDesign(), self.first)
        self.assertEqual("first", self.first.getName())
        self.assertEqual("second", self.second.getName())
        self.assertEqual(1, len(list(self.first.getInstances())))
        self.assertEqual(1, len(list(self.second.getInstances())))
        self.assertTrue(self.model.hasSequentialModel())

    def test_default_then_enabled_then_default_do_not_leak(self):
        default = verify_designs(self.first, self.second, options=self.options(False))
        self.assertEqual(0, default.covered_outputs)
        self.unchanged()
        enabled = verify_designs(self.first, self.second, options=self.options())
        self.assertEqual(VerificationStatus.EQUIVALENT, enabled.status)
        self.assertEqual(1, enabled.covered_outputs)
        self.assertEqual(1, enabled.total_outputs)
        self.unchanged()
        repeated = verify_designs(self.first, self.second, options=self.options(False))
        self.assertEqual(0, repeated.covered_outputs)
        self.unchanged()

    def test_handles_and_raw_designs_share_the_event_contract(self):
        for left, right in ((from_najaeda(self.first), from_najaeda(self.second)),
                            (self.first, self.second), (self.second, self.first)):
            result = verify_designs(left, right, options=self.options())
            self.assertEqual(VerificationStatus.EQUIVALENT, result.status)
            self.assertEqual(1, result.covered_outputs)
            self.unchanged()

    def test_inverted_latch_is_different_on_either_side(self):
        different = self.make_top("different", self.make_model("inverted", invert=True))
        for left, right in ((self.first, different), (different, self.first)):
            result = verify_designs(left, right, options=self.options())
            self.assertEqual(VerificationStatus.DIFFERENT, result.status)
            self.unchanged()

    def test_any_change_race_remains_opaque(self):
        result = verify_designs(self.first, self.second, options=self.options(latch_input_changes="any"))
        self.assertEqual(0, result.covered_outputs)
        self.unchanged()
        strict = verify_designs(self.first, self.second,
                                options=self.options(latch_input_changes="any", error_on_opaque=True))
        self.assertEqual(VerificationStatus.UNSUPPORTED, strict.status)
        self.unchanged()

    def test_explicit_one_initialization_is_supported(self):
        result = verify_designs(self.first, self.second,
                                options=self.options(latch_initial_inputs=1, latch_initial_storage=1))
        self.assertEqual(VerificationStatus.EQUIVALENT, result.status)
        self.assertEqual(1, result.covered_outputs)
        self.unchanged()

    def test_invalid_contract_preserves_live_designs(self):
        for changes in (dict(latch_initial_storage=None), dict(latch_support=False),
                        dict(mode="lec", sec_engine=None, sec_encoding=None),
                        dict(set_as_boundary=(("latch", "latch"),))):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                verify_designs(self.first, self.second, options=self.options(**changes))
            self.unchanged()

    def test_strict_opaque_policy_checks_unused_cells_on_either_side(self):
        unknown = naja.SNLDesign.createPrimitive(self.primitives, "unknown")
        output = naja.SNLScalarTerm.create(unknown, naja.SNLTerm.Direction.Output, "Y")
        cell = naja.SNLInstance.create(self.second, unknown, "unused_opaque")
        cell.getInstTerm(output).setNet(naja.SNLScalarNet.create(self.second, "unused_net"))
        for left, right in ((self.first, self.second), (self.second, self.first)):
            result = verify_designs(left, right, options=self.options(error_on_opaque=True))
            self.assertEqual(VerificationStatus.UNSUPPORTED, result.status)
            self.assertIn("unused_opaque", result.reason)
        relaxed = verify_designs(self.first, self.second, options=self.options())
        self.assertEqual(VerificationStatus.EQUIVALENT, relaxed.status)
        self.assertEqual(1, relaxed.covered_outputs)
        self.assertEqual(2, len(list(self.second.getInstances())))


if __name__ == "__main__":
    unittest.main()
