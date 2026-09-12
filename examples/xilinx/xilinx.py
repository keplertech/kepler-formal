# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Formal models for the Xilinx primitives used by these examples."""

import naja


def _construct_gate(lib, name, inputs, output, truth_table):
    primitive = naja.SNLDesign.createPrimitive(lib, name)
    for input_name in inputs:
        naja.SNLScalarTerm.create(
            primitive, naja.SNLTerm.Direction.Input, input_name
        )
    naja.SNLScalarTerm.create(
        primitive, naja.SNLTerm.Direction.Output, output
    )
    primitive.setTruthTable(truth_table)


def _construct_lut(lib, input_count):
    primitive = naja.SNLDesign.createPrimitive(lib, f"LUT{input_count}")
    inputs = [
        naja.SNLScalarTerm.create(
            primitive, naja.SNLTerm.Direction.Input, f"I{index}"
        )
        for index in range(input_count)
    ]
    output = naja.SNLScalarTerm.create(
        primitive, naja.SNLTerm.Direction.Output, "O"
    )
    init = naja.SNLParameter.create_binary(
        primitive, "INIT", 1 << input_count, 0
    )
    primitive.setTruthTableFromParameter(output, inputs, init)


def _create_term(primitive, direction, name, width=1):
    if width == 1:
        return naja.SNLScalarTerm.create(primitive, direction, name)
    return naja.SNLBusTerm.create(primitive, direction, width - 1, 0, name)


def _set_role(terms, role, active_level=naja.SNLActiveLevel.NA):
    for term in terms:
        bits = list(term.getBits()) if hasattr(term, "getBits") else [term]
        for bit in bits:
            bit.setRole(role, active_level)


def _construct_carry4(lib):
    primitive = naja.SNLDesign.createPrimitive(lib, "CARRY4")
    inputs = [
        _create_term(primitive, naja.SNLTerm.Direction.Input, "CI"),
        _create_term(primitive, naja.SNLTerm.Direction.Input, "CYINIT"),
        _create_term(primitive, naja.SNLTerm.Direction.Input, "DI", 4),
        _create_term(primitive, naja.SNLTerm.Direction.Input, "S", 4),
    ]
    outputs = [
        _create_term(primitive, naja.SNLTerm.Direction.Output, "CO", 4),
        _create_term(primitive, naja.SNLTerm.Direction.Output, "O", 4),
    ]
    naja.SNLDesign.addCombinatorialArcs(inputs, outputs)


def _construct_dsp48e1(lib):
    primitive = naja.SNLDesign.createPrimitive(lib, "DSP48E1")
    input_specs = (
        ("A", 30),
        ("ACIN", 30),
        ("ALUMODE", 4),
        ("B", 18),
        ("BCIN", 18),
        ("C", 48),
        ("CARRYIN", 1),
        ("CARRYINSEL", 3),
        ("CEM", 1),
        ("CEP", 1),
        ("D", 25),
        ("INMODE", 5),
        ("OPMODE", 7),
        ("PCIN", 48),
        ("RSTM", 1),
        ("RSTP", 1),
    )
    inputs = {
        name: _create_term(
            primitive, naja.SNLTerm.Direction.Input, name, width
        )
        for name, width in input_specs
    }
    clock = _create_term(
        primitive, naja.SNLTerm.Direction.Input, "CLK"
    )
    output = _create_term(
        primitive, naja.SNLTerm.Direction.Output, "P", 48
    )
    naja.SNLDesign.addInputsToClockArcs(list(inputs.values()), clock)
    naja.SNLDesign.addClockToOutputsArcs(clock, output)
    _set_role([clock], naja.SNLTermRole.Clock)
    _set_role(inputs.values(), naja.SNLTermRole.DataInput)
    _set_role([output], naja.SNLTermRole.DataOutput)
    _set_role([inputs["CEM"], inputs["CEP"]], naja.SNLTermRole.Enable)
    _set_role(
        [inputs["RSTM"], inputs["RSTP"]],
        naja.SNLTermRole.SyncReset,
        naja.SNLActiveLevel.High,
    )

    for name in (
        "ACASCREG",
        "ADREG",
        "ALUMODEREG",
        "AREG",
        "BCASCREG",
        "BREG",
        "CARRYINREG",
        "CARRYINSELREG",
        "CREG",
        "DREG",
        "INMODEREG",
        "MREG",
        "OPMODEREG",
        "PREG",
    ):
        naja.SNLParameter.create_decimal(primitive, name, 1)
    naja.SNLParameter.create_string(primitive, "A_INPUT", "DIRECT")
    naja.SNLParameter.create_string(primitive, "B_INPUT", "DIRECT")
    naja.SNLParameter.create_boolean(primitive, "USE_DPORT", False)
    naja.SNLParameter.create_string(primitive, "USE_MULT", "MULTIPLY")
    naja.SNLParameter.create_string(primitive, "USE_SIMD", "ONE48")


def _construct_distributed_ram(lib, name, address_width, data_width):
    primitive = naja.SNLDesign.createPrimitive(lib, name)
    addresses = {}
    data_inputs = {}
    data_outputs = {}
    for suffix in "ABCD":
        addresses[suffix] = _create_term(
            primitive,
            naja.SNLTerm.Direction.Input,
            f"ADDR{suffix}",
            address_width,
        )
        data_inputs[suffix] = _create_term(
            primitive,
            naja.SNLTerm.Direction.Input,
            f"DI{suffix}",
            data_width,
        )
        data_outputs[suffix] = _create_term(
            primitive,
            naja.SNLTerm.Direction.Output,
            f"DO{suffix}",
            data_width,
        )
        naja.SNLDesign.addCombinatorialArcs(
            list(addresses[suffix].getBits()),
            list(data_outputs[suffix].getBits()),
        )

    clock = _create_term(
        primitive, naja.SNLTerm.Direction.Input, "WCLK"
    )
    write_enable = _create_term(
        primitive, naja.SNLTerm.Direction.Input, "WE"
    )
    naja.SNLDesign.addInputsToClockArcs(
        list(data_inputs.values()) + [write_enable], clock
    )
    _set_role([clock], naja.SNLTermRole.Clock)
    _set_role([write_enable], naja.SNLTermRole.MemoryWriteEnable)
    _set_role(addresses.values(), naja.SNLTermRole.MemoryReadAddress)
    _set_role(data_inputs.values(), naja.SNLTermRole.MemoryWriteData)
    _set_role(data_outputs.values(), naja.SNLTermRole.MemoryReadData)
    for suffix in "ABCD":
        naja.SNLParameter.create_binary(
            primitive, f"INIT_{suffix}", 64, 0
        )
    if name == "RAM64M":
        naja.SNLParameter.create_binary(
            primitive, "IS_WCLK_INVERTED", 1, 0
        )


def _construct_block_ram(
    lib,
    name,
    address_width,
    data_width,
    parity_width,
    write_enable_width,
    init_count,
    initp_count,
):
    primitive = naja.SNLDesign.createPrimitive(lib, name)
    a_inputs = []
    b_inputs = []
    a_outputs = []
    b_outputs = []

    for suffix, inputs in (("A", a_inputs), ("B", b_inputs)):
        address = _create_term(
            primitive,
            naja.SNLTerm.Direction.Input,
            "ADDRARDADDR" if suffix == "A" else "ADDRBWRADDR",
            address_width,
        )
        data_input = _create_term(
            primitive,
            naja.SNLTerm.Direction.Input,
            "DIADI" if suffix == "A" else "DIBDI",
            data_width,
        )
        parity_input = _create_term(
            primitive,
            naja.SNLTerm.Direction.Input,
            "DIPADIP" if suffix == "A" else "DIPBDIP",
            parity_width,
        )
        data_output = _create_term(
            primitive,
            naja.SNLTerm.Direction.Output,
            "DOADO" if suffix == "A" else "DOBDO",
            data_width,
        )
        parity_output = _create_term(
            primitive,
            naja.SNLTerm.Direction.Output,
            "DOPADOP" if suffix == "A" else "DOPBDOP",
            parity_width,
        )
        inputs.extend(
            list(address.getBits())
            + list(data_input.getBits())
            + list(parity_input.getBits())
        )
        (a_outputs if suffix == "A" else b_outputs).extend(
            list(data_output.getBits()) + list(parity_output.getBits())
        )
        _set_role([address], naja.SNLTermRole.MemoryReadAddress)
        _set_role(
            [data_input, parity_input], naja.SNLTermRole.MemoryWriteData
        )
        _set_role(
            [data_output, parity_output], naja.SNLTermRole.MemoryReadData
        )

    clock_a = _create_term(
        primitive, naja.SNLTerm.Direction.Input, "CLKARDCLK"
    )
    clock_b = _create_term(
        primitive, naja.SNLTerm.Direction.Input, "CLKBWRCLK"
    )
    controls = {}
    for control_name in (
        "ENARDEN",
        "ENBWREN",
        "REGCEAREGCE",
        "REGCEB",
        "RSTRAMARSTRAM",
        "RSTRAMB",
        "RSTREGARSTREG",
        "RSTREGB",
    ):
        controls[control_name] = _create_term(
            primitive, naja.SNLTerm.Direction.Input, control_name
        )
    wea = _create_term(
        primitive,
        naja.SNLTerm.Direction.Input,
        "WEA",
        write_enable_width,
    )
    webwe = _create_term(
        primitive,
        naja.SNLTerm.Direction.Input,
        "WEBWE",
        write_enable_width * 2,
    )
    a_inputs.extend(
        [
            controls["ENARDEN"],
            controls["REGCEAREGCE"],
            controls["RSTRAMARSTRAM"],
            controls["RSTREGARSTREG"],
        ]
        + list(wea.getBits())
    )
    b_inputs.extend(
        [
            controls["ENBWREN"],
            controls["REGCEB"],
            controls["RSTRAMB"],
            controls["RSTREGB"],
        ]
        + list(webwe.getBits())
    )
    naja.SNLDesign.addInputsToClockArcs(a_inputs, clock_a)
    naja.SNLDesign.addInputsToClockArcs(b_inputs, clock_b)
    naja.SNLDesign.addClockToOutputsArcs(clock_a, a_outputs)
    naja.SNLDesign.addClockToOutputsArcs(clock_b, b_outputs)
    _set_role([clock_a, clock_b], naja.SNLTermRole.Clock)
    _set_role(
        [
            controls["ENARDEN"],
            controls["ENBWREN"],
            controls["REGCEAREGCE"],
            controls["REGCEB"],
        ],
        naja.SNLTermRole.Enable,
    )
    _set_role(
        [
            controls["RSTRAMARSTRAM"],
            controls["RSTRAMB"],
            controls["RSTREGARSTREG"],
            controls["RSTREGB"],
        ],
        naja.SNLTermRole.SyncReset,
        naja.SNLActiveLevel.High,
    )
    _set_role([wea, webwe], naja.SNLTermRole.MemoryWriteEnable)

    naja.SNLParameter.create_decimal(primitive, "DOA_REG", 0)
    naja.SNLParameter.create_decimal(primitive, "DOB_REG", 0)
    naja.SNLParameter.create_binary(
        primitive, "INIT_A", data_width + parity_width, 0
    )
    naja.SNLParameter.create_binary(
        primitive, "INIT_B", data_width + parity_width, 0
    )
    for index in range(init_count):
        naja.SNLParameter.create_binary(
            primitive, f"INIT_{index:02X}", 256, 0
        )
    for index in range(initp_count):
        naja.SNLParameter.create_binary(
            primitive, f"INITP_{index:02X}", 256, 0
        )
    if name == "RAMB36E1":
        naja.SNLParameter.create_string(
            primitive, "RAM_EXTENSION_A", "NONE"
        )
        naja.SNLParameter.create_string(
            primitive, "RAM_EXTENSION_B", "NONE"
        )
    naja.SNLParameter.create_string(primitive, "RAM_MODE", "TDP")
    for suffix in "AB":
        naja.SNLParameter.create_decimal(
            primitive, f"READ_WIDTH_{suffix}", 0
        )
        naja.SNLParameter.create_decimal(
            primitive, f"WRITE_WIDTH_{suffix}", 0
        )
        naja.SNLParameter.create_binary(
            primitive,
            f"SRVAL_{suffix}",
            data_width + parity_width,
            0,
        )
        naja.SNLParameter.create_string(
            primitive, f"WRITE_MODE_{suffix}", "WRITE_FIRST"
        )


def _construct_flip_flop(
    lib,
    name,
    control_name,
    control_role,
    next_state,
    init,
    clear=None,
    preset=None,
):
    primitive = naja.SNLDesign.createPrimitive(lib, name)
    q = naja.SNLScalarTerm.create(
        primitive, naja.SNLTerm.Direction.Output, "Q"
    )
    clock = naja.SNLScalarTerm.create(
        primitive, naja.SNLTerm.Direction.Input, "C"
    )
    enable = naja.SNLScalarTerm.create(
        primitive, naja.SNLTerm.Direction.Input, "CE"
    )
    control = naja.SNLScalarTerm.create(
        primitive, naja.SNLTerm.Direction.Input, control_name
    )
    data = naja.SNLScalarTerm.create(
        primitive, naja.SNLTerm.Direction.Input, "D"
    )
    naja.SNLParameter.create_binary(primitive, "INIT", 1, init)

    naja.SNLDesign.addInputsToClockArcs([data, enable, control], clock)
    naja.SNLDesign.addClockToOutputsArcs(clock, q)
    clock.setRole(naja.SNLTermRole.Clock)
    enable.setRole(naja.SNLTermRole.Enable)
    control.setRole(control_role, naja.SNLActiveLevel.High)
    data.setRole(naja.SNLTermRole.DataInput)
    q.setRole(naja.SNLTermRole.DataOutput)

    state = {"name": "IQ", "next_state": next_state}
    if clear is not None:
        state["clear"] = clear
    if preset is not None:
        state["preset"] = preset
    primitive.setSequentialModel(
        clocked_on="C", states=[state], outputs=[(q, "IQ")]
    )


def constructPrimitives(lib):
    for name in ("IBUF", "OBUF", "BUFG"):
        _construct_gate(lib, name, ["I"], "O", 0b10)
    _construct_gate(lib, "INV", ["I"], "O", 0b01)
    _construct_gate(lib, "XOR2", ["A", "B"], "Y", 0b0110)
    _construct_gate(lib, "AND2", ["A", "B"], "Y", 0b1000)
    for name in ("MUXF7", "MUXF8"):
        _construct_gate(lib, name, ["I0", "I1", "S"], "O", 0xCA)
    for input_count in range(1, 7):
        _construct_lut(lib, input_count)
    _construct_carry4(lib)
    _construct_dsp48e1(lib)
    _construct_distributed_ram(lib, "RAM32M", 5, 2)
    _construct_distributed_ram(lib, "RAM64M", 6, 1)
    _construct_block_ram(lib, "RAMB18E1", 14, 16, 2, 2, 64, 8)
    _construct_block_ram(lib, "RAMB36E1", 16, 32, 4, 4, 128, 16)

    enabled_next_state = "(CE & D) | (!CE & IQ)"
    _construct_flip_flop(
        lib,
        "FDCE",
        "CLR",
        naja.SNLTermRole.AsyncReset,
        enabled_next_state,
        0,
        clear="CLR",
    )
    _construct_flip_flop(
        lib,
        "FDPE",
        "PRE",
        naja.SNLTermRole.AsyncSet,
        enabled_next_state,
        1,
        preset="PRE",
    )
    _construct_flip_flop(
        lib,
        "FDRE",
        "R",
        naja.SNLTermRole.SyncReset,
        "(!R) & ((CE & D) | (!CE & IQ))",
        0,
    )
    _construct_flip_flop(
        lib,
        "FDSE",
        "S",
        naja.SNLTermRole.SyncSet,
        "S | ((CE & D) | (!CE & IQ))",
        0,
    )
