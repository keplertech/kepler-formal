# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

"""Validation of default-on latch support with an explicit event contract."""

import ctypes


def build_latch_options(settings, *, mode: str, has_boundaries: bool) -> dict:
    enabled = settings.latch_support
    if not isinstance(enabled, bool):
        raise TypeError("latch_support must be a bool")
    names = (
        "latch_input_changes", "latch_initial_inputs", "latch_initial_storage",
        "latch_workers", "latch_max_waves", "latch_max_states", "latch_max_transactions",
        "latch_max_symbolic_nodes", "latch_max_sat_conflicts", "latch_max_sat_decisions",
    )
    values = {name: getattr(settings, name) for name in names}
    changes = values["latch_input_changes"]
    if changes is not None:
        if not isinstance(changes, str):
            raise TypeError("latch_input_changes must be a string or None")
        if changes not in ("any", "single"):
            raise ValueError("latch_input_changes must be any or single")
    size_max = (1 << (8 * ctypes.sizeof(ctypes.c_size_t))) - 1
    int_max = (1 << (8 * ctypes.sizeof(ctypes.c_int) - 1)) - 1
    unsigned_max = (1 << (8 * ctypes.sizeof(ctypes.c_uint))) - 1
    for name in names[1:]:
        value = values[name]
        if value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{name} must be an integer or None")
        if value < 0 or value > size_max:
            raise ValueError(f"{name} must fit a nonnegative size_t")
        if name in ("latch_initial_inputs", "latch_initial_storage") and value > 1:
            raise ValueError(f"{name} must explicitly be Boolean 0 or 1")
        if name == "latch_workers" and value > int_max:
            raise ValueError("latch_workers must fit a nonnegative int")
        if name.startswith("latch_max_") and value == 0:
            raise ValueError("latch resource limits must be positive integers")
        if name.startswith("latch_max_sat_") and value > unsigned_max:
            raise ValueError("latch SAT limits must fit a positive unsigned int")
    has_tuning = any(value is not None for value in values.values())
    if not enabled:
        if has_tuning:
            raise ValueError("latch event tuning requires latch_support=True")
    elif not has_tuning:
        pass  # Preserve legacy extraction without inventing initialization.
    elif mode != "sec":
        raise ValueError("latch_support is only supported for SEC")
    elif any(values[name] is None for name in names[:3]):
        raise ValueError(
            "latch_support requires explicit latch_input_changes, latch_initial_inputs and latch_initial_storage"
        )
    elif has_boundaries:
        raise ValueError("latch_support requires the complete top interface, not selected leaf boundaries")
    return {"latch_support": enabled, **values}
