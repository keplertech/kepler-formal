#pragma once

#include <cstdint>

namespace KEPLER_FORMAL {

// Lifetime of the dissolved netlist (DNL). Kepler Formal tears the DNL down
// whenever the netlist it views changes, so every teardown starts a new epoch.
// Caches of per-netlist data keyed by object addresses or iso IDs (truth
// tables of models, model input layouts, term-to-iso maps) include the epoch:
// a design recreated at the address of a destroyed one, or an iso ID reused by
// a rebuilt DNL, never hits a stale entry.
uint64_t dnlEpoch();

// Destroys the DNL singleton and starts a new epoch. Use this instead of
// calling naja::DNL::destroy() directly.
void destroyDnl();

}  // namespace KEPLER_FORMAL
