"""Include flags for compiling Kepler targets against Bazel-fetched Naja.

Bazel's transitive `includes` are emitted as system include paths.  On macOS a
host-installed Naja under /usr/local/include can otherwise win the short
`#include "DNL.h"` lookup before the fetched @naja headers are considered.

The -iquote prefix is derived from @naja's repo mapping instead of a
hardcoded canonical repo name: as the root module the repo materializes at
external/+deps+naja, but in a downstream consumer it is
external/<module>++deps+naja — a hardcoded path dangles there, quoted
includes fall through to the transitive system paths, and headers resolve
through two physical copies (source tree + cmake install), which #pragma
once cannot deduplicate.
"""

_NAJA_ROOT = Label("@naja//:BUILD").workspace_root

_NAJA_QUOTE_INCLUDE_DIRS = [
    "src/src/dnl",
    "src/src/nl/netlist/core",
    "src/src/nl/netlist/snl",
    "src/src/core",
    "src/src/bne",
    "src/src/metrics",
    "src/src/nl/netlist/decorators",
    "src/src/nl/netlist/pnl",
    "src/src/nl/netlist/visual",
    "src/src/nl/netlist/serialization/capnp",
    "src/src/nl/formats/lefdef",
    "src/src/nl/formats/liberty",
    "src/src/nl/formats/systemverilog/frontend",
    "src/src/nl/formats/verilog/backend",
    "src/src/nl/formats/verilog/frontend",
    "src/src/nl/python/pyloader",
    "src/src/optimization",
    "src/thirdparty/cpptrace/include",
    "src/thirdparty/naja-verilog/src",
    "src/thirdparty/yosys-liberty/src",
    "src/thirdparty/lefdef/src/def/def",
]

NAJA_HEADER_COPTS = [
    # Kepler and Naja's headers require C++20. The root module sets this via
    # .bazelrc (--cxxopt=-std=c++20), but that does not reach a downstream
    # consumer, so set it on the targets themselves.
    "-std=c++20",
] + [
    "-iquote%s/%s" % (_NAJA_ROOT, d)
    for d in _NAJA_QUOTE_INCLUDE_DIRS
]
