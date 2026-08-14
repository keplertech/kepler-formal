"""Module extension that fetches third-party dependencies via http_archive.

For BCR compatibility, dependencies are fetched as source archives from
GitHub rather than relying on git submodules. The cmake flow continues
to use submodules directly (thirdparty/).

When a dependency gains native Bazel support or is published to BCR,
replace the corresponding http_archive with a bazel_dep in MODULE.bazel.
"""

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load(
    "//bazel:naja_repositories.bzl",
    "alias_repository",
    "host_prefixes_repository",
    "naja_version_repository",
    "python_repository",
    "system_tool_repository",
)

# Pinned dependency versions.
# To update: change the commit, run `bazel fetch @cadical @glucose @kissat @naja`
# to verify, then update the sha256 hashes.
# Hermetic library dependencies (previously system packages).
# capnproto and oneTBB retain CMake install trees for //:deps; the native
# Naja targets consume BCR capnp-cpp and the same oneTBB CcInfo directly.
# Boost and FlexLexer.h are header-only.
_CAPNPROTO_VERSION = "1.4.0"
_ONETBB_VERSION = "2022.3.0"
_BOOST_VERSION = "1_89_0"
_FLEX_VERSION = "2.6.4"

_CADICAL_COMMIT = "7b99c07f0bcab5824a5a3ce62c7066554017f641"
_GLUCOSE_COMMIT = "7f887abba7cf13636a5ac2d28653668a20a91b25"
_KISSAT_COMMIT = "8af8e56f174b778aef3aa45af9f739b2a5f492c2"
_NAJA_COMMIT = "c96e1b18bf6e411a9ed5273eac0d4cea1c79c73f"
_NAJA_VERILOG_COMMIT = "5da040bb34f0e4e5bb8d67223b999a0132fb401f"
_CPPITERTOOLS_COMMIT = "5a7f4aa357ed9b0ad59823e3d2acd57217d5beaf"
_NAJA_IF_COMMIT = "099677d9f52c0db11b12c08d03e32543eebc7888"
_SLANG_COMMIT = "512c327c209d3043aa98ecfd02d06a1b73fcd5fb"
_TOMLPLUSPLUS_COMMIT = "30172438cee64926dc41fdd9c11fb3ba5b2ba9de"

def _deps_impl(_module_ctx):
    http_archive(
        name = "capnproto",
        url = "https://capnproto.org/capnproto-c++-{}.tar.gz".format(_CAPNPROTO_VERSION),
        sha256 = "fa02378ad522b318916b9ad928d1372fc9abd43dd1f4f0392e50450f5c87828f",
        strip_prefix = "capnproto-c++-{}".format(_CAPNPROTO_VERSION),
        build_file = Label("//bazel:capnproto.BUILD.bazel"),
    )

    http_archive(
        name = "onetbb",
        url = "https://github.com/uxlfoundation/oneTBB/archive/refs/tags/v{}.tar.gz".format(_ONETBB_VERSION),
        sha256 = "01598a46c1162c27253a0de0236f520fd8ee8166e9ebb84a4243574f88e6e50a",
        strip_prefix = "oneTBB-{}".format(_ONETBB_VERSION),
        build_file = Label("//bazel:onetbb.BUILD.bazel"),
    )

    http_archive(
        name = "boost_headers",
        url = "https://archives.boost.io/release/{}/source/boost_{}.tar.gz".format(_BOOST_VERSION.replace("_", "."), _BOOST_VERSION),
        sha256 = "9de758db755e8330a01d995b0a24d09798048400ac25c03fc5ea9be364b13c93",
        strip_prefix = "boost_{}".format(_BOOST_VERSION),
        build_file = Label("//bazel:boost.BUILD.bazel"),
    )

    http_archive(
        name = "flex_src",
        url = "https://github.com/westes/flex/releases/download/v{v}/flex-{v}.tar.gz".format(v = _FLEX_VERSION),
        sha256 = "e87aae032bf07c26f85ac0ed3250998c37621d95f8bd748b31f15b33c45ee995",
        strip_prefix = "flex-{}".format(_FLEX_VERSION),
        build_file = Label("//bazel:flexlexer.BUILD.bazel"),
    )

    http_archive(
        name = "cppitertools",
        url = "https://github.com/ryanhaining/cppitertools/archive/{}.tar.gz".format(_CPPITERTOOLS_COMMIT),
        sha256 = "d61fdb7be3222c7b6c039daafd1f6ff54d7f2b9edb77240b1c34376a3038972e",
        strip_prefix = "cppitertools-{}".format(_CPPITERTOOLS_COMMIT),
        build_file = Label("//bazel:cppitertools.BUILD.bazel"),
    )

    http_archive(
        name = "cadical",
        url = "https://github.com/arminbiere/cadical/archive/{}.tar.gz".format(_CADICAL_COMMIT),
        sha256 = "d89bad4091f2203980ab30fdac14be874a4aca9b716cbcc132f5c7283b6fd987",
        strip_prefix = "cadical-{}".format(_CADICAL_COMMIT),
        build_file = Label("//bazel:cadical.BUILD.bazel"),
    )

    http_archive(
        name = "glucose",
        url = "https://github.com/audemard/glucose/archive/{}.tar.gz".format(_GLUCOSE_COMMIT),
        sha256 = "3033a27047f35653f63559e4f31d664cb8b57a7dcdab9d90233be1d1f52f4eda",
        strip_prefix = "glucose-{}".format(_GLUCOSE_COMMIT),
        build_file = Label("//bazel:glucose.BUILD.bazel"),
    )

    http_archive(
        name = "kissat",
        url = "https://github.com/arminbiere/kissat/archive/{}.tar.gz".format(_KISSAT_COMMIT),
        sha256 = "9268b6daaf76ea34ea9da503338beddc5539eb783d1a83a37a7af2a028f3b236",
        strip_prefix = "kissat-{}".format(_KISSAT_COMMIT),
        build_file = Label("//bazel:kissat.BUILD.bazel"),
    )

    alias_repository(
        name = "boost",
        aliases = {"headers": "@boost_headers//:boost"},
    )

    alias_repository(
        name = "tbb",
        aliases = {
            "tbb": "@onetbb//:onetbb",
            "tbbmalloc": "@onetbb//:onetbb",
        },
    )

    alias_repository(
        name = "flex_headers",
        aliases = {"headers": "@flex_src//:flexlexer"},
    )

    system_tool_repository(name = "bison_tool", tool = "bison")
    system_tool_repository(name = "flex_tool", tool = "flex")
    system_tool_repository(name = "m4_tool", tool = "m4")
    python_repository(name = "python")
    host_prefixes_repository(name = "slang_host_prefixes")
    naja_version_repository(name = "naja_git_version", git_hash = _NAJA_COMMIT[:8])

    http_archive(
        name = "naja-if",
        url = "https://github.com/najaeda/naja-if/archive/{}.tar.gz".format(_NAJA_IF_COMMIT),
        sha256 = "0a980530cc12a19a7df5702012cb096c0186b96fe93328bfb250415834fadd1a",
        strip_prefix = "naja-if-{}".format(_NAJA_IF_COMMIT),
    )

    http_archive(
        name = "naja-verilog",
        url = "https://github.com/najaeda/naja-verilog/archive/{}.tar.gz".format(_NAJA_VERILOG_COMMIT),
        sha256 = "8a0513378c419afc462ffd59c35c7d0362fee7787ab40ef23b2f0a360df0a9df",
        strip_prefix = "naja-verilog-{}".format(_NAJA_VERILOG_COMMIT),
        patch_args = ["-p0", "-f"],
        patches = [Label("//bazel:naja_verilog_bazel9.patch")],
    )

    http_archive(
        name = "tomlplusplus",
        url = "https://github.com/marzer/tomlplusplus/archive/{}.tar.gz".format(_TOMLPLUSPLUS_COMMIT),
        sha256 = "291254ffe7f2433f90deef878d0d9335534a350a958ea23ecf511b7b2277bf7f",
        strip_prefix = "tomlplusplus-{}".format(_TOMLPLUSPLUS_COMMIT),
        build_file = Label("//bazel:tomlplusplus.BUILD.bazel"),
    )

    http_archive(
        name = "slang",
        url = "https://github.com/najaeda/slang/archive/{}.tar.gz".format(_SLANG_COMMIT),
        sha256 = "144054285e246801a579e1365fe50c4d0a04a188025c8cb2bbe2355f653f2cbd",
        strip_prefix = "slang-{}".format(_SLANG_COMMIT),
        build_file = Label("//bazel:slang.BUILD.bazel"),
    )

    http_archive(
        name = "naja",
        url = "https://github.com/nanocoh/naja/archive/{}.tar.gz".format(_NAJA_COMMIT),
        sha256 = "055337bc81b41cf6f10720950d99afed7c4ae7c328f42ec28d9185d392f33f13",
        strip_prefix = "naja-{}".format(_NAJA_COMMIT),
        build_file = Label("//bazel:naja.BUILD.bazel"),
        patch_args = ["-p0", "-f"],
        patches = [Label("//bazel:naja_ff_scan.patch")],
    )

deps = module_extension(
    implementation = _deps_impl,
)
