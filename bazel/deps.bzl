"""Module extension that fetches third-party dependencies via http_archive.

For BCR compatibility, dependencies are fetched as source archives from
GitHub rather than relying on git submodules. The cmake flow continues
to use submodules directly (thirdparty/).

When a dependency gains native Bazel support or is published to BCR,
replace the corresponding http_archive with a bazel_dep in MODULE.bazel.
"""

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def _naja_repo_impl(repo_ctx):
    """Repository rule that assembles naja from multiple archives.

    Naja has a nested submodule (naja-verilog) that must be placed at
    thirdparty/naja-verilog/. The whole tree is then placed under src/
    to match the overlay BUILD file paths used by cmake_submodule_repo.

    Any BUILD/BUILD.bazel files in the archives are excluded so the
    entire tree is a single Bazel package (required by rules_foreign_cc
    cmake()).
    """
    # Download naja main archive
    repo_ctx.download_and_extract(
        url = repo_ctx.attr.naja_url,
        sha256 = repo_ctx.attr.naja_sha256,
        stripPrefix = repo_ctx.attr.naja_strip_prefix,
        output = "src",
    )

    # Download naja-verilog archive into the expected subdirectory
    repo_ctx.download_and_extract(
        url = repo_ctx.attr.naja_verilog_url,
        sha256 = repo_ctx.attr.naja_verilog_sha256,
        stripPrefix = repo_ctx.attr.naja_verilog_strip_prefix,
        output = "src/thirdparty/naja-verilog",
    )

    # Download naja-if (capnp schemas)
    repo_ctx.download_and_extract(
        url = repo_ctx.attr.naja_if_url,
        sha256 = repo_ctx.attr.naja_if_sha256,
        stripPrefix = repo_ctx.attr.naja_if_strip_prefix,
        output = "src/thirdparty/naja-if",
    )

    # Download cpptrace
    repo_ctx.download_and_extract(
        url = repo_ctx.attr.cpptrace_url,
        sha256 = repo_ctx.attr.cpptrace_sha256,
        stripPrefix = repo_ctx.attr.cpptrace_strip_prefix,
        output = "src/thirdparty/cpptrace",
    )

    # Download slang (SystemVerilog parser)
    repo_ctx.download_and_extract(
        url = repo_ctx.attr.slang_url,
        sha256 = repo_ctx.attr.slang_sha256,
        stripPrefix = repo_ctx.attr.slang_strip_prefix,
        output = "src/thirdparty/slang",
    )

    # Download googletest
    repo_ctx.download_and_extract(
        url = repo_ctx.attr.googletest_url,
        sha256 = repo_ctx.attr.googletest_sha256,
        stripPrefix = repo_ctx.attr.googletest_strip_prefix,
        output = "src/thirdparty/googletest",
    )

    # Remove any BUILD files from archives to avoid package boundary issues
    result = repo_ctx.execute(["find", "src", "(", "-name", "BUILD", "-o", "-name", "BUILD.bazel", ")", "-delete"])
    if result.return_code != 0:
        fail("Failed to clean BUILD files: " + result.stderr)

    # Write the overlay BUILD file
    repo_ctx.file("BUILD.bazel", repo_ctx.read(repo_ctx.attr.build_file))

naja_repo = repository_rule(
    implementation = _naja_repo_impl,
    attrs = {
        "naja_url": attr.string(mandatory = True),
        "naja_sha256": attr.string(mandatory = True),
        "naja_strip_prefix": attr.string(mandatory = True),
        "naja_verilog_url": attr.string(mandatory = True),
        "naja_verilog_sha256": attr.string(mandatory = True),
        "naja_verilog_strip_prefix": attr.string(mandatory = True),
        "naja_if_url": attr.string(mandatory = True),
        "naja_if_sha256": attr.string(mandatory = True),
        "naja_if_strip_prefix": attr.string(mandatory = True),
        "cpptrace_url": attr.string(mandatory = True),
        "cpptrace_sha256": attr.string(mandatory = True),
        "cpptrace_strip_prefix": attr.string(mandatory = True),
        "slang_url": attr.string(mandatory = True),
        "slang_sha256": attr.string(mandatory = True),
        "slang_strip_prefix": attr.string(mandatory = True),
        "googletest_url": attr.string(mandatory = True),
        "googletest_sha256": attr.string(mandatory = True),
        "googletest_strip_prefix": attr.string(mandatory = True),
        "build_file": attr.label(mandatory = True),
    },
)

# Pinned dependency versions (commit SHAs from thirdparty/ submodules).
# To update: change the commit, run `bazel fetch @cadical @glucose @kissat @naja`
# to verify, then update the sha256 hashes.
# Hermetic library dependencies (previously system packages).
# capnproto and oneTBB are built with rules_foreign_cc cmake() so that
# naja's config-mode find_package(CapnProto)/find_package(TBB) get real
# CMake package trees (and the capnp compiler binary); boost and
# FlexLexer.h are header-only.
_CAPNPROTO_VERSION = "1.4.0"
_ONETBB_VERSION = "2022.3.0"
_BOOST_VERSION = "1_89_0"
_FLEX_VERSION = "2.6.4"

_CADICAL_COMMIT = "7b99c07f0bcab5824a5a3ce62c7066554017f641"
_GLUCOSE_COMMIT = "7f887abba7cf13636a5ac2d28653668a20a91b25"
_KISSAT_COMMIT = "8af8e56f174b778aef3aa45af9f739b2a5f492c2"
_NAJA_COMMIT = "e1a649e2fce182d8b7ca5c5c80ab5d04aad3ffa3"
_NAJA_VERILOG_COMMIT = "8a13b5986c765035548775808273d61defcaf738"
_NAJA_IF_COMMIT = "8719bf93fdcd65534c75eb7a8a1f69393f74a75a"
_CPPTRACE_COMMIT = "3db8da80111171c219ab5839905771386bee06b3"
_CPPITERTOOLS_COMMIT = "5a7f4aa357ed9b0ad59823e3d2acd57217d5beaf"
_SLANG_COMMIT = "512c327c209d3043aa98ecfd02d06a1b73fcd5fb"
_GOOGLETEST_COMMIT = "52eb8108c5bdec04579160ae17225d66034bd723"

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

    naja_repo(
        name = "naja",
        naja_url = "https://github.com/najaeda/naja/archive/{}.tar.gz".format(_NAJA_COMMIT),
        naja_sha256 = "bedbb745eb5110278c3de9da9538b3e5ef025a736e4862d491bc8209626d3f45",
        naja_strip_prefix = "naja-{}".format(_NAJA_COMMIT),
        naja_verilog_url = "https://github.com/najaeda/naja-verilog/archive/{}.tar.gz".format(_NAJA_VERILOG_COMMIT),
        naja_verilog_sha256 = "e5caf041d7c8867bb0805b8a182cc4330afc40d0dde74760f4ee42d10b70c9cb",
        naja_verilog_strip_prefix = "naja-verilog-{}".format(_NAJA_VERILOG_COMMIT),
        naja_if_url = "https://github.com/najaeda/naja-if/archive/{}.tar.gz".format(_NAJA_IF_COMMIT),
        naja_if_sha256 = "d9ac71c5021b38bde4c5c1e66462e7e52df9f0ffe8739c8f64f8dfbe2cd0b0ea",
        naja_if_strip_prefix = "naja-if-{}".format(_NAJA_IF_COMMIT),
        cpptrace_url = "https://github.com/jeremy-rifkin/cpptrace/archive/{}.tar.gz".format(_CPPTRACE_COMMIT),
        cpptrace_sha256 = "77d689fd7956ff80351a079d83e86a03865dbbe2433b4559cc6cea50bed77390",
        cpptrace_strip_prefix = "cpptrace-{}".format(_CPPTRACE_COMMIT),
        slang_url = "https://github.com/najaeda/slang/archive/{}.tar.gz".format(_SLANG_COMMIT),
        slang_sha256 = "144054285e246801a579e1365fe50c4d0a04a188025c8cb2bbe2355f653f2cbd",
        slang_strip_prefix = "slang-{}".format(_SLANG_COMMIT),
        googletest_url = "https://github.com/google/googletest/archive/{}.tar.gz".format(_GOOGLETEST_COMMIT),
        googletest_sha256 = "745c55415660044610f7fcd3af7a6420d5de16a7dbb9ebfe2e131275676232be",
        googletest_strip_prefix = "googletest-{}".format(_GOOGLETEST_COMMIT),
        build_file = Label("//bazel:naja.BUILD.bazel"),
    )

deps = module_extension(
    implementation = _deps_impl,
)
