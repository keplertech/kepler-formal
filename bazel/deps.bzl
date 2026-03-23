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
# To update: change the commit, run `bazel fetch @glucose @kissat @naja`
# to verify, then update the sha256 hashes.
_GLUCOSE_COMMIT = "7f887abba7cf13636a5ac2d28653668a20a91b25"
_KISSAT_COMMIT = "8af8e56f174b778aef3aa45af9f739b2a5f492c2"
_NAJA_COMMIT = "62e28342f2cf2eb22772689d5cd43a933cc0cebf"
_NAJA_VERILOG_COMMIT = "a377d7f2644bbbf98ff1ae9e8511d52eba8dd6ca"
_NAJA_IF_COMMIT = "27ea776a0c3022fea0c29ccd14a4fb7daea941f0"
_CPPTRACE_COMMIT = "3db8da80111171c219ab5839905771386bee06b3"
_SLANG_COMMIT = "7d36365b4ee32c615f30b971bf2fa69d7d5debc6"
_GOOGLETEST_COMMIT = "52eb8108c5bdec04579160ae17225d66034bd723"

def _deps_impl(_module_ctx):
    http_archive(
        name = "glucose",
        url = "https://github.com/audemard/glucose/archive/{}.tar.gz".format(_GLUCOSE_COMMIT),
        sha256 = "3033a27047f35653f63559e4f31d664cb8b57a7dcdab9d90233be1d1f52f4eda",
        strip_prefix = "glucose-{}".format(_GLUCOSE_COMMIT),
        build_file = "@@//bazel:glucose.BUILD.bazel",
    )

    http_archive(
        name = "kissat",
        url = "https://github.com/arminbiere/kissat/archive/{}.tar.gz".format(_KISSAT_COMMIT),
        sha256 = "9268b6daaf76ea34ea9da503338beddc5539eb783d1a83a37a7af2a028f3b236",
        strip_prefix = "kissat-{}".format(_KISSAT_COMMIT),
        build_file = "@@//bazel:kissat.BUILD.bazel",
    )

    naja_repo(
        name = "naja",
        naja_url = "https://github.com/nanocoh/naja/archive/{}.tar.gz".format(_NAJA_COMMIT),
        naja_sha256 = "f50d7c67291dfc08bcfb4fea08de642880de3837df013d8461d378457ab54777",
        naja_strip_prefix = "naja-{}".format(_NAJA_COMMIT),
        naja_verilog_url = "https://github.com/najaeda/naja-verilog/archive/{}.tar.gz".format(_NAJA_VERILOG_COMMIT),
        naja_verilog_sha256 = "beffd84e14da4b146d9a0dc4699ebf44596cd437cf8157e3b39430a4080c606d",
        naja_verilog_strip_prefix = "naja-verilog-{}".format(_NAJA_VERILOG_COMMIT),
        naja_if_url = "https://github.com/najaeda/naja-if/archive/{}.tar.gz".format(_NAJA_IF_COMMIT),
        naja_if_sha256 = "bfc0b22855aaf332efb86a5b2801363e36646b6fac816e4bcbed9284f627ab52",
        naja_if_strip_prefix = "naja-if-{}".format(_NAJA_IF_COMMIT),
        cpptrace_url = "https://github.com/jeremy-rifkin/cpptrace/archive/{}.tar.gz".format(_CPPTRACE_COMMIT),
        cpptrace_sha256 = "77d689fd7956ff80351a079d83e86a03865dbbe2433b4559cc6cea50bed77390",
        cpptrace_strip_prefix = "cpptrace-{}".format(_CPPTRACE_COMMIT),
        slang_url = "https://github.com/najaeda/slang/archive/{}.tar.gz".format(_SLANG_COMMIT),
        slang_sha256 = "31e36f3b53be8bff17b6d1cbb10e7460b24a7866a68c3aea157c697f329548f3",
        slang_strip_prefix = "slang-{}".format(_SLANG_COMMIT),
        googletest_url = "https://github.com/google/googletest/archive/{}.tar.gz".format(_GOOGLETEST_COMMIT),
        googletest_sha256 = "745c55415660044610f7fcd3af7a6420d5de16a7dbb9ebfe2e131275676232be",
        googletest_strip_prefix = "googletest-{}".format(_GOOGLETEST_COMMIT),
        build_file = "@@//bazel:naja.BUILD.bazel",
    )

deps = module_extension(
    implementation = _deps_impl,
)
