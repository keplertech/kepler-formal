# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Root-owned Z3 repository rule for the embedded XLS C2RTL slice."""

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def _z3_extension_impl(_module_ctx):
    http_archive(
        name = "z3",
        build_file = Label("//bazel:z3.BUILD.bazel"),
        integrity = "sha256-gaAsLGTGTWw98jP1kYa5VieZCtoMTC/JAcnCWnByZyo=",
        strip_prefix = "z3-z3-4.14.1",
        urls = ["https://github.com/Z3Prover/z3/archive/z3-4.14.1.tar.gz"],
    )

z3_extension = module_extension(implementation = _z3_extension_impl)
