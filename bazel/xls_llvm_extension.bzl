# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Root-owned LLVM repository rule for the embedded XLS C2RTL slice."""

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def _llvm_raw_ext_impl(_ctx):
    llvm_commit = "1a5376e062e3f3b99ffce25e24e53182202e06b9"
    llvm_sha256 = "7b5c96a3e90a0c414d2745e1d517a859c9fa018f28a53a94c60510e8a68fc40c"

    http_archive(
        name = "llvm-raw",
        build_file_content = "# empty",
        patches = [
            Label("@com_google_xls//dependency_support/llvm:llvm.patch"),
            Label("@com_google_xls//dependency_support/llvm:zlib-header.patch"),
            Label("@com_google_xls//dependency_support/llvm:run_lit.patch"),
        ],
        patch_args = ["-p1"],
        sha256 = llvm_sha256,
        strip_prefix = "llvm-project-" + llvm_commit,
        urls = ["https://github.com/llvm/llvm-project/archive/{commit}.tar.gz".format(commit = llvm_commit)],
    )

llvm_raw_ext = module_extension(implementation = _llvm_raw_ext_impl)
