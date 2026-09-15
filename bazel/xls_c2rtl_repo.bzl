# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

"""Repository rule that overlays a flat BUILD file on the local XLS tree."""

_GTEST_PROD_HEADER = """#ifndef GTEST_INCLUDE_GTEST_GTEST_PROD_H_
#define GTEST_INCLUDE_GTEST_GTEST_PROD_H_

#define FRIEND_TEST(test_case_name, test_name) \\
  friend class test_case_name##_##test_name##_Test

#endif  // GTEST_INCLUDE_GTEST_GTEST_PROD_H_
"""

def _xls_c2rtl_repository_impl(repo_ctx):
    src_root = repo_ctx.path(repo_ctx.attr.src_marker).dirname
    if repo_ctx.attr.src_subdir:
        src_root = src_root.get_child(repo_ctx.attr.src_subdir)
    repo_ctx.watch_tree(src_root)
    dst_root = repo_ctx.path(".")
    copy_result = repo_ctx.execute([
        "/bin/bash",
        "-c",
        "cd \"$1\" && tar cf - . | (cd \"$2\" && tar xf -)",
        "copy_xls",
        str(src_root),
        str(dst_root),
    ])
    if copy_result.return_code != 0:
        fail("failed to copy XLS source tree: " + copy_result.stderr)

    clean_result = repo_ctx.execute([
        "find",
        ".",
        "(",
        "-name",
        "BUILD",
        "-o",
        "-name",
        "BUILD.bazel",
        "-o",
        "-name",
        "MODULE.bazel",
        "-o",
        "-name",
        "WORKSPACE",
        "-o",
        "-name",
        "WORKSPACE.bazel",
        ")",
        "-delete",
    ])
    if clean_result.return_code != 0:
        fail("failed to clean XLS package markers: " + clean_result.stderr)

    # XLS production headers only use FRIEND_TEST; avoid a non-dev GTest edge.
    repo_ctx.file("gtest/gtest_prod.h", _GTEST_PROD_HEADER)
    repo_ctx.file("BUILD.bazel", repo_ctx.read(repo_ctx.attr.build_file))

xls_c2rtl_repository = repository_rule(
    implementation = _xls_c2rtl_repository_impl,
    attrs = {
        "build_file": attr.label(mandatory = True, allow_single_file = True),
        "src_marker": attr.label(mandatory = True, allow_single_file = True),
        "src_subdir": attr.string(default = ""),
    },
)
