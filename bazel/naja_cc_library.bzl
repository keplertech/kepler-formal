"""C++ rule wrapper that keeps fetched Naja headers ahead of host installs."""

load("@rules_cc//cc:cc_library.bzl", _cc_library = "cc_library")
load("//bazel:naja_includes.bzl", "NAJA_HEADER_COPTS")

def cc_library(**kwargs):
    kwargs["copts"] = kwargs.get("copts", []) + NAJA_HEADER_COPTS
    _cc_library(**kwargs)
