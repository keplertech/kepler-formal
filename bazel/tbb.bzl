# oneTBB is provided hermetically by the cmake-built @onetbb repo on all
# platforms (headers + shared libtbb/libtbbmalloc via CcInfo). No system
# link flags: the libraries travel through Bazel's solib mechanism.
TBB_DEPS = ["@onetbb"]
