"""Repository adapters needed by Naja's native Bazel targets."""

def _find_tool(repository_ctx, tool):
    for prefix in ["/opt/homebrew/opt", "/usr/local/opt"]:
        candidate = repository_ctx.path("{}/{}/bin/{}".format(prefix, tool, tool))
        if candidate.exists:
            return candidate
    return repository_ctx.which(tool)

def _alias_repository_impl(repository_ctx):
    lines = ["package(default_visibility = [\"//visibility:public\"])", ""]
    for name, actual in sorted(repository_ctx.attr.aliases.items()):
        lines.extend([
            "alias(",
            "    name = \"{}\",".format(name),
            "    actual = \"{}\",".format(actual),
            ")",
            "",
        ])
    repository_ctx.file("BUILD.bazel", "\n".join(lines))

alias_repository = repository_rule(
    implementation = _alias_repository_impl,
    attrs = {"aliases": attr.string_dict(mandatory = True)},
)

def _system_tool_repository_impl(repository_ctx):
    tool = repository_ctx.attr.tool
    found = _find_tool(repository_ctx, tool)
    if not found:
        fail("{} not found on PATH; required to build Naja".format(tool))
    repository_ctx.symlink(found, tool)
    repository_ctx.file(
        "BUILD.bazel",
        'exports_files(["{}"], visibility = ["//visibility:public"])\n'.format(tool),
    )

system_tool_repository = repository_rule(
    implementation = _system_tool_repository_impl,
    attrs = {"tool": attr.string(mandatory = True)},
    local = True,
)

def _python_repository_impl(repository_ctx):
    python_config = repository_ctx.which("python3-config")
    if not python_config:
        fail("python3-config not found on PATH; required to build Naja")

    def run(args):
        result = repository_ctx.execute([python_config] + args)
        if result.return_code != 0:
            fail("python3-config {} failed: {}".format(" ".join(args), result.stderr))
        return result.stdout.strip().split(" ")

    include_dirs = [arg[2:] for arg in run(["--includes"]) if arg.startswith("-I")]
    if not include_dirs:
        fail("python3-config --includes returned no include directories")
    for index, include_dir in enumerate(include_dirs):
        repository_ctx.symlink(include_dir, "include{}".format(index))

    # Debian's pyconfig.h forwards to a multiarch include directory that
    # python3-config --includes does not report.
    python = repository_ctx.which("python3")
    if not python:
        fail("python3 not found on PATH; required to build Naja")
    multiarch_result = repository_ctx.execute([
        python,
        "-c",
        "import sysconfig; print(sysconfig.get_config_var('MULTIARCH') or '')",
    ])
    if multiarch_result.return_code != 0:
        fail("python3 sysconfig query failed: {}".format(multiarch_result.stderr))

    platform_include_roots = []
    multiarch = multiarch_result.stdout.strip()
    if multiarch:
        for index, include_dir in enumerate(include_dirs):
            include_path = repository_ctx.path(include_dir)
            platform_root = include_path.dirname.get_child(multiarch)
            platform_include = platform_root.get_child(include_path.basename)
            if platform_include.exists:
                root = "platform_include{}".format(index)
                repository_ctx.symlink(
                    platform_include,
                    "{}/{}/{}".format(root, multiarch, include_path.basename),
                )
                platform_include_roots.append(root)

    ldflags = run(["--ldflags", "--embed"])
    linkopts = [arg for arg in ldflags if arg]
    includes = [
        "include{}".format(index)
        for index in range(len(include_dirs))
    ] + platform_include_roots
    header_globs = [
        "include{}/**/*.h".format(index)
        for index in range(len(include_dirs))
    ] + [
        root + "/**/*.h"
        for root in platform_include_roots
    ]
    repository_ctx.file(
        "BUILD.bazel",
        """\
load("@rules_cc//cc:cc_library.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "headers",
    hdrs = glob({header_globs}),
    includes = {includes},
)

cc_library(
    name = "embed",
    hdrs = glob({header_globs}),
    includes = {includes},
    linkopts = {linkopts},
)
""".format(
            header_globs = repr(header_globs),
            includes = repr(includes),
            linkopts = repr(linkopts),
        ),
    )

python_repository = repository_rule(
    implementation = _python_repository_impl,
    local = True,
)

def _host_prefixes_repository_impl(repository_ctx):
    repository_ctx.file(
        "prefixes.bzl",
        'CMAKE_PREFIX_PATH = ""\n',
    )
    repository_ctx.file("BUILD.bazel", "")

host_prefixes_repository = repository_rule(
    implementation = _host_prefixes_repository_impl,
    local = True,
)

def _naja_version_repository_impl(repository_ctx):
    repository_ctx.file(
        "git_version.bzl",
        'NAJA_GIT_HASH = "{}"\n'.format(repository_ctx.attr.git_hash),
    )
    repository_ctx.file("BUILD.bazel", "")

naja_version_repository = repository_rule(
    implementation = _naja_version_repository_impl,
    attrs = {"git_hash": attr.string(mandatory = True)},
)
