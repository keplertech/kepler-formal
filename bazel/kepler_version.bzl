"""Generate CLI build metadata using Bazel's stable workspace status."""

def _kepler_version_impl(ctx):
    output = ctx.actions.declare_file("KeplerVersion.h")
    ctx.actions.run_shell(
        inputs = [ctx.file.template, ctx.info_file],
        outputs = [output],
        arguments = [
            ctx.file.template.path,
            ctx.info_file.path,
            output.path,
        ],
        command = """
set -eu
git_hash=$(sed -n 's/^STABLE_KEPLER_GIT_HASH //p' "$2")
[ -n "$git_hash" ] || git_hash=unknown
sed "s/@KEPLER_GIT_HASH@/$git_hash/" "$1" > "$3"
""",
        mnemonic = "KeplerVersion",
    )
    return [DefaultInfo(files = depset([output]))]

kepler_version = rule(
    implementation = _kepler_version_impl,
    attrs = {
        "template": attr.label(allow_single_file = True, mandatory = True),
    },
)
