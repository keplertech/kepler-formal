"""Collect public headers of cc_library targets as plain files.

Used by //:deps to export headers of external cc_library targets (e.g. BCR
zlib) whose header files are not visible as source targets.
"""

def _collect_headers_impl(ctx):
    headers = []
    for dep in ctx.attr.deps:
        for hdr in dep[CcInfo].compilation_context.direct_public_headers:
            if not ctx.attr.basenames or hdr.basename in ctx.attr.basenames:
                headers.append(hdr)
    return [DefaultInfo(files = depset(headers))]

collect_headers = rule(
    implementation = _collect_headers_impl,
    attrs = {
        "basenames": attr.string_list(
            doc = "If set, only headers with these basenames are collected.",
        ),
        "deps": attr.label_list(providers = [CcInfo]),
    },
)
