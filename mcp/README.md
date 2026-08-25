# kepler-formal MCP server

This directory is the optional Python add-on for the native checker. Install it
from the same source checkout/tag as the checker:

```sh
python -m pip install ./mcp
kepler-formal-mcp --project-root /absolute/path/to/the/design/project
```

The installed command is a local stdio MCP server with three tools:
`gate_lec`, `gate_sec`, and `rtl_sec`. Each accepts a YAML `config_path` and a
required `timeout_seconds`, spawns one fresh `kepler-formal --config ...`
process, waits for it, and returns a structured result. Relative config paths
are resolved against `--project-root`; absolute paths are used unchanged. While
a check is active, clients that request MCP progress receive periodic
`kepler-formal is still running` updates.

Example client configuration:

```json
{
  "mcpServers": {
    "kepler-formal": {
      "command": "/absolute/path/to/venv/bin/kepler-formal-mcp",
      "args": ["--project-root", "/absolute/path/to/design-project"]
    }
  }
}
```

## Binary resolution

Every tool call resolves the checker in this order:

1. `kepler-formal` on `PATH`.
2. An executable under `~/.cache/kepler-formal/bin/`.
3. The matching asset from the latest GitHub Release. The `.sha256` release
   asset is downloaded and verified before the archive is extracted into the
   cache. There are no silent retries or unverified downloads.

If no build exists for the host platform, or any download/checksum step fails,
the tool returns a specific `setup_error`. The current release workflow
publishes `linux-x86_64`; other platforms should put a locally built binary on
`PATH` until their release asset exists.

## Result taxonomy

`status` is one of `passed`, `check_failed`, `inconclusive`, `setup_error`,
`crash`, or `timeout`. `check_failed` means the checker ran normally and found
a real mismatch; it is returned as ordinary tool data, not an MCP protocol
error. `crash` covers abnormal exits and invalid wrapper invocations. Output is
bounded to stdout/stderr tails so solver logs cannot grow the server without
limit.

## Why native stdio, not Docker

This server deliberately runs the native checker with the host user's
permissions. The agentic editing loop needs immediate access to RTL and
netlists in the working tree, and kepler-formal only reads/checks those design
files—it does not edit them. A container would add volume-mount and path-sync
friction without adding a meaningful isolation benefit for this local,
read-only checker. Container packaging remains appropriate for CI and a future
hosted MCP deployment, but it is not the default local transport.
