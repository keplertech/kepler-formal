"""MCP stdio facade for one-shot kepler-formal checks."""

from __future__ import annotations

import argparse
from pathlib import Path

from mcp.server.mcpserver import Context, MCPServer

from . import __version__
from .models import CheckKind, CheckResultPayload
from .resolver import BinaryResolver
from .runner import CheckRunner


def create_server(
    project_root: Path, cache_dir: Path | None = None
) -> MCPServer:
    root = project_root.expanduser().resolve()
    runner = CheckRunner(root, BinaryResolver(cache_dir=cache_dir))
    server = MCPServer(
        name="kepler-formal",
        version=__version__,
        instructions=(
            "Run local one-shot equivalence checks. Relative config paths are "
            f"resolved against the session project root {root}. A check_failed "
            "result is a normal checker verdict, not a server failure."
        ),
    )

    async def invoke(
        kind: CheckKind,
        config_path: str,
        timeout_seconds: float,
        ctx: Context,
    ) -> CheckResultPayload:
        async def report(progress: float, total: float, message: str) -> None:
            await ctx.report_progress(progress=progress, total=total, message=message)

        result = await runner.run(kind, config_path, timeout_seconds, report)
        return result.to_dict()

    @server.tool()
    async def gate_lec(
        config_path: str, timeout_seconds: float, ctx: Context
    ) -> CheckResultPayload:
        """Run gate-level combinational equivalence checking from a YAML config.

        config_path may be absolute or relative to the server's project root.
        The config must select verification: lec (or omit it for the LEC default).
        timeout_seconds is required and bounds this one native checker process.
        """
        return await invoke(CheckKind.GATE_LEC, config_path, timeout_seconds, ctx)

    @server.tool()
    async def gate_sec(
        config_path: str, timeout_seconds: float, ctx: Context
    ) -> CheckResultPayload:
        """Run gate-level sequential equivalence checking from a YAML config.

        config_path may be absolute or relative to the server's project root.
        The config must select verification: sec.
        timeout_seconds is required and bounds this one native checker process.
        """
        return await invoke(CheckKind.GATE_SEC, config_path, timeout_seconds, ctx)

    @server.tool()
    async def rtl_sec(
        config_path: str, timeout_seconds: float, ctx: Context
    ) -> CheckResultPayload:
        """Run RTL or RTL-to-gate sequential equivalence from a YAML config.

        config_path may be absolute or relative to the server's project root.
        The config must select verification: sec.
        timeout_seconds is required and bounds this one native checker process.
        """
        return await invoke(CheckKind.RTL_SEC, config_path, timeout_seconds, ctx)

    return server


# Enables `mcp run .../server.py:mcp`; the installed command below is preferred
# when an explicit project root is needed.
mcp = create_server(Path.cwd())


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run kepler-formal's local stdio MCP server"
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path.cwd(),
        help="base directory for relative tool-call paths (default: current directory)",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="override the binary cache (default: ~/.cache/kepler-formal)",
    )
    args = parser.parse_args()
    project_root = args.project_root.expanduser().resolve()
    if not project_root.is_dir():
        parser.error(f"project root is not a directory: {project_root}")
    create_server(project_root, args.cache_dir).run("stdio")
