"""Local MCP integration for the native kepler-formal checker."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("kepler-formal-mcp")
except PackageNotFoundError:  # source-tree imports before installation
    __version__ = "development"
