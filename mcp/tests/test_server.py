from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from kepler_formal_mcp.server import create_server


class ServerSchemaTests(unittest.IsolatedAsyncioTestCase):
    async def test_tools_have_required_timeouts_and_structured_results(self):
        with tempfile.TemporaryDirectory() as temporary:
            server = create_server(Path(temporary))
            tools = await server.list_tools()
            self.assertEqual(
                {tool.name for tool in tools}, {"gate_lec", "gate_sec", "rtl_sec"}
            )
            for tool in tools:
                self.assertEqual(
                    set(tool.input_schema["required"]),
                    {"config_path", "timeout_seconds"},
                )
                self.assertIn("status", tool.output_schema["properties"])
                self.assertIn(
                    "check_failed",
                    tool.output_schema["properties"]["status"]["enum"],
                )

    async def test_wrapper_input_failure_is_structured_not_protocol_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            server = create_server(Path(temporary))
            result = await server.call_tool(
                "gate_lec",
                {"config_path": "missing.yaml", "timeout_seconds": 1},
            )
            payload = result.structured_content
            self.assertIsNotNone(payload)
            self.assertEqual(payload["status"], "crash")
            self.assertEqual(payload["error_kind"], "crash")


if __name__ == "__main__":
    unittest.main()
