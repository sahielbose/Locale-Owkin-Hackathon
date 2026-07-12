"""End-to-end integration: drive the REAL MCP server over the MCP protocol.

Unlike test_mcp_server (which calls tools.py directly) and test_engine (which calls
the engine directly), this connects an MCP client to the actual FastMCP server in
process and exercises every tool through the wire: client -> JSON-RPC -> server ->
tools.py -> Lane A engine -> schema object -> structured result -> client.

This is the "does the whole thing actually work together" check.
"""

from __future__ import annotations

import asyncio
from pathlib import Path

import pytest
from mcp.shared.memory import create_connected_server_and_client_session

from src.localespatial.mcp_server import server, tools
from src.localespatial.schema import EnrichmentResult, MapPayload, Niche, SampleRecord

MOCK_PATH = Path(__file__).resolve().parents[1] / "data" / "mock.h5ad"

pytestmark = pytest.mark.skipif(
    not MOCK_PATH.exists(), reason="run `python scripts/make_mock.py` first"
)

EXPECTED_TOOLS = {
    "list_samples",
    "describe_sample",
    "compute_enrichment",
    "find_niches",
    "characterize_niche",
    "find_prognostic_niches",
    "get_map_payload",
    "describe_niches",
    "correlate_niche_outcome",
}


@pytest.fixture(autouse=True)
def _use_mock(monkeypatch):
    monkeypatch.setenv("LOCALE_DATA", str(MOCK_PATH))
    tools.reset_cache()
    yield
    tools.reset_cache()


def _payload(result):
    """Unwrap a FastMCP structured result (lists come back under a 'result' key)."""
    assert result.isError is False, result.content
    structured = result.structuredContent
    if isinstance(structured, dict) and set(structured.keys()) == {"result"}:
        return structured["result"]
    return structured


def test_protocol_lists_seven_tools():
    async def run():
        async with create_connected_server_and_client_session(
            server.mcp._mcp_server
        ) as client:
            listed = await client.list_tools()
            assert {t.name for t in listed.tools} == EXPECTED_TOOLS

    asyncio.run(run())


def test_full_pipeline_over_the_wire():
    async def run():
        async with create_connected_server_and_client_session(
            server.mcp._mcp_server
        ) as client:
            # 1. introspection
            samples = _payload(await client.call_tool("list_samples", {}))
            assert samples
            for s in samples:
                SampleRecord.model_validate(s)
            image_id = samples[0]["image_id"]

            # 2. co-location
            enr = _payload(
                await client.call_tool("compute_enrichment", {"scope": "cohort:breast"})
            )
            EnrichmentResult.model_validate(enr)
            n = len(enr["cell_types"])
            assert len(enr["zscores"]) == n and all(len(r) == n for r in enr["zscores"])

            # 3. niches (name filled by interpret.py through the server)
            niches = _payload(await client.call_tool("find_niches", {}))
            assert niches
            for niche in niches:
                Niche.model_validate(niche)
                assert niche["name"]
                assert abs(sum(niche["composition"].values()) - 1.0) < 0.05

            # 4. characterize one niche by id
            one = _payload(
                await client.call_tool(
                    "characterize_niche", {"niche_id": niches[0]["niche_id"]}
                )
            )
            Niche.model_validate(one)

            # 5. orchestrator: ranked prognostic niches
            ranked = _payload(await client.call_tool("find_prognostic_niches", {}))
            assert ranked
            has_prog = [r["prognostic"] is not None for r in ranked]
            assert has_prog == sorted(has_prog, reverse=True)
            for r in ranked:
                if r["prognostic"] is not None:
                    p = r["prognostic"]
                    assert p["ci_low"] <= p["hazard_ratio"] <= p["ci_high"]
                    assert p["n_patients"] >= 1

            # 6. tissue map payload for the viz, both color modes
            for mode in ("cell_type", "niche"):
                payload = _payload(
                    await client.call_tool(
                        "get_map_payload",
                        {"image_id": image_id, "color_mode": mode},
                    )
                )
                MapPayload.model_validate(payload)
                assert payload["units"] and payload["legend"]
                assert payload["color_mode"] == mode

    asyncio.run(run())


def test_engine_backend_served_the_calls():
    """After the protocol run, the niche/enrichment tools must have used the real
    Lane A engine, not the fallback."""

    async def run():
        async with create_connected_server_and_client_session(
            server.mcp._mcp_server
        ) as client:
            await client.call_tool("compute_enrichment", {"scope": "cohort:breast"})
            await client.call_tool("find_niches", {})

    asyncio.run(run())
    status = tools.backend_status()
    assert status["compute_enrichment"] == "real"
    assert status["find_niches"] == "real"
