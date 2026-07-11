"""Locale MCP server (Lane B).

A remote/HTTP MCP server built on the Python MCP SDK's FastMCP. It exposes the
tools listed in CLAUDE.md, each wrapping src/locale/mcp_server/tools.py and
returning schema objects. Wire it into K Pro / Claude as a custom connector
(remote MCP URL), mirroring Owkin's Pathology Explorer.

    python -m src.locale.mcp_server.server                 # streamable-http (default)
    LOCALE_TRANSPORT=stdio python -m src.locale.mcp_server.server
    LOCALE_HOST=0.0.0.0 LOCALE_PORT=8000 python -m src.locale.mcp_server.server

Data source: data/mock.h5ad by default; set LOCALE_DATA=data/locale.h5ad for real data.
"""

from __future__ import annotations

import os

from mcp.server.fastmcp import FastMCP

from ..schema import EnrichmentResult, MapPayload, Niche, SampleRecord
from . import tools

mcp = FastMCP(
    "locale",
    instructions=(
        "Locale reasons about WHERE cells sit in tumor tissue: it finds recurring "
        "cellular niches (for example an immune-excluded tumor core) and links them "
        "to patient survival. Start with list_samples / describe_sample, then "
        "find_prognostic_niches for the full pipeline, and get_map_payload to render "
        "the tissue map. Underspecified requests should elicit cohort / cell types / "
        "n_niches before running."
    ),
    host=os.environ.get("LOCALE_HOST", "127.0.0.1"),
    port=int(os.environ.get("LOCALE_PORT", "8000")),
)


@mcp.tool()
def list_samples() -> list[SampleRecord]:
    """List every sample (one record per IMC image) in the loaded cohort."""
    return tools.list_samples()


@mcp.tool()
def describe_sample(
    image_id: str | None = None, cohort: str | None = None
) -> SampleRecord:
    """Describe one image (pass image_id) or the whole cohort (pass neither)."""
    return tools.describe_sample(image_id=image_id, cohort=cohort)


@mcp.tool()
def compute_enrichment(scope: str = "cohort:breast") -> EnrichmentResult:
    """Cell-type co-location matrix (who is next to whom) for the given scope."""
    return tools.compute_enrichment(scope)


@mcp.tool()
def find_niches(cohort: str = "breast", n_niches: int | None = None) -> list[Niche]:
    """Find recurring cellular niches across the cohort (composition + marker program)."""
    return tools.find_niches(cohort=cohort, n_niches=n_niches)


@mcp.tool()
def characterize_niche(niche_id: int) -> Niche:
    """Composition, marker program, and name for one niche."""
    return tools.characterize_niche(niche_id)


@mcp.tool()
def find_prognostic_niches(
    cohort: str = "breast", patient_subset: list[str] | None = None
) -> list[Niche]:
    """Full pipeline: niches + survival association + naming, ranked by prognostic strength."""
    return tools.find_prognostic_niches(cohort=cohort, patient_subset=patient_subset)


@mcp.tool()
def get_map_payload(image_id: str, color_mode: str = "cell_type") -> MapPayload:
    """Build the interactive tissue-map payload for one image (cell_type | niche)."""
    return tools.get_map_payload(image_id=image_id, color_mode=color_mode)


def main() -> None:
    transport = os.environ.get("LOCALE_TRANSPORT", "streamable-http")
    mcp.run(transport=transport)


if __name__ == "__main__":
    main()
