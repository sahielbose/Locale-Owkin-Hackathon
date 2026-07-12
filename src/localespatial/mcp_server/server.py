"""Locale MCP server (Lane B).

A remote MCP server built on the Python MCP SDK's FastMCP, served over the
Streamable HTTP transport so it can be added to Claude / K Pro as a REMOTE custom
connector (an https .../mcp URL), mirroring Owkin's Pathology Explorer. Every tool
wraps src/locale/mcp_server/tools.py and returns schema objects.

    python -m localespatial.mcp_server.server                    # streamable-http (default, for K Pro)
    LOCALE_TRANSPORT=stdio python -m localespatial.mcp_server.server   # for Claude Desktop (stdio)
    LOCALE_HOST=0.0.0.0 LOCALE_PORT=8000 python -m localespatial.mcp_server.server

Data source: data/locale.h5ad if present, else data/mock.h5ad (override with LOCALE_DATA).
No authentication for now (add a bearer token / OAuth before exposing publicly).
"""

from __future__ import annotations

import logging
import os
import sys

from mcp.server.fastmcp import Context, FastMCP
from pydantic import BaseModel, Field

from ..schema import (
    EnrichmentResult,
    MapPayload,
    Niche,
    RiskModelCard,
    RiskScore,
    SampleRecord,
)
from . import tools

logger = logging.getLogger("locale.mcp.server")

mcp = FastMCP(
    "locale",
    instructions=(
        "Locale reasons about WHERE cells sit in tumor tissue: it finds recurring "
        "cellular niches (for example an immune-excluded tumor core) and links them "
        "to patient survival. Start with list_samples / describe_sample, then "
        "find_prognostic_niches for the full ranked pipeline, and get_map_payload to "
        "render the tissue map. If a cohort is missing or ambiguous the niche tools "
        "will elicit the cohort and niche count before running."
    ),
    host=os.environ.get("LOCALE_HOST", "127.0.0.1"),
    port=int(os.environ.get("LOCALE_PORT", "8000")),
)


class NicheQuery(BaseModel):
    """Structured input elicited when a niche request omits/ambiguates the cohort."""

    cohort: str = Field(
        default="breast", description="Cohort to analyze, e.g. 'breast'"
    )
    n_niches: int = Field(
        default=6, ge=2, le=20, description="How many niches to detect"
    )


async def _resolve_cohort(
    ctx: Context | None, cohort: str | None, n_niches: int | None
) -> tuple[str, int | None]:
    """Elicit cohort + n_niches when the request is under-specified; degrade gracefully."""
    ambiguous = (not cohort) or (cohort not in tools.known_cohorts())
    if ambiguous and ctx is not None:
        try:
            result = await ctx.elicit(
                message=(
                    "Which cohort should I analyze, and how many niches should I look "
                    f"for? Known cohorts: {tools.known_cohorts()}."
                ),
                schema=NicheQuery,
            )
            if getattr(result, "action", None) == "accept" and result.data is not None:
                cohort = result.data.cohort or cohort
                n_niches = result.data.n_niches or n_niches
        except Exception as exc:  # client does not support elicitation, etc.
            logger.info("elicitation unavailable (%s); using defaults", exc)
    if not cohort or cohort not in tools.known_cohorts():
        if cohort:
            logger.info("cohort %r not known; defaulting to %r", cohort, tools.COHORT)
        cohort = tools.COHORT
    return cohort, n_niches


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
async def find_niches(
    cohort: str | None = None,
    n_niches: int | None = None,
    ctx: Context | None = None,
) -> list[Niche]:
    """Find recurring cellular niches across the cohort (composition + marker program).

    If cohort is omitted or ambiguous, elicits the cohort and niche count first.
    """
    cohort, n_niches = await _resolve_cohort(ctx, cohort, n_niches)
    return tools.find_niches(cohort=cohort, n_niches=n_niches)


@mcp.tool()
def characterize_niche(niche_id: int) -> Niche:
    """Composition, marker program, and name for one niche."""
    return tools.characterize_niche(niche_id)


@mcp.tool()
async def find_prognostic_niches(
    cohort: str | None = None,
    patient_subset: list[str] | None = None,
    ctx: Context | None = None,
) -> list[Niche]:
    """Full pipeline: niches + survival association + naming, ranked worst-survival first.

    If cohort is omitted or ambiguous, elicits the cohort (and niche count) first.
    """
    cohort, _ = await _resolve_cohort(ctx, cohort, None)
    return tools.find_prognostic_niches(cohort=cohort, patient_subset=patient_subset)


@mcp.tool()
def get_map_payload(image_id: str, color_mode: str = "cell_type") -> MapPayload:
    """Build the interactive tissue-map payload for one image (cell_type | niche)."""
    return tools.get_map_payload(image_id=image_id, color_mode=color_mode)


@mcp.tool()
def describe_niches() -> list[dict]:
    """List the discovered cellular niches (name, major-class composition, size)."""
    return tools.describe_niches()


@mcp.tool()
def correlate_niche_outcome(niche_id: int) -> dict:
    """Niche abundance vs overall survival, shipped WITH its statistical context.

    Returns hazard_ratio, ci_95, p_raw, n_hypotheses_tested, q_fdr, p_selection_aware,
    n_events, min_detectable_hr (80% power), and a verdict ('supported' |
    'insufficient evidence'). The context ships unasked so the caller cannot mistake a
    point estimate for a finding.
    """
    return tools.correlate_niche_outcome(niche_id)


# --- risk layer (impact) ----------------------------------------------------------
# Every risk tool returns a RiskEvidence inline. A risk_score whose verdict is
# "insufficient evidence" or "not evaluable" MUST be reported to the user as
# EXPLORATORY and NON-ACTIONABLE, and MUST NOT be presented as a clinical prediction.

_RISK_VERDICT_NOTE = (
    "The returned evidence.verdict governs how you may use risk_score. If it is "
    "'supported', report it as a statistically supported association. If it is "
    "'insufficient evidence' or 'not evaluable', you MUST tell the user the score is "
    "exploratory and non-actionable and MUST NOT present it as a clinical prediction; "
    "state evidence.verdict_reason. Never report risk_score without its verdict."
)

_PREDICT_RISK_DESC = (
    "Patient-level risk from spatial niche composition, WITH its trust verdict. "
    "Returns a risk_score, cohort percentile, tertile risk_group, the per-niche "
    "contribution decomposition (which niches drive this patient's risk), and an "
    "inline RiskEvidence (cross-validated c-index + CI, calibration, power context, "
    "verdict).\n\n" + _RISK_VERDICT_NOTE
)
_RANK_RISK_DESC = (
    "Rank cohort patients by risk (highest first); each carries the same "
    "RiskEvidence.\n\n" + _RISK_VERDICT_NOTE
)
_MODEL_CARD_DESC = (
    "The fitted risk model's provenance and honest evaluation: features (niche ids), "
    "covariates adjusted for, training size, the OUT-OF-FOLD c-index with a bootstrap "
    "CI, the calibration slope, the per-niche coefficients, and the RiskEvidence "
    "verdict.\n\n" + _RISK_VERDICT_NOTE
)


@mcp.tool(description=_PREDICT_RISK_DESC)
def predict_risk(
    patient_id: str | None = None, image_id: str | None = None
) -> RiskScore:
    return tools.predict_risk(patient_id=patient_id, image_id=image_id)


@mcp.tool(description=_RANK_RISK_DESC)
def rank_patients_by_risk(
    cohort: str = "breast", top_n: int | None = None
) -> list[RiskScore]:
    return tools.rank_patients_by_risk(cohort=cohort, top_n=top_n)


@mcp.tool(description=_MODEL_CARD_DESC)
def get_risk_model_card() -> RiskModelCard:
    return tools.get_risk_model_card()


def main() -> None:
    transport = os.environ.get("LOCALE_TRANSPORT", "streamable-http")
    host = os.environ.get("LOCALE_HOST", "127.0.0.1")
    port = os.environ.get("LOCALE_PORT", "8000")

    # The stdio transport (used by desktop clients like Claude Desktop) speaks
    # JSON-RPC over stdout, so a single stray byte of logging on stdout corrupts
    # the stream and the client silently drops the connection. Pin ALL logging to
    # stderr. (streamable-http has its own socket, but we keep logs on stderr for
    # consistency.) `force=True` overrides any stdout handler a dependency installed
    # at import time.
    logging.basicConfig(
        level=os.environ.get("LOCALE_LOG_LEVEL", "INFO"),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
        force=True,
    )
    if transport == "stdio":
        logger.info("starting Locale MCP server (transport=stdio, JSON-RPC on stdout)")
    else:
        logger.info(
            "starting Locale MCP server (transport=%s) on http://%s:%s/mcp",
            transport,
            host,
            port,
        )
    mcp.run(transport=transport)


if __name__ == "__main__":
    main()
