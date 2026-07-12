"""Locale shared data contract (Pydantic v2).

This is THE coordination artifact between the three lanes:
  Lane A (data + engine) produces these objects,
  Lane B (MCP server) returns them from tool calls,
  Lane C (viz) renders MapPayload.

The models below are reproduced verbatim from CLAUDE.md. Never change a
schema field without telling the whole team.
"""

from pydantic import BaseModel


class SampleRecord(BaseModel):
    cohort: str
    patient_id: str | None = None
    image_id: str | None = None
    n_cells: int
    cell_types: list[str]
    has_survival: bool


class EnrichmentResult(BaseModel):
    scope: str  # e.g. "cohort:breast" or "image:<id>"
    cell_types: list[str]
    zscores: list[list[float]]  # cell_type x cell_type
    pvalues: list[list[float]]


class KMCurve(BaseModel):
    time: list[float]
    high: list[float]  # survival prob, high niche-abundance group
    low: list[float]  # survival prob, low group


class Prognostic(BaseModel):
    hazard_ratio: float
    ci_low: float
    ci_high: float
    pvalue: float
    n_patients: int
    km: KMCurve | None = None


class Niche(BaseModel):
    niche_id: int
    name: str  # human-readable, filled by interpret.py
    composition: dict[str, float]  # cell_type -> fraction
    marker_program: list[str]  # top enriched markers
    prognostic: Prognostic | None = None


class MapUnit(BaseModel):
    x: float
    y: float
    cell_type: str
    niche_id: int | None = None


class MapPayload(BaseModel):
    units: list[MapUnit]
    legend: dict[str, str]  # label -> hex color
    color_mode: str  # "cell_type" | "niche"
    image_id: str | None = None
