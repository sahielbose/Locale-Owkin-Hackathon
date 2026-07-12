"""Locale shared data contract (Pydantic v2).

This is THE coordination artifact between the three lanes:
  Lane A (data + engine) produces these objects,
  Lane B (MCP server) returns them from tool calls,
  Lane C (viz) renders MapPayload.

The models below are reproduced verbatim from CLAUDE.md. Never change a
schema field without telling the whole team.
"""

from typing import Literal

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


# --- Risk layer (impact) ----------------------------------------------------------
# The risk score and its trust verdict are ONE inseparable object. A RiskScore can
# never exist without its RiskEvidence: `evidence` has no default, so Pydantic refuses
# to construct the model without it. The product predicts AND says "do not act on this".
# These are additive; nothing above changes. Do not edit the models above.


class RiskEvidence(BaseModel):
    """The honesty payload that ships with every risk prediction.

    A point estimate is meaningless without: how many hypotheses were tested, the
    FDR, the selection-aware p, the event count, the smallest HR the events can
    resolve, the OUT-OF-FOLD discrimination (never in-sample) with a CI, and
    calibration. The verdict is the single bit that says whether to act.
    """

    n_events: int
    n_patients: int
    n_hypotheses_tested: int
    q_fdr_min: float  # smallest BH q across the tested niches
    p_selection_aware: float  # best-of-all-niches p under permuted survival
    min_detectable_hr: float  # smallest per-SD HR resolvable at this event count
    c_index_cv: float  # cross-validated (out-of-fold) Harrell c-index
    c_index_ci_95: list[float]  # [low, high] bootstrap CI on c_index_cv
    calibration_slope: float  # held-out recalibration slope (1.0 = calibrated)
    verdict: Literal["supported", "insufficient evidence", "not evaluable"]
    verdict_reason: str  # one plain-English sentence


class NicheContribution(BaseModel):
    """One niche's push on a patient's risk. contribution = coefficient * z(abundance)."""

    niche_id: int
    name: str
    abundance: float  # this patient's raw niche fraction
    coefficient: float  # model coefficient (per SD of abundance)
    contribution: float  # coefficient * standardized abundance (adds to risk_score)


class RiskScore(BaseModel):
    """A patient-level risk prediction, inseparable from its evidence.

    `evidence` is required (no default): a RiskScore cannot be constructed without
    its verdict, confidence interval, and power context.
    """

    patient_id: str | None = None
    image_id: str | None = None
    risk_score: float  # linear predictor (sum of the niche contributions)
    risk_percentile: float  # 0-100, relative to the cohort
    risk_group: Literal["low", "intermediate", "high"]  # by cohort tertile
    top_contributing_niches: list[NicheContribution]
    evidence: RiskEvidence


class NicheCoefficient(BaseModel):
    coef: float  # log hazard per SD of standardized abundance
    hr_per_sd: float  # exp(coef)
    ci_95: list[float]  # [low, high] on hr_per_sd
    p: float


class RiskModelCard(BaseModel):
    """Provenance for the fitted risk model, evidence attached."""

    features: list[int]  # niche ids used as features
    covariates_adjusted: list[str]  # e.g. ["grade", "clinical_type"]
    n_train_patients: int
    n_events: int
    cv_folds: int
    c_index_cv: float
    c_index_ci_95: list[float]
    c_index_in_sample: (
        float  # full model scored on all training patients: optimistic, biased high
    )
    c_index_out_of_fold: float  # the honest cross-validated number (equals c_index_cv), named for the report
    optimism_gap: float  # c_index_in_sample - c_index_out_of_fold: how much the in-sample number overstates
    calibration_slope: float
    evidence: RiskEvidence
    coefficients: dict[int, NicheCoefficient]  # niche_id -> coefficient block
