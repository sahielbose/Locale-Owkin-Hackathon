"""Risk layer tests: the honesty invariants that the whole thesis rests on.

  - a RiskScore can never exist without its RiskEvidence,
  - the c-index is OUT-OF-FOLD (differs from the optimistic in-sample value),
  - reducing events flips the verdict to "insufficient evidence",
  - per-niche contributions sum to the linear predictor,
  - a degenerate cohort is "not evaluable", never a bare number.

No real Anthropic calls are made.
"""

from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from lifelines.utils import concordance_index

from src.localespatial.engine import risk
from src.localespatial.mcp_server import interpret
from src.localespatial.schema import RiskEvidence, RiskScore

warnings.filterwarnings("ignore")

# Favorable honesty stub so verdict logic is exercised without the slow permutation.
_STUB = {
    "cohort": {
        "n_patients": 200,
        "n_events": 40,
        "n_hypotheses_tested": 6,
        "p_selection_aware": 0.01,
        "min_detectable_hr": 1.2,
        "power": 0.8,
        "alpha": 0.05,
    },
    "niches": {
        i: {"hazard_ratio": 1.9 if i == 0 else 1.03, "q_fdr": 0.01 if i == 0 else 0.7}
        for i in range(6)
    },
}


def _synth(n_patients=200, n_niches=6, signal=0.7, seed=2, cells=30) -> ad.AnnData:
    """Per-cell synthetic cohort with survival driven by niche 0's abundance."""
    rng = np.random.default_rng(seed)
    ab = rng.dirichlet(np.ones(n_niches), size=n_patients)
    z0 = (ab[:, 0] - ab[:, 0].mean()) / ab[:, 0].std()
    rate = 0.02 * np.exp(signal * z0)
    t_event = rng.exponential(1.0 / rate)
    t_cens = rng.exponential(1.0 / (0.02 * 0.9))
    os_month = np.clip(np.minimum(t_event, t_cens), 1, 200)
    os_event = (t_event <= t_cens).astype(int)
    grade = rng.integers(1, 4, n_patients)
    ctype = rng.choice(["LumA", "LumB", "HER2", "TNBC"], n_patients)
    cols = {
        k: []
        for k in (
            "patient_id",
            "niche",
            "os_month",
            "os_event",
            "grade",
            "clinical_type",
            "cell_type",
        )
    }
    for p in range(n_patients):
        for d in rng.choice(n_niches, size=cells, p=ab[p]):
            cols["patient_id"].append(f"P{p:03d}")
            cols["niche"].append(int(d))
            cols["os_month"].append(float(os_month[p]))
            cols["os_event"].append(int(os_event[p]))
            cols["grade"].append(int(grade[p]))
            cols["clinical_type"].append(str(ctype[p]))
            cols["cell_type"].append(f"c{d % 4}")
    obs = pd.DataFrame(cols)
    obs.index = obs.index.astype(str)
    return ad.AnnData(X=rng.normal(size=(len(obs), 4)).astype("float32"), obs=obs)


def _reduce_events(adata: ad.AnnData, keep_events: int) -> ad.AnnData:
    obs = adata.obs.copy()
    per = adata.obs.groupby("patient_id")["os_event"].first()
    keep = set(list(per[per == 1].index)[:keep_events])
    pid = obs["patient_id"].to_numpy()
    obs["os_event"] = np.array([1 if p in keep else 0 for p in pid], dtype=int)
    return ad.AnnData(X=adata.X.copy(), obs=obs)


@pytest.fixture(scope="module")
def cohort() -> ad.AnnData:
    return _synth()


@pytest.fixture(scope="module")
def model(cohort):
    return risk.fit_risk_model(cohort, honesty=_STUB)


# --- 1. a RiskScore can never exist without its evidence --------------------------


def test_riskscore_requires_evidence():
    import pydantic

    with pytest.raises(pydantic.ValidationError):
        RiskScore(
            patient_id="P1",
            risk_score=0.1,
            risk_percentile=50.0,
            risk_group="high",
            top_contributing_niches=[],
        )  # no evidence -> refused by the schema


def test_every_returned_score_carries_evidence(cohort, model):
    for score in risk.score_cohort(cohort, model):
        assert isinstance(score.evidence, RiskEvidence)
        assert score.evidence.verdict in {
            "supported",
            "insufficient evidence",
            "not evaluable",
        }


# --- 2. the c-index is out-of-fold, not in-sample ---------------------------------


def test_c_index_is_out_of_fold(cohort, model):
    ab, pat, covs = risk._prep(cohort, model.covariates_adjusted)
    coeffs, mean, std, _ = risk._fit_coeffs(ab, pat, covs, risk._DEFAULT_PENALIZER)
    lp_insample = risk._linear_predictor(ab, coeffs, mean, std)
    c_in = concordance_index(pat["OSmonth"], -lp_insample, pat["event"])

    c_cv = model.evidence.c_index_cv
    assert 0.0 < c_cv < 1.0
    assert abs(c_in - c_cv) > 1e-6  # the two must differ (CV is not in-sample)
    # the bootstrap CI is a real interval around the CV estimate
    lo, hi = model.evidence.c_index_ci_95
    assert lo <= c_cv <= hi


# --- 3. reducing events flips the verdict to "insufficient evidence" --------------


def test_verdict_flips_to_insufficient_when_events_reduced(cohort):
    full = risk.fit_risk_model(cohort, honesty=None, n_perm=40, seed=2)
    assert full.evidence.verdict == "supported"

    reduced = risk.fit_risk_model(
        _reduce_events(cohort, 12), honesty=None, n_perm=40, seed=2
    )
    assert reduced.evidence.verdict == "insufficient evidence"
    # the honest reason is about power: fewer events cannot resolve the effect
    assert reduced.evidence.min_detectable_hr > full.evidence.min_detectable_hr


# --- 4. per-niche contributions sum to the linear predictor -----------------------


def test_contributions_sum_to_linear_predictor(cohort, model):
    for pid in ["P000", "P001", "P050", "P199"]:
        score = risk.score_patient(cohort, model, pid)
        total = sum(c.contribution for c in score.top_contributing_niches)
        assert abs(total - score.risk_score) < 1e-9
    # they are sorted by absolute contribution (interpretability payoff)
    contribs = risk.score_patient(cohort, model, "P000").top_contributing_niches
    mags = [abs(c.contribution) for c in contribs]
    assert mags == sorted(mags, reverse=True)


# --- 5. a degenerate cohort is "not evaluable" (never a bare number) --------------


def test_tiny_cohort_is_not_evaluable(cohort):
    tiny = _reduce_events(cohort, 1)  # essentially no events -> cannot cross-validate
    card = risk.fit_risk_model(tiny, honesty=None, n_perm=10, seed=2)
    assert card.evidence.verdict in {"not evaluable", "insufficient evidence"}
    scores = risk.score_cohort(tiny, card)
    assert all(isinstance(s.evidence, RiskEvidence) for s in scores)


# --- interpret: risk verbalization always states the verdict (no API calls) --------


def test_summarize_risk_states_verdict_without_api(cohort, model, monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    score = risk.score_patient(cohort, model, "P000")
    text = interpret.summarize_risk(score)
    assert model.evidence.verdict in text.lower()
    if model.evidence.verdict != "supported":
        assert "non-actionable" in text.lower()


def test_summarize_unsupported_is_never_confident(monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    ev = RiskEvidence(
        n_events=79,
        n_patients=281,
        n_hypotheses_tested=12,
        q_fdr_min=0.3,
        p_selection_aware=0.44,
        min_detectable_hr=1.37,
        c_index_cv=0.52,
        c_index_ci_95=[0.44, 0.60],
        calibration_slope=0.4,
        verdict="insufficient evidence",
        verdict_reason="cross-validated c-index CI includes 0.5; 79 events cannot resolve hazard ratios below 1.37.",
    )
    score = RiskScore(
        patient_id="BASEL001",
        risk_score=0.8,
        risk_percentile=92.0,
        risk_group="high",
        top_contributing_niches=[],
        evidence=ev,
    )
    text = interpret.summarize_risk(score).lower()
    assert "insufficient evidence" in text
    assert "non-actionable" in text or "not a clinical prediction" in text


# --- predict_risk must ERROR on an unknown patient, never fabricate a score --------


def test_resolve_patient_rejects_unknown_id():
    from src.localespatial.mcp_server import tools

    obs = pd.DataFrame(
        {"patient_id": ["P000", "P000", "P001"], "image_id": ["A", "A", "B"]}
    )
    adata = ad.AnnData(X=np.zeros((3, 1), dtype="float32"), obs=obs)
    assert tools._resolve_patient(adata, "P001", None) == "P001"
    with pytest.raises(ValueError):
        tools._resolve_patient(adata, "GHOST", None)  # unknown id must error
    with pytest.raises(ValueError):
        tools._resolve_patient(adata, None, "NO_SUCH_IMAGE")


def _supported_score(patient_id="S1") -> RiskScore:
    ev = RiskEvidence(
        n_events=200,
        n_patients=400,
        n_hypotheses_tested=6,
        q_fdr_min=0.01,
        p_selection_aware=0.01,
        min_detectable_hr=1.2,
        c_index_cv=0.72,
        c_index_ci_95=[0.64, 0.80],
        calibration_slope=1.0,
        verdict="supported",
        verdict_reason="cross-validated c-index CI lower bound 0.64 exceeds 0.5.",
    )
    return RiskScore(
        patient_id=patient_id,
        risk_score=0.5,
        risk_percentile=88.0,
        risk_group="high",
        top_contributing_niches=[],
        evidence=ev,
    )


def _fake_anthropic(text: str):
    class _Block:
        type = "text"

        def __init__(self, t):
            self.text = t

    class _Msg:
        def __init__(self, t):
            self.content = [_Block(t)]

    class _Messages:
        def create(self, **kwargs):
            return _Msg(text)

    class _Client:
        def __init__(self, **kwargs):
            self.messages = _Messages()

    return _Client


def test_summarize_supported_uses_compliant_llm_rephrase(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key")
    import anthropic

    monkeypatch.setattr(
        anthropic, "Anthropic", _fake_anthropic("S1 is high risk; this is supported.")
    )
    text = interpret.summarize_risk(_supported_score())
    assert text == "S1 is high risk; this is supported."


def test_summarize_supported_rejects_rephrase_missing_verdict(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key")
    import anthropic

    monkeypatch.setattr(anthropic, "Anthropic", _fake_anthropic("S1 is high risk."))
    text = interpret.summarize_risk(_supported_score())
    assert "supported" in text.lower()  # fell back to the deterministic line


def test_summarize_non_supported_never_calls_llm(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key")
    import anthropic

    def _boom(**kwargs):
        raise AssertionError("LLM must not be called for a non-supported verdict")

    monkeypatch.setattr(anthropic, "Anthropic", _boom)
    ev = RiskEvidence(
        n_events=79,
        n_patients=281,
        n_hypotheses_tested=12,
        q_fdr_min=0.3,
        p_selection_aware=0.44,
        min_detectable_hr=1.37,
        c_index_cv=0.5,
        c_index_ci_95=[0.42, 0.59],
        calibration_slope=0.1,
        verdict="insufficient evidence",
        verdict_reason="79 events cannot resolve hazard ratios below 1.37.",
    )
    score = RiskScore(
        patient_id="X",
        risk_score=0.1,
        risk_percentile=60.0,
        risk_group="intermediate",
        top_contributing_niches=[],
        evidence=ev,
    )
    text = interpret.summarize_risk(score).lower()
    assert "insufficient evidence" in text and "non-actionable" in text
