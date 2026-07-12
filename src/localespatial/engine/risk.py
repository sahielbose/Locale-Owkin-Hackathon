"""Patient-level risk from niche composition (Lane A, impact layer). Backed by lifelines.

This turns the cohort-level niche->survival findings into a PER-PATIENT risk score:
a multivariate Cox model over the standardized niche-abundance matrix, adjusted for
the same covariates as outcome.py, with an L2 penalizer (few events per covariate).

The non-negotiable design principle: the risk score and its trust verdict are ONE
object. Every score returned here carries a RiskEvidence (its confidence interval,
its power context, its verdict). The model must be able to predict AND say "do not
act on this". So:

  - c_index is OUT-OF-FOLD (5-fold, stratified on event) with a bootstrap CI, never
    in-sample. In-sample concordance is optimistic and would lie.
  - the multiplicity + power honesty fields (n_hypotheses_tested, q_fdr,
    p_selection_aware, min_detectable_hr) are REUSED from engine.outcome, not
    recomputed, so the risk layer cannot silently disagree with the niche layer.
  - the verdict is conservative: "supported" only if the cross-validated c-index CI
    lower bound is clearly above 0.5 AND the selection-aware p is small AND the
    observed effects exceed the minimum HR the events can resolve. On Basel (79
    events, selection-aware p = 0.44, min detectable HR 1.37) the honest verdict is
    "insufficient evidence", and this code is not tuned to force "supported".

Pure functions over AnnData. No MCP code here.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from anndata import AnnData
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import StratifiedKFold

from ..schema import (
    NicheCoefficient,
    NicheContribution,
    RiskEvidence,
    RiskModelCard,
    RiskScore,
)
from .outcome import _cohort_frame, cohort_survival, min_detectable_hr

logger = logging.getLogger("locale.engine.risk")

_DEFAULT_PENALIZER = 1.0  # L2 ridge: 79 events / ~15 params is few events per covariate
_CV_FOLDS = 5
_BOOTSTRAP = 500
_ALPHA = 0.05
_C_INDEX_FLOOR = (
    0.55  # "clearly above 0.5": the CV c-index CI lower bound must clear this
)
_DEFAULT_COVARIATES = ("grade", "clinical_type")


# --- fitted-model state passed between fit and evaluate ----------------------------


@dataclass
class _ModelState:
    covariates: list[str]
    penalizer: float
    cv_folds: int
    seed: int
    n_perm: int


def _finite(value: float, default: float) -> float:
    value = float(value)
    return value if np.isfinite(value) else float(default)


# --- data prep + fitting (shared by fit_risk_model and the CV folds) ----------------


def _prep(adata: AnnData, covariates) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Per-patient RAW niche-abundance matrix + survival + present covariates.

    Reuses engine.outcome._cohort_frame (the same abundance construction the niche
    survival layer uses) so the two layers never diverge. Drops zero-variance niches.
    """
    ab, pat, covs = _cohort_frame(adata)
    covs = [c for c in covariates if c in covs]
    variance = ab.std(ddof=0)
    ab = ab.loc[:, variance[variance > 0].index]
    ab.columns = [int(c) for c in ab.columns]
    pat = pat.copy()
    pat["event"] = pat["event"].astype(int)
    return ab, pat, covs


def _dummies(pat: pd.DataFrame, covs: list[str]) -> pd.DataFrame:
    if not covs:
        return pd.DataFrame(index=pat.index)
    dums = pd.concat(
        [pd.get_dummies(pat[c], prefix=c, drop_first=True).astype(float) for c in covs],
        axis=1,
    )
    if not dums.empty:
        dums = dums.loc[:, dums.std(ddof=0) > 0]
    return dums


def _fit_coeffs(
    ab: pd.DataFrame, pat: pd.DataFrame, covs: list[str], penalizer: float
) -> tuple[dict[int, dict], pd.Series, pd.Series, list[int]]:
    """Fit the multivariate Cox on standardized niche abundances (+ covariates).

    Returns per-niche coefficient blocks plus the standardization mean/std so the
    exact same transform can be applied to held-out patients.
    """
    niches = [int(c) for c in ab.columns]
    mean = ab[niches].mean()
    std = ab[niches].std(ddof=0).replace(0, 1.0)
    z = (ab[niches] - mean) / std
    z.columns = [f"n{c}" for c in niches]
    covdum = _dummies(pat, covs)
    design = pd.DataFrame(
        {"OSmonth": pat["OSmonth"].to_numpy(), "event": pat["event"].to_numpy()},
        index=pat.index,
    )
    design = pd.concat([design, z, covdum], axis=1)
    cph = CoxPHFitter(penalizer=penalizer).fit(design, "OSmonth", "event")
    coeffs: dict[int, dict] = {}
    for c in niches:
        r = cph.summary.loc[f"n{c}"]
        coeffs[c] = {
            "coef": float(r["coef"]),
            "hr_per_sd": float(np.exp(r["coef"])),
            "ci_95": [
                float(np.exp(r["coef lower 95%"])),
                float(np.exp(r["coef upper 95%"])),
            ],
            "p": float(r["p"]),
        }
    return coeffs, mean, std, niches


def _linear_predictor(
    ab: pd.DataFrame, coeffs: dict[int, dict], mean: pd.Series, std: pd.Series
) -> pd.Series:
    """Niche-only linear predictor: sum over niches of coef * standardized abundance.

    Deliberately excludes covariates so the reported risk score decomposes EXACTLY
    into per-niche contributions (they sum to this). Covariates are adjusters in the
    fit, not part of the spatial risk we attribute to niches.
    """
    niches = [c for c in coeffs if c in ab.columns]
    z = (ab[niches] - mean[niches]) / std[niches]
    weights = np.array([coeffs[c]["coef"] for c in niches])
    return pd.Series(z.to_numpy() @ weights, index=ab.index)


# --- public: fit ------------------------------------------------------------------


def fit_risk_model(
    adata: AnnData,
    covariates=_DEFAULT_COVARIATES,
    penalizer: float = _DEFAULT_PENALIZER,
    honesty: dict | None = None,
    cv_folds: int = _CV_FOLDS,
    seed: int = 0,
    n_perm: int = 1000,
) -> RiskModelCard:
    """Fit the multivariate niche-abundance Cox risk model and evaluate it honestly.

    Args:
        adata: canonical AnnData with obs['niche'], patient id, survival, covariates.
        covariates: clinical covariates to adjust for (same as outcome.py).
        penalizer: L2 strength for lifelines CoxPHFitter.
        honesty: an existing engine.outcome.cohort_survival summary to reuse for the
            multiplicity/power fields (skip the slow permutation). Computed if None.
        cv_folds: folds for the out-of-fold c-index.
        seed: RNG seed for CV + bootstrap reproducibility.
        n_perm: permutations for the selection-aware p when honesty is computed here.

    Returns:
        RiskModelCard with coefficients and an attached RiskEvidence. Never raises for
        a degenerate cohort; it returns a card whose evidence.verdict is "not evaluable".
    """
    ab, pat, covs = _prep(adata, covariates)
    n_patients = int(ab.shape[0])
    n_events = int(pat["event"].sum())

    coeffs: dict[int, dict] = {}
    niches: list[int] = []
    if n_patients >= 2 and ab.shape[1] >= 1 and n_events >= 1:
        try:
            coeffs, _, _, niches = _fit_coeffs(ab, pat, covs, penalizer)
        except Exception as exc:
            logger.warning("risk model fit failed (%s); model is not evaluable", exc)

    state = _ModelState(covs, penalizer, cv_folds, seed, n_perm)
    evidence = _evaluate(adata, state, honesty=honesty)

    return RiskModelCard(
        features=list(niches),
        covariates_adjusted=list(covs),
        n_train_patients=n_patients,
        n_events=n_events,
        cv_folds=cv_folds,
        c_index_cv=evidence.c_index_cv,
        c_index_ci_95=evidence.c_index_ci_95,
        calibration_slope=evidence.calibration_slope,
        evidence=evidence,
        coefficients={c: NicheCoefficient(**coeffs[c]) for c in coeffs},
    )


# --- public: score ----------------------------------------------------------------


def _derive_niche_names(adata: AnnData) -> dict[int, str]:
    obs = adata.obs
    if "niche" not in obs.columns:
        return {}
    have_types = "cell_type" in obs.columns
    types = obs["cell_type"].astype(str).to_numpy() if have_types else None
    niche_arr = obs["niche"].to_numpy().astype(int)
    names: dict[int, str] = {}
    for nid in sorted(set(niche_arr.tolist())):
        if have_types:
            sub = types[niche_arr == nid]
            if sub.size:
                vals, counts = np.unique(sub, return_counts=True)
                names[nid] = f"niche {nid} ({vals[counts.argmax()]}-dominant)"
                continue
        names[nid] = f"niche {nid}"
    return names


def _not_evaluable_score(
    evidence: RiskEvidence,
    patient_id: str | None = None,
    image_id: str | None = None,
) -> RiskScore:
    """A RiskScore whose number is a placeholder; the verdict says do not act on it."""
    return RiskScore(
        patient_id=patient_id,
        image_id=image_id,
        risk_score=0.0,
        risk_percentile=50.0,
        risk_group="intermediate",
        top_contributing_niches=[],
        evidence=evidence,
    )


def score_cohort(
    adata: AnnData, model: RiskModelCard, niche_names: dict[int, str] | None = None
) -> list[RiskScore]:
    """Score every patient in the cohort, each RiskScore carrying model.evidence."""
    if not model.coefficients:
        # model could not be fit: return one not-evaluable score per patient
        ab, _, _ = _prep(adata, model.covariates_adjusted)
        return [
            _not_evaluable_score(model.evidence, patient_id=str(p)) for p in ab.index
        ]

    names = niche_names or _derive_niche_names(adata)
    ab, pat, _ = _prep(adata, model.covariates_adjusted)
    coeffs = {c: model.coefficients[c].coef for c in model.coefficients}
    niches = [c for c in model.features if c in ab.columns and c in coeffs]
    if not niches:
        return [
            _not_evaluable_score(model.evidence, patient_id=str(p)) for p in ab.index
        ]

    mean = ab[niches].mean()
    std = ab[niches].std(ddof=0).replace(0, 1.0)
    z = (ab[niches] - mean) / std
    weights = np.array([coeffs[c] for c in niches])
    lp = pd.Series(z.to_numpy() @ weights, index=ab.index)

    t1, t2 = float(lp.quantile(1 / 3)), float(lp.quantile(2 / 3))
    scores: list[RiskScore] = []
    for pid in ab.index:
        value = float(lp[pid])
        percentile = float(100.0 * (lp <= value).mean())
        group = "high" if value >= t2 else ("low" if value < t1 else "intermediate")
        contribs = [
            NicheContribution(
                niche_id=int(c),
                name=names.get(int(c), f"niche {c}"),
                abundance=float(ab.loc[pid, c]),
                coefficient=float(coeffs[c]),
                contribution=float(z.loc[pid, c] * coeffs[c]),
            )
            for c in niches
        ]
        contribs.sort(key=lambda nc: abs(nc.contribution), reverse=True)
        scores.append(
            RiskScore(
                patient_id=str(pid),
                image_id=None,
                risk_score=value,
                risk_percentile=percentile,
                risk_group=group,
                top_contributing_niches=contribs,
                evidence=model.evidence,
            )
        )
    return scores


def score_patient(
    adata: AnnData,
    model: RiskModelCard,
    patient_id: str,
    niche_names: dict[int, str] | None = None,
) -> RiskScore:
    """Risk for one patient (percentile + tertile are relative to the whole cohort)."""
    for score in score_cohort(adata, model, niche_names):
        if score.patient_id == str(patient_id):
            return score
    raise ValueError(f"patient_id {patient_id!r} not found in the cohort")


# --- public: evaluate (the anti-overconfidence guard) ------------------------------


def evaluate_model(
    adata: AnnData, model: RiskModelCard, honesty: dict | None = None
) -> RiskEvidence:
    """Recompute the honest evidence for a fitted model on this data.

    Public wrapper over the internal evaluator; uses the model's covariate set and
    fold count. The default L2 penalizer is used for the fold refits.
    """
    state = _ModelState(
        list(model.covariates_adjusted),
        _DEFAULT_PENALIZER,
        int(model.cv_folds),
        0,
        1000,
    )
    return _evaluate(adata, state, honesty=honesty)


def _honesty_summary(adata: AnnData, n_perm: int, seed: int) -> dict:
    return cohort_survival(adata, n_perm=n_perm, seed=seed)


def _evaluate(adata: AnnData, state: _ModelState, honesty: dict | None) -> RiskEvidence:
    ab, pat, covs = _prep(adata, state.covariates)
    n_patients = int(ab.shape[0])
    n_events = int(pat["event"].sum())
    n_niches = int(ab.shape[1])

    # multiplicity + power fields: REUSE the niche layer's honesty machinery
    try:
        summary = honesty or _honesty_summary(adata, state.n_perm, state.seed)
        cohort = summary["cohort"]
        niche_stats = summary["niches"]
    except Exception as exc:  # honesty could not be computed (degenerate cohort)
        logger.warning("honesty summary failed (%s)", exc)
        cohort, niche_stats = {}, {}

    n_hypotheses = int(cohort.get("n_hypotheses_tested", n_niches))
    p_selection = _finite(cohort.get("p_selection_aware", 1.0), 1.0)
    min_hr = _finite(cohort.get("min_detectable_hr", min_detectable_hr(n_events)), 99.0)
    q_values = [
        v["q_fdr"] for v in niche_stats.values() if np.isfinite(v.get("q_fdr", np.nan))
    ]
    q_fdr_min = _finite(min(q_values) if q_values else 1.0, 1.0)
    effects = [
        max(v["hazard_ratio"], 1.0 / v["hazard_ratio"])
        for v in niche_stats.values()
        if np.isfinite(v.get("hazard_ratio", np.nan)) and v.get("hazard_ratio", 0) > 0
    ]
    max_effect = max(effects) if effects else 1.0

    n_pos = int((pat["event"] == 1).sum())
    n_neg = int((pat["event"] == 0).sum())
    evaluable = (
        n_events >= 4
        and n_niches >= 1
        and n_patients >= 2 * state.cv_folds
        and n_pos >= state.cv_folds
        and n_neg >= state.cv_folds
    )

    c_index_cv = 0.5
    c_index_ci = [0.5, 0.5]
    calibration = 0.0
    if evaluable:
        try:
            oof = _oof_lp(ab, pat, covs, state.penalizer, state.cv_folds, state.seed)
            times = pat["OSmonth"].to_numpy(dtype=float)
            events = pat["event"].to_numpy(dtype=int)
            ok = np.isfinite(oof)
            c_index_cv = float(concordance_index(times[ok], -oof[ok], events[ok]))
            c_index_ci = _bootstrap_c_index(times[ok], oof[ok], events[ok], state.seed)
            calibration = _calibration_slope(times[ok], oof[ok], events[ok])
        except Exception as exc:
            logger.warning("cross-validation failed (%s); not evaluable", exc)
            evaluable = False

    verdict, reason = _verdict(
        evaluable=evaluable,
        n_patients=n_patients,
        n_events=n_events,
        n_niches=n_niches,
        cv_folds=state.cv_folds,
        c_index_ci_low=c_index_ci[0],
        p_selection=p_selection,
        max_effect=max_effect,
        min_hr=min_hr,
    )

    return RiskEvidence(
        n_events=n_events,
        n_patients=n_patients,
        n_hypotheses_tested=n_hypotheses,
        q_fdr_min=round(q_fdr_min, 4),
        p_selection_aware=round(p_selection, 4),
        min_detectable_hr=round(min_hr, 3),
        c_index_cv=round(c_index_cv, 4),
        c_index_ci_95=[round(c_index_ci[0], 4), round(c_index_ci[1], 4)],
        calibration_slope=round(calibration, 3),
        verdict=verdict,
        verdict_reason=reason,
    )


def _oof_lp(
    ab: pd.DataFrame,
    pat: pd.DataFrame,
    covs: list[str],
    penalizer: float,
    folds: int,
    seed: int,
) -> np.ndarray:
    """Out-of-fold niche linear predictor. Each fold refits (incl. standardization) on
    its own training patients, so the held-out prediction never sees its own outcome."""
    events = pat["event"].to_numpy(dtype=int)
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    oof = np.full(len(ab), np.nan)
    for train_idx, test_idx in skf.split(np.zeros(len(ab)), events):
        coeffs, mean, std, _ = _fit_coeffs(
            ab.iloc[train_idx], pat.iloc[train_idx], covs, penalizer
        )
        lp = _linear_predictor(ab.iloc[test_idx], coeffs, mean, std)
        oof[test_idx] = lp.to_numpy()
    return oof


def _bootstrap_c_index(
    times: np.ndarray,
    oof: np.ndarray,
    events: np.ndarray,
    seed: int,
    n_boot: int = _BOOTSTRAP,
) -> list[float]:
    rng = np.random.default_rng(seed)
    n = len(times)
    idx = np.arange(n)
    vals: list[float] = []
    for _ in range(n_boot):
        sample = rng.choice(idx, n, replace=True)
        if events[sample].sum() == 0:
            continue
        try:
            vals.append(
                float(concordance_index(times[sample], -oof[sample], events[sample]))
            )
        except Exception:
            continue
    if len(vals) < 10:
        return [float(min(vals)) if vals else 0.0, float(max(vals)) if vals else 1.0]
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


def _calibration_slope(times: np.ndarray, oof: np.ndarray, events: np.ndarray) -> float:
    """Cox slope of the held-out linear predictor: 1.0 means well calibrated."""
    try:
        frame = pd.DataFrame({"OSmonth": times, "event": events.astype(int), "lp": oof})
        frame = frame[np.isfinite(frame["lp"])]
        if frame["lp"].nunique() < 2 or int(frame["event"].sum()) == 0:
            return 0.0
        slope = float(
            CoxPHFitter(penalizer=0.0)
            .fit(frame, "OSmonth", "event")
            .summary.loc["lp", "coef"]
        )
        return _finite(slope, 0.0)
    except Exception:
        return 0.0


def _verdict(
    *,
    evaluable: bool,
    n_patients: int,
    n_events: int,
    n_niches: int,
    cv_folds: int,
    c_index_ci_low: float,
    p_selection: float,
    max_effect: float,
    min_hr: float,
) -> tuple[str, str]:
    """Conservative verdict. 'supported' requires ALL THREE: the CV c-index CI lower
    bound clearly above 0.5, a small selection-aware p, and observed effects that
    exceed what the events can resolve. Otherwise 'insufficient evidence'."""
    if not evaluable:
        return (
            "not evaluable",
            f"{n_patients} patients / {n_events} events cannot support "
            f"{cv_folds}-fold cross-validation of a {n_niches}-niche model.",
        )

    ci_ok = c_index_ci_low > _C_INDEX_FLOOR
    p_ok = p_selection < _ALPHA
    effect_ok = max_effect > min_hr
    if ci_ok and p_ok and effect_ok:
        return (
            "supported",
            f"cross-validated c-index CI lower bound {c_index_ci_low:.2f} exceeds 0.5, "
            f"selection-aware p = {p_selection:.3f}, and the strongest niche effect "
            f"exceeds the minimum detectable HR of {min_hr:.2f}.",
        )

    reasons: list[str] = []
    if not ci_ok:
        reasons.append(
            f"cross-validated c-index CI includes 0.5 (lower bound {c_index_ci_low:.2f})"
        )
    if not p_ok:
        reasons.append(f"selection-aware p = {p_selection:.2f}")
    if not effect_ok:
        reasons.append(
            f"{n_events} events cannot resolve hazard ratios below {min_hr:.2f}"
        )
    return "insufficient evidence", "; ".join(reasons) + "."
