"""Niche -> outcome association (Lane A). Backed by lifelines.

This is where the tool becomes a biomarker engine: does a niche's per-patient
abundance track overall survival? Cox proportional hazards + Kaplan-Meier, with
FDR correction across niches (rank_prognostic_niches).

The hazard ratio is reported per 1 standard deviation of niche abundance, so it is
interpretable rather than a per-unit extrapolation. On a tiny/underpowered cohort
(the 3-patient mock) the fit is regularized and the ratio is clamped to a finite,
plausible band; a real cohort (~281 patients) needs neither.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from lifelines import CoxPHFitter, KaplanMeierFitter

from ..schema import KMCurve, Prognostic

_REQUIRED = ("niche", "patient_id", "os_month", "os_event")


def niche_outcome(
    adata: AnnData, niche_id: int, high_low_quantile: float = 0.5
) -> Prognostic:
    """Associate one niche's abundance with survival (lifelines Cox + KM).

    Args:
        adata: canonical AnnData with obs['niche'], obs['patient_id'],
            obs['os_month'], obs['os_event'].
        niche_id: the niche whose abundance is tested.
        high_low_quantile: split patients into high/low niche-abundance groups at
            this quantile for the KM curve (default median).

    Returns:
        Prognostic(hazard_ratio, ci_low, ci_high, pvalue, n_patients, km=KMCurve).
    """
    per_patient = _per_patient_abundance(adata, int(niche_id))
    n_patients = int(per_patient.shape[0])
    if (
        n_patients < 2
        or int(per_patient["os_event"].sum()) == 0
        or per_patient["abundance"].nunique() < 2
    ):
        raise ValueError(
            "survival association not viable (need >=2 patients, an event, and "
            "abundance variation)"
        )

    abundance = per_patient["abundance"]
    sd = float(abundance.std(ddof=0))
    per_patient["abundance_z"] = (
        (abundance - abundance.mean()) / sd if sd > 0 else abundance - abundance.mean()
    )

    # Underpowered cohorts (the mock) make an unregularized Cox fit blow up; use a
    # heavier ridge when n is tiny. A real cohort hits penalizer 0.1.
    penalizer = 0.1 if n_patients >= 20 else 1.0
    design = per_patient[["os_month", "os_event", "abundance_z"]].rename(
        columns={"abundance_z": "abundance"}
    )
    cph = CoxPHFitter(penalizer=penalizer)
    cph.fit(design, duration_col="os_month", event_col="os_event")
    row = cph.summary.loc["abundance"]

    hazard_ratio, ci_low, ci_high = _sane_hr(
        float(np.exp(row["coef"])),
        float(np.exp(row["coef lower 95%"])),
        float(np.exp(row["coef upper 95%"])),
    )
    return Prognostic(
        hazard_ratio=hazard_ratio,
        ci_low=ci_low,
        ci_high=ci_high,
        pvalue=float(row["p"]),
        n_patients=n_patients,
        km=_km_curve(per_patient, high_low_quantile),
    )


def rank_prognostic_niches(adata: AnnData) -> list[tuple[int, Prognostic]]:
    """Rank every viable niche by prognostic strength, BH-FDR corrected across niches.

    Computes niche_outcome for each niche, Benjamini-Hochberg corrects the p-values,
    and returns (niche_id, Prognostic) for the viable niches sorted by q-value
    (most significant first). Niches whose survival fit is not viable are omitted.
    """
    ids = sorted({int(n) for n in adata.obs["niche"].to_numpy()})
    viable: dict[int, Prognostic] = {}
    for niche_id in ids:
        try:
            viable[niche_id] = niche_outcome(adata, niche_id)
        except Exception:
            continue
    qvalues = _benjamini_hochberg({k: v.pvalue for k, v in viable.items()})
    return sorted(viable.items(), key=lambda kv: qvalues.get(kv[0], 1.0))


# --- helpers ----------------------------------------------------------------------


def _per_patient_abundance(adata: AnnData, niche_id: int) -> pd.DataFrame:
    """Per-patient niche abundance + survival, dropping patients without survival."""
    obs = adata.obs
    missing = [c for c in _REQUIRED if c not in obs.columns]
    if missing:
        raise ValueError(f"missing obs columns {missing}")
    frame = pd.DataFrame(
        {
            "patient_id": obs["patient_id"].astype(str).to_numpy(),
            "in_niche": (obs["niche"].to_numpy().astype(int) == niche_id).astype(float),
            "os_month": obs["os_month"].to_numpy(dtype=float),
            "os_event": obs["os_event"].to_numpy(),
        }
    )
    grouped = frame.groupby("patient_id").agg(
        abundance=("in_niche", "mean"),
        os_month=("os_month", "first"),
        os_event=("os_event", "first"),
    )
    grouped = grouped[np.isfinite(grouped["os_month"].to_numpy(dtype=float))]
    grouped = grouped[grouped["os_month"] > 0]
    grouped["os_event"] = grouped["os_event"].astype(int)
    return grouped


def _sane_hr(hr: float, ci_low: float, ci_high: float) -> tuple[float, float, float]:
    """Clamp HR + CI to a finite, plausible band and guarantee ci_low <= hr <= ci_high."""

    def clamp(value: float, lo: float, hi: float) -> float:
        if not np.isfinite(value):
            return hi if value > 0 else lo
        return float(min(max(value, lo), hi))

    hr = clamp(hr, 1e-2, 1e2)
    ci_low = clamp(ci_low, 1e-3, hr)
    ci_high = clamp(ci_high, hr, 1e3)
    return hr, ci_low, ci_high


def _km_curve(per_patient: pd.DataFrame, quantile: float) -> KMCurve | None:
    """Kaplan-Meier survival for high vs low niche-abundance patient groups."""
    try:
        threshold = float(per_patient["abundance"].quantile(quantile))
        high = per_patient[per_patient["abundance"] > threshold]
        low = per_patient[per_patient["abundance"] <= threshold]
        if high.empty or low.empty:  # heavy ties: fall back to a strict median split
            median = float(per_patient["abundance"].median())
            high = per_patient[per_patient["abundance"] > median]
            low = per_patient[per_patient["abundance"] <= median]
        if high.empty or low.empty:
            return None
        t_max = float(per_patient["os_month"].max())
        grid = [round(t_max * i / 12.0, 3) for i in range(13)]

        def survival(group: pd.DataFrame) -> list[float]:
            kmf = KaplanMeierFitter().fit(group["os_month"], group["os_event"])
            return [float(v) for v in kmf.survival_function_at_times(grid).to_numpy()]

        return KMCurve(time=grid, high=survival(high), low=survival(low))
    except Exception:
        return None


def _benjamini_hochberg(pvalues: dict[int, float]) -> dict[int, float]:
    items = [(k, v) for k, v in pvalues.items() if v is not None and np.isfinite(v)]
    m = len(items)
    if m == 0:
        return {}
    order = sorted(items, key=lambda kv: kv[1])
    qvalues: dict[int, float] = {}
    running_min = 1.0
    for rank in range(m - 1, -1, -1):
        key, pvalue = order[rank]
        running_min = min(running_min, pvalue * m / (rank + 1))
        qvalues[key] = running_min
    return qvalues


# --- honest outcome bundle: multiplicity + power aware ------------------------------
# The middle fields of correlate_niche_outcome are the entire product: a point estimate
# means nothing without how many hypotheses were tested, the FDR, the selection-aware p,
# the event count, and the minimum HR this many events can even resolve. An agent cannot
# compute these; they ship with every finding, whether or not they flatter us.


def min_detectable_hr(n_events: int, power: float = 0.80, alpha: float = 0.05) -> float:
    """Smallest hazard ratio detectable at `power`, two-sided `alpha`, given `n_events`,
    for a per-SD (standardized) covariate. Schoenfeld: events = (z_a2 + z_pow)^2 / beta^2.
    Returns HR > 1; the protective mirror is 1/HR. Effects inside [1/HR, HR] are below
    what this many events can resolve, so a non-significant one there is uninformative.
    """
    from scipy.stats import norm

    if n_events < 4:
        return float("inf")
    beta = (norm.ppf(1 - alpha / 2) + norm.ppf(power)) / np.sqrt(n_events)
    return float(np.exp(beta))


def _cohort_frame(adata: AnnData):
    """Per-patient niche abundance + survival, tolerant of Basel (PID/OSmonth/event)
    or schema (patient_id/os_month/os_event) naming."""
    obs = adata.obs
    pid = next(c for c in ("PID", "patient_id") if c in obs.columns)
    osm = next(c for c in ("OSmonth", "os_month") if c in obs.columns)
    ev = next(c for c in ("event", "os_event") if c in obs.columns)
    ab = pd.crosstab(obs[pid].astype(str), obs["niche"], normalize="index")
    ab.columns = [int(c) for c in ab.columns]
    covs = [c for c in ("grade", "clinical_type") if c in obs.columns]
    dedup = obs.drop_duplicates(pid).copy()
    dedup[pid] = dedup[pid].astype(str)
    pat = (
        dedup.set_index(pid)
        .loc[ab.index, [osm, ev, *covs]]
        .rename(columns={osm: "OSmonth", ev: "event"})
    )
    for c in covs:
        pat[c] = pat[c].astype(str).replace("nan", "unknown")
    return ab, pat, covs


def cohort_survival(
    adata: AnnData,
    n_perm: int = 1000,
    power: float = 0.80,
    alpha: float = 0.05,
    seed: int = 0,
) -> dict:
    """Full multiplicity-aware survival summary for every niche, computed ONCE.

    Per niche: Cox HR/CI/p (unpenalized MLE, adjusted for grade + clinical_type) and
    BH q. Cohort: event count, the minimum detectable HR at 80% power, and the
    selection-aware p from permuting survival across patients and taking the best of
    all niches n_perm times. Slow (permutation) -> precompute and cache it.
    """
    ab, pat, covs = _cohort_frame(adata)
    covdum = (
        pd.concat(
            [
                pd.get_dummies(pat[c], prefix=c, drop_first=True).astype(float)
                for c in covs
            ],
            axis=1,
        )
        if covs
        else pd.DataFrame(index=pat.index)
    )
    if not covdum.empty:
        covdum = covdum.loc[:, covdum.std(ddof=0) > 0]
    n_events = int(pat["event"].sum())

    def cox(vec, penalizer):
        d = pd.DataFrame(
            {
                "OSmonth": pat["OSmonth"].values,
                "event": pat["event"].astype(int).values,
                "abundance": (vec - vec.mean()) / vec.std(ddof=0),
            },
            index=pat.index,
        )
        d = pd.concat([d, covdum], axis=1)
        r = (
            CoxPHFitter(penalizer=penalizer)
            .fit(d, "OSmonth", "event")
            .summary.loc["abundance"]
        )
        return (
            float(np.exp(r["coef"])),
            float(np.exp(r["coef lower 95%"])),
            float(np.exp(r["coef upper 95%"])),
            float(r["p"]),
        )

    def cox_mle(vec):
        for pen in (0.0, 0.1):
            try:
                return cox(vec, pen)
            except Exception:
                continue
        return (float("nan"), float("nan"), float("nan"), float("nan"))

    niches = {}
    for n in ab.columns:
        hr, lo, hi, p = cox_mle(ab[n])
        niches[int(n)] = {"hazard_ratio": hr, "ci_95": [lo, hi], "p_raw": p}
    q = _benjamini_hochberg({n: niches[n]["p_raw"] for n in niches})
    for n in niches:
        niches[n]["q_fdr"] = float(q.get(n, 1.0))

    # selection-aware p: best-of-all-niches under permuted survival (ridge for stability)
    designs = {
        int(n): pd.concat(
            [
                pd.DataFrame(
                    {"abundance": (ab[n] - ab[n].mean()) / ab[n].std(ddof=0)},
                    index=pat.index,
                ),
                covdum,
            ],
            axis=1,
        )
        for n in ab.columns
    }
    outc = pat[["OSmonth", "event"]].astype({"event": int})

    def best_p(oc):
        best = 1.0
        for n in designs:
            dd = designs[n].copy()
            dd["OSmonth"] = oc["OSmonth"].values
            dd["event"] = oc["event"].values
            try:
                best = min(
                    best,
                    CoxPHFitter(penalizer=0.1)
                    .fit(dd, "OSmonth", "event")
                    .summary.loc["abundance", "p"],
                )
            except Exception:
                pass
        return best

    obs_best = best_p(outc)
    rng = np.random.default_rng(seed)
    arr = outc.values
    ge = sum(
        best_p(
            pd.DataFrame(
                arr[rng.permutation(len(arr))],
                index=pat.index,
                columns=["OSmonth", "event"],
            )
        )
        <= obs_best
        for _ in range(n_perm)
    )
    return {
        "cohort": {
            "n_patients": int(len(pat)),
            "n_events": n_events,
            "n_hypotheses_tested": int(len(ab.columns)),
            "observed_best_p": float(obs_best),
            "p_selection_aware": float((ge + 1) / (n_perm + 1)),
            "min_detectable_hr": min_detectable_hr(n_events, power, alpha),
            "power": power,
            "alpha": alpha,
        },
        "niches": niches,
    }


def correlate_niche_outcome(
    cohort_summary: dict, niche_id: int, alpha: float = 0.05
) -> dict:
    """Assemble the honest outcome bundle for one niche from a cohort_survival summary.

    verdict is 'supported' only if the FDR q survives alpha AND the selection-aware p is
    below alpha; otherwise 'insufficient evidence'. The selection-aware p is a GLOBAL
    statement: if the best of all niches under permuted survival is no better than chance
    (p_selection_aware >= alpha), then no single niche can be supported whatever its own
    q. Gating on q alone would stamp a niche 'supported' while the panel as a whole is
    chance, which contradicts the thesis. The min_detectable_hr and p_selection_aware
    ship regardless so the caller sees WHY."""
    c = cohort_summary["cohort"]
    n = cohort_summary["niches"][int(niche_id)]
    supported = (
        np.isfinite(n["q_fdr"])
        and n["q_fdr"] < alpha
        and float(c.get("p_selection_aware", 1.0)) < alpha
    )
    return {
        "niche_id": int(niche_id),
        "hazard_ratio": round(n["hazard_ratio"], 3),
        "ci_95": [round(n["ci_95"][0], 3), round(n["ci_95"][1], 3)],
        "p_raw": round(n["p_raw"], 4),
        "n_hypotheses_tested": c["n_hypotheses_tested"],
        "q_fdr": round(n["q_fdr"], 4),
        "p_selection_aware": round(c["p_selection_aware"], 4),
        "n_events": c["n_events"],
        "min_detectable_hr": round(c["min_detectable_hr"], 3),
        "verdict": "supported" if supported else "insufficient evidence",
    }
