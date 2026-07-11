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
