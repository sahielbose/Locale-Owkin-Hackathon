"""MCP tool functions (Lane B).

Each function returns a schema object from src/locale/schema.py. The core pattern
for every tool is graceful degradation:

    try the real engine function (Lane A);
    if it raises NotImplementedError (not finished) or ANY exception,
    fall back to a value derived from the loaded AnnData,

so every tool ALWAYS returns a valid, schema-correct object today on data/mock.h5ad
and AUTO-UPGRADES to the real analysis as Lane A lands each engine function. The
path taken (real vs fallback) is logged and recorded in backend_status().

Engine modules are imported LAZILY inside each call so that an in-progress edit to
Lane A (a syntax error, a changed signature) can never stop this server from
importing or serving; it just routes that one tool to its fallback.

Scope: this file only WRAPS the engine and reads the AnnData. All real analysis
lives in src/locale/engine/. The one bit of computation here is a clearly labeled
TEMPORARY adjacency z-score used by compute_enrichment until engine.enrichment lands.
"""

from __future__ import annotations

import functools
import json
import logging
import os
from pathlib import Path
from typing import Callable, TypeVar

import anndata as ad
import numpy as np

from ..schema import (
    EnrichmentResult,
    MapPayload,
    Niche,
    Prognostic,
    RiskModelCard,
    RiskScore,
    SampleRecord,
)
from ..viz.payload import build_map_payload
from . import interpret

logger = logging.getLogger("locale.mcp.tools")

COHORT = "breast"
_TOP_MARKERS = 6
_REPO_ROOT = Path(__file__).resolve().parents[3]
_MOCK = _REPO_ROOT / "data" / "mock.h5ad"
_REAL = _REPO_ROOT / "data" / "locale.h5ad"
_FINDINGS = (
    _REPO_ROOT / "demo" / "findings.json"
)  # precomputed demo path (scripts/precompute_findings.py)

# Records which backend served each tool most recently ("real" | "fallback").
_BACKEND: dict[str, str] = {}

T = TypeVar("T")


# --- data loading -----------------------------------------------------------------


def data_path() -> Path:
    """Resolve the AnnData to serve.

    Priority: $LOCALE_DATA, then data/locale.h5ad (the real object, once Lane A
    shares it), then the committed data/mock.h5ad.
    """
    env = os.environ.get("LOCALE_DATA")
    if env:
        return Path(env)
    if _REAL.exists():
        return _REAL
    return _MOCK


# Dataset-native obs columns -> the canonical schema names the tools/engine expect.
# The real Basel object (data/basel_niched.h5ad) carries PID/core/OSmonth/event; the
# committed mock already uses the canonical names, so this is a no-op there. Done in
# memory on load only; the .h5ad on disk is never modified.
_OBS_ALIASES = {
    "image_id": ("core",),
    "patient_id": ("PID",),
    "os_month": ("OSmonth",),
    "os_event": ("event",),
}


def _normalize_obs(adata: ad.AnnData) -> ad.AnnData:
    for canonical, sources in _OBS_ALIASES.items():
        if canonical in adata.obs.columns:
            continue
        for src in sources:
            if src in adata.obs.columns:
                adata.obs[canonical] = adata.obs[src].to_numpy()
                break
    return adata


@functools.lru_cache(maxsize=1)
def _load() -> ad.AnnData:
    path = data_path()
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `python scripts/make_mock.py` or set LOCALE_DATA."
        )
    logger.info("loaded AnnData from %s", path)
    return _normalize_obs(ad.read_h5ad(path))


def reset_cache() -> None:
    """Clear the cached AnnData + risk model (used by tests that switch LOCALE_DATA)."""
    _load.cache_clear()
    _risk_model.cache_clear()


def backend_status() -> dict[str, str]:
    """Return the most recent backend ("real" | "fallback") used per tool."""
    return dict(_BACKEND)


def _run(tool: str, engine: Callable[[], T], fallback: Callable[[], T]) -> T:
    """Try the engine path; on NotImplementedError or any error, use the fallback.

    Records and logs which path served the call. Never raises for engine reasons.
    """
    try:
        result = engine()
        _BACKEND[tool] = "real"
        logger.info("tool %s served by REAL engine", tool)
        return result
    except NotImplementedError:
        logger.info("tool %s: engine not implemented, using FALLBACK", tool)
    except Exception as exc:  # engine present but broke (mid-edit, bad data, etc.)
        logger.warning("tool %s: engine error (%s), using FALLBACK", tool, exc)
    _BACKEND[tool] = "fallback"
    return fallback()


# --- shared derivations -----------------------------------------------------------


def _has_survival(obs) -> bool:
    return "os_month" in obs and bool(
        np.isfinite(obs["os_month"].to_numpy(dtype=float)).any()
    )


def known_cohorts() -> list[str]:
    """Cohorts this server can answer about (single-cohort dataset for now)."""
    return [COHORT]


def _niche_ids(adata: ad.AnnData) -> list[int]:
    return sorted({int(n) for n in adata.obs["niche"].to_numpy()})


def _composition(adata: ad.AnnData, mask: np.ndarray) -> dict[str, float]:
    types = adata.obs["cell_type"].astype(str).to_numpy()[mask]
    if types.size == 0:
        return {}
    values, counts = np.unique(types, return_counts=True)
    total = counts.sum()
    return {str(v): round(float(c) / total, 4) for v, c in zip(values, counts)}


def _marker_program(
    adata: ad.AnnData, mask: np.ndarray, k: int = _TOP_MARKERS
) -> list[str]:
    """Top-k markers by mean (z-scored) intensity inside the mask."""
    x = adata.X[mask]
    x = np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)
    if x.shape[0] == 0:
        return []
    means = x.mean(axis=0)
    order = np.argsort(means)[::-1][:k]
    markers = list(adata.var_names)
    return [str(markers[i]) for i in order]


def _emergency_niches(adata: ad.AnnData) -> np.ndarray:
    """Last-resort niche labels when there is no engine AND no precomputed obs['niche'].

    Crude deterministic 2x2 spatial quadrant binning per image (0..3). This should
    almost never run: the mock ships obs['niche'], and real data arrives with the
    engine. It exists only so the niche tools never crash on an unlabeled object.
    """
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    images = adata.obs["image_id"].astype(str).to_numpy()
    labels = np.zeros(coords.shape[0], dtype=int)
    for img in np.unique(images):
        sel = np.where(images == img)[0]
        xs, ys = coords[sel, 0], coords[sel, 1]
        mx, my = np.median(xs), np.median(ys)
        labels[sel] = (xs > mx).astype(int) + 2 * (ys > my).astype(int)
    logger.warning("using crude quadrant niche fallback (no engine, no obs['niche'])")
    return labels


def _labeled_adata(adata: ad.AnnData, n_niches: int | None) -> tuple[ad.AnnData, str]:
    """Return (adata_with_obs['niche'], backend). Try the engine, else precomputed, else crude."""
    # When the object already carries the real engine's niche labeling (the Basel run
    # writes obs['niche'] = the 12 discovered niches) and the caller did not ask for a
    # specific k, serve those. Re-clustering 755k cells on every call is slow and would
    # diverge from the precomputed niche world correlate_niche_outcome/describe_niches use.
    if n_niches is None and "niche" in adata.obs:
        return adata, "real"
    try:
        from ..engine.niches import find_niches as engine_find_niches

        n = n_niches if n_niches is not None else 6
        labeled = engine_find_niches(adata, n_niches=n)
        if "niche" in labeled.obs:
            return labeled, "real"
    except NotImplementedError:
        pass
    except Exception as exc:
        logger.warning("engine.find_niches error (%s), using precomputed/crude", exc)

    if "niche" in adata.obs:
        return adata, "fallback"
    adata = adata.copy()
    adata.obs["niche"] = _emergency_niches(adata)
    return adata, "fallback"


def _characterize(adata: ad.AnnData, niche_id: int, mask: np.ndarray) -> Niche:
    """Composition + marker program + name for one niche (engine, else derived)."""

    def engine() -> Niche:
        from ..engine.characterize import characterize_niche as fn

        niche = fn(adata, niche_id)
        if not niche.name:
            niche.name = interpret.name_niche(
                niche.composition, niche.marker_program, niche_id
            )
        return niche

    def fallback() -> Niche:
        composition = _composition(adata, mask)
        marker_program = _marker_program(adata, mask)
        name = interpret.name_niche(composition, marker_program, niche_id)
        return Niche(
            niche_id=niche_id,
            name=name,
            composition=composition,
            marker_program=marker_program,
            prognostic=None,
        )

    return _run(f"characterize_niche[{niche_id}]", engine, fallback)


# --- sample introspection ---------------------------------------------------------


def list_samples() -> list[SampleRecord]:
    """List every sample (one record per IMC image) in the loaded cohort."""
    adata = _load()
    obs = adata.obs
    records: list[SampleRecord] = []
    for image_id in sorted(set(obs["image_id"].astype(str))):
        sub = obs[obs["image_id"].astype(str) == image_id]
        patient = str(sub["patient_id"].iloc[0]) if "patient_id" in sub else None
        records.append(
            SampleRecord(
                cohort=COHORT,
                patient_id=patient,
                image_id=image_id,
                n_cells=int(sub.shape[0]),
                cell_types=sorted(set(sub["cell_type"].astype(str))),
                has_survival=_has_survival(sub),
            )
        )
    _BACKEND["list_samples"] = "real"  # pure read, no engine involved
    return records


def describe_sample(
    image_id: str | None = None, cohort: str | None = None
) -> SampleRecord:
    """Describe one image (pass image_id) or the whole cohort (pass neither)."""
    adata = _load()
    obs = adata.obs
    _BACKEND["describe_sample"] = "real"
    if image_id is not None:
        sub = obs[obs["image_id"].astype(str) == str(image_id)]
        if sub.shape[0] == 0:
            raise ValueError(f"image_id {image_id!r} not found")
        patient = str(sub["patient_id"].iloc[0]) if "patient_id" in sub else None
        return SampleRecord(
            cohort=COHORT,
            patient_id=patient,
            image_id=str(image_id),
            n_cells=int(sub.shape[0]),
            cell_types=sorted(set(sub["cell_type"].astype(str))),
            has_survival=_has_survival(sub),
        )
    return SampleRecord(
        cohort=cohort or COHORT,
        patient_id=None,
        image_id=None,
        n_cells=int(obs.shape[0]),
        cell_types=sorted(set(obs["cell_type"].astype(str))),
        has_survival=_has_survival(obs),
    )


# --- neighborhood enrichment ------------------------------------------------------


def compute_enrichment(scope: str = f"cohort:{COHORT}") -> EnrichmentResult:
    """Cell-type co-location matrix. Uses engine.compute_enrichment when available."""
    adata = _load()

    def engine() -> EnrichmentResult:
        from ..engine.enrichment import compute_enrichment as fn

        return fn(adata, scope)

    return _run("compute_enrichment", engine, lambda: _mock_enrichment(adata, scope))


def _mock_enrichment(
    adata: ad.AnnData, scope: str, n_perm: int = 100, k: int = 6
) -> EnrichmentResult:
    """TEMPORARY stopgap co-location, replaced by engine.compute_enrichment (squidpy).

    Per-image kNN graph, cell-type-pair adjacency counts, z-scored against a
    within-image label-permutation null. Deterministic (fixed seed).
    """
    obs = adata.obs
    cell_types = sorted(set(obs["cell_type"].astype(str)))
    idx = {ct: i for i, ct in enumerate(cell_types)}
    labels = obs["cell_type"].astype(str).map(idx).to_numpy()
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    images = obs["image_id"].astype(str).to_numpy()
    n_types = len(cell_types)

    def adjacency_counts(lab: np.ndarray) -> np.ndarray:
        counts = np.zeros((n_types, n_types), dtype=float)
        for img in np.unique(images):
            sel = np.where(images == img)[0]
            if sel.size < 2:
                continue
            pts = coords[sel]
            kk = min(k, sel.size - 1)
            d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1)
            neigh = np.argsort(d2, axis=1)[:, 1 : kk + 1]
            for local_i, nbrs in enumerate(neigh):
                a = lab[sel[local_i]]
                for j in nbrs:
                    counts[a, lab[sel[j]]] += 1.0
        return (counts + counts.T) / 2.0

    observed = adjacency_counts(labels)
    rng = np.random.default_rng(0)
    perms = np.zeros((n_perm, n_types, n_types), dtype=float)
    for p in range(n_perm):
        shuffled = labels.copy()
        for img in np.unique(images):
            sel = np.where(images == img)[0]
            shuffled[sel] = rng.permutation(shuffled[sel])
        perms[p] = adjacency_counts(shuffled)
    mean = perms.mean(0)
    std = perms.std(0)
    std[std == 0] = 1.0
    zscores = (observed - mean) / std
    pvalues = (np.sum(perms >= observed[None], axis=0) + 1.0) / (n_perm + 1.0)

    return EnrichmentResult(
        scope=f"mock:{scope}",
        cell_types=cell_types,
        zscores=zscores.tolist(),
        pvalues=pvalues.tolist(),
    )


# --- niches -----------------------------------------------------------------------


def find_niches(cohort: str = COHORT, n_niches: int | None = None) -> list[Niche]:
    """Find recurring niches (composition + marker program + name) for the cohort."""
    adata = _load()
    labeled, backend = _labeled_adata(adata, n_niches)
    _BACKEND["find_niches"] = backend
    niche_col = labeled.obs["niche"].to_numpy()
    return [
        _characterize(labeled, niche_id, niche_col == niche_id)
        for niche_id in _niche_ids(labeled)
    ]


def characterize_niche(niche_id: int) -> Niche:
    """Composition, marker program, and name for one niche."""
    adata = _load()
    labeled, _ = _labeled_adata(adata, None)
    mask = labeled.obs["niche"].to_numpy() == niche_id
    if not mask.any():
        raise ValueError(f"niche_id {niche_id} not present")
    return _characterize(labeled, niche_id, mask)


# --- survival -> prognostic -------------------------------------------------------


def _prognostic_for_niche(adata: ad.AnnData, niche_id: int) -> Prognostic | None:
    """Prognostic for one niche via the engine, else a lifelines fallback, else None."""

    def engine() -> Prognostic:
        from ..engine.outcome import niche_outcome as fn

        return fn(adata, niche_id)

    return _run(
        f"niche_outcome[{niche_id}]",
        engine,
        lambda: _survival_fallback(adata, niche_id),
    )


def _survival_fallback(adata: ad.AnnData, niche_id: int) -> Prognostic | None:
    """Cox + KM on per-patient niche abundance (lifelines). Returns None if not viable."""
    obs = adata.obs
    needed = {"os_month", "os_event", "patient_id", "niche"}
    if not needed.issubset(set(obs.columns)):
        return None
    try:
        import pandas as pd
        from lifelines import CoxPHFitter, KaplanMeierFitter

        frame = pd.DataFrame(
            {
                "patient_id": obs["patient_id"].astype(str).to_numpy(),
                "in_niche": (obs["niche"].to_numpy() == niche_id).astype(float),
                "os_month": obs["os_month"].to_numpy(dtype=float),
                "os_event": obs["os_event"].to_numpy().astype(int),
            }
        )
        per_patient = frame.groupby("patient_id").agg(
            abundance=("in_niche", "mean"),
            os_month=("os_month", "first"),
            os_event=("os_event", "first"),
        )
        per_patient = per_patient[np.isfinite(per_patient["os_month"])]
        n_patients = int(per_patient.shape[0])
        if (
            n_patients < 2
            or per_patient["os_event"].sum() == 0
            or per_patient["abundance"].nunique() < 2
        ):
            return None

        # Small cohorts (the 3-patient mock) make an unregularized Cox fit blow up
        # into absurd hazard ratios. Regularize harder when n is tiny, then clamp to
        # a plausible finite band so the object is always sane. Real cohorts (~281
        # patients) hit penalizer 0.1 and no clamping.
        penalizer = 0.1 if n_patients >= 20 else 1.0
        cph = CoxPHFitter(penalizer=penalizer)
        cph.fit(per_patient, duration_col="os_month", event_col="os_event")
        row = cph.summary.loc["abundance"]
        hazard_ratio = float(np.exp(row["coef"]))
        ci_low = float(np.exp(row["coef lower 95%"]))
        ci_high = float(np.exp(row["coef upper 95%"]))
        pvalue = float(row["p"])
        hazard_ratio, ci_low, ci_high = _sane_hr(
            hazard_ratio, ci_low, ci_high, n_patients
        )

        km = _km_curve(per_patient)
        return Prognostic(
            hazard_ratio=hazard_ratio,
            ci_low=ci_low,
            ci_high=ci_high,
            pvalue=pvalue,
            n_patients=n_patients,
            km=km,
        )
    except Exception as exc:
        logger.warning("survival fallback for niche %s failed: %s", niche_id, exc)
        return None


def _sane_hr(
    hr: float, ci_low: float, ci_high: float, n_patients: int
) -> tuple[float, float, float]:
    """Clamp a hazard ratio + CI to a finite, plausible band (logs if it was degenerate).

    Direction (HR>1 vs HR<1) is preserved; only extreme/non-finite magnitudes from an
    underpowered fit are reined in. On a well-powered cohort nothing is clamped.
    """
    lo_band, hi_band = 1e-2, 1e2

    def clamp(v: float, lo: float, hi: float) -> float:
        if not np.isfinite(v):
            return hi if v > 0 else lo
        return float(min(max(v, lo), hi))

    raw = (hr, ci_low, ci_high)
    hr = clamp(hr, lo_band, hi_band)
    ci_low = clamp(ci_low, 1e-3, hr)
    ci_high = clamp(ci_high, hr, 1e3)
    if raw != (hr, ci_low, ci_high):
        logger.warning(
            "regularized underpowered survival estimate (n=%d): %s -> %s",
            n_patients,
            raw,
            (hr, ci_low, ci_high),
        )
    return hr, ci_low, ci_high


def _km_curve(per_patient) -> "object | None":
    """Kaplan-Meier survival for high vs low niche-abundance patient groups."""
    try:
        from lifelines import KaplanMeierFitter

        from ..schema import KMCurve

        med = float(per_patient["abundance"].median())
        high = per_patient[per_patient["abundance"] > med]
        low = per_patient[per_patient["abundance"] <= med]
        if high.shape[0] == 0 or low.shape[0] == 0:
            return None
        t_max = float(per_patient["os_month"].max())
        grid = [round(t_max * i / 12.0, 3) for i in range(13)]

        def surv(group) -> list[float]:
            kmf = KaplanMeierFitter().fit(group["os_month"], group["os_event"])
            return [float(v) for v in kmf.survival_function_at_times(grid).to_numpy()]

        return KMCurve(time=grid, high=surv(high), low=surv(low))
    except Exception:
        return None


def find_prognostic_niches(
    cohort: str = COHORT, patient_subset: list[str] | None = None
) -> list[Niche]:
    """Orchestrator: niches + survival association + naming, ranked worst-survival first.

    Chains find_niches -> characterize -> survival -> naming, then ranks by prognostic
    strength (highest hazard ratio first). Never crashes: if survival is unavailable,
    ranks by niche size with prognostic=None; naming always has a deterministic fallback.
    """
    adata = _load()
    if patient_subset:
        keep = (
            adata.obs["patient_id"]
            .astype(str)
            .isin([str(p) for p in patient_subset])
            .to_numpy()
        )
        if not keep.any():
            raise ValueError(f"no cells for patient_subset {patient_subset}")
        adata = adata[keep]

    labeled, backend = _labeled_adata(adata, None)
    _BACKEND["find_prognostic_niches"] = backend
    niche_col = labeled.obs["niche"].to_numpy()
    total = niche_col.shape[0]

    scored: list[tuple[float, float, Niche]] = []
    for niche_id in _niche_ids(labeled):
        mask = niche_col == niche_id
        niche = _characterize(labeled, niche_id, mask)
        niche.prognostic = _prognostic_for_niche(labeled, niche_id)
        abundance = float(mask.sum()) / max(total, 1)
        # rank key: prognostic niches first by descending hazard ratio, then by size
        has_prog = 1.0 if niche.prognostic is not None else 0.0
        hazard = (
            niche.prognostic.hazard_ratio if niche.prognostic is not None else abundance
        )
        scored.append((has_prog, hazard, niche))

    scored.sort(key=lambda t: (t[0], t[1]), reverse=True)
    return [niche for _, _, niche in scored]


# --- tissue map -------------------------------------------------------------------


def get_map_payload(image_id: str, color_mode: str = "cell_type") -> MapPayload:
    """Build the tissue-map render payload for one image (delegates to viz.payload)."""
    adata = _load()
    _BACKEND["get_map_payload"] = "real"
    return build_map_payload(adata, image_id, color_mode)


# --- precomputed findings (the demo path; served instantly from cache) -------------


@functools.lru_cache(maxsize=1)
def _findings() -> dict:
    data = json.loads(_FINDINGS.read_text())
    data["niches"] = {int(k): v for k, v in data["niches"].items()}
    return data


def correlate_niche_outcome(niche_id: int) -> dict:
    """One niche's abundance vs overall survival, with the full honesty bundle."""
    from ..engine.outcome import correlate_niche_outcome as _engine

    return _engine(_findings(), int(niche_id))


def describe_niches() -> list[dict]:
    """The niche catalog (name, major composition, size) from the precomputed cache."""
    return [
        {
            "niche_id": k,
            "name": v["name"],
            "composition": v["composition_major"],
            "n_cells": v["n_cells"],
            "n_cores": v["n_cores"],
        }
        for k, v in sorted(_findings()["niches"].items())
    ]


# --- risk layer (impact): patient-level risk, inseparable from its verdict ----------
# Thin wrappers over engine.risk. No analysis here. Hard guard: a RiskScore is
# physically incapable of existing without its RiskEvidence (the schema requires it),
# and if the model cannot be fit/evaluated these tools return a RiskScore whose
# evidence.verdict is "not evaluable", never a bare number.


@functools.lru_cache(maxsize=1)
def _risk_model() -> RiskModelCard:
    """Fit (and cache) the risk model on the loaded cohort.

    Reuses the precomputed honesty bundle (demo/findings.json) for the multiplicity /
    power fields when it matches the loaded cohort, so no slow permutation runs per
    request. On the tiny mock the model is simply not evaluable.
    """
    from ..engine import risk

    adata = _load()
    honesty = None
    try:
        findings = _findings()
        if int(findings["cohort"]["n_patients"]) == int(
            adata.obs["patient_id"].astype(str).nunique()
        ):
            honesty = findings
    except Exception:
        honesty = None
    return risk.fit_risk_model(adata, honesty=honesty)


def _resolve_patient(adata: ad.AnnData, patient_id, image_id) -> str:
    if patient_id is not None:
        # A nonexistent patient must ERROR, never be scored as a fabricated placeholder
        # (that would let a number the model never computed escape with a verdict).
        if str(patient_id) not in set(adata.obs["patient_id"].astype(str)):
            raise ValueError(f"patient_id {patient_id!r} not found")
        return str(patient_id)
    if image_id is not None:
        sub = adata.obs[adata.obs["image_id"].astype(str) == str(image_id)]
        if sub.shape[0] == 0:
            raise ValueError(f"image_id {image_id!r} not found")
        return str(sub["patient_id"].iloc[0])
    raise ValueError("predict_risk requires patient_id or image_id")


def _not_evaluable_score(
    model: RiskModelCard, patient_id: str | None, image_id: str | None
) -> RiskScore:
    """A RiskScore whose number is a placeholder; evidence.verdict says do not act."""
    return RiskScore(
        patient_id=patient_id,
        image_id=image_id,
        risk_score=0.0,
        risk_percentile=50.0,
        risk_group="intermediate",
        top_contributing_niches=[],
        evidence=model.evidence,
    )


def predict_risk(
    patient_id: str | None = None, image_id: str | None = None
) -> RiskScore:
    """Patient-level risk from niche composition, inseparable from its trust verdict."""
    from ..engine import risk

    adata = _load()
    model = _risk_model()
    pid = _resolve_patient(
        adata, patient_id, image_id
    )  # raises if the patient is unknown
    _BACKEND["predict_risk"] = "real" if model.coefficients else "fallback"
    if not model.coefficients:
        # The MODEL cannot evaluate anyone (e.g. the tiny mock): not-evaluable verdict.
        return _not_evaluable_score(model, pid, image_id)
    # The patient is known and the model is fitted; score_patient raises only if the
    # patient cannot be scored (e.g. no survival), which we surface as an error rather
    # than fabricating a number with the model's cohort verdict.
    score = risk.score_patient(adata, model, pid)
    if image_id is not None:
        score.image_id = str(image_id)
    return score


def rank_patients_by_risk(
    cohort: str = COHORT, top_n: int | None = None
) -> list[RiskScore]:
    """Cohort patients ranked by risk (highest first). Each carries the same evidence."""
    from ..engine import risk

    adata = _load()
    model = _risk_model()
    _BACKEND["rank_patients_by_risk"] = "real" if model.coefficients else "fallback"
    scores = risk.score_cohort(adata, model)
    scores.sort(key=lambda s: s.risk_score, reverse=True)
    return scores[: int(top_n)] if top_n else scores


def get_risk_model_card() -> RiskModelCard:
    """The fitted risk model's provenance + evidence (c-index, calibration, verdict)."""
    model = _risk_model()
    _BACKEND["get_risk_model_card"] = "real" if model.coefficients else "fallback"
    return model
