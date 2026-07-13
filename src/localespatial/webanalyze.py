"""Run the whole Locale engine on an arbitrary AnnData and return the viz bundle.

This is what turns the site into a live analysis engine: a user uploads a spatial
single-cell object, this runs niche detection, characterization, neighborhood
enrichment, niche->survival, and the validation battery, and returns exactly the
JSON shape the report frontend (window.LOCALE_LAB) renders.

Nothing is precomputed here; every number comes from the real engine on the given
object, with graceful fallbacks so odd inputs degrade instead of crashing.
"""

from __future__ import annotations

import warnings

import numpy as np
from anndata import AnnData

warnings.filterwarnings("ignore")


def _bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return []
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out.tolist()


def canonicalize(adata: AnnData) -> AnnData:
    """Shape an uploaded object to the schema the engine expects."""
    obs = adata.obs
    if "cell_type" not in obs:
        for alt in ("cell type", "celltype", "cell_types", "phenotype", "label"):
            if alt in obs:
                adata.obs["cell_type"] = obs[alt].astype("category")
                break
    if "cell_type" not in adata.obs:
        raise ValueError("need obs['cell_type'] (or 'cell type'/'phenotype'/'label')")
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    if "image_id" not in adata.obs:
        for alt in ("ImageId", "image", "core", "sample", "roi", "fov"):
            if alt in adata.obs:
                adata.obs["image_id"] = adata.obs[alt].astype(str)
                break
        else:
            adata.obs["image_id"] = "image_1"
    adata.obs["image_id"] = adata.obs["image_id"].astype(str).astype("category")
    if "spatial" not in adata.obsm:
        for a, b in (
            ("x", "y"),
            ("X", "Y"),
            ("x_centroid", "y_centroid"),
            ("Pos_X", "Pos_Y"),
        ):
            if a in adata.obs and b in adata.obs:
                adata.obsm["spatial"] = adata.obs[[a, b]].to_numpy(dtype=float)
                break
    if "spatial" not in adata.obsm:
        raise ValueError("need obsm['spatial'] (or x/y columns in obs)")
    adata.uns["markers"] = list(adata.var_names)
    return adata


def risk_bundle(card, scores) -> dict:
    """Flatten a RiskModelCard + per-patient RiskScores into the report's JSON shape."""
    names: dict[int, str] = {}
    for s in scores:
        for nc in s.top_contributing_niches:
            names[nc.niche_id] = nc.name
    ev = card.evidence
    coeffs = sorted(
        [
            {
                "niche_id": nid,
                "name": names.get(nid, f"niche {nid}"),
                "coef": float(cf.coef),
                "hr_per_sd": float(cf.hr_per_sd),
                "ci_95": [float(cf.ci_95[0]), float(cf.ci_95[1])],
                "p": float(cf.p),
            }
            for nid, cf in card.coefficients.items()
        ],
        key=lambda d: -abs(d["coef"]),
    )
    groups = {"low": 0, "intermediate": 0, "high": 0}
    for s in scores:
        groups[s.risk_group] = groups.get(s.risk_group, 0) + 1

    def pat(s) -> dict:
        return {
            "patient_id": s.patient_id,
            "image_id": s.image_id,
            "risk_score": float(s.risk_score),
            "risk_percentile": float(s.risk_percentile),
            "risk_group": s.risk_group,
            "top": [
                {"name": nc.name, "contribution": float(nc.contribution)}
                for nc in s.top_contributing_niches[:3]
            ],
        }

    ranked = sorted(scores, key=lambda s: -s.risk_score)
    return {
        "verdict": ev.verdict,
        "verdict_reason": ev.verdict_reason,
        "c_index_cv": float(card.c_index_cv),
        "c_index_ci_95": [float(x) for x in card.c_index_ci_95],
        "calibration_slope": float(card.calibration_slope),
        "n_train_patients": int(card.n_train_patients),
        "n_events": int(card.n_events),
        "cv_folds": int(card.cv_folds),
        "covariates_adjusted": list(card.covariates_adjusted),
        "coefficients": coeffs,
        "groups": groups,
        "n_scored": len(scores),
        "top_patients": [pat(s) for s in ranked[:5]],
        "bottom_patients": [pat(s) for s in ranked[-5:][::-1]],
    }


def analyze(adata: AnnData, n_niches: int = 6) -> dict:
    adata = canonicalize(adata)
    from .engine import characterize as CH
    from .engine import enrichment as EN
    from .engine import niches as NI
    from .engine import outcome as OU
    from .engine import validate as VAL
    from .mcp_server import interpret

    n_niches = int(max(2, min(n_niches, adata.obs["cell_type"].nunique() + 4)))
    labeled = NI.find_niches(adata, n_niches=n_niches)
    ids = sorted({int(x) for x in labeled.obs["niche"].to_numpy()})
    niche_col = labeled.obs["niche"].to_numpy().astype(int)

    cards = CH.characterize_all_niches(labeled)
    for c in cards:
        c.marker_program = [str(m) for m in c.marker_program]
        if not c.name:
            c.name = interpret.name_niche(c.composition, c.marker_program, c.niche_id)

    # per-niche survival (if outcome columns exist)
    prog = {}
    for nid in ids:
        try:
            prog[nid] = OU.niche_outcome(labeled, nid)
        except Exception:
            prog[nid] = None

    pvals = [prog[nid].pvalue if prog.get(nid) else 1.0 for nid in ids]
    qvals = _bh_fdr(pvals)

    findings_niches, survival = {}, []
    for k, nid in enumerate(ids):
        card = next((c for c in cards if c.niche_id == nid), None)
        p = prog.get(nid)
        comp = card.composition if card else {}
        findings_niches[str(nid)] = {
            "name": (card.name if card else f"niche {nid}"),
            "composition_major": comp,
            "top_metaclusters": (card.marker_program[:4] if card else []),
            "n_cells": int((niche_col == nid).sum()),
            "n_cores": int(labeled.obs["image_id"].nunique()),
            "hazard_ratio": float(p.hazard_ratio) if p else 1.0,
            "ci_95": [float(p.ci_low), float(p.ci_high)] if p else [1.0, 1.0],
            "p_raw": float(p.pvalue) if p else 1.0,
            "q_fdr": float(qvals[k]),
        }
        if p and p.km:
            survival.append(
                {
                    "name": card.name if card else f"niche {nid}",
                    "hazard_ratio": float(p.hazard_ratio),
                    "ci_low": float(p.ci_low),
                    "ci_high": float(p.ci_high),
                    "pvalue": float(p.pvalue),
                    "n_patients": int(p.n_patients),
                    "km": {
                        "time": list(p.km.time),
                        "high": list(p.km.high),
                        "low": list(p.km.low),
                    },
                }
            )

    try:
        enr = EN.compute_enrichment(adata, "cohort:uploaded")
        enrichment = {
            "cell_types": enr.cell_types,
            "zscores": enr.zscores,
            "pvalues": enr.pvalues,
        }
    except Exception as exc:  # noqa: BLE001
        enrichment = {"cell_types": [], "zscores": [], "pvalues": [], "error": str(exc)}

    try:
        shuffle = VAL.shuffle_negative_control(labeled, n_permutations=60)
        stability = float(VAL.stability_ari(labeled, n_runs=6))
        markers = []
        for nid in ids:
            m = VAL.marker_validation(labeled, nid)
            ev = m.get("evidence", {})
            markers.append(
                {
                    "niche_id": nid,
                    "passed": bool(m.get("passed")),
                    "corr": float(ev.get("marker_composition_corr"))
                    if isinstance(ev, dict)
                    and ev.get("marker_composition_corr") is not None
                    else None,
                    "top": [str(x) for x in ev.get("top_markers", [])]
                    if isinstance(ev, dict)
                    else [],
                }
            )
        validation = {
            "shuffle": {
                k: (float(v) if not isinstance(v, bool) else v)
                for k, v in shuffle.items()
            },
            "stability_ari": stability,
            "markers": markers,
        }
    except Exception as exc:  # noqa: BLE001
        validation = {"error": str(exc)}

    # patient-level risk model (Cox over niche abundances), inseparable from its verdict
    try:
        from .engine import risk as RISK

        cov = [c for c in ("grade", "clinical_type") if c in adata.obs]
        card = RISK.fit_risk_model(labeled, covariates=cov, n_perm=200)
        scores = RISK.score_cohort(
            labeled, card, niche_names={c.niche_id: c.name for c in cards}
        )
        risk = risk_bundle(card, scores)
    except Exception as exc:  # noqa: BLE001
        risk = {
            "verdict": "not evaluable",
            "verdict_reason": str(exc),
            "coefficients": [],
            "groups": {},
            "top_patients": [],
            "bottom_patients": [],
        }

    obs = adata.obs
    has_surv = "os_month" in obs
    cohort = {
        "n_patients": int(obs["patient_id"].nunique())
        if "patient_id" in obs
        else int(obs["image_id"].nunique()),
        "n_events": int(np.nansum(obs["os_event"].to_numpy(dtype=float)))
        if "os_event" in obs
        else 0,
        "n_hypotheses_tested": len(ids),
        "observed_best_p": float(min(pvals)) if pvals else 1.0,
        "p_selection_aware": float(min(1.0, min(pvals) * len(ids))) if pvals else 1.0,
        "min_detectable_hr": None,
        "power": None,
        "alpha": 0.05,
    }

    return {
        "source": "uploaded",
        "findings": {"cohort": cohort, "niches": findings_niches},
        "engine": {
            "cohort_label": "computed live on your uploaded data",
            "enrichment": enrichment,
            "survival": survival,
            "validation": validation,
        },
        "inputs": {
            "command": "POST /api/analyze",
            "data_file": "uploaded object",
            "schema": "src/localespatial/schema.py (Pydantic v2 contract)",
            "n_cells": int(adata.n_obs),
            "n_markers": int(adata.n_vars),
            "n_images": int(obs["image_id"].nunique()),
            "cell_types": sorted({str(x) for x in obs["cell_type"]}),
            "n_niches": len(ids),
            "has_survival": bool(has_surv),
        },
        "risk": risk,
        "figures": [],
        "k_sweep": None,
        "tests": [],
        "meta": {},
    }
