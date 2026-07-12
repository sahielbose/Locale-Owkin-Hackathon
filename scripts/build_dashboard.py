"""Build the dashboard data bundle from the FROZEN real-cohort outputs.

Reads data/basel_niched.h5ad (frozen niche labels + the cached neighborhood
enrichment) and demo/findings.json (the precomputed per-niche survival plus the
honesty bundle). It does NOT re-run the engine and does NOT re-cluster: the niche
numbering is frozen and is what the report, the transcript, and the PDF refer to.
Re-clustering would renumber the niches and break all three.

Writes src/localespatial/viz/app/dashboard_data.js for the offline dashboard.

    python scripts/build_dashboard.py
    LOCALE_DATA=/abs/path/other.h5ad python scripts/build_dashboard.py
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
APP_DIR = ROOT / "src" / "localespatial" / "viz" / "app"
FINDINGS = ROOT / "demo" / "findings.json"
RISK_CARD = ROOT / "reports" / "risk_model_card.json"

MAP_CORE = "BaselTMA_SP43_144_X15Y1"  # 32% tumor, compartmentalised: shows the niches
MAP_MAX_PER_CORE = 10_000  # never ship a whole core's worth of points to a canvas


def _data_path() -> Path:
    """The real cohort object, or an explicit override. Never the mock.

    Silently falling back to data/mock.h5ad is exactly how the dashboard ended up
    reporting 3 patients, so a missing real file is a hard error here.
    """
    env = os.environ.get("LOCALE_DATA")
    path = Path(env) if env else ROOT / "data" / "basel_niched.h5ad"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. The dashboard is built on the real cohort, not the "
            "mock. Build data/basel_niched.h5ad (scripts/build_basel.py then "
            "scripts/run_basel_niches.py) or set LOCALE_DATA to the real object. "
            "This script never falls back to data/mock.h5ad."
        )
    return path


def _risk() -> dict | None:
    """The risk model card's honest evaluation for the dashboard risk panel.

    Reads reports/risk_model_card.json (scripts/run_basel_risk.py). The panel exists to
    show the in-sample c-index next to the out-of-fold one: the number we declined to
    report, beside the one we did. Returns None if the card has not been generated.
    """
    if not RISK_CARD.exists():
        return None
    card = json.loads(RISK_CARD.read_text())
    ev = card["evidence"]
    return {
        "c_index_in_sample": card["c_index_in_sample"],
        "c_index_out_of_fold": card["c_index_out_of_fold"],
        "c_index_ci_95": card["c_index_ci_95"],
        "optimism_gap": card["optimism_gap"],
        "n_patients": card["n_train_patients"],
        "n_events": card["n_events"],
        "n_features": len(card["features"]),
        "verdict": ev["verdict"],
        "verdict_reason": ev["verdict_reason"],
    }


def _major_enrichment(findings: dict) -> dict:
    """The frozen 4-class neighborhood enrichment, read from demo/findings.json.

    findings['enrichment'] is computed over the 25 unique cell-type names, so the
    tumor vs immune headline reads -32, matching the report, the figures, and the
    PDF. Aggregating the object's cached 27-id matrix here would give -22 and disagree
    with the report, so this is read, not recomputed.
    """
    return findings["enrichment"]["major_blocks"]


def _cards(findings: dict) -> list[dict]:
    """One card per frozen niche, each carrying the full honesty bundle."""
    cohort = findings["cohort"]
    alpha = float(cohort.get("alpha", 0.05))
    cards = []
    for nid_str, n in findings["niches"].items():
        q = float(n["q_fdr"])
        stats = {
            "hazard_ratio": round(float(n["hazard_ratio"]), 3),
            "ci_95": [round(float(n["ci_95"][0]), 3), round(float(n["ci_95"][1]), 3)],
            "p_raw": round(float(n["p_raw"]), 4),
            "q_fdr": round(q, 4),
            "n_hypotheses_tested": int(cohort["n_hypotheses_tested"]),
            "p_selection_aware": round(float(cohort["p_selection_aware"]), 3),
            "n_events": int(cohort["n_events"]),
            "min_detectable_hr": round(float(cohort["min_detectable_hr"]), 3),
            # Same gate as engine.correlate_niche_outcome: q AND the global
            # selection-aware p must both clear alpha, or the verdict is insufficient.
            "verdict": (
                "supported"
                if q < alpha and float(cohort.get("p_selection_aware", 1.0)) < alpha
                else "insufficient evidence"
            ),
        }
        cards.append(
            {
                "niche_id": int(nid_str),
                "name": n["name"],
                "composition": n["composition_major"],
                "dominant": n.get("top_metaclusters", []),
                "n_cells": int(n["n_cells"]),
                "n_cores": int(n["n_cores"]),
                "note": n.get("phenotype_note", ""),
                "stats": stats,
            }
        )
    # Lead with the most significant-looking niche, so the honest verdict on it is
    # the first thing a reader sees.
    cards.sort(key=lambda c: c["stats"]["p_raw"])
    return cards


def _pick_core(a: ad.AnnData) -> str:
    cores = a.obs["core"].astype(str)
    if MAP_CORE in set(cores):
        return MAP_CORE
    core_arr = cores.to_numpy()
    maj = a.obs["major"].astype(str).to_numpy()
    counts = cores.value_counts()
    big = counts[counts > 1500].index
    best, best_d = str(counts.index[0]), 9.0
    for c in big:
        m = core_arr == c
        d = abs(float((maj[m] == "tumor").mean()) - 0.4)
        if d < best_d:
            best, best_d = str(c), d
    return best


def _map(a: ad.AnnData, core: str) -> dict:
    """One core's cells for the tissue map, downsampled to MAP_MAX_PER_CORE points."""
    core_arr = a.obs["core"].astype(str).to_numpy()
    idx = np.where(core_arr == core)[0]
    if idx.size > MAP_MAX_PER_CORE:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(idx, MAP_MAX_PER_CORE, replace=False))
    xy = np.asarray(a.obsm["spatial"], dtype=float)[idx]
    maj = a.obs["major"].astype(str).to_numpy()[idx]
    niche = a.obs["niche"].to_numpy().astype(int)[idx]
    units = [
        {
            "x": round(float(xy[i, 0]), 1),
            "y": round(float(xy[i, 1]), 1),
            "cell_type": str(maj[i]),  # major class, so the map colours cleanly
            "niche_id": int(niche[i]),
        }
        for i in range(idx.size)
    ]
    return {"image_id": core, "units": units, "color_mode": "niche"}


def main() -> None:
    a = ad.read_h5ad(_data_path())
    findings = json.loads(FINDINGS.read_text())
    cohort = findings["cohort"]
    cards = _cards(findings)
    core = _pick_core(a)

    bundle = {
        "cohort": "breast",
        "n_cells": int(a.n_obs),
        "n_patients": (
            int(a.obs["PID"].nunique())
            if "PID" in a.obs.columns
            else int(cohort["n_patients"])
        ),
        "n_images": int(a.obs["core"].nunique()),
        "n_niches": len(cards),
        "n_events": int(cohort["n_events"]),
        "context": {
            "n_hypotheses_tested": int(cohort["n_hypotheses_tested"]),
            "p_selection_aware": round(float(cohort["p_selection_aware"]), 3),
            "n_events": int(cohort["n_events"]),
            "min_detectable_hr": round(float(cohort["min_detectable_hr"]), 3),
        },
        "enrichment": _major_enrichment(findings),
        "niches": cards,
        "map": _map(a, core),
        "risk": _risk(),
    }

    # Guard: never let a mock-sized bundle out the door.
    assert (
        bundle["n_cells"] > 700_000
        and bundle["n_patients"] == 281
        and bundle["n_niches"] == 12
    ), (
        f"refusing to write a dashboard built on the wrong data: "
        f"{bundle['n_cells']} cells, {bundle['n_patients']} patients, "
        f"{bundle['n_niches']} niches"
    )

    APP_DIR.mkdir(parents=True, exist_ok=True)
    out = APP_DIR / "dashboard_data.js"
    out.write_text(
        "// Auto-generated by scripts/build_dashboard.py. Do not edit by hand.\n"
        "window.LOCALE_DASHBOARD = " + json.dumps(bundle) + ";\n"
    )
    print(f"wrote {out}")
    print(
        f"  {bundle['n_cells']:,} cells, {bundle['n_patients']} patients, "
        f"{bundle['n_niches']} niches, {bundle['n_images']} cores, "
        f"{bundle['n_events']} deaths"
    )
    print(
        f"  map core {bundle['map']['image_id']}: {len(bundle['map']['units'])} points"
    )
    print(f"  niche verdicts: {sorted({c['stats']['verdict'] for c in cards})}")


if __name__ == "__main__":
    main()
