"""Fit the patient-level risk model on the real Basel cohort and print the honest
evaluation: the cross-validated c-index with CI, the top and bottom risk patients,
an example patient's per-niche contributions, and the verdict. Saves the model card
to reports/risk_model_card.json.

    python scripts/run_basel_risk.py

Expects the real cohort at data/basel_niched.h5ad (obs['niche'] + survival, shared
out of band, never committed). If it is absent, this runs on a clearly-labeled
SYNTHETIC cohort matched to Basel's dimensions (281 patients, ~79 events, 12 niches
with no planted signal) so the pipeline and the honest verdict are demonstrable. The
synthetic run is NOT a real result and is labeled as such.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Allow `python scripts/run_basel_risk.py` from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.localespatial.engine import risk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASEL = ROOT / "data" / "basel_niched.h5ad"
FINDINGS = ROOT / "demo" / "findings.json"
REPORTS = ROOT / "reports"


def _synthetic_basel(seed: int = 0) -> ad.AnnData:
    """281 patients, 12 niches, ~79 events, NO planted niche->survival signal.

    This reproduces the SHAPE of Basel (so min_detectable_hr = 1.37 at 79 events) and
    the honest verdict (insufficient evidence), without pretending to be real data.
    """
    rng = np.random.default_rng(seed)
    n_patients, n_niches, cells = 281, 12, 30
    abundance = rng.dirichlet(np.ones(n_niches), size=n_patients)
    # survival independent of abundance (no signal)
    base = 0.02
    t_event = rng.exponential(1.0 / base, n_patients)
    t_cens = rng.exponential(1.0 / (base * 2.57), n_patients)
    os_month = np.clip(np.minimum(t_event, t_cens), 1, 200)
    os_event = (t_event <= t_cens).astype(int)
    grade = rng.integers(1, 4, n_patients)
    ctype = rng.choice(["LumA", "LumB", "HER2", "TNBC"], n_patients)

    rows = {
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
        for d in rng.choice(n_niches, size=cells, p=abundance[p]):
            rows["patient_id"].append(f"BASEL{p:03d}")
            rows["niche"].append(int(d))
            rows["os_month"].append(float(os_month[p]))
            rows["os_event"].append(int(os_event[p]))
            rows["grade"].append(int(grade[p]))
            rows["clinical_type"].append(str(ctype[p]))
            rows["cell_type"].append(f"class{d % 5}")
    obs = pd.DataFrame(rows)
    obs.index = obs.index.astype(str)
    return ad.AnnData(X=rng.normal(size=(len(obs), 4)).astype("float32"), obs=obs)


def _load() -> tuple[ad.AnnData, dict | None, bool]:
    if BASEL.exists():
        adata = ad.read_h5ad(BASEL)
        honesty = None
        if FINDINGS.exists():
            data = json.loads(FINDINGS.read_text())
            data["niches"] = {int(k): v for k, v in data["niches"].items()}
            n = adata.obs.get("PID", adata.obs.get("patient_id"))
            if n is not None and int(pd.Series(n).nunique()) == int(
                data["cohort"]["n_patients"]
            ):
                honesty = data
        return adata, honesty, True
    print("=" * 74)
    print("data/basel_niched.h5ad NOT FOUND.")
    print("Running on a SYNTHETIC cohort matched to Basel (281 patients, ~79 events,")
    print("12 niches, no planted signal). This is NOT a real result; drop in the real")
    print("object to get real numbers.")
    print("=" * 74)
    return _synthetic_basel(), None, False


def _print_scores(title: str, scores: list) -> None:
    print(f"\n{title}")
    for s in scores:
        drivers = ", ".join(
            f"{c.name} ({c.contribution:+.2f})" for c in s.top_contributing_niches[:2]
        )
        print(
            f"  {s.patient_id:>10}  risk={s.risk_score:+.3f}  "
            f"pct={s.risk_percentile:5.1f}  {s.risk_group:>12}  <- {drivers}"
        )


def main() -> None:
    adata, honesty, is_real = _load()
    n_perm = 1000 if honesty is not None else 300
    print(
        f"\nfitting risk model ({'REAL Basel' if is_real else 'SYNTHETIC'} cohort)..."
    )
    model = risk.fit_risk_model(adata, honesty=honesty, n_perm=n_perm)
    ev = model.evidence

    print("\n--- risk model card ---")
    print(f"  features (niches):   {model.features}")
    print(f"  covariates adjusted: {model.covariates_adjusted}")
    print(f"  train patients:      {model.n_train_patients}   events: {model.n_events}")
    print(f"  cv folds:            {model.cv_folds}")
    print(f"  c-index (out-of-fold): {ev.c_index_cv:.3f}  95% CI {ev.c_index_ci_95}")
    print(f"  calibration slope:   {ev.calibration_slope:.3f}")
    print(
        f"  n hypotheses tested: {ev.n_hypotheses_tested}   q_fdr_min: {ev.q_fdr_min}"
    )
    print(f"  p_selection_aware:   {ev.p_selection_aware}")
    print(f"  min detectable HR:   {ev.min_detectable_hr}")
    print(f"\n  VERDICT: {ev.verdict.upper()}")
    print(f"  reason:  {ev.verdict_reason}")

    if model.coefficients:
        scores = risk.score_cohort(adata, model)
        scores.sort(key=lambda s: s.risk_score, reverse=True)
        _print_scores("--- highest risk patients ---", scores[:3])
        _print_scores("--- lowest risk patients ---", scores[-3:])

        example = scores[0]
        print(f"\n--- per-niche contributions for {example.patient_id} (top risk) ---")
        for c in example.top_contributing_niches:
            print(
                f"  {c.name:>28}  abundance={c.abundance:.3f}  "
                f"coef={c.coefficient:+.3f}  contribution={c.contribution:+.3f}"
            )
        total = sum(c.contribution for c in example.top_contributing_niches)
        print(
            f"  {'sum of contributions':>28}  = {total:+.3f}  (== risk_score {example.risk_score:+.3f})"
        )

    REPORTS.mkdir(exist_ok=True)
    out = REPORTS / "risk_model_card.json"
    card = model.model_dump()
    card["_source"] = "real_basel" if is_real else "synthetic_demo"
    out.write_text(json.dumps(card, indent=2))
    print(f"\nsaved model card -> {out}")


if __name__ == "__main__":
    main()
