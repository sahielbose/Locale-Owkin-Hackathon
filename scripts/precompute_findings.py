"""Precompute the demo path -> demo/findings.json (read-only over data/).

A four-minute niche discovery or a 1000x permutation running live on stage is how
we lose. So we compute everything ONCE from the frozen data/basel_niched.h5ad
(niche labels + the cached spatial graph are already in it, we do NOT rebuild) and
write a cache the MCP tools serve instantly:

  * the multiplicity + power aware survival bundle for every niche (engine.cohort_survival)
  * each niche's honest name, composition, and size

Writes ONLY demo/findings.json. Never touches data/.

    python scripts/precompute_findings.py            # 1000 permutations (demo)
    python scripts/precompute_findings.py --n-perm 200
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq

from src.localespatial.engine.outcome import cohort_survival
from src.localespatial.metaclusters import METACLUSTERS

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parents[1] / "demo" / "findings.json"

MAJOR_ORDER = ["immune", "endothelial", "stroma", "tumor"]


def enrichment_block(a: ad.AnnData) -> dict:
    """Neighborhood enrichment over the 25 unique cell-type NAMES, computed two ways.

    Granularity: the object caches an enrichment over 27 metacluster ids, which splits
    immune into two 'T cells' and two 'Macrophages' populations; that matrix's
    tumor-immune block averages to -22. Merging the duplicates into 25 real cell types
    (the granularity the report and pitch quote) gives -32. We use the 25 names.

    Permutation scope: squidpy's default permutes cell-type labels GLOBALLY across all
    cores, which conflates within-core spatial organization with between-core
    compositional heterogeneity (a 97%-tumor core, globally shuffled, is compared to a
    cohort-wide label mix). We therefore compute BOTH: the global z (the headline the
    report quotes) and the within-core z (library_key='core', the correct null that
    holds each core's composition fixed). Within core the tumor-immune block is smaller
    (about -10) but still clearly negative: the exclusion survives the correct null.

    Computed on the cached graph (no graph rebuild), deterministic seed.
    """
    b = a.copy()
    b.obs["mc_name"] = pd.Categorical(
        [METACLUSTERS[int(c)][0] for c in b.obs["metacluster_id"]]
    )
    b.obs["core"] = b.obs["core"].astype("category")
    name_to_major = {nm: mj for (nm, mj) in METACLUSTERS.values()}
    names = list(b.obs["mc_name"].cat.categories)
    major = [name_to_major[c] for c in names]
    order = sorted(
        range(len(names)), key=lambda i: (MAJOR_ORDER.index(major[i]), names[i])
    )
    names_o = [names[i] for i in order]
    major_o = [major[i] for i in order]
    m = np.array(major_o)

    def _blocks() -> tuple[list[list[int]], np.ndarray]:
        z = np.asarray(b.uns["mc_name_nhood_enrichment"]["zscore"], dtype=float)
        z_o = z[np.ix_(order, order)]
        blocks = [
            [
                int(round(float(np.nanmean(z_o[np.ix_(m == r, m == c)]))))
                for c in MAJOR_ORDER
            ]
            for r in MAJOR_ORDER
        ]
        return blocks, z_o

    sq.gr.nhood_enrichment(
        b, cluster_key="mc_name", seed=0, n_perms=1000, show_progress_bar=False
    )
    global_blocks, z_global = _blocks()

    sq.gr.nhood_enrichment(
        b,
        cluster_key="mc_name",
        library_key="core",
        seed=0,
        n_perms=1000,
        show_progress_bar=False,
    )
    within_blocks, _ = _blocks()

    ti, im = MAJOR_ORDER.index("tumor"), MAJOR_ORDER.index("immune")
    return {
        "granularity": "25 unique cell-type names (duplicate T-cell and Macrophage metaclusters merged)",
        "permutation": "global (squidpy default); within-core reported alongside",
        "names": names_o,
        "major": major_o,
        "zscores": [[round(float(v), 1) for v in row] for row in z_global],
        "major_blocks": {"cell_types": list(MAJOR_ORDER), "zscores": global_blocks},
        "within_core_major_blocks": {
            "cell_types": list(MAJOR_ORDER),
            "zscores": within_blocks,
        },
        "tumor_immune_global": global_blocks[ti][im],
        "tumor_immune_within_core": within_blocks[ti][im],
    }


def _name(fr: dict[str, float], top: list[str]) -> str:
    # honest, cell-level labels. "tumor-rich, immune-poor" NOT "immune-excluded":
    # a cell-level composition cannot tell exclusion (immune present, locked out) from
    # desert (no immune), and on Basel that phenotype is a continuum, not a clean split.
    if fr["immune"] > 0.3:
        return "immune-rich"
    if fr["tumor"] > 0.6 and fr["immune"] < 0.05:
        return "tumor-rich, immune-poor"
    if fr["tumor"] > 0.35 and fr["immune"] > 0.12:
        return "tumor-immune boundary"
    if fr["stroma"] + fr["endothelial"] > 0.45:
        return "stromal / vascular"
    return f"{top[0]}-mixed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=1000)
    args = ap.parse_args()

    a = ad.read_h5ad(DATA / "basel_niched.h5ad")  # niches + graph already frozen inside
    summary = cohort_survival(a, n_perm=args.n_perm)  # the slow part, done once

    maj = a.obs["major"].astype(str).to_numpy()
    mc = a.obs["metacluster"].astype(str).to_numpy()
    niche = a.obs["niche"].to_numpy()
    cores = a.obs["core"].astype(str).to_numpy()

    niches = {}
    for nid in sorted(set(int(n) for n in niche)):
        m = niche == nid
        fr = {
            c: round(float((maj[m] == c).mean()), 3)
            for c in ["tumor", "immune", "stroma", "endothelial"]
        }
        top = pd.Series(mc[m]).value_counts(normalize=True).head(3)
        entry = {
            "name": _name(fr, list(top.index)),
            "composition_major": fr,
            "top_metaclusters": [f"{i} ({v:.0%})" for i, v in top.items()],
            "n_cells": int(m.sum()),
            "n_cores": int(len(set(cores[m]))),
        }
        entry.update(summary["niches"][nid])  # hazard_ratio, ci_95, p_raw, q_fdr
        niches[nid] = entry

    # niche 7 phenotype caveat (from the excluded-vs-desert check) travels with the finding
    if 7 in niches:
        niches[7]["phenotype_note"] = (
            "cell-level 'tumor-rich, immune-poor'. Excluded-vs-desert is a continuum on "
            "Basel (deserts are real but there is no clean boundary), so this is NOT an "
            "'immune-excluded' claim. Testing that would need pre-registration on Zurich."
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        json.dumps(
            {
                "cohort": summary["cohort"],
                "niches": niches,
                "enrichment": enrichment_block(a),
            },
            indent=2,
        )
    )
    c = summary["cohort"]
    print(f"wrote {OUT}")
    print(
        f"cohort: {c['n_patients']} patients, {c['n_events']} events | "
        f"min detectable HR (80% power) = {c['min_detectable_hr']:.2f} | "
        f"selection-aware p = {c['p_selection_aware']:.3f}"
    )


if __name__ == "__main__":
    main()
