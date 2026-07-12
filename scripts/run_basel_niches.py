"""Enrichment sanity + niche discovery on the real Basel cohort (STEP 3-4).

Reads data/basel_nomarkers.h5ad (scripts/build_basel.py), prints the metacluster
neighborhood-enrichment sanity (immune-immune +, tumor-immune -), then discovers
window-only cellular neighborhoods over metacluster_id, sweeping k by stability ARI.

Window-only (NOT [identity | window]): the one-hot identity block has magnitude 1 in
a single dimension while the window is spread over 27, so in squared Euclidean
distance the identity swamps the window and k-means degenerates to a label
partition. And stability_ari cannot catch that, because a degenerate clustering is
extremely stable. Stability is not validity.

    python scripts/run_basel_niches.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score

from src.locale.engine.niches import find_niches, niche_features
from src.locale.metaclusters import METACLUSTERS

DATA = Path(__file__).resolve().parents[1] / "data"


def enrichment_sanity(a: ad.AnnData) -> None:
    a.obs["metacluster_id"] = a.obs["metacluster_id"].astype("category")
    sq.gr.nhood_enrichment(
        a, cluster_key="metacluster_id", seed=0, n_perms=100, show_progress_bar=False
    )
    z = np.asarray(a.uns["metacluster_id_nhood_enrichment"]["zscore"], dtype=float)
    ids = list(a.obs["metacluster_id"].cat.categories)
    major = np.array([METACLUSTERS[i][1] for i in ids])

    def block(r, c):
        return float(np.nanmean(z[np.ix_(major == r, major == c)]))

    print("=== STEP 3 enrichment sanity (metacluster block mean z) ===")
    print(
        f"  immune-immune {block('immune','immune'):.0f} | endo-endo {block('endothelial','endothelial'):.0f} "
        f"| stroma-stroma {block('stroma','stroma'):.0f} | tumor-tumor {block('tumor','tumor'):.0f}"
    )
    print(
        f"  tumor-immune {block('tumor','immune'):.0f}  (immune exclusion; expect strongly negative)"
    )


def sweep_k(a: ad.AnnData, ks=range(6, 13), seed: int = 0) -> int:
    X, _ = niche_features(a, cluster_key="metacluster_id", include_identity=False)
    cores = a.obs["core"].to_numpy()
    uniq = np.unique(cores)
    rng = np.random.default_rng(seed)
    print("\n=== STEP 4 k-sweep (stability ARI, window-only) ===")
    best_k, best_s = 6, -1.0
    for k in ks:
        runs = []
        for r in range(3):
            keep = set(rng.choice(uniq, size=int(len(uniq) * 0.8), replace=False))
            m = np.array([c in keep for c in cores])
            lab = MiniBatchKMeans(
                k, batch_size=4096, random_state=r, n_init=3
            ).fit_predict(X[m])
            full = np.full(a.n_obs, -1)
            full[m] = lab
            runs.append(full)
        aris = [
            adjusted_rand_score(
                runs[i][(runs[i] >= 0) & (runs[j] >= 0)],
                runs[j][(runs[i] >= 0) & (runs[j] >= 0)],
            )
            for i in range(3)
            for j in range(i + 1, 3)
        ]
        s = float(np.mean(aris))
        print(f"  k={k}: stability ARI {s:.3f}")
        if s > best_s:
            best_k, best_s = k, s
    print(f"  chosen k={best_k} (stability {best_s:.3f})")
    return best_k


def main() -> None:
    a = ad.read_h5ad(DATA / "basel_nomarkers.h5ad")
    enrichment_sanity(a)
    k = sweep_k(a)
    find_niches(a, n_niches=k, cluster_key="metacluster_id", include_identity=False)

    maj = a.obs["major"].astype(str).to_numpy()
    mc = a.obs["metacluster"].astype(str).to_numpy()
    niche = a.obs["niche"].to_numpy()
    cores = a.obs["core"].astype(str).to_numpy()
    print(f"\n=== {k} niches ===")
    for nid in sorted(set(niche)):
        m = niche == nid
        fr = {
            c: float((maj[m] == c).mean())
            for c in ["tumor", "immune", "stroma", "endothelial"]
        }
        top = pd.Series(mc[m]).value_counts(normalize=True).head(3)
        tag = (
            "  <- tumor-rich, immune-poor"
            if fr["tumor"] > 0.6 and fr["immune"] < 0.05
            else "  <- immune-RICH" if fr["immune"] > 0.3 else ""
        )
        print(
            f"niche {nid:>2} n={m.sum():>7,} cores={len(set(cores[m])):>3}  "
            f"T {fr['tumor']:.2f} I {fr['immune']:.2f} S {fr['stroma']:.2f}  "
            f"top {', '.join(top.index)}{tag}"
        )

    a.write_h5ad(DATA / "basel_niched.h5ad")
    print(f"\nwrote {DATA/'basel_niched.h5ad'}")


if __name__ == "__main__":
    main()
