"""Ground-truth comparison against Jackson et al. 2020 (Lane A).

Uniquely, the archive ships the paper's own answers, so we can check our niches
against them without ever fitting to them:

  * community_ari      : Adjusted Rand Index of our per-cell niche vs the paper's
                         published tumor community phenotype.
  * enrichment_vs_published : correlation of our neighborhood-enrichment matrix vs
                         the paper's neighborhood heatmap. Kept for the record only;
                         it is NOT a validation (the two statistics differ, see note).

adjusted_rand_score here is EXTERNAL truth, a different thing from the SELF-stability
ARI in validate.stability_ari.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.metrics import adjusted_rand_score


def community_ari(
    adata: AnnData,
    communtiy_data_tumor_csv: str,
    pg_tumor_communities_csv: str,
    niche_key: str = "niche",
) -> dict:
    """ARI of our niche vs the paper's tumor community phenotype (tumor-restricted).

    Their coarse community phenotype (`cluster`, ~23 recurring types) is reached by
    joining the per-cell fine community in ``Communtiy_data_tumor.csv`` (note the
    archive's typo) to the ``(core, Community) -> cluster`` map in
    ``PG_tumor_communities.csv``. The per-cell fine ``Community`` (~1940) is NOT the
    comparison label; the coarse phenotype is.
    """
    ct = pd.read_csv(communtiy_data_tumor_csv)
    pg = pd.read_csv(pg_tumor_communities_csv)
    ct = ct.merge(
        pg[["core", "Community", "cluster"]], on=["core", "Community"], how="left"
    )
    ct["id"] = ct.core.astype(str) + "_" + ct.CellId.astype(str)
    theirs = (
        ct.dropna(subset=["cluster"])
        .drop_duplicates("id")
        .set_index("id")["cluster"]
        .astype(int)
    )

    has_major = "major" in adata.obs.columns
    ours = pd.DataFrame(
        {
            "niche": adata.obs[niche_key].to_numpy(),
            "major": (
                adata.obs["major"].astype(str).to_numpy() if has_major else "tumor"
            ),
        },
        index=adata.obs_names,
    )
    joined = ours.join(theirs, how="inner").dropna(subset=["cluster"])
    tumor = joined[joined["major"] == "tumor"] if has_major else joined
    return {
        "ari_tumor_only": float(adjusted_rand_score(tumor["niche"], tumor["cluster"])),
        "ari_all_matched": float(
            adjusted_rand_score(joined["niche"], joined["cluster"])
        ),
        "n_tumor": int(len(tumor)),
        "n_matched": int(len(joined)),
        "n_their_communities": int(theirs.nunique()),
        "n_our_niches": int(adata.obs[niche_key].nunique()),
    }


def enrichment_vs_published(
    adata: AnnData, neighborhood_heatmap_csv: str, cluster_key: str = "metacluster_id"
) -> dict:
    """Correlate our ``nhood_enrichment`` z-matrix vs the paper's neighborhood heatmap.

    NOTE: kept for the record, NOT a validation. The diagonal (self-enrichment) and
    the immune block align, but the off-diagonal does not correlate: our permutation
    z shows tumor subtypes SEGREGATING while their Delta statistic POOLS them into one
    tumor community. Different statistics measuring different things.
    """
    from scipy.stats import pearsonr

    key = f"{cluster_key}_nhood_enrichment"
    if key not in adata.uns:
        raise ValueError(
            f"run squidpy.gr.nhood_enrichment(cluster_key={cluster_key!r}) first"
        )
    ourz = np.asarray(adata.uns[key]["zscore"], dtype=float)
    theirs = pd.read_csv(neighborhood_heatmap_csv).to_numpy().astype(float)
    if ourz.shape != theirs.shape:
        return {"aligned": False, "our_shape": ourz.shape, "their_shape": theirs.shape}
    off = ~np.eye(ourz.shape[0], dtype=bool)
    return {
        "aligned": True,
        "pearson_full": float(pearsonr(ourz.ravel(), theirs.ravel())[0]),
        "pearson_offdiag": float(pearsonr(ourz[off], theirs[off])[0]),
        "note": "diagonal aligns; off-diagonal differs (different statistic). Not a validation.",
    }
