"""Neighborhood enrichment / co-location (Lane A). Backed by squidpy.

Answers "who is next to whom": for each cell-type pair, are they spatial neighbors
more or less often than chance? Output is the cell_type x cell_type z-score matrix.
"""

from __future__ import annotations

import numpy as np
import squidpy as sq
from anndata import AnnData
from scipy.stats import norm

from ..schema import EnrichmentResult
from .graph import build_spatial_graph


def compute_enrichment(
    adata: AnnData, scope: str, seed: int = 0, n_perms: int = 1000
) -> EnrichmentResult:
    """Neighborhood enrichment via squidpy.gr.nhood_enrichment (permutation test).

    Args:
        adata: canonical AnnData. If the spatial graph is not present it is built
            per image_id first. obs['cell_type'] is the grouping.
        scope: label describing what was analyzed, e.g. "cohort:breast" or
            "image:<image_id>". If it starts with "image:" the analysis is
            restricted to that one image.
        seed: permutation seed (deterministic output).
        n_perms: permutations for the enrichment null.

    Returns:
        EnrichmentResult with cell_types (row/col order), the z-score matrix, and a
        two-sided p-value matrix.
    """
    if scope.startswith("image:"):
        target = build_spatial_graph(adata, image_id=scope.split(":", 1)[1])
    else:
        target = adata
        if "spatial_connectivities" not in target.obsp:
            build_spatial_graph(target)

    if str(target.obs["cell_type"].dtype) != "category":
        target.obs["cell_type"] = target.obs["cell_type"].astype("category")

    sq.gr.nhood_enrichment(
        target,
        cluster_key="cell_type",
        seed=seed,
        n_perms=n_perms,
        show_progress_bar=False,
    )
    result = target.uns["cell_type_nhood_enrichment"]
    cell_types = [str(c) for c in target.obs["cell_type"].cat.categories]
    zscores = np.asarray(result["zscore"], dtype=float)

    # squidpy stores "zscore" and "count", not a p-value. The z is already the
    # standardized permutation statistic, so a two-sided p follows from the normal.
    pvalues = 2.0 * norm.sf(np.abs(zscores))

    return EnrichmentResult(
        scope=scope,
        cell_types=cell_types,
        zscores=zscores.tolist(),
        pvalues=pvalues.tolist(),
    )
