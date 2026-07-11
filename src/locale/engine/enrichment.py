"""Neighborhood enrichment / co-location (Lane A). Backed by squidpy.

Answers "who is next to whom": for each cell-type pair, are they spatial neighbors
more or less often than chance? Output is the cell_type x cell_type z-score matrix.
"""

from __future__ import annotations

from anndata import AnnData

from ..schema import EnrichmentResult


def compute_enrichment(adata: AnnData, scope: str) -> EnrichmentResult:
    """Neighborhood enrichment via squidpy.gr.nhood_enrichment (permutation test).

    Args:
        adata: canonical AnnData. The spatial graph must already be built
            (see engine.graph.build_spatial_graph); obs['cell_type'] is the grouping.
        scope: label describing what was analyzed, e.g. "cohort:breast" or
            "image:<image_id>". Copied into the returned EnrichmentResult.

    Returns:
        EnrichmentResult with cell_types (row/col order), the z-score matrix, and
        the permutation p-value matrix.

    TODO(Lane A):
        1. ensure the spatial graph exists (build it per image_id if not).
        2. squidpy.gr.nhood_enrichment(adata, cluster_key="cell_type", seed=0).
        3. read adata.uns["cell_type_nhood_enrichment"]; note squidpy stores
           "zscore" and "count" here, NOT a "pvalue" key. Derive the p-values
           yourself from the permutation null (e.g. pass n_perms and compute a
           two-sided empirical p per pair from the permuted count distribution).
        4. pack into EnrichmentResult(scope=scope, cell_types=<category order>,
           zscores=..., pvalues=...).
    """
    raise NotImplementedError(
        "compute_enrichment: wire squidpy.gr.nhood_enrichment -> EnrichmentResult."
    )
