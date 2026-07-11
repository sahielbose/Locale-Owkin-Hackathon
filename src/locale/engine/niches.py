"""Niche discovery (Lane A). Backed by CellCharter (VAE + GMM), kmeans fallback.

A "niche" is a recurring cellular neighborhood that appears across many images.
This writes an integer niche label per cell into adata.obs['niche'].
"""

from __future__ import annotations

from anndata import AnnData


def find_niches(adata: AnnData, n_niches: int = 6) -> AnnData:
    """Detect recurring cellular niches across the cohort.

    Primary method: CellCharter, which aggregates each cell's l-hop neighborhood
    features and clusters them with a GMM, explicitly designed to compare niches
    across many samples. Fallback: k-means on per-cell local cell-type composition
    (the Schurch/Nolan "cellular neighborhoods" approach).

    Args:
        adata: canonical AnnData with the spatial graph already built
            (engine.graph.build_spatial_graph).
        n_niches: number of niches (GMM components / k-means k).

    Returns:
        The AnnData with adata.obs['niche'] (int, categorical-friendly) assigned.

    TODO(Lane A):
        1. build neighborhood features: cellcharter.gr.aggregate_neighbors
           (or windowed local cell-type-composition histograms per cell).
        2. cluster: cellcharter.tl.Cluster (GMM, n_clusters=n_niches) or KMeans.
        3. write adata.obs['niche'] = labels (int).
        4. return adata.
    """
    raise NotImplementedError(
        "find_niches: wire CellCharter aggregate_neighbors + GMM (kmeans fallback)."
    )
