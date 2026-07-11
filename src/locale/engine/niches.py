"""Niche discovery (Lane A). k-means on local cell-type composition.

A "niche" is a recurring cellular neighborhood that appears across many images.
Each cell is described by the cell-type composition of its spatial neighborhood
(the Schurch/Nolan "cellular neighborhoods" approach); clustering those vectors
gives niches that are defined by geography, not by a cell's own label, so the same
cell type can sit in several niches. This writes an integer niche label per cell
into adata.obs['niche'].

CellCharter (VAE + GMM, cross-sample) can be swapped in behind this same signature
later; k-means on the neighborhood composition is the reliable first pass.
"""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from sklearn.cluster import KMeans

from .graph import build_spatial_graph


def neighborhood_composition(
    adata: AnnData, cluster_key: str = "cell_type", include_self: bool = True
) -> tuple[np.ndarray, list[str]]:
    """Per-cell mean cell-type composition of its spatial neighborhood (self incl.)."""
    if "spatial_connectivities" not in adata.obsp:
        build_spatial_graph(adata)
    graph = adata.obsp["spatial_connectivities"]
    cats = list(adata.obs[cluster_key].astype("category").cat.categories)
    codes = adata.obs[cluster_key].astype("category").cat.codes.to_numpy()

    one_hot = np.zeros((adata.n_obs, len(cats)), dtype=float)
    valid = codes >= 0
    one_hot[np.arange(adata.n_obs)[valid], codes[valid]] = 1.0

    aggregated = graph @ one_hot
    degree = np.asarray(graph.sum(axis=1)).ravel()
    if include_self:
        aggregated = aggregated + one_hot
        degree = degree + 1.0
    degree[degree == 0] = 1.0
    return aggregated / degree[:, None], cats


def find_niches(adata: AnnData, n_niches: int = 6, seed: int = 0) -> AnnData:
    """Detect recurring cellular niches across the cohort.

    Method: k-means on each cell's local cell-type composition (the neighborhood
    window), computed over the per-image spatial graph so no window crosses a core.

    Args:
        adata: canonical AnnData; the per-image spatial graph is built if absent.
        n_niches: number of niches (k-means k).
        seed: random_state for reproducible labels.

    Returns:
        The AnnData with adata.obs['niche'] (int) assigned. If a view is passed
        (e.g. a patient subset) a copy is made and returned so squidpy and the obs
        write operate on a real object.
    """
    if adata.is_view:
        adata = adata.copy()
    features, _ = neighborhood_composition(adata)
    labels = KMeans(n_clusters=n_niches, n_init=10, random_state=seed).fit_predict(
        features
    )
    adata.obs["niche"] = labels.astype(np.int64)
    adata.uns["niche_params"] = {"n_niches": int(n_niches), "seed": int(seed)}
    return adata
