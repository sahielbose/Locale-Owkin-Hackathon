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


def _one_hot(adata: AnnData, cluster_key: str) -> tuple[np.ndarray, list[str]]:
    cats = list(adata.obs[cluster_key].astype("category").cat.categories)
    codes = adata.obs[cluster_key].astype("category").cat.codes.to_numpy()
    one_hot = np.zeros((adata.n_obs, len(cats)), dtype=float)
    valid = codes >= 0
    one_hot[np.arange(adata.n_obs)[valid], codes[valid]] = 1.0
    return one_hot, cats


def neighborhood_composition(
    adata: AnnData, cluster_key: str = "cell_type", include_self: bool = True
) -> tuple[np.ndarray, list[str]]:
    """Per-cell mean cell-type composition of its spatial neighborhood (self incl.)."""
    if "spatial_connectivities" not in adata.obsp:
        build_spatial_graph(adata)
    graph = adata.obsp["spatial_connectivities"]
    one_hot, cats = _one_hot(adata, cluster_key)
    aggregated = graph @ one_hot
    degree = np.asarray(graph.sum(axis=1)).ravel()
    if include_self:
        aggregated = aggregated + one_hot
        degree = degree + 1.0
    degree[degree == 0] = 1.0
    return aggregated / degree[:, None], cats


def niche_features(
    adata: AnnData, cluster_key: str = "cell_type", include_identity: bool = False
) -> tuple[np.ndarray, list[str]]:
    """Clustering features: the neighborhood window, optionally prepended with the
    cell's own one-hot identity (the ``[C | window]`` variant => 2*n_types dims)."""
    window, cats = neighborhood_composition(adata, cluster_key=cluster_key)
    if include_identity:
        ident, _ = _one_hot(adata, cluster_key)
        return np.hstack([ident, window]), cats
    return window, cats


def find_niches(
    adata: AnnData,
    n_niches: int = 6,
    seed: int = 0,
    cluster_key: str = "cell_type",
    include_identity: bool = False,
) -> AnnData:
    """Detect recurring cellular niches across the cohort.

    Method: k-means on each cell's local cell-type composition (the neighborhood
    window), computed over the per-image spatial graph so no window crosses a core.
    Set ``include_identity`` to cluster on ``[own one-hot | window]`` (2*n_types dims).
    ``cluster_key`` selects the label basis (e.g. ``metacluster_id``).

    Returns the AnnData with ``obs['niche']`` (int). A view is copied first so the
    squidpy call and obs write operate on a real object.
    """
    if adata.is_view:
        adata = adata.copy()
    features, _ = niche_features(adata, cluster_key=cluster_key, include_identity=include_identity)
    labels = KMeans(n_clusters=n_niches, n_init=10, random_state=seed).fit_predict(
        features
    )
    adata.obs["niche"] = labels.astype(np.int64)
    adata.uns["niche_params"] = {
        "n_niches": int(n_niches), "seed": int(seed),
        "cluster_key": cluster_key, "include_identity": bool(include_identity),
    }
    return adata
