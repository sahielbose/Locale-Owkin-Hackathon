"""Niche validation (Lane A). REQUIRED, not optional.

These three checks are how we prove the niches are real biology and not clustering
artifacts. A niche that fails them must not appear in the demo.
"""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from .graph import build_spatial_graph
from .niches import find_niches


def shuffle_negative_control(
    adata: AnnData, n_permutations: int = 100, seed: int = 0
) -> dict[str, float]:
    """Negative control: destroy spatial structure, confirm the signal collapses.

    Permutes cell-type labels WITHIN each image (breaking spatial arrangement while
    preserving each image's composition), re-clusters, and compares the real niche
    spatial coherence (fraction of a cell's graph neighbors sharing its niche) to
    the permuted null. If the niches survive shuffling they are an artifact.

    Returns:
        {"real_effect", "null_mean", "null_std", "empirical_p", "collapsed"}.
    """
    _ensure_niches(adata)
    graph = adata.obsp["spatial_connectivities"]
    niche_labels = adata.obs["niche"].to_numpy().astype(int)
    k = int(np.unique(niche_labels).size)
    real = _coherence(graph, niche_labels, k)

    cell_type = adata.obs["cell_type"].astype("category")
    type_codes = cell_type.cat.codes.to_numpy()
    n_types = cell_type.cat.categories.size
    images = adata.obs["image_id"].astype(str).to_numpy()
    image_idx = [np.where(images == img)[0] for img in np.unique(images)]

    rng = np.random.default_rng(seed)
    null = np.empty(n_permutations, dtype=float)
    for p in range(n_permutations):
        shuffled = type_codes.copy()
        for idx in image_idx:
            shuffled[idx] = rng.permutation(shuffled[idx])
        labels = _labels_from_type_codes(graph, shuffled, n_types, k, seed + p + 1)
        null[p] = _coherence(graph, labels, k)

    null_mean, null_std = float(null.mean()), float(null.std())
    return {
        "real_effect": real,
        "null_mean": null_mean,
        "null_std": null_std,
        "empirical_p": float((np.sum(null >= real) + 1) / (n_permutations + 1)),
        "collapsed": bool(real > null_mean + 3 * null_std and real > float(null.max())),
    }


def stability_ari(
    adata: AnnData, n_runs: int = 10, subsample_frac: float = 0.8, seed: int = 0
) -> float:
    """Stability: re-cluster on image subsamples, mean ARI vs the full-data labels.

    High mean ARI means the niches are reproducible, not seed/subsample noise.
    """
    _ensure_niches(adata)
    reference = adata.obs["niche"].to_numpy().astype(int)
    k = int(np.unique(reference).size)
    images = np.array(sorted(adata.obs["image_id"].astype(str).unique()))
    if len(images) < 2:
        return float("nan")

    rng = np.random.default_rng(seed)
    aris: list[float] = []
    for _ in range(n_runs):
        size = min(len(images), max(2, int(len(images) * subsample_frac)))
        keep = rng.choice(images, size=size, replace=False)
        mask = adata.obs["image_id"].astype(str).isin(keep).to_numpy()
        sub = adata[mask].copy()
        build_spatial_graph(sub)
        find_niches(sub, n_niches=k, seed=seed)
        aris.append(
            adjusted_rand_score(
                reference[mask], sub.obs["niche"].to_numpy().astype(int)
            )
        )
    return float(np.mean(aris))


def marker_validation(
    adata: AnnData, niche_id: int, threshold: float = 0.5
) -> dict[str, object]:
    """Biological sanity: does the niche's marker profile match its cell composition?

    Correlates the niche's observed mean marker profile with the profile PREDICTED
    from its cell-type composition (composition x per-cell-type mean markers). A high
    correlation means the marker program is explained by who is actually in the niche.

    Returns:
        {"niche_id", "passed": bool, "evidence": {marker_composition_corr, top_markers}}.
    """
    _ensure_niches(adata)
    niche_id = int(niche_id)
    x = adata.X
    x = np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)
    cell_type = adata.obs["cell_type"].astype(str)
    types = sorted(cell_type.unique())
    per_type = np.vstack(
        [
            (
                x[(cell_type == t).to_numpy()].mean(axis=0)
                if (cell_type == t).any()
                else np.zeros(x.shape[1])
            )
            for t in types
        ]
    )

    mask = adata.obs["niche"].to_numpy().astype(int) == niche_id
    if not mask.any():
        raise ValueError(f"niche_id {niche_id} not present")
    observed = x[mask].mean(axis=0)

    niche_types = cell_type[mask]
    composition = np.array(
        [float((niche_types == t).mean()) for t in types], dtype=float
    )
    expected = composition @ per_type

    if np.std(expected) > 0 and np.std(observed) > 0:
        corr = float(np.corrcoef(expected, observed)[0, 1])
    else:
        corr = float("nan")
    top = [str(adata.var_names[i]) for i in np.argsort(observed)[::-1][:6]]
    return {
        "niche_id": niche_id,
        "passed": bool(np.isfinite(corr) and corr > threshold),
        "evidence": {
            "marker_composition_corr": round(corr, 4) if np.isfinite(corr) else None,
            "top_markers": top,
        },
    }


# --- helpers ----------------------------------------------------------------------


def _ensure_niches(adata: AnnData) -> None:
    if "spatial_connectivities" not in adata.obsp:
        build_spatial_graph(adata)
    if "niche" not in adata.obs.columns:
        find_niches(adata)


def _one_hot(codes: np.ndarray, n: int) -> np.ndarray:
    out = np.zeros((len(codes), n), dtype=float)
    valid = codes >= 0
    out[np.arange(len(codes))[valid], codes[valid]] = 1.0
    return out


def _labels_from_type_codes(
    graph, type_codes: np.ndarray, n_types: int, k: int, seed: int
) -> np.ndarray:
    one_hot = _one_hot(type_codes, n_types)
    aggregated = graph @ one_hot + one_hot
    degree = np.asarray(graph.sum(axis=1)).ravel() + 1.0
    degree[degree == 0] = 1.0
    features = aggregated / degree[:, None]
    return KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(features)


def _coherence(graph, niche_labels: np.ndarray, k: int) -> float:
    one_hot = _one_hot(niche_labels, k)
    neighbor_same = np.asarray(graph @ one_hot)[
        np.arange(len(niche_labels)), niche_labels
    ]
    degree = np.asarray(graph.sum(axis=1)).ravel()
    degree = np.where(degree == 0, 1.0, degree)
    return float((neighbor_same / degree).mean())
