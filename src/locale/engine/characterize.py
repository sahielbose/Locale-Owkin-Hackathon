"""Niche characterization (Lane A). Backed by scanpy.

For each niche: its cell-type composition and its top enriched markers (the
"marker program"). It is 35 protein markers, not genes, so this is a marker
program, not a gene program. Returns Niche objects with name left blank for
mcp_server.interpret to fill.
"""

from __future__ import annotations

import numpy as np
import scanpy as sc
from anndata import AnnData

from ..schema import Niche

_TOP_MARKERS = 6
_RANK_KEY = "locale_rank_genes"


def characterize_niche(
    adata: AnnData, niche_id: int, n_markers: int = _TOP_MARKERS
) -> Niche:
    """Composition + marker program for one niche.

    Args:
        adata: canonical AnnData with adata.obs['niche'] assigned
            (engine.niches.find_niches).
        niche_id: the niche to characterize.
        n_markers: how many top markers to report.

    Returns:
        Niche(niche_id, name="", composition=<cell_type->fraction>,
        marker_program=<top enriched markers>, prognostic=None). The name is filled
        later by mcp_server.interpret.
    """
    if "niche" not in adata.obs.columns:
        raise ValueError(
            "adata.obs['niche'] missing; run engine.niches.find_niches first"
        )
    niche_id = int(niche_id)
    mask = adata.obs["niche"].to_numpy().astype(int) == niche_id
    if not mask.any():
        raise ValueError(f"niche_id {niche_id} not present")

    return Niche(
        niche_id=niche_id,
        name="",
        composition=_composition(adata, mask),
        marker_program=_marker_program(adata, niche_id, n_markers),
        prognostic=None,
    )


def characterize_all_niches(adata: AnnData) -> list[Niche]:
    """Characterize every niche present in adata.obs['niche']."""
    if "niche" not in adata.obs.columns:
        raise ValueError(
            "adata.obs['niche'] missing; run engine.niches.find_niches first"
        )
    ids = sorted({int(n) for n in adata.obs["niche"].to_numpy()})
    return [characterize_niche(adata, niche_id) for niche_id in ids]


def _composition(adata: AnnData, mask: np.ndarray) -> dict[str, float]:
    """cell_type -> fraction within the niche."""
    types = adata.obs["cell_type"].astype(str).to_numpy()[mask]
    values, counts = np.unique(types, return_counts=True)
    total = counts.sum()
    return {str(v): round(float(c) / total, 4) for v, c in zip(values, counts)}


def _marker_program(adata: AnnData, niche_id: int, n: int) -> list[str]:
    """Top differential markers for the niche via scanpy.tl.rank_genes_groups.

    Falls back to top mean-intensity markers if the ranking cannot be computed
    (for example a single niche present).
    """
    col = "_niche_str"
    adata.obs[col] = adata.obs["niche"].astype(int).astype(str).astype("category")
    if adata.obs[col].cat.categories.size < 2:
        return _top_mean_markers(adata, niche_id, n)
    try:
        sc.tl.rank_genes_groups(
            adata,
            groupby=col,
            groups=[str(niche_id)],
            reference="rest",
            method="wilcoxon",
            key_added=_RANK_KEY,
        )
        names = adata.uns[_RANK_KEY]["names"][str(niche_id)]
        return [str(x) for x in list(names)[:n]]
    except Exception:
        return _top_mean_markers(adata, niche_id, n)


def _top_mean_markers(adata: AnnData, niche_id: int, n: int) -> list[str]:
    mask = adata.obs["niche"].to_numpy().astype(int) == int(niche_id)
    x = adata.X[mask]
    x = np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)
    if x.shape[0] == 0:
        return []
    order = np.argsort(x.mean(axis=0))[::-1][:n]
    markers = list(adata.var_names)
    return [str(markers[i]) for i in order]
