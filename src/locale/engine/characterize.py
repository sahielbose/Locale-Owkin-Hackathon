"""Niche characterization (Lane A). Backed by scanpy.

For each niche: its cell-type composition and its top enriched markers (the
"marker program"). Returns Niche objects with name left blank for interpret.py.
"""

from __future__ import annotations

from anndata import AnnData

from ..schema import Niche


def characterize_niche(adata: AnnData, niche_id: int) -> Niche:
    """Composition + marker program for one niche.

    Args:
        adata: canonical AnnData with adata.obs['niche'] assigned
            (engine.niches.find_niches).
        niche_id: the niche to characterize.

    Returns:
        Niche(niche_id=niche_id, name="", composition=<cell_type->fraction>,
        marker_program=<top enriched markers>, prognostic=None). The name is
        filled later by mcp_server.interpret.

    TODO(Lane A):
        1. composition: value_counts of obs['cell_type'] within the niche, normalized.
        2. marker program: scanpy.tl.rank_genes_groups(adata, "niche",
           groups=[str(niche_id)], method="wilcoxon") and take the top marker names.
        3. return the Niche (prognostic stays None; outcome.py adds it).
    """
    raise NotImplementedError(
        "characterize_niche: composition + scanpy.rank_genes_groups marker program."
    )


def characterize_all_niches(adata: AnnData) -> list[Niche]:
    """Characterize every niche present in adata.obs['niche'].

    TODO(Lane A): iterate the unique niche ids and call characterize_niche.
    """
    raise NotImplementedError(
        "characterize_all_niches: loop unique obs['niche'] -> characterize_niche."
    )
