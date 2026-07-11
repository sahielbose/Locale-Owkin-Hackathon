"""Spatial graph construction (Lane A). Backed by squidpy.

The spatial graph is built PER image_id: cells are only neighbors within the same
IMC core. Everything downstream (enrichment, niches) reads adata.obsp.
"""

from __future__ import annotations

from anndata import AnnData


def build_spatial_graph(
    adata: AnnData,
    image_id: str | None = None,
    method: str = "delaunay",
    n_neighs: int = 6,
    radius: float | None = None,
) -> AnnData:
    """Build the spatial neighbor graph with squidpy.gr.spatial_neighbors.

    Args:
        adata: canonical AnnData (needs obsm['spatial'] and obs['image_id']).
        image_id: if given, restrict to that one core; else build per-core using
            squidpy's library_key='image_id' so cores are not cross-connected.
        method: "delaunay" (default), "knn", or "radius".
        n_neighs: neighbors for the knn graph.
        radius: cutoff (microns) for the radius graph.

    Returns:
        The AnnData with adata.obsp['spatial_connectivities'] and
        adata.obsp['spatial_distances'] populated.

    TODO(Lane A):
        1. subset to image_id if provided.
        2. call squidpy.gr.spatial_neighbors(adata, coord_type="generic",
           library_key="image_id", delaunay=(method=="delaunay"),
           n_neighs=n_neighs, radius=radius).
        3. return adata.
    """
    raise NotImplementedError(
        "build_spatial_graph: wire squidpy.gr.spatial_neighbors (per image_id)."
    )
