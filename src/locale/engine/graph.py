"""Spatial graph construction (Lane A). Backed by squidpy.

The spatial graph is built PER image_id: cells are only neighbors within the same
IMC core, because each core is its own tissue with its own coordinate system.
Connecting cells across cores fabricates edges and every downstream result is
garbage, and it fails silently, so cross_image_edges() exists as a hard guard.
Everything downstream (enrichment, niches) reads adata.obsp.
"""

from __future__ import annotations

import numpy as np
import squidpy as sq
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
        image_id: if given, restrict to that one core (returns a graph-bearing copy
            of just those cells); else build per-core in place using squidpy's
            library_key='image_id' so cores are never cross-connected.
        method: "delaunay" (default), "knn", or "radius".
        n_neighs: neighbors for the knn graph.
        radius: cutoff (microns) for the radius graph.

    Returns:
        The AnnData with adata.obsp['spatial_connectivities'] and
        adata.obsp['spatial_distances'] populated.
    """
    if "spatial" not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] is required to build the graph")
    if "image_id" not in adata.obs.columns:
        raise ValueError("adata.obs['image_id'] is required (graph is built per image)")

    target = adata
    if image_id is not None:
        sel = adata.obs["image_id"].astype(str) == str(image_id)
        if not sel.any():
            raise ValueError(f"image_id {image_id!r} not found")
        target = adata[sel.to_numpy()].copy()

    kwargs: dict = {"coord_type": "generic", "library_key": "image_id"}
    if method == "delaunay":
        kwargs["delaunay"] = True
    elif method == "knn":
        kwargs["n_neighs"] = n_neighs
    elif method == "radius":
        if radius is None:
            raise ValueError("method='radius' requires a radius (microns)")
        kwargs["radius"] = radius
    else:
        raise ValueError(f"unknown method {method!r} (use delaunay, knn, or radius)")

    sq.gr.spatial_neighbors(target, **kwargs)
    target.uns["graph_params"] = {"method": method, "library_key": "image_id"}

    n_bad = cross_image_edges(target)
    if n_bad:
        raise AssertionError(
            f"{n_bad} edges cross an image boundary; the graph leaked across cores"
        )
    return target


def cross_image_edges(adata: AnnData, image_key: str = "image_id") -> int:
    """Count graph edges that connect cells from two different images. Must be 0."""
    if "spatial_connectivities" not in adata.obsp:
        raise ValueError("build the spatial graph first")
    graph = adata.obsp["spatial_connectivities"].tocoo()
    images = adata.obs[image_key].to_numpy()
    return int((images[graph.row] != images[graph.col]).sum())
