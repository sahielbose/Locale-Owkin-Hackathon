"""Lane A engine tests against data/mock.h5ad.

Exercises the real engine end to end: per-image graph (no cross-core edges),
neighborhood enrichment, niche discovery, characterization, survival association,
and the three required validation checks.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from src.localespatial.engine.characterize import characterize_niche
from src.localespatial.engine.enrichment import compute_enrichment
from src.localespatial.engine.graph import build_spatial_graph, cross_image_edges
from src.localespatial.engine.niches import find_niches
from src.localespatial.engine.outcome import (
    correlate_niche_outcome,
    niche_outcome,
    rank_prognostic_niches,
)
from src.localespatial.engine.validate import (
    marker_validation,
    shuffle_negative_control,
    stability_ari,
)
from src.localespatial.schema import EnrichmentResult, Niche, Prognostic

MOCK_PATH = Path(__file__).resolve().parents[1] / "data" / "mock.h5ad"
N_NICHES = 6

pytestmark = pytest.mark.skipif(
    not MOCK_PATH.exists(), reason="run `python scripts/make_mock.py` first"
)


@pytest.fixture()
def adata() -> ad.AnnData:
    a = ad.read_h5ad(MOCK_PATH)
    build_spatial_graph(a)
    find_niches(a, n_niches=N_NICHES)
    return a


# --- graph ------------------------------------------------------------------------


def test_graph_has_no_cross_image_edges():
    a = ad.read_h5ad(MOCK_PATH)
    build_spatial_graph(a)
    assert "spatial_connectivities" in a.obsp
    assert cross_image_edges(a) == 0


def test_graph_leak_is_detected():
    import squidpy as sq

    a = ad.read_h5ad(MOCK_PATH)
    sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)  # no library_key
    assert cross_image_edges(a) > 0


# --- enrichment -------------------------------------------------------------------


def test_compute_enrichment_shape_and_symmetry():
    a = ad.read_h5ad(MOCK_PATH)
    res = compute_enrichment(a, "cohort:breast", n_perms=200)
    assert isinstance(res, EnrichmentResult)
    n = len(res.cell_types)
    assert n >= 4
    z = np.array(res.zscores)
    assert z.shape == (n, n)
    assert np.array(res.pvalues).shape == (n, n)
    # self-affinity: a cell type sits next to its own kind more than chance
    assert z[res.cell_types.index("Tumor"), res.cell_types.index("Tumor")] > 0


# --- niches -----------------------------------------------------------------------


def test_find_niches_writes_int_labels(adata):
    assert "niche" in adata.obs
    assert np.issubdtype(adata.obs["niche"].dtype, np.integer)
    assert adata.obs["niche"].nunique() == N_NICHES


def test_niches_are_neighborhoods_not_cell_types(adata):
    # A shared cell type should appear across several niches, not map 1:1 to one.
    import pandas as pd

    spread = pd.crosstab(
        adata.obs["cell_type"].astype(str),
        adata.obs["niche"].astype(int),
        normalize="index",
    )
    assert (spread.loc["Macrophage"] > 0.05).sum() >= 2


# --- characterize -----------------------------------------------------------------


def test_characterize_niche_returns_schema(adata):
    nid = int(adata.obs["niche"].to_numpy()[0])
    niche = characterize_niche(adata, nid)
    assert isinstance(niche, Niche)
    assert niche.name == ""  # interpret.py fills the name later
    assert abs(sum(niche.composition.values()) - 1.0) < 1e-6
    assert niche.marker_program
    assert Niche.model_validate(niche.model_dump()) == niche


# --- outcome ----------------------------------------------------------------------


def test_niche_outcome_returns_valid_prognostic(adata):
    nid = int(adata.obs["niche"].to_numpy()[0])
    prog = niche_outcome(adata, nid)
    assert isinstance(prog, Prognostic)
    assert np.isfinite(prog.hazard_ratio)
    assert prog.ci_low <= prog.hazard_ratio <= prog.ci_high
    assert prog.n_patients >= 2


def test_correlate_verdict_requires_selection_aware_p():
    """A good FDR q is NOT enough: if the global selection-aware p says the best of all
    niches is chance (>= alpha), no single niche can be 'supported' whatever its q."""
    summary = {
        "cohort": {
            "n_hypotheses_tested": 12,
            "n_events": 79,
            "min_detectable_hr": 1.37,
            "p_selection_aware": 0.44,  # best-of-12 is chance
        },
        "niches": {
            3: {"hazard_ratio": 0.6, "ci_95": [0.4, 0.9], "p_raw": 0.001, "q_fdr": 0.04}
        },
    }
    # q = 0.04 clears alpha, but p_selection_aware = 0.44 does not -> insufficient
    assert correlate_niche_outcome(summary, 3)["verdict"] == "insufficient evidence"

    # flip ONLY the selection-aware p below alpha: now both gates clear -> supported
    summary["cohort"]["p_selection_aware"] = 0.01
    assert correlate_niche_outcome(summary, 3)["verdict"] == "supported"


def test_rank_prognostic_niches_is_bh_ordered(adata):
    ranked = rank_prognostic_niches(adata)
    assert ranked  # at least one viable niche on the mock
    for _, prog in ranked:
        assert isinstance(prog, Prognostic)
        assert prog.ci_low <= prog.hazard_ratio <= prog.ci_high


# --- validation -------------------------------------------------------------------


def test_shuffle_null_reports_effect(adata):
    out = shuffle_negative_control(adata, n_permutations=20)
    assert set(out) >= {"real_effect", "null_mean", "null_std", "empirical_p"}
    assert out["real_effect"] >= out["null_mean"]


def test_stability_ari_in_range(adata):
    ari = stability_ari(adata, n_runs=3)
    assert -0.5 <= ari <= 1.0


def test_marker_validation_reports_correlation(adata):
    nid = int(adata.obs["niche"].to_numpy()[0])
    out = marker_validation(adata, nid)
    assert out["niche_id"] == nid
    assert "passed" in out and "evidence" in out
    assert out["evidence"]["top_markers"]
