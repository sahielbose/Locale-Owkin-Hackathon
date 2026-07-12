"""Round-trip every Pydantic model, and assert data/mock.h5ad matches the
canonical AnnData schema (CLAUDE.md). This is the contract Lanes B and C build on.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from src.localespatial.schema import (
    EnrichmentResult,
    KMCurve,
    MapPayload,
    MapUnit,
    Niche,
    Prognostic,
    SampleRecord,
)

MOCK_PATH = Path(__file__).resolve().parents[1] / "data" / "mock.h5ad"
N_MARKERS = 35


# --- Pydantic round-trips ---------------------------------------------------------


def _roundtrip(model):
    """Serialize to JSON and parse back; the result must equal the original."""
    restored = type(model).model_validate_json(model.model_dump_json())
    assert restored == model
    return restored


def test_sample_record_roundtrip():
    _roundtrip(
        SampleRecord(
            cohort="breast",
            patient_id="P001",
            image_id="IMG001",
            n_cells=150,
            cell_types=["Tumor", "CD8_T", "Fibroblast"],
            has_survival=True,
        )
    )
    # optional fields default to None
    _roundtrip(
        SampleRecord(cohort="breast", n_cells=0, cell_types=[], has_survival=False)
    )


def test_enrichment_result_roundtrip():
    _roundtrip(
        EnrichmentResult(
            scope="cohort:breast",
            cell_types=["Tumor", "CD8_T"],
            zscores=[[3.1, -1.2], [-1.2, 2.0]],
            pvalues=[[0.001, 0.2], [0.2, 0.01]],
        )
    )


def test_km_curve_roundtrip():
    _roundtrip(
        KMCurve(
            time=[0.0, 12.0, 24.0],
            high=[1.0, 0.7, 0.4],
            low=[1.0, 0.95, 0.9],
        )
    )


def test_prognostic_roundtrip():
    km = KMCurve(time=[0.0, 12.0], high=[1.0, 0.6], low=[1.0, 0.9])
    _roundtrip(
        Prognostic(
            hazard_ratio=2.4,
            ci_low=1.3,
            ci_high=4.5,
            pvalue=0.003,
            n_patients=281,
            km=km,
        )
    )
    # km is optional
    _roundtrip(
        Prognostic(hazard_ratio=1.0, ci_low=0.8, ci_high=1.3, pvalue=0.4, n_patients=10)
    )


def test_niche_roundtrip():
    prog = Prognostic(
        hazard_ratio=2.0, ci_low=1.1, ci_high=3.6, pvalue=0.01, n_patients=100
    )
    _roundtrip(
        Niche(
            niche_id=0,
            name="immune-excluded tumor core",
            composition={"Tumor": 0.65, "Fibroblast": 0.22, "CD8_T": 0.02},
            marker_program=["panCK", "SMA", "CAIX"],
            prognostic=prog,
        )
    )
    _roundtrip(
        Niche(niche_id=1, name="", composition={}, marker_program=[], prognostic=None)
    )


def test_map_unit_roundtrip():
    _roundtrip(MapUnit(x=12.5, y=88.1, cell_type="Tumor", niche_id=0))
    _roundtrip(MapUnit(x=0.0, y=0.0, cell_type="CD8_T", niche_id=None))


def test_map_payload_roundtrip():
    _roundtrip(
        MapPayload(
            units=[
                MapUnit(x=1.0, y=2.0, cell_type="Tumor", niche_id=0),
                MapUnit(x=3.0, y=4.0, cell_type="CD8_T", niche_id=1),
            ],
            legend={"Tumor": "#d62728", "CD8_T": "#1f77b4"},
            color_mode="cell_type",
            image_id="IMG001",
        )
    )


# --- canonical AnnData schema (data/mock.h5ad) ------------------------------------


@pytest.fixture(scope="module")
def adata() -> ad.AnnData:
    if not MOCK_PATH.exists():
        pytest.skip(f"{MOCK_PATH} missing; run `python scripts/make_mock.py`")
    return ad.read_h5ad(MOCK_PATH)


def test_mock_X(adata):
    assert adata.X.dtype == np.float32
    assert adata.X.shape[1] == N_MARKERS
    assert adata.n_obs > 0
    # z-scored per marker
    assert np.allclose(adata.X.mean(axis=0), 0.0, atol=1e-3)
    assert np.allclose(adata.X.std(axis=0), 1.0, atol=1e-3)


def test_mock_var_and_markers(adata):
    assert len(adata.var_names) == N_MARKERS
    assert "markers" in adata.uns
    assert list(adata.uns["markers"]) == list(adata.var_names)
    assert len(adata.uns["markers"]) == N_MARKERS


def test_mock_obs_required_columns(adata):
    obs = adata.obs
    for col in ("cell_type", "patient_id", "image_id", "os_month", "os_event"):
        assert col in obs, f"missing obs column {col}"
    assert str(obs["cell_type"].dtype) == "category"
    assert np.issubdtype(obs["os_month"].dtype, np.floating)
    assert set(np.unique(obs["os_event"].to_numpy())).issubset({0, 1})


def test_mock_obs_niche_precomputed(adata):
    assert "niche" in adata.obs
    assert np.issubdtype(adata.obs["niche"].dtype, np.integer)
    assert adata.obs["niche"].nunique() >= 2


def test_mock_spatial(adata):
    assert "spatial" in adata.obsm
    assert adata.obsm["spatial"].shape == (adata.n_obs, 2)


def test_mock_multi_image_and_types(adata):
    assert adata.obs["image_id"].nunique() >= 2
    assert adata.obs["cell_type"].nunique() >= 4
