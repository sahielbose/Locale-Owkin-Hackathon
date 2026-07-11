"""Smoke tests for the MCP tool layer against data/mock.h5ad.

Bonus (beyond the required test_schema.py): proves Lane B's tools return valid
schema objects on the mock before Lane A's engine exists. If these break, an
integration seam between schema, tools, and the mock has drifted.
"""

from __future__ import annotations

from pathlib import Path

import pytest

MOCK_PATH = Path(__file__).resolve().parents[1] / "data" / "mock.h5ad"

pytestmark = pytest.mark.skipif(
    not MOCK_PATH.exists(), reason="run `python scripts/make_mock.py` first"
)

from src.locale.mcp_server import tools  # noqa: E402
from src.locale.schema import (  # noqa: E402
    EnrichmentResult,
    MapPayload,
    Niche,
    SampleRecord,
)


def test_list_samples():
    samples = tools.list_samples()
    assert samples and all(isinstance(s, SampleRecord) for s in samples)
    assert all(s.n_cells > 0 for s in samples)


def test_describe_sample_cohort_and_image():
    cohort = tools.describe_sample()
    assert isinstance(cohort, SampleRecord)
    assert cohort.image_id is None and cohort.n_cells > 0

    image_id = tools.list_samples()[0].image_id
    one = tools.describe_sample(image_id=image_id)
    assert one.image_id == image_id and one.n_cells > 0


def test_compute_enrichment_shape():
    res = tools.compute_enrichment("cohort:breast")
    assert isinstance(res, EnrichmentResult)
    n = len(res.cell_types)
    assert n >= 4
    assert len(res.zscores) == n and all(len(row) == n for row in res.zscores)
    assert len(res.pvalues) == n and all(len(row) == n for row in res.pvalues)


def test_find_and_characterize_niches():
    niches = tools.find_niches()
    assert niches and all(isinstance(nch, Niche) for nch in niches)
    for nch in niches:
        assert nch.composition and abs(sum(nch.composition.values()) - 1.0) < 0.05
        assert nch.marker_program
        assert nch.name

    one = tools.characterize_niche(niches[0].niche_id)
    assert isinstance(one, Niche) and one.niche_id == niches[0].niche_id


def test_find_prognostic_niches_ranked():
    ranked = tools.find_prognostic_niches()
    assert ranked and all(isinstance(nch, Niche) for nch in ranked)


def test_get_map_payload_both_modes():
    image_id = tools.list_samples()[0].image_id
    for mode in ("cell_type", "niche"):
        payload = tools.get_map_payload(image_id, mode)
        assert isinstance(payload, MapPayload)
        assert payload.color_mode == mode
        assert payload.image_id == image_id
        assert payload.units and payload.legend
