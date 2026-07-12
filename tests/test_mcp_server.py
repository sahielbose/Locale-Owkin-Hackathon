"""Lane B server tests: call each tool directly against data/mock.h5ad and assert
every return validates against its Pydantic schema; check the orchestrator ranking,
the elicitation resolver, and interpret.py's fallback + (monkeypatched) API path.

No real Anthropic calls are made (the client is monkeypatched).
"""

from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

from src.locale.mcp_server import interpret, server, tools
from src.locale.schema import (
    EnrichmentResult,
    MapPayload,
    Niche,
    Prognostic,
    SampleRecord,
)

MOCK_PATH = Path(__file__).resolve().parents[1] / "data" / "mock.h5ad"

pytestmark = pytest.mark.skipif(
    not MOCK_PATH.exists(), reason="run `python scripts/make_mock.py` first"
)


@pytest.fixture(autouse=True)
def _use_mock(monkeypatch):
    """Force the tools to serve the committed mock, hermetically."""
    monkeypatch.setenv("LOCALE_DATA", str(MOCK_PATH))
    tools.reset_cache()
    yield
    tools.reset_cache()


def _revalidates(obj) -> None:
    """A schema object must survive a model_dump -> model_validate round-trip."""
    model = type(obj)
    assert model.model_validate(obj.model_dump()) == obj


# --- individual tools -------------------------------------------------------------


def test_list_samples():
    samples = tools.list_samples()
    assert samples and all(isinstance(s, SampleRecord) for s in samples)
    for s in samples:
        _revalidates(s)
        assert s.n_cells > 0 and s.cohort == "breast"


def test_describe_sample_cohort_and_image():
    cohort = tools.describe_sample()
    assert isinstance(cohort, SampleRecord)
    assert cohort.image_id is None and cohort.n_cells > 0
    _revalidates(cohort)

    image_id = tools.list_samples()[0].image_id
    one = tools.describe_sample(image_id=image_id)
    assert one.image_id == image_id and one.n_cells > 0


def test_compute_enrichment_square_matrix():
    res = tools.compute_enrichment("cohort:breast")
    assert isinstance(res, EnrichmentResult)
    _revalidates(res)
    n = len(res.cell_types)
    assert n >= 4
    assert len(res.zscores) == n and all(len(row) == n for row in res.zscores)
    assert len(res.pvalues) == n and all(len(row) == n for row in res.pvalues)


def test_find_and_characterize_niches():
    niches = tools.find_niches()
    assert niches and all(isinstance(n, Niche) for n in niches)
    for n in niches:
        _revalidates(n)
        assert n.composition and abs(sum(n.composition.values()) - 1.0) < 0.05
        assert n.marker_program and n.name

    one = tools.characterize_niche(niches[0].niche_id)
    assert isinstance(one, Niche) and one.niche_id == niches[0].niche_id


def test_find_prognostic_niches_ranked():
    ranked = tools.find_prognostic_niches()
    assert ranked and all(isinstance(n, Niche) for n in ranked)
    for n in ranked:
        _revalidates(n)

    # niches WITH a prognostic must come before those without
    has_prog = [n.prognostic is not None for n in ranked]
    assert has_prog == sorted(has_prog, reverse=True)

    # among prognostic niches, hazard ratio is non-increasing (worst survival first)
    hrs = [n.prognostic.hazard_ratio for n in ranked if n.prognostic is not None]
    assert hrs == sorted(hrs, reverse=True)

    # any prognostic present must be finite and schema-valid
    for n in ranked:
        if n.prognostic is not None:
            p = n.prognostic
            assert isinstance(p, Prognostic)
            assert p.ci_low <= p.hazard_ratio <= p.ci_high
            assert p.n_patients >= 1


def test_find_prognostic_niches_patient_subset():
    ranked = tools.find_prognostic_niches(patient_subset=["P001"])
    assert ranked and all(isinstance(n, Niche) for n in ranked)


def test_get_map_payload_both_modes():
    image_id = tools.list_samples()[0].image_id
    for mode in ("cell_type", "niche"):
        payload = tools.get_map_payload(image_id, mode)
        assert isinstance(payload, MapPayload)
        _revalidates(payload)
        assert payload.color_mode == mode and payload.image_id == image_id
        assert payload.units and payload.legend


def test_all_tools_report_a_backend():
    # exercise everything, then confirm each tool recorded which path served it
    tools.list_samples()
    tools.compute_enrichment()
    tools.find_niches()
    tools.find_prognostic_niches()
    tools.get_map_payload(tools.list_samples()[0].image_id)
    status = tools.backend_status()
    assert status["list_samples"] == "real"
    # Lane A engine is implemented now, so these run the real engine
    assert status["find_niches"] == "real"
    assert status["compute_enrichment"] == "real"


# --- interpret.py -----------------------------------------------------------------


def test_interpret_fallback_without_key(monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    out = interpret.name_and_describe_niche(
        {"Tumor": 0.7, "Fibroblast": 0.2, "CD8_T": 0.02}, ["panCK", "SMA"]
    )
    assert out["name"] and out["description"]
    assert isinstance(interpret.name_niche({"Tumor": 0.9}, ["panCK"]), str)


def test_interpret_empty_composition_never_crashes(monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    out = interpret.name_and_describe_niche({}, [])
    assert out["name"] and out["description"]


def test_interpret_uses_api_when_key_present(monkeypatch):
    """With a key set, interpret must call the (mocked) client and parse its JSON."""
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key-not-real")

    class _Block:
        type = "text"
        text = '{"name": "test tumor niche", "description": "A mocked description."}'

    class _Msg:
        content = [_Block()]

    class _Messages:
        def create(self, **kwargs):
            # ensure thinking is disabled and no real network is used
            assert kwargs.get("thinking") == {"type": "disabled"}
            return _Msg()

    class _Client:
        def __init__(self, **kwargs):
            self.messages = _Messages()

    import anthropic

    monkeypatch.setattr(anthropic, "Anthropic", _Client)
    out = interpret.name_and_describe_niche({"Tumor": 0.8}, ["panCK"])
    assert out["name"] == "test tumor niche"
    assert out["description"] == "A mocked description."


def test_interpret_falls_back_on_api_error(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key-not-real")

    class _Boom:
        def __init__(self, **kwargs):
            raise RuntimeError("network down")

    import anthropic

    monkeypatch.setattr(anthropic, "Anthropic", _Boom)
    out = interpret.name_and_describe_niche({"Tumor": 0.9}, ["panCK"])
    assert out["name"]  # deterministic fallback, no crash


# --- server elicitation resolver --------------------------------------------------


def test_resolve_cohort_defaults_without_context():
    cohort, n = asyncio.run(server._resolve_cohort(None, None, None))
    assert cohort == "breast"


def test_resolve_cohort_ambiguous_defaults_without_context():
    cohort, n = asyncio.run(server._resolve_cohort(None, "spleen", 4))
    assert cohort == "breast" and n == 4


def test_resolve_cohort_elicits_when_ambiguous():
    class _Result:
        action = "accept"
        data = server.NicheQuery(cohort="breast", n_niches=9)

    class _Ctx:
        async def elicit(self, message, schema):
            return _Result()

    cohort, n = asyncio.run(server._resolve_cohort(_Ctx(), None, None))
    assert cohort == "breast" and n == 9


def test_resolve_cohort_skips_elicit_when_known():
    called = {"n": 0}

    class _Ctx:
        async def elicit(self, message, schema):
            called["n"] += 1
            raise AssertionError("should not elicit for a known cohort")

    cohort, n = asyncio.run(server._resolve_cohort(_Ctx(), "breast", 5))
    assert cohort == "breast" and n == 5 and called["n"] == 0


def test_server_registers_all_seven_tools():
    names = {t.name for t in asyncio.run(server.mcp.list_tools())}
    assert names == {
        "list_samples",
        "describe_sample",
        "compute_enrichment",
        "find_niches",
        "characterize_niche",
        "find_prognostic_niches",
        "get_map_payload",
        "describe_niches",
        "correlate_niche_outcome",
    }
