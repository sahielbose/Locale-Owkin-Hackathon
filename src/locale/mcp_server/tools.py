"""MCP tool functions (Lane B).

Each function returns a schema object from src/locale/schema.py. Tools try the real
engine first and, while the engine still raises NotImplementedError, fall back to
lightweight reads/summaries of the loaded AnnData so the server is demoable on
data/mock.h5ad TODAY. Swap nothing here when Lane A lands; the try/except upgrades
automatically.

The only "analysis" done in this layer is a clearly-marked TEMPORARY neighborhood
co-occurrence used by compute_enrichment until engine.enrichment is wired. All real
analysis belongs in src/locale/engine/.
"""

from __future__ import annotations

import functools
import os
from pathlib import Path

import anndata as ad
import numpy as np

from ..engine import enrichment as _engine_enrichment
from ..engine import outcome as _engine_outcome
from ..schema import EnrichmentResult, MapPayload, Niche, Prognostic, SampleRecord
from ..viz.payload import build_map_payload
from . import interpret

COHORT = "breast"
_DEFAULT_DATA = Path(__file__).resolve().parents[3] / "data" / "mock.h5ad"
_TOP_MARKERS = 6


def data_path() -> Path:
    """Path to the AnnData the server serves (LOCALE_DATA env overrides the mock)."""
    return Path(os.environ.get("LOCALE_DATA", str(_DEFAULT_DATA)))


@functools.lru_cache(maxsize=1)
def _load() -> ad.AnnData:
    path = data_path()
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `python scripts/make_mock.py` or set LOCALE_DATA."
        )
    return ad.read_h5ad(path)


def _has_survival(obs) -> bool:
    return "os_month" in obs and bool(
        np.isfinite(obs["os_month"].to_numpy(dtype=float)).any()
    )


# --- sample introspection ---------------------------------------------------------


def list_samples() -> list[SampleRecord]:
    """List every sample (one record per IMC image) in the loaded cohort."""
    adata = _load()
    obs = adata.obs
    records: list[SampleRecord] = []
    for image_id in sorted(set(obs["image_id"].astype(str))):
        sub = obs[obs["image_id"].astype(str) == image_id]
        patient = str(sub["patient_id"].iloc[0]) if "patient_id" in sub else None
        records.append(
            SampleRecord(
                cohort=COHORT,
                patient_id=patient,
                image_id=image_id,
                n_cells=int(sub.shape[0]),
                cell_types=sorted(set(sub["cell_type"].astype(str))),
                has_survival=_has_survival(sub),
            )
        )
    return records


def describe_sample(
    image_id: str | None = None, cohort: str | None = None
) -> SampleRecord:
    """Describe one image, or the whole cohort when no image_id is given."""
    adata = _load()
    obs = adata.obs
    if image_id is not None:
        sub = obs[obs["image_id"].astype(str) == str(image_id)]
        if sub.shape[0] == 0:
            raise ValueError(f"image_id {image_id!r} not found")
        patient = str(sub["patient_id"].iloc[0]) if "patient_id" in sub else None
        return SampleRecord(
            cohort=COHORT,
            patient_id=patient,
            image_id=str(image_id),
            n_cells=int(sub.shape[0]),
            cell_types=sorted(set(sub["cell_type"].astype(str))),
            has_survival=_has_survival(sub),
        )
    return SampleRecord(
        cohort=cohort or COHORT,
        patient_id=None,
        image_id=None,
        n_cells=int(obs.shape[0]),
        cell_types=sorted(set(obs["cell_type"].astype(str))),
        has_survival=_has_survival(obs),
    )


# --- neighborhood enrichment ------------------------------------------------------


def compute_enrichment(scope: str = f"cohort:{COHORT}") -> EnrichmentResult:
    """Cell-type co-location matrix. Uses engine.compute_enrichment when available."""
    adata = _load()
    try:
        return _engine_enrichment.compute_enrichment(adata, scope)
    except NotImplementedError:
        return _mock_enrichment(adata, scope)


def _mock_enrichment(
    adata: ad.AnnData, scope: str, n_perm: int = 100, k: int = 6
) -> EnrichmentResult:
    """TEMPORARY stopgap co-location, replaced by engine.compute_enrichment (squidpy).

    Builds a per-image kNN graph, counts cell-type-pair adjacencies, and z-scores
    them against a within-image label-permutation null. Deterministic (fixed seed).
    """
    obs = adata.obs
    cell_types = sorted(set(obs["cell_type"].astype(str)))
    idx = {ct: i for i, ct in enumerate(cell_types)}
    labels = obs["cell_type"].astype(str).map(idx).to_numpy()
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    images = obs["image_id"].astype(str).to_numpy()
    n_types = len(cell_types)

    def adjacency_counts(lab: np.ndarray) -> np.ndarray:
        counts = np.zeros((n_types, n_types), dtype=float)
        for img in np.unique(images):
            sel = np.where(images == img)[0]
            if sel.size < 2:
                continue
            pts = coords[sel]
            kk = min(k, sel.size - 1)
            d2 = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1)
            neigh = np.argsort(d2, axis=1)[:, 1 : kk + 1]
            for local_i, nbrs in enumerate(neigh):
                a = lab[sel[local_i]]
                for j in nbrs:
                    b = lab[sel[j]]
                    counts[a, b] += 1.0
        return (counts + counts.T) / 2.0

    observed = adjacency_counts(labels)
    rng = np.random.default_rng(0)
    perms = np.zeros((n_perm, n_types, n_types), dtype=float)
    for p in range(n_perm):
        shuffled = labels.copy()
        for img in np.unique(images):
            sel = np.where(images == img)[0]
            shuffled[sel] = rng.permutation(shuffled[sel])
        perms[p] = adjacency_counts(shuffled)
    mean = perms.mean(0)
    std = perms.std(0)
    std[std == 0] = 1.0
    zscores = (observed - mean) / std
    pvalues = (np.sum(perms >= observed[None], axis=0) + 1.0) / (n_perm + 1.0)

    return EnrichmentResult(
        scope=f"mock:{scope}",
        cell_types=cell_types,
        zscores=zscores.tolist(),
        pvalues=pvalues.tolist(),
    )


# --- niches -----------------------------------------------------------------------


def _niche_ids(adata: ad.AnnData) -> list[int]:
    if "niche" not in adata.obs:
        raise ValueError(
            "adata.obs has no 'niche'; run engine.niches.find_niches first"
        )
    return sorted({int(n) for n in adata.obs["niche"].to_numpy()})


def _composition(adata: ad.AnnData, mask: np.ndarray) -> dict[str, float]:
    types = adata.obs["cell_type"].astype(str).to_numpy()[mask]
    if types.size == 0:
        return {}
    values, counts = np.unique(types, return_counts=True)
    total = counts.sum()
    return {str(v): round(float(c) / total, 4) for v, c in zip(values, counts)}


def _marker_program(
    adata: ad.AnnData, mask: np.ndarray, k: int = _TOP_MARKERS
) -> list[str]:
    """Top-k markers by mean (z-scored) intensity inside the mask."""
    x = adata.X[mask]
    x = np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)
    if x.shape[0] == 0:
        return []
    means = x.mean(axis=0)
    order = np.argsort(means)[::-1][:k]
    markers = list(adata.var_names)
    return [str(markers[i]) for i in order]


def _build_niche(adata: ad.AnnData, niche_id: int, mask: np.ndarray) -> Niche:
    composition = _composition(adata, mask)
    marker_program = _marker_program(adata, mask)
    name = interpret.name_niche(composition, marker_program, niche_id)
    return Niche(
        niche_id=niche_id,
        name=name,
        composition=composition,
        marker_program=marker_program,
        prognostic=None,
    )


def find_niches(cohort: str = COHORT, n_niches: int | None = None) -> list[Niche]:
    """Return the cohort's niches (from precomputed obs['niche'] until Lane A lands).

    n_niches is honored once engine.niches.find_niches can recluster; the mock has a
    fixed precomputed set, so n_niches is currently ignored (documented, not silent
    at the schema level).
    """
    adata = _load()
    niche_col = adata.obs["niche"].to_numpy()
    out: list[Niche] = []
    for niche_id in _niche_ids(adata):
        mask = niche_col == niche_id
        out.append(_build_niche(adata, niche_id, mask))
    return out


def characterize_niche(niche_id: int) -> Niche:
    """Composition + marker program + name for one niche."""
    adata = _load()
    mask = adata.obs["niche"].to_numpy() == niche_id
    if not mask.any():
        raise ValueError(f"niche_id {niche_id} not present")
    return _build_niche(adata, niche_id, mask)


def find_prognostic_niches(
    cohort: str = COHORT, patient_subset: list[str] | None = None
) -> list[Niche]:
    """Orchestrator: niches + survival association + naming, ranked.

    Attaches Prognostic via engine.outcome.niche_outcome when available; until then
    ranks by niche abundance and leaves prognostic=None (survival needs Lane A).
    """
    adata = _load()
    if patient_subset:
        keep = (
            adata.obs["patient_id"]
            .astype(str)
            .isin([str(p) for p in patient_subset])
            .to_numpy()
        )
        if not keep.any():
            raise ValueError(f"no cells for patient_subset {patient_subset}")
        adata = adata[keep]

    niche_col = adata.obs["niche"].to_numpy()
    total = niche_col.shape[0]
    niches: list[tuple[float, Niche]] = []
    for niche_id in _niche_ids(adata):
        mask = niche_col == niche_id
        niche = _build_niche(adata, niche_id, mask)
        try:
            niche.prognostic = _engine_outcome.niche_outcome(adata, niche_id)
        except NotImplementedError:
            niche.prognostic = None
        abundance = float(mask.sum()) / max(total, 1)
        niches.append((abundance, niche))

    def rank_key(item: tuple[float, Niche]) -> float:
        abundance, niche = item
        if niche.prognostic is not None:
            return -abs(np.log(max(niche.prognostic.hazard_ratio, 1e-6)))
        return -abundance  # fallback: most abundant first

    niches.sort(key=rank_key)
    return [niche for _, niche in niches]


# --- tissue map -------------------------------------------------------------------


def get_map_payload(image_id: str, color_mode: str = "cell_type") -> MapPayload:
    """Build the tissue-map render payload for one image (delegates to viz.payload)."""
    adata = _load()
    return build_map_payload(adata, image_id, color_mode)
