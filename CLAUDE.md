# CLAUDE.md: Locale

Context for Claude Code. This is the build source of truth. The pitch/vision lives in `Locale.md`; when the two disagree, this file wins (in particular, the dataset changed, see below).

## What we are building

Locale is an MCP server for K Pro (Owkin's AI Scientist), built for the Rewiring Biology Hackathon. It reasons about *where* cells sit in tumor tissue: it finds recurring cellular niches (for example an immune-excluded tumor core) and links them to patient survival. It gives K Pro a spatial-reasoning capability that its current expression-only tools lack.

Track being targeted: Best AI Scientist MCP.

## Team and parallelism (read this first)

3-person team, three lanes that build in parallel:
- Lane A (data + analysis engine): `src/localespatial/data/` and `src/localespatial/engine/`. Critical path, needs biology judgment.
- Lane B (MCP server): `src/localespatial/mcp_server/`.
- Lane C (visualization): `src/localespatial/viz/`.

The linchpin that lets the three lanes work independently is the shared data contract in `src/localespatial/schema.py` plus a committed tiny `data/mock.h5ad`. Lanes B and C build against the schema + mock data; Lane A produces the real data and analysis; everything integrates through the schema objects. **Never change a schema field without telling the team.**

## Dataset (CHANGED from MOSAIC, important)

We do NOT use MOSAIC (no access). We use a public dataset:

**Jackson & Fischer et al. 2020, "The single-cell pathology landscape of breast cancer," Nature.** Imaging mass cytometry (IMC): 35 protein markers, 720 images, 352 breast cancer patients, survival data for 281 of them, roughly 1.7M cells. Each cell has x/y coordinates (microns), 35 marker intensities, and a cell-type label (PhenoGraph phenotypes). Clinical variables: tumor grade, ER/PR/HER2 status, overall survival, disease-free survival, alive/dead status.

Why it fits Locale: multi-patient (ideal for CellCharter's cross-sample niche detection), true single-cell with coordinates (so NO deconvolution is needed), cell types already annotated, real well-powered survival (281 patients), and it is the field's benchmark dataset (lots of reference code exists).

Consequences vs the old MOSAIC plan (which some docs still reference):
- Single-cell, so there is NO deconvolution step. Delete that from your mental model.
- Cell types are pre-annotated, so no annotation step is strictly required.
- It is protein markers (35), not transcriptome, so niche "programs" are marker-based (call them marker programs, not gene programs).
- Cancer type is breast; cell types are breast-TME (tumor/epithelial subsets, CD8 T, CD4 T, B cells, macrophages, endothelial, fibroblasts/stroma, etc.).

### Getting the data (data is NOT committed to git)

Source: Zenodo record `3518284`, file `SingleCell_and_Metadata.zip` (36.8 GB, do NOT download the whole thing). We fetch only the files we need via `remotezip` (HTTP range requests), namely: the single-cell data table, the cell-type/cluster label files, and the patient-metadata table. Exact internal filenames must be confirmed by listing the archive (`infolist()`) first, then extracting only those. See `scripts/download_data.py`.

Data-not-in-git rules (enforce these):
- `data/` is gitignored. Never commit raw CSVs or `data/locale.h5ad`.
- The ONE committed data file is `data/mock.h5ad` (tiny synthetic object for building against).
- The real processed object `data/locale.h5ad` is shared between teammates out of band (shared Drive), not via git. Also put the extracted CSVs in the shared Drive so nobody re-downloads from Zenodo.

## The canonical AnnData object (Lane A produces this -> `data/locale.h5ad`)

Every downstream component reads this. `src/localespatial/data/build_anndata.py` builds it from the extracted CSVs.

```
adata.X                  float32 [n_cells x 35]   arcsinh(cofactor=5) marker intensities, z-scored per marker
adata.var_names          35 marker names
adata.obs['cell_type']   category    cell phenotype
adata.obs['patient_id']  str
adata.obs['image_id']    str         IMC core / image id (spatial graph is built per image_id)
adata.obs['os_month']    float       overall survival time (months)
adata.obs['os_event']    int {0,1}   death event
adata.obs['dfs_month']   float       disease-free survival (optional)
adata.obs['dfs_event']   int {0,1}   (optional)
adata.obs['grade','er','pr','her2','subtype']   clinical (optional but useful)
adata.obsm['spatial']    float [n_cells x 2]     x, y in microns
adata.uns['markers']     list[str]
```

Preprocessing in `build_anndata.py`: arcsinh transform (cofactor 5), optional clip at 99th percentile, z-score per marker. Join cell-type labels and patient survival/clinical onto the cells. Set `obsm['spatial']` from the coordinate columns.

## The shared data contract (`src/localespatial/schema.py`, Pydantic v2)

This is THE coordination artifact. Implement exactly this.

```python
from pydantic import BaseModel

class SampleRecord(BaseModel):
    cohort: str
    patient_id: str | None = None
    image_id: str | None = None
    n_cells: int
    cell_types: list[str]
    has_survival: bool

class EnrichmentResult(BaseModel):
    scope: str                       # e.g. "cohort:breast" or "image:<id>"
    cell_types: list[str]
    zscores: list[list[float]]       # cell_type x cell_type
    pvalues: list[list[float]]

class KMCurve(BaseModel):
    time: list[float]
    high: list[float]                # survival prob, high niche-abundance group
    low: list[float]                 # survival prob, low group

class Prognostic(BaseModel):
    hazard_ratio: float
    ci_low: float
    ci_high: float
    pvalue: float
    n_patients: int
    km: KMCurve | None = None

class Niche(BaseModel):
    niche_id: int
    name: str                        # human-readable, filled by interpret.py
    composition: dict[str, float]    # cell_type -> fraction
    marker_program: list[str]        # top enriched markers
    prognostic: Prognostic | None = None

class MapUnit(BaseModel):
    x: float
    y: float
    cell_type: str
    niche_id: int | None = None

class MapPayload(BaseModel):
    units: list[MapUnit]
    legend: dict[str, str]           # label -> hex color
    color_mode: str                  # "cell_type" | "niche"
    image_id: str | None = None
```

## Repo layout

```
locale/
  CLAUDE.md
  README.md
  .gitignore
  requirements.txt
  data/
    .gitkeep
    mock.h5ad              # committed, tiny synthetic (build against this)
    raw/                   # extracted CSVs (gitignored)
    locale.h5ad            # real canonical object (gitignored, shared out of band)
  scripts/
    download_data.py       # remotezip: list + extract only needed files
    make_mock.py           # generate committable data/mock.h5ad
  src/localespatial/
    __init__.py
    schema.py              # shared data contract (above)
    data/build_anndata.py  # CSVs -> canonical AnnData -> data/locale.h5ad
    engine/
      graph.py             # spatial graph (squidpy.gr.spatial_neighbors, per image_id)
      enrichment.py        # neighborhood enrichment -> EnrichmentResult (squidpy.gr.nhood_enrichment)
      niches.py            # niche detection -> adata.obs['niche'] (CellCharter VAE+GMM, or kmeans on neighborhood composition)
      characterize.py      # composition + marker program -> list[Niche] (scanpy rank_genes_groups)
      outcome.py           # niche abundance vs survival -> Prognostic (lifelines Cox/KM, FDR corrected)
      validate.py          # shuffle negative control, stability (ARI), marker validation
    mcp_server/
      server.py            # MCP server scaffold (Python MCP SDK, remote/HTTP)
      tools.py             # tool fns wrapping engine, returning schema objects (start on mock.h5ad)
      interpret.py         # niche naming via Anthropic API (server-side sampling)
    viz/
      app/                 # HTML/JS tissue-map widget (deck.gl or plotly scattergl)
      payload.py           # AnnData -> MapPayload helper
  tests/
    test_schema.py         # round-trip every model; load mock.h5ad and assert canonical schema
```

## MCP tools (server exposes these; each returns schema objects)

- `list_samples() -> list[SampleRecord]`
- `describe_sample(image_id | cohort) -> SampleRecord`
- `compute_enrichment(scope) -> EnrichmentResult`
- `find_niches(cohort, n_niches?) -> list[Niche]`
- `characterize_niche(niche_id) -> Niche`
- `find_prognostic_niches(cohort, patient_subset?) -> list[Niche]`   # orchestrator: niches + survival + naming, ranked
- `get_map_payload(image_id, color_mode) -> MapPayload`

Connector setup mirrors Owkin's Pathology Explorer (custom connector, remote MCP URL). Underspecified requests use elicitation (ask cohort / cell types / n_niches).

## Tech stack and conventions

Python 3.11. Libraries: squidpy, cellcharter, scanpy, anndata, lifelines, pydantic (v2), mcp (Python SDK), anthropic, remotezip, numpy, pandas, plotly (map; deck.gl via CDN in the HTML is fine), pytest.

Conventions: type hints everywhere, docstrings that name the underlying tool, black formatting. Engine functions take the canonical AnnData and return schema objects. Do NOT put analysis logic in the MCP layer; it only wraps the engine. Keep the schema stable. The validation functions (shuffle control, stability, marker check) are required, not optional; they are how we prove the niches are real.

Style: no em-dashes in any generated text or docs.

## Commands

```
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python scripts/make_mock.py               # writes data/mock.h5ad (committed)
python scripts/download_data.py           # lists the zip, extracts needed CSVs to data/raw/ (run once)
python -m src.localespatial.data.build_anndata   # builds data/locale.h5ad from data/raw/
pytest -q                                 # schema + mock checks
python -m src.localespatial.mcp_server.server    # run the MCP server
```

## Status

Scaffold, Lane B (MCP server), and Lane C (viz) have landed. Lane A (analysis engine) is now implemented: per-image spatial graph, neighborhood enrichment, niche discovery (neighborhood-composition k-means), characterization + marker programs, niche-to-survival association (Cox/KM, BH-ranked), and the three validation checks all run against `data/mock.h5ad` and return the schema objects. Point `LOCALE_DATA` (or drop in `data/locale.h5ad`) to run the same engine on the real cohort.
