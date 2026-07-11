# Locale: Owkin Rewiring Bio Hackathon

An MCP server for K Pro (Owkin's AI Scientist) that reasons about *where* cells sit
in tumor tissue. It finds recurring cellular niches (for example an immune-excluded
tumor core) and links them to patient survival, giving K Pro a spatial-reasoning
capability its expression-only tools lack.

Built for the Rewiring Biology Hackathon. Track: Best AI Scientist MCP.

Dataset: Jackson & Fischer et al. 2020, "The single-cell pathology landscape of
breast cancer" (Nature). Imaging mass cytometry, 35 protein markers, 720 images,
352 patients, survival for 281, roughly 1.7M cells. See `CLAUDE.md` for the full
build spec and `Locale.md` for the vision/pitch. **CLAUDE.md is the source of truth.**

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python scripts/make_mock.py     # writes the committed data/mock.h5ad
pytest                          # schema + mock + tool smoke tests
```

Lanes B and C need only a light install (`pydantic`, `anndata`, `mcp`, `plotly`);
the heavy spatial stack (`squidpy`, `cellcharter`, `scanpy`, `lifelines`) is Lane A's
and is only required to build/analyze the real data.

Run the pieces:

```bash
python -m src.locale.mcp_server.server        # MCP server (streamable-http on :8000)
python -m src.locale.viz.payload              # regenerate the viz sample from the mock
open src/locale/viz/app/index.html            # interactive tissue map
```

## Data is NOT in git

The real dataset is 36.8 GB. Do not download or commit it.

- `data/` is gitignored. Never commit raw CSVs or `data/locale.h5ad`.
- The ONE committed data file is `data/mock.h5ad`: a tiny synthetic object (600
  cells, 4 images, 6 cell types, 4 niches) that matches the canonical schema, so
  Lanes B and C can build before the real data exists.
- Fetch the real data with `python scripts/download_data.py`, which uses `remotezip`
  HTTP range requests to pull ONLY the needed files from Zenodo record 3518284 (it
  first prints the archive listing so you can confirm the exact filenames). Run it
  once, by hand. It is never run by tests or CI.
- Share the extracted CSVs and the resulting `data/locale.h5ad` out of band (Drive),
  so nobody re-downloads from Zenodo.

## How the three lanes work

The linchpin is the shared data contract in `src/locale/schema.py` (Pydantic v2)
plus the committed `data/mock.h5ad`. All three lanes integrate through the schema
objects and the mock. **Never change a schema field without telling the team.**

### Lane A: data + analysis engine (`src/locale/data/`, `src/locale/engine/`)
Critical path, needs biology judgment. Produces the canonical AnnData
(`data/locale.h5ad`) and the real analysis. Each engine function takes the canonical
AnnData and returns schema objects; they are stubbed with `NotImplementedError` today.

- `data/build_anndata.py`: CSVs -> canonical AnnData (arcsinh cofactor 5, optional
  99th-pct clip, per-marker z-score; join cell types + survival; set `obsm['spatial']`).
  The preprocessing math is implemented; only the real column names are TODO.
- `engine/graph.py`: spatial graph per `image_id` (squidpy).
- `engine/enrichment.py`: neighborhood enrichment -> `EnrichmentResult` (squidpy).
- `engine/niches.py`: niche detection -> `obs['niche']` (CellCharter, kmeans fallback).
- `engine/characterize.py`: composition + marker program -> `Niche` (scanpy).
- `engine/outcome.py`: niche abundance vs survival -> `Prognostic` (lifelines).
- `engine/validate.py`: shuffle negative control, stability (ARI), marker validation.
  These are REQUIRED; they are how we prove niches are real, not artifacts.

### Lane B: MCP server (`src/locale/mcp_server/`)
Wraps the engine and exposes tools to K Pro. Do NOT put analysis logic here.

- `server.py`: FastMCP remote/HTTP server; connect as a custom connector, mirroring
  Owkin's Pathology Explorer.
- `tools.py`: the tool functions, wired to `data/mock.h5ad` so they return real
  schema objects TODAY; each tries the engine and falls back to a mock read until
  Lane A lands, then upgrades automatically.
- `interpret.py`: niche auto-naming via the Anthropic API (server-side sampling),
  with a deterministic keyless fallback.

Tools: `list_samples`, `describe_sample`, `compute_enrichment`, `find_niches`,
`characterize_niche`, `find_prognostic_niches` (orchestrator), `get_map_payload`.

### Lane C: visualization (`src/locale/viz/`)
The interactive tissue map, the demo's visual centerpiece.

- `payload.py`: AnnData + `image_id` -> `MapPayload`; regenerates the sample the
  widget renders from the mock.
- `app/index.html`: self-contained WebGL scatter (deck.gl, with a canvas-2D
  fallback) colored by `color_mode`, with a cell_type/niche toggle.

## Layout

```
src/locale/schema.py        shared data contract (the coordination artifact)
src/locale/data/            build the canonical AnnData (Lane A)
src/locale/engine/          spatial analysis engine (Lane A)
src/locale/mcp_server/      MCP tools (Lane B)
src/locale/viz/             tissue-map payload + widget (Lane C)
scripts/make_mock.py        generate committed data/mock.h5ad
scripts/download_data.py    remotezip fetch of only the needed files (run once)
tests/                      schema round-trip + mock schema + tool smoke tests
```

## Conventions

Python 3.11, type hints everywhere, docstrings that name the underlying tool, black
formatting, no em-dashes. Engine functions take the canonical AnnData and return
schema objects. Keep the schema stable.
