# Handoff

Snapshot of where Locale stands so anyone can keep working. `main` is green and
works end to end on the committed mock. Just pull and go.

## Status

- `main` is green: 44 tests passing.
- Lane A (analysis engine), Lane B (MCP server), and Lane C (viz + dashboard) have
  all landed. v1 runs end to end on `data/mock.h5ad`.
- The MCP server serves all seven tools, backed by the real engine when its deps
  are installed, and falls back gracefully per tool otherwise.

## Get running from scratch

```bash
cd ~/Owkin-Hackathon || (cd ~ && git clone https://github.com/sahielbose/Owkin-Hackathon.git && cd Owkin-Hackathon)
git pull
[ -d .venv ] || python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install pydantic anndata mcp anthropic squidpy scanpy scikit-learn lifelines pytest
python scripts/make_mock.py
pytest -q          # expect 44 passing
```

Install `squidpy scanpy scikit-learn`, NOT `cellcharter` from `requirements.txt`.
cellcharter is heavy, often fails to build, and the engine does not use it (niche
detection is k-means on neighborhood composition).

## Run the pieces

```bash
python -m src.localespatial.mcp_server.server     # MCP server on http://127.0.0.1:8000/mcp
python -m src.localespatial.viz.payload && python scripts/build_dashboard.py   # refresh visuals
open src/localespatial/viz/app/dashboard.html
```

To run on the REAL cohort instead of the mock:

```bash
LOCALE_DATA=data/locale.h5ad python -m src.localespatial.mcp_server.server
```

To demo in Claude / K Pro as a remote custom connector:

```bash
brew install cloudflared
cloudflared tunnel --url http://127.0.0.1:8000     # gives an https URL; connector URL is that + /mcp
```

Then Claude Settings -> Connectors -> Add custom connector -> name `Locale`, paste
`https://<tunnel>/mcp`. Full details in `src/localespatial/mcp_server/README.md`.

## What is left

- Lane A: build the real `data/locale.h5ad`. Run `scripts/download_data.py`
  (confirm the internal filenames against the printed listing), fill the
  column-name TODOs in `src/localespatial/data/build_anndata.py`, then run the validation
  checks (shuffle control, ARI stability, marker check) and lock the demo
  cohort/niche. Share `locale.h5ad` and the extracted CSVs on the Drive, never git.
- Lane B: optional bearer-token auth for the remote connector; set a real
  `ANTHROPIC_API_KEY` to turn on LLM niche naming (works without it via a
  deterministic fallback).
- Lane C: polish the dashboard and hero, wire the map to render live from
  `get_map_payload`, capture demo screenshots.
- Together: the killer demo (find_prognostic_niches -> survival curve -> tissue
  map) and a backup recording.

## Gotchas

- Do NOT change any field in `src/localespatial/schema.py` without telling the whole team.
  It is the shared contract all three lanes depend on. See the proposal below.
- Data is not in git. Only `data/mock.h5ad` is committed. `data/locale.h5ad` and
  raw CSVs are shared out of band on the Drive.
- The cloudflared quick-tunnel URL is random and changes on every restart, so
  re-paste it into Claude if you restart the tunnel.

## Proposed schema change (needs team sign-off, DO NOT apply unilaterally)

`src/localespatial/mcp_server/interpret.py` already produces both a niche `name` AND a
one-line `description`, but the `Niche` schema has no field to carry the
description, so it is currently computed and then dropped before reaching K Pro.

Proposal: add one optional, backwards-compatible field to `Niche` in
`src/localespatial/schema.py`:

```python
class Niche(BaseModel):
    niche_id: int
    name: str
    description: str | None = None   # <-- proposed: one-line human summary from interpret.py
    composition: dict[str, float]
    marker_program: list[str]
    prognostic: Prognostic | None = None
```

It defaults to `None`, so nothing breaks for existing callers. If approved, Lane B
will populate it from `interpret.name_and_describe_niche(...)`. Because the schema
is the shared contract, this needs a quick thumbs-up from all three lanes before
anyone edits `schema.py`.
