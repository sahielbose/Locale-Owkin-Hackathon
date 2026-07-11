# Locale MCP server (Lane B)

A remote MCP server that exposes Locale's spatial-niche tools to Claude / K Pro. It
wraps the Lane A engine and returns the shared schema objects. Every tool works on
`data/mock.h5ad` today and auto-upgrades to the real analysis as Lane A lands each
engine function.

## Tools

| Tool | Returns | Notes |
|---|---|---|
| `list_samples()` | `list[SampleRecord]` | one record per IMC image |
| `describe_sample(image_id?, cohort?)` | `SampleRecord` | image-level or cohort-level |
| `compute_enrichment(scope)` | `EnrichmentResult` | cell-type co-location matrix |
| `find_niches(cohort?, n_niches?)` | `list[Niche]` | elicits cohort if ambiguous |
| `characterize_niche(niche_id)` | `Niche` | composition + marker program + name |
| `find_prognostic_niches(cohort?, patient_subset?)` | `list[Niche]` | orchestrator, ranked worst-survival first |
| `get_map_payload(image_id, color_mode)` | `MapPayload` | tissue-map render payload |

## Graceful degradation (the core pattern)

Every tool tries the real Lane A engine function first; on `NotImplementedError`
(not built yet) or any exception (a mid-edit break, bad data) it falls back to a
value derived from the loaded AnnData, so a tool ALWAYS returns a valid schema
object. The path taken is logged and recorded in `tools.backend_status()`
(`"real"` vs `"fallback"`). Fallbacks:

- enrichment: a temporary within-image kNN adjacency z-score (until squidpy lands)
- niches: the precomputed `obs['niche']` in the object
- survival: a lifelines Cox + Kaplan-Meier on per-patient niche abundance
  (regularized and clamped on tiny cohorts; `prognostic=None` if not viable)
- naming: `interpret.py` uses the Anthropic API if `ANTHROPIC_API_KEY` is set, else a
  deterministic name

## Run it locally

```bash
source .venv/bin/activate
pip install -r requirements.txt          # or the light set: pydantic anndata mcp anthropic lifelines
python -m src.locale.mcp_server.server   # Streamable HTTP on http://127.0.0.1:8000/mcp
```

Environment:

| Var | Default | Meaning |
|---|---|---|
| `LOCALE_HOST` | `127.0.0.1` | bind host (`0.0.0.0` to expose) |
| `LOCALE_PORT` | `8000` | bind port |
| `LOCALE_TRANSPORT` | `streamable-http` | use `sse` only if a client needs it |
| `LOCALE_DATA` | auto | force a specific AnnData (else `data/locale.h5ad`, else `data/mock.h5ad`) |
| `ANTHROPIC_API_KEY` | unset | enables LLM niche naming (optional) |
| `LOCALE_LOG_LEVEL` | `INFO` | logging level |

## Verify with the MCP Inspector first

Before wiring it into Claude, confirm the tools with the Inspector:

```bash
# Option A: let the SDK launch the server object with the Inspector
mcp dev src/locale/mcp_server/server.py

# Option B: run the server, then point the standalone Inspector at the URL
python -m src.locale.mcp_server.server
npx @modelcontextprotocol/inspector        # transport: Streamable HTTP, URL: http://127.0.0.1:8000/mcp
```

Call `list_samples` and `find_prognostic_niches` and confirm you get structured results.

## Expose it as a REMOTE custom connector

Claude custom connectors need a public **https** URL. For the hackathon, tunnel the
local server:

```bash
# cloudflared (no account needed for a quick tunnel)
cloudflared tunnel --url http://127.0.0.1:8000
# -> prints https://<random>.trycloudflare.com  ; your MCP URL is that + /mcp

# or ngrok
ngrok http 8000
# -> https://<random>.ngrok-free.app            ; your MCP URL is that + /mcp
```

Then in **Claude → Settings → Connectors → Add custom connector**:

1. Name: `Locale`
2. Remote MCP server URL: `https://<your-tunnel>/mcp`
3. Save, then enable it in a chat and confirm the seven tools appear.

K Pro uses the same custom-connector pattern (the remote MCP URL), mirroring Owkin's
Pathology Explorer.

## Authentication

None for now. If Claude requires auth for a remote connector, add either:

- a static bearer token: wrap the ASGI app with middleware that checks
  `Authorization: Bearer <token>` (FastMCP exposes `mcp.streamable_http_app()`), or
- OAuth: configure an `OAuthAuthorizationServerProvider` / token verifier on the
  `FastMCP` instance (see the MCP Python SDK auth docs).

Add auth before exposing anything beyond a short-lived hackathon tunnel.

## Tests

```bash
pytest tests/test_mcp_server.py -q
```

Covers: every tool validates against its Pydantic schema, the orchestrator ranking,
the elicitation resolver, and `interpret.py` (fallback + a monkeypatched API path;
no real API calls are made).
