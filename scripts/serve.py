"""Serve Locale as a live analysis engine.

Runs the report frontend AND an /api/analyze endpoint: upload a spatial single-cell
object (.h5ad) and the real engine returns the full analysis (niches, enrichment,
niche->survival, validation) that the report renders.

    python scripts/serve.py            # http://localhost:8000
    # then open the site, click "Analyze your data", drop an .h5ad

Requires the engine env plus fastapi + uvicorn (see requirements-web.txt).
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import anndata as ad
import uvicorn
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

ROOT = Path(__file__).resolve().parents[1]
APP_DIR = ROOT / "src" / "localespatial" / "viz" / "app"
sys.path.insert(0, str(ROOT))

api = FastAPI(title="Locale analysis engine")


@api.get("/api/health")
def health() -> dict:
    return {"status": "ok", "engine": "localespatial"}


@api.get("/api/example")
def example() -> JSONResponse:
    """Run the engine on the committed example object, so the analysis path is demoable."""
    from src.localespatial.webanalyze import analyze

    adata = ad.read_h5ad(ROOT / "data" / "mock.h5ad")
    bundle = analyze(adata)
    bundle["inputs"]["data_file"] = "data/mock.h5ad (example)"
    return JSONResponse(bundle)


@api.post("/api/analyze")
async def analyze_endpoint(file: UploadFile = File(...)) -> JSONResponse:
    name = (file.filename or "").lower()
    if not name.endswith((".h5ad", ".h5")):
        raise HTTPException(400, "upload an AnnData .h5ad file")
    raw = await file.read()
    if len(raw) > 200 * 1024 * 1024:
        raise HTTPException(413, "file too large (200 MB max)")
    try:
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=True) as tmp:
            tmp.write(raw)
            tmp.flush()
            adata = ad.read_h5ad(tmp.name)
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(400, f"could not read AnnData: {exc}")

    from src.localespatial.webanalyze import analyze

    try:
        bundle = analyze(adata)
    except ValueError as exc:
        raise HTTPException(422, str(exc))
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, f"analysis failed: {exc}")
    bundle["inputs"]["data_file"] = file.filename or "uploaded object"
    return JSONResponse(bundle)


# static site last so /api/* wins
api.mount("/", StaticFiles(directory=str(APP_DIR), html=True), name="site")


def main() -> None:
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    print(f"Locale analysis engine on http://localhost:{port}")
    uvicorn.run(api, host="0.0.0.0", port=port, log_level="warning")


if __name__ == "__main__":
    main()
