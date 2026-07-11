"""Fetch the Lane C 3D hero model from Sketchfab (Lane C).

We do NOT generate meshes. This pulls a real, externally-authored, CC-BY 4.0 cell
model into src/locale/viz/app/assets/ and records the required attribution. Every
model id below was verified against the Sketchfab v3 API (license.slug == "by",
isDownloadable == true). See src/locale/viz/app/assets/ASSETS.md for the full
sourcing manifest.

Sketchfab gates downloads behind a free API token (no cost, 2 clicks):
    Settings -> Password & API -> API token.

The Sketchfab download API hands back the glTF as a .zip (scene.gltf + scene.bin +
textures + license.txt). We extract it into assets/model/ and the viewer loads
assets/model/scene.gltf. The credit line is read from the archive's own license.txt.

Usage:
    export SKETCHFAB_TOKEN=...
    python scripts/fetch_assets.py                     # default hero -> assets/model/
    python scripts/fetch_assets.py --alt lymphocyte    # an alternate (see CATALOG)
    python scripts/fetch_assets.py --list              # list verified models
"""

from __future__ import annotations

import argparse
import io
import json
import os
import re
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path

APP_ASSETS = (
    Path(__file__).resolve().parents[1] / "src" / "locale" / "viz" / "app" / "assets"
)

# Verified 2026-07-11 via https://api.sketchfab.com/v3/models/<uid>
# (license.slug == "by", isDownloadable == true). key -> metadata.
CATALOG: dict[str, dict[str, str]] = {
    "prostate_cancer_cells": {
        "uid": "79074dd4d9c64bd9af4c0e34eff4d2b8",
        "name": "Prostate Cancer Cells",
        "author": "GLS",
        "faces": "1255997",
    },
    "tumor_vasculature": {
        "uid": "b6d5c662cf194091a152012deca2a932",
        "name": "Tumor Vasculature And Glycolysis",
        "author": "m-product",
        "faces": "53785",
    },
    "animal_cell": {
        "uid": "737b35f5b779418998d834c28ed15295",
        "name": "Animal Cell",
        "author": "James_Anthony",
        "faces": "280912",
    },
    "macrophage_kurzgesagt": {
        "uid": "0a62ebcb78ef484a92ea12463eb55093",
        "name": "Macrophage (from Kurzgesagt)",
        "author": "Spikefilmer",
        "faces": "92452",
    },
    "macrophage": {
        "uid": "c2af04b09e164e2a9e42ed321161ac43",
        "name": "Macrophage",
        "author": "Fezy",
        "faces": "184598",
    },
    "lymphocyte": {
        "uid": "5736f63ebbd54889b427af1c1dc3778e",
        "name": "Lymphocyte",
        "author": "3dcellstudio",
        "faces": "32416",
    },
    "components_of_blood": {
        "uid": "3ae309d331a049918b5788718ee58f35",
        "name": "Components of blood",
        "author": "arloopa",
        "faces": "78352",
    },
}
DEFAULT = "prostate_cancer_cells"
API = "https://api.sketchfab.com/v3/models"


def _get(url: str, token: str | None = None) -> bytes:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    if token:
        req.add_header("Authorization", f"Token {token}")
    with urllib.request.urlopen(req, timeout=120) as resp:
        return resp.read()


def _download_url(uid: str, token: str) -> str:
    """Resolve the short-lived signed glTF URL for a model via the download API."""
    payload = json.loads(_get(f"{API}/{uid}/download", token))
    gltf = payload.get("gltf") or {}
    url = gltf.get("url")
    if not url:
        raise RuntimeError(
            f"No glTF download URL returned for {uid}. Response keys: {sorted(payload)}"
        )
    return url


def _list() -> None:
    print("Verified CC-BY 4.0, downloadable (key: name -- author, faces):")
    for key, m in CATALOG.items():
        tag = "  [default]" if key == DEFAULT else ""
        print(f"  {key:22} {m['name']} -- {m['author']}, {m['faces']} faces{tag}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alt", default=DEFAULT, help="catalog key (see --list)")
    ap.add_argument(
        "--uid",
        help="fetch any Sketchfab model by uid (overrides --alt). Works for your own "
        "private/unlisted models when SKETCHFAB_TOKEN is yours. Verify the license "
        "before reusing a third-party model.",
    )
    ap.add_argument("--name", help="model name for attribution (with --uid)")
    ap.add_argument("--author", help="author for attribution (with --uid)")
    ap.add_argument(
        "--out",
        default=str(APP_ASSETS / "model"),
        help="output directory for the extracted glTF",
    )
    ap.add_argument("--list", action="store_true", help="list verified models and exit")
    args = ap.parse_args()

    if args.list:
        _list()
        return 0

    if args.uid:
        model = {
            "uid": args.uid,
            "name": args.name or "custom model",
            "author": args.author or "unknown",
            "faces": "?",
            "license": "unverified",
        }
    elif args.alt not in CATALOG:
        print(f"Unknown model {args.alt!r}. Options:", file=sys.stderr)
        _list()
        return 2
    else:
        model = CATALOG[args.alt]

    token = os.environ.get("SKETCHFAB_TOKEN")
    if not token:
        print(
            "SKETCHFAB_TOKEN not set.\n"
            "Get a free token at Sketchfab -> Settings -> Password & API -> API token,\n"
            "then:  export SKETCHFAB_TOKEN=...  and re-run.\n"
            f"Or download manually from https://sketchfab.com/3d-models/{model['uid']}\n"
            f"(Download 3D Model -> glTF) and unzip it into {args.out}/.",
            file=sys.stderr,
        )
        return 1

    out_dir = Path(args.out)
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Resolving download for {model['name']} ({model['author']})...")
    url = _download_url(model["uid"], token)
    print("Downloading glTF archive...")
    with urllib.request.urlopen(url, timeout=300) as resp:
        blob = resp.read()
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        zf.extractall(out_dir)
    gltf = next(out_dir.rglob("*.gltf"), None) or next(out_dir.rglob("*.glb"), None)
    if gltf is None:
        print(f"No .gltf/.glb found in the archive for {model['uid']}", file=sys.stderr)
        return 3
    print(f"extracted {out_dir}/ ({len(blob) / 1e6:.1f} MB, primary: {gltf.name})")

    attribution = _attribution(model, out_dir)
    attr_path = out_dir.parent / "attribution.json"
    attr_path.write_text(json.dumps(attribution, indent=2))
    print(f"wrote {attr_path}")
    print("Open src/locale/viz/app/hero3d.html to view.")
    return 0


def _attribution(model: dict[str, str], out_dir: Path) -> dict[str, str]:
    """Prefer the credit line Sketchfab ships in the archive's license.txt."""
    uid = model["uid"]
    lic_file = next(out_dir.rglob("license.txt"), None)
    if lic_file is not None:
        text = lic_file.read_text(errors="ignore")

        def field(label: str) -> str | None:
            m = re.search(rf"{label}:\s*(.+)", text)
            return m.group(1).strip() if m else None

        credit = None
        m = re.search(r"copy paste this credit[^\n]*\n(.+)", text, re.IGNORECASE)
        if m:
            credit = m.group(1).strip()
        return {
            "name": field("title") or model["name"],
            "author": re.sub(r"\s*\(http.*\)$", "", field("author") or model["author"]),
            "license": "CC-BY 4.0",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "source": field("source") or f"https://sketchfab.com/3d-models/{uid}",
            "credit": credit
            or f'"{model["name"]}" by {model["author"]}, licensed CC-BY 4.0 via Sketchfab.',
        }

    lic = model.get("license", "CC-BY 4.0")
    is_ccby = lic == "CC-BY 4.0"
    return {
        "name": model["name"],
        "author": model["author"],
        "license": lic,
        "license_url": "https://creativecommons.org/licenses/by/4.0/"
        if is_ccby
        else "",
        "source": f"https://sketchfab.com/3d-models/{uid}",
        "credit": (
            f'"{model["name"]}" by {model["author"]}, licensed CC-BY 4.0 via Sketchfab.'
            if is_ccby
            else f'"{model["name"]}" by {model["author"]} (verify license before reuse).'
        ),
    }


if __name__ == "__main__":
    raise SystemExit(main())
