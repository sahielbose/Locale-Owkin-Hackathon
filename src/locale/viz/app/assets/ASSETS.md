# Lane C 3D assets — sourcing manifest

The 3D hero viewer (`../hero3d.html`) renders a **real, externally-sourced** model.
Nothing here is procedurally generated. Every model below was verified against the
Sketchfab v3 API for its license (`license.slug`) and download flag
(`isDownloadable`) on 2026-07-11. All picks are **CC-BY 4.0** and downloadable, so
they may be redistributed **with attribution** (which the viewer renders on screen
and `fetch_assets.py` records to `attribution.json`).

Do not commit large binary `.glb` files casually. Run `scripts/fetch_assets.py`
(see below) to pull the chosen model into this folder as `cell.glb`.

## Recommended hero (default)

**Prostate Cancer Cells — @GLS**
- Why: professional medical-viz studio, ~1.26M faces, the most photoreal cancer-cell
  surface in the set. Cancer cells are the centerpiece of Locale's tumor niches.
- License: CC Attribution 4.0 (CC-BY). Downloadable: yes.
- UID: `79074dd4d9c64bd9af4c0e34eff4d2b8`
- Page: https://sketchfab.com/3d-models/79074dd4d9c64bd9af4c0e34eff4d2b8
- Attribution: "Prostate Cancer Cells" by GLS, licensed CC-BY 4.0 via Sketchfab.
- Note: 1.26M faces is heavy for low-end laptops. For the lighter option use the
  tumor-microenvironment model below.

## Verified alternates (all CC-BY 4.0, downloadable)

Ranked by fit to the breast-TME cell types in `CLAUDE.md`
(Tumor / CD8 T / CD4 T / Macrophage / Fibroblast / Endothelial):

| Model | Author | Faces | Fit | UID |
|---|---|---|---|---|
| Tumor Vasculature & Glycolysis | @m-product | 53,785 | tumor microenvironment, web-light | `b6d5c662cf194091a152012deca2a932` |
| Animal Cell (labeled organelles) | @James_Anthony | 280,912 | clean generic cell, education | `737b35f5b779418998d834c28ed15295` |
| Macrophage (from Kurzgesagt) | @Spikefilmer | 92,452 | TME immune cell | `0a62ebcb78ef484a92ea12463eb55093` |
| Macrophage | @Fezy | 184,598 | TME immune cell | `c2af04b09e164e2a9e42ed321161ac43` |
| Lymphocyte | @3dcellstudio | 32,416 | CD8/CD4 T cell, web-light | `5736f63ebbd54889b427af1c1dc3778e` |
| Components of blood | @arloopa | 78,352 | mixed immune/blood cells | `3ae309d331a049918b5788718ee58f35` |

## Embed-only (do NOT download — All Rights Reserved)

**Cancer cell under immune attack (membrane blebs) — @bblakesley** is the single most
on-message image for Locale (a tumor cell being attacked by immune cells, exactly the
immune-excluded vs immune-infiltrated niche story). But it is **not downloadable and
carries no reuse license**. Use it only as a Sketchfab `<iframe>` embed, never as a
committed asset.
- UID: `97de3b002c5f4fae95dfdf641be8b65a`
- Embed: https://sketchfab.com/models/97de3b002c5f4fae95dfdf641be8b65a/embed

## Other genuinely-good sources considered

- **NIH 3D (3d.nih.gov) / NIH BioArt** — public domain, but the 3D catalog is mostly
  untextured print meshes, not photoreal cells. Good for anatomy, weak for TME cells.
- **Poly Pizza** — CC-BY GLBs with direct download, but the API needs a free key and
  its biology catalog is thin.
- **TurboSquid / CGTrader / Free3D** — the highest-fidelity cancer-cell models, but
  paid and/or login-gated; licenses are per-item and often not redistributable.

## How to get the file into this folder

Downloads require a **free** Sketchfab token (2 clicks, no cost):
1. Sign in at sketchfab.com, open Settings → Password & API, copy your API token.
2. `export SKETCHFAB_TOKEN=...`
3. `python scripts/fetch_assets.py`            # pulls the default hero as assets/cell.glb
   `python scripts/fetch_assets.py --alt tumor_vasculature`   # or any alternate

`fetch_assets.py` writes `cell.glb` plus `attribution.json` (author, license, URL),
which the viewer reads to render the required CC-BY credit line.

To pull a specific model by id (for example your own private/unlisted upload, which
your token can access), bypass the catalog:

    python scripts/fetch_assets.py --uid <sketchfab_uid> --name "..." --author "..."

Verify the license yourself before reusing a third-party model this way; the catalog
entries above are the only ones we license-checked.

Manual fallback: on the model page click **Download 3D Model → glTF (.glb)**, drop the
file here as `cell.glb`, and copy its credit line into `attribution.json`.
