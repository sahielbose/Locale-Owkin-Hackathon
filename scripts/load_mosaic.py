"""Load MOSAIC Window spatial data into Locale's schema, then (optionally) analyze it.

MOSAIC Window is the 60-patient public tier of Owkin's MOSAIC spatial-omics dataset
(Glioblastoma, Bladder, Ovarian, DLBCL, Mesothelioma), delivered through the EGA. Its
spatial-transcriptomics modality is exactly what Locale reasons over: cells with a type,
an (x, y) position, and a sample/region id, optionally with patient survival.

This turns a MOSAIC Window export into the one object every Locale entrypoint consumes
(obs['cell_type'], obs['image_id'], obsm['spatial'], optional os_month/os_event/patient_id):

    # 1) write a canonical .h5ad Locale can run on
    python scripts/load_mosaic.py path/to/mosaic_window_spatial.h5ad -o data/mosaic.h5ad

    # from a cell table (csv/tsv/parquet) plus a clinical table for survival
    python scripts/load_mosaic.py cells.parquet --clinical clinical.csv -o data/mosaic.h5ad

    # 2) run the whole engine on it right now and print the headline result
    python scripts/load_mosaic.py path/to/mosaic_window_spatial.h5ad --analyze

    # 3) serve it as the live report / hand it to K Pro via the MCP server
    LOCALE_DATA=data/mosaic.h5ad python -m localespatial.mcp_server.server

Column names are auto-detected across the layouts spatial-transcriptomics exports use
(Visium/Xenium/CosMx/AnnData); pass --cell-type-col / --x-col / --y-col / --image-col /
--patient-col to override when a MOSAIC release names them differently. Nothing about the
science is MOSAIC-specific: the same path accepts any spatial single-cell object.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Candidate column names, most specific first. Auto-detection walks these in order.
CELL_TYPE = [
    "cell_type",
    "cell type",
    "celltype",
    "cell_types",
    "cell_phenotype",
    "phenotype",
    "cellType",
    "majorType",
    "major_type",
    "cell_annotation",
    "annotation",
    "ct",
    "label",
]
X_COLS = [
    "x",
    "X",
    "x_centroid",
    "Pos_X",
    "CenterX",
    "CentroidX",
    "spatial_x",
    "imagecol",
    "array_col",
    "x_um",
]
Y_COLS = [
    "y",
    "Y",
    "y_centroid",
    "Pos_Y",
    "CenterY",
    "CentroidY",
    "spatial_y",
    "imagerow",
    "array_row",
    "y_um",
]
IMAGE = [
    "image_id",
    "ImageId",
    "sample_id",
    "sample",
    "region_id",
    "region",
    "roi",
    "fov",
    "FOV",
    "core",
    "slide",
    "slide_id",
    "block",
    "spot",
]
PATIENT = [
    "patient_id",
    "patient",
    "case_id",
    "submitter_id",
    "donor_id",
    "subject_id",
    "PatientID",
]
OS_MONTH = [
    "os_month",
    "os_months",
    "OS_MONTHS",
    "overall_survival_months",
    "os_time",
    "survival_months",
    "survival_time",
    "followup_months",
    "time",
]
OS_EVENT = [
    "os_event",
    "OS_STATUS",
    "event",
    "death",
    "deceased",
    "vital_status",
    "censor",
    "status",
]


def _first(cols, candidates):
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand in cols:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def _coerce_event(series: pd.Series) -> np.ndarray:
    """Map an event/vital-status column to 1 = event (death), 0 = censored."""
    s = series
    if s.dtype.kind in "biufc":
        return (s.to_numpy(dtype=float) > 0).astype(int)
    txt = s.astype(str).str.strip().str.lower()
    dead = {
        "1",
        "1:deceased",
        "deceased",
        "dead",
        "death",
        "yes",
        "true",
        "event",
        "progressed",
    }
    return txt.isin(dead).to_numpy().astype(int)


def _read_table(path: Path) -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf in (".parquet", ".pq"):
        return pd.read_parquet(path)
    if suf in (".csv", ".gz"):
        return pd.read_csv(path)
    if suf in (".tsv", ".txt"):
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"unsupported table format: {path.name}")


def _adata_from_table(df: pd.DataFrame, args) -> ad.AnnData:
    cols = list(df.columns)
    ct = args.cell_type_col or _first(cols, CELL_TYPE)
    xc = args.x_col or _first(cols, X_COLS)
    yc = args.y_col or _first(cols, Y_COLS)
    img = args.image_col or _first(cols, IMAGE)
    if ct is None:
        raise ValueError(
            f"could not find a cell-type column; pass --cell-type-col. Saw: {cols[:20]}"
        )
    if xc is None or yc is None:
        raise ValueError(
            f"could not find x/y coordinate columns; pass --x-col/--y-col. Saw: {cols[:20]}"
        )

    reserved = {ct, xc, yc, img, args.patient_col, args.os_month_col, args.os_event_col}
    reserved |= set(PATIENT + OS_MONTH + OS_EVENT)
    markers = [
        c for c in cols if c not in reserved and pd.api.types.is_numeric_dtype(df[c])
    ]
    markers = [c for c in markers if c not in (xc, yc)]
    if not markers:
        # No expression columns (a pure phenotype table): give the engine a 1-D placeholder.
        markers = ["_placeholder"]
        X = np.ones((len(df), 1), dtype=float)
    else:
        X = df[markers].to_numpy(dtype=float)

    obs = pd.DataFrame(index=[str(i) for i in range(len(df))])
    # positional (.to_numpy) assignment: df's index need not match obs's string index
    obs["cell_type"] = pd.Categorical(df[ct].astype(str).to_numpy())
    img_vals = df[img].astype(str).to_numpy() if img else np.repeat("mosaic_1", len(df))
    obs["image_id"] = pd.Categorical(img_vals)
    adata = ad.AnnData(X=X, obs=obs)
    adata.var_names = [str(m) for m in markers]
    adata.obsm["spatial"] = df[[xc, yc]].to_numpy(dtype=float)
    # carry candidate patient / survival columns so _attach_clinical can find them
    for want, cands, override in (
        ("patient_id", PATIENT, args.patient_col),
        ("os_month", OS_MONTH, args.os_month_col),
        ("os_event", OS_EVENT, args.os_event_col),
    ):
        src = override or _first(cols, cands)
        if src and src in df:
            adata.obs[want] = df[src].to_numpy()
    return adata


def _attach_clinical(adata: ad.AnnData, clinical_path: Path, args) -> None:
    """Join a per-patient clinical table (survival) onto the cells by patient/sample key."""
    cdf = _read_table(clinical_path)
    ccols = list(cdf.columns)
    key = args.patient_col or _first(ccols, PATIENT) or _first(ccols, IMAGE)
    if key is None:
        raise ValueError(
            f"clinical table needs a patient/sample key column. Saw: {ccols[:20]}"
        )
    osm = args.os_month_col or _first(ccols, OS_MONTH)
    ose = args.os_event_col or _first(ccols, OS_EVENT)
    if osm is None:
        raise ValueError(
            f"clinical table needs a survival-months column; pass --os-month-col. Saw: {ccols[:20]}"
        )

    cdf = cdf.copy()
    cdf["_key"] = cdf[key].astype(str)
    lut_m = dict(zip(cdf["_key"], pd.to_numeric(cdf[osm], errors="coerce")))
    join_on = "patient_id" if "patient_id" in adata.obs else "image_id"
    keys = adata.obs[join_on].astype(str)
    adata.obs["os_month"] = keys.map(lut_m).to_numpy(dtype=float)
    if ose is not None:
        lut_e = dict(zip(cdf["_key"], _coerce_event(cdf[ose])))
        adata.obs["os_event"] = keys.map(lut_e).fillna(0).to_numpy().astype(int)
    if "patient_id" not in adata.obs:
        adata.obs["patient_id"] = keys.astype("category")


def load(args) -> ad.AnnData:
    src = Path(args.input)
    if not src.exists():
        raise SystemExit(f"input not found: {src}")
    if src.suffix.lower() in (".h5ad", ".h5"):
        adata = ad.read_h5ad(src)
    else:
        adata = _adata_from_table(_read_table(src), args)

    if args.clinical:
        _attach_clinical(adata, Path(args.clinical), args)

    from src.localespatial.webanalyze import canonicalize

    adata = canonicalize(adata)
    return adata


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Load MOSAIC Window spatial data into Locale."
    )
    ap.add_argument(
        "input",
        help="MOSAIC Window spatial export: .h5ad/.h5, or a cell table (.csv/.tsv/.parquet)",
    )
    ap.add_argument(
        "-o",
        "--out",
        default=str(ROOT / "data" / "mosaic.h5ad"),
        help="canonical .h5ad to write",
    )
    ap.add_argument(
        "--clinical",
        help="per-patient clinical table (.csv/.tsv/.parquet) with survival, joined by patient/sample key",
    )
    ap.add_argument(
        "--analyze",
        action="store_true",
        help="run the full Locale engine now and print the headline result",
    )
    ap.add_argument("--n-niches", type=int, default=6)
    ap.add_argument("--cell-type-col")
    ap.add_argument("--x-col")
    ap.add_argument("--y-col")
    ap.add_argument("--image-col")
    ap.add_argument("--patient-col")
    ap.add_argument("--os-month-col")
    ap.add_argument("--os-event-col")
    args = ap.parse_args()

    adata = load(args)
    has_surv = "os_month" in adata.obs and bool(
        np.isfinite(adata.obs["os_month"].to_numpy(dtype=float)).any()
    )
    print(
        f"loaded {adata.n_obs:,} cells x {adata.n_vars} markers, "
        f"{adata.obs['image_id'].nunique()} images, "
        f"{adata.obs['cell_type'].nunique()} cell types, survival={'yes' if has_surv else 'no'}"
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"wrote {out}")

    if args.analyze:
        from src.localespatial.webanalyze import analyze

        bundle = analyze(adata, n_niches=args.n_niches)
        niches = bundle["findings"]["niches"]
        surv = bundle["engine"]["survival"]
        risk = bundle["risk"]
        print("\n== Locale on MOSAIC Window ==")
        print(f"niches detected: {len(niches)}")
        top = sorted(surv, key=lambda s: -abs(np.log(max(s["hazard_ratio"], 1e-6))))[:3]
        for s in top:
            print(
                f"  {s['name']}: HR={s['hazard_ratio']:.2f} [{s['ci_low']:.2f},{s['ci_high']:.2f}] p={s['pvalue']:.3g}"
            )
        print(
            f"risk model verdict: {risk.get('verdict')} ({risk.get('verdict_reason', '')})"
        )

    print("\nrun it:")
    print(
        f"  LOCALE_DATA={out} python -m localespatial.mcp_server.server   # hand to K Pro"
    )
    print(
        f"  LOCALE_DATA={out} python scripts/serve.py                     # live report at :8000"
    )


if __name__ == "__main__":
    main()
