"""Build data/basel_nomarkers.h5ad from the extracted CSVs (real Basel cohort).

The unlock: only characterize_niche needs the 35 markers. Graph, enrichment, niches,
abundance, survival, and the ARI need only cell type + coordinates + survival, which
is ~50 MB. So X is a placeholder here; markers are backfilled from SC_dat.csv later.

Inputs (all in data/, produced earlier, gitignored):
  coords.csv                 x/y per cell, from scripts/extract_coords.py (masks)
  PG_final_k20.csv           PhenoGraph cell type per id
  Basel_PatientMetadata.csv  survival + clinical per core

Steps: join coords + cell type on `id`; join survival on `core`; drop non-tumor
cores; build the 0/1 event from the Patientstatus STRING with an explicit map;
add metaclusters (71 -> 27 + 4 major); build the per-core spatial graph.

    python scripts/build_basel.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq

from src.localespatial.metaclusters import annotate

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = DATA / "basel_nomarkers.h5ad"

# Explicit overall-survival event map from the Patientstatus string. Any death is an
# OS event; both alive categories are censored. Do NOT trust a raw string as 0/1.
DEATH = {"death by primary disease", "death"}
ALIVE = {"alive", "alive w metastases"}


def main() -> None:
    coords = pd.read_csv(DATA / "coords.csv")
    pg = pd.read_csv(DATA / "PG_final_k20.csv")
    meta = pd.read_csv(DATA / "Basel_PatientMetadata.csv")

    df = coords.merge(pg[["id", "PhenoGraphBasel"]], on="id", how="inner")
    df = df.merge(
        meta[
            [
                "core",
                "PID",
                "OSmonth",
                "Patientstatus",
                "grade",
                "clinical_type",
                "diseasestatus",
            ]
        ],
        on="core",
        how="left",
    )
    df = df[df["diseasestatus"] == "tumor"].copy()  # drop non-tumor cores

    df["event"] = df["Patientstatus"].map(
        lambda s: 1 if s in DEATH else (0 if s in ALIVE else np.nan)
    )
    if df["event"].isna().any():
        bad = sorted(df.loc[df["event"].isna(), "Patientstatus"].unique())
        raise ValueError(f"unmapped Patientstatus levels: {bad}")
    df["cell_type"] = "PG" + df["PhenoGraphBasel"].astype(int).astype(str)

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(df["cell_type"].values),
            "core": pd.Categorical(df["core"].astype(str).values),
            "PID": df["PID"].astype(str).values,
            "OSmonth": df["OSmonth"].astype(float).values,
            "event": df["event"].astype(int).values,
            "grade": df["grade"].astype(str).values,
            "clinical_type": df["clinical_type"].astype(str).values,
        },
        index=df["id"].values,
    )
    a = ad.AnnData(X=np.zeros((len(df), 1), dtype=np.float32), obs=obs)
    a.obsm["spatial"] = df[["x", "y"]].to_numpy(dtype=np.float32)
    a.uns["markers_placeholder"] = True

    # named metaclusters (27) + major class (4)
    a.obs["pg_cluster"] = a.obs["cell_type"].astype(str).str[2:].astype(int)
    annotate(a, pg_key="pg_cluster")

    # per-core spatial graph (image = core)
    sq.gr.spatial_neighbors(a, library_key="core", coord_type="generic", delaunay=True)

    a.write_h5ad(OUT)
    pat = a.obs.drop_duplicates("PID")
    print(f"wrote {OUT}")
    print(
        f"n_cells {a.n_obs:,} | n_patients {a.obs.PID.nunique()} | n_cores {a.obs.core.nunique()} "
        f"| deaths {int(pat.event.sum())}"
    )
    print("major class counts:", a.obs["major"].value_counts().to_dict())


if __name__ == "__main__":
    main()
