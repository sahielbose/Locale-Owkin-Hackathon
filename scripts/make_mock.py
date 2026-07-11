"""Generate data/mock.h5ad: a tiny synthetic AnnData that matches the canonical
Locale schema (see CLAUDE.md "The canonical AnnData object").

This file IS committed to git. It exists so Lanes B (MCP server) and C (viz) can
build and test against a realistic object before Lane A finishes the real
data/locale.h5ad. Keep it tiny and deterministic (fixed seed => byte-stable-ish
output => clean diffs).

The mock is deliberately biologically coherent so the demo reads well:
  - 4 IMC images across 3 patients (one patient has 2 cores)
  - 6 breast-TME cell types
  - 4 spatial niches with distinct compositions, including an "immune-excluded
    tumor core" that is enriched in the poor-survival patient
  - 35 fake protein markers whose intensities track cell type

Run:
    python scripts/make_mock.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SEED = 0
OUT_PATH = Path(__file__).resolve().parents[1] / "data" / "mock.h5ad"

# --- 35 fake but plausible breast-IMC protein markers ---------------------------
MARKERS: list[str] = [
    "HistoneH3",
    "panCK",
    "CK5",
    "CK7",
    "CK8_18",
    "CK14",
    "CK19",
    "ECadherin",
    "EpCAM",
    "Vimentin",
    "SMA",
    "Fibronectin",
    "Collagen1",
    "CD31",
    "vWF",
    "CD3",
    "CD4",
    "CD8a",
    "FOXP3",
    "CD20",
    "CD68",
    "CD163",
    "CD11c",
    "CD45",
    "CD45RO",
    "GranzymeB",
    "HLA_DR",
    "Ki67",
    "pS6",
    "cPARP",
    "p53",
    "ER",
    "PR",
    "HER2",
    "CAIX",
]
assert len(MARKERS) == 35, "canonical panel is 35 markers"

# --- 6 breast-TME cell types ----------------------------------------------------
CELL_TYPES: list[str] = [
    "Tumor",
    "CD8_T",
    "CD4_T",
    "Macrophage",
    "Fibroblast",
    "Endothelial",
]

# Markers each cell type expresses strongly (drives the synthetic intensities).
CELLTYPE_SIGNATURE: dict[str, list[str]] = {
    "Tumor": [
        "panCK",
        "CK5",
        "CK7",
        "CK8_18",
        "CK19",
        "ECadherin",
        "EpCAM",
        "Ki67",
        "p53",
        "ER",
        "HER2",
        "CAIX",
    ],
    "CD8_T": ["CD3", "CD8a", "CD45", "CD45RO", "GranzymeB", "HLA_DR"],
    "CD4_T": ["CD3", "CD4", "CD45", "CD45RO", "FOXP3", "HLA_DR"],
    "Macrophage": ["CD68", "CD163", "CD11c", "CD45", "HLA_DR"],
    "Fibroblast": ["SMA", "Vimentin", "Fibronectin", "Collagen1"],
    "Endothelial": ["CD31", "vWF", "Vimentin"],
}

# --- 4 spatial niches: composition over CELL_TYPES (in the order above) ---------
NICHE_NAMES: dict[int, str] = {
    0: "immune-excluded tumor core",
    1: "immune infiltrate",
    2: "vascular stroma",
    3: "proliferative tumor edge",
}
NICHE_COMPOSITION: dict[int, list[float]] = {
    # Tumor, CD8_T, CD4_T, Macrophage, Fibroblast, Endothelial
    0: [0.65, 0.02, 0.03, 0.05, 0.22, 0.03],
    1: [0.15, 0.30, 0.25, 0.20, 0.05, 0.05],
    2: [0.10, 0.05, 0.05, 0.10, 0.40, 0.30],
    3: [0.55, 0.08, 0.07, 0.20, 0.07, 0.03],
}
# Blob center for each niche within a 1000x1000 micron field.
NICHE_CENTER: dict[int, tuple[float, float]] = {
    0: (300.0, 300.0),
    1: (700.0, 300.0),
    2: (300.0, 700.0),
    3: (700.0, 700.0),
}
NICHE_SIGMA = 90.0  # microns

# --- images -> patient, and per-image niche cell counts (sum to 150 each) -------
# IMG001/IMG002 belong to the poor-survival patient (niche-0 heavy).
IMAGE_TO_PATIENT: dict[str, str] = {
    "IMG001": "P001",
    "IMG002": "P001",
    "IMG003": "P002",
    "IMG004": "P003",
}
IMAGE_NICHE_COUNTS: dict[str, list[int]] = {
    # niche 0, 1, 2, 3
    "IMG001": [70, 20, 30, 30],
    "IMG002": [65, 25, 30, 30],
    "IMG003": [25, 60, 35, 30],
    "IMG004": [45, 40, 35, 30],
}

# --- per-patient survival + clinical (immune-excluded tracks worse outcome) -----
PATIENT_CLINICAL: dict[str, dict] = {
    "P001": dict(
        os_month=18.0,
        os_event=1,
        dfs_month=9.0,
        dfs_event=1,
        grade=3,
        er="neg",
        pr="neg",
        her2="neg",
        subtype="TNBC",
    ),
    "P002": dict(
        os_month=60.0,
        os_event=0,
        dfs_month=60.0,
        dfs_event=0,
        grade=1,
        er="pos",
        pr="pos",
        her2="neg",
        subtype="LumA",
    ),
    "P003": dict(
        os_month=42.0,
        os_event=1,
        dfs_month=30.0,
        dfs_event=1,
        grade=2,
        er="pos",
        pr="neg",
        her2="pos",
        subtype="HER2",
    ),
}


def _make_raw_intensities(
    cell_types: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Synthesize non-negative marker intensities that track cell type.

    Baseline gamma noise for every marker, plus a boost on each cell's signature
    markers (and HistoneH3 for all cells, since it is a nuclear stain).
    """
    n = cell_types.shape[0]
    marker_idx = {m: i for i, m in enumerate(MARKERS)}
    raw = rng.gamma(shape=2.0, scale=1.0, size=(n, len(MARKERS))).astype(np.float64)
    raw[:, marker_idx["HistoneH3"]] += rng.gamma(shape=4.0, scale=2.0, size=n)
    for ct, sig in CELLTYPE_SIGNATURE.items():
        rows = np.where(cell_types == ct)[0]
        for m in sig:
            raw[rows, marker_idx[m]] += rng.gamma(
                shape=4.0, scale=2.0, size=rows.shape[0]
            )
    return raw


def build_mock() -> ad.AnnData:
    """Build the tiny canonical AnnData in memory."""
    rng = np.random.default_rng(SEED)

    cell_type_col: list[str] = []
    patient_col: list[str] = []
    image_col: list[str] = []
    niche_col: list[int] = []
    coords: list[tuple[float, float]] = []

    for image_id, counts in IMAGE_NICHE_COUNTS.items():
        patient_id = IMAGE_TO_PATIENT[image_id]
        # small per-image offset so images do not overlap exactly
        offset = rng.uniform(-40.0, 40.0, size=2)
        for niche_id, n_cells in enumerate(counts):
            probs = NICHE_COMPOSITION[niche_id]
            drawn = rng.choice(CELL_TYPES, size=n_cells, p=probs)
            cx, cy = NICHE_CENTER[niche_id]
            xy = rng.normal(
                loc=(cx + offset[0], cy + offset[1]),
                scale=NICHE_SIGMA,
                size=(n_cells, 2),
            )
            cell_type_col.extend(drawn.tolist())
            patient_col.extend([patient_id] * n_cells)
            image_col.extend([image_id] * n_cells)
            niche_col.extend([niche_id] * n_cells)
            coords.extend([tuple(p) for p in xy])

    cell_types = np.array(cell_type_col, dtype=object)
    n = cell_types.shape[0]

    # Same pipeline as build_anndata.preprocess_intensities with its defaults:
    # arcsinh(cofactor 5) -> clip at 99th percentile per marker -> z-score per marker.
    # Kept inline (not imported) so this script stays standalone.
    raw = _make_raw_intensities(cell_types, rng)
    transformed = np.arcsinh(raw / 5.0)
    hi = np.percentile(transformed, 99.0, axis=0, keepdims=True)
    transformed = np.minimum(transformed, hi)
    mean = transformed.mean(axis=0, keepdims=True)
    std = transformed.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    x = ((transformed - mean) / std).astype(np.float32)

    obs = pd.DataFrame(index=[f"cell_{i:04d}" for i in range(n)])
    obs["cell_type"] = pd.Categorical(cell_types, categories=CELL_TYPES)
    obs["patient_id"] = patient_col
    obs["image_id"] = image_col
    obs["niche"] = np.array(niche_col, dtype=np.int64)
    for field in (
        "os_month",
        "os_event",
        "dfs_month",
        "dfs_event",
        "grade",
        "er",
        "pr",
        "her2",
        "subtype",
    ):
        obs[field] = [PATIENT_CLINICAL[p][field] for p in patient_col]
    obs["os_month"] = obs["os_month"].astype(float)
    obs["os_event"] = obs["os_event"].astype(int)
    obs["dfs_month"] = obs["dfs_month"].astype(float)
    obs["dfs_event"] = obs["dfs_event"].astype(int)
    obs["grade"] = obs["grade"].astype(int)

    var = pd.DataFrame(index=pd.Index(MARKERS, name=None))

    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["spatial"] = np.asarray(coords, dtype=np.float32)
    adata.uns["markers"] = list(MARKERS)
    # h5ad mapping keys must be strings.
    adata.uns["niche_names"] = {str(k): v for k, v in NICHE_NAMES.items()}
    adata.uns["mock"] = True
    return adata


def main() -> None:
    adata = build_mock()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT_PATH)
    print(f"wrote {OUT_PATH}")
    print(adata)
    print("\ncell_type counts:\n", adata.obs["cell_type"].value_counts())
    print("\nniche counts:\n", adata.obs["niche"].value_counts().sort_index())
    print(
        "\nimages -> patients:\n",
        adata.obs.groupby("image_id", observed=True)["patient_id"].first(),
    )


if __name__ == "__main__":
    main()
