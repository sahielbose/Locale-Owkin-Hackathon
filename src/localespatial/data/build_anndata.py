"""Build the canonical Locale AnnData (data/locale.h5ad) from the CSVs that
scripts/download_data.py extracts into data/raw/.

This is Lane A's entry point. The preprocessing MATH is implemented and real
(arcsinh cofactor 5 -> optional 99th-percentile clip -> per-marker z-score). The
CSV parsing and column mapping are STUBBED: the real column names in the Jackson
2020 tables must be filled into the *_COL / MARKER_COLUMNS constants below, which
are the only TODOs. Everything downstream reads the object this produces.

Canonical object (see CLAUDE.md):
    X            float32 [n_cells x 35]   arcsinh(cofactor=5), z-scored per marker
    var_names    35 marker names
    obs          cell_type, patient_id, image_id, os_month, os_event,
                 dfs_month?, dfs_event?, grade?, er?, pr?, her2?, subtype?
    obsm         'spatial' float [n_cells x 2] (x, y microns)
    uns          'markers' list[str]

    python -m src.locale.data.build_anndata
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

RAW_DIR = Path(__file__).resolve().parents[3] / "data" / "raw"
OUT_PATH = Path(__file__).resolve().parents[3] / "data" / "locale.h5ad"

# File basenames as flattened by scripts/download_data.py (keep in sync).
SINGLE_CELL_CSV = "SC_dat.csv"
CELL_TYPES_CSV = "Basel_metaclusters.csv"
PATIENT_META_CSV = "Basel_PatientMetadata.csv"

ARCSINH_COFACTOR = 5.0
CLIP_PERCENTILE: float | None = 99.0  # set to None to disable the clip

# --- TODO: fill in the REAL column names from the extracted CSVs -------------------
# Inspect the head of each CSV (pandas .head()) and set these. They are the only
# thing standing between this stub and a real data/locale.h5ad.
#
# The single-cell table in this dataset is often LONG (one row per cell x marker);
# if so, pivot it to wide (cells x markers) in _load_single_cell(). If it is already
# wide, list the 35 marker columns in MARKER_COLUMNS.
CELL_ID_COL = "id"  # TODO: unique cell id, present in both SC + cell-type tables
IMAGE_ID_COL = "core"  # TODO: IMC core / image id column
X_COORD_COL = "Location_Center_X"  # TODO: x coordinate (microns)
Y_COORD_COL = "Location_Center_Y"  # TODO: y coordinate (microns)
MARKER_COLUMNS: list[str] | None = None  # TODO: the 35 marker columns (if wide)

CELLTYPE_LABEL_COL = "cluster"  # TODO: PhenoGraph phenotype / metacluster label
CELLTYPE_ID_COL = "id"  # TODO: cell id column in the cell-type table

PATIENT_ID_COL = "PID"  # TODO: patient id in metadata (and joinable to cells)
CORE_TO_PATIENT_COL = "core"  # TODO: core->patient mapping (often in patient meta)
OS_MONTH_COL = "OSmonth"  # TODO
OS_EVENT_COL = "Patientstatus"  # TODO: map to {0,1} in _map_survival()
DFS_MONTH_COL = "DFSmonth"  # TODO (optional)
DFS_EVENT_COL = "DFSevent"  # TODO (optional)
CLINICAL_COLS = {  # TODO: map canonical name -> real column name (optional block)
    "grade": "grade",
    "er": " ERStatus",
    "pr": "PRStatus",
    "her2": "HER2Status",
    "subtype": "clinical_type",
}
# ---------------------------------------------------------------------------------


def arcsinh_transform(x: np.ndarray, cofactor: float = ARCSINH_COFACTOR) -> np.ndarray:
    """arcsinh(x / cofactor). Standard IMC intensity transform (cofactor 5)."""
    return np.arcsinh(x.astype(np.float64) / cofactor)


def clip_percentile(x: np.ndarray, percentile: float = 99.0) -> np.ndarray:
    """Clip each marker (column) at its given upper percentile to tame outliers."""
    hi = np.percentile(x, percentile, axis=0, keepdims=True)
    return np.minimum(x, hi)


def zscore_per_marker(x: np.ndarray) -> np.ndarray:
    """Standardize each marker (column) to mean 0, unit variance."""
    mean = x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    return (x - mean) / std


def preprocess_intensities(raw: np.ndarray) -> np.ndarray:
    """Full marker pipeline: arcsinh(cofactor 5) -> optional 99th clip -> z-score.

    Returns float32 ready for adata.X.
    """
    x = arcsinh_transform(raw)
    if CLIP_PERCENTILE is not None:
        x = clip_percentile(x, CLIP_PERCENTILE)
    x = zscore_per_marker(x)
    return x.astype(np.float32)


def _load_single_cell() -> tuple[pd.DataFrame, list[str]]:
    """Read the single-cell marker table and return (wide cells x markers, marker_names).

    TODO(Lane A): the Jackson table may be long-format (id, channel, value). If so,
    pivot to wide here: df.pivot(index=CELL_ID_COL, columns=<channel>, values=<value>).
    If it is already wide, just select MARKER_COLUMNS.
    """
    path = RAW_DIR / SINGLE_CELL_CSV
    _ = pd.read_csv(path)  # noqa: F841  (kept so the stub fails loudly if file missing)
    raise NotImplementedError(
        "TODO(Lane A): reshape the single-cell CSV to wide (cells x 35 markers), "
        "keyed by CELL_ID_COL, and return (df, marker_names). Confirm MARKER_COLUMNS "
        "or the long->wide pivot against the real SC_dat.csv header."
    )


def _load_cell_types() -> pd.Series:
    """Read PhenoGraph labels -> Series indexed by cell id. TODO confirm columns."""
    path = RAW_DIR / CELL_TYPES_CSV
    _ = pd.read_csv(path)  # noqa: F841
    raise NotImplementedError(
        "TODO(Lane A): return a Series of cell_type labels indexed by CELL_ID_COL, "
        "using CELLTYPE_ID_COL and CELLTYPE_LABEL_COL. Optionally remap raw "
        "metacluster ids to readable breast-TME names."
    )


def _load_patient_meta() -> pd.DataFrame:
    """Read patient metadata (survival + clinical). TODO confirm columns + status map."""
    path = RAW_DIR / PATIENT_META_CSV
    _ = pd.read_csv(path)  # noqa: F841
    raise NotImplementedError(
        "TODO(Lane A): return a per-cell (or per-core) frame with os_month/os_event "
        "(+ optional dfs, grade, er, pr, her2, subtype). Map OS_EVENT_COL text "
        "(e.g. 'death'/'alive') to {1,0}."
    )


def build_anndata() -> ad.AnnData:
    """Assemble the canonical AnnData from the extracted CSVs.

    Steps (Lane A):
      1. load single-cell wide matrix + marker names
      2. preprocess_intensities() -> X
      3. join cell-type labels (per cell)
      4. join patient survival + clinical (per patient/core -> per cell)
      5. set obsm['spatial'] from X/Y coordinate columns
      6. set uns['markers']
    """
    cells, markers = _load_single_cell()
    cell_types = _load_cell_types()
    patient_meta = _load_patient_meta()

    raw = cells[markers].to_numpy()
    x = preprocess_intensities(raw)

    obs = pd.DataFrame(index=cells[CELL_ID_COL].astype(str).to_numpy())
    obs["cell_type"] = pd.Categorical(cell_types.reindex(obs.index).to_numpy())
    obs["image_id"] = cells[IMAGE_ID_COL].astype(str).to_numpy()
    # TODO(Lane A): broadcast patient-level survival/clinical onto cells via the
    # core/patient join in patient_meta, filling os_month/os_event/etc.
    _ = patient_meta  # noqa: F841

    var = pd.DataFrame(index=pd.Index(markers))
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["spatial"] = cells[[X_COORD_COL, Y_COORD_COL]].to_numpy(dtype=np.float32)
    adata.uns["markers"] = list(markers)
    return adata


def main() -> None:
    adata = build_anndata()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT_PATH)
    print(f"wrote {OUT_PATH}")
    print(adata)


if __name__ == "__main__":
    main()
