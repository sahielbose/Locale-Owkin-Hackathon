"""Fetch ONLY the files Locale needs from the 36.8 GB Zenodo archive, using
HTTP range requests (remotezip). We never download the whole zip.

Source: Zenodo record 3518284, file SingleCell_and_Metadata.zip
  Jackson & Fischer et al. 2020, "The single-cell pathology landscape of breast
  cancer," Nature.

What this does:
  (a) prints the archive's full file listing via RemoteZip.infolist(), so you can
      CONFIRM the exact internal paths of the files we want, then
  (b) extracts only the needed files (single-cell table, cell-type/cluster
      labels, patient metadata) into data/raw/.

Idempotent: a target that already exists in data/raw/ is skipped.

IMPORTANT: run this once, by hand. It is NOT run by CI or by any test. The
extracted CSVs and the resulting data/locale.h5ad are shared out of band (Drive),
never committed to git.

    python scripts/download_data.py            # list + extract needed files
    python scripts/download_data.py --list     # only print the archive listing
"""

from __future__ import annotations

import argparse
from pathlib import Path

ZENODO_URL = (
    "https://zenodo.org/records/3518284/files/" "SingleCell_and_Metadata.zip?download=1"
)

RAW_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"

# --- TODO: CONFIRM these against the printed infolist() output before trusting -----
# The exact internal paths inside SingleCell_and_Metadata.zip must be verified by
# running this script with --list first. The names below are best guesses based on
# the published archive layout (Data_publication/BaselTMA/...). Update them to the
# real paths, then re-run to extract. Each value is the member path INSIDE the zip.
NEEDED_FILES: dict[str, str] = {
    # single-cell marker-intensity table (one row per cell x 35 markers + ids)
    "single_cell": "Data_publication/BaselTMA/SC_dat.csv",  # TODO confirm
    # PhenoGraph cell-type / metacluster labels per cell
    "cell_types": "Cluster_labels/Basel_metaclusters.csv",  # TODO confirm
    # patient-level metadata (survival, grade, ER/PR/HER2, subtype)
    "patient_meta": "Data_publication/BaselTMA/Basel_PatientMetadata.csv",  # TODO confirm
}
# Note: the Zurich cohort has parallel Zurich_* files; add them here if you also
# want the validation cohort. Keep Basel as the primary demo cohort.
# ---------------------------------------------------------------------------------


def print_listing(rz) -> list[str]:
    """Print every member of an open RemoteZip and return their names.

    Reads only the zip central directory over HTTP range requests, so this is
    fast and does not download the 36.8 GB payload.
    """
    names: list[str] = []
    print(f"Listing {ZENODO_URL}\n(reading central directory only, not the payload)\n")
    for info in rz.infolist():
        size_mb = info.file_size / 1e6
        print(f"  {size_mb:10.2f} MB  {info.filename}")
        names.append(info.filename)
    print(f"\n{len(names)} members total.")
    return names


def extract_needed(rz, available: set[str]) -> None:
    """Extract only NEEDED_FILES into data/raw/, skipping any already present."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    missing = {
        key: member
        for key, member in NEEDED_FILES.items()
        if not (RAW_DIR / Path(member).name).exists()
    }
    if not missing:
        print("All needed files already present in data/raw/. Nothing to do.")
        return

    for key, member in missing.items():
        if member not in available:
            print(
                f"[SKIP] {key}: '{member}' not found in the archive. "
                f"Confirm the path against --list and update NEEDED_FILES."
            )
            continue
        print(f"[GET ] {key}: {member}")
        # Read the member bytes and write them straight to a flat data/raw/<basename>
        # so downstream code has stable paths and no empty nested dirs are left behind.
        (RAW_DIR / Path(member).name).write_bytes(rz.read(member))
    print(f"\nDone. Extracted files are in {RAW_DIR}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--list",
        action="store_true",
        help="only print the archive listing (use this to confirm NEEDED_FILES)",
    )
    args = parser.parse_args()

    from remotezip import RemoteZip

    # One RemoteZip context: read the central directory once, list, then extract.
    with RemoteZip(ZENODO_URL) as rz:
        names = print_listing(rz)
        if args.list:
            return
        print()
        extract_needed(rz, set(names))


if __name__ == "__main__":
    main()
