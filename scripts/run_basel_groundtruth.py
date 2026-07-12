"""Ground-truth check: our niches vs Jackson et al.'s published tumor communities (STEP 5).

The headline is the tumor-only Adjusted Rand Index. Two different methods (our k-means
on neighborhood composition vs their graph community detection), so a moderate ARI is
the real result; perfect agreement would be suspicious.

    python scripts/run_basel_groundtruth.py
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad

from src.locale.engine.groundtruth import community_ari

DATA = Path(__file__).resolve().parents[1] / "data"


def main() -> None:
    a = ad.read_h5ad(DATA / "basel_niched.h5ad")
    res = community_ari(
        a,
        str(DATA / "Communtiy_data_tumor.csv"),
        str(DATA / "PG_tumor_communities.csv"),
    )
    print("=== STEP 5: ARI vs published tumor communities ===")
    print(
        f"our niches: {res['n_our_niches']} | their community phenotypes: {res['n_their_communities']}"
    )
    print(f"joined tumor cells: {res['n_tumor']:,}")
    print(f"\nHEADLINE tumor-only ARI = {res['ari_tumor_only']:.3f}")
    print(f"(all matched: {res['ari_all_matched']:.3f})")


if __name__ == "__main__":
    main()
