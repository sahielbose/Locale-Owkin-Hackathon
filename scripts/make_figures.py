"""Regenerate every figure in the README and the technical report from the cached
real-cohort object, so the plots in the repository are reproducible rather than
pasted images.

Reads (all gitignored, produced by the build pipeline):
    data/basel_niched.h5ad      cells, metacluster labels, cached spatial graph,
                                cached neighborhood-enrichment z, niche labels
    demo/findings.json          per-niche hazard ratio + CI + BH q (frozen demo path)
    data/neighborhood_heatmap.csv   the authors' published Delta statistic

Writes PNGs to docs/figures/. Deterministic (fixed seeds); touches no data on disk.

    python scripts/make_figures.py
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from src.localespatial.metaclusters import METACLUSTERS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGS = ROOT / "docs" / "figures"

MAJOR_ORDER = ["immune", "endothelial", "stroma", "tumor"]
MAJOR_COLOR = {
    "immune": "#2E86C1",
    "endothelial": "#8E44AD",
    "stroma": "#27AE60",
    "tumor": "#C0392B",
}

# Four cores that span the tumor-fraction range, from homogeneous (nests) to
# compartmentalised, used for the coordinate-validation hero figure.
HERO_CORES = [
    "BaselTMA_SP43_168_X10Y5",
    "BaselTMA_SP42_127_X2Y8",
    "BaselTMA_SP43_108_X13Y8",
    "BaselTMA_SP43_144_X15Y1",
]

plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    }
)


def _save(fig, name: str) -> None:
    out = FIGS / name
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out.relative_to(ROOT)}")


def fig_enrichment_heatmap(a: ad.AnnData) -> None:
    """Figure 1: neighborhood-enrichment z over the 27 metaclusters, ordered by major
    class. Block-diagonal (like-next-to-like) plus the tumor/immune off-diagonal."""
    z = np.asarray(a.uns["metacluster_id_nhood_enrichment"]["zscore"], dtype=float)
    cats = [int(c) for c in a.obs["metacluster_id"].cat.categories]
    order = sorted(
        range(len(cats)),
        key=lambda i: (MAJOR_ORDER.index(METACLUSTERS[cats[i]][1]), cats[i]),
    )
    zo = z[np.ix_(order, order)]
    labels = [METACLUSTERS[cats[i]][0] for i in order]
    majors = [METACLUSTERS[cats[i]][1] for i in order]

    fig, ax = plt.subplots(figsize=(11, 9.5))
    lim = float(np.nanpercentile(np.abs(zo), 98))
    im = ax.imshow(zo, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticklabels(labels, fontsize=7)
    # major-class group separators
    bounds = [i for i in range(1, len(majors)) if majors[i] != majors[i - 1]]
    for b in bounds:
        ax.axhline(b - 0.5, color="black", lw=1.0)
        ax.axvline(b - 0.5, color="black", lw=1.0)
    ax.set_title(
        "Neighborhood enrichment (permutation z) over 27 metaclusters\n"
        "block-diagonal is real tissue; tumor/immune off-diagonal is immune exclusion"
    )
    fig.colorbar(im, ax=ax, shrink=0.7, label="enrichment z")
    _save(fig, "enrichment_heatmap.png")


def fig_major_blocks(a: ad.AnnData) -> None:
    """The four-by-four major-class mean-z matrix. The headline: tumor vs immune z = -32."""
    z = np.asarray(a.uns["metacluster_id_nhood_enrichment"]["zscore"], dtype=float)
    cats = [int(c) for c in a.obs["metacluster_id"].cat.categories]
    majors = np.array([METACLUSTERS[c][1] for c in cats])
    M = np.zeros((4, 4))
    for i, r in enumerate(MAJOR_ORDER):
        for j, c in enumerate(MAJOR_ORDER):
            M[i, j] = float(np.nanmean(z[np.ix_(majors == r, majors == c)]))

    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    lim = float(np.abs(M).max())
    im = ax.imshow(M, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax.set_xticks(range(4))
    ax.set_yticks(range(4))
    ax.set_xticklabels(MAJOR_ORDER, rotation=30, ha="right")
    ax.set_yticklabels(MAJOR_ORDER)
    for i in range(4):
        for j in range(4):
            ax.text(
                j,
                i,
                f"{M[i, j]:+.1f}",
                ha="center",
                va="center",
                fontsize=12,
                color="white" if abs(M[i, j]) > lim * 0.5 else "black",
                fontweight="bold",
            )
    ti = M[MAJOR_ORDER.index("tumor"), MAJOR_ORDER.index("immune")]
    ax.set_title(
        "Major-class neighborhood enrichment (mean z)\n"
        f"tumor vs immune block mean = {ti:+.1f} (immune exclusion)"
    )
    fig.colorbar(im, ax=ax, shrink=0.8, label="mean enrichment z")
    _save(fig, "enrichment_major_blocks.png")


def fig_niche_composition(a: ad.AnnData) -> None:
    """Figure 4: major-class composition of each of the twelve discovered niches."""
    maj = a.obs["major"].astype(str).to_numpy()
    niche = a.obs["niche"].to_numpy().astype(int)
    ids = sorted(set(niche))
    comp = {
        m: [float((maj[niche == n] == m).mean()) for n in ids] for m in MAJOR_ORDER
    }
    fig, ax = plt.subplots(figsize=(11, 5))
    bottom = np.zeros(len(ids))
    for m in MAJOR_ORDER:
        vals = np.array(comp[m])
        ax.bar(
            [str(n) for n in ids],
            vals,
            bottom=bottom,
            color=MAJOR_COLOR[m],
            label=m,
            edgecolor="white",
            linewidth=0.5,
        )
        bottom += vals
    ax.set_ylim(0, 1)
    ax.set_xlabel("niche")
    ax.set_ylabel("fraction of cells")
    ax.set_title(
        "Major-class composition of the 12 niches "
        "(niche 1 immune-rich; niche 7 tumor-rich, immune-poor)"
    )
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.12), frameon=False)
    _save(fig, "niche_composition.png")


def fig_survival_forest() -> None:
    """Forest plot of the per-niche hazard ratios with 95% CIs from the frozen demo
    path. Nothing crosses into significance after correction; the two confirmatory
    niches (H2 niche 1, H1 niche 7) are marked."""
    f = json.loads((ROOT / "demo" / "findings.json").read_text())
    niches = {int(k): v for k, v in f["niches"].items()}
    ids = sorted(niches, key=lambda n: niches[n]["hazard_ratio"])
    fig, ax = plt.subplots(figsize=(8, 6))
    for row, n in enumerate(ids):
        hr = niches[n]["hazard_ratio"]
        lo, hi = niches[n]["ci_95"]
        conf = n in (1, 7)
        color = "#C0392B" if n == 7 else ("#2E86C1" if n == 1 else "#566573")
        ax.plot([lo, hi], [row, row], color=color, lw=2 if conf else 1.4)
        ax.plot(hr, row, "o", color=color, ms=9 if conf else 6)
        tag = {1: "  H2 (immune-rich)", 7: "  H1 (excluded)"}.get(n, "")
        ax.text(hi, row, f"  HR {hr:.2f}{tag}", va="center", fontsize=9, color=color)
    ax.axvline(1.0, color="black", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xticks([0.5, 0.7, 1.0, 1.4])
    ax.set_xticklabels(["0.5", "0.7", "1.0", "1.4"])
    ax.set_yticks(range(len(ids)))
    ax.set_yticklabels([f"niche {n}" for n in ids])
    ax.set_xlabel("hazard ratio per SD of niche abundance (log scale)")
    ax.set_title(
        "Per-niche survival association (Cox, adjusted)\n"
        "nothing survives BH correction; selection-aware empirical p = 0.44"
    )
    _save(fig, "survival_forest.png")


def fig_km_niche7(a: ad.AnnData) -> None:
    """Figure 6: overall survival by niche-7 abundance, median split. The unadjusted
    split inverts the sign of the adjusted Cox estimate; both are null."""
    obs = a.obs.copy()
    obs["PID"] = obs["PID"].astype(str)
    ab = pd.crosstab(obs["PID"], obs["niche"], normalize="index")
    pat = obs.drop_duplicates("PID").set_index("PID")[["OSmonth", "event"]]
    pat = pat.loc[ab.index]
    a7 = ab[7]
    hi = a7 > a7.median()
    lr = logrank_test(
        pat.OSmonth[hi], pat.OSmonth[~hi], pat.event[hi], pat.event[~hi]
    )
    fig, ax = plt.subplots(figsize=(7, 5))
    KaplanMeierFitter().fit(
        pat.OSmonth[hi], pat.event[hi], label=f"high niche-7 (n={int(hi.sum())})"
    ).plot_survival_function(ax=ax, color="#C0392B", ci_show=True)
    KaplanMeierFitter().fit(
        pat.OSmonth[~hi], pat.event[~hi], label=f"low niche-7 (n={int((~hi).sum())})"
    ).plot_survival_function(ax=ax, color="#2E86C1", ci_show=True)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("months")
    ax.set_ylabel("overall survival probability")
    ax.set_title(
        f"Overall survival by niche-7 abundance (median split)\n"
        f"log-rank p = {lr.p_value:.3f}; unadjusted sign is opposite the adjusted Cox"
    )
    _save(fig, "km_niche7.png")


def fig_cores_spatial(a: ad.AnnData) -> None:
    """Coordinate-validation hero: four Basel cores plotted at their recovered
    coordinates and coloured by major class, spanning tumor fraction from
    homogeneous (nests) to compartmentalised. A correct mask-to-cell join produces
    tissue morphology; a scrambled join would produce uniform confetti."""
    obs = a.obs
    xy = np.asarray(a.obsm["spatial"], dtype=float)
    core_arr = obs["core"].astype(str).to_numpy()
    maj_arr = obs["major"].astype(str).to_numpy()
    order = ["tumor", "stroma", "immune", "endothelial"]  # legend + z-order
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.6))
    for ax, core in zip(axes, HERO_CORES):
        m = core_arr == core
        pct = int(round((maj_arr[m] == "tumor").mean() * 100))
        for cls in order:
            mm = m & (maj_arr == cls)
            ax.scatter(xy[mm, 0], xy[mm, 1], s=2, c=MAJOR_COLOR[cls], label=cls)
        ax.set_title(f"{core}\n{int(m.sum()):,} cells, {pct}% tumor", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
    axes[0].legend(markerscale=4, loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(
        "Basel cores at recovered coordinates, coloured by major class. "
        "A correct join shows tissue nests and voids; a scrambled join would be uniform confetti.",
        y=1.03,
        fontsize=12,
    )
    _save(fig, "cores_spatial.png")


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    print("loading data/basel_niched.h5ad ...")
    a = ad.read_h5ad(DATA / "basel_niched.h5ad")
    fig_enrichment_heatmap(a)
    fig_major_blocks(a)
    fig_niche_composition(a)
    fig_survival_forest()
    fig_km_niche7(a)
    fig_cores_spatial(a)
    print("done.")


if __name__ == "__main__":
    main()
