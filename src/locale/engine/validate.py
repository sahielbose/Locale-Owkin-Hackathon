"""Niche validation (Lane A). REQUIRED, not optional.

These three checks are how we prove the niches are real biology and not clustering
artifacts. A niche that fails them must not appear in the demo.
"""

from __future__ import annotations

from anndata import AnnData


def shuffle_negative_control(
    adata: AnnData, n_permutations: int = 100
) -> dict[str, float]:
    """Negative control: destroy spatial structure, confirm the signal collapses.

    Permute cell-type labels WITHIN each image (breaking spatial arrangement while
    preserving composition), re-run enrichment / niche detection, and confirm the
    real effect sizes are far outside the permuted null. If the "niche" survives
    shuffling, it is an artifact.

    Returns:
        Summary stats, e.g. {"real_effect": ..., "null_mean": ..., "null_std": ...,
        "empirical_p": ...}.

    TODO(Lane A):
        1. for each permutation, shuffle obs['cell_type'] within each image_id.
        2. recompute the statistic under test (enrichment z / niche separation).
        3. compare the real statistic to the null distribution; return an empirical p.
    """
    raise NotImplementedError(
        "shuffle_negative_control: within-image label permutation null."
    )


def stability_ari(
    adata: AnnData, n_runs: int = 10, subsample_frac: float = 0.8
) -> float:
    """Stability: are the niches reproducible across reruns?

    Re-run niche detection over different seeds / subsamples and measure the
    adjusted Rand index (ARI) between assignments on shared cells. High mean ARI
    means the niches are stable, not seed-dependent noise.

    Returns:
        Mean pairwise ARI across runs (1.0 = identical, ~0 = random).

    TODO(Lane A):
        1. run engine.niches.find_niches n_runs times (vary seed / subsample_frac).
        2. sklearn.metrics.adjusted_rand_score between every pair on shared cells.
        3. return the mean.
    """
    raise NotImplementedError(
        "stability_ari: rerun find_niches, mean pairwise adjusted_rand_score."
    )


def marker_validation(adata: AnnData, niche_id: int) -> dict[str, object]:
    """Biological sanity check: does the niche's marker program match known biology?

    Example expectations: an immune-excluded tumor core should be high in tumor /
    stromal markers (panCK, SMA) and low in CD8/GranzymeB; a TLS-like aggregate
    should be high in CD20/CD3. Flags niches whose programs are biologically
    incoherent.

    Returns:
        {"niche_id": niche_id, "passed": bool, "evidence": {...}}.

    TODO(Lane A):
        1. get the niche marker program (engine.characterize.characterize_niche).
        2. check expected markers are enriched / depleted as biology predicts.
        3. return pass/fail with the supporting marker z-scores.
    """
    raise NotImplementedError(
        "marker_validation: check niche marker program against known TME biology."
    )
