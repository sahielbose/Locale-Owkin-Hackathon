"""Niche -> outcome association (Lane A). Backed by lifelines.

This is where the tool becomes a biomarker engine: does a niche's per-patient
abundance track overall survival? Cox proportional hazards + Kaplan-Meier, with
FDR correction across niches.
"""

from __future__ import annotations

from anndata import AnnData

from ..schema import Prognostic


def niche_outcome(
    adata: AnnData, niche_id: int, high_low_quantile: float = 0.5
) -> Prognostic:
    """Associate one niche's abundance with survival (lifelines Cox + KM).

    Args:
        adata: canonical AnnData with obs['niche'], obs['patient_id'],
            obs['os_month'], obs['os_event'].
        niche_id: the niche whose abundance is tested.
        high_low_quantile: split patients into high/low niche-abundance groups at
            this quantile for the KM curve (default median).

    Returns:
        Prognostic(hazard_ratio, ci_low, ci_high, pvalue, n_patients, km=KMCurve).

    TODO(Lane A):
        1. per patient: fraction of that patient's cells in niche_id (abundance).
        2. Cox: lifelines.CoxPHFitter on (abundance -> os_month, os_event); read HR
           and 95% CI. FDR-correct pvalues across all niches at the caller level.
        3. KM: split patients at high_low_quantile, KaplanMeierFitter per group,
           sample survival probs on a shared time grid -> KMCurve.
        4. return Prognostic.
    """
    raise NotImplementedError(
        "niche_outcome: wire lifelines Cox + KM on per-patient niche abundance."
    )


def rank_prognostic_niches(adata: AnnData) -> list[tuple[int, Prognostic]]:
    """Rank every niche by prognostic strength (FDR-corrected across niches).

    TODO(Lane A): compute niche_outcome for each niche, BH-correct the p-values,
    and return (niche_id, Prognostic) sorted by significance / |log HR|.
    """
    raise NotImplementedError(
        "rank_prognostic_niches: niche_outcome per niche + BH-FDR, ranked."
    )
