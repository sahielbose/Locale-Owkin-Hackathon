"""Pre-registered niche -> survival analysis on real Basel (STEP 6).

Hypotheses are frozen in PREREGISTRATION_survival.md (H1: niche 7 worse OS; H2:
niche 1 better OS), stated before this ran. Confirmatory tests are reported first
and separately; everything else is exploratory and corrected for selection.

79 events on 281 patients is thin (~11 events/covariate per model), so read the CIs,
not the point estimates.

    python scripts/run_basel_survival.py            # 1000 permutations
    python scripts/run_basel_survival.py --n-perm 200 --km-out km.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

DATA = Path(__file__).resolve().parents[1] / "data"
H1_NICHE, H2_NICHE = 7, 1  # pre-registered (immune-excluded, immune-rich)
DEATH_ALIVE_NOTE = "confirmatory, pre-registered"


def _benjamini_hochberg(pvals: list[float]) -> np.ndarray:
    ps = np.array(pvals, dtype=float)
    m = len(ps)
    order = np.argsort(ps)
    q = np.empty(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, ps[i] * m / (rank + 1))
        q[i] = prev
    return q


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--km-out", default=str(DATA / "km_niche7.png"))
    args = ap.parse_args()

    a = ad.read_h5ad(DATA / "basel_niched.h5ad")
    obs = a.obs.copy()
    obs["PID"] = obs["PID"].astype(str)

    # 1. per-patient niche abundance (281 x k)
    ab = pd.crosstab(obs["PID"], obs["niche"], normalize="index")
    ab.columns = [int(c) for c in ab.columns]
    pat = obs.drop_duplicates("PID").set_index("PID")[
        ["OSmonth", "event", "grade", "clinical_type"]
    ]
    pat = pat.loc[ab.index].copy()
    pat["clinical_type"] = pat["clinical_type"].astype(str).replace("nan", "unknown")
    pat["grade"] = pat["grade"].astype(str)
    n_events = int(pat["event"].sum())
    print(
        f"cohort: {len(pat)} patients, {n_events} events "
        f"(~{n_events/7:.0f} events/covariate; THIN, read the CIs)"
    )

    covdum = pd.concat(
        [
            pd.get_dummies(pat[c], prefix=c, drop_first=True).astype(float)
            for c in ["grade", "clinical_type"]
        ],
        axis=1,
    )
    covdum = covdum.loc[:, covdum.std(ddof=0) > 0]

    def cox(ab_vec, penalizer=0.0):
        d = pd.DataFrame(
            {
                "OSmonth": pat["OSmonth"].values,
                "event": pat["event"].astype(int).values,
                "abundance": (ab_vec - ab_vec.mean()) / ab_vec.std(ddof=0),
            },
            index=pat.index,
        )
        d = pd.concat([d, covdum], axis=1)
        for pen in ([penalizer, 0.1] if penalizer == 0 else [penalizer]):
            try:
                r = (
                    CoxPHFitter(penalizer=pen)
                    .fit(d, "OSmonth", "event")
                    .summary.loc["abundance"]
                )
                return dict(
                    hr=float(np.exp(r["coef"])),
                    lo=float(np.exp(r["coef lower 95%"])),
                    hi=float(np.exp(r["coef upper 95%"])),
                    p=float(r["p"]),
                )
            except Exception:
                continue
        return dict(hr=np.nan, lo=np.nan, hi=np.nan, p=np.nan)

    # 2-3. per-niche Cox (MLE) + BH
    res = {n: cox(ab[n]) for n in ab.columns}
    qs = _benjamini_hochberg([res[n]["p"] for n in ab.columns])
    for i, n in enumerate(ab.columns):
        res[n]["q"] = float(qs[i])

    print("\n=== CONFIRMATORY (pre-registered) ===")
    for niche, name, want in [
        (H1_NICHE, "H1 niche 7 (tumor-rich, immune-poor) -> WORSE OS", "hr>1"),
        (H2_NICHE, "H2 niche 1 (immune-rich) -> BETTER OS", "hr<1"),
    ]:
        r = res[niche]
        ok = (r["hr"] > 1 if want == "hr>1" else r["hr"] < 1) and r["p"] < 0.05
        print(
            f"{name}: HR/SD {r['hr']:.2f} [{r['lo']:.2f}, {r['hi']:.2f}] p={r['p']:.4f}  "
            f"=> {'CONFIRMED' if ok else 'not confirmed (direction ' + ('right' if (r['hr']>1)==(want=='hr>1') else 'wrong') + ')'}"
        )

    print("\n=== ALL NICHES (others EXPLORATORY) ===")
    print(f"{'niche':6}{'HR/SD':>8}{'95% CI':>20}{'p':>9}{'q(BH)':>9}")
    for n in sorted(ab.columns, key=lambda x: res[x]["p"]):
        r = res[n]
        tag = " *conf" if n in (H1_NICHE, H2_NICHE) else ""
        print(
            f"{n:<6}{r['hr']:>8.2f}   [{r['lo']:>6.2f},{r['hi']:>7.2f}]{r['p']:>9.4f}{r['q']:>9.4f}{tag}"
        )

    # 4. selection-aware permutation p (ridge 0.1 for stability across many fits)
    designs = {
        n: pd.concat(
            [
                pd.DataFrame(
                    {"abundance": (ab[n] - ab[n].mean()) / ab[n].std(ddof=0)},
                    index=pat.index,
                ),
                covdum,
            ],
            axis=1,
        )
        for n in ab.columns
    }
    outcome = pat[["OSmonth", "event"]].astype({"event": int})

    def best_p(oc):
        best = 1.0
        for n in ab.columns:
            dd = designs[n].copy()
            dd["OSmonth"] = oc["OSmonth"].values
            dd["event"] = oc["event"].values
            try:
                best = min(
                    best,
                    CoxPHFitter(penalizer=0.1)
                    .fit(dd, "OSmonth", "event")
                    .summary.loc["abundance", "p"],
                )
            except Exception:
                pass
        return best

    obs_best = best_p(outcome)
    rng = np.random.default_rng(0)
    arr = outcome.values
    ge = sum(
        best_p(
            pd.DataFrame(
                arr[rng.permutation(len(arr))],
                index=pat.index,
                columns=["OSmonth", "event"],
            )
        )
        <= obs_best
        for _ in range(args.n_perm)
    )
    print(
        f"\n=== SELECTION-AWARE p ({args.n_perm} perms, {len(ab.columns)} niches tested) ==="
    )
    print(
        f"observed best p {obs_best:.4f} | empirical p = {(ge+1)/(args.n_perm+1):.4f}"
    )

    # 6. KM niche 7
    med = ab[H1_NICHE].median()
    hi = ab[H1_NICHE] > med
    lr = logrank_test(pat.OSmonth[hi], pat.OSmonth[~hi], pat.event[hi], pat.event[~hi])
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    KaplanMeierFitter().fit(
        pat.OSmonth[hi],
        pat.event[hi],
        label=f"high niche-{H1_NICHE} (n={int(hi.sum())})",
    ).plot_survival_function(ax=ax, color="#c0392b")
    KaplanMeierFitter().fit(
        pat.OSmonth[~hi], pat.event[~hi], label=f"low (n={int((~hi).sum())})"
    ).plot_survival_function(ax=ax, color="#2980b9")
    ax.set_title(
        f"Niche {H1_NICHE} (tumor-rich, immune-poor) vs OS, log-rank p={lr.p_value:.4f}"
    )
    ax.set_xlabel("months")
    ax.set_ylim(0, 1.02)
    fig.tight_layout()
    fig.savefig(args.km_out, dpi=120)
    print(
        f"\n=== KM niche {H1_NICHE} === log-rank p={lr.p_value:.4f} (high n={int(hi.sum())}, low n={int((~hi).sum())}) -> {args.km_out}"
    )

    # 7. leave-one-core-out for H1
    obs["is1"] = (obs["niche"] == H1_NICHE).astype(int)
    obs["core"] = obs["core"].astype(str)
    n1 = obs.groupby(["PID", "core"])["is1"].sum().unstack(fill_value=0)
    tot = obs.groupby(["PID", "core"]).size().unstack(fill_value=0)
    n1p, totp = n1.sum(1), tot.sum(1)
    hrs = []
    for c in tot.columns:
        av = (
            ((n1p - n1[c]) / (totp - tot[c]))
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
        )
        av = av[(totp - tot[c]).reindex(av.index) > 0]
        if av.nunique() < 2:
            continue
        hrs.append(cox(av.reindex(pat.index).fillna(av.mean()))["hr"])
    hrs = np.array(hrs)
    print(f"\n=== LEAVE-ONE-CORE-OUT for H1 ({len(hrs)} cores) ===")
    print(
        f"HR range [{np.nanmin(hrs):.2f}, {np.nanmax(hrs):.2f}] | HR>1 in {int((hrs>1).sum())}/{len(hrs)} "
        f"-> {'direction robust' if (hrs > 1).all() else 'flips on some cores'}"
    )


if __name__ == "__main__":
    main()
