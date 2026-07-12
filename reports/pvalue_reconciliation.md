# P-value reconciliation: niche 3 "0.070" vs selection "observed best p = 0.093"

Scope: read-only audit of `scripts/run_basel_survival.py` against the cached
cohort `data/basel_niched.h5ad` (281 patients, 79 events, k = 12 niches). No files
were changed except this report; the spatial graph was not rebuilt.

## The apparent contradiction

- The exploratory per-niche Cox table prints niche 3 at **p = 0.070**.
- The selection-aware permutation prints **observed best p = 0.093**.

If the permutation searched all 12 niches (it does), the "best" observed p should
be the smallest per-niche p, so it looks impossible for niche 3 to sit at 0.070
while the searched minimum is 0.093.

## Root cause: two different Cox model specs (penalizer 0.0 vs 0.1)

The two blocks fit the *same* design (same per-SD z-scored abundance, same
`grade + clinical_type` dummy covariates, same outcome) but with a **different
penalizer**. That single parameter is the entire discrepancy.

**Exploratory table = unpenalized MLE (penalizer = 0.0).**

```
77  def cox(ab_vec, penalizer=0.0):                 # default is MLE
...
87      for pen in ([penalizer, 0.1] if penalizer == 0 else [penalizer]):
88          try:
89              r = (CoxPHFitter(penalizer=pen).fit(...) ...)
...
105 res = {n: cox(ab[n]) for n in ab.columns}        # table uses cox() with default
```

`cox()` tries `pen = 0.0` first and only falls back to `0.1` if that fit throws.
For this cohort all 12 niches converge under MLE, so **every table p-value is the
unpenalized MLE p** (the `0.1` fallback never fires). Niche 3 MLE p = 0.0700 -> prints "0.070".

**Selection permutation = ridge (penalizer = 0.1), always.**

```
131 # 4. selection-aware permutation p (ridge 0.1 for stability across many fits)
...
147 def best_p(oc):
148     best = 1.0
149     for n in ab.columns:                         # loops ALL 12 niches
...
156             CoxPHFitter(penalizer=0.1)            # always ridge 0.1
157             .fit(dd, "OSmonth", "event")
158             .summary.loc["abundance", "p"]
...
164 obs_best = best_p(outcome)
```

`best_p` deliberately uses `penalizer=0.1` for numerical stability across the
~12 x 1001 fits it performs. Under ridge, niche 3's p = 0.0929 -> prints "0.093".

## Proof by recomputation (both specs, all 12 niches)

Recomputed from `data/basel_niched.h5ad` reproducing the script's preprocessing:

```
niche  MLE p (pen=0.0)   ridge p (pen=0.1)
0          0.6816            0.6667
1          0.0458            0.1503
2          0.3788            0.3652
3          0.0700            0.0929   <- niche 3
4          0.2199            0.5163
5          0.4427            0.6608
6          0.5415            0.7050
7          0.5252            0.6588
8          0.3056            0.2752
9          0.9772            0.9693
10         0.8866            0.9842
11         0.0726            0.1087

Table (MLE-first)  best niche = 1, p = 0.0458
Permutation best_p (ridge)  best niche = 3, observed best p = 0.0929  ->  0.093
```

Both reported numbers reproduce exactly: niche 3 MLE = 0.0700 ("0.070") and the
ridge search minimum = 0.0929 ("0.093"). The permutation searched all 12 niches;
niche 3 is genuinely the ridge-model minimum.

### Candidate causes explicitly ruled out

- **Not niche exclusion.** `best_p` (lines 149-161) iterates every `ab.columns`
  entry; the recompute shows the ridge minimum is niche 3 = 0.0929, i.e. all 12
  were evaluated and none was dropped.
- **Not different covariates or z-scoring.** Both blocks build the identical
  `covdum` (`grade` + `clinical_type` drop-first dummies, lines 68-75) and the
  identical per-SD standardization `(ab[n]-mean)/std(ddof=0)` (lines 82 and 136-139).
- **Not a bug.** The ridge choice is intentional and labeled (line 131). The only
  defect is narrative: the printed exploratory p (MLE) and the printed "observed
  best p" (ridge) come from different model specs, so side by side they read as a
  contradiction.

Side note: the penalty also reorders the argmin. Under MLE the smallest p is
niche 1 (0.046, a pre-registered confirmatory niche) with niche 3 second (0.070);
ridge shrinks niche 1 harder (0.046 -> 0.150), so niche 3 (0.093) becomes the
ridge search minimum. That is why the permutation's "best" happens to land on the
same niche 3 quoted in the exploratory row.

## Which number is right, and the internally consistent statement

Both p-values are correct for their own model; they are simply not the same model.

For the **selection-corrected permutation**, the correct observed statistic is the
**ridge best p = 0.093 (niche 3)**, because the permutation null in the same block
is generated with the identical ridge (penalizer = 0.1) Cox model. Comparing a
ridge observed statistic against a ridge null is the like-for-like test, and that
is what yields the empirical p. The **0.070** is the *unpenalized MLE* p for
niche 3 as printed in the exploratory table; it is not the statistic the
permutation minimized over and should not be read as the "best p" the permutation
searched.

Corrected, internally consistent wording:

> Exploratory per-niche Cox (unpenalized MLE): the two smallest p-values are
> niche 1 (p = 0.046) and niche 3 (p = 0.070). The selection-aware permutation
> uses a ridge Cox (penalizer = 0.1) for numerical stability across the ~12,000
> fits it performs; under that identical model the best of the 12 niches is
> niche 3 at observed p = 0.093, and comparing it to the ridge-model permutation
> null gives a selection-corrected empirical p = 0.44. The niche-3 gap between the
> two blocks (0.070 vs 0.093) is entirely the penalizer (MLE vs ridge 0.1), not
> niche exclusion, a covariate change, or a scaling change; all 12 niches were
> searched.

Note: the selection-corrected conclusion is unaffected: empirical p = 0.44 is
nowhere near significance under either reading (MLE-best 0.046 or ridge-best 0.093
both wash out after correcting for searching all 12 niches).
