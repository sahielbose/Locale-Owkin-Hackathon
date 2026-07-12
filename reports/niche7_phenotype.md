# Niche 7 phenotype check: "immune-excluded" vs "immune desert"

**Status:** READ-ONLY analysis. No data written, no code changed, no graph rebuilt, no
survival re-run. Cohort: `data/basel_niched.h5ad` (755,070 cells, 289 cores, 281 PIDs),
loaded read-only.

## TL;DR

The name **"tumor immune-excluded"** for niche 7 is **overreaching and should be
renamed to "tumor-rich, immune-poor."** Niche 7 is genuinely a tumor-dominated,
immune-sparse *cell-level* neighborhood, but the *core-level* claim smuggled into the
word "excluded" — that immune cells are present but locked out — is not supported. A
meaningful fraction of niche-7-heavy cores are true **immune deserts** (no immune-rich
neighborhoods at all), and across the rest the excluded-vs-desert distinction is a
**smooth continuum with no bimodal break**, not a clean phenotype.

## What niche 7 and niche 1 actually are (composition confirmed)

| niche | tumor | immune | stroma | endo | reading |
|------:|------:|-------:|-------:|-----:|---------|
| **7** | 78.6% | **4.4%** | 15.1% | 1.9% | tumor-rich, immune-poor (currently "immune-excluded") |
| **1** | 18.1% | **59.4%** | 18.5% | 4.0% | immune-rich (the only immune-dominated niche) |

Niche 1 is unambiguously *the* immune-rich neighborhood; the next-most-immune niches
(6, 5, 10, 0) top out at 10–17% immune. So "does this core carry an immune-rich
compartment elsewhere?" is well operationalized as "does it contain niche 1?"

## The test

For each core I computed: total cells, niche-7 count/fraction, niche-1 count/fraction,
and core-level immune-major fraction. A core is **immune-excluded** if it is niche-7-heavy
*and* carries niche-1 elsewhere (immune present but compartmentalized); it is an
**immune desert** if it is niche-7-heavy and carries essentially no niche 1 (immune
absent). "Excluded" and "desert" both produce ~79%-tumor/~4%-immune niche-7 cells, so
the cell-level label alone cannot distinguish them — that is exactly the ambiguity being
tested.

**Primary threshold for "niche-7-heavy": niche-7 cells ≥ 10% of the core → 42 cores.**
(Robustness at ≥20% → 37 cores; ≥5% → 48 cores. All give the same story.)

## Result 1 — the excluded/desert split has no stable boundary

There is no natural threshold for "carries niche 1," so I swept several. The split slides
continuously with wherever the line is drawn (frac7 ≥ 10%, 42 cores):

| "immune present" defined as | EXCLUDED | DESERT |
|---|---:|---:|
| any niche-1 cell (>0) | 38 (90%) | **4 (10%)** |
| ≥ 10 niche-1 cells | 30 (71%) | **12 (29%)** |
| niche-1 ≥ 1% of core | 24 (57%) | **18 (43%)** |
| niche-1 ≥ 5% of core | 8 (19%) | **34 (81%)** |

Robustness at frac7 ≥ 20% (37 cores) is essentially identical: 89% / 70% / 57% / 22%
excluded. **A phenotype that flips from 90%-excluded to 81%-desert as you nudge the
presence threshold is not a clean dichotomy — the threshold-sensitivity itself is the
signature of a continuum.**

## Result 2 — the distribution is a continuum, not bimodal

Sorted niche-1 fraction across the 42 niche-7-heavy cores:

```
0, 0, 0, 0, .0003, .0003, .0006, .0009, .0012, .0017, .002, .0036, .0048, .0053,
.0063, .0076, .0097, .0099, .011, .0117, .0157, .0175, .019, .0216, .0223, .0301,
.0304, .0321, .0336, .0361, .0379, .0421, .0456, .0465, .0517, .0616, .0863, .1307,
.1343, .1703, .2196, .2883
```

Histogram (niche-1 fraction bins): 0=4, (0,0.5%]=9, (0.5,1%]=5, (1,2%]=5, (2,5%]=11,
(5,10%]=3, (10,20%]=3, >20%=2. The core-level **immune-major** fraction tells the same
story: min 0.85%, median 6.75%, max 35%, unimodal and right-skewed with no gap.

**There is no bimodal valley separating an "excluded" cluster from a "desert" cluster.**
Cores fill the whole range from zero immune to immune-rich smoothly.

## Result 3 — the deserts are real, not sampling artifacts

The 4 cores with exactly zero niche-1 cells are **large and overwhelmingly tumor**, so
their emptiness is a real biological absence, not a small-n fluke:

| core | n_cells | niche-7 frac | immune-major frac |
|---|---:|---:|---:|
| BaselTMA_SP41_135_X8Y5 | 3,326 | 0.932 | 1.6% |
| BaselTMA_SP42_25_X3Y2 | 3,866 | 0.974 | 1.0% |
| BaselTMA_SP42_74_X15Y3 | 2,696 | 0.991 | 2.3% |
| BaselTMA_SP43_32_X4Y6 | 3,905 | 0.963 | 0.8% |

These are ~3,000-cell tumors that are 93–99% niche-7 with 1–2% immune and zero
immune-rich neighborhoods. That is an **immune desert**, the biological opposite of
"excluded" (which implies immune cells amassed at a border and held out). Calling these
"immune-excluded" would be flatly wrong for at least ~10% of niche-7-heavy cores — and
for 29–43% under any non-trivial presence threshold (≥10 cells / ≥1%).

## Verdict

- **Excluded vs desert does not separate cleanly** — it is a continuum (Results 1 & 2).
- **A meaningful fraction are genuine deserts, not excluded** — 10% (any-cell) up to
  ~30–43% (≥10 cells / ≥1%) of niche-7-heavy cores, robust to the niche-7 threshold, and
  the extreme deserts are large real tumors (Result 3).

By the decision rule (a meaningful desert fraction **and** a continuum rather than clean
exclusion), **"immune-excluded" is a false/overreaching label.**

## Recommendation (do NOT auto-apply — main thread should perform)

**Rename niche 7 from "tumor immune-excluded" to "tumor-rich, immune-poor."**

Reasoning — why a false label is dangerous here, not just imprecise:

1. **It claims a mechanism the cell-level data cannot see.** "Excluded" (Chen & Mellman,
   *Nature* 2017) is a specific spatial mechanism — immune cells present but barred from
   the tumor bed — and is clinically distinct from "desert" (no immune cells), with
   different immunotherapy implications. Niche 7 is a *cell-level* neighborhood label;
   "excluded" is a *core/tumor-level* phenotype. You cannot infer the core phenotype from
   the cell label, and when we actually tested it, the core phenotype was a continuum that
   included real deserts. "Tumor-rich, immune-poor" states exactly and only what the
   composition supports.

2. **`characterize_niche` hands the name straight to K Pro, which reasons from it and
   cannot re-check it.** The downstream agent has no access to the per-core niche-1
   distribution; it sees the string "immune-excluded" and will treat exclusion as an
   established fact — potentially reasoning toward mechanisms or therapy framings (e.g.
   TGF-β/stromal barrier, anti-PD-L1 rationale) that apply to exclusion but not to
   deserts. A wrong label doesn't just lose precision; it manufactures a false premise the
   next reasoner builds on. An accurate, deliberately modest name removes that failure
   mode.

## Future work (NOT done here; flagged only)

Whether excluded-vs-desert is itself **prognostic** is a real, separate question — but
testing it on Basel would be fishing, because Basel has already been used to define and
name these niches. It would require pre-registration and evaluation on the untouched
**Zurich cohort** (`Data_publication/ZurichTMA/`, ~70 patients). Not performed in this
check.

---
*Methods: per-core niche composition from `obs['niche']`/`obs['core']`; immune fraction
from `obs['major']`. Spatial graph in `obsp` was not rebuilt or used (composition-only).
No survival analysis. `data/` untouched.*
