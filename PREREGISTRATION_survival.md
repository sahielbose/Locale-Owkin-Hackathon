# Pre-registration: niche -> overall-survival hypotheses

Written 2026-07-12T03:43:07Z (UTC), BEFORE computing any survival association (STEP 6).
Stated before seeing the result, so STEP 6 is a test, not a fishing expedition.

Cohort: Basel breast IMC, n = 281 patients, 79 death events.
Niches: window-only cellular neighborhoods over metacluster_id, k = 12, frozen in
data/basel_niched.h5ad (KMeans seed=0). Niche ids are stable in that file.

H1: niche 7 (79% tumor, 4% immune; dominant metacluster "HR low CK") -- the
    tumor-rich immune-EXCLUDED niche -- abundance predicts WORSE overall survival
    (Cox hazard ratio > 1).

H2: niche 1 (immune-rich / infiltrated; 59% immune, T cells dominant) abundance
    predicts BETTER overall survival (Cox hazard ratio < 1).

Analysis (pre-specified): per-patient niche abundance -> CoxPHFitter on OSmonth +
event, adjusted for grade + clinical_type; permutation test over the full
best-niche selection (1000x) for an honest p; Benjamini-Hochberg across all k niches.
