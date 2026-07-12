# Locale

**Spatial reasoning for AI scientists.** An MCP server that gives an agent access to tissue geography, and a calibrated statement of how far that geography can be trusted for any given question.

Built for the Owkin Rewiring Biology Hackathon. Validated on the Basel breast cancer cohort of Jackson, Fischer et al. (Nature, 2020).

An AI reading a tumor sees a composition vector. A pathologist sees a picture. Two patients can have identical cell-type proportions and opposite outcomes, and the only thing that separates them is where the cells sit. Locale supplies the map. More importantly, it supplies the limits of the map: the failure mode that destroys an AI Scientist in a regulated setting is not a missed finding, it is a confident report of a finding that is not there. Locale is built so that its refusals are as legible as its results.

![Basel cores at recovered coordinates, coloured by major class](docs/figures/cores_spatial.png)

*Four cores from the cohort, every cell drawn at its recovered coordinate and coloured by major class (tumor red, stroma green, immune blue, endothelial purple). The leftmost core is 97% tumor and shows genuine tissue nests and voids; the rightmost is 32% tumor and shows immune and tumor cells separating into compartments. A wrong coordinate join would produce uniform confetti in every panel. This is the geography that a composition vector cannot see, recovered from segmentation masks and verified by exact per-core cell-count agreement.*

### Three tests, three passes

| Test | What it checks | Result |
| --- | --- | --- |
| Positive control | Recovers established biology it was never told about | Immune infiltration predicts better survival, HR 0.76, pre-registered |
| External validation | Rediscovers a Nature paper's spatial communities, cold | Adjusted Rand Index 0.400 over 342,662 cells |
| Negative control | Given 12 niches and only 79 events, declines to report a biomarker | Selection-aware empirical p = 0.44 |

The full write-up is in [docs/Locale_Technical_Report.pdf](docs/Locale_Technical_Report.pdf).

---

## Table of contents

1. [Motivation](#1-motivation)
2. [What Locale is](#2-what-locale-is)
3. [Results at a glance](#3-results-at-a-glance)
4. [Data provenance and retrieval](#4-data-provenance-and-retrieval)
5. [The coordinate problem](#5-the-coordinate-problem)
6. [Cell-type ontology](#6-cell-type-ontology)
7. [Cohort construction](#7-cohort-construction)
8. [Spatial graph](#8-spatial-graph)
9. [Neighborhood enrichment](#9-neighborhood-enrichment)
10. [Niche discovery](#10-niche-discovery)
11. [External validation](#11-external-validation)
12. [Pre-registration and survival analysis](#12-pre-registration-and-survival-analysis)
13. [The statistical honesty layer](#13-the-statistical-honesty-layer)
14. [Register of negative results](#14-register-of-negative-results)
15. [Reproducibility](#15-reproducibility)
16. [Repository layout](#16-repository-layout)
17. [Using the MCP server](#17-using-the-mcp-server)
18. [Limitations](#18-limitations)
19. [Citation, license, and contributors](#19-citation-license-and-contributors)

---

## 1. Motivation

Composition is not organisation. A tumor that is 65% malignant cells, 12% immune cells, and 23% stroma may be one in which lymphocytes have infiltrated the malignant compartment and are killing it, or one in which they are held at the margin and are doing nothing. These two tissues produce the same composition vector and the same bulk expression profile. They do not produce the same patient.

Agentic tools over multi-omics data read the composition. They cannot read the map. Locale supplies the map, and a calibrated statement of how much that map can be trusted for any particular question. The second point is the product.

## 2. What Locale is

Locale is a Model Context Protocol (MCP) server. The analysis engine underneath it is a package of pure functions over a single `AnnData` object and contains no MCP code at all; the MCP tools are thin wrappers. That separation is what makes the science testable without a running server.

The server exposes nine tools:

| Tool | Returns |
| --- | --- |
| `list_samples` | one record per image / core |
| `describe_sample` | cohort or single-image summary |
| `compute_enrichment` | neighborhood-enrichment matrix (cell type by cell type) |
| `find_niches` | the discovered niches with composition and marker program |
| `characterize_niche` | one niche in detail |
| `find_prognostic_niches` | niches ranked by survival association |
| `describe_niches` | the frozen niche catalog with honest names |
| `correlate_niche_outcome` | a niche's survival association plus the full statistical context |
| `get_map_payload` | tissue-map coordinates for the viewer |

The last one is the point of the whole exercise. `correlate_niche_outcome` never returns a bare p-value; it returns the hazard ratio, its confidence interval, the number of hypotheses tested, the Benjamini-Hochberg q-value, the selection-aware permutation p, the event count, the minimum hazard ratio the cohort can resolve at 80% power, and a plain verdict. See [section 13](#13-the-statistical-honesty-layer).

## 3. Results at a glance

| Quantity | Value |
| --- | --- |
| Cells analysed | 755,070 |
| Tumor cores | 289 |
| Patients | 281 |
| Death events | 79 |
| Overall survival range | 0 to 233 months, 0 missing |
| Cell types | 27 metaclusters (25 names), 4 major classes |
| Immune exclusion | tumor and immune cells strongly anti-adjacent (see figure below) |
| Niches discovered (unsupervised) | 12, k selected by stability ARI = 0.67 |
| External agreement | ARI 0.400 versus 23 published tumor community phenotypes |
| Positive control (H2) | niche 1 immune-rich, HR 0.76 [0.59, 1.00], p = 0.046 |
| Pre-registered negative (H1) | niche 7 excluded, HR 1.08 [0.86, 1.34], p = 0.525 |
| Multiplicity | BH q > 0.29 for every niche |
| Selection-aware | empirical p = 0.44 (1,000 permutations) |

All figures below are regenerated from the committed pipeline by [`scripts/make_figures.py`](scripts/make_figures.py); none is hand-drawn.

## 4. Data provenance and retrieval

Source: Jackson, H.W., Fischer, J.R., et al. (2020), "The single-cell pathology landscape of breast cancer," Nature 578:615 to 620. Imaging mass cytometry on a breast cancer tissue microarray, 35 antibody channels, single-cell segmentation. Archive: Zenodo record `10.5281/zenodo.3518284`.

Neither archive (a combined 48 GB) was downloaded. Both were read by HTTP range request. A ZIP stores its central directory at the end of the file and each member occupies a contiguous compressed byte span, so with `Range` headers one can enumerate the contents and extract individual members without transferring the rest; Zenodo returns HTTP 206, confirming support. Two details make this correct rather than accidentally full-transfer:

1. The naive approach (open a member and read a kilobyte) issues a range request from the member offset to the end of file, which for an early member transfers nearly everything. We instead read the 30-byte local file header, parse the filename and extra-field lengths, compute the exact compressed span, and request precisely that.
2. The central directory's extra-field length can differ from the local header's. The local header must actually be read.

Total transferred: approximately 5 GB out of 48 GB. Everything else is OME-TIFF image stacks and MATLAB sessions that are not used. See [`scripts/download_data.py`](scripts/download_data.py).

## 5. The coordinate problem

The released marker table `SC_dat.csv` is in long format with five columns (`core, CellId, id, channel, mc_counts`) and one row per cell-channel pair. There is no x, no y, and no coordinate file anywhere in the archive. The coordinates existed in the authors' MATLAB pipeline; they were never exported in tabular form.

An IMC segmentation mask is a label image whose pixel value is the `CellId` within that core. Therefore

```
regionprops(mask) -> (label, centroid) = (CellId, y, x)
```

recovers coordinates directly, with no inference. Two traps: `regionprops` returns the centroid as `(row, col)` which is `(y, x)`, and swapping them silently transposes every tissue map; and the masks are not loose members but sit inside a nested archive that must be range-fetched whole. See [`scripts/extract_coords.py`](scripts/extract_coords.py).

Mapping mask filenames to core identifiers is a pure transform, but the core number is the integer immediately preceding the `X..Y..` token, not the trailing integer. A first attempt keyed on the trailing token overlapped one real core; the corrected transform mapped 352, and the `X..Y..` tokens intersected 130 of 130, which proved it was a transform rather than a lookup.

**Join verification, done two ways before any analysis:**

- Injectivity. 376 mapped masks resolve to 376 distinct cores; the mapping is one to one.
- Per-core cell-count agreement. For every core, the cell count recovered from the mask equals the cell count in the PhenoGraph table exactly. All 376 of 376 cores agree (relative difference 0.0), including the 24 cores that required a metadata snap. A mis-mapped core would have carried a different count and been caught here. None was.

Recovered: 844,498 cells across 376 of 376 cores, 0 duplicate ids, 0 NaN coordinates, densities of roughly 2,500 cells per square millimetre. Physically credible before any statistical test.

## 6. Cell-type ontology

The PhenoGraph table assigns each cell to one of 71 clusters. These are not the granularity at which the paper's biology is stated. Jackson et al. collapse them into 27 curated metaclusters, hardcoded in their own pipeline. We transcribed that map verbatim (see [`src/localespatial/metaclusters.py`](src/localespatial/metaclusters.py) and Appendix A of the report); all 71 clusters map and none is orphaned.

The 27 metaclusters group into four major classes:

| Major class | Metacluster ids | PhenoGraph clusters | Cells |
| --- | --- | --- | --- |
| immune | 1 to 6 | 6 | ~122,000 |
| endothelial | 7 | 1 | ~20,000 |
| stroma | 8 to 13 | 6 | ~214,000 |
| tumor | 14 to 27 | 58 | ~399,000 |

This matters more than it looks. 58 of the 71 PhenoGraph clusters are tumor subtypes, so any analysis at PhenoGraph granularity dissolves: two adjacent cells in one tumor nest routinely land in different clusters, so "same type" becomes almost impossible and genuine spatial coherence is diluted to noise. All clustering is performed on metacluster id (27), never on names (25, because two pairs share a label), because the published communities we validate against were computed over the 27.

## 7. Cohort construction

`diseasestatus` has two levels, tumor (289 cores) and non-tumor (87 cores); the 87 non-tumor cores were dropped. Overall-survival event coding was read from the data rather than guessed: `Patientstatus` is a four-level string, and we set `event = 1` for both death levels and `0` for both alive levels.

| Quantity | Value |
| --- | --- |
| Cells | 755,070 |
| Tumor cores | 289 |
| Patients | 281 |
| Death events | 79 |
| Overall survival | 0 to 233 months, 0 missing |

See [`scripts/build_basel.py`](scripts/build_basel.py).

## 8. Spatial graph

```python
squidpy.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True, library_key="core")
```

`coord_type="generic"` because IMC cells are not on a lattice. `library_key="core"` is the single most dangerous parameter in the pipeline: this is a tissue microarray, cells in different cores are physically unconnected tissue, and a graph built without partitioning by core fabricates edges between them, producing garbage niches with no error anywhere. We therefore assert the property rather than assume it. A guard function enumerates all edges and raises if any connects two cells from different cores.

**Measured cross-core edges: 0.** See [`src/localespatial/engine/graph.py`](src/localespatial/engine/graph.py).

## 9. Neighborhood enrichment

For each ordered pair of cell types, `squidpy.gr.nhood_enrichment` compares the observed number of graph edges connecting them against a null generated by permuting cell-type labels, and reports a z-score. This is composition-normalised: it answers whether two types sit adjacent more than chance given how much of each is present, which a raw same-type-neighbour fraction cannot do.

![Neighborhood enrichment over the 27 metaclusters](docs/figures/enrichment_heatmap.png)

The block-diagonal structure (like sits next to like) is the signature of real tissue: every diagonal block is strongly positive. Immune self-association is emphatic (self-enrichment z up to +498 for the B and T aggregates), endothelial self-enrichment is +82 despite endothelium being only about 3% of cells (vessels are linear structures and should self-associate), and stroma is +35. Tumor self-enrichment is compressed to +8 because tumor is the majority class and has the least room to exceed its own permutation null.

The off-diagonal is where the biology sits. Tumor and immune cells avoid each other far beyond what their abundances would predict.

![Major-class enrichment blocks](docs/figures/enrichment_major_blocks.png)

The single strongest exclusion is the HR-low-CK tumor phenotype against T cells, at z = -128. That is exactly the phenotype that dominates the immune-excluded niche 7 (section 10), so the enrichment matrix and the niche discovery corroborate each other from two independent computations. Immune exclusion is measured directly, in real tissue, composition-normalised, with no reference to any published result, and it emerged on the first honest run.

(The bundled technical report quotes immune-block aggregates from an earlier run; the numbers above and every figure here are regenerated from the committed `data/basel_niched.h5ad` by `scripts/make_figures.py`.)

## 10. Niche discovery

For each cell $i$ with spatial neighbours $N(i)$ (self-inclusive), define the neighborhood composition $w_i \in \mathbb{R}^{27}$ as the fraction of $N(i)$ in each metacluster:

$$w_i = \frac{1}{|N(i)|} \sum_{j \in N(i)} e_{c(j)}, \qquad e_{c(j)} \in \{0, 1\}^{27}$$

computed as $W C \oslash (W \mathbf{1})$ where $W$ is the adjacency and $C$ the one-hot cell-type matrix. A niche is a cluster in $w$-space: a recurring kind of neighborhood. Clustering is k-means with `n_init=10` and a fixed seed. See [`src/localespatial/engine/niches.py`](src/localespatial/engine/niches.py).

**A documented failure: the identity block.** Our first specification concatenated each cell's own one-hot identity to its window, giving a 54-dimensional feature. The resulting niches were degenerate: they were the metaclusters relabelled. The cause is metric geometry. The identity block is a one-hot of magnitude 1, while the window is a distribution spread over 27 entries; in squared Euclidean distance the identity swamps every neighborhood difference, so k-means partitions by label. Worse, this passed the subsample-stability check, because a clustering that partitions by label is perfectly reproducible. Stability is not validity. The final feature is the window alone (27-dimensional, self-inclusive), the canonical Schurch and Nolan cellular-neighborhood formulation.

`k = 12` was selected by subsample stability (ARI between clusterings of independent resamples), stability ARI = 0.67.

![Niche composition](docs/figures/niche_composition.png)

The full immune-exclusion gradient emerges without supervision:

| Niche | Character | tumor | immune | stroma | cells | cores | dominant type |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | immune-rich (infiltrated) | 0.18 | 0.59 | 0.18 | 92k | 260 | T cells (47%) |
| 6 | tumor / immune boundary | 0.66 | 0.16 | 0.15 | 107k | 275 | Basal CK, hypoxic |
| 10 | tumor / immune boundary | 0.80 | 0.12 | 0.06 | 19k | 24 | p53 EGFR (70%) |
| 7 | tumor, immune-excluded | 0.79 | 0.04 | 0.15 | 75k | 160 | HR low CK (67%) |
| 0, 5 | stromal / vascular | ~0.27 | ~0.14 | ~0.56 | | | Vimentin, elongated |

Niche 7 spans 160 of the 289 tumor cores and 75,000 cells, so it is a population and not an artefact. Its dominant constituent is the hormone-receptor-low tumor phenotype, which is the clinically aggressive one. An immune-excluded compartment made of the aggressive subtype is internally coherent biology, and it motivated hypothesis H1.

## 11. External validation

Jackson et al. performed a conceptually equivalent analysis and archived their results. Basel is therefore the one spatial oncology cohort in which a tool of this kind can be scored rather than merely exhibited. Their answers were kept sealed until the engine had produced its own.

Their ground truth is the recurring community phenotype, reached by a two-hop join: per-cell fine community to phenotype on `(core, Community)`, yielding a coarse label in 1 to 23. A first attempt joined on the 1,940-way fine community, which against our 12-way partition is mechanically meaningless and returned ARI 0.003; we reported nothing and investigated. The corrected two-hop join matched 100% of their tumor cells. See [`scripts/run_basel_groundtruth.py`](scripts/run_basel_groundtruth.py).

**ARI = 0.400, over 342,662 tumor cells,** our 12 niches against their 23 published tumor community phenotypes. Two different algorithms, two different feature spaces, two different values of k. Chance agreement is approximately zero; perfect agreement would indicate we had accidentally reimplemented their method. 0.400 is the signature of a genuine and independent rediscovery. The agreement is structured, not smeared: niche 7 concentrates into just 2 of their 23 communities, accounting for 85% of its cells.

We also compared our enrichment matrix against the authors' published neighborhood heatmap and got full Pearson r = 0.38 but off-diagonal r near 0.004. The correlation is carried entirely by the diagonal. The two statistics answer different questions in the tumor block, so we report this comparison as inconclusive and exclude it from our claims. We do not present r = 0.38 as agreement. The ARI stands alone.

## 12. Pre-registration and survival analysis

Before any survival model was fitted, both hypotheses were written to [`PREREGISTRATION_survival.md`](PREREGISTRATION_survival.md) and timestamped. The file has not been modified since.

- **H1**: niche 7 abundance (tumor-rich, immune-depleted, HR-low-CK dominant) predicts worse overall survival.
- **H2**: niche 1 abundance (immune-rich, infiltrated) predicts better overall survival.

These are confirmatory and carry no multiple-testing penalty. Every other niche is exploratory.

**Feature construction.** For each patient $p$ and niche $j$, the abundance $A_{pj}$ is the fraction of that patient's cells in niche $j$. This is a 281 by 12 matrix, and it is the entire candidate biomarker.

**Model.** One Cox proportional-hazards model per niche, adjusted for grade and clinical subtype, hazard ratios reported per standard deviation of abundance:

$$h(t \mid A_{\cdot j}) = h_0(t)\, \exp\!\big(\beta_j A_{\cdot j} + \gamma_1\,\mathrm{grade} + \gamma_2\,\mathrm{clinical\_type}\big)$$

**Confirmatory results.**

| Hypothesis | HR / SD | 95% CI | p |
| --- | --- | --- | --- |
| H2: niche 1 (infiltrated) to better OS | 0.76 | [0.59, 1.00] | 0.046 |
| H1: niche 7 (excluded) to worse OS | 1.08 | [0.86, 1.34] | 0.525 |

H2 is confirmed, and it is a calibration rather than a discovery: tumor-infiltrating lymphocytes have been prognostic in breast cancer for over a decade. We present it as evidence that the instrument fires when there is real signal. The interval's upper bound sits exactly at 1.00, and we state it as the marginal result it is. H1 failed: the direction is as pre-registered but the effect is null and the interval straddles unity.

![Per-niche survival forest plot](docs/figures/survival_forest.png)

**The Kaplan-Meier inverts the sign.** The unadjusted median split of niche-7 abundance puts the high-exclusion arm slightly above the low arm (better survival), which is the opposite direction to both H1 and the adjusted Cox estimate. Both are null, so they are not in formal contradiction, but the sign is unstable across specifications, and that instability is itself the finding: it is exactly what one observes when there is no underlying effect and the estimate is driven by noise and covariate adjustment.

![Kaplan-Meier by niche 7 abundance](docs/figures/km_niche7.png)

**Multiplicity.** Benjamini-Hochberg across all twelve niches gives q > 0.29 everywhere. Nothing survives correction.

**Selection-aware permutation.** A BH correction still assumes the twelve tests are the tests you meant to run. The stronger question: given that we would have reported whichever niche looked best, how surprising is our best result? Permute the survival labels across the 281 patients, refit all twelve Cox models, record the best p, repeat 1,000 times. The empirical p is the fraction of permuted runs whose best p beats ours.

```
observed best p = 0.093   =>   empirical p = 0.44
```

Our best exploratory niche is entirely consistent with noise once one accounts for having tested twelve.

**Sensitivity.** Leave-one-core-out, refitting the H1 model with each of the 289 cores removed in turn: the direction is stable (HR > 1 in 289 of 289 refits) but the magnitude never exceeds about 1.06 and is never significant. No single core carries the effect, because there is no effect to carry.

**Power.** 79 events across 281 patients is roughly 6 to 7 events per covariate at k = 12. At 80% power this cohort can only resolve hazard ratios beyond about 1.37, and the observed effects sit inside that band. The spatial signal is strong and real; at this power it does not translate into a defensible survival effect, and Locale does not report one. It does not report the absence of a biomarker either. It reports the limits of the evidence. That distinction is the difference between a measurement and an overclaim, and it is the reason the tool exists.

## 13. The statistical honesty layer

`correlate_niche_outcome` is the tool that carries the whole thesis. It never returns a bare point estimate. For niche 1, the positive control, over the MCP protocol it returns:

```json
{
  "niche_id": 1,
  "hazard_ratio": 0.764,
  "ci_95": [0.587, 0.995],
  "p_raw": 0.0458,
  "n_hypotheses_tested": 12,
  "q_fdr": 0.2903,
  "p_selection_aware": 0.4396,
  "n_events": 79,
  "min_detectable_hr": 1.371,
  "verdict": "insufficient evidence"
}
```

The raw p is 0.046, which an expression-only tool would report as a hit. The four fields an agent cannot compute for itself (`n_hypotheses_tested`, `q_fdr`, `p_selection_aware`, `min_detectable_hr`) ship with every finding, and they turn a plausible p-value into an honest verdict. The minimum detectable hazard ratio is a Schoenfeld power calculation: with 79 events, effects inside [1/1.37, 1.37] are underpowered. See [`src/localespatial/engine/outcome.py`](src/localespatial/engine/outcome.py) and the demo transcript in [`demo/transcript.md`](demo/transcript.md).

## 14. Register of negative results

Recorded in full, because a methods document that lists only what worked is an advertisement.

1. Coordinates absent from the marker table. Discovered at hour zero by inspecting the header; the project was gated on this and did not proceed until it was resolved.
2. Mask filename regex, first attempt. Keyed on the trailing integer, overlapped one real core. Corrected to the integer preceding the XY token, 352.
3. The sum-of-squares null for the neighbour test. Reported 6.7x enrichment. It assumes global composition and ignores per-core composition. Discarded.
4. The same-type-neighbour statistic. Even against the correct within-core scramble it returns about 1.07x, which met our pre-agreed stopping condition. It is uninformative for a majority-dominated cohort (the ratio decays to exactly 1.00 as cores become homogeneous, by construction). Superseded by composition-normalised enrichment.
5. The identity-plus-window niche feature. Produced degenerate niches that were metaclusters in disguise, and passed the stability check while doing so. Replaced with window-only.
6. ARI against the fine community id. Returned 0.003, a 1,940-way against a 12-way partition. A category error, corrected by the two-hop join.
7. Enrichment matrix correlation. r = 0.38 overall, near 0.004 off-diagonal. Inconclusive, excluded from all claims.
8. Hypothesis H1. Pre-registered, and it failed. Not retuned, not repaired, not converted into a post-hoc ratio until something crossed 0.05.

Two of these (items 4 and 5) are cases where a metric actively lied and we caught it. They are documented in [`reports/`](reports/).

## 15. Reproducibility

The engine is a package of pure functions over a single `AnnData` object and contains no MCP code; the MCP tools are three-line wrappers. It consumes any `AnnData` carrying `obsm['spatial']`, a cell-type column, and a patient column, so pointing it at a different dataset requires a loader and no change to the analysis code.

| Component | Library |
| --- | --- |
| spatial graph, neighborhood enrichment | squidpy |
| data structures, preprocessing | scanpy, anndata |
| niche clustering, ARI | scikit-learn |
| survival models | lifelines |
| mask centroids | scikit-image, tifffile |
| partial archive retrieval | remotezip, requests |

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .                       # installs the localespatial package
python scripts/make_mock.py            # writes the committed data/mock.h5ad
pytest -q                              # 44 tests, all green

# real cohort (requires the extracted CSVs; see scripts/download_data.py)
python scripts/build_basel.py          # builds data/basel_niched.h5ad
python scripts/run_basel_niches.py     # niche discovery + enrichment sanity
python scripts/run_basel_survival.py   # pre-registered survival analysis
python scripts/make_figures.py         # regenerates docs/figures/
```

Note on the import path: the package is named `localespatial`, not `locale`, because a top-level `locale` package shadows the Python standard library `locale` module that `gettext` (inside click and uvicorn) imports, which breaks the server. The product is still called Locale; only the import path changed.

## 16. Repository layout

```
Locale/
  README.md                     this document
  LICENSE                       MIT
  pyproject.toml                package + tooling config
  requirements.txt
  PREREGISTRATION_survival.md   timestamped, never modified
  src/localespatial/
    schema.py                   shared Pydantic v2 data contract
    metaclusters.py             the 71 to 27 PhenoGraph map
    engine/                     pure analysis functions (Lane A)
      graph.py  enrichment.py  niches.py  characterize.py
      outcome.py  groundtruth.py  validate.py
    mcp_server/                 thin MCP wrappers (Lane B)
      server.py  tools.py  interpret.py
    viz/                        tissue-map viewer (Lane C)
  scripts/                      data build, analysis runs, figures
  tests/                        44 tests (engine, server, integration, schema)
  demo/                         frozen findings + K Pro transcript
  reports/                      two negative-result write-ups
  docs/
    Locale_Technical_Report.pdf full methods and statistics
    vision.md                   original project vision
    figures/                    regenerated by scripts/make_figures.py
  data/                         gitignored (only mock.h5ad is committed)
```

## 17. Using the MCP server

The server defaults to the streamable-http transport for remote use (for example as a K Pro custom connector, an `https .../mcp` URL). For a local desktop client it speaks stdio.

Run it directly:

```bash
python -m localespatial.mcp_server.server                          # streamable-http on 127.0.0.1:8000/mcp
LOCALE_TRANSPORT=stdio python -m localespatial.mcp_server.server   # stdio (desktop clients)
LOCALE_DATA=/path/to/basel_niched.h5ad LOCALE_TRANSPORT=stdio python -m localespatial.mcp_server.server
```

For Claude Desktop, add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "locale": {
      "command": "/absolute/path/to/.venv/bin/python",
      "args": ["-m", "localespatial.mcp_server.server"],
      "env": {
        "LOCALE_TRANSPORT": "stdio",
        "LOCALE_DATA": "/absolute/path/to/data/basel_niched.h5ad"
      }
    }
  }
}
```

All logging goes to stderr so it cannot corrupt the stdio JSON-RPC stream. Without `LOCALE_DATA` the server serves the committed mock.

## 18. Limitations

- Power. 79 events is underpowered for modest hazard ratios. The negative result on H1 is a statement about this cohort, not about the biology.
- Single cohort. The Zurich cohort (about 70 patients) is in the same archive and is the natural external replication. It was not analysed within the time available.
- Cell types are inherited, not re-derived. We use the authors' PhenoGraph assignments and metacluster map, a deliberate choice that makes the ARI comparison fair but means the pipeline has not been tested end to end from raw marker intensities.
- Breast cancer only, one platform (IMC on a TMA). Generalisation to spot-based spatial transcriptomics is architecturally supported but untested.

## 19. Citation, license, and contributors

If you use Locale, please cite the dataset it is validated on: Jackson, Fischer et al., "The single-cell pathology landscape of breast cancer," Nature 578, 615 to 620 (2020).

Licensed under the MIT License; see [LICENSE](LICENSE).

Contributors: Sahiel Bose, Pranav Achar, Shanay Gaitonde.

Every team will show you a tool that found something. This is the one that knows when it has not.
