# SPATIAL NICHE EXPLORER
### An MCP for K Pro that reasons about *where* cells sit in tumor tissue, the cellular neighborhoods and immune niches that predict outcomes
**"Not what's expressed, but who's next to whom, and whether that geography decides how the patient does."**

**Track:** Best AI Scientist MCP (Rewiring Biology Hackathon, Owkin)
**Team shape:** agentic-workflow builders + at least one computational biologist (this idea is bio-load-bearing, see §11)
*Name is a placeholder, alternatives: **Niche**, **Terrain**, **Neighborhood**.*

---

## 1. The one-paragraph pitch

Every tool in K Pro today treats a tumor sample as a bag of cells: it tells you *what* genes are on and *how much*. But cancer outcomes are often decided by tissue **geography**: whether cytotoxic T-cells are physically infiltrating the tumor or locked out at the margin, whether cancer-associated fibroblasts form a barrier, whether immune cells organize into tertiary lymphoid structures. Spatial Niche Explorer is an MCP server that lets K Pro reason about that geography. It takes MOSAIC's spatial data, builds the spatial graph of a tumor, and answers questions no expression-only tool can: which cell types co-localize more than chance, what recurring **cellular niches** exist across a cohort, what molecular programs define them, and, the payoff, which niches track patient survival. It's the difference between "GENE X is high in this tumor" and "GENE X marks an immune-excluded niche that predicts poor response to checkpoint blockade."

## 2. What it actually does (the capability ladder)

Four capabilities, each strictly deeper than K Pro's current expression-only view:

1. **Co-location**: for any pair of cell types, are they neighbors more often than chance? (e.g. "are CD8 T-cells adjacent to malignant cells, or avoiding them?")
2. **Niche discovery**: find the recurring cellular neighborhoods across a cohort (e.g. an "immune-excluded tumor core," a "TLS-like immune aggregate," a "CAF barrier interface").
3. **Niche characterization**: what cell composition and gene/pathway programs define each niche.
4. **Niche → outcome**: does a niche's abundance track survival or treatment response? This is where it becomes a *biomarker engine*, not just a description.

## 3. Why this wins Track 1 (mapped to the six official criteria)

| Criterion | How Spatial Niche Explorer scores |
|---|---|
| **Integrability** | Standard remote MCP server; connects like Owkin's own Pathology Explorer connector. Tools are independently callable and composable into K Pro skills. |
| **Performance** | Produces a real, quantitative spatial finding, a prognostic niche with a survival curve and a rendered map, not a qualitative gesture. |
| **Depth of biomedical reasoning** | Spatial neighborhood + niche analysis is genuinely advanced tumor-microenvironment biology; reasoning about *tissue architecture* is deeper than differential expression. |
| **Extension of K Pro** | This is the strongest axis: K Pro's Multi-Omics tool is expression/survival-oriented and treats samples non-spatially. Spatial reasoning is a capability it *entirely lacks*, and it's MOSAIC's whole reason to exist. |
| **Usefulness & reusability** | Works for any cohort/cell-type/niche question; spatial biomarkers are exactly what pharma wants for patient stratification. |
| **Technical quality & UX** | The interactive tissue map (MCP-app UI) is a visual centerpiece no text tool can match, the single most "wow" artifact at the event. |

**Themes hit:** *multimodal reasoning over spatial + genomic data* (dead center), and it touches *AI-driven hypothesis generation* (niches suggest testable mechanisms).

**Strategic read:** this is the **highest-ceiling** idea because it leans directly on MOSAIC's differentiator. But it is also **higher-variance**: its whole feasibility rests on one data-access question (§10). Choose it only with a bio person on the team and a plan to resolve that question in hour 0.

## 4. The analysis pipeline (the technical core)

The methods here are standard and well-supported in the `scverse` ecosystem, you are orchestrating known tools, not inventing algorithms. That's a strength: it's reliable and defensible.

**Step 0: Get cell-type-resolved spatial data.** This is the fork that determines everything (see §10):
- *If true single-cell spatial (e.g. CosMx):* you have per-cell x/y coordinates + cell-type labels directly. Proceed to Step 1.
- *If Visium (spot-based, MOSAIC's primary spatial platform):* each ~55µm spot contains several cells, so you first **deconvolve** each spot into cell-type proportions using the paired single-cell reference, `cell2location`, `RCTD`, or `Tangram`. Your "cells" become spot-level compositions. Coarser, but valid and publishable.

**Step 1: Build the spatial graph.** `squidpy.gr.spatial_neighbors()`, a k-nearest-neighbor, Delaunay, or fixed-radius graph over cell/spot coordinates.

**Step 2: Co-location (neighborhood enrichment).** `squidpy.gr.nhood_enrichment()`, a permutation test asking whether each cell-type pair is spatially adjacent more/less than chance. Output: a cell-type × cell-type z-score matrix (the "who's next to whom" heatmap).

**Step 3: Niche discovery.** Identify recurring cellular neighborhoods across the cohort. Best fit for multi-sample data: **CellCharter** (aggregates each cell's l-hop neighborhood features, clusters with a GMM, and is explicitly built to compare niches across many samples and technologies). Alternatives: classic windowed k-means on local cell-type composition (the Schürch/Nolan "cellular neighborhoods" method), `scNiche`, `UTAG`, or `BANKSY`.

**Step 4: Characterize niches.** For each niche: its cell-type composition, and the enriched genes/pathways (spatially variable genes via `squidpy` Moran's I; pathway enrichment on the niche's expression).

**Step 5: Niche → outcome.** Correlate each niche's per-patient abundance with survival/response (Cox model / Kaplan-Meier), using MOSAIC Window's clinical annotations. Rank niches by prognostic strength.

## 5. The biology it unlocks (the "why it matters")

This is what makes the demo land with a science-literate judge. Spatial organization of the tumor microenvironment is repeatedly, clinically decisive:
- **Immune exclusion vs. infiltration.** Whether T-cells penetrate the tumor or are stuck at the stromal margin is roughly how checkpoint-blockade response is judged, an "immune-excluded" niche is a poor-prognosis signature bulk data can't see.
- **Tertiary lymphoid structures (TLS).** Organized immune aggregates in tumors are associated with *better* immunotherapy response, a niche worth detecting.
- **CAF barriers & tumor-immune interfaces.** Cancer-associated fibroblast niches can physically wall off tumors from immune attack.
- **Documented prognostic niches** (precedent for your demo): a CAF+endothelial neighborhood tracked recurrence in salivary duct carcinoma; low stromal-TIL and low CD4-stromal interaction predicted recurrence in colorectal cancer; CellCharter found a neutrophil + hypoxic-tumor niche in NSCLC. Your tool reproduces this *class* of finding on demand.

## 6. Frontier MCP features (they said these impress them)

- **MCP-app UI rendering** → the **interactive tissue map** is the star: plot cells/spots colored by cell type, by niche, or by a gene, highlight the prognostic niche, optionally overlay on the H&E image. This is the visual nobody else will have.
- **Server-side sampling** → auto-interpret/name each niche: an LLM step turns "niche 4 = 60% malignant, 25% CAF, 3% CD8, TGF-β program up" into "immune-excluded fibrotic tumor core." Turns opaque cluster IDs into biology.
- **Elicitation** → underspecified request ("find interesting niches", in which cohort? which cell types matter? niche resolution?) → prompt for structured input before running.

## 7. MCP tool surface (the product)

**Orchestrator**
- `find_prognostic_niches(cohort, patient_subset?, n_niches?) -> RankedNiches`, runs the full pipeline (graph → niches → characterization → outcome), returns niches ranked by prognostic strength + the tissue map. Elicits missing inputs.

**Spatial analysis (also standalone → inspectable)**
- `list_cohorts()` / `list_cell_types(cohort)`
- `describe_spatial_sample(cohort, patient)`, platform, #cells/spots, resolution, what's available (so the agent knows what it's working with)
- `compute_neighborhood_enrichment(cohort, patient?)`, cell-type co-location matrix (Step 2)
- `identify_niches(cohort, n_niches?)`, recurring cellular niches (Step 3)
- `characterize_niche(cohort, niche_id)`, composition + gene/pathway programs (Step 4)
- `correlate_niche_outcome(cohort, niche_id)`, survival/response association (Step 5)
- `compare_niches_across_cohorts(niche_id)`, is this niche conserved across cancer types (CellCharter's cross-sample strength)
- `deconvolve_spots(cohort, patient)`, Visium → per-spot cell-type proportions (utility; may be pre-computed)

**Presentation**
- `render_spatial_map(cohort, patient, color_by="cell_type"|"niche"|"gene")`, interactive tissue map (MCP-app)
- `export_niche_report(cohort, niche_id, format)`, findings + figures

**Coordination note:** the whole idea hinges on the data returned by `describe_spatial_sample` / `deconvolve_spots`. Get those returning real data in hour 0 before building anything downstream.

## 8. How it plugs into K Pro

Remote MCP server (mirror the Pathology Explorer connector pattern: custom connector URL in K Pro/Claude settings). K Pro's orchestrator calls `find_prognostic_niches` as a skill; the tissue map renders in-session; niche findings feed downstream reasoning. Data comes from the provided MOSAIC Window tier. **The "extends K Pro" story is airtight:** K Pro's current Multi-Omics tool analyzes expression and survival but is spatially blind, it cannot tell you what's adjacent to what. Spatial Niche Explorer adds the entire spatial-reasoning dimension of MOSAIC that K Pro doesn't currently expose.

## 9. The killer demo (this is what wins)

**Principle:** produce a complete spatial-biology finding, live, that an expression-only tool physically cannot.

**The demo query:** *"In the MOSAIC bladder cohort, find the tissue niche most associated with poor survival, tell me what defines it, and show me where it is."*

**What K Pro returns:**
1. A ranked list of niches → the top one is, say, an **immune-excluded fibrotic niche** (malignant cells + CAFs, almost no CD8 T-cells).
2. Its defining program (T-cell exclusion / TGF-β signature).
3. A **survival curve**: patients with high abundance of this niche do significantly worse.
4. A **rendered tissue map** of a patient, cells colored by niche, the poor-prognosis niche highlighted, you literally point at the immune desert.

That sequence, discovery → mechanism → clinical association → *visual proof*, is a stronger single artifact than anything a literature-plus-expression tool can produce. It is exactly why MOSAIC's spatial resolution exists.

**Backup demo:** find and show a *good*-prognosis niche (a TLS-like immune aggregate) in a different cohort, demonstrates the tool generalizes and isn't a one-trick result.

## 10. The hardcore bottleneck (read this before committing)

**Single point of failure: what spatial granularity does MOSAIC Window actually expose?** This is not a detail, it's the whole bet.

The facts: MOSAIC's primary spatial platform is **10x Visium, which is spot-based, not single-cell**: each spot (~55µm) captures several cells, ~1,700 spots per sample, up to ~29,000 genes. MOSAIC also has paired **10x Chromium** single-cell (dissociated, no coordinates) and, for some samples, possibly **NanoString CosMx** (true single-cell spatial). The free MOSAIC Window tier's exact contents/granularity must be confirmed.

The three cases:
- **Best case: CosMx / true single-cell spatial is in MOSAIC Window** → per-cell coordinates + labels → run the pipeline directly. Ideal.
- **Likely case: Visium spots with coordinates** → you must add a **deconvolution** step (cell2location / RCTD / Tangram) to get per-spot cell-type composition, then do spot-level niche analysis. Fully valid, but coarser and adds a real chunk of work; "cell A touches cell B" becomes "spot enriched in A is adjacent to spot enriched in B."
- **Worst case: only aggregated/processed spatial outputs, no spot coordinates** → the core idea is blocked.

**Hour-0 test (do this before writing any pipeline code):** call for one MOSAIC Window sample and check, do you get (a) per-cell or per-spot x/y coordinates, and (b) cell-type labels, or expression you can deconvolve against the single-cell reference? The answer decides the idea's fate immediately.

**Fallback ladder (honest):**
1. True single-cell niche analysis (CosMx), best.
2. Visium + deconvolution → spot-niche analysis, solid, coarser.
3. H&E-morphology spatial features via a pathology model (ties to Pathology Explorer's cell-type maps), a last-resort, weaker version if transcriptomic spatial data is unusable.
4. Below that, no viable version, unlike Consilience, this idea does *not* fully survive a total spatial-data failure.

That shorter, front-loaded fallback ladder is the core tradeoff versus Consilience: **higher ceiling, thinner safety net.**

## 11. Team split (this idea is bio-load-bearing)

- **Computational biologist (required, ideally spatial/single-cell experience):** owns Steps 0-5, data access, deconvolution if needed, squidpy graph + enrichment, CellCharter niches, survival association, and (critically) sanity-checking that the niches are biologically real, not artifacts. This is the critical path.
- **Agentic engineer(s) / you:** the MCP server scaffold, orchestrator, elicitation, server-side niche auto-naming, and the interactive tissue-map UI (the visual centerpiece).
- **Floating:** pathway/enrichment characterization, demo polish, pitch, backup video.

If you have no bio person, do **not** pick this idea, the analysis validity and the hour-0 data triage both need domain judgment. (This is the main reason Consilience is the safer default for an agentic-heavy team.)

## 12. 48-hour plan

| Window | Goal |
|---|---|
| **0-3h** | **Resolve the bottleneck first.** Confirm MOSAIC Window spatial granularity (single-cell vs Visium spots vs aggregated). Decide the pipeline branch. Scaffold MCP server + connector; confirm it appears in Claude. |
| **3-8h** | Data loading + (if Visium) deconvolution working on ONE cohort (bladder or lung, most patients). Spatial graph + neighborhood enrichment returning real numbers. |
| **8-20h** | Niche discovery (CellCharter) + characterization. Lock the demo cohort and confirm a real prognostic niche exists. |
| **20-30h** | Niche→survival association. Build the interactive tissue map (the centerpiece). Niche auto-naming via server-side sampling. |
| **30-40h** | Orchestrator (`find_prognostic_niches`) end-to-end. Elicitation. Backup cohort/niche. Export. **Feature freeze.** |
| **40-48h** | Rehearse pitch. Buffer. **Record a backup demo video** against a known-good niche so a live failure can't sink you. |

## 13. Risk management

- **Resolve the data-granularity bottleneck in hours 0-3.** Everything depends on it; find out before you build.
- **Deconvolution is the likely hidden cost**: budget for it if the data is Visium, and consider whether MOSAIC Window ships pre-computed cell-type annotations you can use instead.
- **Validate niches are real, not clustering artifacts**: this is where the bio person earns their seat; a nonsense niche in the demo is fatal.
- **Backup demo video** against a confirmed prognostic niche.
- Unlike Consilience, there is no full graceful-degradation path, the mitigation is front-loading the data check, not a fallback architecture.

## 14. Pre-event checklist

- [ ] **Confirm (or plan to confirm hour-0) MOSAIC Window spatial granularity**: single-cell, Visium spots, or aggregated. This is the single most important prep item.
- [ ] Line up a **computational biologist** on the team, non-negotiable for this idea.
- [ ] Get `squidpy` + `CellCharter` (+ a deconvolution tool like `cell2location`) installed and run their tutorials so the pipeline is muscle memory.
- [ ] Get a **hello-world MCP server** connecting to Claude as a custom connector (mirror Pathology Explorer).
- [ ] Pre-read a couple of TME spatial papers (immune exclusion, TLS, prognostic niches) so you can frame the demo finding in real biology.
- [ ] Pick the demo cohort in advance (bladder / lung = most patients in MOSAIC Window).

---

*One-line version for the group chat:* **Spatial Niche Explorer** = an MCP that lets K Pro reason about *where* cells sit in MOSAIC tumor tissue, finding the cellular niches (like immune-excluded tumor cores) that predict survival and rendering them on an interactive tissue map, a spatial-reasoning capability K Pro completely lacks. Highest-ceiling idea, but the whole thing rests on one hour-0 question: does MOSAIC Window expose cell/spot-level spatial coordinates? Needs a bio person; thinner safety net than Consilience.
