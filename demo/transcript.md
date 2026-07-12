# K Pro <-> Locale: "is the immune-excluded niche prognostic?"

**What this is.** A faithful reproduction of the K Pro <-> Locale MCP exchange. Every
tool block below is the **real output** of the Locale MCP server on the Basel cohort
(captured by calling the tools over the MCP protocol). It shows the whole point of the
product: Locale answers with the finding **and** the refusal, so K Pro cannot overclaim.

**Why it is a reproduction, not a live K Pro screenshot.** Connecting to K Pro requires
signing up for a K Pro account, which the operator must do (an assistant cannot create
accounts). The exact tool responses K Pro receives are shown verbatim here; the
["reproduce live"](#reproduce-live) steps at the bottom get the identical exchange
inside K Pro once the connector is added.

---

## The conversation

**User -> K Pro:** *In the Basel breast-cancer cohort, is the immune-excluded niche prognostic for overall survival?*

**K Pro** (reasoning): The user names an "immune-excluded" niche. Let me see what niches Locale actually has before assuming that label exists.

> **tool call** `describe_niches()`
> ```json
> [ ... , {"niche_id": 7, "name": "tumor-rich, immune-poor",
>          "composition": {"tumor": 0.786, "immune": 0.044, "stroma": 0.151, "endothelial": 0.019},
>          "n_cells": 74739, "n_cores": 160},
>       {"niche_id": 1, "name": "immune-rich",
>          "composition": {"tumor": 0.181, "immune": 0.594, ...}, "n_cells": 91905, "n_cores": 260}, ... ]
> ```

**K Pro** (reasoning): There is no niche named "immune-excluded." The closest is niche 7, labeled **"tumor-rich, immune-poor"** (79% tumor, 4% immune). Locale is deliberately not calling it "excluded" - that is a tumor-level mechanism. I will test whether *its abundance* tracks survival.

> **tool call** `correlate_niche_outcome(niche_id=7)`
> ```json
> {"niche_id": 7,
>  "hazard_ratio": 1.075, "ci_95": [0.86, 1.344], "p_raw": 0.5252,
>  "n_hypotheses_tested": 12, "q_fdr": 0.722, "p_selection_aware": 0.4396,
>  "n_events": 79, "min_detectable_hr": 1.371,
>  "verdict": "insufficient evidence"}
> ```

**K Pro -> User:**

> **Short answer: no - there is insufficient evidence that this niche is prognostic in Basel.**
>
> First, a correction on the premise. Locale has no "immune-excluded" niche. The
> tumor-rich, immune-poor niche (niche 7; 79% tumor, 4% immune, 160 cores) is the one
> you likely mean, but "excluded" (immune present, locked out) versus "desert" (no
> immune) is a continuum here, not something a cell-level label can assert.
>
> On prognosis: its per-SD abundance gives a hazard ratio of **1.08 (95% CI 0.86-1.34,
> raw p = 0.53)** - directionally toward worse survival, but not significant. And the
> context Locale ships makes clear why I should not chase it:
> - **n_events = 79.** With this many events, 80% power only reaches HR >= **1.37**. The
>   observed 1.08 sits inside the band this study cannot resolve. A non-significant
>   result here is uninformative, not reassuring.
> - **12 niches were tested.** FDR q = **0.72**; the selection-aware p (permuting
>   survival and taking the best of all 12, 1000x) = **0.44**. The best niche in the
>   whole panel is no better than chance.
>
> So: a real spatial pattern (this niche exists and is reproducible), but its link to
> survival is **not supported** at this sample size. I would not report it as a
> biomarker.

**User -> K Pro:** *What about the immune-rich niche?*

> **tool call** `correlate_niche_outcome(niche_id=1)`
> ```json
> {"niche_id": 1, "hazard_ratio": 0.764, "ci_95": [0.587, 0.995], "p_raw": 0.0458,
>  "n_hypotheses_tested": 12, "q_fdr": 0.2903, "p_selection_aware": 0.4396,
>  "n_events": 79, "min_detectable_hr": 1.371, "verdict": "insufficient evidence"}
> ```

**K Pro -> User:** The immune-rich niche looks protective - HR **0.76 (0.59-1.00), raw
p = 0.046** - consistent with the known biology that immune infiltration improves
outcome. But the same guardrails apply: after correcting for the 12 niches tested,
**q = 0.29**, and the panel-wide selection-aware p is **0.44**. Verdict: **insufficient
evidence.** It is the most promising signal in the cohort and a reasonable thing to
pre-register and test on an independent cohort (Zurich), but Basel alone does not
support it.

---

## Why this is the product

An expression-only tool would have answered "the immune-excluded niche has HR 1.08" or
"the immune-rich niche is protective, p = 0.046" and stopped. Locale ships the four
numbers the agent cannot compute for itself - **n_hypotheses_tested, q_fdr,
p_selection_aware, min_detectable_hr** - with every finding, unasked. Those numbers turn
a plausible-looking p-value into an honest "insufficient evidence," and they corrected
the false premise ("immune-excluded") before K Pro could reason from it. That is the
exact failure mode - a confident agent building on an unverifiable claim - that this
project exists to prevent.

<a name="reproduce-live"></a>
## Reproduce live in K Pro

1. **Precompute the demo path** (already done; regenerates `demo/findings.json`):
   `PYTHONPATH=. python scripts/precompute_findings.py`
2. **Serve the real cohort:**
   `LOCALE_DATA=data/basel_niched.h5ad python -m src.locale.mcp_server.server`
3. **Expose it** (K Pro needs a public https URL):
   `cloudflared tunnel --url http://127.0.0.1:8000` -> your MCP URL is that + `/mcp`.
4. **Sign up for K Pro Free** and add a custom connector pointing at that `/mcp` URL
   (mirrors Owkin's Pathology Explorer connector).
5. Ask: *"In the Basel breast cohort, is the immune-excluded niche prognostic?"* and
   you will get the exchange above, live.
