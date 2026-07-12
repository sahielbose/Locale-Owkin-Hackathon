"""Niche auto-naming (Lane B). Server-side sampling via the Anthropic API.

Turns an opaque niche (composition + marker program) into a short human-readable
name plus a one-line description, e.g.
    {"name": "immune-excluded fibrotic tumor core",
     "description": "Tumor and fibroblast dense region with almost no cytotoxic T cells."}

Robustness (naming must NEVER break a tool):
  - API key read from ANTHROPIC_API_KEY; if missing, a deterministic fallback runs.
  - Request uses a timeout and thinking disabled (short label task).
  - Any failure (no key, no package, network, bad JSON) returns the deterministic
    fallback derived from the dominant cell types and markers.
"""

from __future__ import annotations

import json
import os
import re

MODEL = "claude-sonnet-5"
_MAX_TOKENS = 200
_TIMEOUT_S = 12.0

# Rough cell-type buckets used by the deterministic fallback namer.
_IMMUNE = {"CD8_T", "CD4_T", "B_cell", "Macrophage", "NK", "Treg", "DC"}
_STROMAL = {"Fibroblast", "CAF", "Endothelial", "Stroma"}
_TUMOR = {"Tumor", "Epithelial", "Malignant"}


def _dominant(composition: dict[str, float], k: int = 2) -> list[str]:
    return [ct for ct, _ in _sorted(composition)[:k]]


def _sorted(composition: dict[str, float]) -> list[tuple[str, float]]:
    return sorted(composition.items(), key=lambda kv: kv[1], reverse=True)


def _fallback_name(composition: dict[str, float], marker_program: list[str]) -> str:
    """Deterministic placeholder name derived from composition + markers. No network."""
    if not composition:
        return "unlabeled niche"
    top = _dominant(composition, 2)
    frac = composition
    immune = sum(frac.get(c, 0.0) for c in _IMMUNE)
    tumor = sum(frac.get(c, 0.0) for c in _TUMOR)
    stromal = sum(frac.get(c, 0.0) for c in _STROMAL)

    if tumor >= 0.5 and immune < 0.15:
        # A cell-level composition cannot distinguish "immune-excluded" (immune present
        # but locked out of the tumor bed) from "immune desert" (no immune at all); on
        # Basel this phenotype is a continuum. So the honest, checkable label describes
        # what we measured, not a mechanism we cannot see. K Pro reasons from this string.
        base = "tumor-rich, immune-poor"
    elif immune >= 0.5:
        base = "immune-infiltrated niche"
    elif stromal >= 0.5:
        base = "stromal / vascular niche"
    else:
        base = f"{top[0].lower()}-enriched niche"

    tag = f" ({', '.join(marker_program[:2])})" if marker_program else ""
    return base + tag


def _fallback_description(
    composition: dict[str, float], marker_program: list[str]
) -> str:
    """One-line deterministic description of the niche composition + markers."""
    if not composition:
        return "No cells assigned to this niche."
    parts = [f"{ct} {frac:.0%}" for ct, frac in _sorted(composition)[:3]]
    markers = ", ".join(marker_program[:3])
    tail = f"; top markers {markers}" if markers else ""
    return "Dominated by " + ", ".join(parts) + tail + "."


def name_and_describe_niche(
    composition: dict[str, float],
    marker_program: list[str],
    niche_id: int | None = None,
) -> dict[str, str]:
    """Return {"name", "description"} for a niche. Never raises.

    Uses the Anthropic API when ANTHROPIC_API_KEY is set; otherwise (or on any
    failure) returns the deterministic fallback.
    """
    fallback = {
        "name": _fallback_name(composition, marker_program),
        "description": _fallback_description(composition, marker_program),
    }

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return fallback

    try:
        import anthropic

        comp = ", ".join(f"{ct} {frac:.0%}" for ct, frac in _sorted(composition))
        markers = ", ".join(marker_program[:6]) or "n/a"
        prompt = (
            "You are labeling a spatial cellular niche from breast-cancer imaging "
            "mass cytometry. Given its cell-type composition and top enriched protein "
            "markers, respond with ONLY a JSON object of the form "
            '{"name": "...", "description": "..."} where name is at most 6 words and '
            "description is one sentence. Do not use em-dashes.\n\n"
            f"Composition: {comp}\n"
            f"Top markers: {markers}\n"
        )
        client = anthropic.Anthropic(api_key=api_key, timeout=_TIMEOUT_S)
        msg = client.messages.create(
            model=MODEL,
            max_tokens=_MAX_TOKENS,
            # Short deterministic label task: disable extended thinking so it does not
            # eat the max_tokens budget and leave no text block (which would force the
            # fallback even with a valid key).
            thinking={"type": "disabled"},
            messages=[{"role": "user", "content": prompt}],
        )
        text = "".join(b.text for b in msg.content if b.type == "text").strip()
        data = _parse_json(text)
        name = str(data.get("name") or "").strip()
        description = str(data.get("description") or "").strip()
        return {
            "name": name or fallback["name"],
            "description": description or fallback["description"],
        }
    except Exception:
        return fallback


def name_niche(
    composition: dict[str, float],
    marker_program: list[str],
    niche_id: int | None = None,
) -> str:
    """Return just the short name (used for the Niche.name schema field)."""
    return name_and_describe_niche(composition, marker_program, niche_id)["name"]


def _parse_json(text: str) -> dict:
    """Parse a JSON object from model text, tolerating surrounding prose/fences."""
    try:
        return json.loads(text)
    except Exception:
        match = re.search(r"\{.*\}", text, re.DOTALL)
        if match:
            return json.loads(match.group(0))
        raise ValueError("no JSON object in model response")
