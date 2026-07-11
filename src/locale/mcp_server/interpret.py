"""Niche auto-naming (Lane B). Server-side sampling via the Anthropic API.

Turns an opaque niche (composition + marker program) into a human-readable name
like "immune-excluded fibrotic tumor core". If no ANTHROPIC_API_KEY is set (or the
anthropic package is missing), falls back to a DETERMINISTIC rule-based name so the
server works keyless during the hackathon.
"""

from __future__ import annotations

import os

MODEL = "claude-sonnet-5"
_MAX_TOKENS = 60

# Rough cell-type buckets used by the deterministic fallback namer.
_IMMUNE = {"CD8_T", "CD4_T", "B_cell", "Macrophage", "NK", "Treg", "DC"}
_STROMAL = {"Fibroblast", "CAF", "Endothelial", "Stroma"}
_TUMOR = {"Tumor", "Epithelial", "Malignant"}


def _dominant(composition: dict[str, float], k: int = 2) -> list[str]:
    return [
        ct
        for ct, _ in sorted(composition.items(), key=lambda kv: kv[1], reverse=True)[:k]
    ]


def _fallback_name(composition: dict[str, float], marker_program: list[str]) -> str:
    """Deterministic placeholder name derived from composition + markers.

    No network, no key. Good enough to demo; the Anthropic path replaces it when a
    key is present.
    """
    if not composition:
        return "unlabeled niche"
    top = _dominant(composition, 2)
    frac = composition
    immune = sum(frac.get(c, 0.0) for c in _IMMUNE)
    tumor = sum(frac.get(c, 0.0) for c in _TUMOR)
    stromal = sum(frac.get(c, 0.0) for c in _STROMAL)

    if tumor >= 0.5 and immune < 0.15:
        base = "immune-excluded tumor core"
    elif immune >= 0.5:
        base = "immune-infiltrated niche"
    elif stromal >= 0.5:
        base = "stromal / vascular niche"
    else:
        base = f"{top[0].lower()}-enriched niche"

    tag = f" ({', '.join(marker_program[:2])})" if marker_program else ""
    return base + tag


def name_niche(
    composition: dict[str, float],
    marker_program: list[str],
    niche_id: int | None = None,
) -> str:
    """Return a short human-readable name for a niche.

    Uses the Anthropic API when ANTHROPIC_API_KEY is set; otherwise a deterministic
    rule-based fallback. Any API error also falls back, so callers never crash.
    """
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return _fallback_name(composition, marker_program)

    try:
        import anthropic

        comp = ", ".join(f"{ct} {frac:.0%}" for ct, frac in _dominant_full(composition))
        markers = ", ".join(marker_program[:6]) or "n/a"
        prompt = (
            "You are labeling a spatial cellular niche from breast-cancer imaging "
            "mass cytometry. Given its cell-type composition and top enriched protein "
            "markers, return ONLY a concise (<= 6 words) descriptive name, no quotes.\n\n"
            f"Composition: {comp}\n"
            f"Top markers: {markers}\n\n"
            "Name:"
        )
        client = anthropic.Anthropic(api_key=api_key)
        msg = client.messages.create(
            model=MODEL,
            max_tokens=_MAX_TOKENS,
            # Naming a niche is a short deterministic label task; disable extended
            # thinking so it does not eat the small max_tokens budget and leave the
            # response with no text block (which would silently force the fallback).
            thinking={"type": "disabled"},
            messages=[{"role": "user", "content": prompt}],
        )
        text = "".join(
            block.text for block in msg.content if block.type == "text"
        ).strip()
        return text or _fallback_name(composition, marker_program)
    except Exception:
        # Keyless-safe: any failure (no package, bad key, network) => deterministic name.
        return _fallback_name(composition, marker_program)


def _dominant_full(composition: dict[str, float]) -> list[tuple[str, float]]:
    return sorted(composition.items(), key=lambda kv: kv[1], reverse=True)
