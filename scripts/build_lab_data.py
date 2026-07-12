"""Build viz/app/lab_data.js for the engine-and-tests frontend (Lane C).

Bundles two real artifacts into window.LOCALE_LAB:
- the engine's findings (demo/findings.json): Basel-cohort niches + survival stats;
- the test-suite result: run `pytest --json-report` and pass its json path here.

    python -m pytest --json-report --json-report-file=report.json
    python scripts/build_lab_data.py report.json

If no report path is given it looks for report.json in the repo root; the findings
are always read from demo/findings.json.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "src" / "localespatial" / "viz" / "app" / "lab_data.js"


def _engine_block() -> dict | None:
    """Compute enrichment, KM survival, and validation live on the committed cohort.

    These are the graphs the Basel findings.json doesn't carry. Runs the real engine on
    data/mock.h5ad (the one object available offline). Returns None if the heavy engine
    deps (squidpy/lifelines/scanpy) aren't installed, so the core build never breaks.
    """
    import os
    import warnings

    warnings.filterwarnings("ignore")
    try:
        import anndata as ad

        os.environ["LOCALE_DATA"] = str((ROOT / "data" / "mock.h5ad").resolve())
        sys.path.insert(0, str(ROOT))
        from src.localespatial.engine import niches as NI
        from src.localespatial.engine import validate as VAL
        from src.localespatial.mcp_server import tools

        tools.reset_cache()
        enr = tools.compute_enrichment()
        prog = tools.find_prognostic_niches()
        survival = [
            {
                "name": n.name,
                "hazard_ratio": n.prognostic.hazard_ratio,
                "ci_low": n.prognostic.ci_low,
                "ci_high": n.prognostic.ci_high,
                "pvalue": n.prognostic.pvalue,
                "n_patients": n.prognostic.n_patients,
                "km": {
                    "time": n.prognostic.km.time,
                    "high": n.prognostic.km.high,
                    "low": n.prognostic.km.low,
                },
            }
            for n in prog
            if n.prognostic and n.prognostic.km
        ]
        adata = ad.read_h5ad(ROOT / "data" / "mock.h5ad")
        lab = NI.find_niches(adata, n_niches=6)
        ids = sorted({int(x) for x in lab.obs["niche"].to_numpy()})
        shuffle = VAL.shuffle_negative_control(lab, n_permutations=80)
        stability = VAL.stability_ari(lab, n_runs=8)
        markers = [VAL.marker_validation(lab, i) for i in ids]
        return {
            "cohort_label": "computed live on the committed cohort (data/mock.h5ad)",
            "enrichment": {
                "cell_types": enr.cell_types,
                "zscores": enr.zscores,
                "pvalues": enr.pvalues,
            },
            "survival": survival,
            "validation": {
                "shuffle": {
                    k: float(v) if not isinstance(v, bool) else v
                    for k, v in shuffle.items()
                },
                "stability_ari": float(stability),
                "markers": [
                    {
                        "niche_id": m["niche_id"],
                        "passed": bool(m["passed"]),
                        "corr": (
                            float(m["evidence"].get("marker_composition_corr"))
                            if isinstance(m.get("evidence"), dict)
                            and m["evidence"].get("marker_composition_corr") is not None
                            else None
                        ),
                        "top": (
                            [str(x) for x in m["evidence"].get("top_markers", [])]
                            if isinstance(m.get("evidence"), dict)
                            else []
                        ),
                    }
                    for m in markers
                ],
            },
        }
    except Exception as exc:  # noqa: BLE001
        print(f"(engine block skipped: {exc})", file=sys.stderr)
        return None


def _humanize(name: str) -> str:
    n = re.sub(r"^test_", "", name).replace("_", " ")
    return n[:1].upper() + n[1:]


def main() -> None:
    findings = json.loads((ROOT / "demo" / "findings.json").read_text())

    report_path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "report.json"
    tests: list[dict] = []
    meta: dict = {
        "total": 0,
        "passed": 0,
        "failed": 0,
        "skipped": 0,
        "duration_s": 0,
        "files": [],
    }
    if report_path.exists():
        rep = json.loads(report_path.read_text())
        for t in rep.get("tests", []):
            file, name = t["nodeid"].split("::", 1)
            dur = (t.get("call") or {}).get("duration") or 0
            tests.append(
                {
                    "file": file.replace("tests/", ""),
                    "name": name,
                    "label": _humanize(name),
                    "outcome": t["outcome"],
                    "ms": round(dur * 1000, 1),
                }
            )
        s = rep.get("summary", {})
        meta = {
            "total": s.get("total", 0),
            "passed": s.get("passed", 0),
            "failed": s.get("failed", 0),
            "skipped": s.get("skipped", 0),
            "duration_s": round(rep.get("duration", 0), 2),
            "files": sorted({t["file"] for t in tests}),
        }
    else:
        print(f"(no {report_path.name}; tests tab will be empty)", file=sys.stderr)

    lab = {
        "findings": findings,
        "tests": tests,
        "meta": meta,
        "engine": _engine_block(),
    }
    OUT.write_text(
        "// Auto-generated by scripts/build_lab_data.py: real engine findings "
        "(demo/findings.json) + pytest results.\n"
        "window.LOCALE_LAB = " + json.dumps(lab) + ";\n"
    )
    print(
        f"wrote {OUT}: {len(findings['niches'])} niches, "
        f"{meta['passed']}/{meta['total']} tests passing"
    )


if __name__ == "__main__":
    main()
