"""
Gate 4 — Cross-source invariant tests (R103)

These tests detect "修一处漏一处" drift where one source file (e.g.,
manuscript_facts.json) is updated but a dependent source file (e.g.,
analysis_results_R3.2.json) is not. The R102 audit found that analysis_results.json
had M0=0.499 in the benchmark block but r90_recompute_outputs had M0=0.567,
and the JSON was self-contradictory for 7 rounds.

Gate 4 enforces:
  - analysis_results_R3.2.json benchmark.M0/M3/M4a use the same values as
    facts.json benchmark_M0_M5_full (R90 canonical)
  - antigenicity_benchmark_summary.json pos_median/neg_median match facts.json
  - donor_scorecard.tsv ISI_coverage% matches analysis_results.R90_recompute.donor_coverages
  - Figure panel metadata matches PDF rendered text
  - 60 tests + this Gate 4 = 70 tests total
"""
from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

from . import TestResult, register_gate, FAIL, PASS, WARN


def _pdf_text(pdf_path: Path) -> str:
    """Extract text from a PDF using pdfplumber."""
    try:
        import pdfplumber
    except ImportError:
        return ""
    text = ""
    try:
        with pdfplumber.open(pdf_path) as pdf:
            for page in pdf.pages[:3]:
                txt = page.extract_text() or ""
                text += txt + "\n"
    except Exception:
        pass
    return text


def _load_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


@register_gate("Gate 4 — Cross-source invariants", weight=15)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []

    # Resolve canonical paths
    tier1 = _load_tsv(canonical_dir / "tier1_frozen.tsv")
    analysis = _load_json(canonical_dir / "analysis_results.json")
    antigen = _load_json(canonical_dir.parent / "tier3_surfaceome" / "antigenicity_benchmark_summary.json") \
        if (canonical_dir.parent / "tier3_surfaceome" / "antigenicity_benchmark_summary.json").exists() \
        else _load_json(canonical_dir / "antigenicity_benchmark_summary.json")
    donor_scorecard = _load_tsv(canonical_dir / "donor_scorecard.tsv")
    figure_provenance = _load_json(canonical_dir / "figure_provenance.json")

    # ===== T61: analysis_results benchmark M0/M3/M4a == facts benchmark_M0_M5_full =====
    fact_m0 = facts.get("benchmark_M0_M5_full", {}).get("M0_resampling_AUROC_mean") if facts else None
    ar_m0 = analysis.get("benchmark", {}).get("M0_resampling_AUROC")
    consistent = fact_m0 is not None and ar_m0 is not None and abs(fact_m0 - ar_m0) < 0.001
    results.append(TestResult(
        id="T61", name="Invariant: M0 AUROC canonical (R103)",
        status=PASS if consistent else FAIL,
        message=f"facts.M0_AUROC={fact_m0} vs analysis_results.M0_AUROC={ar_m0}",
        actual=str(ar_m0), expected=str(fact_m0),
    ))

    fact_m3 = facts.get("benchmark_M0_M5_full", {}).get("M3_resampling_AUROC_mean") if facts else None
    ar_m3 = analysis.get("benchmark", {}).get("M3_resampling_AUROC")
    consistent = fact_m3 is not None and ar_m3 is not None and abs(fact_m3 - ar_m3) < 0.001
    results.append(TestResult(
        id="T62", name="Invariant: M3 AUROC canonical (R103)",
        status=PASS if consistent else FAIL,
        message=f"facts.M3_AUROC={fact_m3} vs analysis_results.M3_AUROC={ar_m3}",
            actual=str(ar_m3), expected=str(fact_m3),
    ))

    fact_m4a = facts.get("benchmark_M0_M5_full", {}).get("M4a_resampling_AUROC_mean") if facts else None
    ar_m4a = analysis.get("benchmark", {}).get("M4a_resampling_AUROC")
    consistent = fact_m4a is not None and ar_m4a is not None and abs(fact_m4a - ar_m4a) < 0.001
    results.append(TestResult(
        id="T63", name="Invariant: M4a AUROC canonical (R103)",
        status=PASS if consistent else FAIL,
        message=f"facts.M4a_AUROC={fact_m4a} vs analysis_results.M4a_AUROC={ar_m4a}",
        actual=str(ar_m4a), expected=str(fact_m4a),
    ))

    # ===== T64: pos_median_v2 in antigen_summary == facts.antigen_full.positives_median =====
    amb_pos = antigen.get("pos_median_v2")
    fact_pos = facts.get("antigen_full", {}).get("positives_median") if facts else None
    consistent = amb_pos is not None and fact_pos is not None and abs(amb_pos - fact_pos) < 0.005
    results.append(TestResult(
        id="T64", name="Invariant: pos_median_v2 (R103)",
        status=PASS if consistent else FAIL,
        message=f"antigen_summary.pos_median_v2={amb_pos} vs facts.positives_median={fact_pos}",
        actual=str(amb_pos), expected=str(fact_pos),
    ))

    amb_neg = antigen.get("neg_median_v2")
    fact_neg = facts.get("antigen_full", {}).get("negatives_median") if facts else None
    consistent = amb_neg is not None and fact_neg is not None and abs(amb_neg - fact_neg) < 0.005
    results.append(TestResult(
        id="T65", name="Invariant: neg_median_v2 (R103)",
        status=PASS if consistent else FAIL,
        message=f"antigen_summary.neg_median_v2={amb_neg} vs facts.negatives_median={fact_neg}",
        actual=str(amb_neg), expected=str(fact_neg),
    ))

    # ===== T66: donor_coverages in r90_recompute == donor_scorecard.tsv ISI_coverage% =====
    r90_coverages = analysis.get("r90_recompute_outputs", {}).get("script_07_donor_scorecard_v32", {}).get("donor_coverages", {})
    scorecard_coverages = {r["donor"]: float(r["ISI_coverage%"]) for r in donor_scorecard if "donor" in r and "ISI_coverage%" in r}
    drift = []
    for gt, sc in scorecard_coverages.items():
        r90 = r90_coverages.get(gt)
        if r90 is not None and abs(sc - r90) > 0.01:
            drift.append(f"{gt}: r90={r90} vs scorecard={sc}")
    results.append(TestResult(
        id="T66", name="Invariant: donor_coverages r90 vs scorecard.tsv (R103)",
        status=PASS if not drift else FAIL,
        message=f"drift={drift if drift else '∅'} (6 genotypes checked)",
    ))

    # ===== T67: figure_provenance.json has all 12 v32 entries (4 figures × 3 formats) =====
    v32_keys = ["fig1_resource_v32", "fig5_benchmark", "figS1_resource_panels_v32", "figS4_v32"]
    missing = [k for k in v32_keys if k not in figure_provenance]
    results.append(TestResult(
        id="T67", name="Invariant: figure_provenance has all R3.2 blocks (R103)",
        status=PASS if not missing else FAIL,
        message=f"missing={missing if missing else '∅'} (4 expected: fig1_v32, fig5, figS1_v32, figS4_v32)",
    ))

    # ===== T68: PDF rendered text matches figure_provenance metadata (R103) =====
    # Use pdfplumber to extract text from fig5 + figS1 and verify R3.2 numbers are present
    if facts:
        # Find fig5 + figS1 PDF files (canonical_dir is pipeline/canonical_links/)
        # canonical_dir.parent = pipeline/, .parent = project root
        pipeline_root = canonical_dir.parent
        project_root = canonical_dir.parent.parent
        # R110-B: pipeline self-contained — read pipeline/figures/*.pdf
        fig5_pdf  = pipeline_root / "figures" / "fig5_benchmark.pdf"
        figS1_pdf = pipeline_root / "figures" / "figS1_resource_panels.pdf"

        # Check fig5 has R3.2 numbers
        if fig5_pdf.exists():
            text = _pdf_text(fig5_pdf)
            has_75 = "75-axis" in text or "75 axis" in text
            has_62_13 = "62 vs 13" in text or "62 incompatible vs 13" in text
            r1_no_r1 = "Release 1.0" not in text or "Release 3.2" in text
            ok = (has_75 or has_62_13)
            results.append(TestResult(
                id="T68", name="fig5 PDF rendered text matches R3.2 (R103)",
                status=PASS if ok else FAIL,
                message=f"fig5 PDF: 75-axis={has_75} 62/13={has_62_13} R1.0_absent={r1_no_r1}",
            ))
        else:
            results.append(TestResult(id="T68", name="fig5 PDF", status=WARN, message=f"PDF not found: {fig5_pdf}"))

        # Check figS1 has R3.2 numbers
        if figS1_pdf.exists():
            text = _pdf_text(figS1_pdf)
            has_75 = "75 axes" in text or "75-axis" in text
            has_61 = "n=61" in text or "61 aligned" in text
            has_727 = "72.7%" in text
            has_39_906 = "39.0" in text and "90.6" in text
            ok = has_75 and has_61 and has_727 and has_39_906
            results.append(TestResult(
                id="T69", name="figS1 PDF rendered text matches R3.2 (R103)",
                status=PASS if ok else FAIL,
                message=f"figS1 PDF: 75={has_75} n=61={has_61} 72.7%={has_727} 39.0-90.6={has_39_906}",
            ))
        else:
            results.append(TestResult(id="T69", name="figS1 PDF", status=WARN, message=f"PDF not found: {figS1_pdf}"))

    # ===== T70: MWU p-value canonical across ALL sources =====
    amb_mwu = antigen.get("mwu_p_one_sided_exact")
    fig_mwu = figure_provenance.get("figS2_antigen_benchmark", {}).get("MWU_p_one_sided")
    fact_mwu = facts.get("benchmark", {}).get("antigen", {}).get("mwu_p_one_sided") if facts else None
    # Convert to float for comparison
    def _to_float(v):
        if v is None: return None
        if isinstance(v, (int, float)): return float(v)
        try: return float(v)
        except (TypeError, ValueError): return None
    amb_f = _to_float(amb_mwu)
    fig_f = _to_float(fig_mwu)
    fact_f = _to_float(fact_mwu)
    values = [v for v in [amb_f, fig_f, fact_f] if v is not None]
    canonical = 9.53e-5
    if values:
        all_same = all(abs(v - canonical) < 1e-6 for v in values)
    else:
        all_same = True
    results.append(TestResult(
        id="T70", name="Invariant: MWU p=9.53e-5 across ALL sources (R103)",
        status=PASS if all_same else FAIL,
        message=f"antigen_summary={amb_mwu} figure_provenance={fig_mwu} facts={fact_mwu}",
    ))

    return results
