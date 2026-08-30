"""
Gate 2 — Analytical integrity tests (Tests 16-30)

These tests verify that the derived analytical outputs (ISI, benchmark, scorecard,
identity, sensitivity) match the frozen expected values. They re-compute key
metrics from the canonical TSV and compare against analysis_results_R3.2.json.
"""
from __future__ import annotations

import csv
import json
import statistics
from collections import Counter
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


@register_gate("Gate 2 — Analytical integrity", weight=15)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []

    tier1 = _load_tsv(canonical_dir / "tier1_frozen.tsv")
    tier2 = _load_tsv(canonical_dir / "tier2_atlas.tsv")
    surfaceome = _load_tsv(canonical_dir / "surfaceome.tsv")
    analysis = _load_json(canonical_dir / "analysis_results.json")
    antigen = _load_json(canonical_dir / "antigenicity_benchmark_summary.json")

    # ===== Test 16: Tier-2 row count = 1,407 =====
    actual = len(tier2)
    results.append(TestResult(
        id="T16", name="Tier-2 row count = 1,407",
        status=PASS if actual == 1407 else FAIL,
        message=f"actual={actual}", expected="1407",
    ))

    # ===== Test 17: Tier-2 raw pairs = 5,835 (OmniPath + CellPhoneDB dedup) =====
    # R100 fix: read from facts.json (L4 canonical) — should be exact value
    raw_omnipath = facts.get("tier2_raw_pairs", None) if facts else None
    results.append(TestResult(
        id="T17", name="Tier-2 raw pairs = 5,835",
        status=PASS if raw_omnipath == 5835 else FAIL,
        message=f"raw_pairs={raw_omnipath} (from facts.json tier2_raw_pairs)",
        actual=str(raw_omnipath), expected="5835",
    ))

    # ===== Test 18: Tier-2 filtered = 4,689 (after membrane filter; final 1,407 after ortholog resolution) =====
    # R100 fix: read from facts.json — should be exact value
    filtered = facts.get("tier2_filtered_pairs", None) if facts else None
    results.append(TestResult(
        id="T18", name="Tier-2 filtered pairs = 4,689",
        status=PASS if filtered == 4689 else FAIL,
        message=f"filtered={filtered} (from facts.json tier2_filtered_pairs)",
        actual=str(filtered), expected="4689",
    ))

    # ===== Test 19: Median sequence identity from TSV ≈ 72.7 =====
    identities = []
    for r in tier1:
        v = r.get("identity", "").strip()
        if v and v != "n/a":
            try:
                identities.append(float(v))
            except ValueError:
                continue
    actual_median = round(statistics.median(identities), 1)
    expected = 72.7
    ok = abs(actual_median - expected) <= 0.1
    results.append(TestResult(
        id="T19", name=f"Median identity ≈ {expected}",
        status=PASS if ok else FAIL,
        message=f"actual={actual_median} (n={len(identities)})",
        actual=str(actual_median), expected=str(expected),
    ))

    # ===== Test 20: Range 39.0-90.6 =====
    actual_min = round(min(identities), 1)
    actual_max = round(max(identities), 1)
    expected_min, expected_max = 39.0, 90.6
    ok = abs(actual_min - expected_min) <= 0.1 and abs(actual_max - expected_max) <= 0.1
    results.append(TestResult(
        id="T20", name=f"Identity range = {expected_min}-{expected_max}",
        status=PASS if ok else FAIL,
        message=f"actual range = {actual_min}-{actual_max}",
        actual=f"{actual_min}-{actual_max}", expected=f"{expected_min}-{expected_max}",
    ))

    # ===== Test 21: n_scoreable pairs = 61 =====
    actual_n = len(identities)
    expected_n = 61
    results.append(TestResult(
        id="T21", name="n_scoreable pairs = 61",
        status=PASS if actual_n == expected_n else FAIL,
        message=f"actual={actual_n}", expected=str(expected_n),
    ))

    # ===== Test 22: Antigen AUROC = 0.958 =====
    auroc = antigen.get("v2", {}).get("auroc") if isinstance(antigen.get("v2"), dict) else None
    ok = auroc is not None and abs(auroc - 0.958) <= 0.001
    results.append(TestResult(
        id="T22", name="Antigen AUROC = 0.958",
        status=PASS if ok else FAIL,
        message=f"actual={auroc}", expected="0.958",
    ))

    # ===== Test 23: Antigen AUPRC = 0.942 =====
    auprc = antigen.get("v2", {}).get("auprc") if isinstance(antigen.get("v2"), dict) else None
    ok = auprc is not None and abs(auprc - 0.942) <= 0.001
    results.append(TestResult(
        id="T23", name="Antigen AUPRC = 0.942",
        status=PASS if ok else FAIL,
        message=f"actual={auprc}", expected="0.942",
    ))

    # ===== Test 24: MWU p one-sided matches canonical (R100: 9.53e-5 scipy exact) =====
    # R100 fix: accept ONLY the scipy exact value 9.53e-5. The previous 2.38e-4 was
    # an outdated cached value in antigenicity_benchmark_summary.json; it has been
    # corrected. T24 now asserts EQUALITY (not "in {a, b}") to enforce Single Source of Truth.
    mwu = antigen.get("mwu_p_one_sided_exact")
    canonical = 9.53e-5
    ok = mwu is not None and abs(mwu - canonical) < 1e-7
    results.append(TestResult(
        id="T24", name=f"MWU p one-sided = {canonical} (scipy exact; R100 canonical)",
        status=PASS if ok else FAIL,
        message=f"actual={mwu} expected={canonical}",
        actual=str(mwu), expected=str(canonical),
    ))

    # ===== Test 25: ISI_v2 = ISI × expr_weight for ALL 75 axes (R100: no sampling) =====
    # R100 fix: with only 75 axes, sampling n=30 leaves 45 axes unverified. Test ALL.
    bad = []
    non_control_rows = [r for r in tier1 if r.get("evidence") != "control"]
    for r in non_control_rows:
        try:
            isi = float(r["ISI"])
            ew = float(r["expr_weight"])
            expected_isi_v2 = isi * ew
            actual_isi_v2 = float(r["ISI_v2"])
            # Allow tolerance 0.5 (manual rounding allowed)
            if abs(expected_isi_v2 - actual_isi_v2) > 0.5:
                bad.append((r["pair_id"], expected_isi_v2, actual_isi_v2))
        except (ValueError, KeyError):
            continue
    results.append(TestResult(
        id="T25", name=f"ISI_v2 = ISI × expr_weight (ALL {len(non_control_rows)} non-control axes)",
        status=PASS if not bad else FAIL,
        message=f"violations={len(bad)} of {len(non_control_rows)} (sample: {bad[:3]})",
    ))

    # ===== Test 26: 13 controls have func_sev = 0 =====
    controls = [r for r in tier1 if r.get("evidence") == "control"]
    bad = [(r["pair_id"], r["func_sev"]) for r in controls if r["func_sev"] not in ("0", "0.0")]
    results.append(TestResult(
        id="T26", name="13 controls have func_sev = 0",
        status=PASS if len(controls) == 13 and not bad else FAIL,
        message=f"n_controls={len(controls)} violations={bad if bad else '∅'}",
    ))

    # ===== Test 27: 62 incompatible pairs have func_sev > 0 =====
    incompatible = [r for r in tier1 if r.get("evidence") in ("direct", "disputed")]
    bad = [(r["pair_id"], r["func_sev"]) for r in incompatible
           if r["func_sev"] in ("0", "0.0", "n/a", "")]
    results.append(TestResult(
        id="T27", name="62 incompatible pairs have func_sev > 0",
        status=PASS if len(incompatible) == 62 and not bad else FAIL,
        message=f"n_incompatible={len(incompatible)} violations={len(bad)}",
    ))

    # ===== Test 28: Donor scorecard has 6 genotypes =====
    scorecard = _load_tsv(canonical_dir / "donor_scorecard.tsv")
    n_genotypes = len(scorecard)
    results.append(TestResult(
        id="T28", name="Donor scorecard has 6 genotypes",
        status=PASS if n_genotypes == 6 else FAIL,
        message=f"actual={n_genotypes}",
    ))

    # ===== Test 29: Gene counts in pathway = 8-32 range =====
    # R100 fix: read from facts.json (L4 derived) — should be canonical values
    if facts and "clinical_pathway_counts" in facts:
        clinical = facts["clinical_pathway_counts"]
        min_c = facts.get("clinical_pathway_min_count", min(clinical.values()))
        max_c = facts.get("clinical_pathway_max_count", max(clinical.values()))
        ok = 8 <= min_c and max_c <= 32
        results.append(TestResult(
            id="T29", name="Pathway gene counts in 8-32 range",
            status=PASS if ok else FAIL,
            message=f"min={min_c} max={max_c} (from facts.json clinical_pathway_counts)",
            actual=f"{min_c}-{max_c}", expected="8-32",
        ))
    else:
        results.append(TestResult(id="T29", name="Pathway gene counts", status=FAIL,
                                  message="clinical_pathway_counts missing from facts.json"))

    # ===== Test 30: Sensitivity Spearman ρ >= 0.92 across variants =====
    # R100 fix: canonical from facts.json sensitivity.min_spearman_overall = 0.92
    if facts and "sensitivity" in facts:
        sens = facts["sensitivity"]
        spearman_min = sens.get("min_spearman_overall", sens.get("min_spearman"))
        threshold = 0.92
        if spearman_min is not None:
            results.append(TestResult(
                id="T30", name=f"Sensitivity Spearman ρ >= {threshold} across variants",
                status=PASS if spearman_min >= threshold else FAIL,
                message=f"min_spearman={spearman_min} (from facts.json sensitivity.min_spearman_overall)",
                actual=str(spearman_min), expected=str(threshold),
            ))
        else:
            results.append(TestResult(id="T30", name="Sensitivity Spearman", status=FAIL,
                                      message="sensitivity.min_spearman_overall missing"))
    else:
        results.append(TestResult(id="T30", name="Sensitivity Spearman", status=FAIL,
                                  message="sensitivity block missing from facts.json"))

    return results


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
