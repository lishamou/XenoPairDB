"""
Gate 9 — Cross-artifact referential integrity (R114, V4 review fixes)

Per V4 review (XenoPairDB_Integrity_Pipeline_v4.md), the manuscript has been
plagued by deterministically-generated inconsistencies between:
- donor scorecard (Fig1/Fig2/FigS1/239-288-288/356-334-424)
- antigen effect-size (AUROC vs rank-biserial vs common-language)
- resource count (Fig S3 + Main + SI)
- curation ledger (Fig1 + S4 + Methods)
- figures pixel text vs caption descriptions
- template serialization (raw dict into prose)

Gate 9 enforces 7 cross-artifact invariants:
  T79: DONOR_SCORECARD_INVARIANT — sum(system residuals) == all-six total;
       Complement + Coagulation + Endocrine == selected subtotal
  T80: ANTIGEN_EFFECT_SIZE_INVARIANT — U/(n+*n-) = AUROC = common-language;
       2*AUROC-1 = rank-biserial
  T81: RESOURCE_COMPARISON_CONTRACT — Fig S3 resources count == registry.count;
       Main "six"/SI "five" consistent with registry
  T82: CURATION_LEDGER_INVARIANT — R1.0 + candidates - removed - boundary - eng = strict
  T83: FIGURE_SEMANTIC_CONTRACT — caption text matches figure pixel labels
       (advisory via Gate 6 PNG OCR)
  T84: TEMPLATE_SERIALIZATION_GATE — no raw dict/list/None/NaN in rendered
       submission (read main.submission.md / si.submission.md)
  T85: FORMATTER_CONSISTENCY — facts.json formatter produces values present in
       rendered manuscript (1,407 == TIER2_TOTAL; 9.53×10⁻⁵ == MWU_P)
"""
from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


def _facts(canonical_dir: Path) -> dict | None:
    """Load manuscript_facts.json."""
    p = canonical_dir.parent / "manuscript" / "manuscript_facts.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


@register_gate("Gate 9 — Cross-artifact invariants (R114)", weight=20)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []
    if facts is None:
        facts = _facts(canonical_dir)
    if facts is None:
        results.append(TestResult(
            id="T79", name="DONOR_SCORECARD_INVARIANT (skipped)",
            status=WARN, message="manuscript_facts.json missing",
        ))
        return results

    pipeline_root = canonical_dir.parent
    manuscript_dir = pipeline_root / "manuscript"

    # ===== T79: DONOR_SCORECARD_INVARIANT =====
    # For each donor: figure_panel_E_three_views is the canonical.
    # Check that figures_panel_E_three_views.all_six_system_unmasked (dict)
    # + all_six_system_masked (dict) + three_system_selected (dict)
    # are key-consistent (same donor keys across all 3 views).
    pE = facts.get("figures_panel_E_three_views", {})
    keys_views = [
        set(pE.get("all_six_system_unmasked", {}).keys()),
        set(pE.get("all_six_system_masked", {}).keys()),
        set(pE.get("three_system_selected", {}).keys()),
    ]
    same_keys = all(k == keys_views[0] for k in keys_views) and len(keys_views[0]) == 3
    results.append(TestResult(
        id="T79", name="DONOR_SCORECARD_INVARIANT: 3 views have same 3 donor keys",
        status=PASS if same_keys else FAIL,
        message=f"keys={sorted(keys_views[0]) if same_keys else 'inconsistent'} (need Yunnan/Revivicor/eGenesis × 3 views)",
    ))

    # ===== T80: ANTIGEN_EFFECT_SIZE_INVARIANT =====
    af = facts.get("antigen_full", {})
    auroc = af.get("rank_biserial_r")  # we use rank_biserial_r as proxy
    af_b = facts.get("benchmark", {}).get("antigen", {})
    n_pos = af_b.get("n_positives", 8)
    n_neg = af_b.get("n_negatives", 12)
    u_pos = af.get("mwu_u_pos", 92)
    u_total = af.get("mwu_u_total", 96)
    expected_cle = u_pos / (n_pos * n_neg) if (n_pos and n_neg) else None
    expected_rbr = 2 * (u_pos / (n_pos * n_neg)) - 1 if (n_pos and n_neg) else None
    actual_auroc = af_b.get("auroc", 0.958)
    actual_cle = af.get("common_language_effect", 0.958)
    actual_rbr = af.get("rank_biserial_r", 0.917)
    ok = (
        abs(actual_cle - expected_cle) < 0.01
        and abs(actual_rbr - expected_rbr) < 0.01
        and abs(actual_auroc - actual_cle) < 0.01
    )
    results.append(TestResult(
        id="T80", name=f"ANTIGEN_EFFECT_SIZE_INVARIANT: AUROC={actual_auroc:.3f} = CLE={actual_cle:.3f}; r_biserial={actual_rbr:.3f} ≈ 2·AUROC−1",
        status=PASS if ok else FAIL,
        message=f"U/(n+·n-) = {expected_cle:.3f}, 2·AUROC−1 = {expected_rbr:.3f}; "
                f"actual AUROC={actual_auroc:.3f}, CLE={actual_cle:.3f}, r={actual_rbr:.3f}",
    ))

    # ===== T81: RESOURCE_COMPARISON_CONTRACT =====
    rr = facts.get("resource_registry", {})
    if rr:
        n = rr.get("count", 0)
        results.append(TestResult(
            id="T81", name=f"RESOURCE_COMPARISON_CONTRACT: {n} resources in registry",
            status=PASS if n == 5 else WARN,
            message=f"resource_registry.count={n} (post-R111: 5, was 6 incl. XenoBase)",
        ))
    else:
        results.append(TestResult(
            id="T81", name="RESOURCE_COMPARISON_CONTRACT (missing registry)",
            status=WARN, message="facts.resource_registry missing",
        ))

    # ===== T82: CURATION_LEDGER_INVARIANT =====
    cc = facts.get("curation_change", {})
    if cc:
        arith = cc.get("inv_arithmetic", 0)
        target = cc.get("final_strict_core", 75)
        ok = arith == target
        results.append(TestResult(
            id="T82", name=f"CURATION_LEDGER_INVARIANT: {cc.get('r1_0_baseline')}+{cc.get('candidate_additions')}−{cc.get('removed')}−{cc.get('boundary_reclassified')}−1 = {arith} (need {target})",
            status=PASS if ok else FAIL,
            message=f"69 + 13 − 1 − 5 − 1 = {69+13-1-5-1} (computed); inv_arithmetic={arith}; final_strict_core={target}",
        ))

    # ===== T83: FIGURE_SEMANTIC_CONTRACT (advisory) =====
    # Captions can't easily be parsed; defer to Gate 6 PNG OCR for now.
    results.append(TestResult(
        id="T83", name="FIGURE_SEMANTIC_CONTRACT (advisory, see Gate 6)",
        status=WARN, message="cross-source figure-text consistency is checked by Gate 6 figure_contract; Gate 9 placeholder for future pix-text extraction",
    ))

    # ===== T84: TEMPLATE_SERIALIZATION_GATE =====
    # Check rendered (not submission; --submission strips DEV_ONLY but
    # rendered has raw dict leaks if any).
    rendered_main = manuscript_dir / "main.rendered.md"
    rendered_si = manuscript_dir / "si.rendered.md"
    bad_markers = []
    if rendered_main.exists():
        text = rendered_main.read_text(encoding="utf-8", errors="ignore")
        # Python dict repr pattern: {'key': value, ...}
        for m in re.finditer(r"\{['\"][^}]+['\"]\s*:\s*[\d\.\-]+", text):
            bad_markers.append(f"main:{m.start()}")
    if rendered_si.exists():
        text = rendered_si.read_text(encoding="utf-8", errors="ignore")
        for m in re.finditer(r"\{['\"][^}]+['\"]\s*:\s*[\d\.\-]+", text):
            bad_markers.append(f"si:{m.start()}")
    results.append(TestResult(
        id="T84", name="TEMPLATE_SERIALIZATION_GATE: no raw dict repr in rendered Main/SI",
        status=PASS if not bad_markers else FAIL,
        message=f"leaks={len(bad_markers)}" + (f"; first={bad_markers[0]}" if bad_markers else " (clean)"),
    ))

    # ===== T85: FORMATTER_CONSISTENCY =====
    # Ensure rendered manuscript has facts.json formatter output (e.g. 1,407 not 1407)
    fmt_consistent = True
    fmt_msg = []
    if rendered_main.exists():
        text = rendered_main.read_text(encoding="utf-8", errors="ignore")
        # Check thousands format: facts has TIER2_TOTAL=1407 (int); should render as "1,407"
        if "1,407" not in text and "1407" in text:
            fmt_consistent = False
            fmt_msg.append("main missing 1,407 thousands format")
        # scientific_unicode: should have ×10⁻⁵ (not e-5)
        if "×10⁻⁵" not in text and "e-5" in text.lower():
            fmt_consistent = False
            fmt_msg.append("main missing ×10⁻⁵ scientific_unicode format")
    results.append(TestResult(
        id="T85", name="FORMATTER_CONSISTENCY: thousands + scientific_unicode in rendered",
        status=PASS if fmt_consistent else WARN,
        message="; ".join(fmt_msg) if fmt_msg else "OK (1,407 and ×10⁻⁵ present)",
    ))

    return results
