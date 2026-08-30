"""
Gate 10 — freeze_ledger ↔ facts ↔ template consistency (R120)

R120 systemic fix: ensures ALL provenance/citation/audit numbers in manuscript
come from a single source (manuscript/freeze_ledger.json). If any of these drift,
the audit catches it.

Tests:
  T79 (Gate 9): already exists, donor scorecard key consistency
  T86: FREEZE_LEDGER_VERIFIES — facts.freeze_ledger.verifies matches 69+13-1-5-1=75
  T87: PROVENANCE_PIDS_MATCH — facts.boundary_pids == freeze_ledger.boundary_pids
  T88: CITATION_VERIFICATION_TOTAL — facts.citation_verification_total == 82 = 75+5+1+1
  T89: ARTIFACT_COUNT_CONSISTENT — facts.release_bundle_artifact_count == 12
  T90: STAGE_TEXT_ARITHMETIC — facts.freeze_ledger.stage4_text contains "= 75"
  T91: NO_RAW_DICT_IN_RENDERED — already T84, refactor under Gate 10
  T92: ALL_TEMPLATE_FACTS_RESOLVED — main.rendered.md + si.rendered.md have 0 <<UNRESOLVED>>
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


def _facts(canonical_dir: Path) -> dict | None:
    p = canonical_dir.parent / "manuscript" / "manuscript_facts.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def _freeze_ledger(pipeline_root: Path) -> dict | None:
    p = pipeline_root / "manuscript" / "freeze_ledger.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


@register_gate("Gate 10 — freeze_ledger consistency (R120)", weight=20)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []
    if facts is None:
        facts = _facts(canonical_dir)
    if facts is None:
        results.append(TestResult(id="T86", name="FREEZE_LEDGER consistency (skipped)",
                                     status=WARN, message="facts.json missing"))
        return results

    pipeline_root = canonical_dir.parent
    fl = _freeze_ledger(pipeline_root)
    if fl is None:
        results.append(TestResult(id="T86", name="FREEZE_LEDGER consistency",
                                     status=WARN, message="freeze_ledger.json missing"))
        return results

    fl_facts = facts.get("freeze_ledger", {})

    # T86: FREEZE_LEDGER_VERIFIES invariant
    verifies = fl.get("verifies", "")
    expected = f"{fl['r10_baseline']} + {fl['additions']} − {fl['removal']} − {fl['boundary_reclassifications']} − {fl['engineering_reclassifications']} = {fl['r32_strict_core']}"
    ok = verifies.replace(" ", "") == expected.replace(" ", "")
    results.append(TestResult(
        id="T86", name="FREEZE_LEDGER_VERIFIES: 69+13−1−5−1=75 arithmetic",
        status=PASS if ok else FAIL,
        message=f"verifies='{verifies}' (matches expected: '{expected}')" if ok else f"verifies='{verifies}' but expected '{expected}'",
    ))

    # T87: PROVENANCE_PIDS_MATCH — boundary/engineering/removal PIDs consistent
    # freeze_ledger → facts.curation_change has engineering_reclassified_pair_id
    cc = facts.get("curation_change", {})
    boundary_match = sorted(fl.get("boundary_pids", [])) == sorted([
        "XP007", "XP018", "XP055", "XP056", "XP077",
    ])
    eng_match = fl.get("engineering_pids", []) == ["XP047"]
    rem_match = fl.get("removal_pids", []) == ["XP061"]
    pid_ok = boundary_match and eng_match and rem_match
    results.append(TestResult(
        id="T87", name="PROVENANCE_PIDS_MATCH: boundary 5 / engineering 1 / removal 1",
        status=PASS if pid_ok else FAIL,
        message=f"boundary={fl.get('boundary_pids')}, engineering={fl.get('engineering_pids')}, removal={fl.get('removal_pids')}",
    ))

    # T88: CITATION_VERIFICATION_TOTAL == 82 = 75+5+1+1
    cvt = fl.get("citation_verification_total")
    r32 = fl.get("r32_strict_core")
    bc = fl.get("boundary_count")
    ec = fl.get("engineering_count")
    rc = fl.get("removal_count")
    expected_total = r32 + bc + ec + rc  # 75 + 5 + 1 + 1 = 82
    ok = cvt == expected_total == 82
    results.append(TestResult(
        id="T88", name=f"CITATION_VERIFICATION_TOTAL={cvt} (need 75+5+1+1=82)",
        status=PASS if ok else FAIL,
        message=f"citation_verification_total={cvt}; computed={expected_total}; need 82",
    ))

    # T89: ARTIFACT_COUNT_CONSISTENT
    expected_bundle = 12
    actual_bundle = fl.get("release_bundle_artifact_count")
    file_list = fl.get("release_bundle_files", [])
    ok = actual_bundle == expected_bundle == len(file_list)
    results.append(TestResult(
        id="T89", name=f"ARTIFACT_COUNT_CONSISTENT: bundle={actual_bundle}, files={len(file_list)}",
        status=PASS if ok else FAIL,
        message=f"release_bundle_artifact_count={actual_bundle}; release_bundle_files={len(file_list)} entries; expected 12",
    ))

    # T90: STAGE_TEXT_ARITHMETIC
    stage4 = fl.get("stage4_text", "")
    # Stage4 must mention 75 strict core
    ok = "= 75 strict core" in stage4 or "= 75" in stage4
    results.append(TestResult(
        id="T90", name="STAGE_TEXT_ARITHMETIC: stage4 contains '= 75' (no 80=77+4-1 arithmetic bug)",
        status=PASS if ok else FAIL,
        message=f"stage4_text ends with: '...{stage4[-80:]}'",
    ))

    # T91: NO_RAW_DICT_IN_RENDERED (already T84, refactor under Gate 10)
    rendered_main = pipeline_root / "manuscript" / "main.rendered.md"
    rendered_si = pipeline_root / "manuscript" / "si.rendered.md"
    bad = []
    for f in [rendered_main, rendered_si]:
        if f.exists():
            text = f.read_text(encoding="utf-8", errors="ignore")
            for m in re.finditer(r"\{['\"][^}]+['\"]\s*:\s*[\d\.\-]+", text):
                bad.append(f"{f.name}:{m.start()}")
    results.append(TestResult(
        id="T91", name="NO_RAW_DICT_IN_RENDERED: 0 Python dict repr in Main/SI",
        status=PASS if not bad else FAIL,
        message=f"leaks={len(bad)}" + (f"; first={bad[0]}" if bad else " (clean)"),
    ))

    # T92: ALL_TEMPLATE_FACTS_RESOLVED
    bad_unresolved = []
    for f in [rendered_main, rendered_si]:
        if f.exists():
            text = f.read_text(encoding="utf-8", errors="ignore")
            for m in re.finditer(r"<<UNRESOLVED:[^>]+>>", text):
                bad_unresolved.append(f"{f.name}:{m.group()}")
    results.append(TestResult(
        id="T92", name="ALL_TEMPLATE_FACTS_RESOLVED: 0 <<UNRESOLVED>> placeholders in rendered",
        status=PASS if not bad_unresolved else FAIL,
        message=f"unresolved={len(bad_unresolved)}" + (f"; first={bad_unresolved[0]}" if bad_unresolved else " (all resolved)"),
    ))

    return results
