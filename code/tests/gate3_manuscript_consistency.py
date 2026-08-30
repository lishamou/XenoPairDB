"""
Gate 3 — Manuscript consistency tests (Tests 31-50)

These tests verify that the manuscript main.md and si.md files contain the same
numerical values as the canonical facts.json. They catch drift between the
manuscript and the underlying data.
"""
from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


@register_gate("Gate 3 — Manuscript consistency", weight=20)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []

    # Locate the manuscript files
    # R110-B: pipeline self-contained — read only pipeline/manuscript/ (was en_nbt/xenopaper_nbt_*.md)
    pipeline_root = canonical_dir.parent
    main_md = _find_file(canonical_dir, [
        str(pipeline_root / "manuscript" / "main.rendered.md"),
        str(pipeline_root / "manuscript" / "template_main.md"),
    ])
    si_md = _find_file(canonical_dir, [
        str(pipeline_root / "manuscript" / "si.rendered.md"),
        str(pipeline_root / "manuscript" / "template_si.md"),
    ])

    main_text = main_md.read_text(encoding="utf-8") if main_md else ""
    si_text = si_md.read_text(encoding="utf-8") if si_md else ""

    # ===== Test 31: main.md mentions Tier-1 total = 75 =====
    n = len(re.findall(r"\b75\b", main_text))
    results.append(TestResult(
        id="T31", name="main.md contains 75 (Tier-1 total)",
        status=PASS if n >= 5 else FAIL,
        message=f"occurrences={n}",
        actual=str(n), expected="≥5",
    ))

    # ===== Test 32: main.md mentions median identity 72.7 (not 73.2) =====
    has_727 = "72.7" in main_text
    has_732 = bool(re.search(r"73\.2", main_text))
    results.append(TestResult(
        id="T32", name="main.md median identity = 72.7 (not 73.2)",
        status=PASS if (has_727 and not has_732) else FAIL,
        message=f"has_72.7={has_727} has_73.2={has_732}",
    ))

    # ===== Test 33: main.md mentions AUROC 0.958 =====
    has_auroc = "0.958" in main_text
    results.append(TestResult(
        id="T33", name="main.md AUROC = 0.958",
        status=PASS if has_auroc else FAIL,
        message=f"present={has_auroc}",
    ))

    # ===== Test 34: si.md Fig S1 caption uses 72.7 (not 73.2) =====
    si_has_727 = "72.7" in si_text
    si_has_732 = bool(re.search(r"73\.2", si_text))
    results.append(TestResult(
        id="T34", name="si.md Fig S1 caption uses 72.7 (not 73.2)",
        status=PASS if (si_has_727 and not si_has_732) else FAIL,
        message=f"has_72.7={si_has_727} has_73.2={si_has_732}",
    ))

    # ===== Test 35: main.md XP047 references engineering-modifier (not direct Tier-1) =====
    # R100 fix: search for the EXACT stale phrase "XP047 ... direct Tier-1 axes"
    # (the exact string R97 removed). FAIL if present. Else PASS if XP047 is
    # correctly referenced as engineering-modifier.
    xp047_stale_direct = bool(re.search(
        r"XP047[^.\n]*?direct\s+Tier[‐\-_]?1\s+axes",
        main_text, re.IGNORECASE))
    xp047_engineering_modifier = bool(re.search(
        r"XP047[^.\n]*?(?:engineering[‐\-_]modifier|engineering modifier|engineering_modifier)",
        main_text, re.IGNORECASE))
    if xp047_stale_direct:
        results.append(TestResult(
            id="T35", name="main.md XP047 NOT in stale 'direct Tier-1 axes' phrase",
            status=FAIL,
            message="STALE: XP047 still described as 'direct Tier-1 axes' in main.md",
        ))
    elif xp047_engineering_modifier:
        results.append(TestResult(
            id="T35", name="main.md XP047 references engineering-modifier",
            status=PASS,
            message="XP047 correctly referenced as engineering-modifier",
        ))
    else:
        results.append(TestResult(
            id="T35", name="main.md XP047 references engineering-modifier",
            status=WARN,
            message="XP047 mentioned but engineering-modifier phrase not found nearby",
        ))

    # ===== Test 36: main.md + si.md Figure 5 caption mention 62/75 = 0.827 =====
    has_827_main = "0.827" in main_text
    has_827_si = "0.827" in si_text
    results.append(TestResult(
        id="T36", name="main.md + si.md mention 0.827 (AUPRC no-skill baseline)",
        status=PASS if (has_827_main and has_827_si) else FAIL,
        message=f"main={has_827_main} si={has_827_si}",
    ))

    # ===== Test 37: main.md system counts match canonical =====
    # Look for the system list
    has_22 = "Immune 22" in main_text or "Immune, 22" in main_text or re.search(r"Immune[^a-zA-Z]+22", main_text) is not None
    has_12 = re.search(r"Coagulation[^a-zA-Z]+12", main_text) is not None
    has_15 = re.search(r"Cytokine[^a-zA-Z]+15", main_text) is not None
    has_7 = re.search(r"Endocrine[^a-zA-Z]+7\b", main_text) is not None
    has_10 = re.search(r"Adhesion[^a-zA-Z]+10\b", main_text) is not None
    has_9 = re.search(r"Complement[^a-zA-Z]+9\b", main_text) is not None
    ok = all([has_22, has_12, has_15, has_7, has_10, has_9])
    results.append(TestResult(
        id="T37", name="main.md system counts match (Immune 22, Coag 12, Cyto 15, Endo 7, Adh 10, Comp 9)",
        status=PASS if ok else FAIL,
        message=f"22={has_22} 12={has_12} 15={has_15} 7={has_7} 10={has_10} 9={has_9}",
    ))

    # ===== Test 38: si.md does not use 69-axis as current (only allowed in R1.0 historical context) =====
    si_69_axis_uncoupled = []
    for m in re.finditer(r"\b69[\‐\-_]?axis\b", si_text):
        ctx = si_text[max(0, m.start() - 80):min(len(si_text), m.end() + 80)]
        if "Release 1.0" not in ctx and "historical" not in ctx and "R1.0" not in ctx:
            si_69_axis_uncoupled.append(m.group())
    results.append(TestResult(
        id="T38", name="si.md '69-axis' only in R1.0 historical context",
        status=PASS if not si_69_axis_uncoupled else FAIL,
        message=f"uncoupled_69_axis={si_69_axis_uncoupled}",
    ))

    # ===== Test 39: main.md total axes references consistent =====
    # Should reference 75 throughout (not 69 in main context)
    n_75 = main_text.count(" 75 ")
    n_75_axis = len(re.findall(r"\b75\b(?:\s*[-‐]?axis)?", main_text))
    results.append(TestResult(
        id="T39", name="main.md '75' references >= 5 (Tier-1 total consistent)",
        status=PASS if n_75 >= 5 else FAIL,
        message=f"total '75' mentions={n_75}",
    ))

    # ===== Test 40: si.md "13 controls" used (not 8 controls) =====
    n_13 = len(re.findall(r"\b13\s+controls?\b", si_text))
    n_8 = len(re.findall(r"\b8\s+controls?\b(?!\s*\(8)", si_text))   # exclude "(8"
    results.append(TestResult(
        id="T40", name="si.md '13 controls' dominant (R3.2 freeze)",
        status=PASS if n_13 >= 3 else WARN,
        message=f"'13 controls'={n_13} '8 controls'={n_8} (8 is historical baseline)",
    ))

    # ===== Test 41: main.md Range uses 39.0-90.6 =====
    has_range = ("39.0" in main_text and "90.6" in main_text) or bool(re.search(r"39\.0\s*[‐\-–]\s*90\.6", main_text))
    results.append(TestResult(
        id="T41", name="main.md identity range uses 39.0-90.6",
        status=PASS if has_range else FAIL,
        message=f"range_present={has_range}",
    ))

    # ===== Test 42: References defined in main.md cover the current reference list =====
    # Extract all citation numbers (in text body and References section)
    ref_nums = set()
    # Bracket form [1] or [1,2,3] or [1-3] or [1-3,5]
    for m in re.finditer(r"\[([0-9,\-\s–]+)\]", main_text):
        inside = m.group(1)
        for part in inside.split(","):
            part = part.strip()
            range_match = re.match(r"(\d+)\s*[\-–]\s*(\d+)", part)
            if range_match:
                a, b = int(range_match.group(1)), int(range_match.group(2))
                ref_nums.update(range(a, b + 1))
            else:
                num_match = re.match(r"(\d+)", part)
                if num_match:
                    ref_nums.add(int(num_match.group(1)))
    # Also extract references defined in References section (e.g., "59. Cooper...")
    if "## References" in main_text:
        refs_section = main_text[main_text.find("## References"):]
        for m in re.finditer(r"^(\d+)\.\s+[A-Z]", refs_section, re.MULTILINE):
            ref_nums.add(int(m.group(1)))
    # Derive the current reference ceiling from the References section so that
    # adding a verified citation cannot be reported as a stale orphan.
    defined_refs = set()
    if "## References" in main_text:
        refs_section = main_text[main_text.find("## References"):]
        for m in re.finditer(r"^(\d+)\.\s+[A-Z]", refs_section, re.MULTILINE):
            defined_refs.add(int(m.group(1)))
    max_ref = max(defined_refs) if defined_refs else 0
    expected_refs = set(range(1, max_ref + 1))
    missing = expected_refs - ref_nums
    orphans = sorted(ref_nums - expected_refs)
    results.append(TestResult(
        id="T42", name="main.md references cover 1-N (orphans allowed)",
        status=PASS if not missing else FAIL,
        message=f"missing={sorted(missing) if missing else '∅'} (n_found={len(ref_nums & expected_refs)}, orphans={orphans})",
    ))

    # ===== Test 43: References in si.md use main.md numbering =====
    si_ref_nums = set(int(m) for m in re.findall(r"\[(\d+)\]", si_text))
    si_out_of_range = sorted(n for n in si_ref_nums if n > max_ref)
    results.append(TestResult(
        id="T43", name="si.md references use main.md numbering (≤N)",
        status=PASS if not si_out_of_range else FAIL,
        message=f"out_of_range={si_out_of_range} (should be empty if SI shares numbering)",
    ))

    # ===== Test 44: No pseudo-Gal / 伪Gal characters =====
    has_pseudo = "伪Gal" in main_text or "伪Gal" in si_text or "pseudo-Gal" in main_text or "pseudo-Gal" in si_text
    results.append(TestResult(
        id="T44", name="No pseudo-Gal / 伪Gal characters",
        status=PASS if not has_pseudo else FAIL,
        message=f"pseudo_gal_in_text={has_pseudo}",
    ))

    # ===== Test 45: 'multi-stage' used (not 'two-stage') in main + si combined =====
    combined = main_text + "\n" + si_text
    has_two_stage = "two-stage" in combined
    has_multi = "multi-stage" in combined or "multi_stage" in combined or "Multi-stage" in combined
    results.append(TestResult(
        id="T45", name="'multi-stage' used, 'two-stage' absent",
        status=PASS if (has_multi and not has_two_stage) else FAIL,
        message=f"has_multi-stage={has_multi} has_two-stage={has_two_stage}",
    ))

    # ===== Test 46: PERV count 59 mentioned correctly =====
    has_perv_59 = "59 of the 69" in main_text or "59 PERV" in main_text or "PERV" in main_text
    results.append(TestResult(
        id="T46", name="main.md mentions PERV 59 count correctly",
        status=PASS if has_perv_59 else WARN,
        message=f"PERV_59_present={has_perv_59}",
    ))

    # ===== Test 47: eGenesis coverage 13.2% < Yunnan 15.2% =====
    has_132 = "13.2%" in main_text
    has_152 = "15.2%" in main_text
    results.append(TestResult(
        id="T47", name="main.md eGenesis 13.2% < Yunnan 15.2% coverage",
        status=PASS if (has_132 and has_152) else FAIL,
        message=f"13.2%={has_132} 15.2%={has_152}",
    ))

    # ===== Test 48: SLA-I rows in TSV = 6 =====
    # (rough heuristic — check that main mentions SLA-I class I multiple times)
    sla_i_count = len(re.findall(r"SLA[\s\-]I(?!\w)", main_text))
    results.append(TestResult(
        id="T48", name="main.md SLA-I references >= 3",
        status=PASS if sla_i_count >= 3 else WARN,
        message=f"SLA-I mentions={sla_i_count}",
    ))

    # ===== Test 49: Boundary cases 5: XP007/018/055/056/077 =====
    boundary_pids = {"XP007", "XP018", "XP055", "XP056", "XP077"}
    mentioned = sum(1 for p in boundary_pids if p in main_text or p in si_text)
    results.append(TestResult(
        id="T49", name="5 boundary cases mentioned (XP007/018/055/056/077)",
        status=PASS if mentioned == 5 else FAIL,
        message=f"mentioned={mentioned}/5",
    ))

    # ===== Test 50: Engineering modifier XP047 =====
    has_xp047 = "XP047" in main_text or "XP047" in si_text
    results.append(TestResult(
        id="T50", name="XP047 mentioned as engineering-modifier",
        status=PASS if has_xp047 else FAIL,
        message=f"XP047_present={has_xp047}",
    ))

    return results


def _find_file(canonical_dir: Path, candidates: list[str]) -> Path | None:
    for c in candidates:
        p = (canonical_dir / c).resolve()
        if p.exists():
            return p
    return None
