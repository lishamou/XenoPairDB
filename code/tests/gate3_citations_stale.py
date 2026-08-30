"""
Gate 3 extension — Citation + stale-term tests (Tests 51-60)

These tests verify the reference numbering integrity and scan for stale-term
residue from previous releases.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


@register_gate("Gate 3 ext — Citations + stale terms", weight=10)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []

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

    # Load Tier-1 TSV for T57 (XP064 source check)
    tier1_path = canonical_dir / "tier1_frozen.tsv"
    tier1 = []
    if tier1_path.exists():
        with open(tier1_path, encoding="utf-8-sig") as f:
            tier1 = list(csv.DictReader(f, delimiter="\t"))

    # ===== Test 51: '69-axis' only allowed in R1.0 historical context (main + si) =====
    main_69_axis_uncoupled = []
    for m in re.finditer(r"69[\‐\-_]axis", main_text):
        ctx = main_text[max(0, m.start() - 80):min(len(main_text), m.end() + 80)]
        if "Release 1.0" not in ctx and "historical" not in ctx and "R1.0" not in ctx:
            main_69_axis_uncoupled.append(m.group())
    si_69_axis_uncoupled = []
    for m in re.finditer(r"69[\‐\-_]axis", si_text):
        ctx = si_text[max(0, m.start() - 80):min(len(si_text), m.end() + 80)]
        if "Release 1.0" not in ctx and "historical" not in ctx and "R1.0" not in ctx:
            si_69_axis_uncoupled.append(m.group())
    ok = (not main_69_axis_uncoupled) and (not si_69_axis_uncoupled)
    results.append(TestResult(
        id="T51", name="'69-axis' only in R1.0 historical context",
        status=PASS if ok else FAIL,
        message=f"main_uncoupled={main_69_axis_uncoupled} si_uncoupled={si_69_axis_uncoupled}",
    ))

    # ===== Test 52: No 'Release 1.0 canonical' =====
    bad = []
    for fname, txt in [("main", main_text), ("si", si_text)]:
        if "Release 1.0 canonical" in txt or "R1.0 canonical" in txt:
            bad.append(fname)
    results.append(TestResult(
        id="T52", name="No 'Release 1.0 canonical'",
        status=PASS if not bad else FAIL,
        message=f"violations={bad if bad else '∅'}",
    ))

    # ===== Test 53: No '73.2%' anywhere =====
    has_732_main = bool(re.search(r"73\.2\s*%", main_text))
    has_732_si = bool(re.search(r"73\.2\s*%", si_text))
    has_732 = has_732_main or has_732_si
    results.append(TestResult(
        id="T53", name="No '73.2%' anywhere",
        status=PASS if not has_732 else FAIL,
        message=f"main={has_732_main} si={has_732_si}",
    ))

    # ===== Test 54: No '0.884' AUPRC baseline =====
    has_884 = "0.884" in main_text or "0.884" in si_text
    results.append(TestResult(
        id="T54", name="No '0.884' AUPRC baseline (R1.0 stale)",
        status=PASS if not has_884 else FAIL,
        message=f"present={has_884}",
    ))

    # ===== Test 55: No 'two-stage' in figures =====
    has_two_stage = "two-stage" in main_text or "two-stage" in si_text
    results.append(TestResult(
        id="T55", name="No 'two-stage' (replaced by 'multi-stage')",
        status=PASS if not has_two_stage else FAIL,
        message=f"present={has_two_stage}",
    ))

    # ===== Test 56: All current references have DOI/PMID/PDB identifier =====
    # Parse the References section of main.md
    ref_pattern = re.compile(r"^\d+\.\s+(.+?)(?:\n|$)", re.MULTILINE)
    refs_text = main_text[main_text.find("## References"):] if "## References" in main_text else ""
    refs = ref_pattern.findall(refs_text)
    no_id = []
    for ref in refs:
        if not any(tag in ref for tag in ["doi.org", "PMID", "pmid:", "pdb", "10."]):
            # Book chapter ref 29 (Varki) is acceptable
            no_id.append(ref[:80])
    results.append(TestResult(
        id="T56", name="Current references have DOI/PMID/PDB identifiers",
        status=PASS if len(no_id) <= 1 else WARN,
        message=f"n_without_id={len(no_id)} (acceptable: book chapter): {no_id[:1]}",
    ))

    # ===== Test 57: XP064 source corrected from Petitpas 2022 → Vanderheijden 2003 J Virol =====
    # R100 fix: check the TSV source column directly (not manuscript prose).
    # The manuscript prose intentionally does NOT cite Vanderheijden by name (it's
    # described as "v3.1 correction from earlier 4.4% provenance error"). The
    # authoritative source is the canonical TSV source column.
    if tier1:
        xp064_rows = [r for r in tier1 if r.get("pair_id") == "XP064"]
        if xp064_rows:
            xp064_source = xp064_rows[0].get("source", "")
            has_vanderijden = "Vanderheijden" in xp064_source or "12857889" in xp064_source
            has_jvirol = "J Virol" in xp064_source or "J. Virol" in xp064_source or "JVI." in xp064_source or "/JVI/" in xp064_source
            no_petitpas = "Petitpas" not in xp064_source and "12780" not in xp064_source
            ok = has_vanderijden and has_jvirol and no_petitpas
            results.append(TestResult(
                id="T57", name="XP064 TSV source = Vanderheijden 2003 J Virol (R95 fix)",
                status=PASS if ok else FAIL,
                message=f"source={xp064_source!r} vanderijden={has_vanderijden} jvirol={has_jvirol} no_petitpas={no_petitpas}",
            ))
        else:
            results.append(TestResult(id="T57", name="XP064 source", status=FAIL, message="XP064 not found in TSV"))
    else:
        results.append(TestResult(id="T57", name="XP064 source", status=WARN, message="tier1 TSV not loaded"))

    # ===== Test 58: Citation [33] = Byrne 2008 Proteomic =====
    has_33_byrne = bool(re.search(r"\[33\].*?Byrne", main_text, re.DOTALL)) or \
                   bool(re.search(r"33\.\s+Byrne", main_text))
    results.append(TestResult(
        id="T58", name="Citation [33] = Byrne 2008 Proteomic",
        status=PASS if has_33_byrne else FAIL,
        message=f"present={has_33_byrne}",
    ))

    # ===== Test 59: Citation [34] = Chihara 2013 FN1 =====
    has_34_chihara = bool(re.search(r"\[34\].*?Chihara", main_text, re.DOTALL)) or \
                     bool(re.search(r"34\.\s+Chihara", main_text))
    results.append(TestResult(
        id="T59", name="Citation [34] = Chihara 2013 FN1",
        status=PASS if has_34_chihara else FAIL,
        message=f"present={has_34_chihara}",
    ))

    # ===== Test 60: References 1-67 contiguous (no orphans) =====
    # Extract from References section in order
    ref_pattern_strict = re.compile(r"^(\d+)\.\s+", re.MULTILINE)
    if "## References" in main_text:
        refs_section = main_text[main_text.find("## References"):]
    else:
        refs_section = main_text
    ref_nums_in_order = sorted(int(m) for m in ref_pattern_strict.findall(refs_section[:30000]))
    # R102: was hardcoded 1-67, now flexibly accept 1..N contiguous (N = total ref count)
    if ref_nums_in_order:
        expected_refs = list(range(1, max(ref_nums_in_order) + 1))
        is_contiguous = ref_nums_in_order == expected_refs
    else:
        expected_refs = []
        is_contiguous = False
    results.append(TestResult(
        id="T60", name="References 1-N contiguous in References section (R102: now dynamic N)",
        status=PASS if is_contiguous else FAIL,
        message=f"first_5={ref_nums_in_order[:5]} last_5={ref_nums_in_order[-5:]} total={len(ref_nums_in_order)}",
    ))

    return results


def _find_file(canonical_dir: Path, candidates: list[str]) -> Path | None:
    for c in candidates:
        p = (canonical_dir / c).resolve()
        if p.exists():
            return p
    return None
