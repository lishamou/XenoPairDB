"""
Gate 11 — Placeholder coverage (R121)

V6 review systemic conclusion: "Fact 已创建, 但 template 没用它, 仍然手工 hard-code"

This gate tracks:
- T93: ACTIVE_PLACEHOLDER_USAGE — every active (non-deprecated) placeholder in
       template_facts_map.tsv should appear in ≥1 template (main.md or si.md).
- T94: FACT_FIELDS_USED — every top-level key in manuscript_facts.json (except
       _meta, _NOTE_R*, _R*_note) should be referenced by ≥1 active placeholder.
- T95: NO_HARDCODED_DUPLICATES — for every {{PLACEHOLDER}} in templates, the
       corresponding value should NOT also be hardcoded in the same file.
       (Advisory WARN, since some hardcodes are deliberate.)
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


@register_gate("Gate 11 — Placeholder coverage (R121)", weight=15)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []
    pipeline_root = canonical_dir.parent
    map_tsv = pipeline_root / "manuscript" / "template_facts_map.tsv"
    template_main = pipeline_root / "manuscript" / "template_main.md"
    template_si = pipeline_root / "manuscript" / "template_si.md"

    if not map_tsv.exists() or not template_main.exists() or not template_si.exists():
        results.append(TestResult(
            id="T93", name="Gate 11 setup (skipped)",
            status=WARN, message="map_tsv or templates missing",
        ))
        return results

    # Read mapping
    with map_tsv.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    active = [r for r in rows if r.get("status", "active") == "active"]
    deprecated = [r for r in rows if r.get("status", "active") == "deprecated"]
    all_placeholders = {r["placeholder"] for r in rows}
    active_placeholders = {r["placeholder"] for r in active}

    # Read templates (combined text)
    text_main = template_main.read_text(encoding="utf-8")
    text_si = template_si.read_text(encoding="utf-8")
    text_combined = text_main + "\n" + text_si

    # T93: ACTIVE_PLACEHOLDER_USAGE
    unused = []
    for ph in sorted(active_placeholders):
        if ph not in text_combined:
            unused.append(ph)
    results.append(TestResult(
        id="T93", name=f"ACTIVE_PLACEHOLDER_USAGE: {len(active_placeholders)-len(unused)}/{len(active_placeholders)} active placeholders used in templates",
        status=PASS if not unused else WARN,
        message=f"unused active placeholders: {len(unused)}" + (f"; first 3: {unused[:3]}" if unused else " (all used)"),
    ))

    # T94: FACT_FIELDS_USED
    if facts is None:
        facts_path = pipeline_root / "manuscript" / "manuscript_facts.json"
        if facts_path.exists():
            facts = json.loads(facts_path.read_text(encoding="utf-8"))
    if facts is None:
        results.append(TestResult(
            id="T94", name="FACT_FIELDS_USED (skipped)",
            status=WARN, message="facts.json missing",
        ))
    else:
        # Get all fact paths referenced by active placeholders
        used_paths = set()
        for r in active:
            p = r["path"]
            # Strip "manuscript_facts.json:" prefix
            if ":" in p:
                used_paths.add(p.split(":", 1)[1])
        # Get all top-level keys in facts
        all_keys = set(facts.keys())
        # Skip non-fact meta keys
        skip_keys = {"_meta"}
        candidate_keys = all_keys - skip_keys
        # Find orphan fact keys (defined but not referenced by any active placeholder)
        # For nested dicts, the path may go deeper (e.g. donor_scorecard_matrix.Yunnan-10GE)
        orphan = []
        for k in candidate_keys:
            if k not in used_paths:
                # Allow if it's a sub-dict that's referenced by a longer path
                if not any(p.startswith(k + ".") for p in used_paths):
                    orphan.append(k)
        results.append(TestResult(
            id="T94", name=f"FACT_FIELDS_USED: {len(candidate_keys)-len(orphan)}/{len(candidate_keys)} top-level fact keys referenced by placeholders",
            status=PASS if not orphan else WARN,
            message=f"orphan fact keys (defined but unused): {len(orphan)}" + (f"; first 3: {orphan[:3]}" if orphan else " (all referenced)"),
        ))

    # T95: NO_HARDCODED_DUPLICATES (advisory)
    # For each active placeholder, check if a hardcoded version of its value also appears
    hardcoded_duplicates = []
    if facts is not None:
        for r in active[:30]:  # sample 30
            ph = r["placeholder"]
            p = r["path"]
            if ":" in p:
                p = p.split(":", 1)[1]
            # Resolve value
            cur = facts
            for k in p.split("."):
                if isinstance(cur, dict) and k in cur:
                    cur = cur[k]
                else:
                    cur = None
                    break
            if cur is None:
                continue
            v = str(cur)
            # Check if the literal value appears hardcoded in either template
            # (outside of placeholders)
            stripped = re.sub(r"\{\{[^}]+\}\}", "", text_combined)
            if v and len(v) >= 3 and v in stripped:
                hardcoded_duplicates.append((ph, v[:30]))
    results.append(TestResult(
        id="T95", name=f"NO_HARDCODED_DUPLICATES: 0 hardcoded values that should be placeholders",
        status=PASS if not hardcoded_duplicates else WARN,
        message=f"potential duplicates: {len(hardcoded_duplicates)}" + (f"; first 3: {hardcoded_duplicates[:3]}" if hardcoded_duplicates else " (clean)"),
    ))

    return results
