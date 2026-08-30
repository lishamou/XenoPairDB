"""
Gate 12 — Orphan hardcoded numbers (R124)

Detects numerical claims in template_main.md + template_si.md that should be
replaced by an existing active placeholder (i.e., the value exists in facts
but the template uses the literal value).

Heuristic: for each active placeholder, check if its canonical value appears
as a literal in the template outside {{ }}. Report as WARN (advisory) —
authors should decide whether to wire the placeholder or accept the hardcode
as deliberate (e.g., for narrative flow).
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


@register_gate("Gate 12 — Orphan hardcodes (R124)", weight=10)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []
    if facts is None:
        p = canonical_dir.parent / "manuscript" / "manuscript_facts.json"
        if p.exists():
            facts = json.loads(p.read_text(encoding="utf-8"))

    pipeline_root = canonical_dir.parent
    map_tsv = pipeline_root / "manuscript" / "template_facts_map.tsv"
    template_main = pipeline_root / "manuscript" / "template_main.md"
    template_si = pipeline_root / "manuscript" / "template_si.md"

    if not map_tsv.exists() or not template_main.exists() or not template_si.exists():
        return [TestResult(id="T93", name="Gate 12 setup (skipped)",
                             status=WARN, message="files missing")]

    # Load mapping
    with map_tsv.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    active = [r for r in rows if r.get("status", "active") == "active"]

    def resolve(path):
        cur = facts
        for k in path.split("."):
            if isinstance(cur, dict) and k in cur:
                cur = cur[k]
            else:
                return None
        return cur

    # T93: scan for orphan hardcodes
    text_main = template_main.read_text(encoding="utf-8")
    text_si = template_si.read_text(encoding="utf-8")
    masked_main = re.sub(r"\{\{[^}]+\}\}", "{{PH}}", text_main)
    masked_si = re.sub(r"\{\{[^}]+\}\}", "{{PH}}", text_si)

    orphan_count = 0
    orphan_samples = []
    for r in active:
        p = r["path"]
        if ":" in p:
            p = p.split(":", 1)[1]
        val = resolve(p)
        if val is None or isinstance(val, bool):
            continue
        if isinstance(val, (int, float)):
            v = int(val) if isinstance(val, int) or (isinstance(val, float) and val.is_integer()) else round(val, 4)
            if isinstance(v, int) and v in (0, 1):
                continue  # 0 and 1 are too common to be meaningful
            # Count literal occurrences in masked text
            n_main = len(re.findall(rf"\b{v}\b", masked_main))
            n_si = len(re.findall(rf"\b{v}\b", masked_si))
            if (n_main + n_si) >= 3:  # threshold
                orphan_count += 1
                if len(orphan_samples) < 5:
                    orphan_samples.append(
                        f"{r['placeholder']}={v} (main={n_main},si={n_si})"
                    )

    results.append(TestResult(
        id="T93", name=f"ORPHAN_HARDCODES: {orphan_count} active placeholders have ≥3 literal occurrences in templates",
        status=PASS if orphan_count == 0 else WARN,
        message=f"{'; '.join(orphan_samples) if orphan_samples else 'all canonical values are placeholder-resolved'}",
    ))

    return results
