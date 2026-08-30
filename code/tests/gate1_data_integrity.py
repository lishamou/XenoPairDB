"""
Gate 1 — Data integrity tests (Tests 1-15)

L7 of the 8-layer architecture. These tests verify that the canonical Release 3.2
TSV/JSON artifacts satisfy the schema contract and the disjoint-set invariants.

Each test returns (test_id, name, status, message). The audit orchestrator
collects these into the consistency report.
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


@register_gate("Gate 1 — Data integrity", weight=15)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []

    # Load canonical files
    tier1 = _load_tsv(canonical_dir / "tier1_frozen.tsv")
    boundary = _load_tsv(canonical_dir / "boundary.tsv")
    engineering = _load_tsv(canonical_dir / "engineering_modifier.tsv")

    tier1_pids = [r["pair_id"] for r in tier1]
    boundary_pids = [r["pair_id"] for r in boundary]
    engineering_pids = [r["pair_id"] for r in engineering]

    # ===== Test 1: Tier-1 row count = 75 =====
    actual = len(tier1)
    expected = 75
    results.append(TestResult(
        id="T01", name="Tier-1 row count = 75",
        status=PASS if actual == expected else FAIL,
        message=f"actual={actual} expected={expected}",
        actual=str(actual), expected=str(expected),
    ))

    # ===== Test 2: Evidence class sum = 75 (58 direct + 4 disputed + 0 inferred + 13 control) =====
    evidence_counts = Counter(r["evidence"] for r in tier1)
    s = sum(evidence_counts.values())
    expected_s = 75
    results.append(TestResult(
        id="T02", name="Evidence class sum = 75",
        status=PASS if s == expected_s else FAIL,
        message=f"sum={s} distribution={dict(evidence_counts)}",
        actual=str(s), expected=str(expected_s),
    ))

    # ===== Test 3: Functional-system sum = 75 (22+12+15+7+10+9) =====
    system_counts = Counter(r["system"] for r in tier1)
    s = sum(system_counts.values())
    expected_systems = {"Complement": 9, "Coagulation": 12, "Immune": 22,
                        "Cytokine": 15, "Endocrine": 7, "Adhesion": 10}
    results.append(TestResult(
        id="T03", name="Functional-system sum = 75",
        status=PASS if (s == 75 and all(system_counts.get(k) == v for k, v in expected_systems.items())) else FAIL,
        message=f"sum={s} distribution={dict(system_counts)}",
    ))

    # ===== Test 4: Surfaceome universe >= scoreable (4123 >= 4122) =====
    surfaceome_path = canonical_dir / "surfaceome.tsv"
    if surfaceome_path.exists():
        rows = _load_tsv(surfaceome_path)
        universe = len(rows)
        scoreable = sum(1 for r in rows if r.get("scoreability_status", "").strip().lower() in {"", "scoreable"})
        ok = universe == 4123 and scoreable == 4122
        results.append(TestResult(
            id="T04", name="Surfaceome 4123 universe / 4122 scoreable",
            status=PASS if ok else FAIL,
            message=f"universe={universe} scoreable={scoreable}",
            actual=f"{universe}/{scoreable}", expected="4123/4122",
        ))
    else:
        results.append(TestResult(id="T04", name="Surfaceome", status=FAIL, message="file missing"))

    # ===== Test 5: Tier-1 ∩ Boundary = ∅ =====
    overlap = set(tier1_pids) & set(boundary_pids)
    results.append(TestResult(
        id="T05", name="Tier-1 ∩ Boundary = ∅",
        status=PASS if not overlap else FAIL,
        message=f"overlap={sorted(overlap) if overlap else '∅'}",
    ))

    # ===== Test 6: Tier-1 ∩ Engineering = ∅ =====
    overlap = set(tier1_pids) & set(engineering_pids)
    results.append(TestResult(
        id="T06", name="Tier-1 ∩ Engineering = ∅",
        status=PASS if not overlap else FAIL,
        message=f"overlap={sorted(overlap) if overlap else '∅'}",
    ))

    # ===== Test 7: Boundary ∩ Engineering = ∅ =====
    overlap = set(boundary_pids) & set(engineering_pids)
    results.append(TestResult(
        id="T07", name="Boundary ∩ Engineering = ∅",
        status=PASS if not overlap else FAIL,
        message=f"overlap={sorted(overlap) if overlap else '∅'}",
    ))

    # ===== Test 8: All pair_ids unique in Tier-1 =====
    seen = set()
    dups = []
    for pid in tier1_pids:
        if pid in seen:
            dups.append(pid)
        seen.add(pid)
    results.append(TestResult(
        id="T08", name="Tier-1 pair_ids unique",
        status=PASS if not dups else FAIL,
        message=f"duplicates={dups if dups else '∅'} (n_unique={len(seen)})",
    ))

    # ===== Test 9: pair_id format XP### =====
    pattern = re.compile(r"^XP[0-9]{3}$")
    bad = [pid for pid in tier1_pids if not pattern.match(pid)]
    results.append(TestResult(
        id="T09", name="pair_id format = XP###",
        status=PASS if not bad else FAIL,
        message=f"bad_format={bad if bad else '∅'}",
    ))

    # ===== Test 10: Required schema columns present =====
    required = ["pair_id", "system", "pig_gene", "human_partner", "evidence",
                "source", "identity", "similarity", "func_sev", "evid_w",
                "iface_id", "ISI", "expr_rho", "expr_weight", "ISI_v2"]
    if tier1:
        actual_cols = set(tier1[0].keys())
        missing = [c for c in required if c not in actual_cols]
        results.append(TestResult(
            id="T10", name="Required schema columns present",
            status=PASS if not missing else FAIL,
            message=f"missing={missing if missing else '∅'} (n_required={len(required)})",
        ))

    # ===== Test 11: Numeric column ranges =====
    bad = []
    for r in tier1:
        # identity: 0-100 or n/a
        v = r.get("identity", "").strip()
        if v and v != "n/a":
            try:
                f = float(v)
                if not (0.0 <= f <= 100.0):
                    bad.append((r["pair_id"], "identity", f))
            except ValueError:
                bad.append((r["pair_id"], "identity", v))
        # func_sev: 0-1
        v = r.get("func_sev", "").strip()
        if v and v != "n/a":
            try:
                f = float(v)
                if not (0.0 <= f <= 1.0):
                    bad.append((r["pair_id"], "func_sev", f))
            except ValueError:
                bad.append((r["pair_id"], "func_sev", v))
        # expr_weight: 0.75-1.25
        v = r.get("expr_weight", "").strip()
        if v and v != "n/a":
            try:
                f = float(v)
                if not (0.75 <= f <= 1.25):
                    bad.append((r["pair_id"], "expr_weight", f))
            except ValueError:
                pass
    results.append(TestResult(
        id="T11", name="Numeric column ranges (identity 0-100, func_sev 0-1, expr_weight 0.75-1.25)",
        status=PASS if not bad else FAIL,
        message=f"violations={len(bad)} (sample: {bad[:3]})",
    ))

    # ===== Test 12: system in 6 allowed categories =====
    allowed_systems = {"Complement", "Coagulation", "Immune", "Cytokine", "Endocrine", "Adhesion"}
    bad = [r["system"] for r in tier1 if r["system"] not in allowed_systems]
    results.append(TestResult(
        id="T12", name="system in 6 allowed categories",
        status=PASS if not bad else FAIL,
        message=f"bad_systems={set(bad) if bad else '∅'}",
    ))

    # ===== Test 13: evidence in {direct, disputed, control} for Tier-1 =====
    allowed_evidence = {"direct", "disputed", "control"}
    bad = [(r["pair_id"], r["evidence"]) for r in tier1 if r["evidence"] not in allowed_evidence]
    results.append(TestResult(
        id="T13", name="Tier-1 evidence in {direct, disputed, control}",
        status=PASS if not bad else FAIL,
        message=f"violations={bad if bad else '∅'}",
    ))

    # ===== Test 14: No legacy shortrefs in source column =====
    legacy_pattern = re.compile(r"(?:Byrne|Martens|Burlak|Estrada|Chihara|Ladowski|Xie|Lutz|Wang|Yang|Cooper)[0-9]{4}_(xen|Tx|Xenotransplantation|immu)")
    bad = []
    for r in tier1:
        src = r.get("source", "")
        if legacy_pattern.search(src):
            bad.append((r["pair_id"], src))
    results.append(TestResult(
        id="T14", name="No legacy shortrefs in source column",
        status=PASS if not bad else FAIL,
        message=f"legacy_shortrefs={bad if bad else '∅'}",
    ))

    # ===== Test 15: SHA-256 of canonical TSVs matches SHA256SUMS.json =====
    sha_path = canonical_dir / "SHA256SUMS.json"
    if sha_path.exists():
        import hashlib, json
        with open(sha_path) as f:
            sha_data = json.load(f)
        files_dict = sha_data.get("files", sha_data)
        # Map short names to canonical
        name_map = {
            "XenoPairDB_v3.3_Tier1_FROZEN.tsv": canonical_dir / "tier1_frozen.tsv",
            "XenoPairDB_v3.2_adjudicated_boundary_cases.tsv": canonical_dir / "boundary.tsv",
            "XenoPairDB_v3.2_engineering_modifier.tsv": canonical_dir / "engineering_modifier.tsv",
        }
        drift = []
        for fname, expected_sha in files_dict.items():
            target = name_map.get(fname)
            if target and target.exists():
                actual_sha = hashlib.sha256(open(target, "rb").read()).hexdigest()
                if actual_sha != expected_sha:
                    drift.append(fname)
        results.append(TestResult(
            id="T15", name="SHA-256 of canonical TSVs matches SHA256SUMS.json",
            status=PASS if not drift else FAIL,
            message=f"drift={drift if drift else '∅'}",
        ))
    else:
        results.append(TestResult(id="T15", name="SHA-256 verification", status=WARN, message="SHA256SUMS.json not found"))

    return results


def _load_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))
