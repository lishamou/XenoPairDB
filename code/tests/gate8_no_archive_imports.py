"""
Gate 8 — No archive imports (R110-B: pipeline self-contained)

Per 方法学_一致性问题.md line 53:
  "submission pipeline 永远禁止读取 archive。"
  (en_nbt/ and _scripts/ under Kimi_Agent_XenoPairDB方案二执行与核实202608141004/
   are archive; pipeline MUST NOT read them.)

Tests:
  T76: No pipeline .py file references Kimi_Agent_XenoPairDB方案二执行与核实202608141004/ or en_nbt/
       except in change_requests/, _archive/, README notes explaining WHY.
  T77: No canonical_links/ symlink points to Kimi_Agent.../ or another legacy workspace.
  T78: pipeline/src/_scripts/ scripts all use parameterized ROOT (XENOPAIR_PROJECT_ROOT)
"""
from __future__ import annotations

import re
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN

ARCHIVE_DIR = "Kimi_Agent_XenoPairDB方案二执行与核实202608141004"
ARCHIVE_PATTERNS = [
    re.compile(re.escape(ARCHIVE_DIR)),       # bare name
    re.compile(r'"en_nbt"'),                  # en_nbt literal in code
    re.compile(r"Kimi_Agent.*_scripts"),      # legacy _scripts path
]


def _scan_file_for_archive_refs(p: Path) -> list[tuple[int, str]]:
    """Return [(line_no, line_text), ...] for each archive reference in p."""
    hits = []
    try:
        text = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return hits
    for i, ln in enumerate(text.splitlines(), start=1):
        for pat in ARCHIVE_PATTERNS:
            if pat.search(ln):
                # Allow if line is a comment explaining WHY (starts with # or //)
                stripped = ln.strip()
                if stripped.startswith("#") or stripped.startswith("//") or "R110-" in ln or "no longer" in ln.lower():
                    continue
                hits.append((i, ln.rstrip()))
                break
    return hits


@register_gate("Gate 8 — No archive imports (R110-B)", weight=20)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []
    pipeline_root = canonical_dir.parent
    project_root = pipeline_root.parent

    # ===== T76: pipeline src has no archive references =====
    pipeline_src = pipeline_root / "src"
    archive_hits: list[str] = []
    if pipeline_src.exists():
        for py_file in pipeline_src.rglob("*.py"):
            # Skip _archive/ and tests/__pycache__/
            rel = py_file.relative_to(pipeline_src)
            if str(rel).startswith("_archive") or "__pycache__" in str(rel):
                continue
            # Skip gate8 itself: it MUST contain the archive name to scan for it.
            if rel.name == "gate8_no_archive_imports.py":
                continue
            for line_no, line_text in _scan_file_for_archive_refs(py_file):
                archive_hits.append(f"{rel}:{line_no}: {line_text.strip()[:80]}")
    results.append(TestResult(
        id="T76", name="T76: pipeline/src/*.py has no archive references (Kimi_Agent.../ or en_nbt)",
        status=PASS if not archive_hits else FAIL,
        message=f"hits={len(archive_hits)}" + (f"; first={archive_hits[0]}" if archive_hits else " (clean)"),
    ))

    # ===== T77: canonical_links symlinks do not point to legacy workspaces =====
    symlink_bad = []
    canonical_links = pipeline_root / "canonical_links"
    if canonical_links.exists():
        for f in canonical_links.iterdir():
            if f.is_symlink():
                target = str(f.resolve())
                if ARCHIVE_DIR in target:
                    symlink_bad.append(f"{f.name} → {target}")
    results.append(TestResult(
        id="T77", name="T77: canonical_links/ symlinks do NOT point to Kimi_Agent.../",
        status=PASS if not symlink_bad else FAIL,
        message=f"bad_symlinks={len(symlink_bad)}" + (f"; first={symlink_bad[0]}" if symlink_bad else " (all point to pipeline releases)"),
    ))

    # ===== T78: pipeline/src/_scripts/ ROOT is parameterized (XENOPAIR_PROJECT_ROOT) =====
    scripts_dir = pipeline_root / "src" / "_scripts"
    unparameterized = []
    if scripts_dir.exists():
        for py_file in scripts_dir.glob("*.py"):
            text = py_file.read_text(encoding="utf-8", errors="ignore")
            if "XENOPAIR_PROJECT_ROOT" not in text:
                unparameterized.append(py_file.name)
    results.append(TestResult(
        id="T78", name="T78: pipeline/src/_scripts/*.py use XENOPAIR_PROJECT_ROOT (parameterized)",
        status=PASS if not unparameterized else WARN,
        message=f"unparameterized={len(unparameterized)}" + (f"; first={unparameterized[0]}" if unparameterized else " (all parameterized)"),
    ))

    return results
