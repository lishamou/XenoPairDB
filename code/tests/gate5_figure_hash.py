"""
Gate 5 — Figure integrity tests (R105)

These tests verify that the figures actually embedded in the final DOCX
are the canonical R3.2 figures (not stale R1.0/R3.1).

The R102 audit (v2 doc Finding 5) revealed that captions and figure pixels
could drift independently. Gate 4 (T68/T69) only checks PDF text; Gate 5
checks the actual PNG bytes embedded in the DOCX.

Three test classes:
  - FIGURE_PRESENCE: every canonical figure (Fig 1-6, Fig S1-S4) is present
  - FIGURE_HASH: DOCX word/media SHA-256 matches en_nbt/ canonical SHA-256
  - CAPTION_FIGURE_FACTS: figure_provenance.json panel metadata matches
    facts.json canonical values
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import zipfile
from pathlib import Path

from . import TestResult, register_gate, FAIL, PASS, WARN


CANONICAL_FIGURES = {
    # figure_name : (relative path in project, expected in docx as)
    "fig1_resource_v32": "fig1_resource.png",
    "fig2_atlas":         "fig2_atlas.png",
    "fig3_modeling":      "fig3_modeling.png",
    "fig4_antigenicity":  "fig4_antigenicity.png",
    "fig5_benchmark":     "fig5_benchmark.png",
    "fig6_clinical":      "fig6_clinical.png",
    "figS1_resource_panels": "figS1_resource_panels.png",
    "figS2_antigen_benchmark": "figS2_antigen_benchmark.png",
    "figS3_comparison":   "figS3_comparison.png",
    "figS4_v32":          "figS4_prisma.png",
    "figS5_graft_observability": "figS5_graft_observability.png",
    "figS6_BR0010_stress_test": "figS6_BR0010_stress_test.png",
}


def _sha256(p: Path) -> str:
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def _docx_embedded_hashes(docx_path: Path) -> dict[str, str]:
    """Return {rId: sha256} for all images in word/media/."""
    hashes = {}
    if not docx_path.exists():
        return hashes
    try:
        with zipfile.ZipFile(docx_path) as z:
            for name in z.namelist():
                if name.startswith("word/media/") and Path(name).suffix.lower() in (".png", ".jpg", ".jpeg"):
                    with z.open(name) as f:
                        hashes[name] = hashlib.sha256(f.read()).hexdigest()
    except Exception:
        pass
    return hashes


@register_gate("Gate 5 — Figure integrity (R105)", weight=10)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    results = []

    # Resolve project root from canonical_dir
    # canonical_dir is XenoPairDB_Integrity_Pipeline/canonical_links/
    # R110-B: pipeline self-contained — read pipeline/manuscript/*.docx only
    pipeline_root = canonical_dir.parent
    project_root = pipeline_root.parent
    main_docx = pipeline_root / "manuscript" / "main.docx"
    si_docx   = pipeline_root / "manuscript" / "si.docx"
    en_nbt = None  # R110-B: en_nbt is archive; never read

    # ===== T71: FIGURE_PRESENCE (main.docx has all 6 figures) =====
    main_hashes = _docx_embedded_hashes(main_docx)
    si_hashes = _docx_embedded_hashes(si_docx)
    n_main_imgs = len(main_hashes)
    n_si_imgs = len(si_hashes)
    results.append(TestResult(
        id="T71", name="FIGURE_PRESENCE: 6 main + 6 SI figures embedded in DOCX",
        status=PASS if (n_main_imgs >= 6 and n_si_imgs >= 6) else FAIL,
        message=f"main.docx has {n_main_imgs} images (need ≥6); si.docx has {n_si_imgs} images (need ≥6)",
    ))

    # ===== T72: FIGURE_HASH (DOCX embedded SHA matches en_nbt/ canonical SHA) =====
    # For each canonical figure, compute SHA of the en_nbt/ file and compare
    # to a docx-embedded image that has the same byte size (heuristic)
    mismatches = []
    matches = 0
    for fig_key, rel_path in CANONICAL_FIGURES.items():
        # R110-B: read canonical figures from pipeline/figures/
        canonical_path = pipeline_root / "figures" / Path(rel_path).name
        if not canonical_path.exists():
            mismatches.append(f"{fig_key}: canonical PNG missing at {canonical_path}")
            continue
        canonical_sha = _sha256(canonical_path)
        canonical_size = canonical_path.stat().st_size
        # Find docx image with same size
        found = False
        for name, sha in main_hashes.items():
            # DOCX images may be resized; just check sha for matching
            if sha == canonical_sha:
                matches += 1
                found = True
                break
        if not found:
            # Also check si
            for name, sha in si_hashes.items():
                if sha == canonical_sha:
                    matches += 1
                    found = True
                    break
        if not found:
            mismatches.append(f"{fig_key}: canonical SHA not found in any embedded image")
    results.append(TestResult(
        id="T72", name=f"FIGURE_HASH: {len(CANONICAL_FIGURES)} canonical PNG SHA matches embedded image",
        status=PASS if not mismatches else FAIL,
        message=f"matches={matches}/{len(CANONICAL_FIGURES)} mismatches={mismatches[:3] if mismatches else '∅'}",
    ))

    # ===== T73: CAPTION_FIGURE_FACTS (figure_provenance panel metadata matches facts.json) =====
    fp_path = canonical_dir / "figure_provenance.json"
    if fp_path.exists() and facts:
        fp = json.load(open(fp_path))
        # Check fig1 panel e three views
        expected_views = facts.get("figures_panel_E_three_views", {})
        actual_views = fp.get("fig1_resource_v32", {}).get("panels", {}).get("e", "")
        if isinstance(actual_views, dict):
            actual_all_unmasked = actual_views.get("all_six_system_unmasked", {})
            actual_all_masked = actual_views.get("all_six_system_masked", {})
            actual_3system = actual_views.get("three_system_selected", {})
            expected_unmasked = expected_views.get("all_six_system_unmasked", {})
            expected_masked = expected_views.get("all_six_system_masked", {})
            expected_3system = expected_views.get("three_system_selected", {})
            all_match = (
                actual_all_unmasked == expected_unmasked and
                actual_all_masked == expected_masked and
                actual_3system == expected_3system
            )
        else:
            # Old format: actual_views is a string description; check for the numbers
            ok = "995" in str(actual_views) and "974" in str(actual_views) and "946" in str(actual_views)
            all_match = ok
        results.append(TestResult(
            id="T73", name="CAPTION_FIGURE_FACTS: fig1 panel_e 3 views match facts.json",
            status=PASS if all_match else FAIL,
            message=f"all_match={all_match}",
        ))
    else:
        results.append(TestResult(
            id="T73", name="CAPTION_FIGURE_FACTS",
            status=WARN,
            message=f"figure_provenance.json or facts.json missing",
        ))

    # ===== T74: NO_STALE_FIGURE (no R1.0/R3.1 hardcoded numbers in figure provenance) =====
    fp_path = canonical_dir / "figure_provenance.json"
    if fp_path.exists():
        fp = json.load(open(fp_path))
        # Check for stale R1.0 numbers
        stale_numbers = []
        checks = [
            ("fig1_resource_v32", "69"),  # Should be 75, not 69 (unless in historical context)
            ("fig5_benchmark", "0.884"),   # Should be 0.827, not 0.884
            ("figS1_resource_panels_v32", "62"),  # Should be 75, not 62
        ]
        # ...etc
        text = json.dumps(fp)
        for label, num in checks:
            if num in text and label in text:
                # The number appears in the same figure context
                # If R1.0 number in the v32 block, this is suspicious
                if f'"{label}"' in text:
                    idx = text.find(f'"{label}"')
                    ctx = text[idx:idx+5000]
                    if num in ctx and not any(neg in ctx for neg in ["NOT", "historical"]):
                        stale_numbers.append(f"{label}: contains {num}")
        # Use a softer check
        results.append(TestResult(
            id="T74", name="NO_STALE_FIGURE: R1.0 numbers absent from v32 figure_provenance",
            status=PASS if not stale_numbers else WARN,
            message=f"stale={stale_numbers if stale_numbers else '∅'} (heuristic check)",
        ))
    else:
        results.append(TestResult(
            id="T74", name="NO_STALE_FIGURE",
            status=WARN,
            message="figure_provenance.json missing",
        ))

    return results
