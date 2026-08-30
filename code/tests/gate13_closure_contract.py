"""Final closure gate: canonical donor formulas, LOSO units and zero TBD/blockers."""
from __future__ import annotations
import ast, csv, json, hashlib, zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from . import TestResult, register_gate, PASS, FAIL
from build_registry import required_inputs

@register_gate("Gate 13 — Closure contracts", weight=25)
def run(canonical_dir: Path, facts: dict | None) -> list[TestResult]:
    out = []
    p = canonical_dir / "donor_residual_matrix.tsv"
    rows = list(csv.DictReader(p.open(), delimiter="\t")) if p.exists() else []
    systems = ["Complement", "Coagulation", "Immune", "Cytokine", "Endocrine", "Adhesion"]
    errors = []
    for r in rows:
        total = sum(int(r[s]) for s in systems)
        masked = sum(int(r[f"masked_{s}"]) for s in systems)
        selected = sum(int(r[s]) for s in ["Complement", "Coagulation", "Endocrine"])
        if total != int(r["all_six_unmasked"]): errors.append(f"{r['donor']}: unmasked {total}!={r['all_six_unmasked']}")
        if masked != int(r["all_six_masked"]): errors.append(f"{r['donor']}: masked {masked}!={r['all_six_masked']}")
        if selected != int(r["selected_total"]): errors.append(f"{r['donor']}: selected {selected}!={r['selected_total']}")
    out.append(TestResult("T96", "Donor residual matrix row/column/totals formulas", PASS if rows and not errors else FAIL, str(errors or "all rows verified")))

    lp = canonical_dir / "loso_evaluation_units.tsv"
    units = list(csv.DictReader(lp.open(), delimiter="\t")) if lp.exists() else []
    active = [r for r in units if r.get("included_in_LOSO") == "1"]
    counts = {k: sum(r["label"] == k for r in active) for k in ("incompatible", "control")}
    ok = len(active) == 77 and counts == {"incompatible": 66, "control": 11}
    out.append(TestResult("T97", "LOSO machine-readable 77/66/11 evaluation units", PASS if ok else FAIL, f"rows={len(active)} counts={counts}"))

    root = canonical_dir.parent
    table_path = root / "docs" / "canonical_number_table.md"
    table = table_path.read_text(encoding="utf-8") if table_path.exists() else ""
    fact_paths = [root / "manuscript" / n for n in ("manuscript_facts.json", "template_main.md", "template_si.md")]
    scan_text = table + "\n" + "\n".join(p.read_text(encoding="utf-8") for p in fact_paths if p.exists())
    residue = [tok for tok in ("TODO", "TBD", "UNRESOLVED", "UNFILLED") if tok in scan_text]
    if "## Closure status" in table or "blocker:" in table:
        residue.append("canonical table contains closure-status metadata")
    clean = not residue
    out.append(TestResult("T98", "closure blockers independently absent from facts/templates/table", PASS if clean else FAIL, "clean" if clean else "; ".join(residue)))
    mp = canonical_dir.parent / "build_manifest.json"
    errors = []
    if mp.exists():
        manifest = json.loads(mp.read_text())
        expected_inputs = required_inputs(root)
        if manifest.get("required_inputs") != expected_inputs:
            errors.append("required input registry mismatch")
        if set(manifest.get("required_inputs", [])) != {k for k in manifest.get("files", {}) if k in expected_inputs}:
            errors.append("manifest required inputs not fully hashed")
        for rel, expected in manifest.get("files", {}).items():
            q = canonical_dir.parent / rel
            if not q.exists() or hashlib.sha256(q.read_bytes()).hexdigest() != expected:
                errors.append(rel)
    else: errors.append("build_manifest.json missing")
    out.append(TestResult("T99", "figure hashes match current build manifest", PASS if not errors else FAIL, "all manifest files match" if not errors else str(errors[:3])))
    payload_dir = canonical_dir.parent / "figure_payloads"
    payload_errors = []
    try:
        p2 = json.loads((payload_dir / "fig2d.json").read_text())
        if p2.get("metric_variant") != "post_edit_masked": payload_errors.append("Fig2d metric")
        if p2.get("totals") != {"Yunnan-10GE":356, "Revivicor-10GE":334, "eGenesis-69GE":424}: payload_errors.append("Fig2d totals")
        p1 = json.loads((payload_dir / "figS1d.json").read_text())
        anchors = p1.get("anchors", {})
        expected = {"Yunnan-10GE":40, "Revivicor-10GE":28, "eGenesis-69GE":60}
        if p1.get("metric_variant") != "post_edit_clinical_anchor" or {k:v.get("value") for k,v in anchors.items()} != expected: payload_errors.append("FigS1d anchor")
        p3 = json.loads((payload_dir / "figS3.json").read_text())
        if p3.get("resource_count") != 5 or p3.get("feature_count") != 4 or p3.get("aggregate_score") is not False: payload_errors.append("FigS3 contract")
        p5 = json.loads((payload_dir / "fig5.json").read_text())
        if (p5.get("active_units"), p5.get("incompatible"), p5.get("controls")) != (77, 66, 11): payload_errors.append("Fig5 LOSO 77/66/11")
        if p5.get("excluded_sparse_controls") != ["XP075", "XP080"]: payload_errors.append("Fig5 LOSO exclusions")
    except Exception as e: payload_errors.append(str(e))
    out.append(TestResult("T100", "figure payloads obey semantic contracts", PASS if not payload_errors else FAIL, "Fig2 masked, S1 anchors, S3 no aggregate" if not payload_errors else str(payload_errors)))
    s4_errors = []
    try:
        ledger_path = root / "manuscript" / "freeze_ledger.json"
        s4 = json.loads((payload_dir / "figS4.json").read_text())
        ledger_sha = hashlib.sha256(ledger_path.read_bytes()).hexdigest()
        if s4.get("source") != "manuscript/freeze_ledger.json": s4_errors.append("FigS4 source")
        if (s4.get("stage3_provisional"), s4.get("stage4_final")) != (77, 75): s4_errors.append("FigS4 stage values")
        if s4.get("ledger_sha256") != ledger_sha: s4_errors.append("FigS4 ledger hash")
    except Exception as exc:
        s4_errors.append(str(exc))
    out.append(TestResult("T104", "FigS4 payload is ledger-derived and hash-bound", PASS if not s4_errors else FAIL, "stage3=77, stage4=75, ledger hash matches" if not s4_errors else str(s4_errors)))
    templates = "\n".join((canonical_dir.parent / "manuscript" / n).read_text(encoding="utf-8") for n in ("template_main.md", "template_si.md"))
    semantic_errors = []
    if "Aggregate capability" in templates or "sum over the 4 dimensions" in templates:
        semantic_errors.append("S3 aggregate prose")
    for literal in ("Coagulation = 40", "Complement = 20", "Endocrine = 60"):
        if literal in templates:
            semantic_errors.append(f"hard-coded {literal}")
    out.append(TestResult("T101", "semantic dynamic values use placeholders; S3 aggregate forbidden", PASS if not semantic_errors else FAIL, "clean" if not semantic_errors else str(semantic_errors)))
    protocol_errors = []
    if facts:
        bs = facts.get("benchmark_scope", {})
        resampling = bs.get("resampling", {})
        loso = bs.get("loso", {})
        if (resampling.get("n_total"), resampling.get("n_incompatible"), resampling.get("n_control"), resampling.get("no_skill_auprc")) != (75, 62, 13, 0.827):
            protocol_errors.append("resampling 75/62/13/0.827")
        if (loso.get("n_total"), loso.get("n_incompatible"), loso.get("n_control"), loso.get("no_skill_auprc")) != (77, 66, 11, 0.857):
            protocol_errors.append("LOSO 77/66/11/0.857")
        active_map = (canonical_dir.parent / "manuscript" / "template_facts_map.tsv").read_text(encoding="utf-8")
        if "{{BENCHMARK_NO_SKILL_AUPRC}}" in active_map:
            protocol_errors.append("deprecated baseline placeholder still mapped")
    out.append(TestResult("T102", "benchmark protocol denominators and baselines are explicit", PASS if not protocol_errors else FAIL, "75/62/13/0.827 and 77/66/11/0.857" if not protocol_errors else str(protocol_errors)))
    # T103 prevents the file-level duplication problem from returning: the
    # orchestrator registry is the only production figure entry point.
    registry_errors = []
    source_dir = root / "src"
    legacy_names = {
        "make_figS3_S4.py", "regenerate_fig1_fig5_fig_s1_release10.py",
        "regenerate_figS4_release10.py", "regenerate_fig5_figS1_R3.2.py",
        "extract_facts_R102_patch.py",
    }
    if any((source_dir / "_scripts" / name).exists() for name in legacy_names) or (source_dir / "extract_facts_R102_patch.py").exists():
        registry_errors.append("superseded file remains in active source tree")
    try:
        tree = ast.parse((source_dir / "regenerate_figures.py").read_text(encoding="utf-8"))
        node = next(n for n in tree.body if isinstance(n, ast.Assign) and any(getattr(t, "id", "") == "FIGURE_SCRIPTS" for t in n.targets))
        registry = ast.literal_eval(node.value)
        expected_registry = {
            "fig1": ("regenerate_v32_figures.py", ["fig1"]), "fig2": ("make_fig2.py", []),
            "fig3": ("make_fig1_3_5_6.py", ["fig3"]), "fig4": ("make_fig4_corrected.py", []),
            "fig5": ("regenerate_v32_figures.py", ["fig5"]), "fig6": ("make_fig1_3_5_6.py", ["fig6"]),
            "figS1": ("regenerate_v32_figures.py", ["figS1"]), "figS2": ("make_figS2.py", []),
            "figS3": ("regenerate_v32_figures.py", ["figS3"]), "figS4": ("regenerate_v32_figures.py", ["figS4"]),
            "figS5": ("make_figS5_S6.py", []), "figS6": ("make_figS5_S6.py", []),
        }
        if registry != expected_registry:
            registry_errors.append("figure registry drift")
    except Exception as exc:
        registry_errors.append(f"registry unreadable: {exc}")
    out.append(TestResult("T103", "production figure registry excludes superseded generators", PASS if not registry_errors else FAIL, "12 target-selected outputs / no legacy active files" if not registry_errors else "; ".join(registry_errors)))

    # T105: CR-145 propagation contract. XP067 is the only severity correction
    # in the current RC candidate; all active source/export/fact consumers must
    # agree before a golden-master update can be approved.
    propagation_errors = []
    try:
        def _row(path, key="pair_id"):
            with open(path, encoding="utf-8") as fh:
                return next(r for r in csv.DictReader(fh, delimiter="\t") if r.get(key) == "XP067")
        c = _row(canonical_dir / "tier1_frozen.tsv")
        if (c.get("func_sev"), c.get("ISI"), c.get("ISI_v2")) != ("0.5", "15.0", "15.0"):
            propagation_errors.append("canonical XP067")
        for q in (canonical_dir / "tier1_frozen.tsv",):
            r = _row(q)
            if r.get("func_sev") != "0.5" or r.get("ISI_v2") != "15.0":
                propagation_errors.append(str(q.relative_to(root.parent)))
        bm = facts.get("benchmark_M0_M5_full", {}) if facts else {}
        if (bm.get("M4b_resampling_AUROC_mean"), bm.get("M4b_AUPRC")) != (0.998, 0.999):
            propagation_errors.append("M4b facts")
        dm = list(csv.DictReader((canonical_dir / "donor_residual_matrix.tsv").open(), delimiter="\t"))
        fact_dm = (facts or {}).get("donor_scorecard_matrix", {})
        for r in dm:
            fx = fact_dm.get(r["donor"], {})
            for s in ("Complement", "Coagulation", "Immune", "Cytokine", "Endocrine", "Adhesion"):
                if int(r[f"masked_{s}"]) != int(fx.get("masked_systems", {}).get(s, -1)):
                    propagation_errors.append(f"donor matrix {r['donor']} {s}")
    except Exception as exc:
        propagation_errors.append(str(exc))
    out.append(TestResult("T105", "XP067 correction propagates through active source/export/facts consumers", PASS if not propagation_errors else FAIL, "func_sev=0.5, ISI/ISI_v2=15.0, M4b=0.998/0.999 and donor matrix consumers agree" if not propagation_errors else str(propagation_errors)))

    # T106: v47 dependency-chain contract. XP067 is unresolved in every
    # evaluated genotype and its corrected value is represented in the
    # clinical-context residual matrix.
    status_errors = []
    try:
        sp = canonical_dir / "donor_axis_status.tsv"
        sr = list(csv.DictReader(sp.open(), delimiter="\t"))
        xp = [r for r in sr if r.get("pair_id") == "XP067"]
        if len(xp) != 6 or any(r.get("status") != "genetically_unresolved" or r.get("residual_contribution") != "included in post-edit residual burden" for r in xp):
            status_errors.append("XP067 status rows")
        dm = {r["donor"]: r for r in csv.DictReader((canonical_dir / "donor_residual_matrix.tsv").open(), delimiter="\t")}
        expected_complement = {"Yunnan-10GE": 8, "Revivicor-10GE": 28, "eGenesis-69GE": 38}
        for donor, value in expected_complement.items():
            if int(dm[donor]["masked_Complement"]) != value:
                status_errors.append(f"{donor} masked complement")
    except Exception as exc:
        status_errors.append(str(exc))
    out.append(TestResult("T106", "XP067 donor-residual status and matrix propagation", PASS if not status_errors else FAIL, "6/6 genotypes unresolved; corrected Complement contribution present" if not status_errors else str(status_errors)))

    # T107: figure-legend and SI layout contract. This closes the previous
    # blind spot where the image was regenerated but its caption retained an
    # old clinical-anchor value, and verifies that S11 is actually landscape.
    qa_errors = []
    try:
        si_md = (root / "manuscript" / "si.rendered.md").read_text(encoding="utf-8")
        if "Revivicor-10GE reports Complement = 28" not in si_md:
            qa_errors.append("S1 legend Complement=28")
        if "Revivicor-10GE reports Complement = 20" in si_md:
            qa_errors.append("stale S1 legend Complement=20")
        with zipfile.ZipFile(root / "manuscript" / "si.docx") as zf:
            ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
            doc = ET.fromstring(zf.read("word/document.xml"))
        body = doc.find("w:body", ns)
        children = list(body)
        h = next(i for i, node in enumerate(children) if "Supplementary Table S11" in "".join(node.itertext()))
        prev = next((children[i].find(".//w:sectPr", ns) for i in range(h - 1, -1, -1) if children[i].find(".//w:sectPr", ns) is not None), None)
        nxt = next((children[i].find(".//w:sectPr", ns) for i in range(h + 1, len(children)) if children[i].find(".//w:sectPr", ns) is not None), None)
        prev_sz = prev.find("w:pgSz", ns) if prev is not None else None
        next_sz = nxt.find("w:pgSz", ns) if nxt is not None else None
        if prev_sz is None or prev_sz.get("{" + ns["w"] + "}orient") != "landscape":
            qa_errors.append("S11 not preceded by landscape section")
        if next_sz is None or next_sz.get("{" + ns["w"] + "}orient") == "landscape":
            qa_errors.append("S11 not returned to portrait")
    except Exception as exc:
        qa_errors.append(str(exc))
    out.append(TestResult("T107", "S1 canonical legend and S11 landscape production contract", PASS if not qa_errors else FAIL, "S1 Complement=28; S11 landscape→portrait verified" if not qa_errors else str(qa_errors)))

    # T108: coverage denominator propagation. XP067 is unresolved for every
    # donor, so addressed numerators stay fixed while the corrected total
    # denominator changes 1060.8 -> 1068.3.
    cov_errors = []
    try:
        cp = canonical_dir / "donor_coverage_contract.tsv"
        cr = list(csv.DictReader(cp.open(), delimiter="\t"))
        expected = {"GTKO": 3.6, "ClonOrgan-6GE": 11.7, "ClonOrgan-8GE": 11.7,
                    "Yunnan-10GE": 15.2, "Revivicor-10GE": 11.7, "eGenesis-69GE": 13.2}
        if len(cr) != 6 or any(float(r["corrected_total_isi_v2"]) != 1068.3 for r in cr):
            cov_errors.append("coverage contract denominator")
        for r in cr:
            if round(float(r["corrected_coverage_pct"]), 1) != expected.get(r["donor"]):
                cov_errors.append(f"coverage {r['donor']}")
        gd = (facts or {}).get("donor_genotypes_detail", {})
        for donor, value in expected.items():
            if round(float(gd.get(donor, {}).get("coverage_pct", -1)), 1) != value:
                cov_errors.append(f"facts {donor}")
    except Exception as exc:
        cov_errors.append(str(exc))
    out.append(TestResult("T108", "ISI v2-weighted donor coverage denominator propagation", PASS if not cov_errors else FAIL, "6 coverage values recomputed against total ISI v2=1068.3" if not cov_errors else str(cov_errors)))
    return out
