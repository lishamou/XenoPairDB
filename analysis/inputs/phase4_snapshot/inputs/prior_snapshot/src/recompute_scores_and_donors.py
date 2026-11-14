"""Replay legacy functions, then compute explicitly provisional scoring/mapping scenarios.

Nothing here replaces the approved scientific scoring contract. Missing scores remain
NA; sums and coverages explicitly refer to the numeric subset. No target coverage,
stored residual, facts or display-rounded coverage feeds the new computation.
"""
from collections import defaultdict
import importlib.util
import math
from common import *

VARIANTS = ["legacy_stored", "formula_sequence_available", "formula_nominal50_sensitivity", "without_divergence_sensitivity"]
SYSTEMS = ["Complement", "Coagulation", "Immune", "Cytokine", "Endocrine", "Adhesion"]

def divergence(identity, interface=None, unit="percent", missing="strict"):
    if unit not in {"percent", "fraction"}:
        raise ValueError("unit must be explicit: percent or fraction")
    scale = 100 if unit == "percent" else 1
    for value in [identity, interface]:
        if value is not None and (not math.isfinite(value) or not 0 <= value <= scale):
            raise ValueError("Divergence input outside declared scale")
    if missing not in {"strict", "nominal50"}:
        raise ValueError("Unknown missing policy")
    if identity is None:
        return 50.0 if missing == "nominal50" else None
    seq = 100 * (1 - identity / scale)
    return seq if interface is None else .5 * seq + .5 * 100 * (1 - interface / scale)

def scores(r):
    e, s, w = [number(r[k]) for k in ["evid_w", "func_sev", "expr_weight"]]
    if e is not None and not 0 <= e <= 1:
        raise ValueError("Evidence weight outside [0,1]")
    if s is not None and not 0 <= s <= 1:
        raise ValueError("Severity outside [0,1]")
    if w is not None and not .75 <= w <= 1.25:
        raise ValueError("Expression weight outside canonical interval")
    base = None if any(x is None for x in [e, s, w]) else e * s * w
    ident, iface = number(r["identity"]), number(r["iface_id"])
    d = divergence(ident, iface, missing="strict")
    n = divergence(ident, iface, missing="nominal50")
    return {"legacy_stored": number(r["ISI_v2"]),
            "formula_sequence_available": None if base is None or d is None else base * d,
            "formula_nominal50_sensitivity": None if base is None else base * n,
            "without_divergence_sensitivity": base}

def aggregate(contributions):
    """Return totals only over numeric rows with known masks; never coerce missing."""
    unique(contributions)
    numeric, unknown = [], []
    for r in contributions:
        s, m = number(r["score"]), number(r["mask"])
        if s is not None and s < 0:
            raise ValueError("Negative burden")
        if m not in {0, 1, None}:
            raise ValueError("Mask must be 0, 1 or explicit NA")
        if s is not None and m is not None:
            numeric.append((s, m))
        elif s is not None:
            unknown.append(s)
    t = math.fsum(s for s, _ in numeric)
    a = math.fsum(s * m for s, m in numeric)
    residual = math.fsum(s * (1 - m) for s, m in numeric)
    return dict(T=t, A=a, R=residual, coverage_fraction=a / t if t else None,
                numeric_known_mask_n=len(numeric), unscored_n=sum(number(r["score"]) is None for r in contributions),
                scored_unknown_mask_n=len(unknown), scored_unknown_mask_burden=math.fsum(unknown),
                full_scope_T=t if len(numeric) == len(contributions) else None)

def direct_masks(edits, axes):
    """Read canonical reported mapping, checking identity before accepting each edge.

    m=0 means no direct reported edge in this scenario, not absence of indirect
    effects. m=1 means genetically addressed under complete-masking assumption,
    not measured functional rescue. Q02 remains open.
    """
    byid = {r["pair_id"]: r for r in axes}
    if len(byid) != len(axes):
        raise ValueError("Duplicate axis ID")
    mappings, mapped = [], defaultdict(set)
    for index, edit in enumerate(edits, 2):
        pid, gene = edit["tier1_pair"], edit["gene"]
        target = byid.get(pid)
        if pid.lower() not in NA and target is None:
            raise ValueError(f"Unknown mapped axis {pid}")
        if edit["edit_type"] not in {"KO", "KI"}:
            raise ValueError("Unknown edit type")
        aliases = {"hEPCR": "PROCR", "hENTPD1": "ENTPD1"}
        symbol = aliases.get(gene, gene[1:] if gene.startswith("h") else gene)
        if target and target["gene_key"] != symbol:
            raise ValueError(f"Gene/entity mismatch: {gene} -> {pid} ({target['gene_key']})")
        if target:
            mapped[edit["genotype"]].add(pid)
        mappings.append(dict(donor=edit["genotype"], edit_type=edit["edit_type"], gene=gene,
                             candidate_pair_id=pid if target else None,
                             stable_axis_id=stable_id(target) if target else None,
                             pig_gene=target["pig_gene"] if target else None,
                             human_partner=target["human_partner"] if target else None,
                             direction=target["direction"] if target else None,
                             functional_endpoint=target["phenotype"] if target else None,
                             mapping_status="REPORTED_DIRECT_EDGE_CANDIDATE_Q02" if target else "MODIFIER_NOT_DISTRIBUTED",
                             source=f"canonical/genotype_definitions.tsv:row{index}",
                             literature_basis_reported=edit["literature_basis"],
                             primary_source_reverified="False", assumption="complete genetic masking for scenario; functional rescue unknown"))
    return mappings, mapped

def replay_legacy(rows):
    # These two source files were inspected: pure functions/constants; no write entry point.
    import pandas as pd
    modules = {}
    for name in ["scoring", "genotypes"]:
        spec = importlib.util.spec_from_file_location("phase1_legacy_" + name,
            SNAP / P / "src/xenopairdb" / (name + ".py"))
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        modules[name] = m
    t = pd.read_csv(canonical("tier1_frozen.tsv"), sep="\t")
    frame = modules["genotypes"].scorecard(t)
    frame.to_csv(REV / "results/legacy_code_replay_scorecard.tsv", sep="\t", index=False)
    modules["genotypes"].residual_isi_by_system(t).to_csv(
        REV / "results/legacy_code_replay_residual_by_system.tsv", sep="\t", index_label="donor")
    legacy_scores = []
    for r in rows:
        e, s, w = (number(r[k]) for k in ["evid_w", "func_sev", "expr_weight"])
        d = modules["scoring"].molecular_divergence(number(r["identity"]), number(r["iface_id"]))
        value = float(modules["scoring"].isi_v2(e, s, d, w))
        legacy_scores.append(dict(pair_id=r["pair_id"], molecular_divergence=d,
                                  replay_formula_ISI_v2=value, stored_ISI_v2=number(r["ISI_v2"]),
                                  replay_matches_new_nominal50=math.isclose(value, scores(r)["formula_nominal50_sensitivity"], abs_tol=1e-8, rel_tol=1e-8)))
    write_tsv(REV / "results/legacy_code_replay_axis_scores.tsv", legacy_scores)
    return modules["genotypes"], frame.to_dict("records")

def main():
    axes = read_tsv(canonical("tier1_frozen.tsv"))
    unique(axes)
    legacy, replay = replay_legacy(axes)
    edits = read_tsv(canonical("genotype_definitions.tsv"))
    mappings, mapped = direct_masks(edits, axes)
    donors = list(dict.fromkeys(e["genotype"] for e in edits))
    write_tsv(REV / "results/edit_axis_mapping.tsv", mappings)
    mappingdiff = []
    for m in mappings:
        oldgene = "GHRHR" if m["gene"] == "GHR" else m["gene"]
        oldpid = legacy.GENE_TO_PAIR.get(oldgene)
        if oldpid != m["candidate_pair_id"] or oldgene != m["gene"]:
            mappingdiff.append(dict(donor=m["donor"], legacy_gene=oldgene, candidate_gene=m["gene"],
                                    legacy_pair_id=oldpid, candidate_pair_id=m["candidate_pair_id"],
                                    status="CANONICAL_REPORTED_MAPPING_CANDIDATE_Q02", basis=m["source"]))
    write_tsv(REV / "results/edit_mapping_diff.tsv", mappingdiff)
    write_tsv(REV / "results/edit_related_axes_pending.tsv", [
        dict(edit="GGTA1 KO", related_axis="XP074", primary_axis="XP004", status="NEEDS_DECISION_Q02",
             reason="XP074 is alphaGal–LGALS3 adhesion, distinct from XP004 antibody/complement. Canonical mapping only lists XP004; no automatic propagation."),
        dict(edit="hCD47 KI", related_axis="XP068", primary_axis="XP020", status="NO_AUTOMATIC_REVERSE_MAPPING",
             reason="XP068 has pig SIRPA/human CD47 species roles and is already control; does not equal XP020"),
        dict(edit="GHR KO", related_axis="XP053", primary_axis="XP052", status="NO_AUTOMATIC_REVERSE_MAPPING",
             reason="XP053 is pig GH1 on human GHR, not donor GHR; no shared-ID or gene-name-only masking"),
        dict(edit="donor ABO phenotype", related_axis="XP079", primary_axis="NA", status="UNKNOWN_APPLICABILITY",
             reason="Axis is phenotype conditional; donor ABO phenotype absent. Kept unchanged, not implicitly absent.")])
    reconciliation, score_rows = [], []
    for r in axes:
        sc = scores(r)
        d = divergence(number(r["identity"]), number(r["iface_id"]))
        base = number(r["evid_w"]) * number(r["func_sev"])
        v1 = base * d if d is not None else None
        stored_v1, stored_v2 = number(r["ISI"]), number(r["ISI_v2"])
        # Staged vs unstaged arithmetic kept separate; never use rounded intermediate for recomputation.
        product_stored_v1 = stored_v1 * number(r["expr_weight"]) if stored_v1 is not None else None
        rho = number(r["expr_rho"])
        rho_weight = min(1.25, max(.75, .75 + .5 * max(0, rho))) if rho is not None else None
        reconciliation.append(dict(pair_id=r["pair_id"], stable_axis_id=stable_id(r), pig_gene=r["pig_gene"],
            human_partner=r["human_partner"], direction=r["direction"], endpoint=r["phenotype"],
            evidence=r["evidence"], source=r["source"], identity_pct=number(r["identity"]), iface_identity_pct=number(r["iface_id"]),
            divergence_pct=d, missing_interface_proxy="sequence_only" if number(r["iface_id"]) is None and d is not None else "no",
            evid_w=number(r["evid_w"]), func_sev=number(r["func_sev"]), expr_weight=number(r["expr_weight"]),
            expr_rho=rho, weight_from_rho_no_override=rho_weight,
            expr_override_or_rounding_flag=(abs(rho_weight - number(r["expr_weight"])) > .00051) if rho_weight is not None else "rho_missing_stored_weight_only",
            stored_ISI=stored_v1, formula_ISI=v1, stored_ISI_v2=stored_v2,
            formula_ISI_v2=sc["formula_sequence_available"], nominal50_ISI_v2=sc["formula_nominal50_sensitivity"],
            without_divergence=sc["without_divergence_sensitivity"], staged_ISI_times_expr=product_stored_v1,
            formula_minus_stored=(sc["formula_sequence_available"] - stored_v2) if sc["formula_sequence_available"] is not None and stored_v2 is not None else None,
            within_single_final_decimal_rounding=(abs(sc["formula_sequence_available"] - stored_v2) <= .05000000001) if sc["formula_sequence_available"] is not None and stored_v2 is not None else None,
            status="XP075_FORMULA_EXCEPTION_Q01" if r["pair_id"] == "XP075" else "MISSING_SEQUENCE_OR_SCORE" if stored_v2 is None or d is None else "COMPARE_UNROUNDED_VS_LEGACY_ROUNDED",
            approval="SCENARIO_ONLY_NOT_OFFICIAL"))
        for variant in VARIANTS:
            score_rows.append(dict(pair_id=r["pair_id"], stable_axis_id=stable_id(r), variant=variant,
                                   score=sc[variant], units="percent_divergence_weighted" if variant != VARIANTS[-1] else "evidence_severity_expression_product",
                                   score_status="MISSING_NOT_ZERO" if sc[variant] is None else "STORED_LEGACY" if variant == VARIANTS[0] else "CANDIDATE_NOT_APPROVED"))
    write_tsv(REV / "results/axis_score_reconciliation.tsv", reconciliation)
    write_tsv(REV / "results/axis_scores_long.tsv", score_rows)
    masks, contributions, cards, systems, missing = [], [], [], [], []
    for donor in donors:
        editrows = [e for e in mappings if e["donor"] == donor]
        for r in axes:
            edges = [e for e in editrows if e["candidate_pair_id"] == r["pair_id"]]
            masks.append(dict(donor=donor, pair_id=r["pair_id"], stable_axis_id=stable_id(r),
                mask=int(r["pair_id"] in mapped[donor]), status="GENETICALLY_ADDRESSED_SCENARIO" if edges else "NO_DIRECT_EDGE_IN_REPORTED_CONTRACT",
                source=";".join(e["source"] for e in edges) if edges else "complete canonical genotype_definitions.tsv; no reported direct edge",
                edits=";".join(e["gene"] for e in edges), functional_rescue="UNKNOWN",
                applicability="UNKNOWN_ABO_PHENOTYPE" if r["pair_id"] == "XP079" else "AS_LEGACY_SCOPE_NOT_RECURATED",
                related_mapping_pending="Q02:GGTA1_also_affects_alphaGal_endpoint?" if r["pair_id"] == "XP074" else "none",
                mapping_contract="reported_direct_only_candidate_Q02"))
        for variant in VARIANTS:
            cohort = []
            for r in axes:
                s = scores(r)[variant]
                mask = int(r["pair_id"] in mapped[donor])
                item = dict(donor=donor, variant=variant, pair_id=r["pair_id"], stable_axis_id=stable_id(r),
                    pig_gene=r["pig_gene"], human_partner=r["human_partner"], direction=r["direction"], endpoint=r["phenotype"],
                    system=r["system"], evidence=r["evidence"], score=s, mask=mask,
                    addressed_contribution=s * mask if s is not None else None,
                    residual_contribution=s * (1 - mask) if s is not None else None,
                    score_status="MISSING_NOT_ZERO" if s is None else "NUMERIC",
                    source="snapshot:canonical/tier1_frozen.tsv:" + r["pair_id"],
                    mask_source="results/axis_donor_masks.tsv", status="CANDIDATE_NOT_OFFICIAL")
                cohort.append(item)
                contributions.append(item)
                if s is None:
                    missing.append(dict(donor=donor, variant=variant, pair_id=r["pair_id"], mask=mask,
                                        identity_status="missing" if number(r["identity"]) is None else "available",
                                        missing_score_burden=None, status="NOT_ASSUMED_ZERO_COMPATIBLE_OR_RESOLVED"))
            nc = [r for r in cohort if r["evidence"] != "control"]
            card = dict(donor=donor, variant=variant, mapping_contract="reported_direct_only_candidate_Q02",
                resource_axes=len(axes), non_control_axes=len(nc),
                addressed_resource_axes=len(mapped[donor]), addressed_non_control_axes=sum(r["mask"] for r in nc),
                axis_count_coverage_fraction=sum(r["mask"] for r in nc) / len(nc) if nc else None,
                addressed_unscored_non_control_axes=sum(r["mask"] and r["score"] is None for r in nc),
                unaddressed_unscored_non_control_axes=sum(not r["mask"] and r["score"] is None for r in nc),
                score_scope="NUMERIC_NON_CONTROL_SUBSET" if any(r["score"] is None for r in nc) else "ALL_NON_CONTROL_ASSUMPTION_SENSITIVE",
                **aggregate(nc), status="CANDIDATE_NOT_OFFICIAL_Q01_Q02",
                clinical_interpretation="genetic targeting only; no functional rescue or clinical risk calibration")
            cards.append(card)
            for system in SYSTEMS:
                subset = [r for r in nc if r["system"] == system]
                systems.append(dict(donor=donor, variant=variant, system=system, axes=len(subset), **aggregate(subset)))
    for name, data in [("axis_donor_masks.tsv", masks), ("axis_donor_contributions.tsv", contributions),
                       ("donor_scorecard.tsv", cards), ("donor_residual_by_system.tsv", systems),
                       ("missing_score_impact.tsv", missing)]:
        write_tsv(REV / "results" / name, data)
    # Attribute historical conflicts without treating a target residual as an input to new scores.
    oldcard = {short_donor(x["donor"]): x for x in read_tsv(canonical("donor_scorecard.tsv"))}
    contracts = {x["donor"]: x for x in read_tsv(canonical("donor_coverage_contract.tsv"))}
    matrices = {x["donor"]: x for x in read_tsv(canonical("donor_residual_matrix.tsv"))}
    replay = {short_donor(x["donor"]): x for x in replay}
    conflicts, syshist = [], []
    for card in cards:
        if card["variant"] != "legacy_stored":
            continue
        d = short_donor(card["donor"])
        h, matrix = contracts[d], matrices.get(d)
        old_T, old_R = number(h["corrected_total_isi_v2"]), number(h["residual_isi_v2"])
        matrix_R = math.fsum(number(matrix["masked_" + x]) for x in SYSTEMS) if matrix else None
        conflicts.append(dict(donor=d, stored_scorecard_coverage_pct=number(oldcard[d]["ISI_coverage%"]),
            stored_contract_coverage_pct=number(h["coverage_percent"]), stored_contract_T=old_T, stored_contract_R=old_R,
            algebra_C_from_contract_T_R_pct=(1-old_R/old_T)*100 if old_T else None,
            legacy_code_replay_coverage_pct=replay[d]["ISI_coverage%"],
            candidate_numeric_subset_T=card["T"], candidate_numeric_subset_A=card["A"], candidate_numeric_subset_R=card["R"],
            candidate_numeric_subset_C_pct=card["coverage_fraction"] * 100 if card["coverage_fraction"] is not None else None,
            candidate_minus_historical_R=card["R"]-old_R,
            saved_matrix_system_sum=matrix_R, historical_contract_R_minus_matrix_sum=old_R-matrix_R if matrix else None,
            saved_matrix_available=bool(matrix), unscored_non_control_axes=card["unscored_n"],
            historical_path="masked matrix sum + historical CFH 0.3 patch" if matrix else "historical worklog explicitly describes reverse inference from displayed coverage",
            attribution_status="HISTORY_TRACED_BUT_NOT_INDEPENDENT_AXIS_RECONSTRUCTION",
            remaining_gap="Saved system cells lack complete axis provenance; old coverage not reproduced; candidate scope is numeric subset only"))
        if matrix:
            for system in SYSTEMS:
                s = next(x for x in systems if x["donor"] == card["donor"] and x["variant"] == "legacy_stored" and x["system"] == system)
                syshist.append(dict(donor=d, system=system, legacy_matrix_unmasked=number(matrix[system]),
                    legacy_matrix_masked=number(matrix["masked_"+system]), candidate_T=s["T"], candidate_A=s["A"], candidate_R=s["R"],
                    legacy_masked_minus_candidate_R=number(matrix["masked_"+system])-s["R"],
                    status="NOT_SAME_COMPUTATION_NO_AXIS_LINEAGE_FOR_LEGACY_MATRIX"))
    write_tsv(REV / "results/coverage_conflict_attribution.tsv", conflicts)
    write_tsv(REV / "results/legacy_matrix_system_reconciliation.tsv", syshist)
    write_json(REV / "audit/scoring_invariants.json", dict(
        tolerance=dict(abs=1e-8, rel=1e-8), mask_rows=len(masks), contribution_rows=len(contributions),
        scorecards=len(cards), variants=VARIANTS,
        T_equals_A_plus_R=all(math.isclose(c["T"], c["A"]+c["R"], abs_tol=1e-8, rel_tol=1e-8) for c in cards),
        coverage_algebra=all(c["coverage_fraction"] is None if c["T"] == 0 else
                            math.isclose(c["coverage_fraction"], 1-c["R"]/c["T"], abs_tol=1e-8, rel_tol=1e-8) for c in cards),
        system_sum_closes=all(all(math.isclose(c[k], math.fsum(s[k] for s in systems if s["donor"] == c["donor"] and s["variant"] == c["variant"]), abs_tol=1e-8, rel_tol=1e-8) for k in ["T", "A", "R"]) for c in cards),
        legacy_formula_replay_matches_new_nominal50=all(math.isclose(scores(r)["formula_nominal50_sensitivity"], float(
            x["replay_formula_ISI_v2"]), abs_tol=1e-8, rel_tol=1e-8) for r, x in zip(axes, read_tsv(REV / "results/legacy_code_replay_axis_scores.tsv"))),
        missing_coerced_to_zero=False, final_release_eligible=False,
        unresolved=["Q01 scoring/severity/CFH exception", "Q02 direct targeting vs related axes and functional rescue",
                    "12 legacy non-control scores missing", "historical matrix has no complete axis-level provenance"] ))

if __name__ == "__main__":
    main()
