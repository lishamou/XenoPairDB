"""Explicit donor-genotype definitions and gene-to-pair mapping.

Compositions are taken verbatim from main-text Table 1 with literature
bases recorded in data/genotype_definitions.tsv. The mapping replaces the
hand-curated mapping of data release v1.0 (see README erratum).
"""

GENOTYPE_DEFINITIONS = {
    # genotype: (knockouts, knock-ins)
    "GTKO (NYU 2021-23)": (["GGTA1"], []),
    "ClonOrgan-6GE (Xijing auxiliary liver/kidney 2024-25)": (
        ["GGTA1", "CMAH", "B4GALNT2"], ["hCD46", "hCD55", "hTHBD"]),
    "ClonOrgan-8GE (Xijing kidney #2 2026)": (
        ["GGTA1", "CMAH", "B4GALNT2"],
        ["hCD46", "hCD55", "hTHBD", "hCD47", "hHMOX1"]),
    "Yunnan-10GE (Anhui auxiliary liver 2024)": (
        ["GGTA1", "CMAH", "B4GALNT2"],
        ["hCD46", "hCD55", "hCD59", "hTHBD", "hEPCR", "hENTPD1", "hCD47"]),
    "Revivicor-10GE (Maryland heart/NYU kidney)": (
        ["GGTA1", "CMAH", "B4GALNT2", "GHRHR"],
        ["hCD46", "hCD55", "hCD47", "hTHBD", "hEPCR", "hHMOX1"]),  # no hCD59
    "eGenesis-69GE (MGH kidney 2024-25)": (
        ["GGTA1", "CMAH", "B4GALNT2"],
        ["hCD46", "hCD55", "hTHBD", "hEPCR", "hCD47", "hCD200",
         "hTNFAIP3"]),  # + 59 PERV inactivations (no Tier-1 pair)
}

# edited gene -> Tier-1 pair_id (explicit alias map; None = no Tier-1 pair)
GENE_TO_PAIR = {
    "GGTA1": "XP004", "CMAH": "XP005", "B4GALNT2": None,
    "GHRHR": None,  # distinct from GHR; GH1/GHR pair XP053 not targeted
    "hCD46": "XP001", "hCD55": "XP002", "hCD59": "XP003",
    "hTHBD": "XP008", "hEPCR": "XP009", "hENTPD1": "XP013",
    "hCD47": "XP020", "hCD200": "XP021", "hHMOX1": None,
    "hTNFAIP3": None,
}

# corrected scorecard values (data release v1.1)
DELIVERED = {
    "GTKO (NYU 2021-23)": (1, 3.6),
    "ClonOrgan-6GE (Xijing auxiliary liver/kidney 2024-25)": (5, 11.8),
    "ClonOrgan-8GE (Xijing kidney #2 2026)": (6, 11.8),
    "Yunnan-10GE (Anhui auxiliary liver 2024)": (9, 15.3),
    "Revivicor-10GE (Maryland heart/NYU kidney)": (7, 11.8),
    "eGenesis-69GE (MGH kidney 2024-25)": (8, 13.3),
}


def fixed_pair_ids(genotype: str) -> set:
    """Tier-1 pair ids fixed by a genotype (KO of pig gene or KI of the
    human counterpart)."""
    kos, kis = GENOTYPE_DEFINITIONS[genotype]
    return {GENE_TO_PAIR[g] for g in kos + kis if GENE_TO_PAIR.get(g)}


def scorecard(t1):
    """Recompute the donor scorecard from a Tier-1 dataframe.

    Returns a dataframe with columns donor / pairs_fixed / ISI_coverage% /
    top_residual. ISI coverage = sum of ISI v2 over fixed incompatible
    pairs / sum of ISI v2 over all 56 incompatible pairs (pairs with
    undefined ISI_v2 contribute 0 to the numerator but count as fixed).
    """
    import pandas as pd

    t1 = t1.copy()
    t1["ISI_v2_num"] = pd.to_numeric(t1["ISI_v2"], errors="coerce")
    inc = t1[t1["evidence"] != "control"]
    total_isi = inc["ISI_v2_num"].sum()
    rows = []
    for gt in GENOTYPE_DEFINITIONS:
        ids = fixed_pair_ids(gt)
        fixed = t1[t1["pair_id"].isin(ids)]
        cov = fixed["ISI_v2_num"].fillna(0).sum() / total_isi * 100
        resid = inc[~inc["pair_id"].isin(ids)].sort_values(
            "ISI_v2_num", ascending=False)
        top = "; ".join(f"{r.pig_gene}/{r.human_partner}"
                        for r in resid.head(3).itertuples())
        rows.append(dict(donor=gt, pairs_fixed=len(fixed),
                         **{"ISI_coverage%": round(cov, 1)}, top_residual=top))
    return pd.DataFrame(rows)


def residual_isi_by_system(t1):
    """Per-system residual ISI v2 sums for each genotype (Fig. 2d heatmap).

    Uses ISI v2 (the published prioritized score, bounded [0, 125]) rather
    than raw ISI, consistent with the scorecard text. For non-1:1-ortholog
    pairs (e.g. GGTA1/CMAH) ISI_v2 is the literature-anchored literature
    evidence score carried by the curated Tier-1 row.
    """
    import pandas as pd

    t1 = t1.copy()
    t1["ISI_v2_num"] = pd.to_numeric(t1["ISI_v2"],
                                     errors="coerce").fillna(0)
    inc = t1[t1["evidence"] != "control"]
    systems = ["Complement", "Coagulation", "Immune", "Cytokine",
               "Endocrine", "Adhesion"]
    mat = {}
    for gt in GENOTYPE_DEFINITIONS:
        ids = fixed_pair_ids(gt)
        resid = inc[~inc["pair_id"].isin(ids)]
        by = resid.groupby("system")["ISI_v2_num"].sum()
        mat[gt.split(" (")[0]] = [int(round(by.get(s, 0))) for s in systems]
    return pd.DataFrame(mat, index=systems).T
