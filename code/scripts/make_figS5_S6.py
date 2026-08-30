#!/usr/bin/env python3
"""Produce Supplementary Figs. S5--S6 and deposit their source-data files.

The panels deliberately distinguish donor--recipient transcript co-observation
from functional interaction and a one-recipient annotation-conditional stress
test from clinical validation.  All values are read from frozen analysis outputs.
"""
from __future__ import annotations

import json
import os
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PIPELINE = Path(__file__).resolve().parents[2]
PROJECT = Path(os.environ.get("XENOPAIR_PROJECT_ROOT", str(PIPELINE.parent)))
OUT = PIPELINE / "figures"
RESULTS = PIPELINE / "results_r172_v37"
CLINICAL = RESULTS / "graft_observability"
BR = RESULTS / "br0010_bulk_proteomics"
ZENODO = PIPELINE / "release_data" / "clinical"
TIER1 = PIPELINE / "canonical_links" / "tier1_frozen.tsv"
TPM = RESULTS / "input_audit" / "BR0010_xeno_whole_blood_bulk_tpm.csv"

plt.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.7, "axes.spines.top": False, "axes.spines.right": False,
})
COL = {"addressed": "#3C78A8", "unresolved": "#C46A3A", "not_evaluable": "#9A9A9A",
       "dark": "#263C54", "neutral": "#E7E4DE", "signal": "#587C9A"}


def save(fig: plt.Figure, stem: str, conclusion: str, inputs: list[str]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=600 if ext == "png" else 300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    provenance = {
        "figure": f"{stem}.png", "produced_by": "src/_scripts/make_figS5_S6.py",
        "command": "python src/_scripts/make_figS5_S6.py",
        "one_line": conclusion, "core_conclusion": conclusion, "inputs": inputs,
        "env": {"python": platform.python_version(), "matplotlib": matplotlib.__version__,
                "numpy": np.__version__, "pandas": pd.__version__}, "backend": "python",
        "seed": 172, "git_commit": None,
        "conversation_ref": "worklog/2026-08-20-R176-v39-molecular-context-manuscript-sync.md",
        "created": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / f"{stem}.prov.json").write_text(json.dumps(provenance, indent=2) + "\n")


def label(ax, letter: str) -> None:
    ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontweight="bold", fontsize=10,
            va="bottom", ha="left")


def deposit_sources() -> None:
    files = [
        CLINICAL / "graft_axis_observability.tsv", CLINICAL / "pig_gene_celltype_time.tsv",
        CLINICAL / "human_gene_celltype_time.tsv", CLINICAL / "observability_threshold_sensitivity.tsv",
        CLINICAL / "clinical_pathway_scoring_sensitivity_timeseries.tsv",
        CLINICAL / "clinical_pathway_scoring_sensitivity_summary.tsv",
        CLINICAL / "clinical_pathway_scoring_sensitivity_metadata.json",
        CLINICAL / "partner_parsing_funnel.tsv", BR / "br0010_axis_measurement_map.tsv",
        BR / "br0010_proteome_measurement_map.tsv", BR / "br0010_residual_vs_addressed_conditional_randomization.tsv",
        BR / "br0010_system_conditioned_randomization.tsv", BR / "br0010_system_conditioned_randomization_null_draws.tsv",
        BR / "br0010_system_score_timeseries.tsv", BR / "br0010_scoring_sensitivity_timeseries.tsv", BR / "br0010_transcript_protein_exact_gene_map.tsv",
        BR / "br0010_sample_contract.tsv",
    ]
    ZENODO.mkdir(parents=True, exist_ok=True)
    for source in files:
        shutil.copy2(source, ZENODO / source.name)


def fig_s5() -> None:
    funnel = pd.read_csv(CLINICAL / "partner_parsing_funnel.tsv", sep="\t")
    threshold = pd.read_csv(CLINICAL / "observability_threshold_sensitivity.tsv", sep="\t")
    obs = pd.read_csv(CLINICAL / "graft_axis_observability.tsv", sep="\t")
    preferred_axes = ["XP002", "XP038", "XP050", "XP067"]
    axes = [x for x in preferred_axes if x in set(obs.pair_id)]
    times = ["D0", "D10", "D14", "D33", "D49", "D61"]
    names = {"XP002": "XP002  CD55–C3/C5", "XP038": "XP038  IL1R1–IL1B",
             "XP050": "XP050  REN–AGT", "XP067": "XP067  C5AR1–C5"}
    fig = plt.figure(figsize=(15.0, 5.6))
    ax1 = fig.add_axes([0.05, 0.18, 0.24, 0.68])
    vals = funnel.n_axes.tolist()
    # The parsing/evaluability funnel stops at 17 before detection filtering.
    # Keep the nominal four-axis endpoint visible, then show the three-axis
    # minimum-cell-qualified endpoint used in panel c and the manuscript.
    vals.append(obs.pair_id.nunique())
    stages = ["Tier-1\naxes", "binary HGNC\nresolvable", "pig gene\npresent", "human partner\npresent", "both sides\nobservable", "cell-count-\nqualified"]
    y = np.arange(len(vals))[::-1]
    ax1.barh(y, vals, color=["#4E6A83", "#6687A3", "#83A0B8", "#A3B7C5", "#D08B5D", "#B56F46"], height=.68)
    for yy, value in zip(y, vals): ax1.text(value + 1.0, yy, str(value), va="center", fontweight="bold")
    ax1.set_yticks(y, stages); ax1.set_xlabel("axes", fontsize=7); ax1.set_xlim(0, 82); label(ax1, "a")
    ax1.tick_params(axis="both", labelsize=7)
    ax1.set_title("Strict evaluability funnel", loc="left", fontsize=9, pad=8)

    ax2 = fig.add_axes([0.36, 0.18, 0.17, 0.68])
    threshold["percent"] = threshold.detection_threshold.mul(100).astype(int).astype(str) + "%"
    colors = [COL["neutral"], "#587C9A", "#A7B7C3"]
    ax2.bar(threshold.percent, threshold.n_distinct_axes, color=colors, width=.62)
    for i, value in enumerate(threshold.n_distinct_axes): ax2.text(i, value + .12, str(value), ha="center", fontweight="bold")
    ax2.axvspan(.5, 1.5, color="#DCE8F1", zorder=-1)
    ax2.text(1, .35, "primary", ha="center", fontsize=6, color=COL["dark"])
    ax2.set_ylim(0, 4.8); ax2.set_ylabel("observable axes", fontsize=7); label(ax2, "b")
    ax2.tick_params(axis="both", labelsize=7)
    ax2.set_title("Detection-threshold sensitivity", loc="left", fontsize=9, pad=8)

    ax3 = fig.add_axes([0.61, 0.18, 0.35, 0.68])
    grid = np.zeros((len(axes), len(times)))
    labels = [["—" for _ in times] for _ in axes]
    for i, pair in enumerate(axes):
        for j, time in enumerate(times):
            sub = obs[(obs.pair_id == pair) & (obs.time_point == time)].copy()
            if sub.empty: continue
            grid[i, j] = len(sub)
            top = sub.assign(score=sub.pig_detection_fraction * sub.human_detection_fraction).sort_values("score", ascending=False).iloc[0]
            labels[i][j] = f"{top.pig_cell_type} →\n{top.human_cell_type}"
    cmap = LinearSegmentedColormap.from_list("obs", ["#F5F4F0", "#9BB5C8", "#3E6F91"])
    im = ax3.imshow(grid, aspect="auto", cmap=cmap, vmin=0, vmax=max(1, grid.max()))
    ax3.set_xticks(range(len(times)), times); ax3.set_yticks(range(len(axes)), [names[x] for x in axes])
    ax3.tick_params(axis="x", length=0); ax3.tick_params(axis="y", length=0)
    for i in range(len(axes)):
        for j in range(len(times)):
            ax3.text(j, i, labels[i][j], ha="center", va="center", fontsize=6.2,
                     color="white" if grid[i, j] > grid.max() * .55 else "#263C54")
    cb = fig.colorbar(im, ax=ax3, fraction=.035, pad=.02); cb.set_label("passing donor–recipient\ncell-type pairs", fontsize=7)
    ax3.tick_params(axis="both", labelsize=7)
    label(ax3, "c"); ax3.set_title("Representative passing cell-type pair at each retained time point", loc="left", fontsize=9, pad=8)
    save(fig, "figS5_graft_observability",
         "Three Tier-1 axes meet a cell-count-qualified donor–recipient transcript co-observation criterion in the GTKO graft.",
         [str(x.relative_to(PIPELINE)) for x in [CLINICAL / "graft_axis_observability.tsv", CLINICAL / "observability_threshold_sensitivity.tsv", CLINICAL / "partner_parsing_funnel.tsv"]])


def _null_distributions() -> tuple[np.ndarray, np.ndarray]:
    raw = pd.read_csv(TPM, header=None)
    genes = raw.iloc[0, 1:-3].astype(str).tolist()
    samples = raw.iloc[1:, 0].astype(str).tolist()
    tpm = pd.DataFrame(raw.iloc[1:, 1:-3].to_numpy(dtype=float), index=samples, columns=genes)
    log = np.log2(tpm + .5); baseline = log.loc[["RGc21_S1", "RGc22_S2"]].mean(axis=0)
    tier = pd.read_csv(TIER1, sep="\t")
    mapping = pd.read_csv(BR / "br0010_axis_measurement_map.tsv", sep="\t")
    pairs = mapping[mapping.analysis_status.isin(["genetically_addressed", "genetically_unresolved"])]
    group = {k: sorted({gene for genestring in x.measured_human_genes.fillna("") for gene in genestring.split(";") if gene}) for k, x in pairs.groupby("analysis_status")}
    event = (log.loc["RGc24_S4"] - baseline).abs()
    background = pd.DataFrame({"gene": genes, "value": event.values, "baseline": baseline.values, "detect": (tpm > 0).sum(axis=0).values})
    background["bin"] = pd.qcut(background.baseline.rank(method="first"), 5, labels=False)
    lookup = background.set_index("gene")
    by_bin = background.groupby(["bin", "detect"])
    pools = {}
    for gene in set(group["genetically_addressed"] + group["genetically_unresolved"]):
        row = lookup.loc[gene]; key = (row.bin, row.detect)
        pools[gene] = by_bin.get_group(key).gene.to_numpy() if key in by_bin.groups else background[background.detect.eq(row.detect)].gene.to_numpy()
    rng = np.random.default_rng(172)
    def draw(metric):
        def sample(gset):
            result = np.zeros(10000)
            for gene in gset: result += rng.choice(metric.loc[pools[gene]].to_numpy(), size=10000, replace=True)
            return result / len(gset)
        return sample(group["genetically_unresolved"]) - sample(group["genetically_addressed"])
    return draw(event), draw(event.rank(method="average", pct=True))


def fig_s6() -> None:
    mapping = pd.read_csv(BR / "br0010_axis_measurement_map.tsv", sep="\t")
    randomization = pd.read_csv(BR / "br0010_residual_vs_addressed_conditional_randomization.tsv", sep="\t")
    system_conditioned = pd.read_csv(BR / "br0010_system_conditioned_randomization.tsv", sep="\t")
    system_conditioned_null = pd.read_csv(BR / "br0010_system_conditioned_randomization_null_draws.tsv", sep="\t")
    protein = pd.read_csv(BR / "br0010_transcript_protein_exact_gene_map.tsv", sep="\t")
    null_mag, null_rank = _null_distributions()
    fig = plt.figure(figsize=(16.0, 6.2))
    ax1 = fig.add_axes([.05, .56, .22, .32]); ax2 = fig.add_axes([.35, .56, .18, .32]); ax3 = fig.add_axes([.62, .56, .32, .32])
    ax4 = fig.add_axes([.05, .10, .42, .30]); ax5 = fig.add_axes([.61, .10, .33, .30])
    noncontrol = mapping[mapping.analysis_status != "control_excluded"]
    counts = noncontrol.analysis_status.value_counts().reindex(["genetically_addressed", "genetically_unresolved", "not_modality_evaluable"])
    ax1.bar(range(3), counts, color=[COL["addressed"], COL["unresolved"], COL["not_evaluable"]])
    ax1.set_xticks(range(3), ["addressed", "unresolved", "not\nevaluable"], fontsize=7); ax1.set_ylabel("non-control axes", fontsize=7)
    for i, v in enumerate(counts): ax1.text(i, v + 1.3, str(v), ha="center", fontweight="bold")
    ax1.set_ylim(0, 58); ax1.tick_params(axis="y", labelsize=7); label(ax1, "a"); ax1.set_title("Prespecified genotype mapping", loc="left", fontsize=9, pad=8)
    usable = mapping[(mapping.analysis_status.isin(["genetically_addressed", "genetically_unresolved"])) & (mapping.n_measured_human_genes > 0)]
    denominators = [7, 50]; values = [usable.analysis_status.eq("genetically_addressed").sum(), usable.analysis_status.eq("genetically_unresolved").sum()]
    ax2.bar(range(2), values, color=[COL["addressed"], COL["unresolved"]]); ax2.set_xticks(range(2), ["addressed", "unresolved"], fontsize=7)
    ax2.set_ylabel("RNA-evaluable axes", fontsize=7); ax2.set_ylim(0, 38); ax2.tick_params(axis="y", labelsize=7)
    for i, (v, d) in enumerate(zip(values, denominators)): ax2.text(i, v + 1, f"{v}/{d}", ha="center", fontweight="bold")
    label(ax2, "b"); ax2.set_title("Exact RNA-feature coverage", loc="left", fontsize=9, pad=8)
    # Panel c must use the primary gene-wise-z output, not the separate
    # baseline-referenced log2(TPM) system-score file.  The latter is retained
    # for the magnitude diagnostic in panel d and has a different estimand.
    system = pd.read_csv(BR / "br0010_scoring_sensitivity_timeseries.tsv", sep="\t")
    system = system[(system.analysis_status == "genetically_unresolved") & (system.timepoint == "D7") & (system.method == "gene_wise_z_mean_primary")]
    system = system.set_index("system").reindex(["Adhesion", "Coagulation", "Complement", "Cytokine", "Immune"])
    ax3.axhline(0, color="#777", lw=.7); ax3.bar(range(len(system)), system.score, color=[COL["unresolved"] if x >= 0 else "#7895AA" for x in system.score])
    ax3.set_xticks(range(len(system)), [x.replace("Coagulation", "Coag.") for x in system.index], rotation=20, ha="right", fontsize=7)
    ax3.set_ylabel("D7 primary gene-wise z mean", fontsize=7); ax3.tick_params(axis="y", labelsize=7); label(ax3, "c"); ax3.set_title("Heterogeneous unresolved-system response", loc="left", fontsize=9, pad=8)
    ax4a = ax4.inset_axes([0.00, .0, .47, 1]); ax4b = ax4.inset_axes([.53, .0, .47, 1]); ax4.set_axis_off()
    for sub, null, metric, title in [(ax4a, null_mag, "mean_abs_log2fc_D7_vs_mean_pre_tx", "Magnitude diagnostic"), (ax4b, null_rank, "mean_event_rank_percentile_D7", "Rank diagnostic")]:
        observed = randomization[(randomization.metric == metric) & (randomization.group == "unresolved_minus_addressed")].observed_value.iloc[0]
        primary_p = randomization[(randomization.metric == metric) & (randomization.group == "unresolved_minus_addressed")].conditional_randomization_p_one_sided.iloc[0]
        cond = system_conditioned[(system_conditioned.metric == metric) & (system_conditioned.group == "unresolved_minus_addressed")].iloc[0]
        cond_null = system_conditioned_null[system_conditioned_null.metric == metric].null_difference.to_numpy()
        sub.hist(null, bins=35, color="#C8D6E0", edgecolor="white", linewidth=.2, alpha=.75, label="abundance + detectability")
        sub.hist(cond_null, bins=35, color="#6F8FA8", edgecolor="white", linewidth=.2, alpha=.55, label="+ system composition")
        sub.axvline(observed, color=COL["unresolved"], lw=1.5)
        n_addr = int(randomization[(randomization.metric == metric) & (randomization.group == "genetically_addressed")].n_genes.iloc[0])
        n_unres = int(randomization[(randomization.metric == metric) & (randomization.group == "genetically_unresolved")].n_genes.iloc[0])
        sub.text(.97, .92, f"n = {n_addr}/{n_unres} genes\nΔ = {observed:.3f}\nPtail {primary_p:.3f} → {cond.system_conditioned_randomization_p_one_sided:.3f}", transform=sub.transAxes, ha="right", va="top", fontsize=7)
        sub.set_title(title, loc="left", fontsize=8); sub.set_xlabel("matched null difference", fontsize=7); sub.set_ylabel("draws", fontsize=7); sub.tick_params(labelsize=7)
        if metric == "mean_abs_log2fc_D7_vs_mean_pre_tx": sub.legend(frameon=False, fontsize=6.2, loc="upper left")
    label(ax4a, "d")
    protein = protein.dropna(subset=["rna_log2fc_D7_vs_D26", "protein_log2fc_D7_vs_D26"])
    colors = protein.analysis_status.map({"genetically_addressed": COL["addressed"], "genetically_unresolved": COL["unresolved"], "control_excluded": "#A3A3A3"}).fillna("#A3A3A3")
    ax5.scatter(protein.rna_log2fc_D7_vs_D26, protein.protein_log2fc_D7_vs_D26, c=colors, s=14, alpha=.75, linewidths=0)
    lim = max(abs(protein.rna_log2fc_D7_vs_D26).max(), abs(protein.protein_log2fc_D7_vs_D26).max()) * 1.05
    ax5.plot([-lim, lim], [-lim, lim], color="#777", lw=.7, ls="--"); ax5.set_xlim(-lim, lim); ax5.set_ylim(-lim, lim)
    ax5.set_xlabel("RNA log2 fold change (D7 vs D26)", fontsize=7); ax5.set_ylabel("protein log2 fold change (D7 vs D26)", fontsize=7)
    ax5.tick_params(labelsize=7)
    ax5.text(.97, .93, "n = 54\nSpearman ρ = 0.129\n55.6% sign agreement", transform=ax5.transAxes, ha="right", va="top", fontsize=7)
    label(ax5, "e"); ax5.set_title("Exact RNA–protein matches", loc="left", fontsize=9, pad=8)
    save(fig, "figS6_BR0010_stress_test",
         "The BR0010 single-recipient independent descriptive stress test is directional but system-composition-sensitive, with limited RNA–protein concordance.",
         [str(x.relative_to(PIPELINE)) for x in [BR / "br0010_axis_measurement_map.tsv", BR / "br0010_residual_vs_addressed_conditional_randomization.tsv", BR / "br0010_system_conditioned_randomization.tsv", BR / "br0010_system_conditioned_randomization_null_draws.tsv", BR / "br0010_scoring_sensitivity_timeseries.tsv", BR / "br0010_system_score_timeseries.tsv", BR / "br0010_transcript_protein_exact_gene_map.tsv"]])


if __name__ == "__main__":
    deposit_sources()
    fig_s5()
    fig_s6()
    print("wrote Supplementary Figs. S5 and S6 plus deposited source-data files")
