#!/usr/bin/env python3
"""R23 Figure regeneration:
- Figure 1: extend from 3-panel to 5-panel (A/B/C + D + E)
- Figure 3: add pipeline schematic panel (a)
- Figure 5: re-layout as 2x2 (a/b on row 1, c/d on row 2)
- Figure 6: improved color palette
"""
import os
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
import numpy as np

# R110-A: pipeline self-contained — derive paths from script location
_PIPELINE_ROOT = Path(__file__).resolve().parents[2]   # pipeline/
_PROJECT_ROOT  = Path(__file__).resolve().parents[3]   # project root (28XenoPairDB/)
ROOT = Path(os.environ.get("XENOPAIR_PROJECT_ROOT", str(_PROJECT_ROOT)))
DATA = _PIPELINE_ROOT / "release_data"
OUT_NBT = _PIPELINE_ROOT / "figures"   # R110-A: write to pipeline/figures/

# Consistent color palette (NBT Resource style)
C = dict(
    direct="#c0392b", disputed="#e67e22", inferred="#b79580", control="#7d8c5a",
    immune="#c0392b", coagulation="#d35400", cytokine="#e67e22",
    endocrine="#2980b9", adhesion="#7d8c5a", complement="#16a085",
    accent="#34495e", panel="#5d6d7e",
    bg_light="#fafafa", bg_box="#f8f4eb", bg_box2="#fbf3df",
    bg_box3="#eef0e6", bg_box4="#fde9d9", bg_grey="#f2f2f2",
    scatter_alpha=0.75,
    # Fig 6 improved palette
    f6_comp1="#7b1f1f",  # dark red - complement cascade
    f6_comp2="#a93226",  # red - complement regulation
    f6_comp3="#cb6e6e",  # light red - complement receptors
    f6_coag1="#1a3a5c",  # dark navy - coagulation
    f6_coag2="#2e5d8a",  # navy - anticoag_fibrinolysis
    f6_adh="#5d6d4f",    # olive - adhesion LFA1
    f6_inn="#8b6f47",    # brown - innate myeloid
    f6_ada="#a05a8e",    # mauve - adaptive infiltrate
    f6_cd47="#7a8b7a",   # sage - checkpoint CD47
)


def save(fig, name, formats=("png", "pdf", "svg")):
    for ext in formats:
        out = OUT_NBT / f"{name}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {name}: {', '.join(formats)}")


# ===================== FIGURE 1: 5-panel =====================
def fig1():
    """Figure 1: Resource architecture (5 panels).
    Panel A: Tiered architecture + validation loop
    Panel B: Tier-1 69 pairs by system & evidence
    Panel C: ISI v2 violin by system
    Panel D: Identity vs ISI scatter (per pair)
    Panel E: Predicted residual vs observed clinical outcomes (3 cases)
    """
    # Load Tier-1
    with open(_PIPELINE_ROOT / "canonical_links" / "tier1_frozen.tsv") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        t1 = list(rdr)

    systems = ["Immune", "Coagulation", "Cytokine", "Endocrine", "Adhesion",
               "Complement"]
    system_colors = {"Immune": C["immune"], "Coagulation": C["coagulation"],
                     "Cytokine": C["cytokine"], "Endocrine": C["endocrine"],
                     "Adhesion": C["adhesion"], "Complement": C["complement"]}

    fig = plt.figure(figsize=(15.0, 10.0))

    # ===== Panel A (top-left, large): tiered architecture =====
    axA = fig.add_axes([0.03, 0.55, 0.45, 0.40])
    axA.axis("off")
    # Keep a small drawing margin so the right-hand frame is never clipped by
    # the axes boundary in PDF/SVG export.
    axA.set_xlim(-0.02, 1.02); axA.set_ylim(0, 1)
    axA.text(0.0, 0.95, "a.", fontsize=12, fontweight="bold")

    boxes = [
        (0.73, C["bg_box"], C["direct"],
         "Tier 1 · Evidence-curated strict core",
         "75 pig–human directional axes",
         "direct 58 · disputed 4 · inferred 0 · control 13"),
        (0.445, C["bg_box2"], C["disputed"],
         "Tier 2 · Computational expansion",
         "1,407 computational candidate pairs",
         "sequence divergence × N-glycosylation sequon shift"),
        (0.16, C["bg_box3"], C["control"],
         "Tier 3 · Surfaceome antigenicity screen",
         "4,123-protein membrane universe (4,122 scoreable)",
         ""),
    ]
    for y, fc, ec, l1, l2, l3 in boxes:
        axA.add_patch(Rectangle((0.02, y), 0.72, 0.22, fc=fc, ec=ec, lw=1.4))
        axA.text(0.38, y + 0.16, l1, ha="center", fontsize=10,
                 fontweight="bold")
        axA.text(0.38, y + 0.10, l2, ha="center", fontsize=9.5)
        axA.text(0.38, y + 0.04, l3, ha="center", fontsize=8)
    for y0, y1 in ((0.73, 0.665), (0.445, 0.375)):
        axA.annotate("", xy=(0.38, y1), xytext=(0.38, y0),
                     arrowprops=dict(arrowstyle="-|>", color="#666", lw=1.2))

    # Validation loop on right
    axA.add_patch(Rectangle((0.79, 0.43), 0.18, 0.32, fc=C["bg_grey"],
                            ec=C["accent"], lw=1.4))
    axA.text(0.88, 0.70, "Robustness &", ha="center", fontsize=10,
             fontweight="bold")
    axA.text(0.88, 0.65, "evaluation loop", ha="center", fontsize=10,
             fontweight="bold")
    axA.text(0.88, 0.50,
             "500× 80/20 +\nLOSO + antigen\nbenchmark +\n61-day graft\nconcordance",
             ha="center", fontsize=7.5)
    axA.annotate("", xy=(0.79, 0.59), xytext=(0.74, 0.47),
                 arrowprops=dict(arrowstyle="-|>", color="#666", lw=1.2))
    axA.annotate("", xy=(0.74, 0.18), xytext=(0.86, 0.43),
                 arrowprops=dict(arrowstyle="-|>", color=C["control"],
                                 linestyle="--", lw=1.2))
    axA.text(0.80, 0.30, "feedback", fontsize=8, style="italic",
             color=C["control"])

    # Outputs (bottom)
    for x, txt in ((0.05, "ISI v2 prioritization"),
                   (0.32, "donor-pig engineering\nscorecard (6 genotypes)"),
                   (0.60, "XenoPairDB portal\nv1.1 · JSON API")):
        axA.add_patch(Rectangle((x, 0.00), 0.22, 0.10, fc=C["bg_box4"],
                                ec="#999", lw=1.0))
        axA.text(x + 0.11, 0.05, txt, ha="center", fontsize=8.5)
    for x0, x1 in ((0.27, 0.32), (0.54, 0.60)):
        axA.annotate("", xy=(x1, 0.04), xytext=(x0, 0.04),
                     arrowprops=dict(arrowstyle="-|>", color="#666", lw=1.0))

    # ===== Panel B (top-right): 69 pairs stacked bar =====
    axB = fig.add_axes([0.55, 0.55, 0.42, 0.40])
    grades = ["direct", "disputed", "inferred", "control"]
    counts = {g: [sum(1 for x in t1 if x["system"] == s and x["evidence"] == g)
                  for s in systems] for g in grades}
    x = np.arange(len(systems))
    bottom = np.zeros(len(systems))
    for g in grades:
        axB.bar(x, counts[g], bottom=bottom, color=C[g], label=g,
                width=0.65, edgecolor="white", linewidth=0.8)
        bottom += np.array(counts[g])
    for i, s in enumerate(systems):
        axB.text(i, bottom[i] + 0.4, str(int(bottom[i])), ha="center",
                 fontweight="bold", fontsize=10)
    axB.set_xticks(x)
    axB.set_xticklabels(systems, rotation=25, ha="right", fontsize=9)
    axB.set_ylabel("pairs")
    axB.set_title("b.", loc="left", fontsize=11, fontweight="bold")
    axB.legend(frameon=False, ncol=4, loc="upper center",
               bbox_to_anchor=(0.5, -0.30), fontsize=9)
    axB.spines[["top", "right"]].set_visible(False)
    axB.set_ylim(0, max(bottom) * 1.15)

    # ===== Panel C (middle-left): ISI v2 violin by system =====
    axC = fig.add_axes([0.03, 0.06, 0.32, 0.42])
    data = []
    for s in systems:
        vals = [float(x["ISI_v2"]) for x in t1
                if x["system"] == s and x["ISI_v2"] not in ("n/a", "")]
        data.append(vals)
    parts = axC.violinplot(data, showmedians=True)
    for pc, s in zip(parts["bodies"], systems):
        pc.set_facecolor(system_colors[s])
        pc.set_alpha(0.85)
        pc.set_edgecolor(C["accent"])
    for i, d in enumerate(data):
        axC.scatter(np.full(len(d), i + 1), d, s=15, color="#333", zorder=3,
                    alpha=0.7)
    axC.set_xticks(range(1, 7))
    axC.set_xticklabels(systems, rotation=25, ha="right", fontsize=9)
    axC.set_ylabel("Incompatibility Severity Index (ISI v2)")
    axC.set_title("c.", loc="left", fontsize=11, fontweight="bold")
    axC.spines[["top", "right"]].set_visible(False)
    axC.axhline(0, color="#666", lw=0.5, ls="--", alpha=0.5)

    # ===== Panel D (middle-center): identity vs ISI scatter =====
    axD = fig.add_axes([0.40, 0.06, 0.32, 0.42])
    for s in systems:
        xs = []  # identity
        ys = []  # ISI v2
        for x in t1:
            if x["system"] == s and x["identity"] and x["ISI_v2"]:
                try:
                    xs.append(float(x["identity"]))
                    ys.append(float(x["ISI_v2"]))
                except ValueError:
                    pass
        axD.scatter(xs, ys, s=50, color=system_colors[s], alpha=0.75,
                    edgecolor="white", linewidth=0.8, label=s, zorder=3)
    # annotate landmarks
    landmarks = {"CD58": (40.2, 49.1), "IFNG": (43.0, 44.0),
                 "THBD": (62.0, 33.0), "REN": (76.0, 28.0),
                 "VWF": (84.5, 19.0), "EPO": (83.3, 12.0)}
    for name, (xi, yi) in landmarks.items():
        axD.annotate(name, (xi, yi), xytext=(5, 5), textcoords="offset points",
                     fontsize=8.5, fontweight="bold", color=C["accent"])
    axD.set_xlabel("human–pig global sequence identity (%)")
    axD.set_ylabel("ISI v2")
    axD.set_title("d.", loc="left", fontsize=11, fontweight="bold")
    axD.legend(frameon=False, fontsize=8, loc="upper right",
               ncol=2, scatterpoints=1)
    axD.spines[["top", "right"]].set_visible(False)
    axD.grid(axis="y", linestyle=":", alpha=0.4)

    # ===== Panel E (middle-right): clinical anchoring stacked bars =====
    # Read the shared canonical donor residual matrix used by Fig. 2d and
    # Fig. S1d; the masked columns are the submission metric.
    axE = fig.add_axes([0.74, 0.06, 0.23, 0.42])
    import pandas as pd
    matrix_path = _PIPELINE_ROOT / "canonical_links" / "donor_residual_matrix.tsv"
    donor_matrix = pd.read_csv(matrix_path, sep="\t").set_index("donor")
    mat = donor_matrix[[f"masked_{s}" for s in systems]].copy()
    mat.columns = systems
    # Get 3 representative genotypes
    cases = [
        ("Yunnan-10GE\nliver context\nobserved: xTMA, day 38",
         "Yunnan-10GE"),
        ("Revivicor-10GE\nheart context\nobserved: AMR × 2",
         "Revivicor-10GE"),
        ("eGenesis-69GE\nkidney context\nobserved: chronic decline",
         "eGenesis-69GE"),
    ]
    # Get the row for each genotype (short name)
    case_colors = system_colors
    totals = []
    for i, (label, gt_long) in enumerate(cases):
        # match short name
        short = gt_long.split(" (")[0]
        if short not in mat.index:
            short = next(idx for idx in mat.index if idx.startswith(short))
        row = mat.loc[short]
        # Show the complete six-system residual sum used by the canonical
        # donor matrix; restricting this panel to three systems would create a
        # misleading subtotal labelled as a donor total.
        vals = {s: int(row[s]) for s in systems}
        total = sum(vals.values())
        totals.append(total)
        bottom = 0
        for cat, v in vals.items():
            axE.bar(i, v, bottom=bottom, color=case_colors[cat],
                    edgecolor="white", linewidth=1, width=0.7)
            if v > 25:
                axE.text(i, bottom + v / 2, str(v), ha="center", va="center",
                         color="white", fontsize=10, fontweight="bold")
            bottom += v
        axE.text(i, bottom + 8, f"total {total}", ha="center",
                 fontweight="bold", fontsize=9)
    axE.set_xticks(range(3))
    axE.set_xticklabels(["Yunnan\n10GE\nxTMA", "Revivicor\n10GE\nAMR", "eGenesis\n69GE\nrenal decline"], fontsize=8)
    axE.set_ylabel("residual ISI v2 sum")
    axE.set_title("e.", loc="left", fontsize=11, fontweight="bold")
    axE.set_ylim(0, max(500, max(totals) * 1.16))
    legend_elems = [Rectangle((0, 0), 1, 1, color=system_colors[s], label=s)
                    for s in systems]
    axE.legend(handles=legend_elems, frameon=False, fontsize=7,
               loc="upper left", ncol=2)
    axE.spines[["top", "right"]].set_visible(False)

    save(fig, "fig1_resource")
def fig3():
    """Figure 3: Structural rationale for cross-species interface.
    Panel (a): template-based complex-modeling pipeline schematic
    Panel (b): bar chart of identity/contacts per case + ICAM1 anchor list
    """
    fig = plt.figure(figsize=(12.5, 7.5))

    # ===== Panel a: pipeline schematic =====
    axA = fig.add_axes([0.05, 0.54, 0.90, 0.34])
    axA.axis("off")
    axA.set_xlim(0, 1); axA.set_ylim(0, 1)
    axA.text(0.0, 0.97, "a.", fontsize=11, fontweight="bold")

    # Step boxes (all kept within axA xlim [0, 1])
    box_w = 0.145
    # Deliberate white gutters keep arrows between boxes rather than inside
    # their borders; the last box remains comfortably inside the frame.
    box_xs = [0.03, 0.225, 0.42, 0.615, 0.81]
    steps = [
        ("Porcine\nAlphaFold\nmonomer", C["bg_box3"], C["control"]),
        ("Sequence\nalignment\n(BLOSUM62)\n→ Cα mapping", C["bg_box2"],
         C["disputed"]),
        ("Rigid\nsuperposition\non human\ncomplex template\n(PDB)",
         C["bg_box"], C["direct"]),
        ("Iterative\noutlier removal\n(2.5× median)", C["bg_box4"],
         C["accent"]),
        ("Contact count\n(<8 Å) +\nclash count\n(<3 Å)", "#f0e6d6",
         "#8b4513"),
    ]
    for x, (label, fc, ec) in zip(box_xs, steps):
        axA.add_patch(FancyBboxPatch((x, 0.46), box_w, 0.22,
                                     boxstyle="round,pad=0.02",
                                     fc=fc, ec=ec, lw=1.4))
        axA.text(x + box_w / 2, 0.57, label, ha="center", va="center",
                 fontsize=7.5)
    # arrows between boxes (4 arrows for 5 boxes)
    arrow_xs = [box_xs[i] + box_w for i in range(len(steps) - 1)]
    for x in arrow_xs:
        axA.annotate("", xy=(x + 0.015, 0.57), xytext=(x, 0.57),
                     arrowprops=dict(arrowstyle="-|>", color="#666", lw=1.4))

    # Sub-caption
    axA.text(0.5, 0.24,
             "Rigid-body reconstruction from AlphaFold/PDB templates",
             ha="center", fontsize=8.5, style="italic", color="#555")

    # ===== Panel b: bar chart =====
    axB = fig.add_axes([0.08, 0.10, 0.50, 0.34])
    cases = ["ICAM1–LFA-1", "ULBP1–NKG2D", "CD86–CTLA4"]
    identity = [25.0, 38.5, 50.0]
    contacts = [23, 27, 20]
    x = np.arange(len(cases))
    width = 0.38
    bars1 = axB.bar(x - width/2, identity, width, color="#a8c69f",
                    edgecolor=C["control"], lw=1.2,
                    label="Interface identity (%)")
    bars2 = axB.bar(x + width/2, contacts, width, color=C["control"],
                    edgecolor="#3d4a2e", lw=1.2,
                    label="Modeled contacts (<8 Å)")
    for bar, v in zip(bars1, identity):
        axB.text(bar.get_x() + bar.get_width()/2, v + 1.2, f"{v}",
                 ha="center", fontsize=9.5, fontweight="bold")
    for bar, v in zip(bars2, contacts):
        axB.text(bar.get_x() + bar.get_width()/2, v + 1.2, f"{v}",
                 ha="center", fontsize=9.5, fontweight="bold")
    axB.set_xticks(x)
    axB.set_xticklabels(cases, fontsize=10)
    axB.set_ylim(0, 60)
    axB.set_ylabel("Interface identity / contact count")
    axB.set_title("b.", loc="left", fontsize=11, fontweight="bold")
    axB.legend(frameon=False, fontsize=8, loc="upper left",
               bbox_to_anchor=(0.0, 1.0))
    axB.spines[["top", "right"]].set_visible(False)
    axB.grid(axis="y", linestyle=":", alpha=0.4)

    # ===== Panel c: ICAM1 anchors (visual contact table) =====
    # Match panel b's vertical footprint; the table rows fill the same data
    # height instead of floating in the upper half of the panel.
    axC = fig.add_axes([0.64, 0.10, 0.31, 0.34])
    axC.axis("off")
    axC.set_xlim(0, 1); axC.set_ylim(0, 1)
    axC.text(0.0, 0.97, "c.", fontsize=10.5,
             fontweight="bold", color="#222")
    axC.text(0.0, 0.91, "contact residues (<5 Å; 1MQ8)",
             fontsize=8.5, color="#666")
    anchors = [
        ("Glu34", "E", "identical", True),
        ("Thr35", "T", "identical", True),
        ("Pro36", "H", "substituted", False),
        ("Gln62", "T", "substituted", False),
        ("Pro63", "L", "substituted", False),
        ("Tyr66", "H", "substituted", False),
    ]
    headers = [(0.04, "pig"), (0.35, "human"), (0.62, "class")]
    for xh, label in headers:
        axC.text(xh, 0.85, label, fontsize=8, fontweight="bold", color="#555")
    for i, (pig, human, status, identical) in enumerate(anchors):
        y = 0.77 - i * 0.11
        fc = "#e8eef5" if identical else "#f1f1f1"
        axC.add_patch(Rectangle((0.01, y - 0.040), 0.96, 0.08,
                                fc=fc, ec="white", lw=0.8,
                                transform=axC.transAxes))
        axC.text(0.04, y, pig, fontsize=8.5, va="center", fontweight="bold")
        axC.text(0.35, y, human, fontsize=9, va="center", fontweight="bold")
        axC.scatter(0.62, y, s=38, color=C["accent"] if identical else "#999",
                    zorder=3)
        axC.text(0.69, y, status, fontsize=7.8, va="center", color="#444")

    save(fig, "fig3_modeling")


# ===================== FIGURE 5: 2x2 layout =====================
def fig5():
    """Figure 5: 4-panel, 2x2 layout."""
    # data
    t2 = list(csv.DictReader(
        open(DATA / "benchmark" / "tier2_weight_sensitivity.tsv"),
        delimiter="\t"))
    isi = list(csv.DictReader(
        open(DATA / "benchmark" / "isi_weight_sensitivity.tsv"),
        delimiter="\t"))
    res = list(csv.DictReader(
        open(DATA / "benchmark" / "resampling_benchmark_v2.tsv"),
        delimiter="\t"))
    loso = {r["model"]: r for r in csv.DictReader(
        open(DATA / "benchmark" / "loso_benchmark.tsv"), delimiter="\t")}
    with open(DATA / "tier3_surfaceome" / "antigenicity_benchmark_summary.json") as f:
        ag = json.load(f)

    fig = plt.figure(figsize=(15.5, 9.0))

    # ===== Panel A (top-left): Tier-2 weighting =====
    axA = fig.add_axes([0.04, 0.50, 0.46, 0.40])
    labels = [f'{r["w_divergence"]}/{r["w_glyco"]}' for r in t2]
    x = np.arange(len(labels))
    sp = [float(r["spearman_vs_0.6/0.4"]) for r in t2]
    t10 = [float(r["top10_overlap"]) for r in t2]
    t20 = [float(r["top20_overlap"]) for r in t2]
    axA.plot(x, sp, "o-", color="#4a6fa5", lw=2.0, ms=8,
             label="Spearman vs 0.6/0.4")
    axA.plot(x, t10, "s-", color="#7a5195", lw=2.0, ms=8,
             label="Top-10 overlap")
    axA.plot(x, t20, "^-", color="#bc5090", lw=2.0, ms=8,
             label="Top-20 overlap")
    axA.set_xticks(x)
    axA.set_xticklabels(labels, fontsize=9)
    axA.set_xlabel("weight (divergence / glycosylation)")
    axA.set_ylabel("stability")
    axA.set_ylim(-0.05, 1.18)
    axA.set_title("a.", loc="left", fontweight="bold")
    axA.legend(loc="lower left", frameon=False, fontsize=9)
    axA.grid(axis="y", linestyle=":", alpha=0.4)
    axA.spines[["top", "right"]].set_visible(False)

    # ===== Panel B (top-right): ISI weight variants =====
    axB = fig.add_axes([0.54, 0.50, 0.42, 0.40])
    labels_b = [r["variant"] for r in isi]
    sp_b = [float(r["spearman_vs_ISIv2"]) for r in isi]
    t10_b = [float(r["top10_overlap"]) for r in isi]
    xb = np.arange(len(labels_b))
    w = 0.38
    axB.bar(xb - w/2, sp_b, w, color="#4a6fa5",
            label="62-pair rank Spearman")
    axB.bar(xb + w/2, t10_b, w, color="#7a5195",
            label="priority-list Spearman")
    short_b = ["canonical\n0.75-1.25", "wide\n0.5-1.5", "fixed\n1.0",
               "no\nexpression", "no\nfunc_sev", "seq-only\n(no iface)"]
    axB.set_xticks(xb)
    axB.set_xticklabels(short_b, fontsize=9)
    axB.set_ylabel("stability")
    axB.set_ylim(0, 1.18)
    axB.set_title("b.", loc="left", fontweight="bold")
    # Reserve the empty strip above the bars for the key, keeping it clear of
    # the bar tops and annotations.
    axB.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
               frameon=False, fontsize=8.5, borderaxespad=0.0)
    axB.grid(axis="y", linestyle=":", alpha=0.4)
    axB.spines[["top", "right"]].set_visible(False)

    # ===== Panel C (bottom-left): M0-M5 benchmark AUROC =====
    axC = fig.add_axes([0.04, 0.05, 0.46, 0.40])
    models = [r["model"] for r in res]
    labels_c = ["M0\nseq", "M1\n+glyco", "M2\n+iface", "M3\n+gly+iface",
                "M4a\nmolecular\nonly", "M4b\nfunc_sev", "M5\nfunc_sev\nonly"]
    auroc_res = [float(r["auroc_mean"]) for r in res]
    auroc_res_sd = [float(r["auroc_sd"]) for r in res]
    auroc_loso = [float(loso[m]["LOSO_AUROC"]) for m in models]
    xc = np.arange(len(models))
    axC.bar(xc - 0.18, auroc_res, width=0.36, yerr=auroc_res_sd, capsize=4,
            color="#4a6fa5", label="repeated 80/20 resampling",
            error_kw={"lw": 1.2})
    axC.bar(xc + 0.18, auroc_loso, width=0.36, color="#7a5195",
            label="leave-one-system-out")
    axC.axhline(0.5, ls=":", color="black", lw=0.8, alpha=0.6)
    axC.set_xticks(xc)
    axC.set_xticklabels(labels_c, fontsize=8.5)
    axC.set_ylabel("AUROC")
    axC.set_ylim(0, 1.18)
    axC.set_title("c.", loc="left", fontweight="bold")
    axC.legend(loc="upper left", frameon=False, fontsize=9)
    axC.grid(axis="y", linestyle=":", alpha=0.4)
    axC.spines[["top", "right"]].set_visible(False)

    # ===== Panel D (bottom-right): antigenicity benchmark =====
    axD = fig.add_axes([0.54, 0.05, 0.42, 0.40])
    pos = [(k, v["score"], v.get("rank"))
           for k, v in ag["covered_positives"].items()]
    pos.sort(key=lambda x: -x[1])
    neg = [(k, v["score"]) for k, v in ag["putative_negatives"].items()]
    neg.sort(key=lambda x: -x[1])
    pos_scores = [p[1] for p in pos]
    neg_scores = [n[1] for n in neg]

    np.random.seed(42)
    pos_x = [0 + np.random.uniform(-0.07, 0.07) for _ in pos_scores]
    neg_x = [1 + np.random.uniform(-0.07, 0.07) for _ in neg_scores]
    axD.scatter(pos_x, pos_scores, s=85, color="#bc5090", alpha=0.78,
                edgecolor="white", lw=1.0, zorder=3, label="positive (n=7)")
    axD.scatter(neg_x, neg_scores, s=85, color="#4a6fa5", alpha=0.78,
                edgecolor="white", lw=1.0, zorder=3,
                label="putative negative (n=12)")
    axD.boxplot([pos_scores, neg_scores], positions=[0, 1], widths=0.42,
                patch_artist=True, showfliers=False,
                boxprops=dict(facecolor="#dddddd", alpha=0.45,
                              edgecolor="#666"),
                medianprops=dict(color="#222", lw=1.6))
    # staggered annotations (cleaner than R22)
    label_pos = [
        ("CD46",   9,    0.985, 0.985),
        ("CD59",   15,   0.978, 0.91),
        ("PROCR",  41,   0.954, 0.83),
        ("FAM234A", 158, 0.881, 0.75),
    ]
    for name, rank, score, dy in label_pos:
        axD.annotate(f"{name} (#{rank})", (0.05, score),
                     xytext=(0.45, dy), fontsize=8.5,
                     arrowprops=dict(arrowstyle="-", color="gray", lw=0.6),
                     color="#222")
    axD.text(0.55, 0.65, f"AUROC = {ag['v2']['auroc']:.3f}",
             transform=axD.transAxes, fontsize=9,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                       edgecolor="#bbb", alpha=0.92))
    axD.text(0.55, 0.50, f"95% bootstrap CI [0.833, 1.000]",
             transform=axD.transAxes, fontsize=8.5,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                       edgecolor="#bbb", alpha=0.92))
    axD.text(0.55, 0.35, f"MWU p (1-sided) = {ag['mwu_p_one_sided_exact']:.4g}",
             transform=axD.transAxes, fontsize=8.5,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                       edgecolor="#bbb", alpha=0.92))
    axD.set_xticks([0, 1])
    axD.set_xticklabels(["positives\n(literature)", "putative\nnegatives"],
                        fontsize=10)
    axD.set_ylabel("antigenicity v2 score")
    axD.set_ylim(0, 1.05)
    axD.set_title("d.", loc="left", fontweight="bold")
    axD.legend(loc="upper left", frameon=False, fontsize=9)
    axD.grid(axis="y", linestyle=":", alpha=0.4)
    axD.spines[["top", "right"]].set_visible(False)

    save(fig, "fig5_benchmark")


# ===================== FIGURE 6: improved color palette =====================
def fig6():
    """Figure 6: Clinical pathway time series with improved palette.
    Panel A: heatmap of per-pathway row-z-scores (cross-pathway comparable)
    Panel B: line plot of top-4 pathway raw composite dynamics (within-pathway)
    """
    with open(DATA / "clinical" / "matrix_validation_pathway_timeseries.tsv") as f:
        rdr = csv.reader(f, delimiter="\t")
        header = next(rdr)
        rows = list(rdr)
    timepoints = header[1:]
    pathways = [r[0] for r in rows]
    mat = np.array([[float(v) for v in r[1:]] for r in rows])

    # Compute row z-scores for cross-pathway comparable heatmap (Panel A)
    z = (mat - mat.mean(axis=1, keepdims=True)) / mat.std(axis=1, keepdims=True)

    fig = plt.figure(figsize=(12.5, 7.0))

    # ===== Panel A: row-z-score heatmap =====
    axA = fig.add_axes([0.07, 0.14, 0.46, 0.72])
    vmax = 1.5
    vmin = -1.5
    cmap = plt.get_cmap("RdBu_r")
    im = axA.imshow(z, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            v = z[i, j]
            color = "white" if abs(v) > 0.9 else "#222"
            axA.text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=8,
                     color=color)
    axA.set_xticks(range(len(timepoints)))
    axA.set_xticklabels(timepoints, fontsize=9.5)
    axA.set_yticks(range(len(pathways)))
    short_p = [p.replace("_", " ").replace("infiltrate", "infil.").replace(
        "Anticoag fibrinolysis", "Anticoag.").replace("Adhesion LFA1 axis",
        "Adhesion LFA1").replace("Checkpoint CD47", "CD47 checkpoint")
        for p in pathways]
    axA.set_yticklabels(short_p, fontsize=10)
    axA.set_title("a.", loc="left", fontweight="bold", fontsize=11)
    cbar = fig.colorbar(im, ax=axA, shrink=0.85, pad=0.02)
    cbar.set_label("row z-score", fontsize=8.5)

    # ===== Panel B: line plot showing top 4 pathway raw composite dynamics =====
    axB = fig.add_axes([0.62, 0.14, 0.35, 0.72])
    highlight_paths = ["Complement_cascade", "Adaptive_infiltrate",
                       "Innate_myeloid", "Coagulation"]
    line_colors = {"Complement_cascade": C["f6_comp1"],
                   "Adaptive_infiltrate": C["f6_ada"],
                   "Innate_myeloid": C["f6_inn"],
                   "Coagulation": C["f6_coag1"]}
    x = np.arange(len(timepoints))
    for path in highlight_paths:
        if path in pathways:
            i = pathways.index(path)
            axB.plot(x, mat[i], "o-", lw=2, ms=7,
                     color=line_colors[path],
                     label=path.replace("_", " "))
    axB.axvline(4, ls="--", color="#999", lw=0.8, alpha=0.6)
    axB.text(4, 2.4, "T45\nrejection\npeak", ha="center", fontsize=8,
             color="#666")
    axB.set_xticks(x[::2])
    axB.set_xticklabels([timepoints[k] for k in range(0, len(timepoints), 2)],
                        fontsize=8.5)
    axB.set_ylabel("Raw composite score")
    axB.set_title("b.", loc="left", fontweight="bold", fontsize=11)
    axB.legend(loc="lower left", frameon=False, fontsize=8)
    axB.spines[["top", "right"]].set_visible(False)
    axB.axhline(0, color="#999", lw=0.5, ls="-", alpha=0.4)
    axB.grid(axis="y", linestyle=":", alpha=0.4)

    save(fig, "fig6_clinical")


if __name__ == "__main__":
    import sys
    target = sys.argv[1] if len(sys.argv) > 1 else "all"
    if target in ("fig1", "all"):
        print("=== Figure 1 (5-panel) ===")
        fig1()
    if target in ("fig3", "all"):
        print("=== Figure 3 (a+b panels) ===")
        fig3()
    if target in ("fig5", "all"):
        print("=== Figure 5 (2x2 layout) ===")
        fig5()
    if target in ("fig6", "all"):
        print("=== Figure 6 (improved palette) ===")
        fig6()
    print("Done.")
