#!/usr/bin/env python3
"""Regenerate Figure 4 (R29 v20-corrected) from the *corrected* surfaceome TSV
(`pig_surfaceome_antigenicity_v2_glyco_corrected.tsv`, 4,123 rows including
CD46/CD55/CD59/SLA-1/SLA-DRA).
"""
import os
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# R110-A: pipeline self-contained — derive paths from script location
_PIPELINE_ROOT = Path(__file__).resolve().parents[2]   # pipeline/
_PROJECT_ROOT  = Path(__file__).resolve().parents[3]   # project root (28XenoPairDB/)
ROOT = Path(os.environ.get("XENOPAIR_PROJECT_ROOT", str(_PROJECT_ROOT)))
DATA = _PIPELINE_ROOT / "release_data"
OUT_NBT = _PIPELINE_ROOT / "figures"   # R110-A: write to pipeline/figures/

C = dict(
    scatter="#9dbb8a",
    scatter_alpha=0.55,
    icam1="#3d6b5e",
    cd86="#4a7a6a",
    slamf7="#5b8a72",
    fam234a="#c96f2a",
    accent="#34495e",
    complement="#8a9b6e",
    cd46="#3a7ca5",
    cd55="#5b8a72",
    cd59="#a35d6a",
    sla1="#7a6a8a",
)


def save(fig, name, formats=("png", "pdf", "svg")):
    for ext in formats:
        out = OUT_NBT / f"{name}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {name}: {', '.join(formats)}")


def fig4():
    """Regenerate Figure 4 from R31 v2 corrected surfaceome."""
    with open(DATA / "tier3_surfaceome" /
              "pig_surfaceome_antigenicity_v2_glyco_corrected.tsv") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        rows = list(rdr)
    n_total = len(rows)
    n_scoreable = sum(1 for r in rows if r["rank_v2"] != "")
    print(f"Membrane universe: {n_total}, scoreable subset: {n_scoreable}")

    rank_lookup = {r["gene"]: int(r["rank_v2"])
                   for r in rows if r["rank_v2"] != ""}

    fig = plt.figure(figsize=(12.5, 7.0))

    # ===== Main scatter panel =====
    ax = fig.add_axes([0.08, 0.12, 0.84, 0.76])
    xs, ys = [], []
    for r in rows:
        try:
            xs.append(float(r["divergence"]))
            ys.append(float(r["p_sequon_density"]))
        except ValueError:
            continue
    ax.scatter(xs, ys, s=10, color=C["scatter"], alpha=C["scatter_alpha"],
               edgecolor="white", linewidth=0.3,
               label=f"sequence-scoreable proteins (n={n_scoreable:,})",
               zorder=2)

    # Highlight key candidates — ranks in the scoreable subset
    marks = {
        "ICAM1": ("icam1",
                  "ICAM1\n#{} overall, #1 among\nadhesion molecules".format(
                      rank_lookup.get("ICAM1", "?"))),
        "FAM234A": ("fam234a",
                    "FAM234A\nrank {} / {} scoreable (top {:.1f}%)".format(
                        rank_lookup.get("FAM234A", "?"), n_scoreable,
                        rank_lookup.get("FAM234A", 1) / n_scoreable * 100)),
        "CD86": ("cd86",
                 "CD86\nrank {} (control)".format(
                     rank_lookup.get("CD86", "?"))),
        "SLAMF7": ("slamf7",
                   "SLAMF7\nrank {} / scoreable".format(
                       rank_lookup.get("SLAMF7", "?"))),
    }
    label_offsets = {
        "ICAM1": (12, 8),
        "CD86": (-70, -48),
        "SLAMF7": (14, -4),
        "FAM234A": (12, 8),
    }
    pos = {g: next(((float(r["divergence"]), float(r["p_sequon_density"]))
                    for r in rows if r["gene"] == g), (None, None))
           for g in marks}
    for gene, (ckey, label) in marks.items():
        xi, yi = pos[gene]
        if xi is None:
            continue
        ax.scatter(xi, yi, s=190, color=C[ckey], edgecolor="#222",
                   linewidth=1.2, zorder=5)
        ax.annotate(label, (xi, yi), xytext=label_offsets.get(gene, (12, 8)),
                    textcoords="offset points", fontsize=8.8,
                    fontweight="bold", color=C["accent"])

    # Newly added canonical markers — coloured diamonds with prominent labels
    new_marks = [
        ("CD46", C["cd46"],
         "CD46\n#{} overall".format(rank_lookup.get("CD46", "?"))),
        ("CD55", C["cd55"],
         "CD55\n#{} overall".format(rank_lookup.get("CD55", "?"))),
        ("CD59", C["cd59"],
         "CD59\n#{} overall".format(rank_lookup.get("CD59", "?"))),
    ]
    for gene, color, label in new_marks:
        xi, yi = pos.get(gene, (None, None))
        if xi is None:
            continue
        ax.scatter(xi, yi, s=190, marker="D", color=color,
                   edgecolor="#222", linewidth=1.2, zorder=6, alpha=0.95)
        ax.annotate(label, (xi, yi), xytext=(14, -22),
                    textcoords="offset points", fontsize=8.5,
                    color=color, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.25",
                              fc="white", ec=color, lw=0.8, alpha=0.9))

    ax.set_xlabel("Sequence divergence (%)", fontsize=11)
    ax.set_ylabel("N-glycan sequon density", fontsize=11)
    ax.text(0.0, 1.02, "a.", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linestyle=":", alpha=0.3)

    save(fig, "fig4_antigenicity")


if __name__ == "__main__":
    fig4()
