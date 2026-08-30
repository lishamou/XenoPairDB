"""R107 Generate Supplementary Figure S2 — antigenicity benchmark (Release 3.2).

3-panel layout (12 × 5 inch, healthy aspect ratio ~2.4:1):
- panel a: Benchmark construction workflow (text + light flow)
- panel b: 8 positives on the 4,122-protein antigenicity v2 distribution
  (background distribution with positive labels)
- panel c: Boxplots — literature positives (n=8) vs high-conservation reference set (n=12)
"""
import os
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# R110-A: pipeline self-contained — derive paths from script location
_PIPELINE_ROOT = Path(__file__).resolve().parents[2]   # pipeline/
_PROJECT_ROOT  = Path(__file__).resolve().parents[3]   # project root (28XenoPairDB/)
ROOT = Path(os.environ.get("XENOPAIR_PROJECT_ROOT", str(_PROJECT_ROOT)))
OUT_EN = _PIPELINE_ROOT / "figures"     # R110-A: write to pipeline/figures/

with open(_PIPELINE_ROOT / "canonical_links" / "antigenicity_benchmark_summary.json") as f:
    bench = json.load(f)

pos = bench["covered_positives"]
neg = bench["putative_negatives"]
pos_items = sorted(pos.items(), key=lambda kv: -kv[1]["score"])
neg_items = sorted(neg.items(), key=lambda kv: -kv[1]["score"])
pos_scores = np.array([v["score"] for v in pos.values()])
neg_scores = np.array([v["score"] for v in neg.values()])
pos_median = float(np.median(pos_scores))
neg_median = float(np.median(neg_scores))
auroc = bench["v2"]["auroc"]
auprc = bench["v2"]["auprc"]
mwu = bench["mwu_p_one_sided_exact"]
ci_95 = bench.get("CI_95_bootstrap", "[0.846, 1.000]")

# Figure — 12 × 5 inch, 3 panels (1 row, 3 cols) with width ratios 1.3 / 1.3 / 1.4
fig, axes = plt.subplots(
    1, 3, figsize=(12, 5.6),
    gridspec_kw={"width_ratios": [1.05, 1.5, 1.45], "wspace": 0.32}
)

# =================== Panel a: construction workflow ===================
ax_a = axes[0]
ax_a.axis("off")
ax_a.set_xlim(0, 1)
ax_a.set_ylim(0, 1)
ax_a.text(0.04, 0.96, "a.",
          fontsize=11, fontweight="bold", va="top")

# Boxes
box_x, box_w, box_h = 0.05, 0.90, 0.11
box_color = "#fbf3df"
ec = "#888"
box_specs = [
    (0.84, "22 publications →\n13 literature-positive antigens", "#fbf3df"),
    (0.69, "13 literature-positive antigens", "#fbf3df"),
    (0.54, "8 scoreable in 4,122-protein\nsubset of the 4,123-protein\nmembrane universe", "#fbf3df"),
    (0.39, "5 excluded:\nSLA-I, SLA-II, ENO1, ANXA2, VIM", "#fce4ec"),
    (0.24, "12 high-conservation reference\nproteins (≥95% identity)", "#fde9d9"),
]
for y, text, color in box_specs:
    ax_a.add_patch(FancyBboxPatch((box_x, y - box_h / 2), box_w, box_h,
                                   boxstyle="round,pad=0.02",
                                   fc=color, ec=ec, lw=0.8))
    ax_a.text(box_x + box_w / 2, y, text, ha="center", va="center",
              fontsize=8.4)
# Arrows
for y_top, y_bot in zip([x[0] for x in box_specs[:-1]],
                         [x[0] for x in box_specs[1:]]):
    ax_a.annotate("", xy=(box_x + box_w / 2, y_bot + box_h / 2 + 0.006),
                  xytext=(box_x + box_w / 2, y_top - box_h / 2 - 0.006),
                  arrowprops=dict(arrowstyle="-|>", color="#444", lw=1.0))


# =================== Panel b: 8 positives on 4,122 distribution ===================
ax_b = axes[1]
surfaceome = pd.read_csv(_PIPELINE_ROOT / "canonical_links" / "surfaceome.tsv", sep="\t")
full_bg = surfaceome.loc[surfaceome["rank_v2"].notna(), "antigenicity_v2"].astype(float).to_numpy()
if len(full_bg) != 4122:
    raise ValueError(f"Expected 4,122 scoreable surfaceome scores, got {len(full_bg)}")

ax_b.hist(full_bg, bins=50, color="#cfd8dc", edgecolor="white",
          linewidth=0.4, label=f"{len(full_bg):,}-protein distribution")

# Mark the 8 positive positions as vertical lines (rank-ordered)
# Sort by descending score for clarity
pos_sorted = sorted(pos.items(), key=lambda kv: -kv[1]["score"])
for i, (gene, info) in enumerate(pos_sorted):
    s = info["score"]
    ax_b.axvline(s, ymin=0.85, ymax=0.99,
                 color="#c0392b", lw=1.6, alpha=0.85,
                 )
# The red lines encode all eight positive references; omit gene labels here
# because the right-tail scores are too close for readable annotations.
ax_b.text(0.98, 0.72, "8 literature-positive\nantigens",
          transform=ax_b.transAxes, ha="right", va="top", fontsize=7.5,
          color="#7a1f12", fontweight="bold")

# 12 high-conservation reference proteins as small dots near the median
ax_b.scatter(neg_scores, [0.02] * len(neg_scores),
             s=22, color="#1f77b4", alpha=0.65,
             edgecolor="#0b3d62", linewidth=0.5,
             label=f"12 high-conservation references (median {neg_median:.3f})", zorder=3)

ax_b.set_xlim(0, 1.02)
ax_b.set_xlabel("Antigenicity v2 score", fontsize=9)
ax_b.set_ylabel("Protein count", fontsize=9)
ax_b.set_title("b.", fontsize=10, fontweight="bold", loc="left")
# The caption defines the two distributions; omitting a separate legend keeps
# the top-tail positive markers and panel label unobstructed.
ax_b.spines[["top", "right"]].set_visible(False)

# =================== Panel c: boxplot literature positives vs reference set ===================
ax_c = axes[2]

# Boxplot
bp = ax_c.boxplot([neg_scores, pos_scores],
                  labels=[f"High-conservation\nreference set\n(n=12)",
                          f"Literature positives\n(n=8)"],
                  patch_artist=True, widths=0.55,
                  medianprops={"color": "#222", "linewidth": 2},
                  whiskerprops={"color": "#444"},
                  capprops={"color": "#444"},
                  flierprops={"marker": "o", "markerfacecolor": "#c0392b",
                              "markeredgecolor": "#7a1f12", "markersize": 5})
for patch, color in zip(bp["boxes"], ["#cfe2f3", "#f9d6cf"]):
    patch.set_facecolor(color)
    patch.set_edgecolor("#444")

# Overlay jittered points
rng = np.random.default_rng(42)
for i, (vals, color) in enumerate(
        [(neg_scores, "#1f77b4"), (pos_scores, "#c0392b")]):
    x = rng.normal(loc=i + 1, scale=0.06, size=len(vals))
    ax_c.scatter(x, vals, s=35, alpha=0.8, color=color,
                 edgecolor="white", linewidth=0.6, zorder=4)

ax_c.set_ylabel("Antigenicity v2 score", fontsize=9)
ax_c.set_title("c.", fontsize=10, fontweight="bold", loc="left")
ax_c.set_ylim(0.0, 1.05)
ax_c.spines[["top", "right"]].set_visible(False)
ax_c.grid(axis="y", linestyle=":", alpha=0.4)

# Annotate medians
ax_c.text(1, 0.36, f"median {neg_median:.3f}",
          ha="center", va="bottom", fontsize=8, color="#0b3d62")
ax_c.text(2, 1.01, f"median {pos_median:.3f}",
          ha="center", va="bottom", fontsize=8, color="#7a1f12",
          fontweight="bold")

plt.tight_layout()
OUT_EN.mkdir(parents=True, exist_ok=True)
for ext in ["png", "pdf", "svg"]:
    out = OUT_EN / f"figS2_antigen_benchmark.{ext}"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  wrote {out}")

plt.close(fig)
print(f"figS2: AUROC={auroc}, AUPRC={auprc}, n+={len(pos_scores)}, n-={len(neg_scores)}")
print(f"Pos median={pos_median:.3f}, Neg median={neg_median:.3f}")
print(f"Ranks: CD46 #1, CD55 #4, CD59 #17, PROCR #44, FAM234A #161")
