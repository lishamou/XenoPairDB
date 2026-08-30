#!/usr/bin/env python3
"""Regenerate Figure 2 (4 panels a-d) from frozen TSVs only — no recompute.

R29 v20-review fixes:
- Restore panels a, b, c from frozen `XenoPairDB_v5_atlas_1407pairs.tsv`
- Panel d uses the frozen donor residual matrix and local genotype helpers.
- Identity units verified: columns are 0-100 (% identity); divergence
  = 100 - mean identity (matches TSV `pair_divergence` directly)
"""
import os
import csv
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# R180: all production inputs are owned by the pipeline; submission/ is read-only history.
_PIPELINE_ROOT = Path(__file__).resolve().parents[2]   # pipeline/
_PROJECT_ROOT  = Path(__file__).resolve().parents[3]   # project root (28XenoPairDB/)
ROOT = Path(os.environ.get("XENOPAIR_PROJECT_ROOT", str(_PROJECT_ROOT)))
DATA = _PIPELINE_ROOT / "release_data"
OUT_NBT = _PIPELINE_ROOT / "figures"   # R110-A: write to pipeline/figures/
PAYLOAD_DIR = _PIPELINE_ROOT / "figure_payloads"
sys.path.insert(0, str(_PIPELINE_ROOT / "src"))

C = dict(
    bar="#c96f4a",
    blue="#2f4858",
    cmap="YlOrBr",
    scatter="#9bb5b1",
    accent="#34495e",
    highlight="#3d6b5e",
)


def save(fig, name, formats=("png", "pdf", "svg")):
    for ext in formats:
        out = OUT_NBT / f"{name}.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {name}: {', '.join(formats)}")


def load_atlas():
    rows = []
    with open(_PIPELINE_ROOT / "canonical_links" / "tier2_atlas.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                r["ligand_identity_pct"] = float(r["ligand_identity_pct"])
                r["receptor_identity_pct"] = float(r["receptor_identity_pct"])
                r["pair_divergence"] = float(r["pair_divergence"])
                r["tier2_score"] = float(r["tier2_score"])
                r["glyco_delta_pig_minus_human"] = float(
                    r["glyco_delta_pig_minus_human"])
                rows.append(r)
            except (ValueError, KeyError):
                continue
    return rows


def load_tier1():
    import pandas as pd
    df = pd.read_csv(_PIPELINE_ROOT / "canonical_links" / "tier1_frozen.tsv",
                     sep="\t")
    df["ISI_num"] = pd.to_numeric(df["ISI"], errors="coerce").fillna(0)
    return df


def fig2():
    atlas = load_atlas()
    t1 = load_tier1()

    fig = plt.figure(figsize=(14.6, 8.6))
    # Equal-size 2 × 2 layout keeps the information-dense a–c panels readable
    # while giving the compact residual matrix only one panel footprint.
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.42,
                          wspace=0.34)

    # ============ Panel a: identity distribution (frozen) ============
    axa = fig.add_subplot(gs[0, 0])
    min_ids = [min(r["ligand_identity_pct"], r["receptor_identity_pct"])
               for r in atlas]
    axa.hist(min_ids, bins=60, color=C["bar"], alpha=0.85,
             edgecolor="white", linewidth=0.4)
    axa.axvline(70, ls="--", color=C["blue"], lw=1.4,
                label="70% identity threshold")
    axa.set_xlim(0, 100)
    axa.set_xlabel("min(partner identity %) per pair", fontsize=10)
    axa.set_ylabel("pairs (n = {:,})".format(len(atlas)), fontsize=10)
    axa.set_title("a.", fontsize=11, fontweight="bold", loc="left")
    axa.legend(loc="upper right", frameon=False, fontsize=9)
    axa.spines[["top", "right"]].set_visible(False)
    axa.grid(axis="y", linestyle=":", alpha=0.3)

    # ============ Panel b: divergence × glyco shift (frozen) ============
    axb = fig.add_subplot(gs[0, 1])
    div = [r["pair_divergence"] for r in atlas]
    gly = [r["glyco_delta_pig_minus_human"] for r in atlas]
    sc = axb.scatter(div, gly, c=[r["tier2_score"] for r in atlas],
                     s=10, cmap="viridis", alpha=0.8,
                     edgecolor="white", linewidth=0.2)
    fig.colorbar(sc, ax=axb, label="Tier-2 score", shrink=0.85)
    axb.set_xlabel("pair divergence (%)", fontsize=10)
    axb.set_ylabel("Δ N-glyco sequons (pig − human)", fontsize=10)
    axb.set_title("b.", fontsize=11, fontweight="bold", loc="left")
    axb.spines[["top", "right"]].set_visible(False)
    axb.grid(linestyle=":", alpha=0.3)

    # ============ Panel c: Top-15 (frozen) ============
    axc = fig.add_subplot(gs[1, 0])
    top = sorted(atlas, key=lambda r: -r["tier2_score"])[:15][::-1]
    labels = [f'{r["ligand"]}–{r["receptor"]}' for r in top]
    scores = [r["tier2_score"] for r in top]
    ypos = np.arange(len(top))
    axc.barh(ypos, scores, color=C["bar"], edgecolor="white", linewidth=0.5)
    for i, v in enumerate(scores):
        axc.text(v + 0.4, i, f"{v:.1f}", va="center", fontsize=8)
    axc.set_yticks(ypos)
    axc.set_yticklabels(labels, fontsize=8.5)
    axc.set_xlabel("Tier-2 candidate score", fontsize=10)
    axc.set_title("c.", fontsize=11, fontweight="bold", loc="left")
    axc.spines[["top", "right"]].set_visible(False)
    axc.grid(axis="x", linestyle=":", alpha=0.3)

    # ============ Panel d: residual ISI matrix by genotype (R121.8) ============
    # R121.8 (V6 P0-1): reads from facts.donor_scorecard_matrix (keyed, single source)
    # Shows canonical post-edit masked per-system residual ISI v2
    # (heatmap of 6 systems × the 3 clinical-context donors in the frozen matrix).
    # Masked totals (356/334/424) are shown in Fig 1e and SI Table S1.
    import json as _json
    _facts_path = Path(__file__).resolve().parents[2] / "manuscript" / "manuscript_facts.json"
    _facts = _json.loads(_facts_path.read_text(encoding="utf-8")) if _facts_path.exists() else {}
    # Single canonical donor matrix shared with Fig1 and FigS1.
    import pandas as pd
    matrix_path = _PIPELINE_ROOT / "canonical_links" / "donor_residual_matrix.tsv"
    systems = ["Complement", "Coagulation", "Immune", "Cytokine", "Endocrine", "Adhesion"]
    mat = pd.read_csv(matrix_path, sep="\t").set_index("donor")[[f"masked_{s}" for s in systems]]
    mat.columns = systems
    PAYLOAD_DIR.mkdir(exist_ok=True)
    (PAYLOAD_DIR / "fig2d.json").write_text(__import__('json').dumps({
        "figure": "Fig2d", "metric_variant": "post_edit_masked",
        "systems": systems, "matrix": mat.astype(int).to_dict(orient="index"),
        "totals": mat.sum(axis=1).astype(int).to_dict()}, indent=2) + "\n")

    axd = fig.add_subplot(gs[1, 1])
    im = axd.imshow(mat.values, cmap=C["cmap"], aspect="auto",
                    vmin=0, vmax=520)
    axd.set_xticks(range(mat.shape[1]))
    axd.set_xticklabels(["Comp.", "Coag.", "Immune", "Cytokine", "Endo.", "Adhesion"],
                        fontsize=9, rotation=25, ha="right")
    axd.set_yticks(range(mat.shape[0]))
    axd.set_yticklabels(mat.index, fontsize=10)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            axd.text(j, i, str(v), ha="center", va="center",
                     fontsize=10,
                     color="white" if v > 320 else "#222222",
                     fontweight="bold")
    fig.colorbar(im, ax=axd, label="masked residual ISI v2", shrink=0.85)
    axd.set_title("d.", fontsize=10.5, fontweight="bold", loc="left")
    axd.set_xlabel("functional system", fontsize=10)
    axd.set_ylabel("donor genotype", fontsize=10)
    save(fig, "fig2_atlas")


if __name__ == "__main__":
    fig2()
