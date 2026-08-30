"""Surface-protein antigenicity v2 screen and external literature benchmark.

The antigenicity layer scores all 4,118 porcine surface proteins as

    A_v2 = 0.5 * Pct(d_sequence) + 0.5 * Pct(glyco_density)

(percentiles within the surfaceome).  The external benchmark evaluates this
score against 7 scoreable literature-validated non-Gal antigens (five present
in the surfaceome - CD9, PROCR, FAM234A, CFTR, CDH1 - plus CD46 and CD59,
which are absent from the porcine UniProt surfaceome owing to gene-annotation
gaps and are rescued by back-computing their percentiles with the same
scoring function) versus 12 putative negatives (essential membrane proteins
with human-pig identity >= 95%).

Published results (manuscript Fig. 5c / Supplementary Fig. S2):
AUROC 0.952, AUPRC 0.924, Mann-Whitney P = 2.4e-4 (one-sided exact;
two-sided exact P = 4.8e-4), positive median 0.881, negative median 0.280.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .benchmark import auroc, auprc
from .scoring import antigenicity_v2, backcompute_v2

# 7 scoreable literature-validated non-Gal antigens (positives)
SURFACEOME_POSITIVES = ("CD9", "PROCR", "FAM234A", "CFTR", "CDH1")
RESCUED_POSITIVES = ("CD46", "CD59")  # rescued via Tier-1 back-computation

# 12 putative negatives: essential membrane proteins, human-pig identity
# >= 95% (CLDN1 borderline at 94.8%; AQP1/ATP1B1 included as reference
# borderline entries), all present in the 4,118-protein surfaceome.
PUTATIVE_NEGATIVES = (
    "ATP1A1", "ATP1A3", "ATP2B1", "KCNJ10", "SLC2A1", "GJA1",
    "GJB1", "SLC9A1", "SLC12A2", "CLDN1", "AQP1", "ATP1B1",
)


def load_surfaceome(path: str) -> pd.DataFrame:
    """Load the surfaceome table and recompute antigenicity v2 from raw
    features (sequence divergence and glycoepitope density)."""
    s = pd.read_csv(path, sep="\t")
    s["v2_calc"] = antigenicity_v2(s.divergence.values, s.glyco_density.values)
    # rank 1 = highest antigenicity (min-rank convention for ties)
    s["rank_calc"] = pd.Series(-s.v2_calc).rank(method="min").astype(int)
    return s


def rescue_antigens(surfaceome: pd.DataFrame, rescue_json: str) -> pd.DataFrame:
    """Re-derive the antigenicity v2 scores of CD46/CD59 from their Tier-1
    identity and sequon counts, back-computed against the surfaceome
    distributions, and return the augmented table with new ranks."""
    with open(rescue_json) as fh:
        rescue = json.load(fh)
    rows = []
    for gene, rec in rescue.items():
        divergence = 100.0 - rec["identity"]
        glyco_density = rec["glyco_density"]  # sites per 100 aa
        v2 = backcompute_v2(surfaceome.divergence.values,
                            surfaceome.glyco_density.values,
                            divergence, glyco_density)
        # published ranks ("rank k of 4,118") are against the surfaceome
        # itself, before any other rescued antigen is merged in
        rank_rescue = int(1 + np.sum(surfaceome.v2_calc.values > v2))
        rows.append({"gene": gene, "divergence": divergence,
                     "glyco_density": glyco_density, "v2_calc": v2,
                     "rank_rescue": rank_rescue, "rescued": True})
    out = surfaceome.copy()
    out["rescued"] = False
    out["rank_rescue"] = np.nan
    out = pd.concat([out, pd.DataFrame(rows)], ignore_index=True)
    out["rank_aug"] = pd.Series(-out.v2_calc).rank(method="min").astype(int)
    return out


def run_antigen_benchmark(augmented: pd.DataFrame) -> dict:
    """Evaluate the 7-positive / 12-putative-negative literature benchmark."""
    pos_genes = list(SURFACEOME_POSITIVES) + list(RESCUED_POSITIVES)
    pos = augmented[augmented.gene.isin(pos_genes)]
    neg = augmented[augmented.gene.isin(PUTATIVE_NEGATIVES)]
    assert len(pos) == 7, f"expected 7 positives, got {len(pos)}"
    assert len(neg) == 12, f"expected 12 negatives, got {len(neg)}"

    pos_scores = pos.v2_calc.values
    neg_scores = neg.v2_calc.values
    y = np.r_[np.ones(len(pos_scores)), np.zeros(len(neg_scores))]
    s = np.r_[pos_scores, neg_scores]

    u_two = mannwhitneyu(pos_scores, neg_scores,
                         alternative="two-sided", method="exact")
    u_one = mannwhitneyu(pos_scores, neg_scores,
                         alternative="greater", method="exact")
    return {
        "positives": pos[["gene", "v2_calc", "rank_aug"]].sort_values("rank_aug"),
        "negatives": neg[["gene", "v2_calc", "rank_aug"]].sort_values("rank_aug"),
        "auroc": auroc(y, s),
        "auprc": auprc(y, s),
        "mwu_p_two_sided_exact": u_two.pvalue,
        "mwu_p_one_sided_exact": u_one.pvalue,
        "pos_median": float(np.median(pos_scores)),
        "neg_median": float(np.median(neg_scores)),
    }
