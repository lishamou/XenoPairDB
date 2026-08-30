"""Benchmarking framework for the 62 Tier-1 gold-standard pairs.

Implements the seven nested/ablated model configurations (M0-M5) and the two
evaluation protocols of the manuscript (Section 4.9.2):

* Protocol A - 500 repeated stratified 80/20 draws (random seed 42).  In each
  iteration percentile normalisation is computed within the reference fold
  only and test points are mapped onto the reference distribution via
  ``searchsorted`` (:func:`xenopairdb.scoring.percentile_map`).
* Protocol B - leave-one-system-out (LOSO) over the six functional systems;
  out-of-fold predictions are pooled before computing AUROC / AUPRC.

Model ladder (all scores are higher = more incompatible):

* M0  - sequence divergence only:                Pct(d_seq)
* M1  - sequence + glycosylation shift:          Pct(0.5*d_seq + 0.5*dglyco)
* M2  - sequence + interface divergence:         Pct(0.5*d_seq + 0.5*d_iface)
* M3  - sequence + interface + glycosylation:    Pct((d_seq+d_iface+dglyco)/3)
* M4a - molecular-feature layer:                 M3 * w_expr
* M4b - annotation-integrated framework:         M4a * s_f
        (evidence weight w_e deliberately excluded - it is derived from the
        curated evidence tier and would constitute label leakage)
* M5  - annotation-concordance control:          s_f alone

Feature imputation rules (manuscript Section 4.9.2):

* pairs without an interface structure: d_interface := d_sequence;
* glycan-antigen entries (GGTA1/CMAH) and other non-alignable non-SLA entries
  (CD99, UOX): nominal sequence/interface divergence of 50;
* SLA-family interactions lack any orthologous sequence and are excluded from
  the sequence-based models M0-M4a; they enter M4b/M5 only through their
  curated functional severity;
* labels: 1 = incompatible pair (56), 0 = compatibility control (6);
  positive rate 56/62 = 0.903 = AUPRC no-skill baseline.

Note on functional severity: the curated s_f values are used as catalogued
(including non-zero severities of retained-control pairs such as
CD80/CD86-CTLA4 and CD274-PDCD1); this reproduces the published benchmark
exactly (M5 AUROC 0.961 / 0.960 in LOSO / resampling).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .scoring import percentile_map

SLA_GENES = ("SLA-I", "SLA-II")
NO_SKILL_BASELINE = 56 / 62  # 0.903 - positive rate of the benchmark set


# ---------------------------------------------------------------------------
# metrics (dependency-free equivalents of sklearn roc_auc / average_precision)
# ---------------------------------------------------------------------------

def auroc(y_true, scores) -> float:
    """AUROC via its Mann-Whitney rank-sum definition (ties count 1/2)."""
    y = np.asarray(y_true, dtype=float)
    s = np.asarray(scores, dtype=float)
    pos = s[y == 1]
    neg = s[y == 0]
    concordant = 0.0
    for sp in pos:
        concordant += np.sum(sp > neg) + 0.5 * np.sum(sp == neg)
    return concordant / (len(pos) * len(neg))


def auprc(y_true, scores) -> float:
    """Average precision (sklearn ``average_precision_score`` semantics).

    Samples are sorted by decreasing score; precision/recall are evaluated at
    the thresholds where the score changes, and AP is the sum of
    (recall_i - recall_{i-1}) * precision_i.
    """
    y = np.asarray(y_true, dtype=float)
    s = np.asarray(scores, dtype=float)
    order = np.argsort(-s, kind="stable")
    y = y[order]
    s = s[order]
    # evaluate at distinct score thresholds (all ties included together)
    distinct = np.where(np.diff(s) != 0)[0] + 1
    thr_idx = np.r_[distinct, len(y)]
    tp = np.cumsum(y)[thr_idx - 1]
    fp = (thr_idx - tp)
    precision = tp / (tp + fp)
    recall = tp / tp[-1]
    ap = 0.0
    prev_recall = 0.0
    for p, r in zip(precision, recall):
        ap += (r - prev_recall) * p
        prev_recall = r
    return ap


# ---------------------------------------------------------------------------
# feature frame construction
# ---------------------------------------------------------------------------

def _tokens(name: str):
    """Split a composite gene field ('CD80/CD86', 'IFNG(pig)') into symbols."""
    name = re.sub(r"\(.*?\)", "", str(name))
    return [t.strip() for t in re.split(r"[/+,;]", name) if t.strip()]


def load_interface_divergence(path: str) -> dict:
    """Return {gene_symbol: 100 - interface_identity} from the 18-pair
    interface JSON (flagship + extended sections)."""
    with open(path) as fh:
        data = json.load(fh)
    out = {}
    for section in data.values():
        for gene, rec in section.items():
            if isinstance(rec, dict) and "iface_identity" in rec:
                out[gene] = 100.0 - rec["iface_identity"]
    return out


def build_feature_frame(v3_path: str, glyco_path: str,
                        interface_path: str) -> pd.DataFrame:
    """Assemble the per-pair benchmark feature table from delivered data.

    Parameters
    ----------
    v3_path : XenoPairDB_v3_ISIv2.tsv (62 gold-standard pairs)
    glyco_path : XenoPairDB_62pairs_glycosylation.tsv (per-gene sequon shifts)
    interface_path : interface_all_18pairs.json (interface divergences)
    """
    df = pd.read_csv(v3_path, sep="\t")
    gly = pd.read_csv(glyco_path, sep="\t")
    delta = dict(zip(gly.gene, gly.delta_sequon))
    iface = load_interface_divergence(interface_path)

    sla = df.gene_key.isin(SLA_GENES)

    # sequence divergence; nominal 50 for non-alignable non-SLA entries
    df["d_seq"] = 100.0 - df.identity
    df.loc[df.d_seq.isna() & ~sla, "d_seq"] = 50.0

    # interface divergence mapped at the *interaction* level (any resolved
    # partner token), then neutrally imputed with d_sequence
    def _map_iface(row):
        for tok in _tokens(row.pig_gene) + _tokens(row.human_partner):
            if tok in iface:
                return iface[tok]
        return np.nan

    df["d_iface"] = df.apply(_map_iface, axis=1)
    df["d_iface"] = df.d_iface.where(df.d_iface.notna(), df.d_seq)

    # glycosylation shift of the porcine donor-side gene
    df["glyco"] = df.gene_key.map(delta).fillna(0.0)

    # curated functional severity as catalogued (controls keep their
    # catalogued values; CD47-SIRPA has no severity annotation -> 0)
    df["fs"] = df.func_sev.fillna(0.0)
    df["ew"] = df.expr_weight
    df["label"] = (df.evidence != "control").astype(int)
    df["mol"] = 0.5 * df.d_seq + 0.5 * df.d_iface
    df["m3_raw"] = (df.d_seq + df.d_iface + df.glyco) / 3.0
    df["is_sla"] = sla
    return df


MODELS = ("M0", "M1", "M2", "M3", "M4a", "M4b", "M5")


@dataclass
class FoldScores:
    """Scores of the seven models for one reference/test split."""

    labels: dict  # model -> array of test labels
    scores: dict  # model -> array of test scores


def score_fold(ref: pd.DataFrame, test: pd.DataFrame) -> FoldScores:
    """Compute all seven model scores for one reference/test split.

    Percentile normalisation is computed within the reference fold only;
    SLA rows are excluded from M0-M4a and contribute their curated s_f to
    M4b/M5 only.
    """
    ref_ns = ref[~ref.is_sla]
    test_ns = test[~test.is_sla]
    test_sla = test[test.is_sla]

    labels, scores = {}, {}
    yt = test_ns.label.values

    p_seq = percentile_map(ref_ns.d_seq.values, test_ns.d_seq.values)
    m1 = percentile_map(0.5 * ref_ns.d_seq.values + 0.5 * ref_ns.glyco.values,
                        0.5 * test_ns.d_seq.values + 0.5 * test_ns.glyco.values)
    p_mol = percentile_map(ref_ns.mol.values, test_ns.mol.values)
    p_m3 = percentile_map(ref_ns.m3_raw.values, test_ns.m3_raw.values)
    m4a = p_m3 * test_ns.ew.values
    m4b = m4a * test_ns.fs.values

    labels["M0"], scores["M0"] = yt, p_seq
    labels["M1"], scores["M1"] = yt, m1
    labels["M2"], scores["M2"] = yt, p_mol
    labels["M3"], scores["M3"] = yt, p_m3
    labels["M4a"], scores["M4a"] = yt, m4a
    labels["M4b"] = np.r_[yt, test_sla.label.values]
    scores["M4b"] = np.r_[m4b, test_sla.fs.values]
    labels["M5"] = test.label.values
    scores["M5"] = percentile_map(ref.fs.values, test.fs.values)
    return FoldScores(labels, scores)


# ---------------------------------------------------------------------------
# protocol A: repeated stratified 80/20 resampling
# ---------------------------------------------------------------------------

def run_resampling(df: pd.DataFrame, n_iter: int = 500, seed: int = 42,
                   progress: bool = False) -> pd.DataFrame:
    """500x stratified 80/20 resampling benchmark (Protocol A).

    Uses ``sklearn.model_selection.train_test_split`` with a single
    ``numpy.random.RandomState(seed)`` consumed sequentially, which fixes the
    draws exactly and reproduces the published numbers.
    """
    from sklearn.model_selection import train_test_split

    idx_all = np.arange(len(df))
    y_all = df.label.values
    idx_ns = idx_all[~df.is_sla.values]
    y_ns = df.label.values[idx_ns]

    acc = {m: [] for m in MODELS}
    accp = {m: [] for m in MODELS}
    rngs = {m: np.random.RandomState(seed) for m in MODELS}
    for m in MODELS:
        rs = rngs[m]
        for _ in range(n_iter):
            if m in ("M4b", "M5"):
                tr, te = train_test_split(idx_all, test_size=0.2,
                                          stratify=y_all, random_state=rs)
            else:
                tr, te = train_test_split(idx_ns, test_size=0.2,
                                          stratify=y_ns, random_state=rs)
            fold = score_fold(df.iloc[tr], df.iloc[te])
            acc[m].append(auroc(fold.labels[m], fold.scores[m]))
            accp[m].append(auprc(fold.labels[m], fold.scores[m]))
        if progress:
            print(f"  {m} done")
    rows = []
    for m in MODELS:
        rows.append({
            "model": m,
            "auroc_mean": round(float(np.mean(acc[m])), 3),
            "auroc_sd": round(float(np.std(acc[m])), 3),
            "auprc_mean": round(float(np.mean(accp[m])), 3),
            "auprc_sd": round(float(np.std(accp[m])), 3),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# protocol B: leave-one-system-out
# ---------------------------------------------------------------------------

def run_loso(df: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-system-out benchmark (Protocol B).

    Each of the six functional systems is held out in turn; out-of-fold
    predictions are pooled before computing AUROC / AUPRC.
    """
    pooled_y = {m: [] for m in MODELS}
    pooled_s = {m: [] for m in MODELS}
    for system in sorted(df.system.unique()):
        ref = df[df.system != system]
        test = df[df.system == system]
        fold = score_fold(ref, test)
        for m in MODELS:
            pooled_y[m].extend(fold.labels[m])
            pooled_s[m].extend(fold.scores[m])
    rows = []
    for m in MODELS:
        rows.append({
            "model": m,
            "LOSO_AUROC": round(auroc(pooled_y[m], pooled_s[m]), 3),
            "LOSO_AUPRC": round(auprc(pooled_y[m], pooled_s[m]), 3),
        })
    return pd.DataFrame(rows)
