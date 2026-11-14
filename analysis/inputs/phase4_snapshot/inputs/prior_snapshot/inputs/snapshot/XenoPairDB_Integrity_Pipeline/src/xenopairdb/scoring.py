"""Core scoring functions of XenoPairDB.

All formulas are pure functions with no fitted parameters, exactly as defined
in the manuscript (Sections 2.2 / 4.4–4.6):

* ISI v1 / v2 (Tier-1, 62 evidence-curated pairs):

      ISI_v1 = w_e * s_f * mol_div
      ISI_v2 = ISI_v1 * w_expr

  with mol_div = 0.5 * d_interface + 0.5 * d_sequence (interface divergence
  replaced by d_sequence when no interface structure is available), a nominal
  mol_div = 50 for glycan/SLA-class entries with no alignable sequence, and

      w_expr = clip(0.75 + 0.5 * max(0, rho), 0.75, 1.25)

  where rho is the cross-species expression Spearman correlation.  For
  low-count inducible genes (e.g. SELE, IL1B-pathway inducible entries,
  ULBP1, IFNG) rho is not interpretable and a neutral w_expr = 1.0 is used.

* Tier-2 atlas score (1,407 predicted ligand-receptor pairs):

      Tier2 = (0.6 * Pct(d_pair) + 0.4 * Pct(delta_glyco)) * 100

  percentiles computed within the atlas (0-1).

* Surface-protein antigenicity v2 (4,118 porcine surface proteins):

      A_v2 = 0.5 * Pct(d_sequence) + 0.5 * Pct(glyco_density)

  percentiles computed within the surfaceome (0-1, reported as fraction).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# generic percentile helpers
# ---------------------------------------------------------------------------

def percentile_rank(values) -> np.ndarray:
    """Percentile rank (0-1] of *values* within themselves.

    Equivalent to ``pandas.Series(values).rank(pct=True)`` (mean ranks for
    ties), which is the convention used for all in-atlas / in-surfaceome
    percentiles in the manuscript.
    """
    return pd.Series(np.asarray(values, dtype=float)).rank(pct=True).values


def percentile_map(reference, x) -> np.ndarray:
    """Map query points *x* onto the percentile scale of a reference sample.

    Percentiles are computed **within the reference fold only**; test points
    are mapped back onto the reference distribution with
    ``numpy.searchsorted`` ("left" convention: fraction of reference values
    strictly smaller than x).  This is the leakage-free normalisation used in
    the 80/20-resampling and LOSO benchmark protocols.
    """
    ref = np.sort(np.asarray(reference, dtype=float))
    return np.searchsorted(ref, np.asarray(x, dtype=float), side="left") / len(ref)


# ---------------------------------------------------------------------------
# Tier-1: ISI v1 / v2
# ---------------------------------------------------------------------------

NOMINAL_DIVERGENCE = 50.0  # nominal mol_div for glycan/SLA-class entries


def sequence_divergence(identity_pct: float) -> float:
    """d_sequence = 100 - identity(%) for a cross-species alignment."""
    return 100.0 - identity_pct


def molecular_divergence(identity_pct, iface_identity_pct=None) -> float:
    """Composite sequence-structure divergence (mol_div).

    mol_div = 0.5 * d_interface + 0.5 * d_sequence  when an interface
    structure is available, otherwise d_sequence alone (neutral setting).
    Entries with no alignable sequence (glycan antigens, SLA class) receive
    the nominal divergence of 50.
    """
    if identity_pct is None or (isinstance(identity_pct, float) and np.isnan(identity_pct)):
        return NOMINAL_DIVERGENCE
    d_seq = sequence_divergence(identity_pct)
    if iface_identity_pct is None or (isinstance(iface_identity_pct, float)
                                      and np.isnan(iface_identity_pct)):
        return d_seq
    return 0.5 * (100.0 - iface_identity_pct) + 0.5 * d_seq


def expression_weight(rho, inducible: bool = False,
                      lo: float = 0.75, hi: float = 1.25) -> float:
    """Expression-concordance weight.

        w_expr = clip(lo + (hi - lo) * max(0, rho), lo, hi)

    The canonical interval is [0.75, 1.25]; the wide sensitivity variant uses
    [0.5, 1.5].  Genes with extremely low constitutive expression in both
    species that are induced under physiological conditions (``inducible``)
    receive the neutral weight 1.0, as do entries without an interpretable
    rho (NaN).
    """
    if inducible or rho is None or (isinstance(rho, float) and np.isnan(rho)):
        return 1.0
    w = lo + (hi - lo) * max(0.0, rho)
    return float(np.clip(w, lo, hi))


def isi_v2(evid_w: float, func_sev: float, mol_div: float,
           expr_weight: float = 1.0) -> float:
    """Incompatibility Severity Index v2.

        ISI_v2 = w_e * s_f * mol_div * w_expr

    w_e  : evidence weight (1.0 direct / 0.6 disputed / 0.4 inferred /
           0.0 control)
    s_f  : functional severity in [0, 1] (0 for compatibility controls)
    mol_div : molecular divergence from :func:`molecular_divergence`
    expr_weight : from :func:`expression_weight`

    Accepts scalars or numpy/pandas vectors (elementwise).
    """
    return np.asarray(evid_w) * np.asarray(func_sev) * np.asarray(mol_div) \
        * np.asarray(expr_weight)


def isi_v1(evid_w: float, func_sev: float, mol_div: float) -> float:
    """Incompatibility Severity Index v1 (legacy, without expression weight).

        ISI_v1 = w_e * s_f * mol_div

    v1 is the v2 formula without the bounded expression-relevance
    multiplier (default w_expr = 1.0). Maintained for backward
    compatibility with the donor scorecard's residual coverage
    definition (manuscript Methods: coverage = sum-of-fixed-ISI / total-ISI,
    both v1).

    Accepts scalars or numpy/pandas vectors (elementwise).
    """
    return np.asarray(evid_w) * np.asarray(func_sev) * np.asarray(mol_div)


# ---------------------------------------------------------------------------
# Tier-2 atlas score
# ---------------------------------------------------------------------------

def tier2_score(pair_divergence, glyco_delta,
                w_div: float = 0.6, w_glyco: float = 0.4) -> np.ndarray:
    """Tier-2 atlas score (0-100).

        Tier2 = (w_div * Pct(d_pair) + w_glyco * Pct(delta_glyco)) * 100

    Percentiles are computed within the atlas (all 1,407 pairs).
    """
    p_div = percentile_rank(pair_divergence)
    p_gly = percentile_rank(glyco_delta)
    return (w_div * p_div + w_glyco * p_gly) * 100.0


# ---------------------------------------------------------------------------
# Tier-3: surface-protein antigenicity v2
# ---------------------------------------------------------------------------

def antigenicity_v2(divergence, glyco_density) -> np.ndarray:
    """Antigenicity v2 score in [0, 1].

        A_v2 = 0.5 * Pct(d_sequence) + 0.5 * Pct(glyco_density)

    Percentiles are computed within the 4,118-protein porcine surfaceome.
    """
    return 0.5 * percentile_rank(divergence) + 0.5 * percentile_rank(glyco_density)


def backcompute_v2(divergence_ref, glyco_density_ref,
                   divergence_x: float, glyco_density_x: float) -> float:
    """Back-compute the antigenicity v2 score of a protein that is absent
    from the surfaceome table, using the surfaceome distributions as the
    reference (mean-rank percentile, matching :func:`percentile_rank`).

    Used to rescue literature-validated antigens (CD46, CD59) that are
    missing from the porcine UniProt surfaceome owing to gene-annotation
    gaps.
    """
    def _pct(ref, x):
        # mean-rank percentile consistent with pandas rank(pct=True): an
        # element of the reference set itself scores (less + (eq+1)/2)/n
        ref = np.asarray(ref, dtype=float)
        less = np.sum(ref < x)
        equal = np.sum(ref == x)
        return (less + 0.5 * (equal + 1)) / len(ref)

    return 0.5 * _pct(divergence_ref, divergence_x) + \
        0.5 * _pct(glyco_density_ref, glyco_density_x)
