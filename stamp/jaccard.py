"""Jaccard similarity between tissues based on switching gene sets.

Two metrics, both defined in Kahveci et al. 2025:

  J_age(A, B):  mean over the 5 switching age brackets of J(S_A,i, S_B,i)
                where S_X,i is the set of genes switching in tissue X at
                bracket i. Captures "do these two tissues switch the
                same genes at the same ages?".

  J_life(A, B): J(union_i S_A,i, union_i S_B,i). Captures "do these two
                tissues share a similar overall pool of aging-related
                switching genes?".

Memory considerations
---------------------
Both matrices are 54x54 float32 (or smaller). They are tiny; we save
them as plain CSV with float values rounded to 6 decimals.

Edge cases
----------
- Jaccard of two empty sets is conventionally defined as 0 here
  (a tissue with no switching genes at a bracket cannot share with
  any other). This is consistent with how the paper visualizes
  these matrices (empty rows = pure black in heatmaps).
- Tissues with missing brackets (Cervix Ectocervix, Cervix
  Endocervix, Fallopian Tube, Kidney Medulla) contribute empty
  sets for the missing brackets. J_age averages over the brackets
  that BOTH tissues have populated; if a bracket is missing in
  either, it is skipped (not counted as 0). This avoids
  artificially deflating the score for incomplete tissues.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from stamp.config import SWITCHING_BRACKETS


# ---------------------------------------------------------------------------
# Core: pairwise Jaccard on sets
# ---------------------------------------------------------------------------

def jaccard(set_a: set[str], set_b: set[str]) -> float:
    """Standard Jaccard index. Returns 0.0 for two empty sets."""
    if not set_a and not set_b:
        return 0.0
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union > 0 else 0.0


# ---------------------------------------------------------------------------
# Per-pair similarity metrics
# ---------------------------------------------------------------------------

def jaccard_age_pair(
    sets_a: dict[str, list[str]],
    sets_b: dict[str, list[str]],
) -> float:
    """Mean Jaccard across the switching brackets present in both tissues.

    Parameters
    ----------
    sets_a, sets_b : dict
        Keys are bracket labels (subset of SWITCHING_BRACKETS), values are
        gene_id lists. Missing keys are treated as "bracket not measured".

    Returns
    -------
    float in [0, 1]
        Mean Jaccard over the brackets common to both tissues. If no
        bracket is common, returns 0.0.
    """
    common = [b for b in SWITCHING_BRACKETS if b in sets_a and b in sets_b]
    if not common:
        return 0.0
    values = [jaccard(set(sets_a[b]), set(sets_b[b])) for b in common]
    return float(np.mean(values))


def jaccard_life_pair(
    sets_a: dict[str, list[str]],
    sets_b: dict[str, list[str]],
) -> float:
    """Jaccard on the union of switching genes across all brackets."""
    union_a: set[str] = set()
    for genes in sets_a.values():
        union_a.update(genes)
    union_b: set[str] = set()
    for genes in sets_b.values():
        union_b.update(genes)
    return jaccard(union_a, union_b)


# ---------------------------------------------------------------------------
# Full similarity matrices
# ---------------------------------------------------------------------------

def tissue_similarity_age(
    sets_by_tissue: dict[str, dict[str, list[str]]],
) -> pd.DataFrame:
    """Compute the J_age similarity matrix.

    Parameters
    ----------
    sets_by_tissue : dict
        Keys = tissue names. Values = dict mapping bracket label to
        list of gene_ids switching there.

    Returns
    -------
    DataFrame
        Symmetric matrix indexed by tissue name on both axes, dtype
        float32. Diagonal is 1.0 (a tissue is identical to itself).
    """
    return _build_matrix(sets_by_tissue, jaccard_age_pair)


def tissue_similarity_life(
    sets_by_tissue: dict[str, dict[str, list[str]]],
) -> pd.DataFrame:
    """Compute the J_life similarity matrix.

    See `tissue_similarity_age` for parameters / returns.
    """
    return _build_matrix(sets_by_tissue, jaccard_life_pair)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_matrix(
    sets_by_tissue: dict[str, dict[str, list[str]]],
    pair_metric,
) -> pd.DataFrame:
    """Build a symmetric similarity matrix using `pair_metric`."""
    tissues = sorted(sets_by_tissue.keys())
    n = len(tissues)
    mat = np.zeros((n, n), dtype=np.float32)

    # Pre-convert lists to sets once (avoids re-doing it n*(n-1)/2 times)
    pre = {
        t: {b: set(genes) for b, genes in sets.items()}
        for t, sets in sets_by_tissue.items()
    }

    for i, t_i in enumerate(tissues):
        mat[i, i] = 1.0
        for j in range(i + 1, n):
            t_j = tissues[j]
            v = pair_metric(pre[t_i], pre[t_j])
            mat[i, j] = v
            mat[j, i] = v  # symmetric

    return pd.DataFrame(mat, index=tissues, columns=tissues, dtype=np.float32)