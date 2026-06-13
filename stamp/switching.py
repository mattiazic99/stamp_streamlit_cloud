"""Switching gene detection (Definition 1 in Kahveci et al. 2025).

A gene is "switching" in a tissue if its binary expression vector across
age brackets shows exactly one permanent transition: e.g. [0,0,1,1,1,1]
or [1,1,0,0,0,0]. Vectors with multiple transitions or constant values
are NOT switching.

Formal definition (Definition 1 of the paper)
---------------------------------------------
Given a time-ordered sequence of n binary states [t_1, t_2, ..., t_n],
the gene is switching if there exists an integer k such that:
    1. 1 < k <= n
    2. for all i, j < k: t_i == t_j         (the prefix is constant)
    3. for all i, j >= k: t_i == t_j        (the suffix is constant)
    4. there is no (i, j) with i < k and j >= k and t_i == t_j
       (the prefix value and suffix value differ)

This implementation is more rigorous than the original R code, which
in `auto_sets.R` accepts vectors with multiple transitions (only the
first is detected). Divergence is documented in docs/divergences_from_R.md.

Bracket convention
------------------
Per `stamp.config.AGE_BRACKETS`, vectors have length 6:
    [20-29, 30-39, 40-49, 50-59, 60-69, 70-79]

The switching bracket is the bracket where the new state is first
observed. By Definition 1 condition (1), k must be > 1, so the switch
cannot be assigned to bracket index 0 (i.e. [20-29]). Valid switching
bracket indices are therefore {1, 2, 3, 4, 5}, corresponding to the
five `SWITCHING_BRACKETS` defined in `stamp.config`.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from stamp.config import AGE_BRACKETS, SWITCHING_BRACKETS


# ---------------------------------------------------------------------------
# Core predicates: pure functions on a single binary vector
# ---------------------------------------------------------------------------

def is_switching_gene(binary_vector: Sequence[int]) -> bool:
    """Check whether a binary vector qualifies as switching per Definition 1.

    Parameters
    ----------
    binary_vector : sequence of int
        Values must be 0 or 1. Length must be >= 2.

    Returns
    -------
    bool
        True iff there exists exactly one k > 0 such that all elements
        before k are equal, all elements from k onward are equal, and
        the two values differ.

    Raises
    ------
    ValueError
        If the vector is shorter than 2 elements, or contains values
        other than 0 and 1.

    Examples
    --------
    >>> is_switching_gene([0, 0, 1, 1, 1, 1])    # switch at index 2
    True
    >>> is_switching_gene([1, 1, 1, 0, 0, 0])    # switch at index 3
    True
    >>> is_switching_gene([0, 0, 0, 0, 0, 0])    # constant
    False
    >>> is_switching_gene([0, 1, 0, 1, 0, 1])    # multiple transitions
    False
    """
    v = _validate_binary_vector(binary_vector)

    # A vector is switching iff there is exactly one position where v[i] != v[i-1].
    # This is the most direct translation of Definition 1: a single permanent
    # transition is equivalent to a single index where consecutive elements differ.
    transitions = np.diff(v)
    n_transitions = int(np.count_nonzero(transitions))
    return n_transitions == 1


def get_switching_bracket(binary_vector: Sequence[int]) -> Optional[int]:
    """Return the 0-based bracket index where the switch occurs.

    By Definition 1, this is the smallest k > 0 such that v[k] != v[k-1],
    provided the vector is switching. Otherwise None is returned.

    Parameters
    ----------
    binary_vector : sequence of int
        Values must be 0 or 1.

    Returns
    -------
    int or None
        Index k in [1, len(vector) - 1] for switching vectors, else None.

    Examples
    --------
    >>> get_switching_bracket([0, 0, 1, 1, 1, 1])    # switches at index 2
    2
    >>> get_switching_bracket([1, 1, 1, 0, 0, 0])    # switches at index 3
    3
    >>> get_switching_bracket([0, 0, 0, 0, 0, 0])    # not switching
    >>> get_switching_bracket([0, 1, 0, 1, 0, 1])    # multiple transitions
    """
    if not is_switching_gene(binary_vector):
        return None
    v = np.asarray(binary_vector, dtype=np.int8)
    # First (and only) position where consecutive elements differ.
    # np.diff returns differences of length n-1; the +1 maps back to v's index.
    return int(np.argmax(np.abs(np.diff(v)) != 0)) + 1


# ---------------------------------------------------------------------------
# Pipeline-level function: from a normalized matrix to per-bracket gene sets
# ---------------------------------------------------------------------------

def identify_switching_genes(
    normalized_df: pd.DataFrame,
    threshold: float,
) -> dict[str, list[str]]:
    """Identify switching genes per age bracket for one tissue.

    Each gene's row is binarised with the given threshold, then tested
    against Definition 1. Switching genes are bucketed by the bracket
    at which the transition occurs.

    Parameters
    ----------
    normalized_df : DataFrame
        Rows = gene_id (index), columns = age brackets in chronological
        order. Values must be in [0, 1] (or NaN). The columns must be a
        subset of `AGE_BRACKETS` and contain at least 2 brackets.
    threshold : float
        Binarisation threshold tau in [0, 1]. A value v is mapped to 1
        iff v >= threshold. NaN values propagate (a row containing any
        NaN is excluded from switching identification).

    Returns
    -------
    dict[str, list[str]]
        Keys are the SWITCHING_BRACKETS labels present in the input
        columns (e.g. "30-39", "40-49", ...). Values are the lists of
        gene_ids that switch at that bracket. Every key in
        SWITCHING_BRACKETS that is also a column of normalized_df is
        present in the dict, even if the value is an empty list. Genes
        that switch but whose bracket falls outside SWITCHING_BRACKETS
        (i.e. index 0, [20-29]) are by construction impossible because
        Definition 1 requires k > 0; this function therefore never
        produces a "20-29" key.

    Notes
    -----
    Genes with any NaN in their normalised vector are silently skipped:
    the binarisation of NaN is undefined, and a partial vector cannot
    be evaluated against Definition 1 without making implicit assumptions.
    """
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(f"threshold must be in [0, 1], got {threshold}")

    columns = list(normalized_df.columns)
    if len(columns) < 2:
        raise ValueError(
            f"normalized_df must have at least 2 columns, got {len(columns)}"
        )

    # Verify that columns are valid age brackets and in the correct order.
    valid_columns_in_order = [b for b in AGE_BRACKETS if b in columns]
    if list(columns) != valid_columns_in_order:
        raise ValueError(
            "normalized_df columns must be a subset of AGE_BRACKETS in "
            f"chronological order. Got {columns}, expected an ordered "
            f"subset of {AGE_BRACKETS}."
        )

    # Pre-allocate result with empty lists for every reachable bracket.
    # Brackets reachable in this call are those columns that are in
    # SWITCHING_BRACKETS. Brackets not present as columns are not keys.
    result: dict[str, list[str]] = {
        b: [] for b in columns if b in SWITCHING_BRACKETS
    }

    # Binarise: NaN propagates as NaN through the comparison; we cast to
    # float so that NaN survives. Values >= threshold become 1.0; values
    # below become 0.0; NaN stays NaN.
    values = normalized_df.to_numpy(dtype=np.float64, copy=True)
    binarised = np.where(np.isnan(values), np.nan, (values >= threshold).astype(float))

    # Identify rows free of NaN (the only ones we can evaluate).
    rows_with_nan = np.isnan(binarised).any(axis=1)

    # For NaN-free rows, count transitions vectorially.
    binary_int = np.where(rows_with_nan[:, None], 0, binarised).astype(np.int8)
    diffs = np.diff(binary_int, axis=1)
    n_transitions = np.count_nonzero(diffs, axis=1)

    # A row is switching iff it has no NaN AND exactly one transition.
    is_switching = (~rows_with_nan) & (n_transitions == 1)

    # For switching rows, find the bracket index where the transition occurs.
    # argmax on the boolean mask of "row has nonzero diff at position j"
    # returns the FIRST nonzero column; this is correct because for
    # switching rows there is exactly one nonzero diff.
    diff_nonzero = diffs != 0
    switch_positions = np.argmax(diff_nonzero, axis=1) + 1  # shift back to v index

    gene_ids = normalized_df.index.to_numpy()
    for row_idx in np.flatnonzero(is_switching):
        bracket_idx = int(switch_positions[row_idx])
        bracket_label = columns[bracket_idx]
        # bracket_label cannot be AGE_BRACKETS[0] because diff index 0
        # corresponds to v index 1 at minimum. So it's always in result.
        result[bracket_label].append(str(gene_ids[row_idx]))

    return result


def identify_switching_with_direction(
    normalized_df: pd.DataFrame,
    threshold: float,
) -> dict[str, list[tuple[str, str]]]:
    """Like ``identify_switching_genes`` but also records the switch direction.

    A switching gene's single transition is either low->high ("up") or
    high->low ("down"). The set of switching genes and their assigned brackets
    are IDENTICAL to ``identify_switching_genes`` (same binarisation and
    one-transition rule); this function only adds the direction label.

    Parameters
    ----------
    normalized_df : DataFrame
        Rows = gene_id, columns = age brackets in chronological order,
        values in [0, 1] (or NaN). Same contract as ``identify_switching_genes``.
    threshold : float
        Binarisation threshold tau; v is mapped to 1 iff v >= threshold.

    Returns
    -------
    dict[str, list[tuple[str, str]]]
        Keys are SWITCHING_BRACKETS present in the columns. Values are lists of
        ``(gene_id, direction)`` where direction is ``"up"`` (low->high) or
        ``"down"`` (high->low).
    """
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(f"threshold must be in [0, 1], got {threshold}")

    columns = list(normalized_df.columns)
    if len(columns) < 2:
        raise ValueError(
            f"normalized_df must have at least 2 columns, got {len(columns)}"
        )

    valid_columns_in_order = [b for b in AGE_BRACKETS if b in columns]
    if list(columns) != valid_columns_in_order:
        raise ValueError(
            "normalized_df columns must be a subset of AGE_BRACKETS in "
            f"chronological order. Got {columns}."
        )

    result: dict[str, list[tuple[str, str]]] = {
        b: [] for b in columns if b in SWITCHING_BRACKETS
    }

    values = normalized_df.to_numpy(dtype=np.float64, copy=True)
    binarised = np.where(np.isnan(values), np.nan, (values >= threshold).astype(float))
    rows_with_nan = np.isnan(binarised).any(axis=1)

    binary_int = np.where(rows_with_nan[:, None], 0, binarised).astype(np.int8)
    diffs = np.diff(binary_int, axis=1)
    n_transitions = np.count_nonzero(diffs, axis=1)
    is_switching = (~rows_with_nan) & (n_transitions == 1)

    diff_nonzero = diffs != 0
    switch_positions = np.argmax(diff_nonzero, axis=1)  # index into diffs

    gene_ids = normalized_df.index.to_numpy()
    for row_idx in np.flatnonzero(is_switching):
        pos = int(switch_positions[row_idx])           # diff index
        bracket_label = columns[pos + 1]               # bracket after transition
        direction = "up" if diffs[row_idx, pos] > 0 else "down"
        result[bracket_label].append((str(gene_ids[row_idx]), direction))

    return result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_binary_vector(binary_vector: Sequence[int]) -> np.ndarray:
    """Validate and convert input to a 1-D numpy int8 array."""
    v = np.asarray(binary_vector)
    if v.ndim != 1:
        raise ValueError(f"binary_vector must be 1-D, got shape {v.shape}")
    if v.size < 2:
        raise ValueError(
            f"binary_vector must have length >= 2, got length {v.size}"
        )
    # Values must be exactly 0 or 1.
    unique_values = set(np.unique(v).tolist())
    if not unique_values.issubset({0, 1}):
        raise ValueError(
            f"binary_vector must contain only 0 and 1, got {sorted(unique_values)}"
        )
    return v.astype(np.int8)