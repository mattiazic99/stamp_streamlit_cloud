"""Threshold sweep: switching gene counts as a function of tau in [0, 1].

For each tissue and each value of tau in a grid, we binarize the
normalized expression vector and count switching genes per Definition 1.
This reproduces Figure 6 and 7 of Kahveci et al. 2025 and provides the
robustness analysis requested by reviewers.

The sweep works directly on the normalized Parquet files produced by
01_normalize.py, so it does NOT need Cassandra or the raw TPM matrix.

Memory note
-----------
For each tissue we load one normalized Parquet (~2 MB), apply the sweep
vectorially across all tau values, and discard. Peak RAM is O(n_genes *
n_taus) per tissue, which is negligible.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from stamp.config import SWITCHING_BRACKETS
from stamp.switching import identify_switching_genes


# ---------------------------------------------------------------------------
# Default tau grid (matches Figure 6 of the paper)
# ---------------------------------------------------------------------------

DEFAULT_TAUS: list[float] = [round(t, 3) for t in np.arange(0.05, 1.0, 0.05).tolist()]


# ---------------------------------------------------------------------------
# Core sweep function
# ---------------------------------------------------------------------------

def sweep_threshold(
    normalized_df: pd.DataFrame,
    taus: list[float] | None = None,
) -> pd.DataFrame:
    """Count switching genes for each tau value, for one tissue.

    Parameters
    ----------
    normalized_df : DataFrame
        Rows = gene_id (index), columns = age brackets, values in [0, 1].
        Produced by stamp.normalize.normalize_tissue.
    taus : list of float, optional
        Grid of threshold values to test. Default: 0.05, 0.10, ..., 0.95.

    Returns
    -------
    DataFrame
        Rows = tau values. Columns:
            tau            : the threshold value tested
            n_switching    : total switching genes across all brackets
            n_30-39        : switching genes at bracket 30-39
            n_40-49        : switching genes at bracket 40-49
            n_50-59        : switching genes at bracket 50-59
            n_60-69        : switching genes at bracket 60-69
            n_70-79        : switching genes at bracket 70-79
    """
    if taus is None:
        taus = DEFAULT_TAUS

    rows = []
    for tau in taus:
        sets = identify_switching_genes(normalized_df, threshold=tau)
        row: dict[str, int | float] = {"tau": tau}
        total = 0
        for bracket in SWITCHING_BRACKETS:
            n = len(sets.get(bracket, []))
            row[f"n_{bracket}"] = n
            total += n
        row["n_switching"] = total
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Multi-tissue sweep
# ---------------------------------------------------------------------------

def sweep_all_tissues(
    tissues_normalized: dict[str, pd.DataFrame],
    taus: list[float] | None = None,
) -> pd.DataFrame:
    """Run sweep_threshold for multiple tissues and concatenate results.

    Parameters
    ----------
    tissues_normalized : dict
        Keys = tissue name, values = normalized DataFrame.

    Returns
    -------
    DataFrame
        Columns: tissue, tau, n_switching, n_30-39, ..., n_70-79.
        One row per (tissue, tau) combination.
    """
    all_rows = []
    for tissue, df in tissues_normalized.items():
        sweep = sweep_threshold(df, taus=taus)
        sweep.insert(0, "tissue", tissue)
        all_rows.append(sweep)
    return pd.concat(all_rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Stability score per tissue
# ---------------------------------------------------------------------------

def stability_score(sweep_df: pd.DataFrame, tau_range: tuple[float, float] = (0.3, 0.7)) -> pd.Series:
    """Coefficient of variation of n_switching across a tau range per tissue.

    A low CV means the tissue's switching count is stable (robust) to
    threshold choice. A high CV means it is very sensitive.

    Parameters
    ----------
    sweep_df : DataFrame
        Output of sweep_all_tissues.
    tau_range : tuple (min_tau, max_tau)
        Only consider taus within this range. Default (0.3, 0.7) covers
        the "central" region where switching is typically most stable.

    Returns
    -------
    Series
        Index = tissue name, values = CV (std/mean) of n_switching.
        Lower is more stable.
    """
    sub = sweep_df[
        (sweep_df["tau"] >= tau_range[0]) &
        (sweep_df["tau"] <= tau_range[1])
    ]
    grouped = sub.groupby("tissue")["n_switching"]
    mean = grouped.mean()
    std  = grouped.std()
    cv   = (std / mean).fillna(0).rename("stability_cv")
    return cv.sort_values()