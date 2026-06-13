"""Tissue completeness helpers for the "complete age-bin coverage" mode.

A tissue has *complete age-bin coverage* if it has at least one valid sample
in ALL six GTEx age brackets (`stamp.config.AGE_BRACKETS`). Tissues missing
one or more brackets cannot be evaluated for switching across the full
lifespan without treating a missing bracket as "no data" — so the main
analysis can optionally restrict to complete tissues, keeping all-tissue
results as a sensitivity analysis.

These helpers compute completeness *dynamically* from the metadata, so they
stay correct if the underlying GTEx dump changes. The hard-coded reference
sets in `stamp.config` are for documentation/tests only.
"""
from __future__ import annotations

import pandas as pd

from stamp.config import AGE_BRACKETS


def tissues_with_complete_age_bins(metadata: pd.DataFrame) -> list[str]:
    """Return the sorted list of tissues with a valid sample in every bracket.

    Parameters
    ----------
    metadata : DataFrame
        Must contain columns ``tissue`` and ``age_bracket``.

    Returns
    -------
    list[str]
        Tissues (sorted) that have >= 1 sample in EACH of the six
        ``AGE_BRACKETS``. Samples whose ``age_bracket`` is NaN or outside
        ``AGE_BRACKETS`` are ignored when assessing completeness.

    Raises
    ------
    ValueError
        If the required columns are missing.
    """
    required_cols = {"tissue", "age_bracket"}
    missing = required_cols - set(metadata.columns)
    if missing:
        raise ValueError(f"metadata is missing required columns: {sorted(missing)}")

    required_brackets = set(AGE_BRACKETS)

    # Keep only rows with a valid (in-range, non-null) bracket.
    valid = metadata[metadata["age_bracket"].isin(AGE_BRACKETS)]

    present_by_tissue = valid.groupby("tissue")["age_bracket"].agg(
        lambda s: set(s.unique())
    )

    return sorted(
        tissue
        for tissue, present in present_by_tissue.items()
        if required_brackets.issubset(present)
    )


def incomplete_tissues(metadata: pd.DataFrame) -> list[str]:
    """Return the sorted list of tissues missing at least one age bracket."""
    all_tissues = set(metadata["tissue"].dropna().unique())
    complete = set(tissues_with_complete_age_bins(metadata))
    return sorted(all_tissues - complete)


def common_complete_tissues(
    metadata_a: pd.DataFrame,
    metadata_b: pd.DataFrame,
) -> list[str]:
    """Tissues with complete coverage in BOTH versions (for v8-vs-v10).

    This is the like-for-like set used when comparing two GTEx releases on
    identical tissues: a tissue is kept only if it is complete in both
    ``metadata_a`` and ``metadata_b`` (equivalently, excluding the union of
    each version's incomplete tissues).
    """
    a = set(tissues_with_complete_age_bins(metadata_a))
    b = set(tissues_with_complete_age_bins(metadata_b))
    return sorted(a & b)
