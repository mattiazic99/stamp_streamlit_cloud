"""Find tissues with samples in all six GTEx age brackets.

Complete mode excludes tissues missing a bracket; all-tissues mode retains them
for sensitivity analyses. Coverage is computed from metadata. The reference sets
in stamp.config are for documentation and tests, not for filtering.
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
        Tissues (sorted) that have >= 1 sample in each of the six
        ``AGE_BRACKETS``. Samples whose ``age_bracket`` is NaN or outside
        ``AGE_BRACKETS`` are ignored when assessing completeness.
    
    Raises
    ------
    ValueError
        If the required columns are missing."""
    required_cols = {"tissue", "age_bracket"}
    missing = required_cols - set(metadata.columns)
    if missing:
        raise ValueError(f"metadata is missing required columns: {sorted(missing)}")

    required_brackets = set(AGE_BRACKETS)

    # Missing and out-of-range labels cannot satisfy bracket coverage.
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
    """Return sorted tissues with full age-bin coverage in both metadata tables.
    
    This keeps the tissue set identical when comparing GTEx releases.
    """
    a = set(tissues_with_complete_age_bins(metadata_a))
    b = set(tissues_with_complete_age_bins(metadata_b))
    return sorted(a & b)
