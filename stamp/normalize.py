"""Per-gene min-max normalization and age-bracket averaging.

Implements the canonical normalization (x - min) / (max - min), with an
epsilon filter that drops genes whose dynamic range is below threshold
across the tissue's samples.

Divergence from the R reference implementation
----------------------------------------------
The original R code in `auto_normalized.R` divides by (max - min) WITHOUT
subtracting min: a non-standard formula that does not produce values in
[0, 1]. Our Python implementation follows the paper specification
(page 5: "per-gene min-max normalisation"), which is the canonical
formula. See docs/divergences_from_R.md.

Memory considerations
---------------------
All numerical outputs are float32 (sufficient precision for TPM data
and half the memory of float64). The normalized matrix per tissue is
designed to fit comfortably in <100 MB even for the largest GTEx tissues.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from stamp.config import AGE_BRACKETS, EPSILON


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def normalize_tissue(
    tpm_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    tissue: str,
    epsilon: float = EPSILON,
) -> pd.DataFrame:
    """Normalize TPMs for one tissue and average per age bracket.

    Pipeline (per the paper, page 5):
        1. Select samples belonging to the requested tissue.
        2. Drop genes whose (max - min) across these samples is below
           epsilon (constant or near-constant: not informative).
        3. Apply per-gene min-max normalization: (x - min) / (max - min).
        4. Group samples by age bracket; compute the per-gene mean.
        5. Reindex columns to match `AGE_BRACKETS` order; drop empty
           brackets (no samples for this tissue in that age range).

    Parameters
    ----------
    tpm_matrix : DataFrame
        Rows = gene_id (index), columns = sample_id, values = TPM.
        Float dtype; NaN allowed (treated as missing).
    metadata : DataFrame
        Sample metadata with columns 'sample_id', 'tissue', 'age_bracket'
        at minimum. 'sample_id' must match column names of tpm_matrix.
    tissue : str
        Tissue name to filter on. Matched exactly against metadata['tissue'].
    epsilon : float, default 0.01
        Minimum (max - min) for a gene to be retained.

    Returns
    -------
    DataFrame
        Rows = gene_id (subset of input that passed the epsilon filter),
        columns = age brackets present for this tissue (in AGE_BRACKETS
        order), dtype float32, values in [0, 1].

    Raises
    ------
    ValueError
        If tissue has no samples, or required metadata columns are missing.

    Notes
    -----
    A gene with all-NaN values for this tissue is dropped (range is
    undefined). A gene whose values are all equal is also dropped (range
    is zero, fails the epsilon filter).
    """
    _validate_metadata_schema(metadata)

    # Step 1: select the samples for this tissue that will actually contribute
    # to an age-bracket average. A sample is used only if it
    #   (a) belongs to the requested tissue,
    #   (b) has a matching column in the TPM matrix, and
    #   (c) has a valid age_bracket present in AGE_BRACKETS.
    # Condition (c) is essential: a sample without a valid bracket never enters
    # any bracket mean, so it must NOT influence gene_min / gene_max /
    # gene_range, the epsilon filter, or the min-max scale. Filtering it here
    # (rather than after normalization) keeps the [0, 1] scale defined only by
    # the samples that carry signal.
    tissue_meta = metadata.loc[metadata["tissue"] == tissue]
    if tissue_meta.empty:
        raise ValueError(f"No samples found for tissue '{tissue}'")

    # isin() treats NaN and out-of-range labels (e.g. "80-89") as invalid.
    valid_bracket = tissue_meta["age_bracket"].isin(AGE_BRACKETS)
    tissue_samples = [
        s
        for s in tissue_meta.loc[valid_bracket, "sample_id"].tolist()
        if s in tpm_matrix.columns
    ]
    if not tissue_samples:
        raise ValueError(
            f"Tissue '{tissue}' has no samples with both a matching TPM column "
            f"and a valid age_bracket (one of {AGE_BRACKETS})."
        )

    # Sub-matrix for this tissue: copy to avoid mutating caller's data
    sub = tpm_matrix[tissue_samples].astype(np.float32, copy=True)

    # Step 2: epsilon filter
    # NaN-safe min/max; if a gene is all NaN, range will be NaN and filter drops it.
    gene_min = sub.min(axis=1, skipna=True)
    gene_max = sub.max(axis=1, skipna=True)
    gene_range = gene_max - gene_min
    keep = (gene_range > epsilon).fillna(False)

    if not keep.any():
        raise ValueError(
            f"All genes filtered out for tissue '{tissue}' "
            f"(epsilon={epsilon}). Check input data."
        )

    sub = sub.loc[keep]
    gene_min = gene_min.loc[keep]
    gene_range = gene_range.loc[keep]

    # Step 3: per-gene min-max normalization
    # Broadcasting subtracts min row-wise, then divides by range.
    normalized = sub.sub(gene_min, axis=0).div(gene_range, axis=0)

    # Step 4: average per age bracket
    # Build a sample -> bracket map for the tissue's samples
    sample_to_bracket = (
        metadata.loc[metadata["sample_id"].isin(tissue_samples)]
        .set_index("sample_id")["age_bracket"]
    )

    # Pandas groupby on columns: transpose, group, mean, transpose back
    bracket_per_column = normalized.columns.map(sample_to_bracket)

    # Some samples may have a missing/unknown age_bracket: drop those columns.
    valid_cols_mask = bracket_per_column.notna()
    if not valid_cols_mask.any():
        raise ValueError(
            f"No samples with valid age_bracket for tissue '{tissue}'"
        )
    if not valid_cols_mask.all():
        normalized = normalized.loc[:, valid_cols_mask]
        bracket_per_column = bracket_per_column[valid_cols_mask]

    # Group columns by bracket and mean (skipna=True by default)
    averaged = normalized.T.groupby(bracket_per_column, observed=True).mean().T

    # Step 5: reindex columns to AGE_BRACKETS chronological order
    # Keep only brackets actually present (no all-NaN columns inserted)
    ordered = [b for b in AGE_BRACKETS if b in averaged.columns]
    averaged = averaged[ordered]

    # Final dtype: float32 to halve memory vs float64
    averaged = averaged.astype(np.float32)

    # Drop genes that ended up all-NaN after averaging (rare but possible if
    # every kept sample of a gene was NaN in every bracket)
    averaged = averaged.dropna(how="all")

    return averaged


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_metadata_schema(metadata: pd.DataFrame) -> None:
    """Check that metadata has the columns we need."""
    required = {"sample_id", "tissue", "age_bracket"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(
            f"metadata is missing required columns: {sorted(missing)}"
        )