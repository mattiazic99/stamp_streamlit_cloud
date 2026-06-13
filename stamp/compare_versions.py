"""Numerical comparison between GTEx v8 and v10 pipeline outputs.

Used in the paper's Validation section to demonstrate that STAMP results
are qualitatively consistent across GTEx releases.
"""
import pandas as pd


def compare_normalized(v8_df: pd.DataFrame, v10_df: pd.DataFrame) -> pd.DataFrame:
    """Diff normalized matrices: gene-level statistics on overlap and deviation."""
    raise NotImplementedError("TODO")


def compare_switching_sets(
    v8_sets: dict[str, dict[str, list[str]]],
    v10_sets: dict[str, dict[str, list[str]]],
) -> pd.DataFrame:
    """Per-tissue Jaccard of switching gene sets between v8 and v10."""
    raise NotImplementedError("TODO")


def compare_jaccard_matrices(
    v8_mat: pd.DataFrame, v10_mat: pd.DataFrame
) -> pd.DataFrame:
    """Element-wise comparison of Jaccard tissue-similarity matrices."""
    raise NotImplementedError("TODO")
