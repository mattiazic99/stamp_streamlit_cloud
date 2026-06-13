"""Tests for per-gene min-max normalization and age-bracket averaging."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from stamp.normalize import normalize_tissue


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_metadata():
    """Three samples, one tissue, three different age brackets."""
    return pd.DataFrame({
        "sample_id": ["S1", "S2", "S3"],
        "tissue": ["Liver", "Liver", "Liver"],
        "age_bracket": ["20-29", "40-49", "60-69"],
    })


@pytest.fixture
def simple_tpm():
    """Three genes, three samples (S1, S2, S3)."""
    return pd.DataFrame(
        {
            "S1": [10.0, 0.0, 5.0],
            "S2": [20.0, 0.0, 5.0],
            "S3": [40.0, 0.0, 5.0],
        },
        index=["g_dynamic", "g_zero", "g_constant"],
    )


# ---------------------------------------------------------------------------
# Core behaviour
# ---------------------------------------------------------------------------

def test_min_max_canonical_formula(simple_tpm, simple_metadata):
    """The output must be in [0, 1] with min=0 and max=1 per gene."""
    out = normalize_tissue(simple_tpm, simple_metadata, "Liver")
    # g_dynamic: TPM [10, 20, 40] -> normalized [0, 1/3, 1]
    # After bracket averaging (one sample per bracket, no aggregation):
    assert out.loc["g_dynamic", "20-29"] == pytest.approx(0.0)
    assert out.loc["g_dynamic", "40-49"] == pytest.approx(1.0 / 3.0, rel=1e-5)
    assert out.loc["g_dynamic", "60-69"] == pytest.approx(1.0)


def test_epsilon_filter_drops_near_constant_genes(simple_tpm, simple_metadata):
    """Genes with range == 0 (e.g. all zeros, all 5s) must be dropped."""
    out = normalize_tissue(simple_tpm, simple_metadata, "Liver")
    assert "g_zero" not in out.index      # range is 0
    assert "g_constant" not in out.index  # range is 0
    assert "g_dynamic" in out.index


def test_epsilon_filter_custom_threshold():
    """Genes whose range is below custom epsilon must be dropped."""
    tpm = pd.DataFrame(
        {"S1": [1.0, 1.0], "S2": [1.005, 100.0]},
        index=["small_range", "big_range"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "tissue": ["Liver", "Liver"],
        "age_bracket": ["20-29", "30-39"],
    })
    # Default epsilon=0.01: small_range (range=0.005) dropped
    out = normalize_tissue(tpm, md, "Liver", epsilon=0.01)
    assert "small_range" not in out.index
    assert "big_range" in out.index
    # Stricter epsilon: even big_range survives but small_range is gone
    out2 = normalize_tissue(tpm, md, "Liver", epsilon=0.001)
    assert "small_range" in out2.index


def test_nan_handling_partial():
    """A NaN in one sample must not blow up the normalization."""
    tpm = pd.DataFrame(
        {"S1": [10.0, np.nan], "S2": [20.0, 30.0], "S3": [40.0, 60.0]},
        index=["g_full", "g_nan_in_S1"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2", "S3"],
        "tissue": ["Liver", "Liver", "Liver"],
        "age_bracket": ["20-29", "40-49", "60-69"],
    })
    out = normalize_tissue(tpm, md, "Liver")
    # g_full: standard
    assert "g_full" in out.index
    # g_nan_in_S1: min=30, max=60, range=30 -> S1 stays NaN, S2->0, S3->1
    assert "g_nan_in_S1" in out.index
    assert pd.isna(out.loc["g_nan_in_S1", "20-29"])
    assert out.loc["g_nan_in_S1", "40-49"] == pytest.approx(0.0)
    assert out.loc["g_nan_in_S1", "60-69"] == pytest.approx(1.0)


def test_nan_handling_all_nan_dropped():
    """A gene with all-NaN values must be dropped (range undefined)."""
    tpm = pd.DataFrame(
        {"S1": [np.nan, 10.0], "S2": [np.nan, 30.0]},
        index=["all_nan", "g_ok"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "tissue": ["Liver", "Liver"],
        "age_bracket": ["20-29", "40-49"],
    })
    out = normalize_tissue(tpm, md, "Liver")
    assert "all_nan" not in out.index
    assert "g_ok" in out.index


def test_age_bracket_averaging():
    """Multiple samples in the same bracket must be averaged."""
    tpm = pd.DataFrame(
        {"S1": [10.0], "S2": [20.0], "S3": [60.0]},
        index=["g1"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2", "S3"],
        "tissue": ["Liver", "Liver", "Liver"],
        # S1 and S2 in same bracket; S3 in another
        "age_bracket": ["20-29", "20-29", "40-49"],
    })
    out = normalize_tissue(tpm, md, "Liver")
    # min=10, max=60 -> normalized [0, 0.2, 1]
    # 20-29 = mean(0, 0.2) = 0.1
    # 40-49 = 1
    assert out.loc["g1", "20-29"] == pytest.approx(0.1, rel=1e-5)
    assert out.loc["g1", "40-49"] == pytest.approx(1.0, rel=1e-5)


def test_tissue_filter_correctness():
    """Samples from other tissues must be ignored entirely."""
    tpm = pd.DataFrame(
        {"S1": [10.0], "S2": [20.0], "S3": [99999.0]},
        index=["g1"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2", "S3"],
        "tissue": ["Liver", "Liver", "Brain"],  # S3 is in Brain!
        "age_bracket": ["20-29", "40-49", "20-29"],
    })
    out = normalize_tissue(tpm, md, "Liver")
    # Only S1 and S2 used: min=10, max=20, range=10
    # S1 -> 0, S2 -> 1
    assert out.loc["g1", "20-29"] == pytest.approx(0.0)
    assert out.loc["g1", "40-49"] == pytest.approx(1.0)


def test_columns_in_chronological_order():
    """Output columns must follow AGE_BRACKETS order, not arbitrary."""
    tpm = pd.DataFrame(
        {f"S{i}": [float(i)] for i in range(1, 7)},
        index=["g1"],
    )
    md = pd.DataFrame({
        "sample_id": [f"S{i}" for i in range(1, 7)],
        "tissue": ["Liver"] * 6,
        # Provide brackets in REVERSE order in metadata
        "age_bracket": ["70-79", "60-69", "50-59", "40-49", "30-39", "20-29"],
    })
    out = normalize_tissue(tpm, md, "Liver")
    # Output columns must be chronological: 20-29, 30-39, ..., 70-79
    assert list(out.columns) == ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]


def test_missing_brackets_excluded():
    """Brackets with no samples must NOT appear as empty/NaN columns."""
    tpm = pd.DataFrame(
        {"S1": [10.0], "S2": [40.0]},
        index=["g1"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "tissue": ["Liver", "Liver"],
        "age_bracket": ["30-39", "60-69"],  # only two brackets
    })
    out = normalize_tissue(tpm, md, "Liver")
    assert list(out.columns) == ["30-39", "60-69"]
    assert "20-29" not in out.columns
    assert "70-79" not in out.columns


def test_output_dtype_is_float32():
    """Output must be float32 for memory efficiency (Streamlit Cloud)."""
    tpm = pd.DataFrame({"S1": [1.0], "S2": [2.0]}, index=["g1"])
    md = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "tissue": ["Liver", "Liver"],
        "age_bracket": ["20-29", "30-39"],
    })
    out = normalize_tissue(tpm, md, "Liver")
    assert out.dtypes.iloc[0] == np.float32


# ---------------------------------------------------------------------------
# Validation / error cases
# ---------------------------------------------------------------------------

def test_unknown_tissue_raises(simple_tpm, simple_metadata):
    with pytest.raises(ValueError, match="No samples found"):
        normalize_tissue(simple_tpm, simple_metadata, "NonExistentTissue")


def test_metadata_missing_columns_raises(simple_tpm):
    bad_metadata = pd.DataFrame({"sample_id": ["S1"]})
    with pytest.raises(ValueError, match="missing required columns"):
        normalize_tissue(simple_tpm, bad_metadata, "Liver")


def test_all_genes_filtered_raises():
    """If every gene fails epsilon, error out instead of returning empty."""
    tpm = pd.DataFrame(
        {"S1": [1.0, 1.0], "S2": [1.0, 1.0]},  # all constant
        index=["g1", "g2"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "tissue": ["Liver", "Liver"],
        "age_bracket": ["20-29", "30-39"],
    })
    with pytest.raises(ValueError, match="All genes filtered"):
        normalize_tissue(tpm, md, "Liver")


def test_does_not_mutate_input(simple_tpm, simple_metadata):
    """The function must not alter the caller's DataFrames."""
    tpm_before = simple_tpm.copy()
    md_before = simple_metadata.copy()
    _ = normalize_tissue(simple_tpm, simple_metadata, "Liver")
    pd.testing.assert_frame_equal(simple_tpm, tpm_before)
    pd.testing.assert_frame_equal(simple_metadata, md_before)


def test_metadata_extra_columns_ignored(simple_tpm):
    """Extra columns in metadata must not break the function."""
    md = pd.DataFrame({
        "sample_id": ["S1", "S2", "S3"],
        "tissue": ["Liver", "Liver", "Liver"],
        "age_bracket": ["20-29", "40-49", "60-69"],
        "sex": ["M", "F", "M"],
        "rin": [7.5, 8.0, 6.9],
    })
    out = normalize_tissue(simple_tpm, md, "Liver")
    assert not out.empty


def test_sample_with_unknown_bracket_dropped():
    """A sample with NaN/unknown age_bracket must be excluded BEFORE min/max.

    A sample without a valid age_bracket never enters any bracket average, so
    it must not influence the min-max scale either. Here S_unknown (TPM 9999)
    must be ignored entirely: the scale is defined only by S1 and S2.
    """
    tpm = pd.DataFrame(
        {"S1": [10.0], "S2": [20.0], "S_unknown": [9999.0]},
        index=["g1"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2", "S_unknown"],
        "tissue": ["Liver", "Liver", "Liver"],
        "age_bracket": ["20-29", "40-49", None],  # S_unknown has no bracket
    })
    out = normalize_tissue(tpm, md, "Liver")
    # Output should only have 20-29 and 40-49
    assert set(out.columns) == {"20-29", "40-49"}
    # S_unknown (9999) is excluded from the scale: min=10, max=20, range=10.
    # So S1 -> 0.0 and S2 -> 1.0. (Under the old behaviour, max would have been
    # 9999 and S2 would have been ~0.001 — this assertion locks in the fix.)
    assert out.loc["g1", "20-29"] == pytest.approx(0.0)
    assert out.loc["g1", "40-49"] == pytest.approx(1.0)


def test_sample_with_out_of_range_bracket_excluded_from_scale():
    """A sample whose age_bracket is not in AGE_BRACKETS must not affect min/max."""
    tpm = pd.DataFrame(
        {"S1": [10.0], "S2": [20.0], "S_old": [9999.0]},
        index=["g1"],
    )
    md = pd.DataFrame({
        "sample_id": ["S1", "S2", "S_old"],
        "tissue": ["Liver", "Liver", "Liver"],
        # "80-89" is a real GTEx-style label but NOT in AGE_BRACKETS.
        "age_bracket": ["20-29", "40-49", "80-89"],
    })
    out = normalize_tissue(tpm, md, "Liver")
    assert set(out.columns) == {"20-29", "40-49"}
    assert "80-89" not in out.columns
    # Scale defined by S1/S2 only: S2 -> 1.0, not ~0.001.
    assert out.loc["g1", "40-49"] == pytest.approx(1.0)