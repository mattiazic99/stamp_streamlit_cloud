"""Exhaustive tests for Definition 1 of switching genes.

Reference: Kahveci et al. 2025, Definition 1 (page 6).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from stamp.switching import (
    get_switching_bracket,
    identify_switching_genes,
    identify_switching_with_direction,
    is_switching_gene,
)


# ---------------------------------------------------------------------------
# identify_switching_with_direction: coherent with identify_switching_genes
# ---------------------------------------------------------------------------

def _toy_norm():
    # g_up: low->high at 40-49; g_down: high->low at 60-69; g_const: no switch;
    # g_osc: oscillating (excluded)
    return pd.DataFrame(
        {
            "20-29": [0.0, 1.0, 0.9, 0.0],
            "30-39": [0.0, 1.0, 0.9, 1.0],
            "40-49": [1.0, 1.0, 0.9, 0.0],
            "50-59": [1.0, 1.0, 0.9, 1.0],
            "60-69": [1.0, 0.0, 0.9, 0.0],
            "70-79": [1.0, 0.0, 0.9, 1.0],
        },
        index=["g_up", "g_down", "g_const", "g_osc"],
    )


def test_direction_matches_plain_switching_set():
    df = _toy_norm()
    plain = identify_switching_genes(df, 0.5)
    withdir = identify_switching_with_direction(df, 0.5)
    # Same gene assigned to same bracket in both functions.
    for bracket in plain:
        assert set(plain[bracket]) == {g for g, _ in withdir[bracket]}


def test_direction_labels_up_and_down():
    df = _toy_norm()
    withdir = identify_switching_with_direction(df, 0.5)
    assert ("g_up", "up") in withdir["40-49"]
    assert ("g_down", "down") in withdir["60-69"]
    # constant and oscillating genes never appear
    flat = [g for b in withdir for g, _ in withdir[b]]
    assert "g_const" not in flat
    assert "g_osc" not in flat


# ---------------------------------------------------------------------------
# is_switching_gene: low-to-high vectors at every valid k
# ---------------------------------------------------------------------------

def test_low_to_high_at_bracket_2():
    # k=1 in 0-indexed; switch happens between bracket 0 and bracket 1
    assert is_switching_gene([0, 1, 1, 1, 1, 1]) is True


def test_low_to_high_at_bracket_3():
    assert is_switching_gene([0, 0, 1, 1, 1, 1]) is True  # ATP11C example


def test_low_to_high_at_bracket_4():
    assert is_switching_gene([0, 0, 0, 1, 1, 1]) is True


def test_low_to_high_at_bracket_5():
    assert is_switching_gene([0, 0, 0, 0, 1, 1]) is True


def test_low_to_high_at_bracket_6():
    assert is_switching_gene([0, 0, 0, 0, 0, 1]) is True


# ---------------------------------------------------------------------------
# is_switching_gene: high-to-low vectors at every valid k
# ---------------------------------------------------------------------------

def test_high_to_low_at_bracket_2():
    assert is_switching_gene([1, 0, 0, 0, 0, 0]) is True


def test_high_to_low_at_bracket_3():
    assert is_switching_gene([1, 1, 0, 0, 0, 0]) is True


def test_high_to_low_at_bracket_4():
    assert is_switching_gene([1, 1, 1, 0, 0, 0]) is True


def test_high_to_low_at_bracket_5():
    assert is_switching_gene([1, 1, 1, 1, 0, 0]) is True


def test_high_to_low_at_bracket_6():
    assert is_switching_gene([1, 1, 1, 1, 1, 0]) is True


# ---------------------------------------------------------------------------
# is_switching_gene: NOT switching cases
# ---------------------------------------------------------------------------

def test_constant_zero_is_not_switching():
    assert is_switching_gene([0, 0, 0, 0, 0, 0]) is False


def test_constant_one_is_not_switching():
    assert is_switching_gene([1, 1, 1, 1, 1, 1]) is False


def test_oscillating_010101_is_not_switching():
    assert is_switching_gene([0, 1, 0, 1, 0, 1]) is False


def test_oscillating_101010_is_not_switching():
    assert is_switching_gene([1, 0, 1, 0, 1, 0]) is False


def test_two_transitions_is_not_switching():
    # Switches up then back down: not Definition 1
    assert is_switching_gene([0, 0, 1, 1, 0, 0]) is False


def test_three_transitions_is_not_switching():
    assert is_switching_gene([0, 1, 0, 1, 1, 1]) is False


def test_short_vector_edge_case():
    # Length-2 vector is the minimum valid input
    assert is_switching_gene([0, 1]) is True
    assert is_switching_gene([1, 0]) is True
    assert is_switching_gene([0, 0]) is False
    assert is_switching_gene([1, 1]) is False


# ---------------------------------------------------------------------------
# is_switching_gene: input validation
# ---------------------------------------------------------------------------

def test_rejects_too_short():
    with pytest.raises(ValueError, match="length >= 2"):
        is_switching_gene([1])


def test_rejects_non_binary_values():
    with pytest.raises(ValueError, match="only 0 and 1"):
        is_switching_gene([0, 0, 2, 1, 1, 1])


def test_rejects_non_1d_array():
    with pytest.raises(ValueError, match="1-D"):
        is_switching_gene(np.array([[0, 1], [1, 0]]))


# ---------------------------------------------------------------------------
# get_switching_bracket
# ---------------------------------------------------------------------------

def test_get_switching_bracket_low_to_high():
    # ATP11C in artery coronary (paper Figure 1): switches at index 2 (40-49)
    assert get_switching_bracket([0, 0, 1, 1, 1, 1]) == 2


def test_get_switching_bracket_high_to_low():
    # ATP11C in brain amygdala (paper Figure 1): switches at index 3 (50-59)
    assert get_switching_bracket([1, 1, 1, 0, 0, 0]) == 3


def test_get_switching_bracket_returns_none_for_non_switching():
    assert get_switching_bracket([0, 0, 0, 0, 0, 0]) is None
    assert get_switching_bracket([0, 1, 0, 1, 0, 1]) is None
    assert get_switching_bracket([0, 0, 1, 1, 0, 0]) is None


def test_get_switching_bracket_grik2_example():
    # GRIK2 in brain cortex (paper Figure 1): vector [1,1,0,0,0,0]
    assert get_switching_bracket([1, 1, 0, 0, 0, 0]) == 2


def test_get_switching_bracket_min_and_max_positions():
    # Earliest possible switch: index 1
    assert get_switching_bracket([0, 1, 1, 1, 1, 1]) == 1
    # Latest possible switch: index 5
    assert get_switching_bracket([0, 0, 0, 0, 0, 1]) == 5


# ---------------------------------------------------------------------------
# identify_switching_genes: end-to-end on small DataFrames
# ---------------------------------------------------------------------------

def test_identify_switching_genes_basic():
    df = pd.DataFrame(
        {
            "20-29": [0.1, 0.9, 0.5, 0.0],
            "30-39": [0.1, 0.9, 0.5, 0.0],
            "40-49": [0.8, 0.9, 0.5, 0.0],  # gene1 switches here (low->high)
            "50-59": [0.8, 0.1, 0.5, 0.0],  # gene2 switches here (high->low)
            "60-69": [0.8, 0.1, 0.5, 0.0],
            "70-79": [0.8, 0.1, 0.5, 0.0],
        },
        index=["gene1", "gene2", "gene3_const", "gene4_const"],
    )
    result = identify_switching_genes(df, threshold=0.5)
    # Expected: gene1 in 40-49, gene2 in 50-59, gene3 and gene4 nowhere
    assert result["40-49"] == ["gene1"]
    assert result["50-59"] == ["gene2"]
    assert result["30-39"] == []
    assert result["60-69"] == []
    assert result["70-79"] == []


def test_identify_switching_genes_oscillating_excluded():
    df = pd.DataFrame(
        {
            "20-29": [0.1],
            "30-39": [0.9],
            "40-49": [0.1],
            "50-59": [0.9],
            "60-69": [0.1],
            "70-79": [0.9],
        },
        index=["oscillator"],
    )
    result = identify_switching_genes(df, threshold=0.5)
    # No bracket should contain the oscillating gene
    for bracket, genes in result.items():
        assert "oscillator" not in genes


def test_identify_switching_genes_nan_row_skipped():
    df = pd.DataFrame(
        {
            "20-29": [0.0, np.nan],
            "30-39": [0.0, np.nan],
            "40-49": [1.0, 0.5],
            "50-59": [1.0, 0.5],
            "60-69": [1.0, 0.5],
            "70-79": [1.0, 0.5],
        },
        index=["good_gene", "nan_gene"],
    )
    result = identify_switching_genes(df, threshold=0.5)
    assert "good_gene" in result["40-49"]
    # nan_gene must not appear anywhere
    for bracket, genes in result.items():
        assert "nan_gene" not in genes


def test_identify_switching_genes_threshold_validation():
    df = pd.DataFrame(
        {"20-29": [0.5], "30-39": [0.5]},
        index=["g1"],
    )
    with pytest.raises(ValueError, match=r"threshold must be in"):
        identify_switching_genes(df, threshold=1.5)
    with pytest.raises(ValueError, match=r"threshold must be in"):
        identify_switching_genes(df, threshold=-0.1)


def test_identify_switching_genes_returns_only_switching_brackets():
    df = pd.DataFrame(
        {
            "20-29": [0.0],
            "30-39": [1.0],
            "40-49": [1.0],
            "50-59": [1.0],
            "60-69": [1.0],
            "70-79": [1.0],
        },
        index=["early_switch"],
    )
    result = identify_switching_genes(df, threshold=0.5)
    # 20-29 must NOT be a key (only switching brackets are)
    assert "20-29" not in result
    # The gene switches at the first switching bracket
    assert result["30-39"] == ["early_switch"]


def test_identify_switching_genes_all_brackets_keys_present():
    df = pd.DataFrame(
        {
            "20-29": [0.0],
            "30-39": [0.0],
            "40-49": [0.0],
            "50-59": [0.0],
            "60-69": [0.0],
            "70-79": [0.0],
        },
        index=["constant_gene"],
    )
    result = identify_switching_genes(df, threshold=0.5)
    # All five SWITCHING_BRACKETS must be keys, even though all empty
    expected_keys = {"30-39", "40-49", "50-59", "60-69", "70-79"}
    assert set(result.keys()) == expected_keys
    for genes in result.values():
        assert genes == []


def test_identify_switching_genes_subset_of_columns_works():
    # User passes only 4 brackets instead of 6: still valid, just less data
    df = pd.DataFrame(
        {
            "30-39": [0.0],
            "40-49": [0.0],
            "50-59": [1.0],
            "60-69": [1.0],
        },
        index=["g1"],
    )
    result = identify_switching_genes(df, threshold=0.5)
    assert result["50-59"] == ["g1"]
    # 70-79 not in input, must not be in output
    assert "70-79" not in result


def test_identify_switching_genes_columns_must_be_in_order():
    df = pd.DataFrame(
        {
            "50-59": [0.0],
            "30-39": [1.0],  # order is wrong
            "40-49": [0.5],
            "20-29": [0.0],
        },
        index=["g1"],
    )
    with pytest.raises(ValueError, match="chronological order"):
        identify_switching_genes(df, threshold=0.5)


def test_identify_switching_genes_threshold_boundary():
    # Value exactly at threshold maps to 1 (>=, not >)
    df = pd.DataFrame(
        {
            "20-29": [0.0, 0.4],
            "30-39": [0.5, 0.5],
            "40-49": [0.5, 0.5],
            "50-59": [0.5, 0.5],
            "60-69": [0.5, 0.5],
            "70-79": [0.5, 0.5],
        },
        index=["boundary_low", "boundary_below"],
    )
    result = identify_switching_genes(df, threshold=0.5)
    # boundary_low: vector becomes [0,1,1,1,1,1] -> switch at 30-39
    assert "boundary_low" in result["30-39"]
    # boundary_below: vector becomes [0,1,1,1,1,1] (because 0.5 >= 0.5)
    # so it ALSO switches at 30-39
    assert "boundary_below" in result["30-39"]