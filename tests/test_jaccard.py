"""Tests for Jaccard tissue similarity."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from stamp.jaccard import (
    jaccard,
    jaccard_age_pair,
    jaccard_life_pair,
    tissue_similarity_age,
    tissue_similarity_life,
)


# ---------------------------------------------------------------------------
# Low-level jaccard()
# ---------------------------------------------------------------------------

def test_jaccard_identical_sets():
    assert jaccard({"a", "b", "c"}, {"a", "b", "c"}) == 1.0


def test_jaccard_disjoint_sets():
    assert jaccard({"a", "b"}, {"c", "d"}) == 0.0


def test_jaccard_partial_overlap():
    # 1 element in common, 3 in union -> 1/3
    assert jaccard({"a", "b"}, {"a", "c"}) == pytest.approx(1 / 3)


def test_jaccard_empty_sets_returns_zero():
    # By convention here: empty vs empty = 0
    assert jaccard(set(), set()) == 0.0


def test_jaccard_one_empty_set():
    assert jaccard({"a"}, set()) == 0.0
    assert jaccard(set(), {"a"}) == 0.0


def test_jaccard_subset():
    # {a, b} is a subset of {a, b, c, d} -> 2 / 4 = 0.5
    assert jaccard({"a", "b"}, {"a", "b", "c", "d"}) == 0.5


# ---------------------------------------------------------------------------
# jaccard_age_pair: per-pair similarity averaged across brackets
# ---------------------------------------------------------------------------

def test_age_pair_identical_tissues():
    sets = {
        "30-39": ["g1", "g2"],
        "40-49": ["g3"],
        "50-59": ["g4", "g5"],
        "60-69": [],
        "70-79": ["g6"],
    }
    assert jaccard_age_pair(sets, sets) == pytest.approx(4 / 5)
    # The empty bracket gives Jaccard(empty,empty)=0 so result is 4/5


def test_age_pair_disjoint_in_every_bracket():
    a = {"30-39": ["g1"], "40-49": ["g2"], "50-59": ["g3"],
         "60-69": ["g4"], "70-79": ["g5"]}
    b = {"30-39": ["x1"], "40-49": ["x2"], "50-59": ["x3"],
         "60-69": ["x4"], "70-79": ["x5"]}
    assert jaccard_age_pair(a, b) == 0.0


def test_age_pair_partial_overlap_one_bracket():
    a = {"30-39": ["g1", "g2"], "40-49": [], "50-59": [],
         "60-69": [], "70-79": []}
    b = {"30-39": ["g1"], "40-49": [], "50-59": [],
         "60-69": [], "70-79": []}
    # 1/2 in 30-39, then 0+0+0+0 -> mean = (0.5 + 0 + 0 + 0 + 0) / 5
    assert jaccard_age_pair(a, b) == pytest.approx(0.1)


def test_age_pair_missing_bracket_in_one_tissue():
    # Tissue B is missing 70-79 entirely (incomplete sampling).
    # That bracket must be skipped, not counted as 0.
    a = {"30-39": ["g1"], "40-49": [], "50-59": [],
         "60-69": [], "70-79": ["g2"]}
    b = {"30-39": ["g1"], "40-49": [], "50-59": [],
         "60-69": []}  # no 70-79
    # Brackets in common: 30-39 (1.0), 40-49 (0), 50-59 (0), 60-69 (0)
    # Mean over 4 brackets = 0.25
    assert jaccard_age_pair(a, b) == pytest.approx(0.25)


def test_age_pair_no_common_brackets():
    # Pathological edge case: no shared brackets at all
    a = {"30-39": ["g1"]}
    b = {"50-59": ["g1"]}
    # No bracket is in both dicts
    # ... but the way `b` is defined here, only 50-59. The function
    # iterates over SWITCHING_BRACKETS. The intersection of keys is
    # empty -> returns 0.0
    # However in this case both ARE in SWITCHING_BRACKETS and they
    # are different. The function should find no bracket present
    # in BOTH dicts. Let's check.
    result = jaccard_age_pair(a, b)
    assert result == 0.0


# ---------------------------------------------------------------------------
# jaccard_life_pair: union over brackets
# ---------------------------------------------------------------------------

def test_life_pair_identical():
    sets = {"30-39": ["g1"], "40-49": ["g2"], "50-59": ["g3"],
            "60-69": ["g4"], "70-79": ["g5"]}
    assert jaccard_life_pair(sets, sets) == 1.0


def test_life_pair_disjoint():
    a = {"30-39": ["g1"], "40-49": ["g2"]}
    b = {"30-39": ["x1"], "40-49": ["x2"]}
    assert jaccard_life_pair(a, b) == 0.0


def test_life_pair_partial():
    # Union A = {g1, g2, g3}, Union B = {g2, g3, g4}, intersection = 2, union = 4
    a = {"30-39": ["g1", "g2"], "40-49": ["g3"]}
    b = {"30-39": ["g2"], "40-49": ["g3", "g4"]}
    assert jaccard_life_pair(a, b) == pytest.approx(2 / 4)


def test_life_pair_empty_tissues():
    a = {"30-39": [], "40-49": []}
    b = {"30-39": [], "40-49": []}
    assert jaccard_life_pair(a, b) == 0.0


# ---------------------------------------------------------------------------
# Full matrices
# ---------------------------------------------------------------------------

def test_similarity_age_matrix_shape_and_diagonal():
    sets_by_tissue = {
        "Liver":          {"30-39": ["g1"], "40-49": ["g2"], "50-59": [], "60-69": [], "70-79": []},
        "Brain_Cortex":   {"30-39": ["g1"], "40-49": ["g3"], "50-59": [], "60-69": [], "70-79": []},
        "Whole_Blood":    {"30-39": ["g4"], "40-49": [], "50-59": [], "60-69": [], "70-79": []},
    }
    mat = tissue_similarity_age(sets_by_tissue)
    # Shape 3x3
    assert mat.shape == (3, 3)
    # Diagonal is 1
    for t in mat.index:
        assert mat.loc[t, t] == pytest.approx(1.0)
    # Symmetric
    for i in mat.index:
        for j in mat.columns:
            assert mat.loc[i, j] == pytest.approx(mat.loc[j, i])


def test_similarity_age_known_value():
    sets_by_tissue = {
        "T1": {"30-39": ["g1"], "40-49": [], "50-59": [],
               "60-69": [], "70-79": []},
        "T2": {"30-39": ["g1"], "40-49": [], "50-59": [],
               "60-69": [], "70-79": []},
    }
    mat = tissue_similarity_age(sets_by_tissue)
    # Both tissues identical -> J_age(T1, T2) = 1/5 = 0.2
    # because 30-39 = 1.0 and the other 4 brackets are 0 (empty vs empty)
    assert mat.loc["T1", "T2"] == pytest.approx(0.2)


def test_similarity_life_matrix_shape_and_diagonal():
    sets_by_tissue = {
        "Liver":        {"30-39": ["g1"], "40-49": ["g2"]},
        "Brain_Cortex": {"30-39": ["g1"], "40-49": ["g3"]},
    }
    mat = tissue_similarity_life(sets_by_tissue)
    assert mat.shape == (2, 2)
    assert mat.loc["Liver", "Liver"] == pytest.approx(1.0)
    assert mat.loc["Brain_Cortex", "Brain_Cortex"] == pytest.approx(1.0)
    # Union(Liver) = {g1, g2}, Union(Brain_Cortex) = {g1, g3}
    # Intersection = {g1}, union = {g1, g2, g3} -> 1/3
    assert mat.loc["Liver", "Brain_Cortex"] == pytest.approx(1 / 3)


def test_matrix_is_sorted_alphabetically():
    sets_by_tissue = {
        "Zebra":  {"30-39": ["g1"]},
        "Alpha":  {"30-39": ["g2"]},
        "Mango":  {"30-39": ["g3"]},
    }
    mat = tissue_similarity_age(sets_by_tissue)
    assert list(mat.index) == ["Alpha", "Mango", "Zebra"]
    assert list(mat.columns) == ["Alpha", "Mango", "Zebra"]


def test_matrix_dtype_is_float32():
    sets_by_tissue = {
        "T1": {"30-39": ["g1"]},
        "T2": {"30-39": ["g2"]},
    }
    mat = tissue_similarity_age(sets_by_tissue)
    assert mat.dtypes.iloc[0] == np.float32