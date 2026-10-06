from __future__ import annotations

import numpy as np
import pytest

from statistical_validation_cli.common import (
    age_means_from_codes,
    bh_qvalues,
    empirical_right_tail_p,
    enforce_final_scope,
    normalize_sample_values,
    switching_state,
)
from statistical_validation_cli.smote import _strategy_and_neighbors


def test_final_scope_rejects_partial_tissue_sets() -> None:
    with pytest.raises(RuntimeError, match="protected final run"):
        enforce_final_scope(["Liver"], ["Liver", "Lung"], "final", "Test")
    enforce_final_scope(["Liver"], ["Liver", "Lung"], "pilot", "Test")


def test_empirical_right_tail_p_uses_plus_one() -> None:
    null = np.array([1, 2, 3, 4])
    assert empirical_right_tail_p(5, null) == 1 / 5
    assert empirical_right_tail_p(3, null) == 3 / 5


def test_bh_qvalues_are_aligned_and_monotone_by_rank() -> None:
    p = np.array([0.032, 0.001, 0.20, 0.006, 0.012])
    q = bh_qvalues(p)
    expected = np.array([0.040, 0.005, 0.20, 0.015, 0.020])
    np.testing.assert_allclose(q, expected)
    order = np.argsort(p)
    assert np.all(np.diff(q[order]) >= 0)


def test_age_means_match_nan_aware_group_means() -> None:
    values = np.array(
        [
            [1.0, np.nan, 3.0, 4.0, 5.0, 6.0, 7.0],
            [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0],
        ],
        dtype=np.float32,
    )
    codes = np.array([0, 0, 1, 2, 3, 4, 5], dtype=np.int8)
    expected = np.array(
        [
            [1.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            [3.0, 6.0, 8.0, 10.0, 12.0, 14.0],
        ],
        dtype=np.float32,
    )
    np.testing.assert_allclose(
        age_means_from_codes(values, codes), expected, equal_nan=True
    )


def test_switching_state_detects_one_transition_and_direction() -> None:
    means = np.array(
        [
            [0.1, 0.2, 0.8, 0.9, 0.9, 1.0],
            [0.9, 0.8, 0.7, 0.2, 0.1, 0.0],
            [0.1, 0.8, 0.2, 0.9, 0.1, 0.8],
            [0.1, 0.2, 0.3, 0.4, 0.4, 0.4],
        ],
        dtype=np.float32,
    )
    state = switching_state(means, tau=0.5)
    np.testing.assert_array_equal(state.mask, [True, True, False, False])
    np.testing.assert_array_equal(state.position[:2], [2, 3])
    np.testing.assert_array_equal(state.direction[:2], [1, -1])


def test_normalization_applies_epsilon_and_minmax() -> None:
    raw = np.array(
        [
            [1.0, 2.0, 3.0],
            [5.0, 5.0, 5.0],
            [2.000, 2.005, 2.009],
        ],
        dtype=np.float32,
    )
    normalized, kept = normalize_sample_values(raw, epsilon=0.01)
    np.testing.assert_array_equal(kept, [0])
    np.testing.assert_allclose(normalized[0], [0.0, 0.5, 1.0])


def test_smote_strategy_is_prespecified_and_singletons_are_untouched() -> None:
    counts = np.array([1, 2, 5, 19, 20, 100])
    strategy, k, singletons = _strategy_and_neighbors(counts, target_n=20)
    assert strategy == {1: 20, 2: 20, 3: 20}
    assert k == 1
    assert singletons == ["20-29"]


def test_adaptive_smote_raises_isolated_minimum_to_second_smallest() -> None:
    counts = np.array([40, 42, 63, 121, 139, 14])
    strategy, k, singletons = _strategy_and_neighbors(
        counts,
        target_n=20,
        target_policy="adaptive-second-smallest",
    )
    assert strategy == {5: 40}
    assert k == 5
    assert singletons == []


def test_adaptive_smote_uses_fixed_floor_when_two_bins_are_small() -> None:
    counts = np.array([8, 18, 38, 99, 93, 6])
    strategy, k, singletons = _strategy_and_neighbors(
        counts,
        target_n=20,
        target_policy="adaptive-second-smallest",
    )
    assert strategy == {0: 20, 1: 20, 5: 20}
    assert k == 5
    assert singletons == []


def test_adaptive_smote_is_noop_when_all_bins_reach_floor() -> None:
    counts = np.array([26, 66, 70, 116, 256, 269])
    strategy, k, singletons = _strategy_and_neighbors(
        counts,
        target_n=20,
        target_policy="adaptive-second-smallest",
    )
    assert strategy == {}
    assert k is None
    assert singletons == []


def test_adaptive_smote_does_not_interpolate_isolated_singleton() -> None:
    counts = np.array([1, 22, 35, 41, 56, 70])
    strategy, k, singletons = _strategy_and_neighbors(
        counts,
        target_n=20,
        target_policy="adaptive-second-smallest",
    )
    assert strategy == {}
    assert k is None
    assert singletons == ["20-29"]
