"""Tests for complete age-bin coverage helpers and path routing."""
from __future__ import annotations

import pandas as pd
import pytest

from stamp.config import AGE_BRACKETS, paths_for
from stamp.tissues import (
    tissues_with_complete_age_bins,
    incomplete_tissues,
    common_complete_tissues,
)


def _meta(rows: list[tuple[str, str]]) -> pd.DataFrame:
    """Build a metadata DataFrame from (tissue, age_bracket) tuples."""
    return pd.DataFrame(
        [
            {"sample_id": f"S{i}", "tissue": t, "age_bracket": b}
            for i, (t, b) in enumerate(rows)
        ]
    )


# ---------------------------------------------------------------------------
# tissues_with_complete_age_bins
# ---------------------------------------------------------------------------

def test_includes_tissue_with_all_six_brackets():
    rows = [("Liver", b) for b in AGE_BRACKETS]            # all 6
    rows += [("Liver", "40-49")]                            # an extra dup is fine
    out = tissues_with_complete_age_bins(_meta(rows))
    assert out == ["Liver"]


def test_excludes_tissue_missing_one_bracket():
    # Bladder missing 70-79
    rows = [("Bladder", b) for b in AGE_BRACKETS if b != "70-79"]
    out = tissues_with_complete_age_bins(_meta(rows))
    assert "Bladder" not in out
    assert out == []


def test_ignores_age_bracket_outside_AGE_BRACKETS():
    # A tissue with the 6 valid brackets PLUS an out-of-range "80-89" is complete;
    # a tissue whose only "coverage" of a bracket is an out-of-range label is not.
    rows = [("Lung", b) for b in AGE_BRACKETS] + [("Lung", "80-89")]
    rows += [("Weird", b) for b in AGE_BRACKETS if b != "70-79"]
    rows += [("Weird", "80-89")]  # does NOT count as 70-79
    out = tissues_with_complete_age_bins(_meta(rows))
    assert "Lung" in out
    assert "Weird" not in out


def test_ignores_null_age_bracket():
    rows = [("Ovary", b) for b in AGE_BRACKETS if b != "60-69"]
    rows += [("Ovary", None)]  # null does not fill 60-69
    out = tissues_with_complete_age_bins(_meta(rows))
    assert "Ovary" not in out


def test_multiple_tissues_mixed():
    rows = []
    rows += [("A", b) for b in AGE_BRACKETS]                       # complete
    rows += [("B", b) for b in AGE_BRACKETS if b != "30-39"]       # incomplete
    rows += [("C", b) for b in AGE_BRACKETS]                       # complete
    out = tissues_with_complete_age_bins(_meta(rows))
    assert out == ["A", "C"]
    assert incomplete_tissues(_meta(rows)) == ["B"]


def test_missing_columns_raises():
    with pytest.raises(ValueError, match="missing required columns"):
        tissues_with_complete_age_bins(pd.DataFrame({"tissue": ["X"]}))


# ---------------------------------------------------------------------------
# common_complete_tissues (v8 vs v10 like-for-like)
# ---------------------------------------------------------------------------

def test_common_complete_is_intersection():
    # v8: A,B complete (C missing a bracket).  v10: A,C complete (B missing).
    m8 = _meta(
        [("A", b) for b in AGE_BRACKETS]
        + [("B", b) for b in AGE_BRACKETS]
        + [("C", b) for b in AGE_BRACKETS if b != "70-79"]
    )
    m10 = _meta(
        [("A", b) for b in AGE_BRACKETS]
        + [("C", b) for b in AGE_BRACKETS]
        + [("B", b) for b in AGE_BRACKETS if b != "70-79"]
    )
    assert common_complete_tissues(m8, m10) == ["A"]


# ---------------------------------------------------------------------------
# paths_for routing
# ---------------------------------------------------------------------------

def test_paths_for_complete_routes_to_separate_dir():
    normal = paths_for("v10", complete=False)
    comp = paths_for("v10", complete=True)
    # Derived outputs go to a separate *_complete directory.
    for key in ("normalized", "sets", "jaccard", "migration", "threshold_sweep"):
        assert normal[key].parent.name == "v10"
        assert comp[key].parent.name == "v10_complete"


def test_paths_for_complete_shares_raw_inputs():
    normal = paths_for("v8", complete=False)
    comp = paths_for("v8", complete=True)
    # Raw inputs are shared between modes (same parquet under data/parquet/v8).
    assert normal["tpm_parquet"] == comp["tpm_parquet"]
    assert normal["metadata_parquet"] == comp["metadata_parquet"]


def test_paths_for_default_is_all_tissues():
    # Default (no flag) must be byte-identical to the historic behaviour.
    assert paths_for("v10")["sets"].parent.name == "v10"
