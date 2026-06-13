"""Integration checks on the generated complete-age-bins outputs.

These tests read the pre-computed CSV/sets outputs directly (pandas only, no
Streamlit) and verify the complete-mode artifacts are well-formed:

  * the complete Jaccard matrices are square and exclude incomplete tissues;
  * v8 has 49 complete tissues, v10 has 50;
  * the like-for-like v8-vs-v10 comparison set is 49 tissues;
  * complete Jaccard values equal the all-tissues matrix sliced to the
    complete tissues (per-pair Jaccard is independent of the tissue set).

If the complete outputs have not been generated yet (output/{version}_complete
missing), the tests are skipped.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "output"

# Incomplete tissues (filesystem-safe names) expected to be excluded.
INCOMPLETE_SAFE = {
    "v8": {"Bladder", "Cervix_Ectocervix", "Cervix_Endocervix",
           "Fallopian_Tube", "Kidney_Medulla"},
    "v10": {"Cervix_Ectocervix", "Cervix_Endocervix",
            "Fallopian_Tube", "Kidney_Medulla"},
}
EXPECTED_N_COMPLETE = {"v8": 49, "v10": 50}


def _complete_jaccard(version: str, metric: str = "life") -> pd.DataFrame:
    path = OUTPUT / f"{version}_complete" / "jaccard" / f"jaccard_{metric}.csv"
    if not path.exists():
        pytest.skip(f"complete outputs not generated: {path} missing")
    return pd.read_csv(path, index_col=0)


@pytest.mark.parametrize("version", ["v8", "v10"])
def test_complete_jaccard_is_square_expected_size(version):
    m = _complete_jaccard(version)
    assert m.shape[0] == m.shape[1]
    assert m.shape[0] == EXPECTED_N_COMPLETE[version]


@pytest.mark.parametrize("version", ["v8", "v10"])
def test_complete_jaccard_excludes_incomplete_tissues(version):
    m = _complete_jaccard(version)
    assert INCOMPLETE_SAFE[version].isdisjoint(set(m.index))
    assert INCOMPLETE_SAFE[version].isdisjoint(set(m.columns))


def test_comparison_set_is_49_common_tissues():
    m8 = _complete_jaccard("v8")
    m10 = _complete_jaccard("v10")
    common = set(m8.index) & set(m10.index)
    assert len(common) == 49
    # Bladder is complete in v10 but not v8, so it is NOT in the comparison set.
    assert "Bladder" not in common


@pytest.mark.parametrize("version", ["v8", "v10"])
@pytest.mark.parametrize("metric", ["life", "age"])
def test_complete_matches_full_sliced(version, metric):
    full_path = OUTPUT / version / "jaccard" / f"jaccard_{metric}.csv"
    if not full_path.exists():
        pytest.skip("all-tissues matrix missing")
    comp = _complete_jaccard(version, metric)
    full = pd.read_csv(full_path, index_col=0)
    sub = full.loc[comp.index, comp.columns]
    # Per-pair Jaccard does not depend on which other tissues are present.
    assert np.allclose(sub.values, comp.values, atol=1e-6)
