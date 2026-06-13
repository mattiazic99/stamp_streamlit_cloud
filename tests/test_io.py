"""Tests for the I/O module: round-trip Parquet and sets-txt."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from stamp.config import SWITCHING_BRACKETS
from stamp.io import (
    _safe_filename,
    load_normalized_tissue,
    load_sets_txt,
    save_normalized_tissue,
    save_sets_txt,
)


# ---------------------------------------------------------------------------
# Filename sanitisation
# ---------------------------------------------------------------------------

def test_safe_filename_simple():
    assert _safe_filename("Liver") == "Liver"


def test_safe_filename_spaces_and_dashes():
    assert _safe_filename("Brain - Substantia nigra") == "Brain_Substantia_nigra"


def test_safe_filename_parentheses():
    assert _safe_filename("Skin - Sun Exposed (Lower leg)") == "Skin_Sun_Exposed_Lower_leg"


def test_safe_filename_collapses_underscores():
    assert _safe_filename("a -- b") == "a_b"


def test_safe_filename_strips_leading_trailing():
    assert _safe_filename(" -- xyz -- ") == "xyz"


# ---------------------------------------------------------------------------
# Sets file: write + read round-trip
# ---------------------------------------------------------------------------

def test_save_and_load_sets_full(monkeypatch, tmp_path):
    """Round-trip: write a sets file, read it back, check identity."""
    _redirect_paths_to_tmp(monkeypatch, tmp_path)

    sets = {
        "30-39": ["g1", "g2", "g3"],
        "40-49": [],
        "50-59": ["g4"],
        "60-69": ["g5", "g6"],
        "70-79": [],
    }
    written = save_sets_txt("v10", "Liver", sets)
    assert written.exists()

    loaded = load_sets_txt("v10", "Liver")
    assert loaded == sets


def test_save_sets_with_partial_dict(monkeypatch, tmp_path):
    """Missing brackets in input dict produce empty lines on disk."""
    _redirect_paths_to_tmp(monkeypatch, tmp_path)

    # Only 2 of 5 brackets provided
    sets = {"30-39": ["g1"], "60-69": ["g2"]}
    save_sets_txt("v10", "Liver", sets)

    loaded = load_sets_txt("v10", "Liver")
    assert loaded["30-39"] == ["g1"]
    assert loaded["40-49"] == []
    assert loaded["50-59"] == []
    assert loaded["60-69"] == ["g2"]
    assert loaded["70-79"] == []


def test_save_sets_rejects_unknown_bracket(monkeypatch, tmp_path):
    _redirect_paths_to_tmp(monkeypatch, tmp_path)
    with pytest.raises(ValueError, match="unknown brackets"):
        save_sets_txt("v10", "Liver", {"99-99": ["g1"]})

def test_save_sets_writes_exactly_5_lines(monkeypatch, tmp_path):
    """The sets file MUST have 5 lines, one per SWITCHING_BRACKET."""
    _redirect_paths_to_tmp(monkeypatch, tmp_path)
    save_sets_txt("v10", "Liver", {"30-39": ["g1"]})

    from stamp.config import paths_for
    path = paths_for("v10")["sets"] / "Liver_sets.txt"
    text = path.read_text(encoding="utf-8")
    # Convention: the file ends with a single trailing newline.
    # Stripping that trailing newline must leave exactly 5 lines (some empty).
    assert text.endswith("\n")
    lines = text[:-1].split("\n")
    assert len(lines) == len(SWITCHING_BRACKETS)
    # First line should contain g1, the rest should be empty
    assert lines[0] == "g1"
    assert all(line == "" for line in lines[1:])


def test_load_sets_rejects_wrong_line_count(monkeypatch, tmp_path):
    _redirect_paths_to_tmp(monkeypatch, tmp_path)
    from stamp.config import paths_for

    out_dir = paths_for("v10")["sets"]
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "Liver_sets.txt"
    # Write only 3 lines instead of 5
    path.write_text("a b\nc\nd\n", encoding="utf-8")

    with pytest.raises(ValueError, match="expected"):
        load_sets_txt("v10", "Liver")


# ---------------------------------------------------------------------------
# Normalized matrix: write + read round-trip
# ---------------------------------------------------------------------------
def test_save_and_load_normalized_round_trip(monkeypatch, tmp_path):
    _redirect_paths_to_tmp(monkeypatch, tmp_path)

    df = pd.DataFrame(
        {
            "20-29": [0.0, 0.5],
            "30-39": [1.0, 0.5],
            "40-49": [1.0, 0.5],
        },
        index=["gene_A", "gene_B"],
    ).astype(np.float32)
    df.index.name = "gene_id"

    save_normalized_tissue("v10", "Liver", df)
    loaded = load_normalized_tissue("v10", "Liver")

    # load_normalized_tissue returns gene_id as index, so direct compare
    pd.testing.assert_frame_equal(
        loaded,
        df,
        check_exact=False,
        rtol=1e-5,
    )


def test_load_normalized_missing_file_raises(monkeypatch, tmp_path):
    _redirect_paths_to_tmp(monkeypatch, tmp_path)
    with pytest.raises(FileNotFoundError, match="Normalized Parquet"):
        load_normalized_tissue("v10", "NonExistentTissue")


def test_load_tpm_missing_file_raises(monkeypatch, tmp_path):
    _redirect_paths_to_tmp(monkeypatch, tmp_path)
    from stamp.io import load_tpm_matrix
    with pytest.raises(FileNotFoundError, match="TPM Parquet"):
        load_tpm_matrix("v10")


def test_load_metadata_missing_file_raises(monkeypatch, tmp_path):
    _redirect_paths_to_tmp(monkeypatch, tmp_path)
    from stamp.io import load_metadata
    with pytest.raises(FileNotFoundError, match="Metadata Parquet"):
        load_metadata("v10")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _redirect_paths_to_tmp(monkeypatch, tmp_path):
    """Redirect stamp.config paths to a pytest tmp directory for isolation.

    This prevents tests from polluting the real output/ folder.
    """
    import stamp.config as cfg

    monkeypatch.setattr(cfg, "ROOT", tmp_path)
    monkeypatch.setattr(cfg, "DATA_DIR", tmp_path / "data")
    monkeypatch.setattr(cfg, "PARQUET_DIR", tmp_path / "data" / "parquet")
    monkeypatch.setattr(cfg, "EXTERNAL_DIR", tmp_path / "data" / "external")
    monkeypatch.setattr(cfg, "OUTPUT_DIR", tmp_path / "output")
    monkeypatch.setattr(cfg, "COMPARE_DIR", tmp_path / "output" / "compare")