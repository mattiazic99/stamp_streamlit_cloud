"""Smoke tests for data_loader.py — run with:

    python -m pytest gui/test_data_loader.py -v

or simply:

    python gui/test_data_loader.py
"""
import sys
from pathlib import Path

# Ensure stamp package is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

# Bypass Streamlit caching for tests: call __wrapped__ on cached fns
from data_loader import (
    get_available_tissues,
    get_available_tissues_safe,
    get_jaccard_matrix,
    get_sets_for_tissue,
    get_threshold_sweep,
    get_stability_scores,
    get_migration_paths,
    sets_to_gui_format,
    sets_to_gene_sets,
    safe_to_display,
    display_to_safe,
)

# ── Helpers ───────────────────────────────────────────────────────
# @st.cache_data wraps functions; __wrapped__ gives the original fn
_tissues = get_available_tissues.__wrapped__
_tissues_safe = get_available_tissues_safe.__wrapped__
_jaccard = get_jaccard_matrix.__wrapped__
_sets = get_sets_for_tissue.__wrapped__
_sweep = get_threshold_sweep.__wrapped__
_stability = get_stability_scores.__wrapped__
_migration = get_migration_paths.__wrapped__


def test_get_available_tissues_v10():
    tissues = _tissues("v10")
    assert len(tissues) == 54, f"Expected 54 tissues for v10, got {len(tissues)}"
    # Should be display names (spaces/hyphens)
    assert any(" " in t or "-" in t for t in tissues), "Names should be display-formatted"


def test_get_available_tissues_v8():
    tissues = _tissues("v8")
    assert len(tissues) == 54, f"Expected 54 tissues for v8, got {len(tissues)}"


def test_jaccard_matrix_shape_v10():
    df = _jaccard("v10", "life")
    assert df.shape == (54, 54), f"Expected (54, 54), got {df.shape}"
    # Diagonal must be 1.0 (self-similarity)
    import numpy as np
    assert np.allclose(df.values.diagonal(), 1.0), "Diagonal should be 1.0"


def test_jaccard_matrix_shape_v8():
    df = _jaccard("v8", "life")
    assert df.shape == (54, 54), f"Expected (54, 54), got {df.shape}"


def test_sets_for_tissue():
    tissues = _tissues("v10")
    first = tissues[0]
    sets_dict = _sets("v10", first)
    assert isinstance(sets_dict, dict), "Should return a dict"
    # Should have exactly 5 brackets
    assert len(sets_dict) == 5, f"Expected 5 brackets, got {len(sets_dict)}"
    from stamp.config import SWITCHING_BRACKETS
    for b in SWITCHING_BRACKETS:
        assert b in sets_dict, f"Missing bracket {b}"


def test_sets_to_gui_format():
    tissues = _tissues("v10")
    first = tissues[0]
    sets_dict = _sets("v10", first)
    age_groups, counts, df = sets_to_gui_format(sets_dict)
    assert len(age_groups) == 5
    assert len(counts) == 5
    assert "Age" in df.columns
    assert "Gene" in df.columns
    # En-dash format
    assert "–" in age_groups[0], f"Expected en-dash, got {age_groups[0]}"


def test_sets_to_gene_sets():
    tissues = _tissues("v10")
    first = tissues[0]
    sets_dict = _sets("v10", first)
    gene_sets = sets_to_gene_sets(sets_dict)
    assert len(gene_sets) == 5
    assert all(isinstance(s, set) for s in gene_sets)


def test_safe_to_display_roundtrip():
    assert safe_to_display("Brain_Cortex") == "Brain - Cortex"
    assert display_to_safe("Brain - Cortex") == "Brain_Cortex"
    # Roundtrip
    for safe in _tissues_safe("v10"):
        disp = safe_to_display(safe)
        assert display_to_safe(disp) == safe, f"Roundtrip failed: {safe} → {disp} → {display_to_safe(disp)}"


def test_threshold_sweep():
    df = _sweep("v10")
    assert "tissue" in df.columns
    assert "tau" in df.columns
    assert "n_switching" in df.columns
    assert len(df) > 0


def test_stability_scores():
    df = _stability("v10")
    assert "stability_cv" in df.columns
    assert len(df) == 54


def test_migration_paths():
    df = _migration("v10")
    assert "source_tissue" in df.columns
    assert "target_tissue" in df.columns
    assert "n_genes" in df.columns


# ── Run as script ─────────────────────────────────────────────────
if __name__ == "__main__":
    tests = [
        test_get_available_tissues_v10,
        test_get_available_tissues_v8,
        test_jaccard_matrix_shape_v10,
        test_jaccard_matrix_shape_v8,
        test_sets_for_tissue,
        test_sets_to_gui_format,
        test_sets_to_gene_sets,
        test_safe_to_display_roundtrip,
        test_threshold_sweep,
        test_stability_scores,
        test_migration_paths,
    ]
    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            print(f"  PASS {test.__name__}")
            passed += 1
        except Exception as e:
            print(f"  FAIL {test.__name__}: {e}")
            failed += 1

    print(f"\n{'='*50}")
    print(f"  Results: {passed} passed, {failed} failed, {len(tests)} total")
    if failed == 0:
        print("  All tests passed!")
    else:
        print("  Some tests failed.")
