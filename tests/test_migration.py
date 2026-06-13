"""Tests for migration path detection."""
from __future__ import annotations

import pandas as pd
import pytest

from stamp.migration import find_migrations, migration_paths, migration_summary


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def simple_sets() -> dict[str, dict[str, list[str]]]:
    """Two tissues sharing gene g1, switching at different brackets."""
    return {
        "Liver": {
            "30-39": ["g1", "g2"],
            "40-49": [],
            "50-59": [],
            "60-69": [],
            "70-79": [],
        },
        "Brain_Cortex": {
            "30-39": [],
            "40-49": ["g1"],    # g1 switches LATER in Brain than in Liver
            "50-59": ["g3"],
            "60-69": [],
            "70-79": [],
        },
        "Pancreas": {
            "30-39": [],
            "40-49": [],
            "50-59": ["g1"],    # g1 switches even later in Pancreas
            "60-69": [],
            "70-79": [],
        },
    }


# ---------------------------------------------------------------------------
# find_migrations
# ---------------------------------------------------------------------------

def test_migration_requires_earlier_source_age():
    sets = simple_sets()
    mig = find_migrations(sets)
    # g1 in Liver@30-39 -> Brain_Cortex@40-49: valid (30-39 < 40-49)
    liver_to_brain = mig[
        (mig["source_tissue"] == "Liver") &
        (mig["target_tissue"] == "Brain_Cortex") &
        (mig["gene_id"] == "g1")
    ]
    assert len(liver_to_brain) == 1
    # g1 in Brain_Cortex@40-49 -> Liver@30-39 must NOT appear (wrong order)
    brain_to_liver = mig[
        (mig["source_tissue"] == "Brain_Cortex") &
        (mig["target_tissue"] == "Liver") &
        (mig["gene_id"] == "g1")
    ]
    assert len(brain_to_liver) == 0


def test_migration_no_self_loops():
    sets = {
        "Liver": {"30-39": ["g1"], "40-49": ["g1"],
                  "50-59": [], "60-69": [], "70-79": []},
    }
    # g1 switches in Liver at two brackets but same tissue -> no migration
    mig = find_migrations(sets)
    self_loops = mig[mig["source_tissue"] == mig["target_tissue"]]
    assert len(self_loops) == 0


def test_migration_empty_when_no_overlap():
    sets = {
        "Liver":        {"30-39": ["g1"], "40-49": [], "50-59": [],
                         "60-69": [], "70-79": []},
        "Brain_Cortex": {"30-39": ["g2"], "40-49": [], "50-59": [],
                         "60-69": [], "70-79": []},
    }
    # g1 and g2 never appear together -> no migrations
    mig = find_migrations(sets)
    assert len(mig) == 0


def test_migration_counts_genes_correctly():
    sets = simple_sets()
    mig = find_migrations(sets)
    # g2 only in Liver -> no migration partner -> not in migrations
    g2_mig = mig[mig["gene_id"] == "g2"]
    assert len(g2_mig) == 0
    # g1 migrates: Liver@30-39 -> Brain@40-49, Liver@30-39 -> Pancreas@50-59,
    #              Brain@40-49 -> Pancreas@50-59
    g1_mig = mig[mig["gene_id"] == "g1"]
    assert len(g1_mig) == 3


def test_migration_requires_earlier_source_not_same_bracket():
    # Same bracket in two tissues: no migration (need i < j strictly)
    sets = {
        "Liver":        {"30-39": ["g1"], "40-49": [], "50-59": [],
                         "60-69": [], "70-79": []},
        "Brain_Cortex": {"30-39": ["g1"], "40-49": [], "50-59": [],
                         "60-69": [], "70-79": []},
    }
    mig = find_migrations(sets)
    # Same bracket = order_a == order_b -> skipped
    assert len(mig) == 0


def test_migration_empty_sets_no_crash():
    sets = {
        "Liver":        {"30-39": [], "40-49": [], "50-59": [],
                         "60-69": [], "70-79": []},
        "Brain_Cortex": {"30-39": [], "40-49": [], "50-59": [],
                         "60-69": [], "70-79": []},
    }
    mig = find_migrations(sets)
    assert len(mig) == 0
    assert list(mig.columns) == [
        "source_tissue", "source_bracket",
        "target_tissue", "target_bracket", "gene_id"
    ]


def test_migration_output_columns():
    mig = find_migrations(simple_sets())
    assert set(mig.columns) == {
        "source_tissue", "source_bracket",
        "target_tissue", "target_bracket", "gene_id"
    }


# ---------------------------------------------------------------------------
# migration_paths
# ---------------------------------------------------------------------------

def test_migration_aggregation_top_n():
    sets = simple_sets()
    mig = find_migrations(sets)
    paths = migration_paths(mig, top_n=2)
    assert len(paths) <= 2
    # top path should have the highest n_genes
    assert paths["n_genes"].is_monotonic_decreasing


def test_migration_paths_columns():
    mig = find_migrations(simple_sets())
    paths = migration_paths(mig)
    assert "n_genes" in paths.columns
    assert "genes" in paths.columns


def test_migration_paths_empty_input():
    empty = pd.DataFrame(columns=[
        "source_tissue", "source_bracket",
        "target_tissue", "target_bracket", "gene_id"
    ])
    paths = migration_paths(empty)
    assert len(paths) == 0


# ---------------------------------------------------------------------------
# migration_summary
# ---------------------------------------------------------------------------

def test_migration_summary_basic():
    mig = find_migrations(simple_sets())
    summary = migration_summary(mig)
    assert summary["n_migrating_genes"] == 1  # only g1 migrates
    assert summary["n_migration_events"] == 3  # 3 pairs involving g1
    assert summary["n_source_tissues"] >= 1
    assert summary["n_target_tissues"] >= 1


def test_migration_summary_empty():
    empty = pd.DataFrame(columns=[
        "source_tissue", "source_bracket",
        "target_tissue", "target_bracket", "gene_id"
    ])
    summary = migration_summary(empty)
    assert summary["n_migration_events"] == 0
    assert summary["n_migrating_genes"] == 0