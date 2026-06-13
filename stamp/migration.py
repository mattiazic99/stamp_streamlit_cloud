"""Cross-tissue migration of switching genes across age brackets.

A gene g "migrates" from tissue A at bracket i to tissue B at bracket j
if:
  - g switches in tissue A at bracket i (g ∈ S_{A,i})
  - g switches in tissue B at bracket j (g ∈ S_{B,j})
  - i < j  (the switch in A precedes the switch in B)

This models the observation that a gene's switching event propagates
across tissues as the organism ages, which the paper terms "gene
switching migration" (Figure 13 of Kahveci et al. 2025).

Output tables
-------------
migrations_raw.csv :  one row per (source_tissue, source_bracket,
                       target_tissue, target_bracket, gene_id)
migration_paths.csv : aggregated by (source_tissue, source_bracket,
                       target_tissue, target_bracket), sorted by
                       n_genes descending.

Note on the paper's Figure 13
------------------------------
The paper shows the top migration paths aggregated across ALL bracket
pairs (i, j) with i < j. The `migration_paths` function does the same.
"""
from __future__ import annotations

import pandas as pd

from stamp.config import SWITCHING_BRACKETS


# ---------------------------------------------------------------------------
# Core migration detection
# ---------------------------------------------------------------------------

def find_migrations(
    sets_by_tissue: dict[str, dict[str, list[str]]],
) -> pd.DataFrame:
    """Find all cross-tissue gene migration events.

    Parameters
    ----------
    sets_by_tissue : dict
        Keys = tissue name.
        Values = dict mapping bracket label -> list of gene_ids switching
        there. Produced by stamp.io.load_sets_txt.

    Returns
    -------
    DataFrame with columns:
        source_tissue   : tissue where the gene switches first
        source_bracket  : bracket where the gene switches in source_tissue
        target_tissue   : tissue where the gene switches later
        target_bracket  : bracket where the gene switches in target_tissue
        gene_id         : the gene undergoing migration
    """
    # Build a reverse index: gene_id -> list of (tissue, bracket) where it switches
    gene_index: dict[str, list[tuple[str, str]]] = {}
    for tissue, sets in sets_by_tissue.items():
        for bracket, genes in sets.items():
            for gene in genes:
                if gene not in gene_index:
                    gene_index[gene] = []
                gene_index[gene].append((tissue, bracket))

    # Bracket ordering for "i < j" comparison
    bracket_order = {b: i for i, b in enumerate(SWITCHING_BRACKETS)}

    rows = []
    for gene, occurrences in gene_index.items():
        # Only genes that switch in more than one tissue can migrate
        if len(occurrences) < 2:
            continue

        # Sort occurrences by bracket order
        occurrences_sorted = sorted(occurrences, key=lambda x: bracket_order.get(x[1], 999))

        # For every pair (source, target) where source bracket < target bracket
        # and source tissue != target tissue
        for idx_a, (tissue_a, bracket_a) in enumerate(occurrences_sorted):
            for tissue_b, bracket_b in occurrences_sorted[idx_a + 1:]:
                order_a = bracket_order.get(bracket_a, 999)
                order_b = bracket_order.get(bracket_b, 999)
                if order_a >= order_b:
                    continue
                if tissue_a == tissue_b:
                    continue
                rows.append({
                    "source_tissue":  tissue_a,
                    "source_bracket": bracket_a,
                    "target_tissue":  tissue_b,
                    "target_bracket": bracket_b,
                    "gene_id":        gene,
                })

    if not rows:
        return pd.DataFrame(columns=[
            "source_tissue", "source_bracket",
            "target_tissue", "target_bracket", "gene_id"
        ])

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def migration_paths(
    migrations: pd.DataFrame,
    top_n: int = 100,
) -> pd.DataFrame:
    """Aggregate migrations into top-N paths by gene count.

    Parameters
    ----------
    migrations : DataFrame
        Output of find_migrations.
    top_n : int
        Return only the top N paths. Default 100 (same as paper Figure 13).

    Returns
    -------
    DataFrame with columns:
        source_tissue, source_bracket, target_tissue, target_bracket,
        n_genes, genes (space-separated list of gene_ids)
    """
    if migrations.empty:
        return pd.DataFrame(columns=[
            "source_tissue", "source_bracket",
            "target_tissue", "target_bracket",
            "n_genes", "genes"
        ])

    grp = migrations.groupby(
        ["source_tissue", "source_bracket", "target_tissue", "target_bracket"],
        sort=False,
    )
    agg = grp["gene_id"].agg(
        n_genes="count",
        genes=lambda x: " ".join(sorted(x))
    ).reset_index()

    return (
        agg.sort_values("n_genes", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def migration_summary(migrations: pd.DataFrame) -> dict:
    """Return a dict of summary statistics on the migration DataFrame."""
    if migrations.empty:
        return {"n_migration_events": 0, "n_migrating_genes": 0,
                "n_source_tissues": 0, "n_target_tissues": 0}
    return {
        "n_migration_events":  len(migrations),
        "n_migrating_genes":   migrations["gene_id"].nunique(),
        "n_source_tissues":    migrations["source_tissue"].nunique(),
        "n_target_tissues":    migrations["target_tissue"].nunique(),
        "top_source_tissue":   migrations["source_tissue"].value_counts().idxmax(),
        "top_target_tissue":   migrations["target_tissue"].value_counts().idxmax(),
    }