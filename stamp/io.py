"""Loading and dumping of TPM matrices, metadata, and pipeline outputs.

Two distinct phases:
    1. One-shot ingestion from Cassandra to Parquet (cassandra_to_parquet).
       Done once per GTEx version on a machine that has access to the
       Cassandra cluster. Produces compact Parquet files in float32 with
       gene_id index and original (with-dashes) sample_ids as columns.
    2. Runtime loading from Parquet, used by every step of the pipeline.

Memory considerations
---------------------
- TPM matrix on disk: float32, snappy-compressed Parquet (~1-2 GB).
- Metadata: pyarrow with proper categorical dtypes for tissue/bracket.
- Normalized per-tissue outputs: float32 Parquet, <100 MB each.
- Sets files: plain text, 5 lines, gene_ids space-separated. Tiny.

Streamlit Cloud constraint
--------------------------
The TPM Parquet (1-2 GB) is NOT meant to be deployed on Streamlit Cloud.
Only pre-computed pipeline outputs (sets, Jaccard, PPI, migration) are.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from stamp.config import (
    AGE_BRACKETS,
    SWITCHING_BRACKETS,
    GtexVersion,
    paths_for,
)


# ===========================================================================
# Cassandra -> Parquet (one-shot ingestion)
# ===========================================================================
# Implementation of the Cassandra dump is in scripts/00_dump_cassandra.py
# rather than here, because:
#   - it requires the cassandra-driver and a live Cassandra connection,
#     which we do not want to import or mock in the io module
#   - it is a one-shot operation, not a runtime function
#
# This module exposes only the loaders that the pipeline needs at runtime.
# ===========================================================================


# ===========================================================================
# Parquet loaders (runtime)
# ===========================================================================

def load_tpm_matrix(
    version: GtexVersion,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """Load the TPM matrix from Parquet.

    Parameters
    ----------
    version : "v8" or "v10"
        GTEx release.
    columns : list of str, optional
        If provided, load only these sample columns (plus the gene index).
        Useful to load a single tissue without pulling all 19k samples
        into memory.

    Returns
    -------
    DataFrame
        Index = gene_id (str), columns = sample_id (str), dtype float32.

    Raises
    ------
    FileNotFoundError
        If the Parquet file does not exist; user must run the dump script.
    """
    path = paths_for(version)["tpm_parquet"]
    if not path.exists():
        raise FileNotFoundError(
            f"TPM Parquet for {version} not found at {path}. "
            f"Run scripts/00_dump_cassandra.py --version {version} first."
        )

    if columns is not None:
        # Always include the gene_id column; user provides sample columns.
        cols_to_read = ["gene_id", *columns]
        df = pd.read_parquet(path, columns=cols_to_read)
    else:
        df = pd.read_parquet(path)

    df = df.set_index("gene_id")
    return df


def load_metadata(version: GtexVersion) -> pd.DataFrame:
    """Load sample metadata from Parquet.

    Returns
    -------
    DataFrame with columns at minimum: sample_id (str), tissue (category),
    age_bracket (category), subject_id (str), sex (category).
    """
    path = paths_for(version)["metadata_parquet"]
    if not path.exists():
        raise FileNotFoundError(
            f"Metadata Parquet for {version} not found at {path}. "
            f"Run scripts/00_dump_cassandra.py --version {version} first."
        )
    return pd.read_parquet(path)
    
def load_normalized_tissue(
    version: GtexVersion, tissue: str, complete: bool = False
) -> pd.DataFrame:
    """Load the normalized matrix for a single tissue.

    Returns DataFrame with rows = gene_id, columns = age brackets, float32.
    """
    path = paths_for(version, complete)["normalized"] / f"{_safe_filename(tissue)}.parquet"
    if not path.exists():
        raise FileNotFoundError(
            f"Normalized Parquet for tissue '{tissue}' ({version}) not found "
            f"at {path}. Run scripts/01_normalize.py --version {version} first."
        )
    df = pd.read_parquet(path)
    if "gene_id" in df.columns:
        df = df.set_index("gene_id")
    return df

def load_sets_txt(
    version: GtexVersion, tissue: str, complete: bool = False
) -> dict[str, list[str]]:
    """Load the 5-line switching gene sets file for a tissue.

    Returns
    -------
    dict mapping bracket label -> list of gene_ids that switch there.
    Empty lines map to empty lists. Keys are SWITCHING_BRACKETS in order.
    """
    path = paths_for(version, complete)["sets"] / f"{_safe_filename(tissue)}_sets.txt"
    if not path.exists():
        raise FileNotFoundError(
            f"Sets file for tissue '{tissue}' ({version}) not found at {path}."
        )

    with path.open("r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    if len(lines) != len(SWITCHING_BRACKETS):
        raise ValueError(
            f"{path} has {len(lines)} lines; expected {len(SWITCHING_BRACKETS)} "
            f"(one per switching bracket: {SWITCHING_BRACKETS})."
        )

    return {
        bracket: line.split() if line else []
        for bracket, line in zip(SWITCHING_BRACKETS, lines)
    }


# ===========================================================================
# Parquet writers (runtime)
# ===========================================================================

def save_normalized_tissue(
    version: GtexVersion,
    tissue: str,
    normalized_df: pd.DataFrame,
    complete: bool = False,
) -> Path:
    """Save a tissue's normalized matrix to Parquet.

    The output is float32 with snappy compression. The gene_id index is
    materialized as a column on disk for compatibility with all readers.

    Returns
    -------
    Path of the written file.
    """
    out_dir = paths_for(version, complete)["normalized"]
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{_safe_filename(tissue)}.parquet"

    df_to_write = normalized_df.reset_index()
    # Ensure the index column is named 'gene_id' on disk
    df_to_write = df_to_write.rename(columns={df_to_write.columns[0]: "gene_id"})

    df_to_write.to_parquet(
        path,
        engine="pyarrow",
        compression="snappy",
        index=False,
    )
    return path


def save_sets_txt(
    version: GtexVersion,
    tissue: str,
    sets_by_bracket: dict[str, list[str]],
    complete: bool = False,
) -> Path:
    """Save switching gene sets to a text file with exactly 5 lines.

    The file format is one line per SWITCHING_BRACKET, in chronological
    order:
        line 1 -> 30-39
        line 2 -> 40-49
        line 3 -> 50-59
        line 4 -> 60-69
        line 5 -> 70-79
    Each line contains the space-separated gene_ids switching at that
    bracket. Empty brackets produce empty lines.

    Parameters
    ----------
    sets_by_bracket : dict
        Keys must be a subset of SWITCHING_BRACKETS. Missing brackets
        yield empty lines. Extra (unknown) keys raise ValueError.

    Returns
    -------
    Path of the written file.
    """
    unknown = set(sets_by_bracket) - set(SWITCHING_BRACKETS)
    if unknown:
        raise ValueError(
            f"sets_by_bracket has unknown brackets: {sorted(unknown)}. "
            f"Allowed: {SWITCHING_BRACKETS}"
        )

    out_dir = paths_for(version, complete)["sets"]
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{_safe_filename(tissue)}_sets.txt"

    lines = [
        " ".join(sets_by_bracket.get(bracket, []))
        for bracket in SWITCHING_BRACKETS
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ===========================================================================
# Helpers
# ===========================================================================

def _safe_filename(tissue: str) -> str:
    """Convert a tissue display name into a filesystem-safe filename.

    GTEx tissues contain spaces, dashes, parentheses, and slashes. We map
    them all to underscores. This matches the convention used by the R
    pipeline outputs so files are interchangeable.

    Examples
    --------
    "Brain - Substantia nigra" -> "Brain_Substantia_nigra"
    "Cells - EBV-transformed lymphocytes" -> "Cells_EBV_transformed_lymphocytes"
    """
    bad_chars = " -()/,"
    out = tissue
    for c in bad_chars:
        out = out.replace(c, "_")
    # Collapse runs of underscores
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def list_available_tissues(version: GtexVersion, complete: bool = False) -> list[str]:
    """List tissues that have a normalized Parquet output already."""
    out_dir = paths_for(version, complete)["normalized"]
    if not out_dir.exists():
        return []
    return sorted(p.stem for p in out_dir.glob("*.parquet"))