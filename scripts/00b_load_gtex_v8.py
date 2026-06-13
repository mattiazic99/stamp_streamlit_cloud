"""Load GTEx v8 local files and produce Parquet outputs for STAMP pipeline.

This script is the v8 equivalent of 00_dump_cassandra.py (which reads v10
from Cassandra). Both scripts produce the same Parquet schema so that the
rest of the pipeline (01_normalize, 02_switching, 03_jaccard, ...) can run
identically on both versions with just --version v8 / --version v10.

Input files (must exist in data/external/v8/):
    gene_tpm.gct.gz         TPM matrix in GCT format (gzip-compressed)
    SampleAttributesDS.txt  Per-sample metadata (tissue, etc.)
    SubjectPhenotypesDS.txt Per-subject phenotypes (age bracket, sex)

Output files:
    data/parquet/v8/tpm_matrix.parquet
    data/parquet/v8/metadata.parquet

GCT format overview
-------------------
Line 0: "#1.2"
Line 1: "<n_genes>\t<n_samples>"
Line 2: header row: "Name\tDescription\t<sample_id_1>\t<sample_id_2>..."
Lines 3+: one gene per row: "<gene_id>\t<symbol>\t<tpm_1>\t<tpm_2>..."

Memory strategy
---------------
Reading 56200 x 17382 float values at once would require ~4 GB RAM.
We read the GCT in chunks of GENE_BATCH_SIZE lines, convert each chunk
to float32, and flush to Parquet via streaming writer. Peak RAM stays
bounded by the batch size (~100 MB per batch of 1000 genes).

Age bracket mapping
-------------------
GTEx v8 uses the same 6 age brackets as v10:
    20-29, 30-39, 40-49, 50-59, 60-69, 70-79
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

from stamp.config import paths_for

# -----------------------------------------------------------------
# Config
# -----------------------------------------------------------------
V8_DIR = Path("data/external/v8")
GCT_FILE = V8_DIR / "gene_tpm.gct.gz"
SAMPLE_ATTRS = V8_DIR / "SampleAttributesDS.txt"
SUBJECT_PHENO = V8_DIR / "SubjectPhenotypesDS.txt"

GENE_BATCH_SIZE = 1000          # genes per Parquet row group
FLOAT_TYPE = pa.float32()       # half the memory of float64


# -----------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------

def extract_subject_id(sample_id: str) -> str:
    """'GTEX-1117F-0226-SM-5GZZ7' -> 'GTEX-1117F'."""
    parts = sample_id.split("-")
    return "-".join(parts[:2])


def load_subject_age(pheno_path: Path) -> dict[str, str]:
    """Return {subject_id: age_bracket}."""
    df = pd.read_csv(pheno_path, sep="\t", usecols=["SUBJID", "AGE"])
    return dict(zip(df["SUBJID"], df["AGE"]))


def load_subject_sex(pheno_path: Path) -> dict[str, str]:
    """Return {subject_id: 'male'|'female'}."""
    df = pd.read_csv(pheno_path, sep="\t", usecols=["SUBJID", "SEX"])
    sex_map = {1: "male", 2: "female"}
    return {row.SUBJID: sex_map.get(row.SEX, "unknown") for row in df.itertuples()}


def load_sample_attrs(attrs_path: Path) -> pd.DataFrame:
    """Return sample attributes with columns: sample_id, tissue, tissue_id."""
    df = pd.read_csv(attrs_path, sep="\t", usecols=["SAMPID", "SMTS", "SMTSD"])
    df = df.rename(columns={
        "SAMPID": "sample_id",
        "SMTS":   "tissue",        # broad tissue (e.g. "Brain")
        "SMTSD":  "tissue_id",     # specific tissue (e.g. "Brain - Cortex")
    })
    return df.dropna(subset=["sample_id"])


def read_gct_header(gct_path: Path) -> tuple[list[str], int]:
    """Parse GCT header; return (sample_ids, n_genes).

    Reads only the first 3 lines (fast), returns the list of sample IDs
    from the column headers and the declared gene count.
    """
    with gzip.open(gct_path, "rt", encoding="utf-8") as f:
        f.readline()                          # "#1.2"
        dims = f.readline().strip().split()   # "56200\t17382"
        n_genes = int(dims[0])
        header = f.readline().strip().split("\t")
    # header[0] = "Name", header[1] = "Description", header[2:] = sample IDs
    sample_ids = header[2:]
    return sample_ids, n_genes


# -----------------------------------------------------------------
# Metadata dump
# -----------------------------------------------------------------

def dump_metadata(version: str = "v8") -> None:
    """Build and save metadata.parquet for v8."""
    print("  Loading sample attributes...")
    # SMTSD = specific tissue display name (e.g. "Brain - Cortex")
    # SMTS  = broad tissue (e.g. "Brain") - not used downstream
    attrs = pd.read_csv(
        SAMPLE_ATTRS,
        sep="\t",
        usecols=["SAMPID", "SMTSD"],
    ).rename(columns={"SAMPID": "sample_id", "SMTSD": "tissue"})
    attrs = attrs.dropna(subset=["sample_id", "tissue"])

    print("  Loading subject phenotypes...")
    age_map = load_subject_age(SUBJECT_PHENO)
    sex_map = load_subject_sex(SUBJECT_PHENO)

    print("  Reading GCT header to get sample list...")
    sample_ids_in_tpm, _ = read_gct_header(GCT_FILE)
    sample_ids_set = set(sample_ids_in_tpm)

    # Keep only samples present in the TPM matrix
    attrs = attrs[attrs["sample_id"].isin(sample_ids_set)].copy()
    print(f"  Samples in TPM matrix:         {len(sample_ids_in_tpm)}")
    print(f"  Samples matched to attrs:      {len(attrs)}")

    # Derive subject_id from sample_id
    attrs["subject_id"] = attrs["sample_id"].apply(extract_subject_id)

    # Add age_bracket and sex from subject phenotypes
    attrs["age_bracket"] = attrs["subject_id"].map(age_map)
    attrs["sex"] = attrs["subject_id"].map(sex_map)

    # tissue_id: filesystem-safe slug of the tissue display name
    attrs["tissue_id"] = (
        attrs["tissue"]
        .str.replace(" - ", "_", regex=False)
        .str.replace(" ", "_", regex=False)
        .str.replace("(", "", regex=False)
        .str.replace(")", "", regex=False)
    )

    # Drop samples without age (some donors may not have phenotype data)
    n_before = len(attrs)
    meta = attrs[["sample_id", "subject_id", "tissue", "tissue_id",
                  "age_bracket", "sex"]].dropna(subset=["age_bracket"]).reset_index(drop=True)
    n_dropped = n_before - len(meta)
    if n_dropped > 0:
        print(f"  Dropped {n_dropped} samples with missing age bracket.")

    print(f"  Final metadata: {len(meta)} samples, "
          f"{meta['tissue'].nunique()} tissues, "
          f"{meta['age_bracket'].nunique()} age brackets")
    print(f"  Age distribution:\n{meta['age_bracket'].value_counts().sort_index().to_string()}")

    out_path = paths_for(version)["metadata_parquet"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    meta.to_parquet(out_path, index=False, compression="snappy")
    size_kb = out_path.stat().st_size / 1024
    print(f"  Written: {out_path} ({size_kb:.1f} KB)")
    return meta
# -----------------------------------------------------------------
# TPM dump (streaming)
# -----------------------------------------------------------------

def dump_tpm(sample_ids: list[str], n_genes: int, version: str = "v8") -> None:
    """Stream GCT → Parquet in batches of GENE_BATCH_SIZE."""
    out_path = paths_for(version)["tpm_parquet"]
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists():
        print(f"  [skip] {out_path.name} already exists. Delete to regenerate.")
        return

    # Build pyarrow schema: gene_id (str), description (str), then one
    # float32 column per sample_id (using original dash-separated IDs).
    schema = pa.schema(
        [pa.field("gene_id", pa.string()), pa.field("description", pa.string())]
        + [pa.field(sid, FLOAT_TYPE) for sid in sample_ids]
    )

    writer = pq.ParquetWriter(out_path, schema, compression="snappy")

    batch_names: list[str] = []
    batch_descs: list[str] = []
    batch_values: list[np.ndarray] = []  # each item: float32 array of len(sample_ids)

    skipped = 0
    t0 = time.time()

    with gzip.open(GCT_FILE, "rt", encoding="utf-8") as f:
        f.readline()   # "#1.2"
        f.readline()   # dimensions
        f.readline()   # header

        pbar = tqdm(total=n_genes, desc="Genes", unit="gene")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                skipped += 1
                pbar.update(1)
                continue

            gene_id = parts[0]
            description = parts[1]
            tpm_values = np.array(parts[2:], dtype=np.float32)

            # Safety: some rows may have fewer values than expected
            if len(tpm_values) != len(sample_ids):
                skipped += 1
                pbar.update(1)
                continue

            batch_names.append(gene_id)
            batch_descs.append(description)
            batch_values.append(tpm_values)
            pbar.update(1)

            if len(batch_names) >= GENE_BATCH_SIZE:
                _flush_batch(writer, schema, batch_names, batch_descs,
                             batch_values, sample_ids)
                batch_names.clear()
                batch_descs.clear()
                batch_values.clear()

        # Final flush
        if batch_names:
            _flush_batch(writer, schema, batch_names, batch_descs,
                         batch_values, sample_ids)

        pbar.close()

    writer.close()
    elapsed = time.time() - t0
    size_mb = out_path.stat().st_size / 1e6
    print(f"\n  Genes skipped (malformed rows): {skipped}")
    print(f"  Time: {elapsed/60:.1f} min")
    print(f"  Written: {out_path} ({size_mb:.1f} MB)")


def _flush_batch(
    writer: pq.ParquetWriter,
    schema: pa.Schema,
    names: list[str],
    descs: list[str],
    values: list[np.ndarray],
    sample_ids: list[str],
) -> None:
    """Write one row group to the Parquet file."""
    arrays: list[pa.Array] = [
        pa.array(names, type=pa.string()),
        pa.array(descs, type=pa.string()),
    ]
    matrix = np.stack(values, axis=0)  # (n_genes_in_batch, n_samples)
    for j in range(matrix.shape[1]):
        arrays.append(pa.array(matrix[:, j], type=FLOAT_TYPE))
    table = pa.Table.from_arrays(arrays, schema=schema)
    writer.write_table(table)


# -----------------------------------------------------------------
# Main
# -----------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-tpm", action="store_true",
                        help="Skip the (slow) TPM dump; only do metadata.")
    parser.add_argument("--skip-metadata", action="store_true",
                        help="Skip metadata dump; only do TPM.")
    args = parser.parse_args()

    # Verify files exist
    for f in [GCT_FILE, SAMPLE_ATTRS, SUBJECT_PHENO]:
        if not f.exists():
            print(f"[ERROR] Missing file: {f}")
            print("Run download_gtex_v8.py first.")
            sys.exit(1)

    VERSION = "v8"
    print(f"=== GTEx v8 loading ===")

    if not args.skip_metadata:
        print("\n--- Metadata ---")
        dump_metadata(VERSION)

    if not args.skip_tpm:
        print("\n--- TPM matrix ---")
        print("  Reading GCT header...")
        sample_ids, n_genes = read_gct_header(GCT_FILE)
        print(f"  Genes declared in GCT: {n_genes}")
        print(f"  Samples declared in GCT: {len(sample_ids)}")
        dump_tpm(sample_ids, n_genes, VERSION)

    print("\nAll done.")


if __name__ == "__main__":
    main()