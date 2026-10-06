"""Build the STAMP Parquet inputs from the public GTEx release.

Stage 0 of the pipeline, and the only step that touches raw data. It reads the
three open-access files GTEx publishes for a release and writes the two Parquet
files every later stage consumes:

    data/external/{version}/gene_tpm.gct.gz          TPM matrix, GCT (gzipped)
    data/external/{version}/SampleAttributesDS.txt   per-sample annotation
    data/external/{version}/SubjectPhenotypesDS.txt  per-subject phenotypes
        |
        v
    data/parquet/{version}/tpm_matrix.parquet        genes x samples, float32
    data/parquet/{version}/metadata.parquet          sample -> tissue, donor, age

Both GTEx releases publish the same three file types, so the same code builds
v8 and v10. Download them from the GTEx portal Downloads page
(https://gtexportal.org/home/downloads/adult-gtex), put them in
``data/external/{version}/`` under the names above -- rename them if the
release ships longer file names -- then run:

    python scripts/00_build_parquet_from_gtex.py --version v10
    python scripts/00_build_parquet_from_gtex.py --version v8

Only gene-level TPM and the open-access annotations are used. Individual
genotypes, which require a dbGaP application, are not needed anywhere in STAMP.

GCT format
----------
Line 0: "#1.2"
Line 1: "<n_genes>\t<n_samples>"
Line 2: "Name\tDescription\t<sample_id_1>\t<sample_id_2>..."
Line 3+: one gene per row.

Memory strategy
---------------
Reading the whole matrix at once would need several GB. The GCT is streamed in
chunks of GENE_BATCH_SIZE rows, converted to float32 and flushed through a
streaming Parquet writer, so peak memory stays bounded by the batch size
(~100 MB per 1000 genes) whatever the size of the release.

Age brackets
------------
Both releases use the same six: 20-29, 30-39, 40-49, 50-59, 60-69, 70-79.
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

# Support direct execution without installing stamp.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from stamp.config import EXTERNAL_DIR, SUPPORTED_VERSIONS, paths_for

# Config
GENE_BATCH_SIZE = 1000          # genes per Parquet row group
FLOAT_TYPE = pa.float32()       # half the memory of float64

# The three input files, resolved once --version is known.
GCT_FILE: Path
SAMPLE_ATTRS: Path
SUBJECT_PHENO: Path


def set_input_paths(version: str) -> list[Path]:
    """Point the module at data/external/{version}/ and return the input files."""
    global GCT_FILE, SAMPLE_ATTRS, SUBJECT_PHENO
    base = EXTERNAL_DIR / version
    GCT_FILE = base / "gene_tpm.gct.gz"
    SAMPLE_ATTRS = base / "SampleAttributesDS.txt"
    SUBJECT_PHENO = base / "SubjectPhenotypesDS.txt"
    return [GCT_FILE, SAMPLE_ATTRS, SUBJECT_PHENO]


# Helpers

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


# Metadata dump

def dump_metadata(version: str) -> pd.DataFrame:
    """Build and save metadata.parquet for one GTEx release."""
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
# TPM dump (streaming)

def dump_tpm(sample_ids: list[str], n_genes: int, version: str) -> None:
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


# Main

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build data/parquet/{version}/ from the public GTEx files."
    )
    parser.add_argument("--version", choices=list(SUPPORTED_VERSIONS), required=True)
    parser.add_argument("--skip-tpm", action="store_true",
                        help="Skip the (slow) TPM dump; only build metadata.")
    parser.add_argument("--skip-metadata", action="store_true",
                        help="Skip metadata; only build the TPM matrix.")
    args = parser.parse_args()

    inputs = set_input_paths(args.version)
    missing = [f for f in inputs if not f.exists()]
    if missing:
        print(f"[ERROR] Missing input files for GTEx {args.version}:")
        for f in missing:
            print(f"          {f}")
        print("\n        Download the gene-level TPM matrix (GCT) and the two")
        print("        annotation files from the GTEx portal Downloads page:")
        print("          https://gtexportal.org/home/downloads/adult-gtex")
        print(f"        and place them in {EXTERNAL_DIR / args.version}/ as:")
        print("          gene_tpm.gct.gz")
        print("          SampleAttributesDS.txt")
        print("          SubjectPhenotypesDS.txt")
        sys.exit(1)

    print(f"=== Building Parquet inputs for GTEx {args.version} ===")

    if not args.skip_metadata:
        print("\n--- Metadata ---")
        dump_metadata(args.version)

    if not args.skip_tpm:
        print("\n--- TPM matrix ---")
        print("  Reading GCT header...")
        sample_ids, n_genes = read_gct_header(GCT_FILE)
        print(f"  Genes declared in GCT:   {n_genes}")
        print(f"  Samples declared in GCT: {len(sample_ids)}")
        dump_tpm(sample_ids, n_genes, args.version)

    print("\nAll done. Next:")
    print(f"  python scripts/01_normalize.py --version {args.version} --complete-age-bins")


if __name__ == "__main__":
    main()