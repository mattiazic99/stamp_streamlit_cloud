"""Dump TPM matrix and sample metadata from Cassandra to Parquet.

Usage
-----
    python scripts/00_dump_cassandra.py --version v10
    python scripts/00_dump_cassandra.py --version v10 --skip-tpm
    python scripts/00_dump_cassandra.py --version v10 --no-resume

Strategy
--------
Cassandra struggles with thousands of single-row PK lookups on a
super-wide table (gene_tpm_full has 19,618 columns). Instead we issue
ONE full-table scan with paginated fetch_size and stream rows out as
they arrive. ~10x faster than gene-by-gene queries on this schema.

Memory notes
------------
- Output dtype: float32. Missing values: NaN (NOT 0.0).
- Peak RAM is bounded by GENE_BATCH_SIZE.

Cassandra column-name quirk (important!)
----------------------------------------
The original GTEx sample IDs contain dashes ('GTEX-1117F-0005-SM-HL9SH').
When stored as Cassandra column names (quoted), the dashes are preserved
in `system_schema.columns`. However, the Python driver replaces dashes
with underscores when exposing rows as named tuples / dicts (because
identifiers in Python attributes cannot contain dashes).

So:
    system_schema column_name : 'GTEX-1117F-0005-SM-HL9SH'  (true name)
    row._asdict() key         : 'GTEX_1117F_0005_SM_HL9SH'  (driver-mangled)

This script keeps both forms and uses each in the right place:
- TRUE names (with dashes) are written as Parquet column headers (so the
  Parquet preserves the original GTEx sample IDs).
- DICT keys (with underscores) are used for lookups inside row._asdict().
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from cassandra.cluster import Cluster, Session
from cassandra.query import SimpleStatement
from tqdm import tqdm

from stamp.config import (
    CASSANDRA_HOSTS,
    CASSANDRA_PORT,
    KEYSPACE,
    METADATA_TABLE,
    TPM_TABLE,
    GtexVersion,
    paths_for,
)


GENE_BATCH_SIZE = 1000
CASSANDRA_FETCH_SIZE = 200


# ---------------------------------------------------------------------------
# Cassandra helpers
# ---------------------------------------------------------------------------

def connect_to_cassandra() -> tuple[Session, Cluster]:
    cluster = Cluster(list(CASSANDRA_HOSTS), port=CASSANDRA_PORT)
    session = cluster.connect(KEYSPACE)
    session.default_timeout = 120
    return session, cluster


def safe_float(x: object) -> np.float32:
    """Cassandra stores TPMs as text. Parse to float32; NaN on failure."""
    if x is None:
        return np.float32("nan")
    try:
        return np.float32(x)
    except (ValueError, TypeError):
        return np.float32("nan")


def get_sample_column_mapping(session: Session) -> tuple[list[str], list[str]]:
    """Return (true_names, dict_keys) for the sample columns.

    true_names: with dashes, as stored in Cassandra and as we want them
                in the Parquet header (= original GTEx sample IDs).
    dict_keys:  with underscores, as exposed by row._asdict() at lookup time.

    Both lists are aligned: dict_keys[i] is the dict key for true_names[i].
    """
    rows = session.execute(
        """
        SELECT column_name FROM system_schema.columns
        WHERE keyspace_name = %s AND table_name = %s
        """,
        [KEYSPACE, TPM_TABLE],
    )
    raw_names = [
        r.column_name for r in rows
        if r.column_name not in ("Name", "Description")
    ]
    raw_names.sort()
    # The Python driver mangles non-identifier characters to '_'.
    # For GTEx ids, only dashes need replacing.
    dict_keys = [n.replace("-", "_") for n in raw_names]
    return raw_names, dict_keys


# ---------------------------------------------------------------------------
# TPM dump (full-scan strategy)
# ---------------------------------------------------------------------------

def dump_tpm(session: Session, version: GtexVersion, resume: bool = True) -> None:
    out_path = paths_for(version)["tpm_parquet"]
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sample_true_names, sample_dict_keys = get_sample_column_mapping(session)
    print(f"  Total sample columns: {len(sample_true_names)}")

    # Resume: detect already-dumped genes
    already_done: set[str] = set()
    partial_path: Path | None = None
    if resume and out_path.exists():
        already_done = _read_existing_gene_ids(out_path)
        if already_done:
            print(f"  Resume: {len(already_done)} genes already dumped, will skip them.")
            partial_path = out_path.with_suffix(".previous.parquet")
            if partial_path.exists():
                partial_path.unlink()
            out_path.rename(partial_path)
            print(f"  Existing file renamed to {partial_path.name} (will merge at end).")
    elif not resume and out_path.exists():
        out_path.unlink()
        print(f"  Removed existing file (--no-resume).")

    schema = _build_tpm_schema(sample_true_names)
    writer = pq.ParquetWriter(out_path, schema, compression="snappy")

    stmt = SimpleStatement(f"SELECT * FROM {TPM_TABLE}", fetch_size=CASSANDRA_FETCH_SIZE)

    skipped = 0
    written = 0
    t0 = time.time()

    batch_names: list[str] = []
    batch_descs: list[str] = []
    batch_values: list[list[np.float32]] = []

    try:
        total_estimate = _try_get_total_count(session)
        pbar = tqdm(
            total=total_estimate,
            desc="Genes",
            unit="gene",
            initial=len(already_done),
        )

        for row in session.execute(stmt):
            d = row._asdict()
            gene_id = d.get("Name")
            if gene_id is None:
                continue

            if gene_id in already_done:
                continue

            description = d.get("Description")
            # IMPORTANT: lookup uses dict keys (with underscores)
            values = [safe_float(d.get(k)) for k in sample_dict_keys]

            if description is None and all(np.isnan(v) for v in values):
                skipped += 1
                pbar.update(1)
                continue

            batch_names.append(gene_id)
            batch_descs.append(description if description is not None else "")
            batch_values.append(values)
            pbar.update(1)

            if len(batch_names) >= GENE_BATCH_SIZE:
                _flush_tpm_batch(writer, schema, batch_names, batch_descs, batch_values)
                written += len(batch_names)
                batch_names.clear()
                batch_descs.clear()
                batch_values.clear()

        if batch_names:
            _flush_tpm_batch(writer, schema, batch_names, batch_descs, batch_values)
            written += len(batch_names)

        pbar.close()
    finally:
        writer.close()

    if partial_path is not None and partial_path.exists():
        print(f"\n  Merging with previous dump ({partial_path.name})...")
        _merge_parquet(out_path, partial_path)
        partial_path.unlink()

    elapsed = time.time() - t0
    size_mb = out_path.stat().st_size / 1e6
    print(f"\n  Done.")
    print(f"  Genes written this run: {written}")
    print(f"  Genes skipped (empty): {skipped}")
    print(f"  Time elapsed: {elapsed/60:.2f} min")
    print(f"  File: {out_path} ({size_mb:.1f} MB)")


def _try_get_total_count(session: Session) -> int | None:
    try:
        rows = list(session.execute("SELECT COUNT(*) FROM gene_mapping_light"))
        return rows[0].count if rows else None
    except Exception:
        return None


def _build_tpm_schema(sample_true_names: list[str]) -> pa.Schema:
    """Pyarrow schema: gene_id, description, then one column per sample."""
    return pa.schema(
        [pa.field("gene_id", pa.string()), pa.field("description", pa.string())]
        + [pa.field(sid, pa.float32()) for sid in sample_true_names]
    )


def _flush_tpm_batch(
    writer: pq.ParquetWriter,
    schema: pa.Schema,
    names: list[str],
    descs: list[str],
    values: list[list[np.float32]],
) -> None:
    arrays: list[pa.Array] = [
        pa.array(names, type=pa.string()),
        pa.array(descs, type=pa.string()),
    ]
    matrix = np.array(values, dtype=np.float32)
    for j in range(matrix.shape[1]):
        arrays.append(pa.array(matrix[:, j], type=pa.float32()))
    table = pa.Table.from_arrays(arrays, schema=schema)
    writer.write_table(table)


def _read_existing_gene_ids(parquet_path: Path) -> set[str]:
    try:
        t = pq.read_table(parquet_path, columns=["gene_id"])
        return set(t.column("gene_id").to_pylist())
    except Exception as e:
        print(f"  Warning: could not read existing Parquet ({e}); starting over.")
        return set()


def _merge_parquet(new_path: Path, previous_path: Path) -> None:
    prev = pq.read_table(previous_path)
    new = pq.read_table(new_path)
    merged = pa.concat_tables([prev, new])
    pq.write_table(merged, new_path, compression="snappy")
    print(f"  Total genes after merge: {merged.num_rows}")


# ---------------------------------------------------------------------------
# Metadata dump (unchanged)
# ---------------------------------------------------------------------------

def dump_metadata(session: Session, version: GtexVersion) -> None:
    out_path = paths_for(version)["metadata_parquet"]
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("  Reading metadata from Cassandra...")
    rows = list(session.execute(f"SELECT * FROM {METADATA_TABLE}"))
    print(f"  Read {len(rows)} rows.")

    base_records: list[dict] = []
    pathology_records: list[dict] = []
    all_pathologies: set[str] = set()

    for row in rows:
        d = row._asdict()
        base = {
            "sample_id": d.get("sample_id"),
            "subject_id": d.get("subject_id"),
            "tissue": d.get("tissue"),
            "tissue_id": d.get("tissue_id"),
            "age_bracket": d.get("age_bracket"),
            "sex": d.get("sex"),
            "rin": d.get("rin"),
            "hardy_scale": d.get("hardy_scale"),
            "ischemic_time": d.get("ischemic_time"),
            "data_type": d.get("data_type"),
        }
        base_records.append(base)

        path_dict: dict[str, bool] = {}
        raw = d.get("pathology_annotations")
        if raw:
            try:
                parsed = json.loads(raw) if isinstance(raw, str) else raw
                for k, v in parsed.items():
                    flag = bool(v)
                    path_dict[k.lower()] = flag
                    all_pathologies.add(k.lower())
            except (json.JSONDecodeError, AttributeError):
                pass
        pathology_records.append(path_dict)

    df = pd.DataFrame(base_records)
    sorted_pathologies = sorted(all_pathologies)
    for path_name in sorted_pathologies:
        col = [rec.get(path_name, False) for rec in pathology_records]
        df[f"pat_{path_name}"] = pd.array(col, dtype="boolean")

    df = df.dropna(subset=["sample_id"]).reset_index(drop=True)

    print(f"  Pathologies detected: {len(sorted_pathologies)}")
    print(f"  Final shape: {df.shape}")

    df.to_parquet(out_path, engine="pyarrow", compression="snappy", index=False)
    size_kb = out_path.stat().st_size / 1024
    print(f"  Written: {out_path} ({size_kb:.1f} KB)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    parser.add_argument("--skip-tpm", action="store_true")
    parser.add_argument("--skip-metadata", action="store_true")
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()

    print(f"Connecting to Cassandra at {CASSANDRA_HOSTS}:{CASSANDRA_PORT}...")
    session, cluster = connect_to_cassandra()

    try:
        if not args.skip_metadata:
            print("\n=== Metadata dump ===")
            dump_metadata(session, args.version)

        if not args.skip_tpm:
            print("\n=== TPM matrix dump (full scan) ===")
            dump_tpm(session, args.version, resume=not args.no_resume)

        print("\nAll done.")
    finally:
        cluster.shutdown()


if __name__ == "__main__":
    main()