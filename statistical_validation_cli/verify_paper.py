"""Check calculated counts, distributions and gene tables against manuscript references."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def compare_frames(actual: pd.DataFrame, expected: pd.DataFrame, keys: list[str]) -> None:
    if actual.duplicated(keys).any() or expected.duplicated(keys).any():
        raise AssertionError(f"Duplicate result keys: {keys}")
    actual = actual.set_index(keys).sort_index()
    expected = expected.set_index(keys).sort_index()
    pd.testing.assert_index_equal(actual.index, expected.index, check_names=True)
    missing = set(expected.columns) - set(actual.columns)
    if missing:
        raise AssertionError(f"Missing columns: {sorted(missing)}")
    for column in expected:
        a, e = actual[column], expected[column]
        np.testing.assert_array_equal(a.isna(), e.isna(), err_msg=column + " missingness")
        if pd.api.types.is_numeric_dtype(e.dtype):
            # Match counts, seeds and decisions exactly. Allow only tiny floating-point
            # serialization differences in continuous summaries.
            atol = 0.0
            rtol = 1e-12 if pd.api.types.is_float_dtype(e.dtype) else 0.0
            np.testing.assert_allclose(a.to_numpy(dtype=float), e.to_numpy(dtype=float),
                                       rtol=rtol, atol=atol, equal_nan=True, err_msg=column)
        else:
            pd.testing.assert_series_equal(a, e, check_dtype=False, check_names=False)


def gene_digest(frame: pd.DataFrame) -> str:
    """Hash gene-level frequencies using a platform-independent row encoding."""
    columns = ["gene_id", "observed_switching", "bootstrap_switch_frequency",
               "same_bracket_and_direction_frequency"]
    if frame.gene_id.duplicated().any():
        raise AssertionError("Duplicate gene IDs")
    digest = hashlib.sha256()
    for row in frame.sort_values("gene_id")[columns].itertuples(index=False, name=None):
        payload = [str(row[0]), bool(row[1]), float(row[2]).hex(), float(row[3]).hex()]
        digest.update((json.dumps(payload, separators=(",", ":")) + "\n").encode())
    return digest.hexdigest()


def gene_table_digest(frame: pd.DataFrame) -> str:
    """Hash gene IDs and every column of an observed or SMOTE gene table."""
    if frame.gene_id.duplicated().any():
        raise AssertionError("Duplicate gene IDs")
    columns = sorted(frame.columns)
    digest = hashlib.sha256(json.dumps(columns).encode())
    ordered = frame.sort_values("gene_id")[columns]
    for row in ordered.itertuples(index=False, name=None):
        payload = []
        for column, value in zip(columns, row):
            dtype = ordered[column].dtype
            if pd.isna(value):
                payload.append(None)
            elif pd.api.types.is_bool_dtype(dtype):
                payload.append(bool(value))
            elif pd.api.types.is_float_dtype(dtype):
                payload.append(float(value).hex())
            elif pd.api.types.is_integer_dtype(dtype):
                payload.append(int(value))
            else:
                payload.append(str(value))
        digest.update((json.dumps(payload, separators=(",", ":")) + "\n").encode())
    return digest.hexdigest()


def manuscript_claims(directory: Path) -> dict:
    primary = pd.read_csv(directory / "primary.csv")
    bootstrap = pd.read_csv(directory / "bootstrap.csv")
    fixed = pd.read_csv(directory / "smote_fixed.csv")
    adaptive = pd.read_csv(directory / "smote_adaptive.csv")
    primary_set = set(primary.loc[primary.q_value_bh <= 0.05, "tissue"])
    exploratory = set(primary.loc[primary.q_value_bh <= 0.10, "tissue"])
    fixed_05 = set(fixed.loc[fixed.smote_q_value_bh <= 0.05, "tissue"])
    adaptive_05 = set(adaptive.loc[adaptive.smote_q_value_bh <= 0.05, "tissue"])
    fixed_10 = set(fixed.loc[fixed.smote_q_value_bh <= 0.10, "tissue"])
    adaptive_10 = set(adaptive.loc[adaptive.smote_q_value_bh <= 0.10, "tissue"])
    thresholds = pd.read_csv(directory / "threshold_summary.csv")
    return {"n_tissues": len(primary), "primary_q_0_01": int((primary.q_value_bh <= 0.01).sum()),
            "primary_q_0_05": len(primary_set), "primary_q_0_10": len(exploratory),
            "primary_tissues": sorted(primary_set),
            "smote_modified": int(fixed.smote_applied.sum()),
            "smote_unchanged": int((~fixed.smote_applied).sum()),
            "shared_primary_smote_q_0_05": sorted(primary_set & fixed_05 & adaptive_05),
            "shared_exploratory_smote_q_0_10": sorted(exploratory & fixed_10 & adaptive_10),
            "adaptive_only_exploratory_q_0_10": sorted((exploratory & adaptive_10) - fixed_10),
            "unchanged_exploratory": sorted(set(fixed.loc[~fixed.smote_applied, "tissue"]) & exploratory),
            "observed_primary_associations": int(primary.loc[primary.tissue.isin(primary_set), "observed_switching"].sum()),
            "stable_primary_associations_0_80": int(bootstrap.loc[bootstrap.tissue.isin(primary_set), "observed_genes_stable_0_80"].sum()),
            "threshold_medians": [{"tau": float(r.tau), "gene": f"{r.median_jaccard_gene:.2f}",
                                   "gene_bracket": f"{r.median_jaccard_gene_bracket:.2f}"}
                                  for r in thresholds.itertuples()]}


def validate_reference_statistics(reference: Path) -> None:
    """Recalculate empirical p-values, BH q-values and null summaries from saved draws."""
    from .common import safe_filename
    for filename, label, stage, array_folder, observed, p_column, q_column, prefix in (
        ("primary.csv", "paper", "02_primary", "null_counts", "observed_switching", "empirical_p_value", "q_value_bh", "null"),
        ("smote_fixed.csv", "paper-smote-fixed", "03_smote", "null_counts", "smote_test_observed", "smote_empirical_p_value", "smote_q_value_bh", "smote_null"),
        ("smote_adaptive.csv", "paper-smote-adaptive", "03_smote", "null_counts", "smote_test_observed", "smote_empirical_p_value", "smote_q_value_bh", "smote_null"),
    ):
        frame = pd.read_csv(reference / filename)
        eligible = frame["smote_bh_eligible"].astype(bool) if "smote_bh_eligible" in frame else pd.Series(True, index=frame.index)
        pvalues = []
        for row in frame.loc[eligible].to_dict("records"):
            path = reference / "monte_carlo" / label / stage / array_folder / (safe_filename(row["tissue"]) + ".npz")
            with np.load(path) as arrays:
                counts = arrays[arrays.files[0]]
                if counts.size != 9999:
                    raise AssertionError("Expected 9999 null counts: " + str(path))
                # Recalculate independently of the inference helpers.
                p = (1 + int(np.count_nonzero(counts >= row[observed]))) / (counts.size + 1)
                np.testing.assert_allclose(p, row[p_column], rtol=1e-12, atol=0)
                for suffix, value in (("mean", np.mean(counts)), ("median", np.median(counts)),
                                      ("min", np.min(counts)), ("max", np.max(counts))):
                    np.testing.assert_allclose(value, row[prefix + "_" + suffix], rtol=1e-12, atol=0)
                pvalues.append(p)
        order = np.argsort(pvalues, kind="stable")
        adjusted = np.asarray(pvalues)[order] * len(pvalues) / np.arange(1, len(pvalues) + 1)
        adjusted = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
        qvalues = np.empty(len(pvalues))
        qvalues[order] = adjusted
        np.testing.assert_allclose(qvalues, frame.loc[eligible, q_column], rtol=1e-12, atol=0)
        if (~eligible).any() and frame.loc[~eligible, [p_column, q_column]].notna().any().any():
            raise AssertionError("Unmodified SMOTE tissues must not have compound p/q values")
    primary = pd.read_csv(reference / "primary.csv").set_index("tissue")
    table = pd.read_csv(reference / "manuscript_table.csv", dtype={"empirical_p_value": str, "q_value_bh": str})
    if set(table.tissue) != set(primary.loc[primary.q_value_bh <= .10].index):
        raise AssertionError("Published table tissue identities differ")
    for row in table.to_dict("records"):
        expected = primary.loc[row["tissue"]]
        for column in ("observed_switching", "observed_up", "observed_down", "null_median"):
            if expected[column] != row[column]:
                raise AssertionError("Published table mismatch: " + row["tissue"] + " " + column)
        if f"{expected.empirical_p_value:.4f}" != row["empirical_p_value"] or f"{expected.q_value_bh:.3f}" != row["q_value_bh"]:
            raise AssertionError("Published rounded p/q mismatch: " + row["tissue"])


def verify_outputs(output: Path, reference: Path) -> None:
    runs = output / "v10/statistical_validation_cli/runs"
    paths = {"primary.csv": runs / "paper/02_primary/primary_permutation_results.csv",
             "bootstrap.csv": runs / "paper/04_bootstrap/bootstrap_stability_results.csv",
             "smote_fixed.csv": runs / "paper-smote-fixed/03_smote/smote_sensitivity_results.csv",
             "smote_adaptive.csv": runs / "paper-smote-adaptive/03_smote/smote_sensitivity_results.csv",
             "threshold_per_tissue.csv": output / "paper-threshold/per_tissue.csv",
             "threshold_summary.csv": output / "paper-threshold/summary.csv"}
    for filename, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Incomplete reproduction: {path}")
        expected = pd.read_csv(reference / filename)
        actual = pd.read_csv(path)
        keys = ["tau"] if filename == "threshold_summary.csv" else ["tissue"]
        if filename == "threshold_per_tissue.csv":
            keys.append("tau")
        compare_frames(actual, expected, keys)
        print(filename + ": all numerical columns and tissue identities match.", flush=True)
    from .common import safe_filename
    gene_hashes = json.loads((reference / "bootstrap_gene_sha256.json").read_text())
    for tissue, digest in gene_hashes.items():
        path = runs / "paper/04_bootstrap/gene_stability" / (safe_filename(tissue) + ".parquet")
        if gene_digest(pd.read_parquet(path)) != digest:
            raise AssertionError("Bootstrap gene-level frequencies differ: " + tissue)
    gene_tables = json.loads((reference / "observed_and_smote_gene_sha256.json").read_text())
    for relative, digest in gene_tables.items():
        if gene_table_digest(pd.read_parquet(runs / relative)) != digest:
            raise AssertionError("Observed calls or SMOTE gene-level frequencies differ: " + relative)
    for policy, label in (("fixed", "paper-smote-fixed"), ("adaptive", "paper-smote-adaptive")):
        for path in sorted((reference / "seed_runs" / policy).glob("*.csv")):
            actual = runs / label / "03_smote/seed_runs" / path.name
            compare_frames(pd.read_csv(actual), pd.read_csv(path), ["seed_index"])
    # Compare individual Monte Carlo draws as well as summary p-values.
    null_hashes = json.loads((reference / "monte_carlo_sha256.json").read_text())
    for relative, entry in null_hashes.items():
        with np.load(runs / relative) as data:
            for name, digest in entry.items():
                values = np.asarray(data[name], dtype="<i4")
                if hashlib.sha256(values.tobytes()).hexdigest() != digest:
                    raise AssertionError("Monte Carlo sequence differs: " + relative + ":" + name)
    print("All full Monte Carlo sequences, observed gene/bracket/direction calls and bootstrap/SMOTE gene frequencies match.")
