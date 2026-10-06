"""Explore sensitivity to SMOTE age-bin balancing.

These results are separate from the primary real-donor analysis. The optional
compound permutation test reruns SMOTE inside every permutation before applying
STAMP, rather than testing a single fixed synthetic dataset.
"""

from __future__ import annotations

import os
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from stamp.config import AGE_BRACKETS
from stamp.io import load_metadata

from .common import (
    age_means_from_codes,
    base_parser,
    bh_qvalues,
    bin_count_columns,
    config_matches,
    checkpoint_provenance,
    record_checkpoint_execution,
    empirical_right_tail_p,
    enforce_final_scope,
    load_tissue_data,
    parse_tissue_arguments,
    read_json,
    safe_filename,
    select_tissues,
    stable_seed,
    summary_statistics,
    switching_state,
    validation_root,
    write_json,
)


def _load_smote_class():
    """Load SMOTE and report incompatible or missing dependencies."""
    os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")
    try:
        from imblearn.over_sampling import SMOTE
    except Exception as exc:  # An internal dependency can raise ImportError during the import.
        try:
            installed = package_version("imbalanced-learn")
        except PackageNotFoundError:
            installed = "not installed"
        raise RuntimeError(
            "SMOTE is unavailable. Installed imbalanced-learn version: "
            f"{installed}. Create the validation virtual environment and install "
            "`statistical_validation_cli/requirements.txt`. Original error: "
            f"{exc}"
        ) from exc
    return SMOTE


def _dependency_versions() -> tuple[str, str]:
    """Record package versions for the SMOTE run."""
    try:
        imblearn_version = package_version("imbalanced-learn")
    except PackageNotFoundError:
        imblearn_version = "not installed"
    try:
        sklearn_version = package_version("scikit-learn")
    except PackageNotFoundError:
        sklearn_version = "not installed"
    return imblearn_version, sklearn_version

def _strategy_and_neighbors(
    age_counts: np.ndarray,
    target_n: int,
    target_policy: str = "fixed-floor",
) -> tuple[dict[int, int], int | None, list[str]]:
    """Create a prespecified per-tissue SMOTE sampling rule.

    ``fixed-floor`` raises every bin with 2..target_n-1 real donors to
    ``target_n``.

    ``adaptive-second-smallest`` first orders the six bin sizes. If only the
    smallest bin is below ``target_n``, that bin is raised to the
    second-smallest real bin size. If at least two bins are below
    ``target_n``, all bins with 2..target_n-1 donors are raised to
    ``target_n``.

    Under both policies, singleton bins remain untouched because SMOTE
    interpolation is mathematically impossible with one real donor.
    """
    counts = np.asarray(age_counts, dtype=np.int32)
    allowed_policies = {"fixed-floor", "adaptive-second-smallest"}
    if target_policy not in allowed_policies:
        raise ValueError(
            f"Unknown SMOTE target policy {target_policy!r}; "
            f"choose one of {sorted(allowed_policies)}"
        )

    singleton_bins = [
        AGE_BRACKETS[code]
        for code, count in enumerate(counts)
        if int(count) == 1
    ]

    if target_policy == "fixed-floor":
        strategy = {
            code: int(target_n)
            for code, count in enumerate(counts)
            if 2 <= int(count) < target_n
        }
    else:
        ordered_counts = np.sort(counts)
        minimum_n = int(ordered_counts[0])
        second_minimum_n = int(ordered_counts[1])
        if minimum_n >= target_n:
            strategy = {}
        elif second_minimum_n < target_n:
            strategy = {
                code: int(target_n)
                for code, count in enumerate(counts)
                if 2 <= int(count) < target_n
            }
        elif minimum_n >= 2:
            strategy = {
                code: second_minimum_n
                for code, count in enumerate(counts)
                if int(count) == minimum_n
            }
        else:
            strategy = {}

    if not strategy:
        return strategy, None, singleton_bins
    minimum_targeted_n = min(int(counts[code]) for code in strategy)
    k_neighbors = min(5, minimum_targeted_n - 1)
    return strategy, k_neighbors, singleton_bins


def _smote_once(
    normalized_values: np.ndarray,
    age_codes: np.ndarray,
    strategy: dict[int, int],
    k_neighbors: int | None,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return an augmented gene-by-sample matrix and age codes."""
    if not strategy:
        return normalized_values, age_codes
    SMOTE = _load_smote_class()
    sampler = SMOTE(
        sampling_strategy=strategy,
        k_neighbors=k_neighbors,
        random_state=random_state,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Could not find the number of physical cores.*",
            category=UserWarning,
        )
        samples_by_genes, resampled_codes = sampler.fit_resample(
            normalized_values.T,
            age_codes,
        )
    return (
        np.asarray(samples_by_genes, dtype=np.float32).T,
        np.asarray(resampled_codes, dtype=np.int8),
    )


def _jaccard(mask_a: np.ndarray, mask_b: np.ndarray) -> float:
    union = int(np.count_nonzero(mask_a | mask_b))
    if union == 0:
        return 1.0
    return float(np.count_nonzero(mask_a & mask_b) / union)


def _process_tissue(job: dict[str, Any]) -> dict[str, Any]:
    tissue = job["tissue"]
    out_dir = Path(job["out_dir"])
    stem = safe_filename(tissue)
    checkpoint = out_dir / "checkpoints" / f"{stem}.json"
    seed_path = out_dir / "seed_runs" / f"{stem}.csv"
    gene_path = out_dir / "gene_stability" / f"{stem}.parquet"
    null_path = out_dir / "null_counts" / f"{stem}.npz"
    config = {
        "provenance": job.get("provenance"),
        "version": job["version"],
        "tau": job["tau"],
        "epsilon": job["epsilon"],
        "seed": job["seed"],
        "n_seeds": job["n_seeds"],
        "target_n": job["target_n"],
        "target_policy": job["target_policy"],
        "n_permutations": job["n_permutations"],
    }
    required_paths = [checkpoint, seed_path, gene_path]
    if not job["force"] and all(path.exists() for path in required_paths):
        cached = read_json(checkpoint)
        if config_matches(cached, config):
            result = cached["result"]
            compound_complete = not (
                job["n_permutations"] > 0
                and result.get("smote_applied", False)
                and not null_path.exists()
            )
            if compound_complete:
                record_checkpoint_execution(out_dir, stem, config, "reused")
                return result

    data = load_tissue_data(job["version"], tissue, epsilon=job["epsilon"])
    real_means = age_means_from_codes(data.normalized_values, data.age_codes)
    real_state = switching_state(real_means, tau=job["tau"])
    strategy, k_neighbors, singleton_bins = _strategy_and_neighbors(
        data.age_counts,
        job["target_n"],
        job["target_policy"],
    )

    frequency = np.zeros(data.gene_ids.size, dtype=np.int32)
    exact_real_call = np.zeros(data.gene_ids.size, dtype=np.int32)
    seed_rows: list[dict[str, Any]] = []
    for seed_index in range(job["n_seeds"]):
        current_seed = stable_seed(
            job["seed"],
            "smote_seed_run",
            tissue,
            seed_index,
        )
        augmented, resampled_codes = _smote_once(
            data.normalized_values,
            data.age_codes,
            strategy,
            k_neighbors,
            current_seed,
        )
        smote_state = switching_state(
            age_means_from_codes(augmented, resampled_codes),
            tau=job["tau"],
        )
        frequency += smote_state.mask.astype(np.int32)
        exact_real_call += (
            smote_state.mask
            & real_state.mask
            & (smote_state.position == real_state.position)
            & (smote_state.direction == real_state.direction)
        ).astype(np.int32)
        seed_rows.append(
            {
                "tissue": tissue,
                "seed_index": seed_index,
                "random_state": current_seed,
                "real_switching": real_state.count,
                "smote_switching": smote_state.count,
                "delta_switching": smote_state.count - real_state.count,
                "ratio_switching": (
                    smote_state.count / real_state.count
                    if real_state.count > 0
                    else None
                ),
                "jaccard_real_vs_smote": _jaccard(
                    real_state.mask,
                    smote_state.mask,
                ),
                "n_samples_real": int(data.sample_ids.size),
                "n_samples_after_smote": int(augmented.shape[1]),
                "n_synthetic": int(augmented.shape[1] - data.sample_ids.size),
                "k_neighbors": k_neighbors,
            }
        )

    seed_frame = pd.DataFrame(seed_rows)
    seed_path.parent.mkdir(parents=True, exist_ok=True)
    seed_frame.to_csv(seed_path, index=False)

    stability_frame = pd.DataFrame(
        {
            "gene_id": data.gene_ids,
            "real_switching": real_state.mask,
            "smote_switch_frequency": frequency / job["n_seeds"],
            "same_bracket_and_direction_frequency": exact_real_call / job["n_seeds"],
        }
    )
    stability_frame = stability_frame.loc[
        stability_frame["real_switching"]
        | (stability_frame["smote_switch_frequency"] > 0)
    ]
    gene_path.parent.mkdir(parents=True, exist_ok=True)
    stability_frame.to_parquet(gene_path, index=False)

    smote_counts = seed_frame["smote_switching"].to_numpy(dtype=np.int32)
    result: dict[str, Any] = {
        "version": job["version"],
        "tissue": tissue,
        "tau": job["tau"],
        "epsilon": job["epsilon"],
        "n_samples_real": int(data.sample_ids.size),
        "n_genes_after_epsilon": int(data.gene_ids.size),
        "n_min_age_bin": int(np.min(data.age_counts)),
        "real_switching": real_state.count,
        "smote_applied": bool(strategy),
        "smote_target_policy": job["target_policy"],
        "smote_target_n": int(job["target_n"]),
        "smote_targeted_bins": ",".join(
            AGE_BRACKETS[code] for code in strategy
        ),
        "smote_sampling_strategy": ",".join(
            f"{AGE_BRACKETS[code]}:{strategy[code]}"
            for code in sorted(strategy)
        ),
        "singleton_bins_not_oversampled": ",".join(singleton_bins),
        "k_neighbors": k_neighbors,
        "n_smote_seeds": int(job["n_seeds"]),
        "jaccard_mean": float(seed_frame["jaccard_real_vs_smote"].mean()),
        "jaccard_min": float(seed_frame["jaccard_real_vs_smote"].min()),
        "jaccard_max": float(seed_frame["jaccard_real_vs_smote"].max()),
        "smote_empirical_p_value": None,
        "smote_test_observed": None,
        "n_smote_permutations": int(
            job["n_permutations"] if strategy else 0
        ),
        **bin_count_columns(data.age_counts),
        **summary_statistics(smote_counts, "smote"),
    }

    if job["n_permutations"] > 0 and strategy:
        test_seed = stable_seed(job["seed"], "smote_compound_test", tissue)
        observed_augmented, observed_codes = _smote_once(
            data.normalized_values,
            data.age_codes,
            strategy,
            k_neighbors,
            test_seed,
        )
        observed_test_state = switching_state(
            age_means_from_codes(observed_augmented, observed_codes),
            tau=job["tau"],
        )

        rng = np.random.default_rng(
            stable_seed(job["seed"], "smote_label_permutations", tissue)
        )
        null_counts = np.empty(job["n_permutations"], dtype=np.int32)
        for permutation_index in range(job["n_permutations"]):
            permuted_codes = rng.permutation(data.age_codes)
            permuted_augmented, permuted_resampled_codes = _smote_once(
                data.normalized_values,
                permuted_codes,
                strategy,
                k_neighbors,
                test_seed,
            )
            null_counts[permutation_index] = switching_state(
                age_means_from_codes(
                    permuted_augmented,
                    permuted_resampled_codes,
                ),
                tau=job["tau"],
            ).count

        null_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(null_path, null_counts=null_counts)
        result["smote_test_observed"] = observed_test_state.count
        result["smote_empirical_p_value"] = empirical_right_tail_p(
            observed_test_state.count,
            null_counts,
        )
        result.update(summary_statistics(null_counts, "smote_null"))

    write_json(checkpoint, {"config": config, "result": result})
    record_checkpoint_execution(out_dir, stem, config, "recomputed")
    return result


def run(args) -> Path:
    if args.seeds < 1:
        raise ValueError("--seeds must be at least 1")
    if args.run_label == "final" and args.seeds < 20:
        raise ValueError(
            "The protected final SMOTE sensitivity requires at least 20 "
            "seeds. Use a non-final run label for shorter runs."
        )
    if args.target_n < 2:
        raise ValueError("--target-n must be at least 2")
    if args.permutations < 0:
        raise ValueError("--permutations cannot be negative")
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")

    _load_smote_class()
    imblearn_version, sklearn_version = _dependency_versions()
    print(
        "SMOTE environment: imbalanced-learn "
        f"{imblearn_version}; scikit-learn {sklearn_version}"
    )
    metadata = load_metadata(args.version)
    all_complete = select_tissues(metadata)
    tissues = select_tissues(
        metadata,
        requested=parse_tissue_arguments(args.tissue),
    )
    enforce_final_scope(tissues, all_complete, args.run_label, "SMOTE")
    if (
        args.permutations > 1_000
        and len(tissues) > 4
        and not bool(getattr(args, "allow_expensive_test", False))
    ):
        raise RuntimeError(
            "The compound SMOTE+STAMP test is intentionally limited to at "
            "most 4 tissues when permutations exceed 1000. It is exploratory "
            "and computationally very expensive. Reduce the scope/permutations "
            "or pass --allow-expensive-test after explicitly accepting the cost."
        )
    out_dir = validation_root(
        args.version, args.output_root, args.run_label
    ) / "03_smote"
    out_dir.mkdir(parents=True, exist_ok=True)
    provenance = checkpoint_provenance(args.version)
    jobs = [
        {
            "provenance": provenance,
            "version": args.version,
            "tissue": tissue,
            "tau": args.tau,
            "epsilon": args.epsilon,
            "seed": args.seed,
            "n_seeds": args.seeds,
            "target_n": args.target_n,
            "target_policy": args.target_policy,
            "n_permutations": args.permutations,
            "force": args.force,
            "out_dir": str(out_dir),
        }
        for tissue in tissues
    ]

    results: list[dict[str, Any]] = []
    if args.jobs == 1:
        for index, job in enumerate(jobs, start=1):
            print(f"[{index}/{len(jobs)}] SMOTE sensitivity: {job['tissue']}")
            results.append(_process_tissue(job))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(_process_tissue, job): job["tissue"]
                for job in jobs
            }
            completed = 0
            for future in as_completed(futures):
                completed += 1
                tissue = futures[future]
                results.append(future.result())
                print(f"[{completed}/{len(jobs)}] Completed: {tissue}")

    frame = pd.DataFrame(results).sort_values("tissue", kind="stable")
    frame["smote_q_value_bh"] = np.nan
    frame["smote_significant_bh_0_05"] = False
    frame["smote_bh_eligible"] = False
    frame["smote_bh_family_size"] = (
        0 if args.permutations == 0 else np.nan
    )
    if args.permutations > 0:
        eligible = frame["smote_applied"].astype(bool)
        frame.loc[~eligible, "smote_empirical_p_value"] = np.nan
        frame.loc[eligible, "smote_bh_eligible"] = True
        if eligible.any():
            q_values = bh_qvalues(
                frame.loc[eligible, "smote_empirical_p_value"]
            )
            frame.loc[eligible, "smote_q_value_bh"] = q_values
            frame.loc[eligible, "smote_significant_bh_0_05"] = (
                q_values <= 0.05
            )
            frame.loc[eligible, "smote_bh_family_size"] = int(
                eligible.sum()
            )

    result_path = out_dir / "smote_sensitivity_results.csv"
    frame.to_csv(result_path, index=False)
    print(f"SMOTE sensitivity results written to: {result_path}")
    if args.permutations == 0:
        print(
            "No SMOTE p/q values were requested. This output is descriptive "
            "(seed stability, switching deltas and Jaccard overlap)."
        )
    return result_path


def build_parser():
    parser = base_parser(
        "Run exploratory SMOTE sensitivity and an optional compound permutation test."
    )
    parser.add_argument(
        "--seeds",
        type=int,
        default=20,
        help="Independent SMOTE random seeds per tissue.",
    )
    parser.add_argument(
        "--target-n",
        type=int,
        default=20,
        help=(
            "Minimum donor count used by the selected target policy "
            "(default: 20)."
        ),
    )
    parser.add_argument(
        "--target-policy",
        choices=["fixed-floor", "adaptive-second-smallest"],
        default="fixed-floor",
        help=(
            "fixed-floor raises every 2..target_n-1 bin to target_n; "
            "adaptive-second-smallest raises an isolated smallest bin to the "
            "second-smallest size, otherwise applies the fixed floor."
        ),
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=0,
        help=(
            "Optional permutations for the compound SMOTE+STAMP test. "
            "SMOTE is recomputed inside every permutation; 0 disables it."
        ),
    )
    parser.add_argument(
        "--allow-expensive-test",
        action="store_true",
        help="Allow >1000 compound permutations across more than 4 tissues.",
    )
    return parser


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
