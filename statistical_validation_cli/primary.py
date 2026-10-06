"""Test for excess switching genes using one-sided age-label permutations on real donors."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

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
    switching_calls_frame,
    switching_state,
    validation_root,
    write_json,
)


def _process_tissue(job: dict[str, Any]) -> dict[str, Any]:
    tissue = job["tissue"]
    out_dir = Path(job["out_dir"])
    stem = safe_filename(tissue)
    checkpoint = out_dir / "checkpoints" / f"{stem}.json"
    calls_path = out_dir / "observed_calls" / f"{stem}.parquet"
    null_path = out_dir / "null_counts" / f"{stem}.npz"
    config = {
        "provenance": job.get("provenance"),
        "version": job["version"],
        "tau": job["tau"],
        "epsilon": job["epsilon"],
        "seed": job["seed"],
        "n_permutations": job["n_permutations"],
    }

    if (
        not job["force"]
        and checkpoint.exists()
        and calls_path.exists()
        and null_path.exists()
    ):
        cached = read_json(checkpoint)
        if config_matches(cached, config):
            record_checkpoint_execution(out_dir, stem, config, "reused")
            return cached["result"]

    data = load_tissue_data(job["version"], tissue, epsilon=job["epsilon"])
    observed_means = age_means_from_codes(data.normalized_values, data.age_codes)
    observed_state = switching_state(observed_means, tau=job["tau"])
    observed_calls = switching_calls_frame(data.gene_ids, observed_state)

    calls_path.parent.mkdir(parents=True, exist_ok=True)
    observed_calls.to_parquet(calls_path, index=False)

    rng = np.random.default_rng(stable_seed(job["seed"], "primary", tissue))
    null_counts = np.empty(job["n_permutations"], dtype=np.int32)
    for permutation_index in range(job["n_permutations"]):
        permuted_codes = rng.permutation(data.age_codes)
        permuted_means = age_means_from_codes(
            data.normalized_values,
            permuted_codes,
        )
        null_counts[permutation_index] = switching_state(
            permuted_means,
            tau=job["tau"],
        ).count

    null_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(null_path, null_counts=null_counts)

    null_std = float(np.std(null_counts, ddof=1)) if null_counts.size > 1 else 0.0
    z_score = (
        float((observed_state.count - np.mean(null_counts)) / null_std)
        if null_std > 0
        else None
    )
    result: dict[str, Any] = {
        "version": job["version"],
        "tissue": tissue,
        "tau": job["tau"],
        "epsilon": job["epsilon"],
        "n_samples": int(data.sample_ids.size),
        "n_genes_input": int(data.raw_gene_ids.size),
        "n_genes_after_epsilon": int(data.gene_ids.size),
        "n_age_brackets": int(np.count_nonzero(data.age_counts)),
        "n_min_age_bin": int(np.min(data.age_counts)),
        "observed_switching": observed_state.count,
        "observed_up": int(np.count_nonzero(observed_state.direction == 1)),
        "observed_down": int(np.count_nonzero(observed_state.direction == -1)),
        "empirical_p_value": empirical_right_tail_p(
            observed_state.count,
            null_counts,
        ),
        "z_score_descriptive": z_score,
        "n_permutations": int(null_counts.size),
        **bin_count_columns(data.age_counts),
        **summary_statistics(null_counts, "null"),
    }
    write_json(checkpoint, {"config": config, "result": result})
    record_checkpoint_execution(out_dir, stem, config, "recomputed")
    return result


def run(args) -> Path:
    if args.permutations < 1:
        raise ValueError("--permutations must be at least 1")
    if args.run_label == "final" and args.permutations < 9_999:
        raise ValueError(
            "The protected final run requires at least 9999 permutations. "
            "Use --run-label smoke or pilot for shorter runs."
        )
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")

    metadata = load_metadata(args.version)
    all_complete_tissues = select_tissues(metadata)
    tissues = select_tissues(
        metadata,
        requested=parse_tissue_arguments(args.tissue),
    )
    full_family = set(tissues) == set(all_complete_tissues)
    enforce_final_scope(
        tissues, all_complete_tissues, args.run_label, "Primary permutation"
    )
    out_dir = validation_root(
        args.version, args.output_root, args.run_label
    ) / "02_primary"
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
            "n_permutations": args.permutations,
            "force": args.force,
            "out_dir": str(out_dir),
        }
        for tissue in tissues
    ]

    results: list[dict[str, Any]] = []
    if args.jobs == 1:
        for index, job in enumerate(jobs, start=1):
            print(f"[{index}/{len(jobs)}] Primary permutation: {job['tissue']}")
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
    allow_partial = bool(getattr(args, "allow_partial_family", False))
    if full_family or allow_partial:
        frame["q_value_bh"] = bh_qvalues(frame["empirical_p_value"])
        frame["significant_bh_0_05"] = frame["q_value_bh"] <= 0.05
        frame["bh_family_size"] = int(frame.shape[0])
        frame["bh_applied"] = True
        frame["bh_family_scope"] = (
            "all_complete_tissues" if full_family else "explicit_partial_family"
        )
    else:
        frame["q_value_bh"] = np.nan
        frame["significant_bh_0_05"] = False
        frame["bh_family_size"] = 0
        frame["bh_applied"] = False
        frame["bh_family_scope"] = "partial_run_no_bh"
        print(
            "Partial tissue run: p-values were computed, but BH was not "
            "applied. Use all complete tissues for the final family or pass "
            "--allow-partial-family for an explicitly defined pilot family."
        )
    result_path = out_dir / "primary_permutation_results.csv"
    frame.to_csv(result_path, index=False)

    ranked = frame.sort_values(
        ["q_value_bh", "empirical_p_value", "tissue"],
        kind="stable",
    )
    ranked.to_csv(out_dir / "primary_permutation_results_ranked.csv", index=False)

    print(f"Primary results written to: {result_path}")
    if frame["bh_applied"].all():
        print(
            "BH-significant tissues: "
            f"{int(frame['significant_bh_0_05'].sum())}/{frame.shape[0]}"
        )
    else:
        print("BH not applied to this partial run.")
    return result_path


def build_parser():
    parser = base_parser(
        "Run the primary tissue-level age-label permutation test on real donors."
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=9_999,
        help="Monte Carlo permutations per tissue (final analysis: 9999).",
    )
    parser.add_argument(
        "--allow-partial-family",
        action="store_true",
        help="Explicitly apply BH to only the selected pilot tissues.",
    )
    return parser


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
