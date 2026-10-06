"""Estimate switching-call stability with a donor bootstrap stratified by age bin."""

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
    bin_count_columns,
    config_matches,
    checkpoint_provenance,
    record_checkpoint_execution,
    enforce_final_scope,
    load_tissue_data,
    normalize_sample_values,
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


def _bootstrap_draws(
    age_codes: np.ndarray,
    rng: np.random.Generator,
) -> list[np.ndarray]:
    """Resample donors within each age bin, retaining the original bin sizes."""
    draws: list[np.ndarray] = []
    for code in range(6):
        candidates = np.flatnonzero(age_codes == code)
        draws.append(rng.choice(candidates, size=candidates.size, replace=True))
    return draws


def _fixed_normalization_state(
    normalized_values: np.ndarray,
    draws: list[np.ndarray],
    tau: float,
):
    selected = np.concatenate(draws)
    bootstrap_codes = np.concatenate(
        [
            np.full(draw.size, code, dtype=np.int8)
            for code, draw in enumerate(draws)
        ]
    )
    means = age_means_from_codes(
        normalized_values[:, selected],
        bootstrap_codes,
    )
    return switching_state(means, tau=tau)


def _refit_normalization_state(
    raw_values: np.ndarray,
    draws: list[np.ndarray],
    epsilon: float,
    tau: float,
):
    selected = np.concatenate(draws)
    bootstrap_raw = raw_values[:, selected]
    normalized, kept_indices = normalize_sample_values(
        bootstrap_raw,
        epsilon=epsilon,
    )
    bootstrap_codes = np.concatenate(
        [
            np.full(draw.size, code, dtype=np.int8)
            for code, draw in enumerate(draws)
        ]
    )
    state = switching_state(
        age_means_from_codes(normalized, bootstrap_codes),
        tau=tau,
    )
    return state, kept_indices


def _process_tissue(job: dict[str, Any]) -> dict[str, Any]:
    tissue = job["tissue"]
    out_dir = Path(job["out_dir"])
    stem = safe_filename(tissue)
    checkpoint = out_dir / "checkpoints" / f"{stem}.json"
    gene_path = out_dir / "gene_stability" / f"{stem}.parquet"
    count_path = out_dir / "bootstrap_counts" / f"{stem}.npz"
    config = {
        "provenance": job.get("provenance"),
        "version": job["version"],
        "tau": job["tau"],
        "epsilon": job["epsilon"],
        "seed": job["seed"],
        "n_bootstrap": job["n_bootstrap"],
        "normalization_mode": job["normalization_mode"],
    }
    if (
        not job["force"]
        and checkpoint.exists()
        and gene_path.exists()
        and count_path.exists()
    ):
        cached = read_json(checkpoint)
        if config_matches(cached, config):
            record_checkpoint_execution(out_dir, stem, config, "reused")
            return cached["result"]

    data = load_tissue_data(job["version"], tissue, epsilon=job["epsilon"])
    observed_state_kept = switching_state(
        age_means_from_codes(data.normalized_values, data.age_codes),
        tau=job["tau"],
    )

    if job["normalization_mode"] == "fixed":
        base_gene_ids = data.gene_ids
        observed_mask = observed_state_kept.mask
        observed_position = observed_state_kept.position
        observed_direction = observed_state_kept.direction
    else:
        base_gene_ids = data.raw_gene_ids
        observed_mask = np.zeros(data.raw_gene_ids.size, dtype=bool)
        observed_position = np.full(data.raw_gene_ids.size, -1, dtype=np.int8)
        observed_direction = np.zeros(data.raw_gene_ids.size, dtype=np.int8)
        observed_mask[data.kept_indices] = observed_state_kept.mask
        observed_position[data.kept_indices] = observed_state_kept.position
        observed_direction[data.kept_indices] = observed_state_kept.direction

    switch_frequency = np.zeros(base_gene_ids.size, dtype=np.int32)
    exact_call_frequency = np.zeros(base_gene_ids.size, dtype=np.int32)
    bootstrap_counts = np.empty(job["n_bootstrap"], dtype=np.int32)
    rng = np.random.default_rng(stable_seed(job["seed"], "bootstrap", tissue))

    for bootstrap_index in range(job["n_bootstrap"]):
        draws = _bootstrap_draws(data.age_codes, rng)
        if job["normalization_mode"] == "fixed":
            state = _fixed_normalization_state(
                data.normalized_values,
                draws,
                tau=job["tau"],
            )
            switch_frequency += state.mask.astype(np.int32)
            exact_call_frequency += (
                state.mask
                & observed_mask
                & (state.position == observed_position)
                & (state.direction == observed_direction)
            ).astype(np.int32)
        else:
            state, kept_indices = _refit_normalization_state(
                data.raw_values,
                draws,
                epsilon=job["epsilon"],
                tau=job["tau"],
            )
            switch_frequency[kept_indices] += state.mask.astype(np.int32)
            exact = (
                state.mask
                & observed_mask[kept_indices]
                & (state.position == observed_position[kept_indices])
                & (state.direction == observed_direction[kept_indices])
            )
            exact_call_frequency[kept_indices] += exact.astype(np.int32)
        bootstrap_counts[bootstrap_index] = state.count

    count_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(count_path, bootstrap_counts=bootstrap_counts)

    gene_frame = pd.DataFrame(
        {
            "gene_id": base_gene_ids,
            "observed_switching": observed_mask,
            "bootstrap_switch_frequency": (
                switch_frequency / job["n_bootstrap"]
            ),
            "same_bracket_and_direction_frequency": (
                exact_call_frequency / job["n_bootstrap"]
            ),
        }
    )
    gene_frame = gene_frame.loc[
        gene_frame["observed_switching"]
        | (gene_frame["bootstrap_switch_frequency"] > 0)
    ]
    gene_path.parent.mkdir(parents=True, exist_ok=True)
    gene_frame.to_parquet(gene_path, index=False)

    observed_gene_stability = gene_frame.loc[
        gene_frame["observed_switching"],
        "same_bracket_and_direction_frequency",
    ]
    result: dict[str, Any] = {
        "version": job["version"],
        "tissue": tissue,
        "tau": job["tau"],
        "epsilon": job["epsilon"],
        "normalization_mode": job["normalization_mode"],
        "n_samples": int(data.sample_ids.size),
        "n_min_age_bin": int(np.min(data.age_counts)),
        "observed_switching": observed_state_kept.count,
        "n_bootstrap": int(job["n_bootstrap"]),
        "observed_gene_exact_stability_mean": (
            float(observed_gene_stability.mean())
            if not observed_gene_stability.empty
            else None
        ),
        "observed_genes_stable_0_80": int(
            np.count_nonzero(observed_gene_stability >= 0.80)
        ),
        "observed_genes_stable_0_90": int(
            np.count_nonzero(observed_gene_stability >= 0.90)
        ),
        **bin_count_columns(data.age_counts),
        **summary_statistics(bootstrap_counts, "bootstrap"),
    }
    write_json(checkpoint, {"config": config, "result": result})
    record_checkpoint_execution(out_dir, stem, config, "recomputed")
    return result


def run(args) -> Path:
    if args.replicates < 1:
        raise ValueError("--replicates must be at least 1")
    if args.run_label == "final" and args.replicates < 1_000:
        raise ValueError(
            "The protected final bootstrap requires at least 1000 replicates. "
            "Use a non-final run label for shorter runs."
        )
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")

    metadata = load_metadata(args.version)
    all_complete = select_tissues(metadata)
    tissues = select_tissues(
        metadata,
        requested=parse_tissue_arguments(args.tissue),
    )
    enforce_final_scope(tissues, all_complete, args.run_label, "Bootstrap")
    out_dir = validation_root(
        args.version, args.output_root, args.run_label
    ) / "04_bootstrap"
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
            "n_bootstrap": args.replicates,
            "normalization_mode": args.normalization,
            "force": args.force,
            "out_dir": str(out_dir),
        }
        for tissue in tissues
    ]

    results: list[dict[str, Any]] = []
    if args.jobs == 1:
        for index, job in enumerate(jobs, start=1):
            print(f"[{index}/{len(jobs)}] Bootstrap: {job['tissue']}")
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
    result_path = out_dir / "bootstrap_stability_results.csv"
    frame.to_csv(result_path, index=False)
    print(f"Bootstrap stability results written to: {result_path}")
    print(
        "Bootstrap outputs are stability estimates, not null-hypothesis "
        "p-values; BH is therefore not applied here."
    )
    return result_path


def build_parser():
    parser = base_parser(
        "Run donor-level stratified bootstrap stability analysis."
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=1_000,
        help="Bootstrap resamples per tissue.",
    )
    parser.add_argument(
        "--normalization",
        choices=["fixed", "refit"],
        default="fixed",
        help=(
            "`fixed` conditions on the original STAMP scale and is faster; "
            "`refit` re-estimates epsilon filtering/min-max normalization in "
            "every bootstrap and measures full-pipeline instability."
        ),
    )
    return parser


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
