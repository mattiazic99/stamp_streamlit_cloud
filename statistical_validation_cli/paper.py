"""Run and verify the manuscript's statistical experiments.

check-reference validates the saved tables; replay recomputes short random-number
prefixes. A full run recomputes all experiments, including both compound SMOTE tests.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import importlib
import runpy
from importlib.metadata import version

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "reproducibility" / "reference"
LABELS = {"primary": "paper", "bootstrap": "paper", "smote-fixed": "paper-smote-fixed",
          "smote-adaptive": "paper-smote-adaptive"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def environment(strict: bool = True) -> dict:
    expected = read_json(ROOT / "reproducibility/environment.json")
    packages = {name: version(name) for name in expected["packages"]}
    actual = {"python": platform.python_version(), "platform": platform.platform(),
              "packages": packages}
    differences = [f"{name}: expected {value}, installed {packages[name]}"
                   for name, value in expected["packages"].items() if packages[name] != value]
    if sys.version_info[:2] != (3, 11):
        differences.append("Python 3.11 is required by the reference environment")
    if differences and strict:
        raise RuntimeError("Environment differs from the lock file:\n" + "\n".join(differences))
    return actual


def check_reference() -> dict:
    manifest = read_json(REFERENCE / "manifest.json")
    for relative, expected in manifest["files"].items():
        path = REFERENCE / relative
        if not path.is_file() or sha256(path) != expected:
            raise RuntimeError(f"Reference file missing or modified: {relative}")
    from .verify_paper import manuscript_claims, validate_reference_statistics
    validate_reference_statistics(REFERENCE)
    claims = manuscript_claims(REFERENCE)
    if claims != manifest["claims"]:
        raise RuntimeError("Reference tables do not support the manuscript claims")
    print("Reference tables support all recorded manuscript claims (no experiments rerun).")
    return claims


def input_directory(value: Path | None) -> Path:
    return (value or ROOT / "data/parquet/v10").resolve()


def check_inputs(directory: Path, demo: bool = False) -> dict:
    expected = read_json(ROOT / "reproducibility/inputs.json")
    expected = expected["demo" if demo else "full"]
    observed = {}
    for name, entry in expected["files"].items():
        path = directory / name
        print(f"Checking input {name} ...", flush=True)
        if not path.is_file() or path.stat().st_size != entry["bytes"]:
            raise RuntimeError(f"Missing or different input: {path}")
        observed[name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
        if observed[name]["sha256"] != entry["sha256"]:
            raise RuntimeError(f"Input SHA-256 mismatch: {path}")
    return observed


def configure_inputs(directory: Path) -> None:
    # Pass the input location to subprocesses and Windows workers.
    os.environ["STAMP_VALIDATION_INPUT_DIR"] = str(directory)
    os.environ["MPLBACKEND"] = "Agg"
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
        os.environ[name] = "1"


def replay(directory: Path, tissues: list[str], demo: bool) -> None:
    """Compare freshly recomputed prefixes with the historical full-run results."""
    configure_inputs(directory)
    import numpy as np
    import pandas as pd
    from threadpoolctl import threadpool_limits
    from .common import (load_tissue_data, switching_state, age_means_from_codes,
                         stable_seed, safe_filename)
    from .bootstrap import _bootstrap_draws, _fixed_normalization_state
    from .smote import _smote_once, _strategy_and_neighbors
    from .verify_paper import compare_frames

    configure_inputs(directory)
    primary = pd.read_csv(REFERENCE / "primary.csv").set_index("tissue")
    with threadpool_limits(limits=1):
        for tissue in tissues:
            data = load_tissue_data("v10", tissue, epsilon=0.01)
            observed = switching_state(age_means_from_codes(data.normalized_values, data.age_codes))
            expected = primary.loc[tissue]
            if (observed.count != int(expected.observed_switching)
                or int(np.count_nonzero(observed.direction == 1)) != int(expected.observed_up)
                or int(np.count_nonzero(observed.direction == -1)) != int(expected.observed_down)):
                raise AssertionError("Observed switching/direction counts differ: " + tissue)
            fixture = np.load(REFERENCE / "replay" / (safe_filename(tissue) + ".npz"))
            rng = np.random.default_rng(stable_seed(42, "primary", tissue))
            actual = [switching_state(age_means_from_codes(data.normalized_values,
                      rng.permutation(data.age_codes))).count for _ in fixture["primary"]]
            np.testing.assert_array_equal(actual, fixture["primary"], err_msg=tissue + " primary")
            rng = np.random.default_rng(stable_seed(42, "bootstrap", tissue))
            actual = [_fixed_normalization_state(data.normalized_values,
                      _bootstrap_draws(data.age_codes, rng), 0.5).count for _ in fixture["bootstrap"]]
            np.testing.assert_array_equal(actual, fixture["bootstrap"], err_msg=tissue + " bootstrap")
            for policy, key in (("fixed-floor", "fixed"), ("adaptive-second-smallest", "adaptive")):
                strategy, neighbors, _ = _strategy_and_neighbors(data.age_counts, 20, policy)
                rows = []
                for index in range(2):
                    seed = stable_seed(42, "smote_seed_run", tissue, index)
                    augmented, codes = _smote_once(data.normalized_values, data.age_codes,
                                                  strategy, neighbors, seed)
                    state = switching_state(age_means_from_codes(augmented, codes))
                    union = np.count_nonzero(state.mask | observed.mask)
                    rows.append({"seed_index": index, "random_state": seed,
                                 "smote_switching": state.count,
                                 "jaccard_real_vs_smote": np.count_nonzero(state.mask & observed.mask) / union if union else 1.0})
                seed_reference = pd.read_csv(REFERENCE / "seed_runs" / key / (safe_filename(tissue) + ".csv"))
                compare_frames(pd.DataFrame(rows), seed_reference.head(2)[list(rows[0])], ["seed_index"])
                if strategy:
                    seed = stable_seed(42, "smote_compound_test", tissue)
                    augmented, codes = _smote_once(data.normalized_values, data.age_codes, strategy, neighbors, seed)
                    observed_compound = switching_state(age_means_from_codes(augmented, codes)).count
                    if observed_compound != int(fixture[key + "_observed"]):
                        raise AssertionError("Compound observed count differs: " + tissue + " " + policy)
                    rng = np.random.default_rng(stable_seed(42, "smote_label_permutations", tissue))
                    actual = []
                    for _ in fixture[key]:
                        augmented, codes = _smote_once(data.normalized_values, rng.permutation(data.age_codes), strategy, neighbors, seed)
                        actual.append(switching_state(age_means_from_codes(augmented, codes)).count)
                    np.testing.assert_array_equal(actual, fixture[key], err_msg=tissue + " " + policy)
            print(f"{tissue}: fresh primary, bootstrap and both SMOTE prefixes match.", flush=True)
    print("Replay passed. Short prefixes do not replace the complete 9,999-permutation run.")


def full_command_plan(output: Path, jobs: int, force: bool) -> list[tuple[str, list[str]]]:
    common = ["--version", "v10", "--tau", "0.5", "--epsilon", "0.01", "--seed", "42",
              "--jobs", str(jobs), "--output-root", str(output)]
    if force:
        common.append("--force")
    steps = [
        ("audit", ["--run-label", "paper"]),
        ("regression", ["--run-label", "paper"]),
        ("primary", ["--run-label", "paper", "--permutations", "9999"]),
        ("bootstrap", ["--run-label", "paper", "--replicates", "1000", "--normalization", "fixed"]),
    ]
    for policy, label in (("fixed-floor", "paper-smote-fixed"), ("adaptive-second-smallest", "paper-smote-adaptive")):
        steps.append(("smote", ["--run-label", label, "--target-policy", policy, "--seeds", "20",
                               "--target-n", "20", "--permutations", "9999", "--allow-expensive-test"]))
    steps.append(("report", ["--run-label", "paper"]))
    return [(module, [*common, *extra]) for module, extra in steps]


def validate_command_plan(plan) -> None:
    for module, arguments in plan:
        importlib.import_module("statistical_validation_cli." + module).build_parser().parse_args(arguments)


def threshold_arguments(output: Path) -> list[str]:
    return ["--version", "v10", "--sample-level", "--primary-summary",
            str(output / "v10/statistical_validation_cli/runs/paper/02_primary/primary_permutation_results.csv"),
            "--output-dir", str(output / "paper-threshold")]


def validate_threshold_command(output: Path) -> None:
    script = runpy.run_path(str(ROOT / "scripts/04_threshold_sensitivity.py"))
    script["build_parser"]().parse_args(threshold_arguments(output))


def prepare_output(output: Path, provenance: dict, force: bool) -> Path:
    """Require matching provenance before reusing a non-empty output directory."""
    provenance_path = output / "paper_provenance.json"
    if output.exists() and any(output.iterdir()) and not provenance_path.is_file() and not force:
        raise RuntimeError("Non-empty output has no provenance. Choose a new output directory or --force.")
    if provenance_path.exists() and read_json(provenance_path) != provenance and not force:
        raise RuntimeError("Inputs, environment or source changed. Choose a new output directory or --force.")
    output.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    (output / "paper_verified.json").unlink(missing_ok=True)
    return provenance_path


def execution_summary(output: Path) -> dict:
    summaries = {}
    for label, step in (("paper", "02_primary"), ("paper", "04_bootstrap"),
                        ("paper-smote-fixed", "03_smote"), ("paper-smote-adaptive", "03_smote")):
        folder = output / "v10/statistical_validation_cli/runs" / label / step / "execution"
        events = [read_json(path) for path in sorted(folder.glob("*.json"))]
        if len(events) != 50 or any(e["status"] not in ("recomputed", "reused") for e in events):
            raise RuntimeError("Missing execution records: " + str(folder))
        summaries[label + "/" + step] = {status: [e["tissue"] for e in events if e["status"] == status]
                                          for status in ("recomputed", "reused")}
    return summaries


def run_full(args) -> None:
    directory = input_directory(args.input_dir)
    configure_inputs(directory)
    output = args.output_root.resolve()
    plan = full_command_plan(output, args.jobs, args.force)
    # Validate child command arguments before launching the experiments.
    validate_command_plan(plan)
    validate_threshold_command(output)
    inputs = check_inputs(directory)
    env = environment()
    check_reference()
    configure_inputs(directory)
    from .common import select_tissues
    from stamp.io import load_metadata
    if len(select_tissues(load_metadata("v10"))) != 50:
        raise RuntimeError("The manuscript reproduction requires all 50 complete v10 tissues")
    provenance = {"inputs": inputs, "environment": env, "seed": 42,
                  "reference_manifest_sha256": sha256(REFERENCE / "manifest.json"),
                  "source_sha256": {str(p.relative_to(ROOT)): sha256(p)
                                    for folder in (ROOT / "stamp", ROOT / "statistical_validation_cli")
                                    for p in sorted(folder.glob("*.py"))}}
    threshold_script = ROOT / "scripts/04_threshold_sensitivity.py"
    provenance["source_sha256"][str(threshold_script.relative_to(ROOT))] = sha256(threshold_script)
    # Include the atlas in provenance because regression reads it.
    provenance["atlas_sha256"] = {p.relative_to(ROOT).as_posix(): sha256(p)
                                  for sub in ("normalized", "sets")
                                  for p in sorted((ROOT / "output/v10_complete" / sub).glob("*")) if p.is_file()}
    provenance_path = prepare_output(output, provenance, args.force)

    def command(module, arguments):
        print("Running " + module + " " + " ".join(arguments), flush=True)
        subprocess.run([sys.executable, "-m", "statistical_validation_cli." + module,
                        *arguments], cwd=ROOT, check=True)

    for module, arguments in plan[:-1]:
        command(module, arguments)
    subprocess.run([sys.executable, "scripts/04_threshold_sensitivity.py",
                    *threshold_arguments(output)], cwd=ROOT, check=True)
    from .verify_paper import verify_outputs
    verify_outputs(output, REFERENCE)
    command(*plan[-1])
    (output / "paper_verified.json").write_text(json.dumps({"status": "all manuscript numerical results matched",
        "provenance_sha256": sha256(provenance_path),
        "execution": execution_summary(output)}, indent=2) + "\n", encoding="utf-8")
    print("Full reproduction verified against the manuscript references.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["check-reference", "check-inputs", "replay", "run", "verify"])
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--demo", action="store_true", help="Use the bundled Ovary demo, only for replay/check-inputs.")
    parser.add_argument("--tissue", action="append", help="Replay selected tissues; default Ovary; use --all-tissues for all 50.")
    parser.add_argument("--all-tissues", action="store_true")
    parser.add_argument("--output-root", type=Path, default=ROOT / "output/reproduction")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.demo and args.action not in ("replay", "check-inputs"):
        parser.error("The one-tissue demo cannot reproduce the complete tissue-level inference.")
    if args.action == "check-reference":
        check_reference()
    elif args.action == "verify":
        check_reference()
        from .verify_paper import verify_outputs
        verify_outputs(args.output_root.resolve(), REFERENCE)
    elif args.action == "run":
        run_full(args)
    else:
        directory = (args.input_dir or (ROOT / "data/demo/v10" if args.demo else ROOT / "data/parquet/v10")).resolve()
        check_inputs(directory, args.demo)
        if args.action == "replay":
            check_reference()
            environment()
            tissues = args.tissue or ["Ovary"]
            if args.all_tissues:
                tissues = read_json(REFERENCE / "manifest.json")["tissues"]
            if args.demo and tissues != ["Ovary"]:
                parser.error("The demo contains only Ovary.")
            replay(directory, tissues, args.demo)


if __name__ == "__main__":
    main()
