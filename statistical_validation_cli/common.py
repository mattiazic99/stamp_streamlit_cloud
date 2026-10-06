"""Data loading, normalization and inference helpers shared by validation commands."""

from __future__ import annotations

import hashlib
import json
import re
import warnings
import platform
from importlib.metadata import version as package_version, PackageNotFoundError
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from stamp.config import AGE_BRACKETS, DEFAULT_THRESHOLD, EPSILON, GtexVersion, paths_for
from stamp.io import load_metadata, load_tpm_matrix
from stamp.tissues import tissues_with_complete_age_bins


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "output"
AGE_TO_CODE = {age: index for index, age in enumerate(AGE_BRACKETS)}


@dataclass(slots=True)
class TissueData:
    """Sample-level arrays used by the analyses for one tissue."""

    tissue: str
    sample_ids: np.ndarray
    subject_ids: np.ndarray
    age_labels: np.ndarray
    age_codes: np.ndarray
    age_counts: np.ndarray
    raw_gene_ids: np.ndarray
    raw_values: np.ndarray
    kept_indices: np.ndarray
    gene_ids: np.ndarray
    normalized_values: np.ndarray


@dataclass(slots=True)
class SwitchState:
    """Switching status, bracket and direction for each row of an age-bin matrix."""

    mask: np.ndarray
    position: np.ndarray
    direction: np.ndarray

    @property
    def count(self) -> int:
        return int(np.count_nonzero(self.mask))


def validation_root(
    version: str,
    output_root: str | Path | None = None,
    run_label: str = "final",
) -> Path:
    """Return the output directory for a GTEx validation run."""
    if (
        not run_label
        or run_label in {".", ".."}
        or re.fullmatch(r"[A-Za-z0-9_.-]+", run_label) is None
    ):
        raise ValueError("run_label may contain only letters, numbers, _, - and .")
    base = Path(output_root).resolve() if output_root else DEFAULT_OUTPUT_ROOT
    return base / version / "statistical_validation_cli" / "runs" / run_label


def safe_filename(value: str) -> str:
    """Convert a tissue label to a stable name suitable for output files."""
    out = value
    for character in " -()/,":
        out = out.replace(character, "_")
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def stable_seed(base_seed: int, *parts: object) -> int:
    """Derive a uint32 seed independent of Python's randomized hash()."""
    payload = "|".join([str(base_seed), *(str(part) for part in parts)])
    digest = hashlib.sha256(payload.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], byteorder="little", signed=False)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON through a temporary file, then replace the destination."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False),
        encoding="utf-8",
    )
    temporary.replace(path)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def select_tissues(
    metadata: pd.DataFrame,
    requested: Iterable[str] | None = None,
) -> list[str]:
    """Return complete tissues, optionally restricted to explicit names."""
    complete = tissues_with_complete_age_bins(metadata)
    if requested is None:
        return complete

    requested_list = list(dict.fromkeys(requested))
    missing = sorted(set(requested_list) - set(metadata["tissue"].dropna().unique()))
    if missing:
        raise ValueError(f"Unknown tissue names: {missing}")

    incomplete = sorted(set(requested_list) - set(complete))
    if incomplete:
        raise ValueError(
            "Requested tissues without all six age bins: "
            f"{incomplete}. The primary pipeline uses complete tissues only."
        )
    return requested_list


def enforce_final_scope(
    selected: Iterable[str],
    all_complete: Iterable[str],
    run_label: str,
    stage: str,
) -> None:
    """Keep subset runs from replacing the protected final results."""
    selected_set = set(selected)
    complete_set = set(all_complete)
    if run_label == "final" and selected_set != complete_set:
        raise RuntimeError(
            f"{stage}: the protected final run requires all "
            f"{len(complete_set)} complete tissues; received "
            f"{len(selected_set)}. Use --run-label smoke or pilot for subsets."
        )


def parse_tissue_arguments(values: list[str] | None) -> list[str] | None:
    """Parse repeatable/comma-separated ``--tissue`` values."""
    if not values:
        return None
    parsed: list[str] = []
    for value in values:
        parsed.extend(item.strip() for item in value.split(",") if item.strip())
    return list(dict.fromkeys(parsed))


def load_tissue_data(
    version: GtexVersion,
    tissue: str,
    epsilon: float = EPSILON,
) -> TissueData:
    """Load a tissue while keeping normalized values at sample level.
    
    The epsilon filter and min-max formula match stamp.normalize.normalize_tissue.
    Permutations and SMOTE need these values before averaging by age bin.
    """
    metadata = load_metadata(version)
    tissue_meta = metadata.loc[
        (metadata["tissue"] == tissue)
        & metadata["age_bracket"].isin(AGE_BRACKETS)
    ].copy()
    if tissue_meta.empty:
        raise ValueError(f"No valid samples found for tissue '{tissue}'")

    duplicate_subject = tissue_meta.duplicated("subject_id", keep=False)
    if duplicate_subject.any():
        examples = (
            tissue_meta.loc[duplicate_subject, ["subject_id", "sample_id"]]
            .head(10)
            .to_dict("records")
        )
        raise ValueError(
            f"Tissue '{tissue}' contains multiple samples for the same donor. "
            f"Aggregate or block these samples before analysis. Examples: {examples}"
        )

    tissue_meta = tissue_meta.sort_values("sample_id", kind="stable")
    requested_samples = tissue_meta["sample_id"].astype(str).tolist()
    tpm = load_tpm_matrix(version, columns=requested_samples)
    available = [sample for sample in requested_samples if sample in tpm.columns]
    if len(available) != len(requested_samples):
        missing = sorted(set(requested_samples) - set(available))
        raise ValueError(
            f"Tissue '{tissue}' has {len(missing)} metadata samples absent from TPM: "
            f"{missing[:10]}"
        )

    tissue_meta = tissue_meta.set_index("sample_id").loc[available].reset_index()
    raw = tpm.loc[:, available].to_numpy(dtype=np.float32, copy=True)
    raw_gene_ids = tpm.index.astype(str).to_numpy()
    normalized, kept_indices = normalize_sample_values(raw, epsilon=epsilon)

    age_labels = tissue_meta["age_bracket"].astype(str).to_numpy()
    age_codes = np.asarray([AGE_TO_CODE[label] for label in age_labels], dtype=np.int8)
    age_counts = np.bincount(age_codes, minlength=len(AGE_BRACKETS)).astype(np.int32)
    if np.any(age_counts == 0):
        missing_bins = [
            AGE_BRACKETS[index]
            for index, count in enumerate(age_counts)
            if count == 0
        ]
        raise ValueError(
            f"Tissue '{tissue}' is incomplete; missing age bins: {missing_bins}"
        )

    return TissueData(
        tissue=tissue,
        sample_ids=tissue_meta["sample_id"].astype(str).to_numpy(),
        subject_ids=tissue_meta["subject_id"].astype(str).to_numpy(),
        age_labels=age_labels,
        age_codes=age_codes,
        age_counts=age_counts,
        raw_gene_ids=raw_gene_ids,
        raw_values=raw,
        kept_indices=kept_indices,
        gene_ids=raw_gene_ids[kept_indices],
        normalized_values=normalized,
    )


def normalize_sample_values(
    raw_values: np.ndarray,
    epsilon: float = EPSILON,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the canonical per-gene epsilon filter and min-max normalization."""
    values = np.asarray(raw_values, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        gene_min = np.nanmin(values, axis=1)
        gene_max = np.nanmax(values, axis=1)
    gene_range = gene_max - gene_min
    keep = np.isfinite(gene_range) & (gene_range > epsilon)
    kept_indices = np.flatnonzero(keep)
    if kept_indices.size == 0:
        raise ValueError(
            f"No genes remain after applying epsilon={epsilon} to sample matrix"
        )
    normalized = (
        values[kept_indices] - gene_min[kept_indices, None]
    ) / gene_range[kept_indices, None]
    return normalized.astype(np.float32, copy=False), kept_indices


def age_means_from_codes(
    normalized_values: np.ndarray,
    age_codes: np.ndarray,
) -> np.ndarray:
    """Compute gene means for all six bins in chronological order."""
    values = np.asarray(normalized_values, dtype=np.float32)
    codes = np.asarray(age_codes, dtype=np.int8)
    if values.ndim != 2 or codes.ndim != 1 or values.shape[1] != codes.size:
        raise ValueError("normalized_values and age_codes have incompatible shapes")
    if np.any((codes < 0) | (codes >= len(AGE_BRACKETS))):
        raise ValueError("age_codes contains an invalid age-bin code")

    indicator = np.zeros((codes.size, len(AGE_BRACKETS)), dtype=np.float32)
    indicator[np.arange(codes.size), codes] = 1.0
    if np.isfinite(values).all():
        sums = values @ indicator
        denominators = np.bincount(
            codes,
            minlength=len(AGE_BRACKETS),
        ).astype(np.float32)
    else:
        finite = np.isfinite(values)
        finite_values = np.nan_to_num(
            values, nan=0.0, posinf=0.0, neginf=0.0
        )
        sums = finite_values @ indicator
        denominators = finite.astype(np.float32) @ indicator

    means = np.full(
        (values.shape[0], len(AGE_BRACKETS)),
        np.nan,
        dtype=np.float32,
    )
    denominator_matrix = (
        denominators[None, :] if denominators.ndim == 1 else denominators
    )
    np.divide(
        sums,
        denominator_matrix,
        out=means,
        where=denominator_matrix > 0,
    )
    return means


def switching_state(
    age_means: np.ndarray,
    tau: float = DEFAULT_THRESHOLD,
) -> SwitchState:
    """Find genes with one permanent transition, including bracket and direction."""
    if not 0.0 <= tau <= 1.0:
        raise ValueError(f"tau must be in [0, 1], got {tau}")
    values = np.asarray(age_means)
    if values.ndim != 2 or values.shape[1] != len(AGE_BRACKETS):
        raise ValueError(
            "age_means must have shape (genes, 6); "
            f"received {values.shape}"
        )

    valid = np.isfinite(values).all(axis=1)
    binary = values >= tau
    changes = binary[:, 1:] != binary[:, :-1]
    n_changes = np.count_nonzero(changes, axis=1)
    mask = valid & (n_changes == 1)

    position = np.full(values.shape[0], -1, dtype=np.int8)
    direction = np.zeros(values.shape[0], dtype=np.int8)
    selected = np.flatnonzero(mask)
    if selected.size:
        change_index = np.argmax(changes[selected], axis=1)
        position[selected] = (change_index + 1).astype(np.int8)
        before = binary[selected, change_index]
        after = binary[selected, change_index + 1]
        direction[selected] = np.where((~before) & after, 1, -1).astype(np.int8)

    return SwitchState(mask=mask, position=position, direction=direction)


def switching_calls_frame(gene_ids: np.ndarray, state: SwitchState) -> pd.DataFrame:
    """Return one row for every switching gene."""
    selected = np.flatnonzero(state.mask)
    positions = state.position[selected].astype(int)
    return pd.DataFrame(
        {
            "gene_id": np.asarray(gene_ids)[selected],
            "switch_bracket": [AGE_BRACKETS[position] for position in positions],
            "direction": np.where(state.direction[selected] == 1, "up", "down"),
        }
    )


def empirical_right_tail_p(observed: int | float, null_values: np.ndarray) -> float:
    """Return the right-tail Monte Carlo p-value using the +1 correction."""
    null_array = np.asarray(null_values)
    if null_array.size == 0:
        raise ValueError("At least one permutation is required")
    exceedances = int(np.count_nonzero(null_array >= observed))
    return float((exceedances + 1) / (null_array.size + 1))


def bh_qvalues(p_values: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values aligned to the input order."""
    p = np.asarray(list(p_values), dtype=float)
    if p.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    if p.size == 0:
        return p.copy()
    if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("p_values must all be finite and lie in [0, 1]")

    order = np.argsort(p, kind="stable")
    ranked = p[order]
    m = p.size
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


def bin_count_columns(counts: np.ndarray) -> dict[str, int]:
    """Flatten six age-bin counts into stable CSV column names."""
    return {
        f"n_{age.replace('-', '_')}": int(count)
        for age, count in zip(AGE_BRACKETS, counts)
    }


def summary_statistics(values: np.ndarray, prefix: str) -> dict[str, float | int]:
    """Common descriptive summaries for permutation/bootstrap distributions."""
    array = np.asarray(values, dtype=float)
    return {
        f"{prefix}_mean": float(np.mean(array)),
        f"{prefix}_std": float(np.std(array, ddof=1)) if array.size > 1 else 0.0,
        f"{prefix}_median": float(np.median(array)),
        f"{prefix}_q025": float(np.quantile(array, 0.025)),
        f"{prefix}_q975": float(np.quantile(array, 0.975)),
        f"{prefix}_min": int(np.min(array)),
        f"{prefix}_max": int(np.max(array)),
    }


def config_matches(payload: dict[str, Any], expected: dict[str, Any]) -> bool:
    """Check that checkpoint settings and provenance match the current run."""
    return bool(expected.get("provenance")) and payload.get("config") == expected


def checkpoint_provenance(version: str) -> dict:
    """Hash inputs and source in the parent process for reuse by workers."""
    def digest(path):
        hasher = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                hasher.update(chunk)
        return hasher.hexdigest()
    paths = paths_for(version)
    packages = {}
    for name in ("numpy", "pandas", "scipy", "pyarrow", "scikit-learn", "imbalanced-learn"):
        try:
            packages[name] = package_version(name)
        except PackageNotFoundError:
            packages[name] = None
    return {
        "inputs": {key: digest(paths[key]) for key in ("tpm_parquet", "metadata_parquet")},
        "source": {p.relative_to(PROJECT_ROOT).as_posix(): digest(p)
                   for folder in (PROJECT_ROOT / "stamp", PROJECT_ROOT / "statistical_validation_cli")
                   for p in sorted(folder.glob("*.py"))},
        "python": platform.python_version(), "platform": platform.platform(),
        "packages": packages,
    }


def record_checkpoint_execution(out_dir: Path, stem: str, config: dict, status: str) -> None:
    write_json(out_dir / "execution" / f"{stem}.json",
               {"tissue": stem, "status": status, "config": config})


def base_parser(description: str):
    """Build shared command-line options without introducing an import cycle."""
    import argparse

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--version", choices=["v8", "v10"], default="v10")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="Base output directory (default: project output/).",
    )
    parser.add_argument(
        "--run-label",
        default="final",
        help="Isolated output label, for example final, smoke or pilot.",
    )
    parser.add_argument(
        "--tissue",
        action="append",
        default=None,
        help="Restrict to a complete tissue; repeat or comma-separate names.",
    )
    parser.add_argument("--tau", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--epsilon", type=float, default=EPSILON)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Ignore compatible checkpoints and recompute.",
    )
    return parser

