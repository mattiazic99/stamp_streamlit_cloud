"""Central configuration: paths, constants, GTEx version handling.

The pipeline is parametric on the GTEx version. Use `paths_for(version)`
to get version-specific paths.
"""
import os
from pathlib import Path
from typing import Literal

GtexVersion = Literal["v8", "v10"]
SUPPORTED_VERSIONS: tuple[GtexVersion, ...] = ("v8", "v10")

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
PARQUET_DIR = DATA_DIR / "parquet"
EXTERNAL_DIR = DATA_DIR / "external"
OUTPUT_DIR = ROOT / "output"
COMPARE_DIR = OUTPUT_DIR / "compare"

# Six age brackets used in the paper.
AGE_BRACKETS = ("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
# Switches are assigned from the second bracket onward.
SWITCHING_BRACKETS = AGE_BRACKETS[1:]

DEFAULT_THRESHOLD = 0.5
EPSILON = 0.01

# Complete age-bin coverage
# Complete tissues have samples in all six age brackets. The all-tissues
# mode keeps incomplete tissues for sensitivity analyses.
# These reference sets document the GTEx dumps and support tests; filtering
# uses stamp.tissues.tissues_with_complete_age_bins(), computed from metadata.
INCOMPLETE_TISSUES_BY_VERSION: dict[str, set[str]] = {
    "v8": {
        "Bladder",
        "Cervix - Ectocervix",
        "Cervix - Endocervix",
        "Fallopian Tube",
        "Kidney - Medulla",
    },
    "v10": {
        "Cervix - Ectocervix",
        "Cervix - Endocervix",
        "Fallopian Tube",
        "Kidney - Medulla",
    },
}

# Union of incomplete tissues across v8 and v10: excluded from BOTH versions
# when a like-for-like v8-vs-v10 comparison on identical tissues is required.
INCOMPLETE_TISSUES_UNION_V8_V10: set[str] = (
    INCOMPLETE_TISSUES_BY_VERSION["v8"] | INCOMPLETE_TISSUES_BY_VERSION["v10"]
)

CASSANDRA_HOSTS = ("127.0.0.1",)
CASSANDRA_PORT = 9042
KEYSPACE = "gtex_keyspace"
TPM_TABLE = "gene_tpm_full"
METADATA_TABLE = "sample_metadata_completa"


def paths_for(version: GtexVersion, complete: bool = False) -> dict[str, Path]:
    """Return all version-specific paths.

    Parameters
    ----------
    version : "v8" or "v10"
    complete : bool, default False
        If True, derived outputs (normalized, sets, jaccard,
        threshold_sensitivity) are routed to ``output/{version}_complete``
        instead of ``output/{version}``. This keeps the complete-age-bins
        analysis fully separate from the all-tissues analysis so neither
        overwrites the other.

        The raw inputs (``tpm_parquet``, ``metadata_parquet``) are shared
        between modes and always live under ``data/parquet/{version}``.
    """
    if version not in SUPPORTED_VERSIONS:
        raise ValueError(f"Unknown GTEx version: {version}")
    out_name = f"{version}_complete" if complete else version
    out = OUTPUT_DIR / out_name
    # Optional exact paper inputs; inherited by validation worker processes.
    input_dir = (Path(os.environ["STAMP_VALIDATION_INPUT_DIR"])
                 if version == "v10" and os.environ.get("STAMP_VALIDATION_INPUT_DIR")
                 else PARQUET_DIR / version)
    return {
        "tpm_parquet": input_dir / "tpm_matrix.parquet",
        "metadata_parquet": input_dir / "metadata.parquet",
        "normalized": out / "normalized",
        "sets": out / "sets",
        "jaccard": out / "jaccard",
        "threshold_sensitivity": out / "threshold_sensitivity",
    }
