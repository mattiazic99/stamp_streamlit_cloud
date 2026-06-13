"""Central configuration: paths, constants, GTEx version handling.

The pipeline is parametric on the GTEx version. Use `paths_for(version)`
to get version-specific paths.
"""
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

# 6 fasce di età come da paper
AGE_BRACKETS = ("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
# Lo switch è valutato dalla seconda fascia in poi
SWITCHING_BRACKETS = AGE_BRACKETS[1:]

DEFAULT_THRESHOLD = 0.5
EPSILON = 0.01

# ---------------------------------------------------------------------------
# Complete age-bin coverage
# ---------------------------------------------------------------------------
# A tissue has "complete age-bin coverage" if it has at least one valid sample
# in ALL six GTEx age brackets. Tissues missing one or more brackets are
# excluded from the main "complete" analysis (see stamp.tissues), but remain
# available as a sensitivity analysis through the default (all-tissues) mode.
#
# The sets below are a documented reference of the incomplete tissues observed
# in the current GTEx dumps (computed dynamically by
# stamp.tissues.tissues_with_complete_age_bins from the metadata; these
# constants are NOT used for filtering, only for documentation / tests).
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
        If True, derived outputs (normalized, sets, jaccard, migration,
        threshold_sweep, plots, ppi) are routed to ``output/{version}_complete``
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
    return {
        "tpm_parquet": PARQUET_DIR / version / "tpm_matrix.parquet",
        "metadata_parquet": PARQUET_DIR / version / "metadata.parquet",
        "normalized": out / "normalized",
        "sets": out / "sets",
        "jaccard": out / "jaccard",
        "ppi": out / "ppi",
        "migration": out / "migration",
        "threshold_sweep": out / "threshold_sweep",
        "plots": out / "plots",
    }
