"""Cached loaders for the switching atlas, normalized matrices and Jaccard tables.

paths_for() selects the GTEx release and complete-age-bin directory. Main pages
use v10; Panel Explorer requests v8 separately for the release comparison.
The repository root is added to sys.path so the loaders work without installing
the stamp package on Streamlit Cloud.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Allow direct Streamlit startup without installing stamp.
_STAMP_ROOT = Path(__file__).resolve().parent.parent
if str(_STAMP_ROOT) not in sys.path:
    sys.path.insert(0, str(_STAMP_ROOT))

import pandas as pd
import streamlit as st

from stamp.config import paths_for, SWITCHING_BRACKETS


# Released configuration

#: Main atlas release. Panel Explorer requests v8 separately for comparison.
GTEX_VERSION = "v10"

#: Require all six age brackets, without imputing missing brackets.
#: Load the atlas from output/{version}_complete/.
COMPLETE_AGE_BINS = True


# Tissue name helpers

# Hardcoded safe → display mapping for all 54 GTEx tissues.
# The key is the filesystem-safe name (underscores), the value is the
# human-readable GTEx display name used in papers and metadata.
_SAFE_TO_DISPLAY: dict[str, str] = {
    "Adipose_Subcutaneous": "Adipose - Subcutaneous",
    "Adipose_Visceral_Omentum": "Adipose - Visceral (Omentum)",
    "Adrenal_Gland": "Adrenal Gland",
    "Artery_Aorta": "Artery - Aorta",
    "Artery_Coronary": "Artery - Coronary",
    "Artery_Tibial": "Artery - Tibial",
    "Bladder": "Bladder",
    "Brain_Amygdala": "Brain - Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24": "Brain - Anterior cingulate cortex (BA24)",
    "Brain_Caudate_basal_ganglia": "Brain - Caudate (basal ganglia)",
    "Brain_Cerebellar_Hemisphere": "Brain - Cerebellar Hemisphere",
    "Brain_Cerebellum": "Brain - Cerebellum",
    "Brain_Cortex": "Brain - Cortex",
    "Brain_Frontal_Cortex_BA9": "Brain - Frontal Cortex (BA9)",
    "Brain_Hippocampus": "Brain - Hippocampus",
    "Brain_Hypothalamus": "Brain - Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia": "Brain - Nucleus accumbens (basal ganglia)",
    "Brain_Putamen_basal_ganglia": "Brain - Putamen (basal ganglia)",
    "Brain_Spinal_cord_cervical_c_1": "Brain - Spinal cord (cervical c-1)",
    "Brain_Substantia_nigra": "Brain - Substantia nigra",
    "Breast_Mammary_Tissue": "Breast - Mammary Tissue",
    "Cells_Cultured_fibroblasts": "Cells - Cultured fibroblasts",
    "Cells_EBV_transformed_lymphocytes": "Cells - EBV-transformed lymphocytes",
    "Cervix_Ectocervix": "Cervix - Ectocervix",
    "Cervix_Endocervix": "Cervix - Endocervix",
    "Colon_Sigmoid": "Colon - Sigmoid",
    "Colon_Transverse": "Colon - Transverse",
    "Esophagus_Gastroesophageal_Junction": "Esophagus - Gastroesophageal Junction",
    "Esophagus_Mucosa": "Esophagus - Mucosa",
    "Esophagus_Muscularis": "Esophagus - Muscularis",
    "Fallopian_Tube": "Fallopian Tube",
    "Heart_Atrial_Appendage": "Heart - Atrial Appendage",
    "Heart_Left_Ventricle": "Heart - Left Ventricle",
    "Kidney_Cortex": "Kidney - Cortex",
    "Kidney_Medulla": "Kidney - Medulla",
    "Liver": "Liver",
    "Lung": "Lung",
    "Minor_Salivary_Gland": "Minor Salivary Gland",
    "Muscle_Skeletal": "Muscle - Skeletal",
    "Nerve_Tibial": "Nerve - Tibial",
    "Ovary": "Ovary",
    "Pancreas": "Pancreas",
    "Pituitary": "Pituitary",
    "Prostate": "Prostate",
    "Skin_Not_Sun_Exposed_Suprapubic": "Skin - Not Sun Exposed (Suprapubic)",
    "Skin_Sun_Exposed_Lower_leg": "Skin - Sun Exposed (Lower leg)",
    "Small_Intestine_Terminal_Ileum": "Small Intestine - Terminal Ileum",
    "Spleen": "Spleen",
    "Stomach": "Stomach",
    "Testis": "Testis",
    "Thyroid": "Thyroid",
    "Uterus": "Uterus",
    "Vagina": "Vagina",
    "Whole_Blood": "Whole Blood",
}

# Reverse mapping: display → safe
_DISPLAY_TO_SAFE: dict[str, str] = {v: k for k, v in _SAFE_TO_DISPLAY.items()}


def safe_to_display(tissue_safe: str) -> str:
    """Convert a filesystem-safe tissue name to its GTEx display name.

    Falls back to replacing underscores with spaces if the tissue is
    not found in the hardcoded mapping.
    """
    return _SAFE_TO_DISPLAY.get(tissue_safe, tissue_safe.replace("_", " "))


def display_to_safe(tissue_display: str) -> str:
    """Convert a GTEx display name back to its filesystem-safe form."""
    return _DISPLAY_TO_SAFE.get(tissue_display, tissue_display.replace(" ", "_"))


# Tissue listing

@st.cache_data(ttl=600)
def get_available_tissues(version: str, complete: bool = False) -> list[str]:
    """Return **display** names of tissues with pre-computed sets.

    The list is sorted alphabetically by display name. When ``complete`` is
    True, reads from ``output/{version}_complete/`` (only tissues with all six
    age brackets).
    """
    sets_dir = paths_for(version, complete)["sets"]
    safe_names = sorted(
        p.stem.replace("_sets", "")
        for p in sets_dir.glob("*_sets.txt")
    )
    return [safe_to_display(s) for s in safe_names]


@st.cache_data(ttl=600)
def get_available_tissues_safe(version: str, complete: bool = False) -> list[str]:
    """Return **filesystem-safe** names (underscores) of available tissues."""
    sets_dir = paths_for(version, complete)["sets"]
    return sorted(
        p.stem.replace("_sets", "")
        for p in sets_dir.glob("*_sets.txt")
    )


# Sets loading

@st.cache_data(ttl=600)
def get_sets_for_tissue(
    version: str, tissue_display: str, complete: bool = False
) -> dict[str, list[str]]:
    """Load the switching gene sets for a tissue (by display name).

    Returns a dict mapping bracket label → list of gene IDs.
    Keys are the SWITCHING_BRACKETS: '30-39', '40-49', '50-59', '60-69', '70-79'.
    """
    from stamp.io import load_sets_txt
    tissue_safe = display_to_safe(tissue_display)
    return load_sets_txt(version, tissue_safe, complete=complete)


# Jaccard matrices

@st.cache_data(ttl=600)
def get_jaccard_matrix(
    version: str, metric: str = "life", complete: bool = False
) -> pd.DataFrame:
    """Load a pre-computed Jaccard similarity matrix.

    Parameters
    ----------
    metric : str
        ``"life"`` for whole-lifespan, ``"age"`` for age-averaged.
    complete : bool
        If True, read the complete-age-bins matrix (49×49 for v8, 50×50 for
        v10) from ``output/{version}_complete/``.

    Returns
    -------
    DataFrame with tissue names as both index and columns.
    """
    path = paths_for(version, complete)["jaccard"] / f"jaccard_{metric}.csv"
    return pd.read_csv(path, index_col=0)


# Normalized matrix (for interactive tau slider)

@st.cache_data(ttl=600)
def get_normalized_tissue(
    version: str, tissue_display: str, complete: bool = False
) -> pd.DataFrame:
    """Load the normalized expression matrix for a single tissue.

    Returns DataFrame with rows = gene_id, columns = age brackets.
    Typically < 2 MB per tissue.
    """
    from stamp.io import load_normalized_tissue
    tissue_safe = display_to_safe(tissue_display)
    return load_normalized_tissue(version, tissue_safe, complete=complete)


# Bridge: convert backend sets dict → GUI format

# The existing GUI pages expect the data in the format produced by
# ``parse_stamp_file()``:
# (used_age_groups, gene_counts_per_age, flattened_dataframe)
# where age groups use en-dashes (30–39) and the DataFrame has columns
# ["Age", "Gene"].

_GUI_AGE_GROUPS = ["30–39", "40–49", "50–59", "60–69", "70–79"]

# Map backend bracket labels (hyphens) → GUI labels (en-dashes)
_BRACKET_TO_GUI = dict(zip(SWITCHING_BRACKETS, _GUI_AGE_GROUPS))


def sets_to_gui_format(
    sets_dict: dict[str, list[str]],
) -> tuple[list[str], list[int], pd.DataFrame]:
    """Convert a ``{bracket: [gene_id, ...]}`` dict to the GUI triple.

    Returns
    -------
    used_age_groups : list[str]
        The 5 age-group labels with en-dashes (e.g. ``"30–39"``).
    gene_counts : list[int]
        Number of genes per age group, in order.
    df : pd.DataFrame
        Flattened DataFrame with columns ``["Age", "Gene"]``.
    """
    rows: list[dict[str, str]] = []
    counts: list[int] = []

    for bracket in SWITCHING_BRACKETS:
        gui_label = _BRACKET_TO_GUI[bracket]
        genes = sets_dict.get(bracket, [])
        counts.append(len(genes))
        for g in genes:
            rows.append({"Age": gui_label, "Gene": g})

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["Age", "Gene"])
    return list(_GUI_AGE_GROUPS), counts, df


def sets_to_gene_sets(
    sets_dict: dict[str, list[str]],
) -> list[set[str]]:
    """Convert a ``{bracket: [gene_id, ...]}`` dict to a list of sets.

    Returns a list of 5 sets (one per age bracket), in the same order
    used by ``parse_multiple_stamp_files`` in ``utils/parsing.py``.
    This is the format expected by multi-tissue analysis functions.
    """
    return [set(sets_dict.get(bracket, [])) for bracket in SWITCHING_BRACKETS]
