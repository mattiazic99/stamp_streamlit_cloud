"""Cloud path adapter for the GTEx v10 SMOTE Panel Explorer.

The offline project stores generated SMOTE sets under ``output/v10``. That
directory is intentionally ignored by the Streamlit Cloud repository, so the
deployable copy stores only the two required mapped datasets under
``gui/data/smote_switching_genes``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from smote_panel_data import (
    SMOTE_POLICIES,
    SWITCHING_BRACKETS,
    PROJECT_ROOT,
    build_panel_results,
    load_interface_alias_map,
    parse_panel,
)
from smote_panel_data import available_tissues as _available_tissues
from smote_panel_data import load_policy_events as _load_policy_events


DEFAULT_EXPORT_ROOT = (
    PROJECT_ROOT / "gui" / "data" / "smote_switching_genes"
)


def available_tissues(
    policies: Iterable[str],
    export_root: Path = DEFAULT_EXPORT_ROOT,
) -> list[str]:
    return _available_tissues(policies, export_root)


def load_policy_events(
    policy: str,
    tissues: Iterable[str],
    export_root: Path = DEFAULT_EXPORT_ROOT,
) -> dict[tuple[str, str], tuple[str, ...]]:
    return _load_policy_events(policy, tissues, export_root)


__all__ = [
    "SMOTE_POLICIES",
    "SWITCHING_BRACKETS",
    "available_tissues",
    "build_panel_results",
    "load_interface_alias_map",
    "load_policy_events",
    "parse_panel",
]
