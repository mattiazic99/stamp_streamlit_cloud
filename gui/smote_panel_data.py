"""Pure data helpers for the GTEx v10 SMOTE Panel Explorer."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable


SWITCHING_BRACKETS = ("30-39", "40-49", "50-59", "60-69", "70-79")
SMOTE_POLICIES = ("fixed-floor", "adaptive-second-smallest")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_EXPORT_ROOT = (
    PROJECT_ROOT
    / "output"
    / "v10"
    / "statistical_validation_cli"
    / "runs"
    / "smote-switching-gene-lists"
    / "smote_switching_genes"
)
DEFAULT_GENE_MAP = PROJECT_ROOT / "gui" / "data" / "all_genes.txt"


def parse_panel(raw: str) -> list[str]:
    """Parse newline/comma/semicolon/space separated genes without duplicates."""
    tokens = re.split(r"[\s,;]+", raw.strip())
    seen: set[str] = set()
    panel: list[str] = []
    for token in tokens:
        cleaned = token.strip().upper()
        if cleaned and cleaned not in seen:
            seen.add(cleaned)
            panel.append(cleaned)
    return panel


def load_interface_alias_map(path: Path = DEFAULT_GENE_MAP) -> dict[str, str]:
    """Return every accepted symbol/Ensembl alias -> mapped interface token."""
    if not path.is_file():
        raise FileNotFoundError(f"Gene map not found: {path}")

    aliases: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if "(" not in line or ")" not in line:
                continue
            versioned_ensembl = line.split("(", 1)[0].strip()
            symbol = line.split("(", 1)[1].replace(")", "").strip()
            if not versioned_ensembl or not symbol:
                continue
            versionless_ensembl = versioned_ensembl.split(".", 1)[0]
            mapped_token = symbol.upper()
            aliases.setdefault(symbol.upper(), mapped_token)
            aliases[versionless_ensembl.upper()] = mapped_token
            aliases[versioned_ensembl.upper()] = mapped_token
    if not aliases:
        raise ValueError(f"No gene aliases parsed from: {path}")
    return aliases


def mapped_policy_directory(
    policy: str,
    export_root: Path = DEFAULT_EXPORT_ROOT,
) -> Path:
    if policy not in SMOTE_POLICIES:
        raise ValueError(
            f"Unknown SMOTE policy {policy!r}; expected {SMOTE_POLICIES}"
        )
    return export_root / f"{policy}-mapped"


def available_tissues(
    policies: Iterable[str],
    export_root: Path = DEFAULT_EXPORT_ROOT,
) -> list[str]:
    """Return safe tissue names present in every selected SMOTE policy."""
    tissue_sets: list[set[str]] = []
    for policy in policies:
        directory = mapped_policy_directory(policy, export_root)
        if not directory.is_dir():
            raise FileNotFoundError(f"Mapped SMOTE directory not found: {directory}")
        tissues = {
            path.name.removesuffix("_sets_stamp_mapped.txt")
            for path in directory.glob("*_sets_stamp_mapped.txt")
        }
        if not tissues:
            raise FileNotFoundError(f"No mapped SMOTE files found in: {directory}")
        tissue_sets.append(tissues)
    return sorted(set.intersection(*tissue_sets)) if tissue_sets else []


def load_policy_events(
    policy: str,
    tissues: Iterable[str],
    export_root: Path = DEFAULT_EXPORT_ROOT,
) -> dict[tuple[str, str], tuple[str, ...]]:
    """Load mapped events as ``(gene, tissue_safe) -> ordered brackets``."""
    directory = mapped_policy_directory(policy, export_root)
    events: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    bracket_order = {bracket: index for index, bracket in enumerate(SWITCHING_BRACKETS)}

    for tissue in tissues:
        path = directory / f"{tissue}_sets_stamp_mapped.txt"
        if not path.is_file():
            raise FileNotFoundError(f"Mapped SMOTE tissue file not found: {path}")
        lines = path.read_text(encoding="utf-8").splitlines()
        if len(lines) != len(SWITCHING_BRACKETS):
            raise ValueError(
                f"{path} has {len(lines)} lines; "
                f"expected {len(SWITCHING_BRACKETS)}"
            )
        for bracket, line in zip(SWITCHING_BRACKETS, lines):
            for gene in line.split():
                events[(gene.upper(), tissue)].add(bracket)

    return {
        key: tuple(sorted(brackets, key=bracket_order.__getitem__))
        for key, brackets in events.items()
    }


def build_panel_results(
    panel: Iterable[str],
    tissues: Iterable[str],
    events_by_policy: dict[
        str,
        dict[tuple[str, str], tuple[str, ...]],
    ],
    aliases: dict[str, str],
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    """Return long event records and one resolution row per queried gene."""
    tissue_list = list(tissues)
    event_tokens = {
        gene
        for policy_events in events_by_policy.values()
        for gene, _tissue in policy_events
    }
    event_records: list[dict[str, str]] = []
    resolution_records: list[dict[str, object]] = []

    for query in panel:
        query_upper = query.upper()
        versionless_query = query_upper.split(".", 1)[0]
        canonical_gene = aliases.get(
            query_upper,
            aliases.get(versionless_query, query_upper),
        )
        annotated = (
            query_upper in aliases
            or versionless_query in aliases
            or canonical_gene in event_tokens
        )
        has_event = False
        for policy, policy_events in events_by_policy.items():
            for tissue in tissue_list:
                for bracket in policy_events.get(
                    (canonical_gene, tissue),
                    (),
                ):
                    has_event = True
                    event_records.append(
                        {
                            "Query": query_upper,
                            "Gene": canonical_gene,
                            "Dataset": policy,
                            "Tissue": tissue,
                            "Age bracket": bracket,
                        }
                    )
        resolution_records.append(
            {
                "Query": query_upper,
                "Gene": canonical_gene,
                "Annotated": annotated,
                "Has event": has_event,
            }
        )

    return event_records, resolution_records
