from __future__ import annotations

from gui.smote_panel_data import (
    available_tissues,
    build_panel_results,
    load_interface_alias_map,
    load_policy_events,
    parse_panel,
)


def _write_sets(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_parse_panel_accepts_common_separators_and_deduplicates() -> None:
    assert parse_panel("APP, mapT\nAPP; ENSG000001.2") == [
        "APP",
        "MAPT",
        "ENSG000001.2",
    ]


def test_available_tissues_returns_policy_intersection(tmp_path) -> None:
    for tissue in ("Liver", "Lung"):
        _write_sets(
            tmp_path
            / "fixed-floor-mapped"
            / f"{tissue}_sets_stamp_mapped.txt",
            ["", "", "", "", ""],
        )
    _write_sets(
        tmp_path
        / "adaptive-second-smallest-mapped"
        / "Liver_sets_stamp_mapped.txt",
        ["", "", "", "", ""],
    )

    assert available_tissues(
        ("fixed-floor", "adaptive-second-smallest"),
        tmp_path,
    ) == ["Liver"]


def test_load_events_preserves_multiple_brackets_for_same_symbol(tmp_path) -> None:
    _write_sets(
        tmp_path
        / "fixed-floor-mapped"
        / "Liver_sets_stamp_mapped.txt",
        ["APP", "MAPT APP", "", "", ""],
    )

    events = load_policy_events("fixed-floor", ["Liver"], tmp_path)
    assert events[("APP", "Liver")] == ("30-39", "40-49")
    assert events[("MAPT", "Liver")] == ("40-49",)


def test_aliases_and_panel_results_support_symbols_and_ensembl(tmp_path) -> None:
    mapping_path = tmp_path / "all_genes.txt"
    mapping_path.write_text(
        "ENSG000001.12 (APP)\n"
        "ENSG000002.5 (MAPT)\n",
        encoding="utf-8",
    )
    aliases = load_interface_alias_map(mapping_path)
    events = {
        "fixed-floor": {
            ("APP", "Liver"): ("30-39",),
            ("MAPT", "Lung"): ("50-59",),
        }
    }

    records, resolution = build_panel_results(
        ["APP", "ENSG000002.9", "UNKNOWN"],
        ["Liver", "Lung"],
        events,
        aliases,
    )

    assert records == [
        {
            "Query": "APP",
            "Gene": "APP",
            "Dataset": "fixed-floor",
            "Tissue": "Liver",
            "Age bracket": "30-39",
        },
        {
            "Query": "ENSG000002.9",
            "Gene": "MAPT",
            "Dataset": "fixed-floor",
            "Tissue": "Lung",
            "Age bracket": "50-59",
        },
    ]
    assert resolution == [
        {
            "Query": "APP",
            "Gene": "APP",
            "Annotated": True,
            "Has event": True,
        },
        {
            "Query": "ENSG000002.9",
            "Gene": "MAPT",
            "Annotated": True,
            "Has event": True,
        },
        {
            "Query": "UNKNOWN",
            "Gene": "UNKNOWN",
            "Annotated": False,
            "Has event": False,
        },
    ]
