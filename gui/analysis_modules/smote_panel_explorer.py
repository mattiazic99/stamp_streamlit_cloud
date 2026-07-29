"""Panel explorer for the two mapped GTEx v10 SMOTE switching datasets."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


GUI_DIR = Path(__file__).resolve().parent.parent
if str(GUI_DIR) not in sys.path:
    sys.path.insert(0, str(GUI_DIR))

from data_loader import safe_to_display  # noqa: E402
from smote_panel_data import (  # noqa: E402
    SMOTE_POLICIES,
    SWITCHING_BRACKETS,
    available_tissues,
    build_panel_results,
    load_interface_alias_map,
    load_policy_events,
    parse_panel,
)


_DEFAULT_PANEL = [
    "APP",
    "PSEN1",
    "MAPT",
    "APOE",
    "TREM2",
    "BIN1",
    "IL6",
    "TNF",
    "NFKB1",
    "STAT3",
]

_MODE_POLICIES = {
    "Compare fixed and adaptive": SMOTE_POLICIES,
    "Fixed floor only": ("fixed-floor",),
    "Adaptive only": ("adaptive-second-smallest",),
}

_POLICY_LABEL = {
    "fixed-floor": "Fixed floor",
    "adaptive-second-smallest": "Adaptive",
}

_NO_EVENT, _FIXED_ONLY, _ADAPTIVE_ONLY, _DIFFERENT, _SAME = range(5)
_STATUS_LABEL = {
    _FIXED_ONLY: "Fixed only",
    _ADAPTIVE_ONLY: "Adaptive only",
    _DIFFERENT: "Different bracket",
    _SAME: "Same event",
}
_STATUS_COLOR = {
    _NO_EVENT: "#eceff1",
    _FIXED_ONLY: "#64b5f6",
    _ADAPTIVE_ONLY: "#ba68c8",
    _DIFFERENT: "#ffb74d",
    _SAME: "#81c784",
}


@st.cache_data(ttl=600)
def _cached_tissues(policies: tuple[str, ...]) -> list[str]:
    return available_tissues(policies)


@st.cache_data(ttl=600)
def _cached_events(
    policy: str,
    tissues: tuple[str, ...],
) -> dict[tuple[str, str], tuple[str, ...]]:
    return load_policy_events(policy, tissues)


@st.cache_data(ttl=3600)
def _cached_aliases() -> dict[str, str]:
    return load_interface_alias_map()


def _plotly_config() -> dict:
    return {
        "toImageButtonOptions": {
            "format": "png",
            "scale": 2,
            "filename": "smote_panel_event_map",
        },
        "displayModeBar": True,
    }


def _comparison_cell(
    gene: str,
    tissue: str,
    events_by_policy: dict[
        str,
        dict[tuple[str, str], tuple[str, ...]],
    ],
) -> tuple[int, str]:
    fixed = events_by_policy.get("fixed-floor", {}).get((gene, tissue), ())
    adaptive = events_by_policy.get(
        "adaptive-second-smallest",
        {},
    ).get((gene, tissue), ())
    if not fixed and not adaptive:
        return _NO_EVENT, ""
    if fixed and adaptive:
        if fixed == adaptive:
            return _SAME, " / ".join(fixed)
        return (
            _DIFFERENT,
            f"Fixed: {' / '.join(fixed)}<br>Adaptive: {' / '.join(adaptive)}",
        )
    if fixed:
        return _FIXED_ONLY, f"Fixed: {' / '.join(fixed)}"
    return _ADAPTIVE_ONLY, f"Adaptive: {' / '.join(adaptive)}"


def _single_cell(
    gene: str,
    tissue: str,
    policy: str,
    events: dict[tuple[str, str], tuple[str, ...]],
) -> tuple[int, str]:
    brackets = events.get((gene, tissue), ())
    if not brackets:
        return 0, ""
    first_position = min(SWITCHING_BRACKETS.index(bracket) for bracket in brackets)
    return first_position + 1, " / ".join(brackets)


def _render_event_map(
    genes: list[str],
    tissues: list[str],
    selected_policies: tuple[str, ...],
    events_by_policy: dict[
        str,
        dict[tuple[str, str], tuple[str, ...]],
    ],
) -> None:
    display_tissues = [safe_to_display(tissue) for tissue in tissues]
    compare = len(selected_policies) == 2

    if compare:
        cells = [
            [_comparison_cell(gene, tissue, events_by_policy) for tissue in tissues]
            for gene in genes
        ]
        z = [[cell[0] for cell in row] for row in cells]
        text = [[cell[1] for cell in row] for row in cells]
        colorscale = [
            [0.0, _STATUS_COLOR[_NO_EVENT]],
            [0.2, _STATUS_COLOR[_NO_EVENT]],
            [0.2, _STATUS_COLOR[_FIXED_ONLY]],
            [0.4, _STATUS_COLOR[_FIXED_ONLY]],
            [0.4, _STATUS_COLOR[_ADAPTIVE_ONLY]],
            [0.6, _STATUS_COLOR[_ADAPTIVE_ONLY]],
            [0.6, _STATUS_COLOR[_DIFFERENT]],
            [0.8, _STATUS_COLOR[_DIFFERENT]],
            [0.8, _STATUS_COLOR[_SAME]],
            [1.0, _STATUS_COLOR[_SAME]],
        ]
        zmin, zmax = -0.5, 4.5
    else:
        policy = selected_policies[0]
        cells = [
            [
                _single_cell(gene, tissue, policy, events_by_policy[policy])
                for tissue in tissues
            ]
            for gene in genes
        ]
        z = [[cell[0] for cell in row] for row in cells]
        text = [[cell[1] for cell in row] for row in cells]
        bracket_colors = [
            "#eceff1",
            "#d7e9ff",
            "#acd2ff",
            "#7ab8ff",
            "#4d9ef5",
            "#1769aa",
        ]
        colorscale = []
        for index, color in enumerate(bracket_colors):
            start = index / len(bracket_colors)
            end = (index + 1) / len(bracket_colors)
            colorscale.extend([[start, color], [end, color]])
        zmin, zmax = -0.5, 5.5

    figure = go.Figure(
        go.Heatmap(
            z=z,
            x=display_tissues,
            y=genes,
            text=text,
            texttemplate="%{text}",
            textfont={"size": 10},
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            showscale=False,
            xgap=2,
            ygap=2,
            hovertemplate=(
                "<b>%{y}</b><br><b>%{x}</b><br>%{text}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        height=max(350, 45 * len(genes) + 150),
        margin={"l": 80, "r": 20, "t": 30, "b": 170},
        plot_bgcolor="white",
        xaxis={"side": "bottom", "tickangle": -45},
        yaxis={"autorange": "reversed"},
    )
    st.plotly_chart(
        figure,
        use_container_width=True,
        key="smote_panel_heatmap",
        config=_plotly_config(),
    )

    if compare:
        legend = " &nbsp; ".join(
            (
                f'<span style="background:{_STATUS_COLOR[code]};'
                'padding:2px 10px;border-radius:4px;">&nbsp;</span> '
                f"{label}"
            )
            for code, label in [
                (_SAME, "Same event"),
                (_DIFFERENT, "Different bracket"),
                (_FIXED_ONLY, "Fixed only"),
                (_ADAPTIVE_ONLY, "Adaptive only"),
                (_NO_EVENT, "No event"),
            ]
        )
        st.markdown(legend, unsafe_allow_html=True)
    else:
        st.caption(
            "Cells report the switching age bracket; empty cells have no "
            "switching event for that gene and tissue."
        )


def show() -> None:
    st.header("🧪 SMOTE Panel Explorer — GTEx v10")
    st.markdown(
        "Query the **mapped switching-gene files generated from the two SMOTE "
        "datasets**. Select all 50 tissues or a custom subset and identify "
        "the tissue and switching age bracket for every panel gene."
    )
    st.info(
        "This page is restricted to **GTEx v10** and uses the fixed, already "
        "generated files. It does not recompute SMOTE or run statistical tests."
    )

    mode = st.radio(
        "Dataset:",
        options=list(_MODE_POLICIES),
        horizontal=True,
        help=(
            "Compare both methods or query only the fixed-floor or adaptive "
            "SMOTE dataset."
        ),
    )
    selected_policies = _MODE_POLICIES[mode]

    try:
        all_tissues_safe = _cached_tissues(selected_policies)
        aliases = _cached_aliases()
    except Exception as error:
        st.error(f"Cannot load the SMOTE v10 datasets: {error}")
        return

    display_to_safe = {
        safe_to_display(tissue): tissue
        for tissue in all_tissues_safe
    }
    all_tissues_display = sorted(display_to_safe)

    c1, c2 = st.columns([1, 1])
    with c1:
        panel_raw = st.text_area(
            "Gene panel (symbols or Ensembl IDs):",
            value="\n".join(_DEFAULT_PANEL),
            height=220,
            help="One gene per line, or use spaces, commas, or semicolons.",
        )
    with c2:
        select_all = st.checkbox(
            f"Use all {len(all_tissues_display)} tissues",
            value=True,
        )
        if select_all:
            selected_tissues_display = all_tissues_display
            st.caption(
                f"All {len(selected_tissues_display)} complete v10 tissues "
                "are selected."
            )
        else:
            selected_tissues_display = st.multiselect(
                "Tissues to query:",
                options=all_tissues_display,
                default=all_tissues_display[:5],
            )

    panel = parse_panel(panel_raw)
    if not panel:
        st.warning("Enter at least one gene.")
        return
    if not selected_tissues_display:
        st.warning("Select at least one tissue.")
        return

    selected_tissues_safe = [
        display_to_safe[tissue]
        for tissue in selected_tissues_display
    ]
    tissue_tuple = tuple(selected_tissues_safe)
    try:
        events_by_policy = {
            policy: _cached_events(policy, tissue_tuple)
            for policy in selected_policies
        }
    except Exception as error:
        st.error(f"Cannot load the selected SMOTE files: {error}")
        return

    records, resolution = build_panel_results(
        panel,
        selected_tissues_safe,
        events_by_policy,
        aliases,
    )
    resolution_frame = pd.DataFrame(resolution)
    event_frame = pd.DataFrame(records)
    if not event_frame.empty:
        event_frame["Tissue"] = event_frame["Tissue"].map(safe_to_display)
        event_frame["Dataset"] = event_frame["Dataset"].map(_POLICY_LABEL)
        event_frame = event_frame.sort_values(
            ["Gene", "Tissue", "Dataset", "Age bracket"],
            kind="stable",
        ).reset_index(drop=True)

    present = resolution_frame.loc[
        resolution_frame["Has event"],
        "Gene",
    ].drop_duplicates().tolist()
    absent = resolution_frame.loc[
        resolution_frame["Annotated"] & ~resolution_frame["Has event"],
        "Gene",
    ].drop_duplicates().tolist()
    unresolved = resolution_frame.loc[
        ~resolution_frame["Annotated"],
        "Query",
    ].drop_duplicates().tolist()

    st.markdown("### Panel summary")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Genes queried", len(panel))
    m2.metric("With switching events", len(present))
    m3.metric("Annotated, no event", len(absent))
    m4.metric("Not in annotation", len(unresolved))

    if absent:
        st.caption(
            "Annotated but without a switching event in the selected "
            f"datasets/tissues: {', '.join(absent)}"
        )
    if unresolved:
        st.warning(
            "Genes not found in `all_genes.txt`: "
            f"{', '.join(unresolved)}"
        )
    if not present:
        st.info(
            "None of the resolved panel genes has a switching event in the "
            "selected datasets and tissues."
        )
        return

    st.markdown("### Tissue × gene event map")
    _render_event_map(
        present,
        selected_tissues_safe,
        selected_policies,
        events_by_policy,
    )

    st.markdown("### Switching events")
    st.dataframe(event_frame, use_container_width=True, hide_index=True)
    st.download_button(
        "⬇️ Download events CSV",
        event_frame.to_csv(index=False).encode("utf-8"),
        "smote_v10_panel_switching_events.csv",
        "text/csv",
        key="smote_panel_download",
    )

    st.info(
        "An absent event means only that the gene is not called switching "
        "under the selected SMOTE dataset, tissues, τ=0.5 and ε=0.01. "
        "SMOTE is an exploratory sensitivity analysis and synthetic samples "
        "are not independent biological donors."
    )
