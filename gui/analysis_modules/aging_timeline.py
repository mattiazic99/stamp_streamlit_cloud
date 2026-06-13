"""Aging Timeline — the "when & where" map of switching genes.

A single cross-tissue view of how many candidate switching genes occur in each
tissue at each age interval: tissues on the rows, the five switching brackets
(30-39 … 70-79) on the columns, colour = number of switching genes in that
(tissue, bracket) cell. This is the temporal architecture of molecular ageing
the application note foregrounds ("when transcriptional state reorganization
becomes detectable across tissues").

All numbers come straight from the pre-computed ``*_sets.txt`` files (the same
counts shown elsewhere); nothing is recomputed.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_loader import get_available_tissues, get_sets_for_tissue, complete_mode

_BRACKETS = ["30-39", "40-49", "50-59", "60-69", "70-79"]


def _plotly_cfg(fn="aging_timeline"):
    return {
        "toImageButtonOptions": {"format": "png", "scale": 2, "filename": fn},
        "displayModeBar": True,
    }


def show():
    st.header("🗺️ Aging Timeline — switching genes across tissues & age")
    st.markdown(
        "When and where do age-related transcriptional shifts appear? Each cell "
        "is the number of candidate switching genes in a tissue at a given age "
        "interval. Counts come directly from the pre-computed `*_sets.txt`."
    )

    version = st.session_state.get("gtex_version", "v10")
    complete = complete_mode()
    if complete:
        st.caption("🧪 Complete age-bins mode (only tissues with all six brackets).")

    tissues = get_available_tissues(version, complete)
    if not tissues:
        st.info("No pre-computed datasets found for this version/mode.")
        return

    selected = st.multiselect(
        "Tissues to display:",
        tissues,
        default=tissues,
        help="By default all available tissues are shown.",
    )
    if not selected:
        st.warning("Select at least one tissue.")
        return

    metric = st.radio(
        "Cell value:",
        ["Absolute counts", "Row-normalised (% of tissue's switching genes)"],
        horizontal=True,
        key="at_metric",
    )

    # Build the tissue × bracket count matrix from the sets files.
    rows = []
    for t in selected:
        sets = get_sets_for_tissue(version, t, complete)
        rows.append([len(sets.get(b, [])) for b in _BRACKETS])
    counts = pd.DataFrame(rows, index=selected, columns=_BRACKETS)

    # Order tissues by their peak bracket then by total, so the "wave" is visible.
    peak_idx = counts.values.argmax(axis=1)
    totals = counts.sum(axis=1).values
    order = np.lexsort((-totals, peak_idx))
    counts = counts.iloc[order]

    if metric.startswith("Row"):
        z = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0) * 100
        ztext = z.round(0).astype(int).astype(str)
        colorbar_title = "% of tissue"
        hover_val = "%{z:.1f}% of tissue"
    else:
        z = counts
        ztext = counts.astype(int).astype(str)
        colorbar_title = "n switching"
        hover_val = "%{z} switching genes"

    height = max(360, 22 * len(counts) + 140)
    fig = go.Figure(data=go.Heatmap(
        z=z.values,
        x=_BRACKETS,
        y=list(counts.index),
        colorscale="YlOrRd",
        text=ztext.values,
        texttemplate="%{text}",
        textfont=dict(size=9),
        hovertemplate="<b>%{y}</b><br>Age %{x}<br>" + hover_val + "<extra></extra>",
        colorbar=dict(title=colorbar_title, thickness=14),
        xgap=1, ygap=1,
    ))
    fig.update_layout(
        title=dict(text=f"Switching genes by tissue and age interval — GTEx {version}",
                   font=dict(size=15), x=0.5),
        xaxis_title="Age interval",
        height=height,
        margin=dict(l=10, r=10, t=60, b=40),
        plot_bgcolor="white",
        yaxis=dict(autorange="reversed"),
    )
    st.plotly_chart(fig, use_container_width=True, key="at_heatmap", config=_plotly_cfg())

    # Aggregate across tissues: which age interval concentrates switching events?
    st.markdown("### 📈 Where the aging 'wave' concentrates")
    per_bracket_total = counts.sum(axis=0)
    bar = go.Figure(go.Bar(
        x=_BRACKETS, y=per_bracket_total.values,
        text=[f"<b>{int(v)}</b>" for v in per_bracket_total.values],
        textposition="outside",
        marker=dict(color=per_bracket_total.values, colorscale="YlOrRd",
                    line=dict(color="rgb(120,40,20)", width=1)),
        hovertemplate="Age %{x}<br>Total switching genes: <b>%{y}</b><extra></extra>",
    ))
    ymax = per_bracket_total.max() if len(per_bracket_total) else 1
    bar.update_layout(
        title=dict(text="Total switching genes per age interval (all selected tissues)",
                   font=dict(size=13), x=0.5),
        xaxis_title="Age interval", yaxis_title="Total switching genes",
        height=320, margin=dict(l=40, r=20, t=50, b=40), plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", range=[0, ymax * 1.2 + 3]),
    )
    st.plotly_chart(bar, use_container_width=True, key="at_bar", config=_plotly_cfg("aging_wave"))

    # Per-tissue peak interval table (candidate "transition window" per tissue).
    peak_bracket = counts.idxmax(axis=1)
    peak_count = counts.max(axis=1)
    peak_df = pd.DataFrame({
        "Tissue": counts.index,
        "Peak age interval": peak_bracket.values,
        "Switching genes at peak": peak_count.values.astype(int),
        "Total switching genes": counts.sum(axis=1).values.astype(int),
    })
    st.markdown("### 🏔️ Peak switching interval per tissue")
    st.caption(
        "The age interval where each tissue shows the most switching events — a "
        "candidate window of transcriptional reorganisation (exploratory)."
    )
    st.dataframe(peak_df, use_container_width=True, hide_index=True)

    st.download_button(
        "⬇️ Download tissue × age count matrix (CSV)",
        counts.to_csv().encode("utf-8"),
        file_name=f"aging_timeline_{version}{'_complete' if complete else ''}.csv",
        mime="text/csv",
    )
