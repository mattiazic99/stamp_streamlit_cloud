"""Switch direction explorer — low→high vs high→low.

The application note states STAMP "records both the direction (low→high or
high→low) and the age interval of the change", but the `*_sets.txt` files store
only which genes switch at each interval, not the direction. This page recovers
the direction from the pre-computed normalized matrices, applying exactly the
same binarisation and one-transition rule as the switching backend
(`stamp.switching.identify_switching_with_direction`), so the set of switching
genes is identical to the rest of the app — only the up/down split is added.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_loader import get_available_tissues, get_normalized_tissue, complete_mode
from stamp.switching import identify_switching_with_direction

_BRACKETS = ["30-39", "40-49", "50-59", "60-69", "70-79"]
_UP = "#d6604d"     # low → high (activation)
_DOWN = "#4393c3"   # high → low (silencing)


def _plotly_cfg(fn="switch_direction"):
    return {
        "toImageButtonOptions": {"format": "png", "scale": 2, "filename": fn},
        "displayModeBar": True,
    }


def show():
    st.header("🔼🔽 Switch direction — activation vs silencing with age")
    st.markdown(
        "Splits candidate switching genes by **direction**: low→high "
        "(*activation*) versus high→low (*silencing*), per age interval. "
        "Derived from the normalized matrices with the same rule as the pipeline."
    )

    version = st.session_state.get("gtex_version", "v10")
    complete = complete_mode()
    if complete:
        st.caption("🧪 Complete age-bins mode.")

    tissues = get_available_tissues(version, complete)
    if not tissues:
        st.info("No pre-computed datasets found for this version/mode.")
        return

    col1, col2 = st.columns([3, 2])
    with col1:
        selected = st.multiselect(
            "Tissues:", tissues, default=tissues[:1] or tissues,
            help="Counts are summed over the selected tissues.",
        )
    with col2:
        tau = st.slider("Threshold τ", 0.0, 1.0, 0.5, 0.05, key="sd_tau",
                        help="Default 0.5, as in the paper.")
    if not selected:
        st.warning("Select at least one tissue.")
        return

    # Aggregate up/down counts per bracket across the selected tissues.
    up = {b: 0 for b in _BRACKETS}
    down = {b: 0 for b in _BRACKETS}
    rows = []
    skipped = []
    for t in selected:
        try:
            df = get_normalized_tissue(version, t, complete)
        except FileNotFoundError:
            skipped.append(t)
            continue
        res = identify_switching_with_direction(df, tau)
        for b in _BRACKETS:
            u = sum(1 for _, d in res.get(b, []) if d == "up")
            dn = sum(1 for _, d in res.get(b, []) if d == "down")
            up[b] += u
            down[b] += dn
            rows.append({"Tissue": t, "Age interval": b, "Up (low→high)": u,
                         "Down (high→low)": dn, "Total": u + dn})

    if skipped:
        st.caption(f"({len(skipped)} tissue(s) without a normalized matrix skipped.)")

    up_v = [up[b] for b in _BRACKETS]
    down_v = [down[b] for b in _BRACKETS]
    tot_up, tot_down = sum(up_v), sum(down_v)

    c1, c2, c3 = st.columns(3)
    c1.markdown(f"<div class='metric-card'><h3>{tot_up}</h3><p>Activation<br>(low→high)</p></div>",
                unsafe_allow_html=True)
    c2.markdown(f"<div class='metric-card'><h3>{tot_down}</h3><p>Silencing<br>(high→low)</p></div>",
                unsafe_allow_html=True)
    ratio = (tot_up / tot_down) if tot_down else float("nan")
    c3.markdown(f"<div class='metric-card'><h3>{ratio:.2f}</h3><p>Up / Down<br>ratio</p></div>",
                unsafe_allow_html=True)

    # Grouped bar: up vs down per age interval.
    fig = go.Figure()
    fig.add_trace(go.Bar(x=_BRACKETS, y=up_v, name="low→high (activation)",
                         marker_color=_UP,
                         text=[f"<b>{v}</b>" for v in up_v], textposition="outside",
                         hovertemplate="Age %{x}<br>Activation: <b>%{y}</b><extra></extra>"))
    fig.add_trace(go.Bar(x=_BRACKETS, y=down_v, name="high→low (silencing)",
                         marker_color=_DOWN,
                         text=[f"<b>{v}</b>" for v in down_v], textposition="outside",
                         hovertemplate="Age %{x}<br>Silencing: <b>%{y}</b><extra></extra>"))
    ymax = max(up_v + down_v) if (up_v + down_v) else 1
    fig.update_layout(
        barmode="group",
        title=dict(text=f"Switch direction by age interval — GTEx {version} (τ={tau:g})",
                   font=dict(size=15), x=0.5),
        xaxis_title="Age interval", yaxis_title="Number of switching genes",
        height=400, margin=dict(l=50, r=20, t=60, b=40), plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", range=[0, ymax * 1.2 + 3]),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    )
    st.plotly_chart(fig, use_container_width=True, key="sd_bar", config=_plotly_cfg())

    if rows:
        det = pd.DataFrame(rows)
        st.markdown("### 📋 Per-tissue breakdown")
        st.dataframe(det, use_container_width=True, hide_index=True)
        st.download_button(
            "⬇️ Download direction breakdown (CSV)",
            det.to_csv(index=False).encode("utf-8"),
            file_name=f"switch_direction_{version}{'_complete' if complete else ''}_tau{tau:g}.csv",
            mime="text/csv",
        )
