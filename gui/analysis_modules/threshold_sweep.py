"""Threshold Sweep Analysis Page.

Visualizes how the number of switching genes varies across different
tau thresholds for each tissue, using pre-computed sweep data.
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from data_loader import (
    get_threshold_sweep,
    get_stability_scores,
    get_available_tissues,
    safe_to_display,
    complete_mode,
)


# ── Plotly chart config ──
def _plotly_cfg(filename="chart"):
    return {
        "toImageButtonOptions": {"format": "png", "scale": 2, "filename": filename.replace(".png", "")},
        "displayModeBar": True,
    }

def _download_plotly_as_png(plotly_fig, filename):
    """Render a working Download PNG button using Plotly.js from CDN."""
    import streamlit.components.v1 as _components
    safe_name = filename.replace(".png", "").replace("'", r"\'")
    fig_json = plotly_fig.to_json().replace("</", r"<\/")
    html = (
        '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
        '<div id="hc" style="position:absolute;left:-9999px;width:1200px;height:600px;"></div>'
        '<div style="background:linear-gradient(90deg,#e8f4fd,#f0f7ff);border:1px solid #b8daff;'
        'border-radius:8px;padding:8px 14px;display:flex;align-items:center;gap:10px;">'
        '<button id="dlbtn" disabled style="background:#0d6efd;color:white;border:none;'
        'border-radius:6px;padding:6px 16px;cursor:pointer;font-size:13px;font-weight:600;'
        'white-space:nowrap;opacity:0.6;">Loading\\u2026</button>'
        '<span style="color:#555;font-size:12px;">' + filename + '</span></div>'
        '<script>'
        'var _fig=' + fig_json + ';'
        'function _init(){if(typeof Plotly==="undefined"){setTimeout(_init,100);return;}'
        'var b=document.getElementById("dlbtn");b.textContent="\\u2b07\\ufe0f Download PNG";'
        'b.disabled=false;b.style.opacity="1";'
        'b.onclick=function(){'
        'b.textContent="Generating\\u2026";b.disabled=true;'
        'Plotly.newPlot("hc",_fig.data,_fig.layout).then(function(){'
        'return Plotly.downloadImage("hc",{format:"png",scale:2,filename:"' + safe_name + '"});'
        '}).then(function(){Plotly.purge("hc");b.textContent="\\u2b07\\ufe0f Download PNG";b.disabled=false;});'
        '};}_init();'
        '</script>'
    )
    _components.html(html, height=50)


def show():
    """Threshold Sweep Analysis Page."""

    version = st.session_state.get("gtex_version", "v10")
    complete = complete_mode()

    st.header("📊 Threshold Sweep Analysis")
    st.markdown(
        f"Explore how switching gene counts change with threshold τ "
        f"across tissues. **GTEx {version}**"
        + (" · _complete age bins_" if complete else "")
    )

    # ── Load data ──────────────────────────────────────────────────
    try:
        sweep_df = get_threshold_sweep(version, complete)
        stability_df = get_stability_scores(version, complete)
    except FileNotFoundError as e:
        st.error(f"❌ Pre-computed data not found for {version}: {e}")
        return

    # Available tissues (from sweep data, display names)
    raw_tissues = sorted(sweep_df["tissue"].unique())
    tissue_display_map = {t: safe_to_display(t) for t in raw_tissues}
    display_names = [tissue_display_map[t] for t in raw_tissues]

    # ── Tissue selector ────────────────────────────────────────────
    st.markdown("""
    <div class="analysis-section">
        <h3>🧪 Select Tissues</h3>
    </div>
    """, unsafe_allow_html=True)

    selected_display = st.multiselect(
        "Select tissues to display:",
        display_names,
        default=display_names[:5],
        key="ts_tissues",
        help="Choose tissues to include in the threshold sweep plot",
    )

    if not selected_display:
        st.warning("⚠️ Please select at least one tissue.")
        return

    # Map back to safe names for filtering
    display_to_safe_map = {v: k for k, v in tissue_display_map.items()}
    selected_safe = [display_to_safe_map[d] for d in selected_display]

    # ── Summary metrics ────────────────────────────────────────────
    default_tau = 0.5
    default_data = sweep_df[
        (sweep_df["tau"] == default_tau) & (sweep_df["tissue"].isin(selected_safe))
    ]

    col1, col2, col3 = st.columns(3)
    with col1:
        total_switching = int(default_data["n_switching"].sum()) if not default_data.empty else 0
        st.markdown(f"""
        <div class="metric-card">
            <h3>{total_switching:,}</h3>
            <p>Total Switching Genes<br>(τ = {default_tau})</p>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(selected_display)}</h3>
            <p>Selected Tissues</p>
        </div>
        """, unsafe_allow_html=True)
    with col3:
        n_taus = sweep_df["tau"].nunique()
        st.markdown(f"""
        <div class="metric-card">
            <h3>{n_taus}</h3>
            <p>Threshold Values</p>
        </div>
        """, unsafe_allow_html=True)

    # ── Line chart: n_switching vs τ ───────────────────────────────
    st.markdown("""
    <div class="analysis-section">
        <h2>📈 Switching Gene Count vs Threshold (τ)</h2>
    </div>
    """, unsafe_allow_html=True)

    filtered = sweep_df[sweep_df["tissue"].isin(selected_safe)].copy()
    filtered["tissue_display"] = filtered["tissue"].map(tissue_display_map)

    # Choose a colorscale
    colors = px.colors.qualitative.Dark24
    if len(selected_safe) > len(colors):
        colors = colors * (len(selected_safe) // len(colors) + 1)

    fig = go.Figure()
    for idx, tissue_safe in enumerate(selected_safe):
        tissue_data = filtered[filtered["tissue"] == tissue_safe].sort_values("tau")
        display_name = tissue_display_map[tissue_safe]
        fig.add_trace(go.Scatter(
            x=tissue_data["tau"],
            y=tissue_data["n_switching"],
            mode="lines+markers",
            name=display_name,
            line=dict(color=colors[idx], width=2),
            marker=dict(size=4),
            hovertemplate=(
                f"<b>{display_name}</b><br>"
                "τ = <b>%{x:.2f}</b><br>"
                "Switching genes: <b>%{y:,}</b>"
                "<extra></extra>"
            ),
        ))

    fig.update_layout(
        title=dict(
            text=f"Switching Gene Count vs Threshold — GTEx {version}",
            font=dict(size=15), x=0.5,
        ),
        xaxis_title="Threshold (τ)",
        yaxis_title="Number of Switching Genes",
        height=500,
        margin=dict(l=60, r=30, t=70, b=60),
        plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False),
        xaxis=dict(gridcolor="rgba(0,0,0,0.08)"),
        legend=dict(
            font=dict(size=10),
            itemsizing="constant",
        ),
    )
    st.plotly_chart(fig, use_container_width=True, key="ts_line_chart", config=_plotly_cfg())
    _download_plotly_as_png(fig, f"threshold_sweep_{version}.png")

    # ── Stability scores table ─────────────────────────────────────
    st.markdown("""
    <div class="analysis-section">
        <h2>🏆 Stability Scores</h2>
        <p>Tissues ranked by stability of switching gene count across τ values.
        Lower CV = more stable.</p>
    </div>
    """, unsafe_allow_html=True)

    # Prepare table
    stab = stability_df.copy()
    stab.index.name = "tissue"
    stab = stab.reset_index()
    stab["display_name"] = stab["tissue"].map(tissue_display_map)
    stab = stab.sort_values("stability_cv", ascending=True)
    stab["rank"] = range(1, len(stab) + 1)

    display_stab = stab[["rank", "display_name", "stability_cv"]].copy()
    display_stab.columns = ["Rank", "Tissue", "Stability CV"]

    st.dataframe(
        display_stab,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Stability CV": st.column_config.NumberColumn(format="%.4f"),
        },
    )

    # Download
    csv_stab = display_stab.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="⬇️ Download Stability Scores (CSV)",
        data=csv_stab,
        file_name=f"stability_scores_{version}.csv",
        mime="text/csv",
        key="dl_stability",
    )

    # ── Per-bracket breakdown ──────────────────────────────────────
    if st.checkbox("🔬 Show per-bracket breakdown", value=False, key="ts_bracket"):
        st.markdown("### 📅 Switching Genes by Age Bracket")

        bracket_cols = [c for c in sweep_df.columns if c.startswith("n_") and c != "n_switching"]
        tau_select = st.slider(
            "Select τ value:",
            float(sweep_df["tau"].min()),
            float(sweep_df["tau"].max()),
            default_tau,
            step=0.05,
            key="ts_tau_slider",
        )

        tau_data = filtered[filtered["tau"].abs().sub(tau_select).lt(0.001)]
        if not tau_data.empty:
            bracket_fig = go.Figure()
            for col in bracket_cols:
                bracket_label = col.replace("n_", "")
                bracket_fig.add_trace(go.Bar(
                    name=bracket_label,
                    x=tau_data["tissue_display"],
                    y=tau_data[col],
                    hovertemplate=(
                        "<b>%{x}</b><br>"
                        f"Age {bracket_label}: <b>%{{y:,}}</b>"
                        "<extra></extra>"
                    ),
                ))

            bracket_fig.update_layout(
                barmode="stack",
                title=dict(
                    text=f"Switching Genes by Bracket (τ = {tau_select:.2f})",
                    font=dict(size=14), x=0.5,
                ),
                xaxis_title="Tissue",
                yaxis_title="Number of Switching Genes",
                height=450,
                margin=dict(l=60, r=30, t=70, b=100),
                plot_bgcolor="white",
                yaxis=dict(gridcolor="rgba(0,0,0,0.08)"),
                xaxis=dict(tickangle=45),
                legend=dict(orientation="h", yanchor="bottom", y=1.02,
                            xanchor="center", x=0.5),
            )
            st.plotly_chart(bracket_fig, use_container_width=True,
                            key="ts_bracket_chart", config=_plotly_cfg())
            _download_plotly_as_png(bracket_fig, f"bracket_breakdown_{version}_tau{tau_select:.2f}.png")

    # ── Download full sweep data ───────────────────────────────────
    st.markdown("""
    <div class="download-section">
        <h4>📥 Download Sweep Data</h4>
    </div>
    """, unsafe_allow_html=True)

    csv_sweep = filtered.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="⬇️ Download Filtered Sweep Results (CSV)",
        data=csv_sweep,
        file_name=f"sweep_results_{version}_filtered.csv",
        mime="text/csv",
        key="dl_sweep",
    )
