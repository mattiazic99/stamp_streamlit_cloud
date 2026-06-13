"""Migration Paths Analysis Page.

Visualizes gene migration patterns between tissues across age brackets
using pre-computed migration path data.
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_loader import get_migration_paths, safe_to_display, complete_mode


def _plotly_cfg(fn="chart"):
    return {"toImageButtonOptions": {"format": "png", "scale": 2, "filename": fn.replace(".png", "")}, "displayModeBar": True}


def show():
    version = st.session_state.get("gtex_version", "v10")
    complete = complete_mode()
    st.header("🔀 Migration Paths Analysis")
    st.markdown(
        f"Explore gene migration patterns between tissues. **GTEx {version}**"
        + (" · _complete age bins_" if complete else "")
    )

    try:
        migration_df = get_migration_paths(version, complete)
    except FileNotFoundError as e:
        st.error(f"❌ Data not found for {version}: {e}")
        return
    if migration_df.empty:
        st.info("ℹ️ No migration paths found.")
        return

    migration_df["source_display"] = migration_df["source_tissue"].apply(safe_to_display)
    migration_df["target_display"] = migration_df["target_tissue"].apply(safe_to_display)

    # Metrics
    st.markdown('<div class="analysis-section"><h2>📊 Migration Overview</h2></div>', unsafe_allow_html=True)
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.markdown(f'<div class="metric-card"><h3>{len(migration_df)}</h3><p>Migration Paths</p></div>', unsafe_allow_html=True)
    with c2:
        st.markdown(f'<div class="metric-card"><h3>{int(migration_df["n_genes"].sum()):,}</h3><p>Total Gene Events</p></div>', unsafe_allow_html=True)
    with c3:
        st.markdown(f'<div class="metric-card"><h3>{migration_df["source_tissue"].nunique()}</h3><p>Source Tissues</p></div>', unsafe_allow_html=True)
    with c4:
        st.markdown(f'<div class="metric-card"><h3>{migration_df["target_tissue"].nunique()}</h3><p>Target Tissues</p></div>', unsafe_allow_html=True)

    # Filters
    st.markdown('<div class="analysis-section"><h3>🔍 Filter Migration Paths</h3></div>', unsafe_allow_html=True)
    src_opts = ["All"] + sorted(migration_df["source_display"].unique())
    tgt_opts = ["All"] + sorted(migration_df["target_display"].unique())
    c1, c2 = st.columns(2)
    with c1:
        sel_src = st.selectbox("🔹 Source Tissue:", src_opts, key="mp_source")
    with c2:
        sel_tgt = st.selectbox("🔸 Target Tissue:", tgt_opts, key="mp_target")

    filtered = migration_df.copy()
    if sel_src != "All":
        filtered = filtered[filtered["source_display"] == sel_src]
    if sel_tgt != "All":
        filtered = filtered[filtered["target_display"] == sel_tgt]

    st.markdown(f"### 📋 Migration Paths ({len(filtered)} results)")
    disp = filtered[["source_display", "source_bracket", "target_display", "target_bracket", "n_genes"]].copy()
    disp.columns = ["Source Tissue", "Source Bracket", "Target Tissue", "Target Bracket", "N Genes"]
    disp = disp.sort_values("N Genes", ascending=False)
    st.dataframe(disp, use_container_width=True, hide_index=True)
    st.download_button("⬇️ Download Filtered Table (CSV)", disp.to_csv(index=False).encode("utf-8"),
                        f"migration_paths_{version}_filtered.csv", "text/csv", key="dl_mig_table")

    # Heatmap
    st.markdown(f'<div class="analysis-section"><h2>🔥 Migration Heatmap</h2>'
                f'<p>Migrating genes per source → target pair (all brackets).</p></div>', unsafe_allow_html=True)
    agg = migration_df.groupby(["source_display", "target_display"])["n_genes"].sum().reset_index()
    all_t = sorted(set(agg["source_display"]) | set(agg["target_display"]))
    pivot = agg.pivot_table(index="source_display", columns="target_display", values="n_genes", fill_value=0)
    pivot = pivot.reindex(index=all_t, columns=all_t, fill_value=0)

    hover = [[f"<b>{s}</b> → <b>{t}</b><br>Genes: {int(pivot.loc[s, t])}" for t in pivot.columns] for s in pivot.index]
    hm = go.Figure(data=go.Heatmap(z=pivot.values, x=list(pivot.columns), y=list(pivot.index),
                                    colorscale="YlOrRd", hovertext=hover,
                                    hovertemplate="%{hovertext}<extra></extra>",
                                    colorbar=dict(title="Genes", thickness=14, len=0.75), xgap=1, ygap=1))
    side = min(700, max(400, len(all_t) * 20 + 150))
    hm.update_layout(title=dict(text=f"Gene Migration Heatmap — GTEx {version}", font=dict(size=14), x=0.5),
                     height=side, margin=dict(l=120, r=40, t=60, b=120),
                     yaxis=dict(autorange="reversed", title="Source"), xaxis=dict(title="Target", tickangle=45),
                     plot_bgcolor="white")
    st.plotly_chart(hm, use_container_width=True, key="mp_heatmap", config=_plotly_cfg())

    # Top pairs
    st.markdown("### 🏆 Top Migration Pairs")
    top = agg.sort_values("n_genes", ascending=False).head(10).rename(
        columns={"source_display": "Source", "target_display": "Target", "n_genes": "Total Migrating Genes"})
    st.dataframe(top, use_container_width=True, hide_index=True)
