import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from utils.parsing import parse_stamp_file, extract_tissue_name  # Import della funzione di pulizia
from components.downloads import create_csv_download, display_download_section

# Canonical age-group order
_AGE_ORDER = ["30–39", "40–49", "50–59", "60–69", "70–79"]


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
    fig_json = plotly_fig.to_json().replace("</", r"<\/")   # prevent HTML injection
    html = (
        '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
        '<div id="hc" style="position:absolute;left:-9999px;width:1200px;height:600px;"></div>'
        '<div style="background:linear-gradient(90deg,#e8f4fd,#f0f7ff);border:1px solid #b8daff;'
        'border-radius:8px;padding:8px 14px;display:flex;align-items:center;gap:10px;">'
        '<button id="dlbtn" disabled style="background:#0d6efd;color:white;border:none;'
        'border-radius:6px;padding:6px 16px;cursor:pointer;font-size:13px;font-weight:600;'
        'white-space:nowrap;opacity:0.6;">Loading\u2026</button>'
        '<span style="color:#555;font-size:12px;">' + filename + '</span></div>'
        '<script>'
        'var _fig=' + fig_json + ';'
        'function _init(){if(typeof Plotly==="undefined"){setTimeout(_init,100);return;}'
        'var b=document.getElementById("dlbtn");b.textContent="\u2b07\ufe0f Download PNG";'
        'b.disabled=false;b.style.opacity="1";'
        'b.onclick=function(){'
        'b.textContent="Generating\u2026";b.disabled=true;'
        'Plotly.newPlot("hc",_fig.data,_fig.layout).then(function(){'
        'return Plotly.downloadImage("hc",{format:"png",scale:2,filename:"' + safe_name + '"});'
        '}).then(function(){Plotly.purge("hc");b.textContent="\u2b07\ufe0f Download PNG";b.disabled=false;});'
        '};}_init();'
        '</script>'
    )
    _components.html(html, height=50)


def show():
    """Tissue Comparison Analysis Page"""
    
    st.header("🔄 Tissue Comparison Analysis")
    st.markdown("Compare gene expression patterns between two specific tissues across age groups.")
    
    age_groups = ["30–39", "40–49", "50–59", "60–69", "70–79"]

    # Track whether we are on pre-computed GTEx data (so we can read the
    # canonical Jaccard matrices instead of recomputing them client-side).
    preloaded = False
    version = None

    # ── Data source toggle ─────────────────────────────────────────
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
    from data_loader import get_available_tissues, get_sets_for_tissue, sets_to_gui_format, complete_mode

    # Pre-loaded (GTEx) mode removed by design: users always upload their files.
    data_source = "📂 Upload files"

    complete = complete_mode()

    if data_source == "📦 Pre-loaded (GTEx)":
        preloaded = True
        version = st.session_state.get("gtex_version", "v10")
        tissues = get_available_tissues(version, complete)

        st.markdown(f"""
        <div class="analysis-section">
            <h3>🧪 Select Two Tissues — GTEx {version}</h3>
        </div>
        """, unsafe_allow_html=True)

        col1, col2 = st.columns(2)
        with col1:
            tissue1 = st.selectbox("📊 First Tissue", tissues, index=0, key="tc_tissue1")
        with col2:
            default_idx = min(1, len(tissues) - 1)
            tissue2 = st.selectbox("📊 Second Tissue", tissues, index=default_idx, key="tc_tissue2")

        sets1 = get_sets_for_tissue(version, tissue1, complete)
        sets2 = get_sets_for_tissue(version, tissue2, complete)
        fasce1, counts1, df1 = sets_to_gui_format(sets1)
        fasce2, counts2, df2 = sets_to_gui_format(sets2)
        tissue1_name = tissue1
        tissue2_name = tissue2

    else:
        st.markdown("""
        <div class="analysis-section">
            <h3>📂 Upload Two Tissues for Comparison</h3>
        </div>
        """, unsafe_allow_html=True)

        col1, col2 = st.columns(2)
        with col1:
            file1 = st.file_uploader(
                "📊 First Tissue", 
                type=["txt"], 
                key="tissue1_comp",
                help="Upload first tissue gene switching data"
            )
        with col2:
            file2 = st.file_uploader(
                "📊 Second Tissue", 
                type=["txt"], 
                key="tissue2_comp",
                help="Upload second tissue gene switching data"
            )

        if not (file1 and file2):
            st.info("👆 Please upload both tissue files to start the comparison.")
            return

        fasce1, counts1, df1 = parse_stamp_file(file1, age_groups)
        fasce2, counts2, df2 = parse_stamp_file(file2, age_groups)
        tissue1_name = extract_tissue_name(file1.name)
        tissue2_name = extract_tissue_name(file2.name)
    
    # Age group selection
    st.markdown("### 🎯 Age Group Selection")
    selected_ages = st.multiselect(
        "Select age groups for comparison:",
        age_groups, 
        default=age_groups,
        help="Choose which age groups to include in the comparison"
    )
    
    if not selected_ages:
        st.warning("⚠️ Please select at least one age group.")
        return
    
    # Filter data
    df1_filtered = df1[df1["Age"].isin(selected_ages)]
    df2_filtered = df2[df2["Age"].isin(selected_ages)]
    
    counts1_filtered = [
        df1_filtered[df1_filtered["Age"] == age].shape[0]
        for age in age_groups if age in selected_ages
    ]
    counts2_filtered = [
        df2_filtered[df2_filtered["Age"] == age].shape[0]
        for age in age_groups if age in selected_ages
    ]
    
    # === COMPARISON METRICS ===
    st.markdown("""
    <div class="analysis-section">
        <h2>📊 Comparison Overview</h2>
    </div>
    """, unsafe_allow_html=True)
    
    # Calculate overall statistics
    total_genes_1 = len(set(df1["Gene"]))
    total_genes_2 = len(set(df2["Gene"]))
    common_genes = len(set(df1["Gene"]) & set(df2["Gene"]))
    unique_genes = len(set(df1["Gene"]) | set(df2["Gene"]))

    # ── Jaccard: prefer the PRE-COMPUTED backend matrices ───────────
    # In pre-loaded (GTEx) mode the canonical J_life / J_age values come
    # from output/{version}/jaccard/jaccard_{life,age}.csv — the same
    # matrices used by the paper and the v8-vs-v10 page. We read them here
    # instead of recomputing, so the GUI never shows a value that drifts
    # from the published one. In upload mode (no matrices) we fall back to
    # the local computation.
    j_life_pre = None
    j_age_pre = None
    if preloaded:
        try:
            from data_loader import get_jaccard_matrix, display_to_safe
            _s1, _s2 = display_to_safe(tissue1_name), display_to_safe(tissue2_name)
            _jl = get_jaccard_matrix(version, "life", complete)
            _ja = get_jaccard_matrix(version, "age", complete)
            if _s1 in _jl.index and _s2 in _jl.columns:
                j_life_pre = float(_jl.loc[_s1, _s2])
            if _s1 in _ja.index and _s2 in _ja.columns:
                j_age_pre = float(_ja.loc[_s1, _s2])
        except Exception:
            pass

    # Headline similarity = pre-computed J_life when available, else local.
    jaccard_sim = (
        j_life_pre
        if j_life_pre is not None
        else (common_genes / unique_genes if unique_genes > 0 else 0)
    )
    jaccard_is_precomputed = j_life_pre is not None
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{total_genes_1}</h3>
            <p>{tissue1_name}<br>Total Genes</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{total_genes_2}</h3>
            <p>{tissue2_name}<br>Total Genes</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{common_genes}</h3>
            <p>Shared<br>Genes</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{jaccard_sim:.2%}</h3>
            <p>Jaccard J<sub>life</sub><br>Similarity</p>
        </div>
        """, unsafe_allow_html=True)

    if jaccard_is_precomputed:
        _cap = f"J<sub>life</sub> = {j_life_pre:.4f}"
        if j_age_pre is not None:
            _cap += f" · J<sub>age</sub> = {j_age_pre:.4f}"
        st.caption(
            f"Pre-computed backend values ({_cap}) from "
            f"`output/{version}/jaccard/`, identical to the v8-vs-v10 page.",
            unsafe_allow_html=True,
        )
    elif not preloaded:
        st.caption(
            "J_life computed from the uploaded files "
            "(no pre-computed matrix available in upload mode)."
        )

    # === SIDE-BY-SIDE COMPARISON CHART (Plotly) ===
    st.markdown("### 📊 Gene Count Comparison by Age Group")
    
    # Tissue toggle
    show_tissue = st.radio(
        "Show:",
        options=["both", tissue1_name, tissue2_name],
        index=0,
        horizontal=True,
        key="tc_tissue_toggle",
    )

    # Ensure canonical order
    paired = list(zip(selected_ages, counts1_filtered, counts2_filtered))
    paired.sort(key=lambda t: _AGE_ORDER.index(t[0]) if t[0] in _AGE_ORDER else 99)
    ages_sorted = [t[0] for t in paired]
    c1_sorted = [t[1] for t in paired]
    c2_sorted = [t[2] for t in paired]

    n = len(ages_sorted)
    is_both = show_tissue == "both"
    bar_w = (0.22 if n >= 4 else 0.16) if is_both else (0.26 if n >= 4 else 0.18)

    comp_fig = go.Figure()

    if show_tissue in ("both", tissue1_name):
        total1 = sum(c1_sorted) or 1
        pct1 = [(c / total1 * 100) for c in c1_sorted]
        comp_fig.add_trace(go.Bar(
            name=tissue1_name,
            x=ages_sorted, y=c1_sorted,
            width=[bar_w] * n,
            marker=dict(color="rgba(100,170,230,0.85)",
                        line=dict(color="rgb(8,48,107)", width=1)),
            text=[f"<b>{c}</b>" for c in c1_sorted],
            textposition="outside",
            hovertemplate=(
                f"<b>{tissue1_name}</b><br>"
                "Age: <b>%{x}</b><br>"
                "Genes: <b>%{y}</b><br>"
                "Tissue %: <b>%{customdata:.1f}%</b>"
                "<extra></extra>"
            ),
            customdata=pct1,
        ))

    if show_tissue in ("both", tissue2_name):
        total2 = sum(c2_sorted) or 1
        pct2 = [(c / total2 * 100) for c in c2_sorted]
        comp_fig.add_trace(go.Bar(
            name=tissue2_name,
            x=ages_sorted, y=c2_sorted,
            width=[bar_w] * n,
            marker=dict(color="rgba(250,128,114,0.85)",
                        line=dict(color="darkred", width=1)),
            text=[f"<b>{c}</b>" for c in c2_sorted],
            textposition="outside",
            hovertemplate=(
                f"<b>{tissue2_name}</b><br>"
                "Age: <b>%{x}</b><br>"
                "Genes: <b>%{y}</b><br>"
                "Tissue %: <b>%{customdata:.1f}%</b>"
                "<extra></extra>"
            ),
            customdata=pct2,
        ))

    all_vals = c1_sorted + c2_sorted
    y_max = max(all_vals) if all_vals else 1
    y_upper = y_max * 1.20 + 3

    comp_fig.update_layout(
        barmode="group",
        title=dict(
            text=f"Gene Expression Comparison: {tissue1_name} vs {tissue2_name}",
            font=dict(size=15), x=0.5,
        ),
        xaxis_title="Age Group",
        yaxis_title="Number of Switching Genes",
        height=380,
        margin=dict(l=50, r=30, t=70, b=50),
        plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False, range=[0, y_upper]),
        xaxis=dict(gridcolor="rgba(0,0,0,0.08)", type="category",
                   categoryorder="array", categoryarray=_AGE_ORDER),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    )
    st.plotly_chart(comp_fig, use_container_width=True, key="tc_bar_comparison", config=_plotly_cfg())
    _download_plotly_as_png(comp_fig, f"comparison_{tissue1_name}_vs_{tissue2_name}.png")
    
    # === SHARED VS EXCLUSIVE ANALYSIS ===
    st.markdown("""
    <div class="analysis-section">
        <h2>🤝 Shared vs Exclusive Gene Analysis</h2>
    </div>
    """, unsafe_allow_html=True)
    
    # Calculate for whole lifespan
    genes_1_all = set(df1["Gene"])
    genes_2_all = set(df2["Gene"])
    shared_genes_all = sorted(genes_1_all & genes_2_all)
    exclusive_1_all = sorted(genes_1_all - genes_2_all)
    exclusive_2_all = sorted(genes_2_all - genes_1_all)
    
    # Pie chart (Plotly)
    sizes = [len(shared_genes_all), len(exclusive_1_all), len(exclusive_2_all)]
    labels_pie = [f'Shared ({len(shared_genes_all)})',
                  f'{tissue1_name} Exclusive ({len(exclusive_1_all)})',
                  f'{tissue2_name} Exclusive ({len(exclusive_2_all)})']
    colors_pie = ['#2ecc71', '#3498db', '#e74c3c']

    pie_fig = go.Figure(data=go.Pie(
        labels=labels_pie, values=sizes,
        marker=dict(colors=colors_pie, line=dict(color='white', width=2)),
        textinfo='percent+label', textposition='inside',
        hovertemplate='<b>%{label}</b><br>Genes: %{value}<br>Percentage: %{percent}<extra></extra>',
        hole=0.0,
    ))
    pie_fig.update_layout(
        title=dict(text=f"Gene Distribution: {tissue1_name} vs {tissue2_name}",
                   font=dict(size=14), x=0.5),
        height=350, margin=dict(l=30, r=30, t=60, b=30),
        showlegend=False,
    )
    _left, _c, _r = st.columns([1, 2, 1])
    with _c:
        st.plotly_chart(pie_fig, use_container_width=True, key="tc_pie_shared", config=_plotly_cfg())
    _download_plotly_as_png(pie_fig, f"pie_shared_exclusive_{tissue1_name}_vs_{tissue2_name}.png")

    # Bar chart (Plotly)
    categories = ['Shared', f'{tissue1_name}\nExclusive', f'{tissue2_name}\nExclusive']
    values = sizes
    bar_colors = ['#2ecc71', '#3498db', '#e74c3c']

    bar_fig = go.Figure()
    bar_fig.add_trace(go.Bar(
        x=categories, y=values,
        marker=dict(color=bar_colors, line=dict(color="black", width=1.2)),
        text=[f"<b>{v}</b>" for v in values],
        textposition="outside",
        hovertemplate="<b>%{x}</b><br>Genes: <b>%{y}</b><extra></extra>",
    ))
    y_max_bar = max(values) if values else 1
    bar_fig.update_layout(
        title=dict(text="Gene Count by Category", font=dict(size=15), x=0.5),
        yaxis_title="Number of Genes",
        height=380,
        margin=dict(l=50, r=30, t=70, b=50),
        plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False,
                   range=[0, y_max_bar * 1.20 + 3]),
    )
    _, col_bar, _ = st.columns([1, 2, 1])
    with col_bar:
        st.plotly_chart(bar_fig, use_container_width=True, key="tc_bar_categories", config=_plotly_cfg())
        _download_plotly_as_png(bar_fig, f"shared_exclusive_{tissue1_name}_vs_{tissue2_name}.png")
    
    # === AGE-SPECIFIC ANALYSIS ===
    st.markdown("### 📅 Age-Specific Shared Genes Analysis")

    # The per-bracket Jaccard values below are the COMPONENTS of the backend
    # J_age metric (which is their mean over the brackets populated in both
    # tissues). The backend stores only the aggregate, so the per-bracket
    # breakdown is computed here; in pre-loaded mode we show the canonical
    # aggregate alongside it for reference.
    if preloaded and j_age_pre is not None:
        st.caption(
            f"Backend J<sub>age</sub> (pre-computed) = **{j_age_pre:.4f}** — "
            "mean of the per-bracket Jaccard values below over the brackets "
            "populated in both tissues.",
            unsafe_allow_html=True,
        )

    # Calculate shared genes for each age group
    age_shared_data = []
    for age in age_groups:
        genes_1_age = set(df1[df1["Age"] == age]["Gene"])
        genes_2_age = set(df2[df2["Age"] == age]["Gene"])
        shared_age = len(genes_1_age & genes_2_age)
        total_1_age = len(genes_1_age)
        total_2_age = len(genes_2_age)
        union_age = len(genes_1_age | genes_2_age)
        jaccard_age = shared_age / union_age if union_age > 0 else 0
        
        age_shared_data.append({
            'Age Group': age,
            f'{tissue1_name} Genes': total_1_age,
            f'{tissue2_name} Genes': total_2_age,
            'Shared Genes': shared_age,
            'Jaccard Similarity': jaccard_age
        })
    
    df_age_analysis = pd.DataFrame(age_shared_data)
    
    # Display table with download / search / fullscreen toolbar
    st.dataframe(
        df_age_analysis,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Jaccard Similarity": st.column_config.NumberColumn(format="%.4f"),
        },
    )

    # Explicit download button for the table
    csv_table = df_age_analysis.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="⬇️ Download Age-Specific Analysis (CSV)",
        data=csv_table,
        file_name=f"age_specific_analysis_{tissue1_name}_vs_{tissue2_name}.csv",
        mime="text/csv",
        key="dl_age_table",
    )
    
    # Line chart of Jaccard similarity across ages (Plotly)
    jaccard_vals = list(df_age_analysis['Jaccard Similarity'])
    age_labels = list(df_age_analysis['Age Group'])
    y_top = max(jaccard_vals) * 1.15 if jaccard_vals else 1

    line_fig = go.Figure()
    line_fig.add_trace(go.Scatter(
        x=age_labels, y=jaccard_vals,
        mode='lines+markers+text',
        line=dict(color='#9b59b6', width=3),
        marker=dict(size=9, color='#9b59b6'),
        fill='tozeroy', fillcolor='rgba(155,89,182,0.18)',
        text=[f'{v:.3f}' for v in jaccard_vals],
        textposition='top center',
        hovertemplate='<b>%{x}</b><br>Jaccard: <b>%{y:.3f}</b><extra></extra>',
    ))
    line_fig.update_layout(
        title=dict(text=f"Jaccard Similarity: {tissue1_name} vs {tissue2_name}",
                   font=dict(size=14), x=0.5),
        xaxis_title='Age Group', yaxis_title='Jaccard Similarity',
        height=350, margin=dict(l=50, r=30, t=60, b=50),
        plot_bgcolor='white',
        yaxis=dict(gridcolor='rgba(0,0,0,0.08)', zeroline=False, range=[0, y_top]),
        xaxis=dict(gridcolor='rgba(0,0,0,0.08)'),
    )
    st.plotly_chart(line_fig, use_container_width=True, key='tc_line_jaccard', config=_plotly_cfg())
    _download_plotly_as_png(line_fig, f"jaccard_similarity_{tissue1_name}_vs_{tissue2_name}.png")
    
    # === DOWNLOAD SECTION ===
    display_download_section("📥 Download Analysis Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("**📊 Data Tables**")
        create_csv_download(df_age_analysis, f"age_analysis_{tissue1_name}_vs_{tissue2_name}.csv", 
                           "⬇️ Age Analysis CSV")
        
        # Combined dataset
        df1_labeled = df1.copy()
        df1_labeled['Tissue'] = tissue1_name
        df2_labeled = df2.copy()
        df2_labeled['Tissue'] = tissue2_name
        combined_df = pd.concat([df1_labeled, df2_labeled], ignore_index=True)
        create_csv_download(combined_df, f"combined_{tissue1_name}_{tissue2_name}.csv", 
                           "⬇️ Combined Dataset CSV")
    
    with col2:
        st.markdown("**🤝 Shared Genes**")
        shared_df = pd.DataFrame(shared_genes_all, columns=['Gene'])
        create_csv_download(shared_df, f"shared_genes_{tissue1_name}_{tissue2_name}.csv", 
                           "⬇️ Shared Genes CSV")
        
        # Show preview
        if len(shared_genes_all) > 0:
            st.text(f"Preview: {', '.join(shared_genes_all[:5])}{'...' if len(shared_genes_all) > 5 else ''}")
    
    with col3:
        st.markdown("**🧬 Exclusive Genes**")
        exclusive_1_df = pd.DataFrame(exclusive_1_all, columns=['Gene'])
        create_csv_download(exclusive_1_df, f"exclusive_{tissue1_name}.csv", 
                           f"⬇️ {tissue1_name} Exclusive CSV")
        
        exclusive_2_df = pd.DataFrame(exclusive_2_all, columns=['Gene'])
        create_csv_download(exclusive_2_df, f"exclusive_{tissue2_name}.csv", 
                           f"⬇️ {tissue2_name} Exclusive CSV")
    
    # === DETAILED GENE LISTS ===
    st.markdown("""
    <div class="analysis-section">
        <h2>📋 Detailed Gene Lists</h2>
    </div>
    """, unsafe_allow_html=True)
    
    tab1, tab2, tab3 = st.tabs(["🤝 Shared Genes", f"🧬 {tissue1_name} Exclusive", f"🧬 {tissue2_name} Exclusive"])
    
    with tab1:
        st.markdown(f"**{len(shared_genes_all)} genes shared between both tissues**")
        if shared_genes_all:
            st.text(", ".join(shared_genes_all))
    
    with tab2:
        st.markdown(f"**{len(exclusive_1_all)} genes exclusive to {tissue1_name}**")
        if exclusive_1_all:
            st.text(", ".join(exclusive_1_all))
    
    with tab3:
        st.markdown(f"**{len(exclusive_2_all)} genes exclusive to {tissue2_name}**")
        if exclusive_2_all:
            st.text(", ".join(exclusive_2_all))