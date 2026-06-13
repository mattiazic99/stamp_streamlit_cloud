import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from utils.parsing import parse_multiple_stamp_files  # USO PARSING CENTRALIZZATO
from components.downloads import create_csv_download, display_download_section, create_multiple_csv_download


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


# ── helper: Plotly heatmap (tissue × tissue, full or masked) ──
def _interactive_heatmap_tissue(matrix, labels, title, colorbar_label,
                                 cmap="rdbu", mask_upper=False,
                                 fmt="float", key=None,
                                 group_separator=None):
    n = matrix.shape[0]
    z_vals = matrix.copy().astype(float)
    text_vals, hover_vals = [], []
    for i in range(n):
        row_t, row_h = [], []
        for j in range(n):
            if mask_upper and j > i:
                z_vals[i, j] = None
                row_t.append("")
                row_h.append("")
            else:
                val = matrix[i, j]
                if fmt == "int":
                    row_t.append(f"{int(val)}")
                    row_h.append(f"<b>{labels[i]}</b> × <b>{labels[j]}</b><br>{colorbar_label}: {int(val)}")
                else:
                    row_t.append(f"{val:.2f}")
                    row_h.append(f"<b>{labels[i]}</b> × <b>{labels[j]}</b><br>{colorbar_label}: {val:.2f}")
        text_vals.append(row_t)
        hover_vals.append(row_h)

    fig = go.Figure(data=go.Heatmap(
        z=z_vals, x=labels, y=labels,
        colorscale=cmap,
        text=text_vals,
        texttemplate="<b>%{text}</b>",
        textfont=dict(size=10),
        hovertext=hover_vals,
        hovertemplate="%{hovertext}<extra></extra>",
        colorbar=dict(title=colorbar_label, thickness=14, len=0.75),
        xgap=2, ygap=2,
    ))

    side = min(560, max(380, n * 65))
    fig.update_layout(
        title=dict(text=title, font=dict(size=13), x=0.5),
        height=side, width=side + 40,
        margin=dict(l=100, r=40, t=60, b=100),
        yaxis=dict(autorange="reversed", scaleanchor="x", scaleratio=1, constrain="domain"),
        xaxis=dict(constrain="domain", tickangle=45),
        plot_bgcolor="white",
    )

    if group_separator is not None:
        fig.add_hline(y=group_separator - 0.5, line_width=3, line_color="black")
        fig.add_vline(x=group_separator - 0.5, line_width=3, line_color="black")

    _left, _centre, _right = st.columns([1, 3, 1])
    with _centre:
        st.plotly_chart(fig, use_container_width=True, key=key, config=_plotly_cfg())
    _download_plotly_as_png(fig, f"{title.replace(' ', '_')[:60]}.png")


def show():
    """Group Comparison Analysis Page"""
    
    st.header("👥 Group Comparison Analysis")
    st.markdown("Compare gene expression patterns between custom tissue groups.")
    
    age_groups = ["30–39", "40–49", "50–59", "60–69", "70–79"]
    
    # ── Data source toggle ─────────────────────────────────────────
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
    from data_loader import get_available_tissues, get_sets_for_tissue, sets_to_gui_format, sets_to_gene_sets, complete_mode

    # Pre-loaded (GTEx) mode removed by design: users always upload their files.
    data_source = "📂 Upload files"

    complete = complete_mode()

    if data_source == "📦 Pre-loaded (GTEx)":
        version = st.session_state.get("gtex_version", "v10")
        all_tissues = get_available_tissues(version, complete)

        st.markdown(f"""
        <div class="analysis-section">
            <h3>🧪 Select Tissues — GTEx {version}</h3>
            <p>Select at least 4 tissues to create and compare custom groups</p>
        </div>
        """, unsafe_allow_html=True)

        selected_tissues = st.multiselect(
            "Select tissues:",
            all_tissues,
            default=all_tissues[:6],
            key="gc_tissues",
        )

        if not selected_tissues or len(selected_tissues) < 4:
            st.info("👆 Please select at least 4 tissues for group comparison analysis.")
            return

        data_parsed = {}
        for tissue in selected_tissues:
            sets_dict = get_sets_for_tissue(version, tissue, complete)
            _, counts, df = sets_to_gui_format(sets_dict)
            gene_sets_list = sets_to_gene_sets(sets_dict)
            data_parsed[tissue] = {
                'gene_sets': gene_sets_list,
                'dataframe': df,
                'counts': counts,
                'total_genes': sum(len(s) for s in gene_sets_list),
                'original_filename': f"{tissue}.txt",
                'clean_name': tissue,
            }
        summary = {
            'total_files': len(selected_tissues),
            'successful_parses': len(selected_tissues),
            'failed_parses': 0,
            'total_unique_genes': len(set().union(*(set().union(*d['gene_sets']) for d in data_parsed.values()))),
            'tissue_names': list(data_parsed.keys()),
            'rejected_files': [],
        }
        st.success(f"✅ {len(selected_tissues)} tissues loaded from GTEx {version}!")

    else:
        st.markdown("""
        <div class="analysis-section">
            <h3>📂 Upload Multiple Tissue Files</h3>
            <p>Upload tissue files to create and compare custom groups</p>
        </div>
        """, unsafe_allow_html=True)

        uploaded_files = st.file_uploader(
            "📂 Upload STAMP .txt files", 
            type=["txt"],
            accept_multiple_files=True, 
            key="group_comparison_files",
            help="Upload multiple tissue gene switching files for group comparison"
        )

        if not uploaded_files or len(uploaded_files) < 4:
            st.info("👆 Please upload at least 4 tissue files to enable group comparison analysis.")

            st.markdown("### 👥 Group Comparison Features")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("""
                **🎯 Group Analysis:**
                - Custom tissue grouping
                - Inter-group comparisons
                - Shared vs exclusive genes
                - Statistical significance testing
                """)
            with col2:
                st.markdown("""
                **📊 Comparison Metrics:**
                - Jaccard similarity
                - Gene overlap percentages
                - Group-specific patterns
                - Multi-group intersections
                """)

            st.markdown("### 💡 Example Group Comparisons")
            examples = [
                "**Organ Systems**: Heart, Liver, Kidney vs Brain, Muscle, Lung",
                "**Metabolic vs Structural**: Liver, Pancreas vs Bone, Cartilage", 
                "**Central vs Peripheral**: Brain, Spinal Cord vs Skin, Muscle",
                "**High vs Low Metabolism**: Heart, Brain, Liver vs Bone, Skin"
            ]
            for example in examples:
                st.markdown(f"- {example}")
            return

        st.success(f"✅ {len(uploaded_files)} tissue files loaded successfully!")

        result = parse_multiple_stamp_files(uploaded_files, age_groups)
        _rejected = result['summary'].get('rejected_files', [])
        if _rejected:
            for _rej in _rejected:
                st.error(
                    f"❌ **{_rej['filename']}** is not a valid STAMP file:\n\n"
                    + "\n".join(f"- {e}" for e in _rej.get('errors', []))
                )
            st.info(
                "ℹ️ **STAMP format**: each file must have exactly **5 lines** "
                "(one per age group), with **space-separated gene names** on each line."
            )

        data_parsed = result['data']
        summary = result['summary']

    # Estrai tessuti e converti formato per compatibilità
    tissues = list(data_parsed.keys())
    data = {}
    for tissue in tissues:
        data[tissue] = data_parsed[tissue]['gene_sets']
    
    # Group creation section
    st.markdown("""
    <div class="analysis-section">
        <h2>👥 Create Tissue Groups</h2>
    </div>
    """, unsafe_allow_html=True)
    
    st.info("💡 Select at least 2 tissues for each group to enable comparison.")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### 🔵 Group 1")
        group1 = st.multiselect(
            "Select tissues for Group 1:",
            tissues,
            key="group1_tissues",
            help="Choose tissues for the first group"
        )
        
        if group1:
            group1_name = st.text_input(
                "Group 1 Name:",
                value="Group 1",
                key="group1_name",
                help="Enter a descriptive name for Group 1"
            )
    
    with col2:
        st.markdown("#### 🔴 Group 2")
        available_for_group2 = [t for t in tissues if t not in group1]
        group2 = st.multiselect(
            "Select tissues for Group 2:",
            available_for_group2,
            key="group2_tissues",
            help="Choose tissues for the second group"
        )
        
        if group2:
            group2_name = st.text_input(
                "Group 2 Name:",
                value="Group 2",
                key="group2_name",
                help="Enter a descriptive name for Group 2"
            )
    
    # Validation
    if len(group1) < 2 or len(group2) < 2:
        st.warning("⚠️ Please select at least 2 tissues for each group to enable comparison.")
        return
    
    # Analysis scope selection
    st.markdown("### ⚙️ Analysis Options")
    col1, col2 = st.columns(2)
    
    with col1:
        analysis_scope = st.selectbox(
            "📊 Analysis Scope:",
            ["Full Lifespan", "Specific Age Group", "All Age Groups Separately"],
            help="Choose the temporal scope for group comparison"
        )
        
        if analysis_scope == "Specific Age Group":
            selected_age = st.selectbox("🎯 Select Age Group:", age_groups, index=2)
    
    with col2:
        comparison_metrics = st.multiselect(
            "📈 Comparison Metrics:",
            ["Jaccard Similarity", "Shared Gene Count", "Overlap Percentage", "Statistical Tests"],
            default=["Jaccard Similarity", "Shared Gene Count"],
            help="Choose which metrics to calculate"
        )
    
    # Group overview
    st.markdown("""
    <div class="analysis-section">
        <h2>📊 Group Overview</h2>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown(f"#### 🔵 {group1_name}")
        st.write(f"**Tissues:** {', '.join(group1)}")
        
        group1_genes = set()
        if analysis_scope == "Full Lifespan":
            group1_genes = set.union(*[set.union(*data[tissue]) for tissue in group1])
        elif analysis_scope == "Specific Age Group":
            selected_idx = age_groups.index(selected_age)
            group1_genes = set.union(*[data[tissue][selected_idx] for tissue in group1])
        else:
            group1_genes = set.union(*[set.union(*data[tissue]) for tissue in group1])
        
        st.write(f"**Total Unique Genes:** {len(group1_genes)}")
        st.write(f"**Average Genes per Tissue:** {len(group1_genes) / len(group1):.1f}")
    
    with col2:
        st.markdown(f"#### 🔴 {group2_name}")
        st.write(f"**Tissues:** {', '.join(group2)}")
        
        group2_genes = set()
        if analysis_scope == "Full Lifespan":
            group2_genes = set.union(*[set.union(*data[tissue]) for tissue in group2])
        elif analysis_scope == "Specific Age Group":
            group2_genes = set.union(*[data[tissue][selected_idx] for tissue in group2])
        else:
            group2_genes = set.union(*[set.union(*data[tissue]) for tissue in group2])
        
        st.write(f"**Total Unique Genes:** {len(group2_genes)}")
        st.write(f"**Average Genes per Tissue:** {len(group2_genes) / len(group2):.1f}")
    
    # === MAIN COMPARISON ANALYSIS ===
    if analysis_scope in ["Full Lifespan", "Specific Age Group"]:
        perform_group_comparison(
            group1, group2, group1_name, group2_name, 
            group1_genes, group2_genes, data, 
            comparison_metrics, analysis_scope,
            selected_age if analysis_scope == "Specific Age Group" else None,
            age_groups
        )
    
    elif analysis_scope == "All Age Groups Separately":
        st.markdown("""
        <div class="analysis-section">
            <h2>📅 Age-by-Age Group Comparison</h2>
        </div>
        """, unsafe_allow_html=True)
        
        # Perform comparison for each age group
        age_comparison_results = []
        
        for age_idx, age in enumerate(age_groups):
            st.markdown(f"### 📊 Age Group: {age}")
            
            group1_age_genes = set.union(*[data[tissue][age_idx] for tissue in group1])
            group2_age_genes = set.union(*[data[tissue][age_idx] for tissue in group2])
            
            shared_genes = group1_age_genes & group2_age_genes
            exclusive1 = group1_age_genes - group2_age_genes
            exclusive2 = group2_age_genes - group1_age_genes
            union_genes = group1_age_genes | group2_age_genes
            
            jaccard_sim = len(shared_genes) / len(union_genes) if len(union_genes) > 0 else 0
            
            age_comparison_results.append({
                'Age Group': age,
                f'{group1_name} Genes': len(group1_age_genes),
                f'{group2_name} Genes': len(group2_age_genes),
                'Shared Genes': len(shared_genes),
                'Jaccard Similarity': jaccard_sim,
                f'{group1_name} Exclusive': len(exclusive1),
                f'{group2_name} Exclusive': len(exclusive2)
            })
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.markdown(f"""
                <div class="metric-card">
                    <h3>{len(shared_genes)}</h3>
                    <p>Shared Genes</p>
                </div>
                """, unsafe_allow_html=True)
            with col2:
                st.markdown(f"""
                <div class="metric-card">
                    <h3>{jaccard_sim:.3f}</h3>
                    <p>Jaccard Similarity</p>
                </div>
                """, unsafe_allow_html=True)
            with col3:
                st.markdown(f"""
                <div class="metric-card">
                    <h3>{len(exclusive1)}</h3>
                    <p>{group1_name} Exclusive</p>
                </div>
                """, unsafe_allow_html=True)
            with col4:
                st.markdown(f"""
                <div class="metric-card">
                    <h3>{len(exclusive2)}</h3>
                    <p>{group2_name} Exclusive</p>
                </div>
                """, unsafe_allow_html=True)
        
        df_age_comparison = pd.DataFrame(age_comparison_results)
        
        st.markdown("### 📋 Age-by-Age Comparison Summary")
        st.dataframe(df_age_comparison, use_container_width=True)
        
        # Plot age-wise trends — split into 4 separate charts
        ages = list(df_age_comparison['Age Group'])

        # 1. Jaccard similarity trend (Plotly line chart)
        jaccard_values = list(df_age_comparison['Jaccard Similarity'])
        y_top_j = max(jaccard_values) * 1.15 if jaccard_values else 1

        line_fig = go.Figure()
        line_fig.add_trace(go.Scatter(
            x=ages, y=jaccard_values,
            mode='lines+markers+text',
            line=dict(color='#9b59b6', width=3),
            marker=dict(size=9, color='#9b59b6'),
            fill='tozeroy', fillcolor='rgba(155,89,182,0.18)',
            text=[f'{v:.3f}' for v in jaccard_values],
            textposition='top center',
            hovertemplate='<b>%{x}</b><br>Jaccard: <b>%{y:.3f}</b><extra></extra>',
        ))
        line_fig.update_layout(
            title=dict(text='Jaccard Similarity Across Age Groups', font=dict(size=14), x=0.5),
            xaxis_title='Age Group', yaxis_title='Jaccard Similarity',
            height=350, margin=dict(l=50, r=30, t=60, b=50),
            plot_bgcolor='white',
            yaxis=dict(gridcolor='rgba(0,0,0,0.08)', zeroline=False, range=[0, y_top_j]),
            xaxis=dict(gridcolor='rgba(0,0,0,0.08)'),
        )
        st.plotly_chart(line_fig, use_container_width=True, key='gc_line_jaccard', config=_plotly_cfg())
        _download_plotly_as_png(line_fig, 'jaccard_similarity_trend.png')

        # 2. Shared genes bar (Plotly)
        shared_values = list(df_age_comparison['Shared Genes'])
        shared_fig = go.Figure()
        shared_fig.add_trace(go.Bar(
            x=ages, y=shared_values,
            marker=dict(color="rgba(46,204,113,0.85)", line=dict(color="black", width=1)),
            text=[f"<b>{v}</b>" for v in shared_values], textposition="outside",
            hovertemplate="<b>%{x}</b><br>Shared Genes: <b>%{y}</b><extra></extra>",
        ))
        y_max_s = max(shared_values) if shared_values else 1
        shared_fig.update_layout(
            title=dict(text="Shared Genes Across Age Groups", font=dict(size=14), x=0.5),
            yaxis_title="Number of Shared Genes", height=350,
            margin=dict(l=50, r=30, t=60, b=50), plot_bgcolor="white",
            yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False, range=[0, y_max_s*1.2+3]),
            xaxis=dict(type="category", categoryorder="array", categoryarray=age_groups),
        )
        st.plotly_chart(shared_fig, use_container_width=True, key="gc_bar_shared_age", config=_plotly_cfg())
        _download_plotly_as_png(shared_fig, 'shared_genes_by_age.png')

        # 3. Group gene counts comparison (Plotly grouped bar)
        group1_values = list(df_age_comparison[f'{group1_name} Genes'])
        group2_values = list(df_age_comparison[f'{group2_name} Genes'])
        n = len(ages)
        bar_w = 0.22 if n >= 4 else 0.16

        grp_fig = go.Figure()
        grp_fig.add_trace(go.Bar(
            name=group1_name, x=ages, y=group1_values, width=[bar_w]*n,
            marker=dict(color="rgba(100,170,230,0.85)", line=dict(color="rgb(8,48,107)", width=1)),
            text=[f"<b>{v}</b>" for v in group1_values], textposition="outside",
            hovertemplate=f"<b>{group1_name}</b><br>Age: <b>%{{x}}</b><br>Genes: <b>%{{y}}</b><extra></extra>",
        ))
        grp_fig.add_trace(go.Bar(
            name=group2_name, x=ages, y=group2_values, width=[bar_w]*n,
            marker=dict(color="rgba(250,128,114,0.85)", line=dict(color="darkred", width=1)),
            text=[f"<b>{v}</b>" for v in group2_values], textposition="outside",
            hovertemplate=f"<b>{group2_name}</b><br>Age: <b>%{{x}}</b><br>Genes: <b>%{{y}}</b><extra></extra>",
        ))
        all_v = group1_values + group2_values
        y_max_g = max(all_v) if all_v else 1
        grp_fig.update_layout(
            barmode="group",
            title=dict(text="Group Gene Counts by Age", font=dict(size=14), x=0.5),
            yaxis_title="Number of Genes", height=380,
            margin=dict(l=50, r=30, t=60, b=50), plot_bgcolor="white",
            yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False, range=[0, y_max_g*1.2+3]),
            xaxis=dict(type="category", categoryorder="array", categoryarray=age_groups),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        )
        st.plotly_chart(grp_fig, use_container_width=True, key="gc_bar_grp_counts", config=_plotly_cfg())
        _download_plotly_as_png(grp_fig, f"group_gene_counts_{group1_name}_vs_{group2_name}.png")

        # 4. Exclusive genes comparison (Plotly grouped bar)
        excl1_values = list(df_age_comparison[f'{group1_name} Exclusive'])
        excl2_values = list(df_age_comparison[f'{group2_name} Exclusive'])

        excl_fig = go.Figure()
        excl_fig.add_trace(go.Bar(
            name=f'{group1_name} Exclusive', x=ages, y=excl1_values, width=[bar_w]*n,
            marker=dict(color="rgba(243,156,18,0.85)", line=dict(color="black", width=1)),
            text=[f"<b>{v}</b>" for v in excl1_values], textposition="outside",
            hovertemplate=f"<b>{group1_name} Exclusive</b><br>Age: <b>%{{x}}</b><br>Genes: <b>%{{y}}</b><extra></extra>",
        ))
        excl_fig.add_trace(go.Bar(
            name=f'{group2_name} Exclusive', x=ages, y=excl2_values, width=[bar_w]*n,
            marker=dict(color="rgba(230,126,34,0.85)", line=dict(color="black", width=1)),
            text=[f"<b>{v}</b>" for v in excl2_values], textposition="outside",
            hovertemplate=f"<b>{group2_name} Exclusive</b><br>Age: <b>%{{x}}</b><br>Genes: <b>%{{y}}</b><extra></extra>",
        ))
        all_e = excl1_values + excl2_values
        y_max_e = max(all_e) if all_e else 1
        excl_fig.update_layout(
            barmode="group",
            title=dict(text="Exclusive Genes by Age", font=dict(size=14), x=0.5),
            yaxis_title="Number of Exclusive Genes", height=380,
            margin=dict(l=50, r=30, t=60, b=50), plot_bgcolor="white",
            yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False, range=[0, y_max_e*1.2+3]),
            xaxis=dict(type="category", categoryorder="array", categoryarray=age_groups),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        )
        st.plotly_chart(excl_fig, use_container_width=True, key="gc_bar_exclusive", config=_plotly_cfg())
        _download_plotly_as_png(excl_fig, f"exclusive_genes_{group1_name}_vs_{group2_name}.png")

        _download_plotly_as_png(grp_fig, f"age_comparison_{group1_name}_vs_{group2_name}.png")
        create_csv_download(df_age_comparison, f"age_comparison_{group1_name}_vs_{group2_name}.csv", 
                           "⬇️ Download Age Comparison CSV")
    
    # === DOWNLOAD SECTION ===
    display_download_section("📥 Download Group Comparison Results")
    
    download_data = {}
    
    group_composition = pd.DataFrame({
        'Group': [group1_name] * len(group1) + [group2_name] * len(group2),
        'Tissue': group1 + group2
    })
    download_data['group_composition'] = group_composition
    
    summary_stats = pd.DataFrame({
        'Metric': ['Group 1 Name', 'Group 2 Name', 'Group 1 Tissues', 'Group 2 Tissues',
                  'Analysis Scope', 'Total Comparisons'],
        'Value': [group1_name, group2_name, ', '.join(group1), ', '.join(group2),
                 analysis_scope, '1' if analysis_scope != "All Age Groups Separately" else str(len(age_groups))]
    })
    download_data['summary_statistics'] = summary_stats
    
    if 'df_age_comparison' in locals():
        download_data['age_by_age_comparison'] = df_age_comparison
    
    create_multiple_csv_download(
        download_data,
        f"group_comparison_{group1_name}_vs_{group2_name}.zip",
        "⬇️ Download Complete Analysis Package (ZIP)"
    )


def perform_group_comparison(group1, group2, group1_name, group2_name, 
                           group1_genes, group2_genes, data, 
                           comparison_metrics, analysis_scope, selected_age=None,
                           age_groups=None):
    """Perform detailed comparison between two tissue groups"""
    
    st.markdown(f"""
    <div class="analysis-section">
        <h2>🔄 {group1_name} vs {group2_name} Comparison</h2>
    </div>
    """, unsafe_allow_html=True)
    
    # Calculate comparison metrics
    shared_genes = group1_genes & group2_genes
    exclusive1 = group1_genes - group2_genes
    exclusive2 = group2_genes - group1_genes
    union_genes = group1_genes | group2_genes
    
    jaccard_similarity = len(shared_genes) / len(union_genes) if len(union_genes) > 0 else 0
    overlap_percentage = len(shared_genes) / min(len(group1_genes), len(group2_genes)) * 100 if min(len(group1_genes), len(group2_genes)) > 0 else 0
    
    # Display main metrics
    col1, col2, col3, col4, col5 = st.columns(5)
    
    with col1:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(shared_genes)}</h3>
            <p>Shared Genes</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{jaccard_similarity:.3f}</h3>
            <p>Jaccard<br>Similarity</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{overlap_percentage:.1f}%</h3>
            <p>Overlap<br>Percentage</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(exclusive1)}</h3>
            <p>{group1_name}<br>Exclusive</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col5:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(exclusive2)}</h3>
            <p>{group2_name}<br>Exclusive</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Visualization section
    st.markdown("### 📊 Visual Comparison")
    
    # 1. Pie chart (Plotly)
    sizes = [len(shared_genes), len(exclusive1), len(exclusive2)]
    labels_pie = [f'Shared ({len(shared_genes)})',
              f'{group1_name} Exclusive ({len(exclusive1)})',
              f'{group2_name} Exclusive ({len(exclusive2)})']
    colors_pie = ['#2ecc71', '#3498db', '#e74c3c']
    
    if sum(sizes) > 0:
        pie_fig = go.Figure(data=go.Pie(
            labels=labels_pie, values=sizes,
            marker=dict(colors=colors_pie, line=dict(color='white', width=2)),
            textinfo='percent+label', textposition='inside',
            hovertemplate='<b>%{label}</b><br>Genes: %{value}<br>Percentage: %{percent}<extra></extra>',
        ))
        pie_fig.update_layout(
            title=dict(text=f"Gene Distribution: {group1_name} vs {group2_name}",
                       font=dict(size=13), x=0.5),
            height=350, margin=dict(l=30, r=30, t=60, b=30),
            showlegend=False,
        )
        _left, _c, _r = st.columns([1, 2, 1])
        with _c:
            st.plotly_chart(pie_fig, use_container_width=True, key='gc_pie_shared', config=_plotly_cfg())
        _download_plotly_as_png(pie_fig, f"pie_{group1_name}_vs_{group2_name}.png")

    # 2. Bar chart comparison (Plotly)
    categories = [f'{group1_name}\nTotal', f'{group2_name}\nTotal', 'Shared', 'Union']
    values = [len(group1_genes), len(group2_genes), len(shared_genes), len(union_genes)]
    bar_colors = ['#3498db', '#e74c3c', '#2ecc71', '#9b59b6']

    bar_fig = go.Figure()
    bar_fig.add_trace(go.Bar(
        x=categories, y=values,
        marker=dict(color=bar_colors, line=dict(color="black", width=1)),
        text=[f"<b>{v}</b>" for v in values], textposition="outside",
        hovertemplate="<b>%{x}</b><br>Genes: <b>%{y}</b><extra></extra>",
    ))
    y_max = max(values) if values else 1
    bar_fig.update_layout(
        title=dict(text="Gene Count Comparison", font=dict(size=14), x=0.5),
        yaxis_title="Number of Genes", height=380,
        margin=dict(l=50, r=30, t=60, b=50), plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False, range=[0, y_max*1.2+3]),
    )
    _, _cb, _ = st.columns([1, 2, 1])
    with _cb:
        st.plotly_chart(bar_fig, use_container_width=True, key="gc_bar_count_comp", config=_plotly_cfg())
        _download_plotly_as_png(bar_fig, f"gene_count_comparison_{group1_name}_vs_{group2_name}.png")

    # 3. Overlap analysis barh (Plotly)
    overlap_data_vals = [len(shared_genes), len(exclusive1), len(exclusive2)]
    overlap_labels = ['Shared', f'{group1_name} Only', f'{group2_name} Only']
    overlap_colors = ['#2ecc71', '#3498db', '#e74c3c']

    barh_fig = go.Figure()
    barh_fig.add_trace(go.Bar(
        y=overlap_labels, x=overlap_data_vals, orientation="h",
        marker=dict(color=overlap_colors, line=dict(color="black", width=1)),
        text=[f"<b>{v}</b>" for v in overlap_data_vals], textposition="outside",
        hovertemplate="<b>%{y}</b><br>Genes: <b>%{x}</b><extra></extra>",
    ))
    barh_fig.update_layout(
        title=dict(text="Gene Overlap Analysis", font=dict(size=14), x=0.5),
        xaxis_title="Number of Genes", height=300,
        margin=dict(l=120, r=60, t=60, b=50), plot_bgcolor="white",
        xaxis=dict(gridcolor="rgba(0,0,0,0.08)"),
    )
    _, _ch, _ = st.columns([1, 2, 1])
    with _ch:
        st.plotly_chart(barh_fig, use_container_width=True, key="gc_barh_overlap", config=_plotly_cfg())
        _download_plotly_as_png(barh_fig, f"overlap_analysis_{group1_name}_vs_{group2_name}.png")

    # 4. Individual tissue contributions Group 1 (Plotly grouped bar)
    st.markdown(f"#### {group1_name} — Individual Tissue Contributions")
    tissue_contributions1 = []
    for tissue in group1:
        if analysis_scope == "Full Lifespan":
            tissue_genes = set.union(*data[tissue])
        else:
            selected_idx = age_groups.index(selected_age)
            tissue_genes = data[tissue][selected_idx]
        contribution = len(tissue_genes & shared_genes)
        total_genes = len(tissue_genes)
        tissue_contributions1.append({
            'Tissue': tissue, 'Shared Contribution': contribution, 'Total Genes': total_genes,
        })

    df_c1 = pd.DataFrame(tissue_contributions1)
    n1 = len(group1)
    bar_w1 = 0.22 if n1 >= 4 else 0.16

    c1_fig = go.Figure()
    c1_fig.add_trace(go.Bar(
        name='Shared Genes', x=df_c1['Tissue'], y=df_c1['Shared Contribution'],
        width=[bar_w1]*n1,
        marker=dict(color="rgba(46,204,113,0.85)", line=dict(color="black", width=1)),
        text=[f"<b>{v}</b>" for v in df_c1['Shared Contribution']], textposition="outside",
        hovertemplate="<b>%{x}</b><br>Shared: <b>%{y}</b><extra></extra>",
    ))
    c1_fig.add_trace(go.Bar(
        name='Total Genes', x=df_c1['Tissue'], y=df_c1['Total Genes'],
        width=[bar_w1]*n1,
        marker=dict(color="rgba(100,170,230,0.85)", line=dict(color="rgb(8,48,107)", width=1)),
        text=[f"<b>{v}</b>" for v in df_c1['Total Genes']], textposition="outside",
        hovertemplate="<b>%{x}</b><br>Total: <b>%{y}</b><extra></extra>",
    ))
    y_max_c1 = max(list(df_c1['Total Genes']) + list(df_c1['Shared Contribution'])) if len(df_c1) > 0 else 1
    c1_fig.update_layout(
        barmode="group",
        title=dict(text=f"{group1_name} — Tissue Contributions", font=dict(size=13), x=0.5),
        yaxis_title="Number of Genes", height=380,
        margin=dict(l=50, r=30, t=60, b=80), plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False, range=[0, y_max_c1*1.2+3]),
        xaxis=dict(tickangle=45),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    )
    st.plotly_chart(c1_fig, use_container_width=True, key="gc_bar_contrib1", config=_plotly_cfg())
    _download_plotly_as_png(c1_fig, f"tissue_contributions_{group1_name}.png")

    # 5. Individual tissue contributions Group 2 (Plotly grouped bar)
    st.markdown(f"#### {group2_name} — Individual Tissue Contributions")
    tissue_contributions2 = []
    for tissue in group2:
        if analysis_scope == "Full Lifespan":
            tissue_genes = set.union(*data[tissue])
        else:
            selected_idx = age_groups.index(selected_age)
            tissue_genes = data[tissue][selected_idx]
        contribution = len(tissue_genes & shared_genes)
        total_genes = len(tissue_genes)
        tissue_contributions2.append({
            'Tissue': tissue, 'Shared Contribution': contribution, 'Total Genes': total_genes,
        })

    df_c2 = pd.DataFrame(tissue_contributions2)
    n2 = len(group2)
    bar_w2 = 0.22 if n2 >= 4 else 0.16

    c2_fig = go.Figure()
    c2_fig.add_trace(go.Bar(
        name='Shared Genes', x=df_c2['Tissue'], y=df_c2['Shared Contribution'],
        width=[bar_w2]*n2,
        marker=dict(color="rgba(46,204,113,0.85)", line=dict(color="black", width=1)),
        text=[f"<b>{v}</b>" for v in df_c2['Shared Contribution']], textposition="outside",
        hovertemplate="<b>%{x}</b><br>Shared: <b>%{y}</b><extra></extra>",
    ))
    c2_fig.add_trace(go.Bar(
        name='Total Genes', x=df_c2['Tissue'], y=df_c2['Total Genes'],
        width=[bar_w2]*n2,
        marker=dict(color="rgba(250,128,114,0.85)", line=dict(color="darkred", width=1)),
        text=[f"<b>{v}</b>" for v in df_c2['Total Genes']], textposition="outside",
        hovertemplate="<b>%{x}</b><br>Total: <b>%{y}</b><extra></extra>",
    ))
    y_max_c2 = max(list(df_c2['Total Genes']) + list(df_c2['Shared Contribution'])) if len(df_c2) > 0 else 1
    c2_fig.update_layout(
        barmode="group",
        title=dict(text=f"{group2_name} — Tissue Contributions", font=dict(size=13), x=0.5),
        yaxis_title="Number of Genes", height=380,
        margin=dict(l=50, r=30, t=60, b=80), plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False, range=[0, y_max_c2*1.2+3]),
        xaxis=dict(tickangle=45),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
    )
    st.plotly_chart(c2_fig, use_container_width=True, key="gc_bar_contrib2", config=_plotly_cfg())
    _download_plotly_as_png(c2_fig, f"tissue_contributions_{group2_name}.png")

    # 6. Similarity heatmap between groups (Plotly)
    st.markdown("#### 🔥 Inter-tissue Similarity (Within and Between Groups)")
    all_group_tissues = group1 + group2
    similarity_matrix = np.zeros((len(all_group_tissues), len(all_group_tissues)))
    
    for i, t1 in enumerate(all_group_tissues):
        for j, t2 in enumerate(all_group_tissues):
            if analysis_scope == "Full Lifespan":
                g1 = set.union(*data[t1])
                g2 = set.union(*data[t2])
            else:
                selected_idx = age_groups.index(selected_age)
                g1 = data[t1][selected_idx]
                g2 = data[t2][selected_idx]
            intersection = len(g1 & g2)
            union = len(g1 | g2)
            similarity_matrix[i, j] = intersection / union if union > 0 else 0

    _interactive_heatmap_tissue(
        similarity_matrix, all_group_tissues,
        "Inter-tissue Similarity (Within and Between Groups)",
        "Jaccard Similarity", cmap="rdbu", mask_upper=False,
        fmt="float", key="gc_hm_similarity",
        group_separator=len(group1)
    )
    
    # Detailed gene lists
    st.markdown("### 📋 Detailed Gene Lists")
    
    tab1, tab2, tab3 = st.tabs([
        f"🤝 Shared Genes ({len(shared_genes)})",
        f"🔵 {group1_name} Exclusive ({len(exclusive1)})",
        f"🔴 {group2_name} Exclusive ({len(exclusive2)})"
    ])
    
    with tab1:
        if shared_genes:
            st.markdown(f"**{len(shared_genes)} genes shared between both groups:**")
            shared_list = sorted(list(shared_genes))
            st.text(", ".join(shared_list))
            shared_df = pd.DataFrame(shared_list, columns=['Gene'])
            create_csv_download(shared_df, f"shared_genes_{group1_name}_vs_{group2_name}.csv", 
                               "⬇️ Download Shared Genes CSV")
        else:
            st.info("No shared genes found between the groups.")
    
    with tab2:
        if exclusive1:
            st.markdown(f"**{len(exclusive1)} genes exclusive to {group1_name}:**")
            exclusive1_list = sorted(list(exclusive1))
            st.text(", ".join(exclusive1_list))
            excl1_df = pd.DataFrame(exclusive1_list, columns=['Gene'])
            create_csv_download(excl1_df, f"exclusive_{group1_name}.csv", 
                               f"⬇️ Download {group1_name} Exclusive CSV")
        else:
            st.info(f"No genes exclusive to {group1_name}.")
    
    with tab3:
        if exclusive2:
            st.markdown(f"**{len(exclusive2)} genes exclusive to {group2_name}:**")
            exclusive2_list = sorted(list(exclusive2))
            st.text(", ".join(exclusive2_list))
            excl2_df = pd.DataFrame(exclusive2_list, columns=['Gene'])
            create_csv_download(excl2_df, f"exclusive_{group2_name}.csv", 
                               f"⬇️ Download {group2_name} Exclusive CSV")
        else:
            st.info(f"No genes exclusive to {group2_name}.")
    
    # Statistical analysis
    if "Statistical Tests" in comparison_metrics:
        st.markdown("""
        <div class="analysis-section">
            <h3>📊 Statistical Analysis</h3>
        </div>
        """, unsafe_allow_html=True)
        
        from scipy.stats import hypergeom
        
        all_genes_universe = set()
        for tissue in list(data.keys()):
            if analysis_scope == "Full Lifespan":
                all_genes_universe.update(set.union(*data[tissue]))
            else:
                selected_idx = age_groups.index(selected_age)
                all_genes_universe.update(data[tissue][selected_idx])
        
        total_genes = len(all_genes_universe)
        overlap_observed = len(shared_genes)
        group1_size = len(group1_genes)
        group2_size = len(group2_genes)
        
        p_value = 1 - hypergeom.cdf(overlap_observed - 1, total_genes, group1_size, group2_size)
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{p_value:.2e}</h3>
                <p>P-value<br>(Hypergeometric)</p>
            </div>
            """, unsafe_allow_html=True)
        with col2:
            significance = "Significant" if p_value < 0.05 else "Not Significant"
            st.markdown(f"""
            <div class="metric-card">
                <h3>{significance}</h3>
                <p>Statistical<br>Significance</p>
            </div>
            """, unsafe_allow_html=True)
        with col3:
            expected_overlap = (group1_size * group2_size) / total_genes
            fold_enrichment = overlap_observed / expected_overlap if expected_overlap > 0 else 0
            st.markdown(f"""
            <div class="metric-card">
                <h3>{fold_enrichment:.2f}</h3>
                <p>Fold<br>Enrichment</p>
            </div>
            """, unsafe_allow_html=True)
        
        if p_value < 0.05:
            st.success(f"✅ The overlap between {group1_name} and {group2_name} is statistically significant (p < 0.05)")
        else:
            st.info(f"ℹ️ The overlap between {group1_name} and {group2_name} is not statistically significant (p ≥ 0.05)")
    
    # Summary statistics table
    st.markdown("### 📈 Summary Statistics")
    
    summary_data = {
        'Metric': [
            f'{group1_name} Total Genes', f'{group2_name} Total Genes',
            'Shared Genes', 'Union Genes', 'Jaccard Similarity',
            'Overlap Percentage', f'{group1_name} Exclusive', f'{group2_name} Exclusive'
        ],
        'Value': [
            len(group1_genes), len(group2_genes), len(shared_genes), len(union_genes),
            f"{jaccard_similarity:.3f}", f"{overlap_percentage:.1f}%",
            len(exclusive1), len(exclusive2)
        ]
    }
    
    summary_df = pd.DataFrame(summary_data)
    st.dataframe(summary_df, use_container_width=True)
    
    create_csv_download(summary_df, f"summary_statistics_{group1_name}_vs_{group2_name}.csv", 
                       "⬇️ Download Summary Statistics CSV")
