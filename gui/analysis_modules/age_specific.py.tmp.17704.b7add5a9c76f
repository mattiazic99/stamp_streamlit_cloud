import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
from scipy.spatial.distance import squareform
from utils.analysis import compute_jaccard_matrix, compute_linkage, compute_common_genes_matrix, compute_percent_overlap_matrix
from utils.parsing import parse_multiple_stamp_files, extract_tissue_name  # AGGIUNTO IMPORT
from components.downloads import create_csv_download, display_download_section


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


# ── helper: Plotly heatmap (tissue × tissue, optional triangular mask) ──
def _interactive_heatmap_tissue(matrix, labels, title, colorbar_label,
                                 cmap="YlGnBu", mask_upper=True,
                                 fmt="float", key=None):
    """Reusable interactive heatmap for tissue-vs-tissue matrices."""
    n = matrix.shape[0]

    # Build display values & hover
    z_vals = matrix.copy().astype(float)
    text_vals = []
    hover_vals = []
    for i in range(n):
        row_t = []
        row_h = []
        for j in range(n):
            if mask_upper and j > i:
                z_vals[i, j] = None
                row_t.append("")
                row_h.append("")
            else:
                val = matrix[i, j]
                if fmt == "int":
                    row_t.append(f"{int(val)}")
                    row_h.append(f"<b>{labels[i]}</b> ∩ <b>{labels[j]}</b><br>{colorbar_label}: {int(val)}")
                else:
                    row_t.append(f"{val:.3f}")
                    row_h.append(f"<b>{labels[i]}</b> ∩ <b>{labels[j]}</b><br>{colorbar_label}: {val:.3f}")
        text_vals.append(row_t)
        hover_vals.append(row_h)

    fig = go.Figure(data=go.Heatmap(
        z=z_vals, x=labels, y=labels,
        colorscale=cmap,
        text=text_vals,
        texttemplate="<b>%{text}</b>",
        textfont=dict(size=11),
        hovertext=hover_vals,
        hovertemplate="%{hovertext}<extra></extra>",
        colorbar=dict(title=colorbar_label, thickness=14, len=0.75),
        xgap=2, ygap=2,
    ))

    side = min(520, max(380, n * 70))
    fig.update_layout(
        title=dict(text=title, font=dict(size=14), x=0.5),
        height=side, width=side + 40,
        margin=dict(l=80, r=40, t=60, b=80),
        yaxis=dict(autorange="reversed", scaleanchor="x", scaleratio=1, constrain="domain"),
        xaxis=dict(constrain="domain", tickangle=45),
        plot_bgcolor="white",
    )

    _left, _centre, _right = st.columns([1, 3, 1])
    with _centre:
        st.plotly_chart(fig, use_container_width=True, key=key, config=_plotly_cfg())
    _download_plotly_as_png(fig, f"{title.replace(' ', '_').replace('—', '_')}.png")


def show():
    """Age-Specific Analysis Page"""
    
    st.header("📅 Age-Specific Analysis")
    st.markdown("Analyze gene expression patterns for specific age groups across multiple tissues.")
    
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
            <p>Select multiple tissues to analyze age-specific patterns</p>
        </div>
        """, unsafe_allow_html=True)

        selected_tissues = st.multiselect(
            "Select tissues:",
            all_tissues,
            default=all_tissues[:5],
            key="as_tissues",
        )

        if not selected_tissues or len(selected_tissues) < 2:
            st.info("👆 Please select at least 2 tissues for age-specific analysis.")
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
        # Original upload mode
        st.markdown("""
        <div class="analysis-section">
            <h3>📂 Upload Multiple Tissue Files</h3>
            <p>Upload multiple tissue files to analyze age-specific patterns</p>
        </div>
        """, unsafe_allow_html=True)

        uploaded_files = st.file_uploader(
            "📂 Upload STAMP .txt files", 
            type=["txt"],
            accept_multiple_files=True, 
            key="age_specific_files",
            help="Upload multiple tissue gene switching files"
        )

        if not uploaded_files or len(uploaded_files) < 2:
            st.info("👆 Please upload at least 2 tissue files for age-specific analysis.")

            # Show age group information
            st.markdown("### 📋 Age Group Information")
            age_info_df = pd.DataFrame({
                'Age Group': age_groups,
                'Age Range': ['30-39 years', '40-49 years', '50-59 years', '60-69 years', '70-79 years'],
                'Life Stage': ['Early Adult', 'Middle Adult', 'Late Middle Age', 'Early Senior', 'Senior'],
                'Description': [
                    'Peak physical performance period',
                    'Career establishment phase',
                    'Pre-retirement transition',
                    'Early retirement phase',
                    'Advanced aging period'
                ]
            })
            st.dataframe(age_info_df, use_container_width=True)
            create_csv_download(age_info_df, "age_groups_info.csv", "⬇️ Download Age Groups Info (CSV)")
            return

        st.success(f"✅ {len(uploaded_files)} tissue files loaded successfully!")

        # Parse all files USANDO LA FUNZIONE CHE GIÀ PULISCE I NOMI
        result = parse_multiple_stamp_files(uploaded_files, age_groups)
        # ── Show rejected files (STAMP validation errors) ──────────────
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

    # Converti in formato compatibile con le funzioni di analisi
    data = {tissue: data_parsed[tissue]['gene_sets'] for tissue in data_parsed.keys()}
    tissues = list(data.keys())  # Nomi già puliti
    
    # Age group selection
    st.markdown("### 🎯 Select Age Group for Analysis")
    selected_age = st.selectbox(
        "Choose age group:",
        age_groups,
        index=2,  # Default to middle age (50-59)
        help="Select the specific age group to analyze across all tissues"
    )
    
    selected_idx = age_groups.index(selected_age)
    
    # Analysis options
    col1, col2 = st.columns(2)
    with col1:
        show_dendrograms = st.checkbox("🌳 Show dendrogram", value=True)
    with col2:
        clustering_method = st.selectbox("🔗 Clustering Method", 
                                       ["average", "complete", "single", "ward"])
        analysis_type = st.selectbox("📊 Analysis Focus", 
                                   ["Similarity", "Gene Counts", "Both"])
    
    # === AGE-SPECIFIC OVERVIEW ===
    st.markdown(f"""
    <div class="analysis-section">
        <h2>📊 Age Group {selected_age} - Tissue Overview</h2>
    </div>
    """, unsafe_allow_html=True)
    
    # Calculate statistics for selected age group
    age_stats = []
    for tissue in tissues:
        genes_in_age = data[tissue][selected_idx]
        total_genes_tissue = len(set.union(*data[tissue]))
        percentage_in_age = (len(genes_in_age) / total_genes_tissue * 100) if total_genes_tissue > 0 else 0
        
        age_stats.append({
            'Tissue': tissue,  # Nome già pulito
            'Genes in Age Group': len(genes_in_age),
            'Total Tissue Genes': total_genes_tissue,
            'Percentage in Age': f"{percentage_in_age:.1f}%"
        })
    
    df_age_stats = pd.DataFrame(age_stats)
    
    # Display metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        total_age_genes = sum([len(data[t][selected_idx]) for t in tissues])
        st.markdown(f"""
        <div class="metric-card">
            <h3>{total_age_genes}</h3>
            <p>Total Genes<br>in Age Group</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        unique_age_genes = len(set.union(*[data[t][selected_idx] for t in tissues]))
        st.markdown(f"""
        <div class="metric-card">
            <h3>{unique_age_genes}</h3>
            <p>Unique Genes<br>in Age Group</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        avg_genes = np.mean([len(data[t][selected_idx]) for t in tissues])
        st.markdown(f"""
        <div class="metric-card">
            <h3>{avg_genes:.1f}</h3>
            <p>Average Genes<br>per Tissue</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        max_genes = max([len(data[t][selected_idx]) for t in tissues])
        st.markdown(f"""
        <div class="metric-card">
            <h3>{max_genes}</h3>
            <p>Max Genes<br>in One Tissue</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Display detailed statistics table
    st.dataframe(df_age_stats, use_container_width=True)
    
    # Bar chart of gene counts per tissue for selected age (Plotly)
    if analysis_type in ["Gene Counts", "Both"]:
        st.markdown(f"### 📊 Gene Counts by Tissue (Age Group {selected_age})")
        
        gene_counts = [len(data[tissue][selected_idx]) for tissue in tissues]
        n_t = len(tissues)
        bar_w = 0.28 if n_t >= 4 else (0.22 if n_t >= 2 else 0.15)
        total = sum(gene_counts) or 1
        pcts = [(c / total * 100) for c in gene_counts]

        import plotly.express as px
        colors_list = px.colors.qualitative.Set3[:n_t]

        bar_fig = go.Figure()
        bar_fig.add_trace(go.Bar(
            x=tissues, y=gene_counts,
            width=[bar_w] * n_t,
            marker=dict(color=colors_list, line=dict(color="black", width=1.2)),
            text=[f"<b>{c}</b>" for c in gene_counts],
            textposition="outside",
            hovertemplate=(
                "<b>%{x}</b><br>"
                "Genes: <b>%{y}</b><br>"
                "Percentage: <b>%{customdata:.1f}%</b>"
                "<extra></extra>"
            ),
            customdata=pcts,
        ))

        y_max = max(gene_counts) if gene_counts else 1
        bar_fig.update_layout(
            title=dict(text=f"Gene Expression by Tissue — Age Group {selected_age}",
                       font=dict(size=15), x=0.5),
            xaxis_title="Tissue", yaxis_title="Number of Switching Genes",
            height=400,
            margin=dict(l=50, r=30, t=70, b=80),
            plot_bgcolor="white",
            yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False,
                       range=[0, y_max * 1.20 + 3]),
            xaxis=dict(tickangle=45),
        )
        st.plotly_chart(bar_fig, use_container_width=True, key="as_bar_gene_counts", config=_plotly_cfg())
        _download_plotly_as_png(bar_fig, f"gene_counts_age_{selected_age.replace('–','_')}.png")
    
    # === SIMILARITY ANALYSIS FOR SELECTED AGE ===
    if analysis_type in ["Similarity", "Both"]:
        st.markdown(f"""
        <div class="analysis-section">
            <h2>🔍 Similarity Analysis - Age Group {selected_age}</h2>
        </div>
        """, unsafe_allow_html=True)
        
        # Compute similarity matrix for selected age
        age_data = {tissue: [data[tissue][selected_idx]] for tissue in tissues}
        matrix_age_specific, tissues_sorted = compute_jaccard_matrix(age_data, mode="life")
        
        # Similarity heatmap (Plotly)
        st.markdown(f"### 📊 Tissue Similarity Matrix (Age {selected_age})")
        _interactive_heatmap_tissue(
            matrix_age_specific, tissues_sorted,
            f"Tissue Similarity — Age Group {selected_age}",
            "Jaccard Similarity", cmap="YlGnBu", mask_upper=True,
            fmt="float", key="as_hm_similarity"
        )
        
        # Hierarchical clustering for selected age (Plotly dendrogram)
        if show_dendrograms:
            st.markdown(f"### 🌳 Hierarchical Clustering (Age {selected_age})")
            
            linkage_age_specific = compute_linkage(matrix_age_specific, method=clustering_method)
            
            if linkage_age_specific is not None:
                # Convert similarity matrix to distance matrix for dendrogram
                dist_matrix = 1 - matrix_age_specific
                np.fill_diagonal(dist_matrix, 0)
                dist_matrix = np.clip(dist_matrix, 0, None)  # ensure non-negative

                dendro_fig = ff.create_dendrogram(
                    dist_matrix, labels=tissues_sorted,
                    linkagefun=lambda x: linkage_age_specific,
                )
                dendro_fig.update_layout(
                    title=dict(text=f"Tissue Clustering — Age {selected_age} ({clustering_method.title()} Linkage)",
                               font=dict(size=14), x=0.5),
                    xaxis_title='Tissue', yaxis_title='Distance (1 − Similarity)',
                    height=380, margin=dict(l=50, r=30, t=60, b=80),
                    plot_bgcolor='white',
                    yaxis=dict(gridcolor='rgba(0,0,0,0.08)'),
                )
                st.plotly_chart(dendro_fig, use_container_width=True, key='as_dendro', config=_plotly_cfg())
                _download_plotly_as_png(dendro_fig, f"dendrogram_age_{selected_age.replace('–','_')}.png")
        
        # Tissue ranking for selected age (Plotly horizontal bar)
        st.markdown(f"### 🏆 Tissue Similarity Ranking (Age {selected_age})")
        
        if len(tissues) >= 3:
            ref_tissue = st.selectbox("🔍 Choose reference tissue:", 
                                     tissues_sorted, key="age_ranking_ref")
            
            ref_idx = tissues_sorted.index(ref_tissue)
            similarities = []
            
            for i, tissue in enumerate(tissues_sorted):
                if tissue != ref_tissue:
                    similarity = matrix_age_specific[ref_idx, i] * 100
                    similarities.append((tissue, similarity))
            
            similarities.sort(key=lambda x: x[1], reverse=True)
            
            tissues_rank = [s[0] for s in similarities]
            sims_rank = [s[1] for s in similarities]

            import plotly.express as px
            colors_rank = px.colors.sequential.Viridis[:len(similarities)]

            rank_fig = go.Figure()
            rank_fig.add_trace(go.Bar(
                y=tissues_rank, x=sims_rank,
                orientation="h",
                marker=dict(color=sims_rank, colorscale="Viridis",
                            line=dict(color="black", width=1)),
                text=[f"<b>{v:.1f}%</b>" for v in sims_rank],
                textposition="outside",
                hovertemplate=(
                    "<b>%{y}</b><br>"
                    "Similarity: <b>%{x:.1f}%</b>"
                    "<extra></extra>"
                ),
            ))
            rank_fig.update_layout(
                title=dict(text=f"Similarity Ranking — Age {selected_age} (Ref: {ref_tissue})",
                           font=dict(size=14), x=0.5),
                xaxis_title="Similarity (%)",
                height=max(300, len(similarities) * 40 + 100),
                margin=dict(l=120, r=60, t=60, b=50),
                plot_bgcolor="white",
                xaxis=dict(gridcolor="rgba(0,0,0,0.08)"),
                yaxis=dict(autorange="reversed"),
            )
            st.plotly_chart(rank_fig, use_container_width=True, key="as_bar_ranking", config=_plotly_cfg())
            _download_plotly_as_png(rank_fig, f"ranking_age_{selected_age.replace('–','_')}_{ref_tissue}.png")
    
    # === AGE COMPARISON ACROSS TISSUES ===
    st.markdown("""
    <div class="analysis-section">
        <h2>📈 Age Group Comparison</h2>
    </div>
    """, unsafe_allow_html=True)
    
    # Compare selected age with other age groups
    st.markdown("### 📊 Gene Count Across All Age Groups")
    
    # Create matrix of gene counts per age per tissue
    age_comparison_data = []
    for tissue in tissues:
        for i, age in enumerate(age_groups):
            gene_count = len(data[tissue][i])
            is_selected = (age == selected_age)
            age_comparison_data.append({
                'Tissue': tissue,  # Nome già pulito
                'Age Group': age,
                'Gene Count': gene_count,
                'Selected': is_selected
            })
    
    df_age_comparison = pd.DataFrame(age_comparison_data)
    
    # Heatmap of gene counts across all ages and tissues (Plotly)
    pivot_data = df_age_comparison.pivot(index='Tissue', columns='Age Group', values='Gene Count')
    pivot_values = pivot_data.values.astype(np.int64)
    tissue_labels = list(pivot_data.index)
    age_labels = list(pivot_data.columns)

    # Build hover
    hover = []
    text_disp = []
    for i, t in enumerate(tissue_labels):
        row_h = []
        row_t = []
        for j, a in enumerate(age_labels):
            val = int(pivot_values[i, j])
            row_t.append(str(val))
            row_h.append(f"<b>{t}</b> × <b>{a}</b><br>Gene Count: {val}")
        hover.append(row_h)
        text_disp.append(row_t)

    hm_fig = go.Figure(data=go.Heatmap(
        z=pivot_values, x=age_labels, y=tissue_labels,
        colorscale="YlOrRd",
        text=text_disp,
        texttemplate="<b>%{text}</b>",
        textfont=dict(size=12),
        hovertext=hover,
        hovertemplate="%{hovertext}<extra></extra>",
        colorbar=dict(title="Gene<br>Count", thickness=14, len=0.75),
        xgap=2, ygap=2,
    ))

    # Highlight selected age column with a shape
    selected_col_idx = age_labels.index(selected_age) if selected_age in age_labels else None

    hm_fig.update_layout(
        title=dict(text=f"Gene Counts Across Age Groups (Highlighted: {selected_age})",
                   font=dict(size=14), x=0.5),
        height=max(400, len(tissue_labels) * 50 + 120),
        margin=dict(l=100, r=40, t=60, b=60),
        yaxis=dict(autorange="reversed"),
        plot_bgcolor="white",
    )

    if selected_col_idx is not None:
        hm_fig.add_shape(
            type="rect",
            x0=selected_col_idx - 0.5, x1=selected_col_idx + 0.5,
            y0=-0.5, y1=len(tissue_labels) - 0.5,
            line=dict(color="red", width=3),
        )

    _left, _centre, _right = st.columns([1, 3, 1])
    with _centre:
        st.plotly_chart(hm_fig, use_container_width=True, key="as_hm_age_comparison", config=_plotly_cfg())
    _download_plotly_as_png(hm_fig, f"age_comparison_heatmap_{selected_age.replace('–','_')}.png")
    
    # === DOWNLOAD SECTION ===
    display_download_section("📥 Download Age-Specific Analysis Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("**📊 Data Tables**")
        create_csv_download(df_age_stats, f"age_stats_{selected_age.replace('–','_')}.csv", 
                           "⬇️ Age Statistics CSV")
        create_csv_download(pivot_data, f"age_comparison_matrix_{selected_age.replace('–','_')}.csv", 
                           "⬇️ Age Comparison Matrix CSV")
    
    with col2:
        st.markdown("**🔍 Similarity Data**")
        if 'matrix_age_specific' in locals():
            df_similarity = pd.DataFrame(matrix_age_specific, 
                                       index=tissues_sorted, columns=tissues_sorted)
            create_csv_download(df_similarity, f"similarity_matrix_{selected_age.replace('–','_')}.csv", 
                               "⬇️ Similarity Matrix CSV")
    
    with col3:
        st.markdown("**🧬 Gene Lists**")
        # Create comprehensive gene list for selected age
        all_age_genes = []
        for tissue in tissues:
            tissue_genes = data[tissue][selected_idx]
            for gene in tissue_genes:
                all_age_genes.append({'Gene': gene, 'Tissue': tissue, 'Age Group': selected_age})
        
        if all_age_genes:
            df_all_genes = pd.DataFrame(all_age_genes)
            create_csv_download(df_all_genes, f"all_genes_{selected_age.replace('–','_')}.csv", 
                               "⬇️ All Genes CSV")