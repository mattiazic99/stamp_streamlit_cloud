"""STAMP Generator page.

Applies the switching rule to one, several or all age-complete tissues at a
user-selected threshold, resolves Ensembl identifiers to gene symbols, and
prepares the resulting five-line STAMP files as a ZIP download in session memory.
The files produced here are the input accepted by every analysis page.
"""
import os
import sys
import zipfile
from io import BytesIO
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from components.downloads import display_download_section  # noqa: E402
from data_loader import COMPLETE_AGE_BINS, GTEX_VERSION  # noqa: E402

STAMP_ROOT = Path(__file__).resolve().parent.parent.parent
if str(STAMP_ROOT) not in sys.path:
    sys.path.insert(0, str(STAMP_ROOT))




def _clear_generated_files():
    """Clear this session before rendering, so its old ZIP is not registered again."""
    st.session_state.processed = False
    for key in ("generator_exports", "generator_zip", "generator_threshold",
                "generator_run_dir", "file_editor"):
        st.session_state.pop(key, None)
    st.session_state.tissues = []
    st.session_state.selected_files = []


def show():
    """STAMP Dataset Generator Page"""
    version = GTEX_VERSION
    complete = COMPLETE_AGE_BINS
    mode_label = "complete age bins" if complete else "all tissues"
    st.header(f"🛠️ STAMP Dataset Generator (GTEx {version})")
    st.markdown(
        f"Generate STAMP-compatible datasets from normalized expression data "
        f"(**{version}**, {mode_label}) for analysis."
    )

    # Paths
    from stamp.config import paths_for, SWITCHING_BRACKETS
    from stamp.io import _safe_filename

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(os.path.dirname(BASE_DIR), "data")
    
    DIR_NORMALIZED = str(paths_for(version, complete)["normalized"])
    
    # Downloads live only in session memory; never write them to server storage.
    GENE_MAP_FILE = os.path.join(DATA_DIR, "all_genes.txt")
    
    # Session state
    if "processed" not in st.session_state:
        st.session_state.processed = False
    if st.session_state.processed and "generator_exports" not in st.session_state:
        # Discard old disk-based state after a development hot reload.
        st.session_state.processed = False
    if "tissues" not in st.session_state:
        st.session_state.tissues = []
    if "selected_files" not in st.session_state:
        st.session_state.selected_files = []
        
    # Info
    st.markdown("""
    <div class="analysis-section">
        <h3>📋 About STAMP Generator</h3>
        <p>This tool processes gene expression data to generate STAMP-compatible files
        for age-based gene switching analysis.</p>
    </div>
    """, unsafe_allow_html=True)
    with st.expander("🔍 View Process Overview"):
        st.markdown("""
        ### 📄 Generation Process:
        1. **🎯 Threshold Application** — Binarize gene expression (uses main pipeline)
        2. **🧬 Gene Set Creation** — Identify switching genes
        3. **🏷️ Symbol Mapping** — Convert Ensembl IDs to gene symbols
        4. **📁 File Export** — Generate downloadable files
        ### 📂 Output Files:
        - `*_sets.txt` — Raw gene sets (Ensembl IDs)
        - `*_sets_stamp_mapped.txt` — Gene sets (Gene Symbols)
        """)
    st.info("Generated files are temporary. Download the ZIP to keep a copy on your "
            "computer. After downloading, use Start New Generation to clear the current "
            "files. Generating again also replaces them.")
    # Prerequisites
    # Show missing prerequisites; stop before processing if any are absent.
    checks = [
        (DIR_NORMALIZED, "Normalized data directory"),
        (GENE_MAP_FILE, "Gene mapping file"),
    ]
    missing = [desc for path, desc in checks if not os.path.exists(path)]
    if missing:
        for desc in missing:
            st.error(f"❌ {desc}: not found")
        st.error("❌ Prerequisites not met. Please ensure all required files and directories are present.")
        st.info("💡 Make sure you have run the data preparation pipeline before using this generator.")
        return
    # Tissue Selection
    st.markdown("""
    <div class="analysis-section"><h3>📂 Tissue Selection</h3></div>
    """, unsafe_allow_html=True)
    
    all_tissues = []
    for f in os.listdir(DIR_NORMALIZED):
        if f.endswith(".parquet"):
            all_tissues.append(f.replace(".parquet", ""))
        elif f.endswith("_normalized.csv"):
            all_tissues.append(f.replace("_normalized.csv", ""))
            
    all_tissues = sorted(list(set(all_tissues)))
    
    if not all_tissues:
        st.error(f"❌ No normalized tissue files found in {DIR_NORMALIZED}.")
        return
    st.info(f"📊 Found {len(all_tissues)} normalized tissues available for analysis.")
    col1, col2 = st.columns([2, 1])
    with col1:
        choice = st.radio(
            "📌 Select tissues to process:",
            ["Single Tissue", "Multiple Tissues", "All Tissues"],
            help="Choose how many tissues to process in this run"
        )
        if choice == "Single Tissue":
            selected_tissues = [st.selectbox("🔬 Choose one tissue:", all_tissues)]
        elif choice == "Multiple Tissues":
            selected_tissues = st.multiselect(
                "🔬 Choose multiple tissues:",
                all_tissues,
                help="Select specific tissues to process"
            )
        else:
            selected_tissues = all_tissues
            st.info(f"🔄 All {len(all_tissues)} tissues will be processed.")
    with col2:
        st.markdown("### 📊 Selection Summary")
        if selected_tissues:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(selected_tissues)}</h3>
                <p>Tissues Selected</p>
            </div>
            """, unsafe_allow_html=True)
            
            if len(selected_tissues) <= 5:
                st.markdown("**Selected tissues:**")
                for tissue in selected_tissues:
                    st.markdown(f"• {tissue}")
            else:
                st.markdown(f"**Selected:** {selected_tissues[0]}, {selected_tissues[1]}, ... and {len(selected_tissues)-2} more")
    
    if not selected_tissues:
        st.warning("⚠️ No tissues selected. Please choose at least one tissue to process.")
        return
    
    # Processing Parameters
    st.markdown("""
    <div class="analysis-section">
        <h3>⚙️ Processing Parameters</h3>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        threshold = st.slider(
            "🎚️ Gene Expression Threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.5,
            step=0.1,
            help="Threshold to consider a gene as 'expressed' (0.0 = all genes, 1.0 = only highly expressed)"
        )
        
        st.markdown(f"""
        **Current threshold:** `{threshold}`
        
        - **Lower values** (0.1-0.3): Include more genes, detect subtle changes
        - **Medium values** (0.4-0.6): Balanced approach (recommended)
        - **Higher values** (0.7-1.0): Only highly expressed genes
        """)
    
    with col2:
        st.markdown("### 🎯 Expected Output")
        st.markdown(f"""
        **For each tissue, you'll get:**
        - 📁 Raw gene sets (Ensembl IDs)
        - 🏷️ Mapped gene sets (Gene symbols)
        - 📊 5 age groups per tissue (30-79 years)
        
        **Total files:** {len(selected_tissues) * 2} files
        """)
    
    # Generation
    st.markdown("""
    <div class="analysis-section"><h3>🚀 File Generation</h3></div>
    """, unsafe_allow_html=True)
    if st.button("🚀 Generate STAMP Files", type="primary", use_container_width=True):
        st.session_state.processed = False
        st.session_state.generator_exports = {}
        st.session_state.pop("generator_zip", None)
        st.session_state.pop("generator_run_dir", None)
        st.session_state.generator_threshold = threshold
        st.session_state.selected_files = []
        st.session_state.pop("file_editor", None)
        st.session_state.tissues = selected_tissues
        progress_bar = st.progress(0)
        status_text = st.empty()
        try:
            with st.spinner("🔄 Generating STAMP files using main pipeline..."):
                progress_bar.progress(5)
                
                # Processing
                status_text.text("⚙️ Processing data (Threshold & Switching)...")

                # Run in this process: Cloud does not install stamp for child processes.
                from stamp.io import load_normalized_tissue
                from stamp.switching import identify_switching_genes

                _n = max(len(selected_tissues), 1)
                _failures = []
                tissue_sets = {}
                for _i, _tissue in enumerate(selected_tissues):
                    try:
                        _df = load_normalized_tissue(version, _tissue, complete=complete)
                        _sets = identify_switching_genes(_df, threshold=threshold)
                        tissue_sets[_tissue] = _sets
                    except Exception as _e:
                        _failures.append(f"{_tissue}: {_e}")
                    progress_bar.progress(5 + (_i + 1) * 55 // _n)

                if _failures:
                    st.error("Some tissues failed during switching:")
                    st.code("\n".join(_failures))
                    raise RuntimeError(f"{len(_failures)} tissue(s) failed")
                progress_bar.progress(60)
                
                # Mapping
                status_text.text("🧬 Mapping gene symbols...")
                gene_map = {}
                with open(GENE_MAP_FILE, encoding="utf-8") as f:
                    for line in f:
                        if "(" in line and ")" in line:
                            ensembl = line.split("(")[0].strip()
                            symbol = line.split("(")[1].replace(")", "").strip()
                            # Match Ensembl IDs without the version suffix; the same gene
                            # can have different suffixes in v8 and v10.
                            gene_map[ensembl.split(".")[0]] = symbol
                exports = {}
                for i, tissue in enumerate(selected_tissues):
                    status_text.text(f"🔬 Mapping symbols: {tissue} ({i+1}/{len(selected_tissues)})")
                    lines = [" ".join(tissue_sets[tissue].get(bracket, []))
                             for bracket in SWITCHING_BRACKETS]
                    mapped_lines = [" ".join(gene_map.get(g.split(".")[0], g)
                                             for g in line.split()) for line in lines]
                    prefix = f"{_safe_filename(tissue)}_{version}_tau_{threshold:g}"
                    for suffix, text_lines, description in (
                        ("sets_stamp.txt", lines, "Raw (Ensembl IDs)"),
                        ("sets_stamp_mapped.txt", mapped_lines, "Mapped (Gene Symbols)"),
                    ):
                        exports[f"{prefix}_{suffix}"] = {
                            "data": ("\n".join(text_lines) + "\n").encode("utf-8"),
                            "type": description,
                            "tissue": tissue,
                        }
                    progress_bar.progress(60 + (i + 1) * 35 // len(selected_tissues))
                st.session_state.generator_exports = exports
                progress_bar.progress(100)
                status_text.text("✅ Generation completed successfully!")
                st.session_state.processed = True
                st.success("🎉 All STAMP files have been generated successfully!")
                st.markdown("### 📊 Generation Summary")
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"""
                    <div class="metric-card">
                        <h3>{len(selected_tissues)}</h3>
                        <p>Tissues Processed</p>
                    </div>
                    """, unsafe_allow_html=True)
                with col2:
                    n_files = len(exports)
                    st.markdown(f"""
                    <div class="metric-card">
                        <h3>{n_files}</h3>
                        <p>Files Generated</p>
                    </div>
                    """, unsafe_allow_html=True)
        except Exception as e:
            st.error(f"❌ Unexpected error: {e}")
    # Download Section
    if st.session_state.processed:
        display_download_section("📥 Download Generated Files")
        exports = st.session_state.generator_exports
        if exports:
            import pandas as pd
            st.markdown("### 📁 Available Files")
            
            all_file_names = list(exports)
            
            file_df_data = []
            for name, export in exports.items():
                file_df_data.append({
                    "Select": name in st.session_state.selected_files,
                    "File Name": name,
                    "Type": export["type"],
                    "Size (KB)": "{:.1f}".format(len(export["data"]) / 1024),
                    "Tissue": export["tissue"]
                })
            
            df_files = pd.DataFrame(file_df_data)
            
            edited_df = st.data_editor(
                df_files,
                column_config={
                    "Select": st.column_config.CheckboxColumn("Select", default=False),
                },
                disabled=["File Name", "Type", "Size (KB)", "Tissue"],
                hide_index=True,
                use_container_width=True,
                key="file_editor"
            )
            
            st.session_state.selected_files = edited_df[edited_df["Select"]]["File Name"].tolist()
            
            prepared = st.session_state.get("generator_zip")
            if prepared and prepared[0] != tuple(st.session_state.selected_files):
                st.session_state.pop("generator_zip", None)

            st.markdown("#### ⚡ Quick Selection")
            col_a, col_b, col_c, col_d = st.columns(4)
            with col_a:
                if st.button("✅ Select All", use_container_width=True):
                    st.session_state.selected_files = all_file_names
                    st.session_state.pop("file_editor", None)
                    st.rerun()
            with col_b:
                if st.button("📊 Only Mapped", use_container_width=True):
                    st.session_state.selected_files = [x for x in all_file_names if "mapped" in x]
                    st.session_state.pop("file_editor", None)
                    st.rerun()
            with col_c:
                if st.button("🧬 Only Raw", use_container_width=True):
                    st.session_state.selected_files = [x for x in all_file_names if "mapped" not in x]
                    st.session_state.pop("file_editor", None)
                    st.rerun()
            with col_d:
                if st.button("❌ Clear All", use_container_width=True):
                    st.session_state.selected_files = []
                    st.session_state.pop("file_editor", None)
                    st.rerun()
            
            if st.session_state.selected_files:
                st.markdown("### ⬇️ Download Options")
                st.info(f"📦 **{len(st.session_state.selected_files)}** files selected")
                selection = tuple(st.session_state.selected_files)
                if st.button("📦 Create ZIP Archive", use_container_width=True):
                    buf = BytesIO()
                    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as archive:
                        for name in selection:
                            archive.writestr(name, exports[name]["data"])
                    st.session_state.generator_zip = (selection, buf.getvalue())
                prepared = st.session_state.get("generator_zip")
                if prepared and prepared[0] == selection:
                    st.download_button(
                        label="📥 Download ZIP Archive",
                        data=prepared[1],
                        file_name="stamp_files_generated.zip",
                        mime="application/zip",
                        use_container_width=True,
                        on_click="ignore",
                        key="download_zip"
                    )
            else:
                st.info("👆 Select files above to enable download options.")
        st.markdown("---")
        st.button("🔄 Start New Generation", type="secondary", use_container_width=True,
                  on_click=_clear_generated_files)
    # Help
    with st.expander("❓ Need Help?"):
        st.markdown("""
        ### 🆘 Troubleshooting
        **Common Issues:**
        1. **Missing prerequisites**: Ensure the normalized atlas and gene mapping file are present
        2. **Processing errors**: Check that the normalized data can be read
        3. **No output files**: Verify that input data is properly formatted
        ### 📚 File Formats
        **Generated STAMP files contain:**
        - 5 lines (one per age group: 30-39, 40-49, 50-59, 60-69, 70-79)
        - Space-separated gene names/IDs per line
        - Compatible with all STAMP analysis modules
        ### 🔧 Parameters
        **Expression Threshold:**
        - Determines which genes are considered "active" in each age group
        - Higher values = more stringent filtering
        - Recommended: 0.4-0.6 for balanced analysis
        """)
