
'''import streamlit as st
import subprocess
import os
import time
import zipfile
from io import BytesIO
from components.downloads import create_download_button, display_download_section
from components.styling import apply_plot_style
def show():
    """STAMP Dataset Generator Page"""
    
    st.header("🛠️ STAMP Dataset Generator")
    st.markdown("Generate STAMP-compatible datasets from raw expression data for analysis.")
    
    # === Configuration ===
    DIR_NORMALIZED = "tpm_normalizzati"
    DIR_SETS = "sets_stamp"
    DIR_SETS_MAPPED = "sets_stamp_symboli"
    GENE_MAP_FILE = "all_genes.txt"
    SCRIPT_STEP2 = "step_2.py"
    SCRIPT_STEP3A = "step_3a_corretto.py"
    
    # === Session state ===
    if "processed" not in st.session_state:
        st.session_state.processed = False
    if "tissues" not in st.session_state:
        st.session_state.tissues = []
    
    # === Info Section ===
    st.markdown("""
    <div class="analysis-section">
        <h3>📋 About STAMP Generator</h3>
        <p>This tool processes raw gene expression data to generate STAMP-compatible files for age-based gene switching analysis.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Show process overview
    with st.expander("🔍 View Process Overview"):
        st.markdown("""
        ### 🔄 Generation Process:
        
        1. **📊 Data Normalization** - Normalize TPM values across age groups
        2. **🎯 Threshold Application** - Apply expression threshold to identify active genes
        3. **🧬 Gene Set Creation** - Generate age-specific gene sets per tissue
        4. **🏷️ Symbol Mapping** - Convert Ensembl IDs to gene symbols
        5. **📁 File Export** - Create downloadable STAMP-compatible files
        
        ### 📂 Output Files:
        - **`*_sets_stamp.txt`** - Raw gene sets with Ensembl IDs
        - **`*_sets_stamp_mapped.txt`** - Gene sets with human-readable symbols
        """)
    
    # === Prerequisites Check ===
    st.markdown("""
    <div class="analysis-section">
        <h3>🔍 Prerequisites Check</h3>
    </div>
    """, unsafe_allow_html=True)
    
    # Check directories and files
    checks = [
        (DIR_NORMALIZED, "Normalized data directory"),
        (GENE_MAP_FILE, "Gene mapping file"),
        (SCRIPT_STEP2, "Step 2 processing script"),
        (SCRIPT_STEP3A, "Step 3A processing script")
    ]
    
    all_checks_passed = True
    for path, description in checks:
        if os.path.exists(path):
            st.success(f"✅ {description}: `{path}` found")
        else:
            st.error(f"❌ {description}: `{path}` not found")
            all_checks_passed = False
    
    if not all_checks_passed:
        st.error("❌ Prerequisites not met. Please ensure all required files and directories are present.")
        st.info("💡 Make sure you have run the data preparation pipeline before using this generator.")
        return
    
    # === Tissue Selection ===
    st.markdown("""
    <div class="analysis-section">
        <h3>📂 Tissue Selection</h3>
    </div>
    """, unsafe_allow_html=True)
    
    # Get available tissues
    file_list = [f for f in os.listdir(DIR_NORMALIZED) if f.endswith("_normalized.csv")]
    all_tissues = sorted([f.replace("_normalized.csv", "") for f in file_list])
    
    if not all_tissues:
        st.error("❌ No normalized tissue files found. Please run the normalization step first.")
        return
    
    st.info(f"📊 Found {len(all_tissues)} normalized tissue datasets available for processing.")
    
    # Selection interface
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
    
    # === Processing Parameters ===
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
    
    # === Generation Process ===
    st.markdown("""
    <div class="analysis-section">
        <h3>🚀 File Generation</h3>
    </div>
    """, unsafe_allow_html=True)
    
    if st.button("🚀 Generate STAMP Files", type="primary", use_container_width=True):
        st.session_state.tissues = selected_tissues
        
        # Progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        try:
            with st.spinner("🔄 Generating STAMP files..."):
                # Step 1: Update threshold in script
                status_text.text("📝 Updating processing parameters...")
                progress_bar.progress(10)
                
                with open(SCRIPT_STEP3A, "r", encoding="utf-8") as f:
                    lines = f.readlines()
                with open(SCRIPT_STEP3A, "w", encoding="utf-8") as f:
                    for line in lines:
                        if line.strip().startswith("THRESHOLD ="):
                            f.write(f"THRESHOLD = {threshold}  # threshold to consider a gene 'expressed'\n")
                        else:
                            f.write(line)
                
                # Step 2: Run processing scripts
                status_text.text("⚙️ Running data processing pipeline...")
                progress_bar.progress(30)
                
                subprocess.run(["python", SCRIPT_STEP2], check=True)
                progress_bar.progress(50)
                
                subprocess.run(["python", SCRIPT_STEP3A], check=True)
                progress_bar.progress(70)
                
                # Step 3: Gene mapping
                status_text.text("🧬 Mapping gene symbols...")
                
                def load_gene_map():
                    gene_map = {}
                    with open(GENE_MAP_FILE, encoding="utf-8") as f:
                        for line in f:
                            if "(" in line and ")" in line:
                                ensembl = line.split("(")[0].strip()
                                symbol = line.split("(")[1].replace(")", "").strip()
                                gene_map[ensembl] = symbol
                    return gene_map
                
                os.makedirs(DIR_SETS_MAPPED, exist_ok=True)
                gene_map = load_gene_map()
                
                # Process each tissue
                for i, tissue in enumerate(selected_tissues):
                    status_text.text(f"🔬 Processing tissue: {tissue} ({i+1}/{len(selected_tissues)})")
                    
                    input_path = os.path.join(DIR_SETS, f"{tissue}_sets_stamp.txt")
                    output_path = os.path.join(DIR_SETS_MAPPED, f"{tissue}_sets_stamp_mapped.txt")
                    
                    if os.path.exists(input_path):
                        with open(input_path, encoding="utf-8") as fin, open(output_path, "w", encoding="utf-8") as fout:
                            for line in fin:
                                mapped_genes = [gene_map.get(g, g) for g in line.strip().split()]
                                fout.write(" ".join(mapped_genes) + "\n")
                    
                    progress_bar.progress(70 + (i + 1) * 25 // len(selected_tissues))
                
                progress_bar.progress(100)
                status_text.text("✅ Generation completed successfully!")
                
                st.session_state.processed = True
                st.success("🎉 All STAMP files have been generated successfully!")
                
                # Show summary
                st.markdown("### 📊 Generation Summary")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown(f"""
                    <div class="metric-card">
                        <h3>{len(selected_tissues)}</h3>
                        <p>Tissues Processed</p>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col2:
                    st.markdown(f"""
                    <div class="metric-card">
                        <h3>{threshold}</h3>
                        <p>Expression Threshold</p>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col3:
                    total_files = len([f for f in os.listdir(DIR_SETS_MAPPED) if f.endswith("_sets_stamp_mapped.txt")])
                    st.markdown(f"""
                    <div class="metric-card">
                        <h3>{total_files}</h3>
                        <p>Files Generated</p>
                    </div>
                    """, unsafe_allow_html=True)
                
        except subprocess.CalledProcessError as e:
            st.error(f"❌ Error during processing: {e}")
            st.error("Please check that all required scripts are present and executable.")
        except Exception as e:
            st.error(f"❌ Unexpected error: {e}")
    
    # === Download Section ===
    if st.session_state.processed:
        display_download_section("📥 Download Generated Files")
        
        # Collect available files
        download_options = []
        for tissue in st.session_state.tissues:
            path1 = os.path.join(DIR_SETS, f"{tissue}_sets_stamp.txt")
            path2 = os.path.join(DIR_SETS_MAPPED, f"{tissue}_sets_stamp_mapped.txt")
            
            if os.path.exists(path1):
                download_options.append((f"{tissue}_sets_stamp.txt", path1, "Raw (Ensembl IDs)"))
            if os.path.exists(path2):
                download_options.append((f"{tissue}_sets_stamp_mapped.txt", path2, "Mapped (Gene Symbols)"))
        
        if download_options:
            st.markdown("### 📁 Available Files")
            
            # Show file list with descriptions
            file_df_data = []
            for name, path, desc in download_options:
                file_size = os.path.getsize(path) / 1024  # KB
                file_df_data.append({
                    "File Name": name,
                    "Type": desc,
                    "Size (KB)": f"{file_size:.1f}",
                    "Tissue": name.split("_")[0]
                })
            
            import pandas as pd
            file_df = pd.DataFrame(file_df_data)
            st.dataframe(file_df, use_container_width=True)
            
            # File selection for download
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("#### 📋 Select Files to Download")
                selected_files = st.multiselect(
                    "Choose files:",
                    [x[0] for x in download_options],
                    help="Select one or more files to download"
                )
            
            with col2:
                st.markdown("#### ⚡ Quick Selection")
                if st.button("📊 Select All Mapped Files", use_container_width=True):
                    selected_files = [x[0] for x in download_options if "mapped" in x[0]]
                    st.experimental_rerun()
                
                if st.button("🧬 Select All Raw Files", use_container_width=True):
                    selected_files = [x[0] for x in download_options if "mapped" not in x[0]]
                    st.experimental_rerun()
                
                if st.button("📁 Select All Files", use_container_width=True):
                    selected_files = [x[0] for x in download_options]
                    st.experimental_rerun()
            
            # Download buttons
            if selected_files:
                st.markdown("### ⬇️ Download Options")
                
                col1, col2 = st.columns(2)
                
                with col1:
                    # Individual file downloads
                    st.markdown("**Individual Downloads:**")
                    for name, path, desc in download_options:
                        if name in selected_files:
                            with open(path, 'rb') as f:
                                st.download_button(
                                    label=f"📄 {name}",
                                    data=f.read(),
                                    file_name=name,
                                    mime='text/plain',
                                    key=f"download_{name}"
                                )
                
                with col2:
                    # ZIP download
                    st.markdown("**Bulk Download:**")
                    if st.button("📦 Create ZIP Archive", use_container_width=True):
                        zip_buffer = BytesIO()
                        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zipf:
                            for name, path, desc in download_options:
                                if name in selected_files:
                                    zipf.write(path, arcname=name)
                        
                        zip_buffer.seek(0)
                        st.download_button(
                            label="📥 Download ZIP Archive",
                            data=zip_buffer.getvalue(),
                            file_name="stamp_files_generated.zip",
                            mime='application/zip',
                            use_container_width=True
                        )
            else:
                st.info("👆 Select files above to enable download options.")
        
        # Reset option
        st.markdown("---")
        if st.button("🔄 Start New Generation", type="secondary", use_container_width=True):
            st.session_state.processed = False
            st.session_state.tissues = []
            st.experimental_rerun()
    
    # === Help Section ===
    with st.expander("❓ Need Help?"):
        st.markdown("""
        ### 🆘 Troubleshooting
        
        **Common Issues:**
        
        1. **Missing prerequisites**: Ensure all required scripts and data directories are present
        2. **Processing errors**: Check that Python scripts have proper permissions
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
'''
''' Versione che funziona
import streamlit as st
import subprocess
import os
import time
import zipfile
from io import BytesIO
from components.downloads import create_download_button, display_download_section
from components.styling import apply_plot_style
def run_script_live(command, env=None):
    """
    Esegue uno script esterno mostrando l'output in tempo reale in Streamlit.
    Ritorna il codice di uscita del processo.
    """
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env
    )
    log_placeholder = st.empty()
    full_log = ""
    for line in process.stdout:
        full_log += line
        log_placeholder.text(full_log)
    return process.wait()
def show():
    """STAMP Dataset Generator Page"""
    st.header("🛠️ STAMP Dataset Generator")
    st.markdown("Generate STAMP-compatible datasets from raw expression data for analysis.")
    # === PATH DINAMICO ALLA ROOT DEL PROGETTO ===
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    ROOT_DIR = os.path.dirname(BASE_DIR)
    # === CARTELLE PRINCIPALI ===
    DATA_DIR = os.path.join(ROOT_DIR, "data")
    SCRIPTS_DIR = os.path.join(ROOT_DIR, "scripts")
    OUTPUT_DIR = os.path.join(ROOT_DIR, "output")
    # === PERCORSI STAMP ===
    DIR_NORMALIZED = os.path.join(OUTPUT_DIR, "tpm_normalizzati")
    DIR_SETS = os.path.join(OUTPUT_DIR, "sets_stamp")
    DIR_SETS_MAPPED = os.path.join(OUTPUT_DIR, "sets_stamp_symboli")
    # === FILE NECESSARI ===
    GENE_MAP_FILE = os.path.join(DATA_DIR, "all_genes.txt")
    # === SCRIPT PYTHON ===
    SCRIPT_STEP2 = os.path.join(SCRIPTS_DIR, "step_2.py")
    SCRIPT_STEP3A = os.path.join(SCRIPTS_DIR, "step_3a_corretto.py")
    # === SESSION STATE ===
    if "processed" not in st.session_state:
        st.session_state.processed = False
    if "tissues" not in st.session_state:
        st.session_state.tissues = []
    if "selected_files" not in st.session_state:
        st.session_state.selected_files = []
    # === INFO SECTION ===
    st.markdown("""
    <div class="analysis-section">
        <h3>📋 About STAMP Generator</h3>
        <p>This tool processes raw gene expression data to generate STAMP-compatible files for age-based gene switching analysis.</p>
    </div>
    """, unsafe_allow_html=True)
    # === OVERVIEW ===
    with st.expander("🔍 View Process Overview"):
        st.markdown("""
        ### 📄 Generation Process:
        
        1. **📊 Data Normalization** (Step 2)
        2. **🎯 Threshold Application** (Step 3A)
        3. **🧬 Gene Set Creation**
        4. **🏷️ Symbol Mapping**
        """)
    # === CHECK PREREQUISITI ===
    st.markdown("""
    <div class="analysis-section">
        <h3>🔍 Prerequisites Check</h3>
    </div>
    """, unsafe_allow_html=True)
    checks = [
        (DIR_NORMALIZED, "Normalized data directory"),
        (GENE_MAP_FILE, "Gene mapping file"),
        (SCRIPT_STEP2, "Step 1 processing script"),
        (SCRIPT_STEP3A, "Step 2 processing script")
    ]
    all_checks_passed = True
    for path, desc in checks:
        if os.path.exists(path):
            st.success(f"✅ {desc} found")
        else:
            st.error(f"❌ {desc} not found")
            all_checks_passed = False
    if not all_checks_passed:
        st.error("❌ Prerequisites not met.")
        return
    # === TISSUE SELECTION ===
    st.markdown("""
    <div class="analysis-section">
        <h3>📂 Tissue Selection</h3>
    </div>
    """, unsafe_allow_html=True)
    file_list = [f for f in os.listdir(DIR_NORMALIZED) if f.endswith("_normalized.csv")]
    all_tissues = sorted([f.replace("_normalized.csv", "") for f in file_list])
    if not all_tissues:
        st.error("❌ No normalized files found.")
        return
    st.info(f"📊 Found {len(all_tissues)} normalized tissues.")
    col1, col2 = st.columns([2, 1])
    with col1:
        choice = st.radio(
            "📌 Select tissues:",
            ["Single Tissue", "Multiple Tissues", "All Tissues"]
        )
        if choice == "Single Tissue":
            selected_tissues = [st.selectbox("🔬 Choose:", all_tissues)]
        elif choice == "Multiple Tissues":
            selected_tissues = st.multiselect("🔬 Choose:", all_tissues)
        else:
            selected_tissues = all_tissues
            st.info(f"📄 All {len(all_tissues)} tissues selected.")
    with col2:
        st.markdown("### 📊 Summary")
        if selected_tissues:
            st.markdown(f"<h3>{len(selected_tissues)}</h3><p>Tissues Selected</p>", unsafe_allow_html=True)
    if not selected_tissues:
        st.warning("⚠️ Select at least one tissue.")
        return
    # === THRESHOLD ===
    st.markdown("""
    <div class="analysis-section">
        <h3>⚙️ Processing Parameters</h3>
    </div>
    """, unsafe_allow_html=True)
    threshold = st.slider("🎚️ Gene Expression Threshold", 0.0, 1.0, 0.5, 0.1)
    # === GENERATION ===
    st.markdown("""
    <div class="analysis-section">
        <h3>🚀 File Generation</h3>
    </div>
    """, unsafe_allow_html=True)
    if st.button("🚀 Generate STAMP Files", use_container_width=True):
        st.session_state.tissues = selected_tissues
        progress = st.progress(0)
        status = st.empty()
        try:
            with st.spinner("🔄 Running pipeline..."):
                # === STEP 2 ===
                status.text("⚙️ Running Step 2 (Normalization)")
                rc = run_script_live(["python", SCRIPT_STEP2])
                if rc != 0:
                    raise Exception("Step 2 failed.")
                progress.progress(40)
                # === STEP 3A (con threshold via variabile ambiente) ===
                status.text("⚙️ Running Step 3A (Threshold & Switching)")
                env = os.environ.copy()
                env["STAMP_THRESHOLD"] = str(threshold)
                rc = run_script_live(["python", SCRIPT_STEP3A], env=env)
                if rc != 0:
                    raise Exception("Step 3A failed.")
                progress.progress(70)
                # === MAPPING ===
                status.text("🧬 Mapping gene symbols...")
                os.makedirs(DIR_SETS_MAPPED, exist_ok=True)
                gene_map = {}
                with open(GENE_MAP_FILE, encoding="utf-8") as f:
                    for line in f:
                        if "(" in line:
                            ensembl = line.split("(")[0].strip()
                            symbol = line.split("(")[1].replace(")", "").strip()
                            gene_map[ensembl] = symbol
                for i, tissue in enumerate(selected_tissues):
                    input_path = os.path.join(DIR_SETS, f"{tissue}_sets_stamp.txt")
                    output_path = os.path.join(DIR_SETS_MAPPED, f"{tissue}_sets_stamp_mapped.txt")
                    if os.path.exists(input_path):
                        with open(input_path) as fin, open(output_path, "w") as fout:
                            for line in fin:
                                mapped = [gene_map.get(g, g) for g in line.strip().split()]
                                fout.write(" ".join(mapped) + "\n")
                    progress.progress(70 + (i+1) * 30 // len(selected_tissues))
                progress.progress(100)
                status.text("✅ Generation completed!")
                st.session_state.processed = True
                st.success("🎉 All STAMP files have been generated.")
        except Exception as e:
            st.error(f"❌ Error: {e}")
            return
    # === DOWNLOAD SECTION ===
    if st.session_state.processed:
        display_download_section("📥 Download Generated Files")
        download_options = []
        for tissue in st.session_state.tissues:
            raw = os.path.join(DIR_SETS, f"{tissue}_sets_stamp.txt")
            mapped = os.path.join(DIR_SETS_MAPPED, f"{tissue}_sets_stamp_mapped.txt")
            if os.path.exists(raw):
                download_options.append((f"{tissue}_sets_stamp.txt", raw, "Raw"))
            if os.path.exists(mapped):
                download_options.append((f"{tissue}_sets_stamp_mapped.txt", mapped, "Mapped"))
        if download_options:
            st.markdown("### 📁 Available Files")
            file_data = []
            for name, path, desc in download_options:
                file_data.append({
                    "File Name": name,
                    "Type": desc,
                    "Size (KB)": f"{os.path.getsize(path)/1024:.1f}",
                    "Tissue": name.split("_")[0]
                })
            import pandas as pd
            st.dataframe(pd.DataFrame(file_data, use_container_width=True), use_container_width=True)
            selected = st.multiselect("Select Files:", [x[0] for x in download_options])
            if selected:
                for name, path, desc in download_options:
                    if name in selected:
                        with open(path, "rb") as f:
                            st.download_button(
                                label=f"📄 {name}",
                                data=f.read(),
                                file_name=name,
                                mime="text/plain"
                            )
    # === HELP ===
    with st.expander("❓ Need Help?"):
        st.write("""
        ### Common Issues
        - Check that all normalized files exist
        - Ensure scripts have permission to run
        """)
'''
import streamlit as st
import subprocess
import os
import zipfile
from io import BytesIO
from components.downloads import create_download_button, display_download_section
from components.styling import apply_plot_style
def run_script_live(command, env=None):
    """
    Esegue uno script esterno mostrando l'output in tempo reale in Streamlit.
    Ritorna il codice di uscita del processo.
    """
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env
    )
    log_placeholder = st.empty()
    full_log = ""
    for line in process.stdout:
        full_log += line
        # Mostra l'output man mano (come se fosse il terminale)
        log_placeholder.text(full_log)
    return process.wait()
import streamlit as st
import subprocess
import sys
import os
import zipfile
from io import BytesIO
from components.downloads import create_download_button, display_download_section
from components.styling import apply_plot_style
def run_script_live(command, env=None):
    """Esegue uno script esterno catturando l'output silenziosamente."""
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env
    )
    # Cattura output senza mostrarlo
    for _ in process.stdout:
        pass
    return process.wait()
def show():
    """STAMP Dataset Generator Page"""
    version = st.session_state.get("gtex_version", "v10")
    complete = st.session_state.get("tissue_mode", "all") == "complete"
    mode_label = "complete age bins" if complete else "all tissues"
    st.header(f"🛠️ STAMP Dataset Generator (GTEx {version})")
    st.markdown(
        f"Generate STAMP-compatible datasets from normalized expression data "
        f"(**{version}**, {mode_label}) for analysis."
    )
    
    # === PATHS ===
    import sys
    from pathlib import Path
    stamp_root = Path(__file__).resolve().parent.parent.parent
    if str(stamp_root) not in sys.path:
        sys.path.insert(0, str(stamp_root))
    from stamp.config import paths_for
    from stamp.io import _safe_filename

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(os.path.dirname(BASE_DIR), "data")
    
    DIR_NORMALIZED = str(paths_for(version, complete)["normalized"])
    DIR_SETS = str(paths_for(version, complete)["sets"])
    
    # We will save mapped sets in a temporary directory for download
    OUTPUT_DIR = os.path.join(os.path.dirname(BASE_DIR), "output")
    mapped_suffix = f"{version}_complete" if complete else version
    DIR_SETS_MAPPED = os.path.join(OUTPUT_DIR, f"sets_stamp_symboli_{mapped_suffix}")
    GENE_MAP_FILE = os.path.join(DATA_DIR, "all_genes.txt")
    SCRIPT_SWITCHING = str(stamp_root / "scripts" / "02_switching.py")
    
    # === Session state ===
    if "processed" not in st.session_state:
        st.session_state.processed = False
    if "tissues" not in st.session_state:
        st.session_state.tissues = []
    if "selected_files" not in st.session_state:
        st.session_state.selected_files = []
        
    # === Info ===
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
    # === Prerequisites Check ===
    st.markdown("""
    <div class="analysis-section"><h3>🔍 Prerequisites Check</h3></div>
    """, unsafe_allow_html=True)
    checks = [
        (DIR_NORMALIZED, "Normalized data directory"),
        (GENE_MAP_FILE, "Gene mapping file"),
        (SCRIPT_SWITCHING, "Main pipeline switching script (02_switching.py)"),
    ]
    all_ok = True
    for path, desc in checks:
        # Green "found" confirmations removed by design; only surface failures.
        if not os.path.exists(path):
            st.error(f"❌ {desc}: not found")
            all_ok = False
    if not all_ok:
        st.error("❌ Prerequisites not met. Please ensure all required files and directories are present.")
        st.info("💡 Make sure you have run the data preparation pipeline before using this generator.")
        return
    # === Tissue Selection ===
    st.markdown("""
    <div class="analysis-section"><h3>📂 Tissue Selection</h3></div>
    """, unsafe_allow_html=True)
    os.makedirs(DIR_NORMALIZED, exist_ok=True)
    
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
    
    # === Processing Parameters ===
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
    
    # === Generation ===
    st.markdown("""
    <div class="analysis-section"><h3>🚀 File Generation</h3></div>
    """, unsafe_allow_html=True)
    if st.button("🚀 Generate STAMP Files", type="primary", use_container_width=True):
        st.session_state.tissues = selected_tissues
        progress_bar = st.progress(0)
        status_text = st.empty()
        try:
            with st.spinner("🔄 Generating STAMP files using main pipeline..."):
                progress_bar.progress(5)
                
                # === Processing ===
                status_text.text("⚙️ Processing data (Threshold & Switching)...")

                # Run switching IN-PROCESS (no external subprocess) so it works
                # reliably on Streamlit Cloud, where the `stamp` package is not
                # pip-installed and a fresh subprocess cannot import it.
                from stamp.io import load_normalized_tissue, save_sets_txt
                from stamp.switching import identify_switching_genes

                _n = max(len(selected_tissues), 1)
                _failures = []
                for _i, _tissue in enumerate(selected_tissues):
                    try:
                        _df = load_normalized_tissue(version, _tissue, complete=complete)
                        _sets = identify_switching_genes(_df, threshold=threshold)
                        save_sets_txt(version, _tissue, _sets, complete=complete)
                    except Exception as _e:
                        _failures.append(f"{_tissue}: {_e}")
                    progress_bar.progress(5 + (_i + 1) * 55 // _n)

                if _failures:
                    st.error("Some tissues failed during switching:")
                    st.code("\n".join(_failures))
                    raise RuntimeError(f"{len(_failures)} tissue(s) failed")
                progress_bar.progress(60)
                
                # === Mapping ===
                status_text.text("🧬 Mapping gene symbols...")
                gene_map = {}
                with open(GENE_MAP_FILE, encoding="utf-8") as f:
                    for line in f:
                        if "(" in line and ")" in line:
                            ensembl = line.split("(")[0].strip()
                            symbol = line.split("(")[1].replace(")", "").strip()
                            # Key by the version-less Ensembl ID (drop the ".N"
                            # suffix) so the mapping works across GTEx releases:
                            # v8 and v10 use different version numbers for the
                            # same gene, and all_genes.txt carries v10 versions.
                            gene_map[ensembl.split(".")[0]] = symbol
                os.makedirs(DIR_SETS_MAPPED, exist_ok=True)
                for i, tissue in enumerate(selected_tissues):
                    status_text.text(f"🔬 Mapping symbols: {tissue} ({i+1}/{len(selected_tissues)})")
                    in_path = os.path.join(DIR_SETS, f"{_safe_filename(tissue)}_sets.txt")
                    out_path = os.path.join(DIR_SETS_MAPPED, f"{tissue}_sets_stamp_mapped.txt")
                    if os.path.exists(in_path):
                        with open(in_path, encoding="utf-8") as fin, open(out_path, "w", encoding="utf-8") as fout:
                            for line in fin:
                                fout.write(" ".join([gene_map.get(g.split(".")[0], g) for g in line.strip().split()]) + "\n")
                    progress_bar.progress(60 + (i + 1) * 35 // len(selected_tissues))
                progress_bar.progress(100)
                status_text.text("✅ Generation completed successfully!")
                st.session_state.processed = True
                st.success("🎉 All STAMP files have been generated successfully!")
                # Summary
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
                    n_files = len([f for f in os.listdir(DIR_SETS_MAPPED) if f.endswith("_sets_stamp_mapped.txt")])
                    st.markdown(f"""
                    <div class="metric-card">
                        <h3>{n_files}</h3>
                        <p>Files Generated</p>
                    </div>
                    """, unsafe_allow_html=True)
        except subprocess.CalledProcessError as e:
            st.error(f"❌ Error during processing: {e}")
            st.error("Please check that all required scripts are present and executable.")
        except Exception as e:
            st.error(f"❌ Unexpected error: {e}")
    # === Download Section ===
    if st.session_state.processed:
        display_download_section("📥 Download Generated Files")
        dl = []
        for tissue in st.session_state.tissues:
            p1 = os.path.join(DIR_SETS, f"{_safe_filename(tissue)}_sets.txt")
            p2 = os.path.join(DIR_SETS_MAPPED, f"{tissue}_sets_stamp_mapped.txt")
            if os.path.exists(p1):
                dl.append((f"{tissue}_sets_stamp.txt", p1, "Raw (Ensembl IDs)"))
            if os.path.exists(p2):
                dl.append((f"{tissue}_sets_stamp_mapped.txt", p2, "Mapped (Gene Symbols)"))
        if dl:
            import pandas as pd
            st.markdown("### 📁 Available Files")
            
            # Build list of all file names for quick selection
            all_file_names = [x[0] for x in dl]
            
            # Interactive table with checkboxes
            file_df_data = []
            for name, path, desc in dl:
                file_df_data.append({
                    "Select": name in st.session_state.selected_files,
                    "File Name": name,
                    "Type": desc,
                    "Size (KB)": f"{os.path.getsize(path)/1024:.1f}",
                    "Tissue": name.split("_")[0]
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
            
            # Update selected files from table
            st.session_state.selected_files = edited_df[edited_df["Select"]]["File Name"].tolist()
            
            # Quick selection buttons
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
                if st.button("📦 Create ZIP Archive", use_container_width=True):
                    buf = BytesIO()
                    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
                        for name, path, desc in dl:
                            if name in st.session_state.selected_files:
                                z.write(path, arcname=name)
                    buf.seek(0)
                    st.download_button(
                        label="📥 Download ZIP Archive",
                        data=buf.getvalue(),
                        file_name="stamp_files_generated.zip",
                        mime="application/zip",
                        use_container_width=True,
                        key="download_zip"
                    )
            else:
                st.info("👆 Select files above to enable download options.")
        st.markdown("---")
        if st.button("🔄 Start New Generation", type="secondary", use_container_width=True):
            st.session_state.processed = False
            st.session_state.tissues = []
            st.session_state.selected_files = []
            st.rerun()
    # === Help ===
    with st.expander("❓ Need Help?"):
        st.markdown("""
        ### 🆘 Troubleshooting
        **Common Issues:**
        1. **Missing prerequisites**: Ensure all required scripts and data directories are present
        2. **Processing errors**: Check that Python scripts have proper permissions
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
