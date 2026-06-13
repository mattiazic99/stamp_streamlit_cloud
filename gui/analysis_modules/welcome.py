import streamlit as st


def show():
    """Welcome page - explains how to use the STAMP application."""
    
    # Hero section
    st.markdown(
        '<div style="'
        'background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);'
        'padding: 40px;'
        'border-radius: 20px;'
        'text-align: center;'
        'color: white;'
        'margin-bottom: 30px;'
        'box-shadow: 0 12px 35px rgba(17, 153, 142, 0.4);'
        '">'
        '<h1 style="margin: 0; font-size: 2.5rem; text-shadow: 2px 2px 6px rgba(0,0,0,0.3);">'
        '&#x1F44B; Welcome to STAMP'
        '</h1>'
        '<div style="font-size: 1.2rem; margin-top: 10px; opacity: 0.95;">'
        'Gene Switching Explorer &#8212; Analyze gene expression patterns in chronic pathologies'
        '</div>'
        '</div>',
        unsafe_allow_html=True
    )

    # How it works section
    st.markdown(
        '<div style="'
        'background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);'
        'padding: 30px;'
        'border-radius: 20px;'
        'margin-bottom: 25px;'
        'box-shadow: 0 8px 25px rgba(0, 0, 0, 0.1);'
        'border: 2px solid rgba(102, 126, 234, 0.15);'
        '">'
        '<h2 style="color: #2c3e50; text-align: center; margin-bottom: 15px; font-size: 1.8rem;">'
        '&#x1F680; How does it work?'
        '</h2>'
        '<div style="color: #555; text-align: center; font-size: 1.05rem; margin-bottom: 10px; line-height: 1.6;">'
        'STAMP works in <strong>3 simple steps</strong>. First, generate your analysis files, '
        'then you can use them as many times as you want across the different analysis pages.'
        '</div>'
        '</div>',
        unsafe_allow_html=True
    )

    # 3 columns for the 3 steps
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown(
            '<div style="'
            'background: linear-gradient(135deg, #ff6b35 0%, #f7931e 100%);'
            'padding: 25px;'
            'border-radius: 18px;'
            'text-align: center;'
            'color: white;'
            'min-height: 280px;'
            'box-shadow: 0 10px 30px rgba(255, 107, 53, 0.35);'
            'border: 2px solid rgba(255, 255, 255, 0.15);'
            '">'
            '<div style="font-size: 3rem; margin-bottom: 10px;">1&#xFE0F;&#x20E3;</div>'
            '<h3 style="margin: 0 0 12px 0; font-size: 1.3rem; text-shadow: 1px 1px 3px rgba(0,0,0,0.2);">'
            '&#x1F6E0;&#xFE0F; Generate Files'
            '</h3>'
            '<div style="font-size: 0.95rem; line-height: 1.6; opacity: 0.95;">'
            'Use the <strong>STAMP Generator</strong> to create analysis files from raw expression data (TPM).'
            '<br><br>'
            '&#x26A0;&#xFE0F; <em>This step must be done <strong>at least once</strong> before you can analyze any data.</em>'
            '</div>'
            '</div>',
            unsafe_allow_html=True
        )

    with col2:
        st.markdown(
            '<div style="'
            'background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);'
            'padding: 25px;'
            'border-radius: 18px;'
            'text-align: center;'
            'color: white;'
            'min-height: 280px;'
            'box-shadow: 0 10px 30px rgba(102, 126, 234, 0.35);'
            'border: 2px solid rgba(255, 255, 255, 0.15);'
            '">'
            '<div style="font-size: 3rem; margin-bottom: 10px;">2&#xFE0F;&#x20E3;</div>'
            '<h3 style="margin: 0 0 12px 0; font-size: 1.3rem; text-shadow: 1px 1px 3px rgba(0,0,0,0.2);">'
            '&#x1F4CA; Analyze Results'
            '</h3>'
            '<div style="font-size: 0.95rem; line-height: 1.6; opacity: 0.95;">'
            'Load the generated files into the <strong>analysis pages</strong> to explore gene switching patterns.'
            '<br><br>'
            '&#x1F504; <em>You can reuse the generated files <strong>as many times as you want</strong>.</em>'
            '</div>'
            '</div>',
            unsafe_allow_html=True
        )

    with col3:
        st.markdown(
            '<div style="'
            'background: linear-gradient(135deg, #28a745 0%, #20c997 100%);'
            'padding: 25px;'
            'border-radius: 18px;'
            'text-align: center;'
            'color: white;'
            'min-height: 280px;'
            'box-shadow: 0 10px 30px rgba(40, 167, 69, 0.35);'
            'border: 2px solid rgba(255, 255, 255, 0.15);'
            '">'
            '<div style="font-size: 3rem; margin-bottom: 10px;">3&#xFE0F;&#x20E3;</div>'
            '<h3 style="margin: 0 0 12px 0; font-size: 1.3rem; text-shadow: 1px 1px 3px rgba(0,0,0,0.2);">'
            '&#x1F4E5; Explore &amp; Download'
            '</h3>'
            '<div style="font-size: 0.95rem; line-height: 1.6; opacity: 0.95;">'
            'View <strong>heatmaps, charts, and tables</strong>. Then download your results as images or CSV files.'
            '<br><br>'
            '&#x1F4C8; <em>Each page offers <strong>interactive visualizations</strong> and downloads.</em>'
            '</div>'
            '</div>',
            unsafe_allow_html=True
        )

    st.markdown("<br>", unsafe_allow_html=True)

    # STAMP file format explanation section
    example_lines = (
        "APOE TP53 BRCA1 EGFR<br>"
        "MYC PTEN RB1 VHL<br>"
        "APC KRAS PIK3CA IDH1<br>"
        "CDKN2A ATM SMAD4<br>"
        "MLH1 MSH2 MSH6 PMS2"
    )
    
    st.markdown(
        '<div style="'
        'background: linear-gradient(135deg, #fff3cd 0%, #ffeeba 100%);'
        'padding: 30px;'
        'border-radius: 20px;'
        'margin-bottom: 25px;'
        'box-shadow: 0 8px 25px rgba(0, 0, 0, 0.1);'
        'border: 2px solid rgba(255, 193, 7, 0.3);'
        '">'
        '<h2 style="color: #856404; text-align: center; margin-bottom: 20px; font-size: 1.8rem;">'
        '&#x1F4C4; What are STAMP files?'
        '</h2>'
        '<div style="color: #856404; font-size: 1.05rem; line-height: 1.7; margin-bottom: 20px; text-align: center;">'
        'The analysis pages <strong>only accept files in STAMP format</strong>. '
        'These are <code>.txt</code> files with a specific structure, generated automatically '
        'by the <strong>STAMP Generator</strong>.'
        '</div>'
        '<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">'
        '<div style="'
        'background: white;'
        'padding: 20px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(133, 100, 4, 0.15);'
        '">'
        '<div style="color: #856404; font-weight: bold; font-size: 1.1rem; margin-bottom: 12px;">'
        '&#x1F4CB; File Structure'
        '</div>'
        '<div style="color: #555; font-size: 0.95rem; line-height: 1.7;">'
        '&#x2022; Each file is a <code>.txt</code> file representing <strong>one tissue</strong><br>'
        '&#x2022; Contains <strong>5 lines</strong>, one for each age group (30-39, 40-49, 50-59, 60-69, 70-79)<br>'
        '&#x2022; Each line contains <strong>space-separated gene names</strong> that switch in that age group<br>'
        '&#x2022; These files are the output of the STAMP Generator'
        '</div>'
        '</div>'
        '<div style="'
        'background: white;'
        'padding: 20px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(133, 100, 4, 0.15);'
        '">'
        '<div style="color: #856404; font-weight: bold; font-size: 1.1rem; margin-bottom: 12px;">'
        '&#x1F4A1; Example (tissue_brain.txt)'
        '</div>'
        '<div style="background: #f8f9fa; padding: 12px; border-radius: 8px; font-size: 0.85rem;'
        ' color: #155724; border: 1px solid #c3e6cb; font-family: Monaco, Consolas, monospace; line-height: 1.8;">'
        + example_lines +
        '</div>'
        '<div style="color: #888; font-size: 0.8rem; margin-top: 8px; font-style: italic;">'
        'Line 1 = Age 30-39, Line 2 = Age 40-49, ... Line 5 = Age 70-79'
        '</div>'
        '</div>'
        '</div>'
        '</div>',
        unsafe_allow_html=True
    )

    # Available pages section
    st.markdown(
        '<div style="'
        'background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);'
        'padding: 30px;'
        'border-radius: 20px;'
        'margin-bottom: 25px;'
        'box-shadow: 0 8px 25px rgba(0, 0, 0, 0.1);'
        'border: 2px solid rgba(102, 126, 234, 0.15);'
        '">'
        '<h2 style="color: #2c3e50; text-align: center; margin-bottom: 25px; font-size: 1.8rem;">'
        '&#x1F4CB; Available Pages'
        '</h2>'
        '<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px;">'
        '<div style="'
        'background: rgba(102, 126, 234, 0.1);'
        'padding: 18px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(102, 126, 234, 0.2);'
        '">'
        '<strong style="color: #667eea; font-size: 1.1rem;">&#x1F6E0;&#xFE0F; STAMP Generator</strong><br>'
        '<span style="font-size: 0.9rem; color: #555;">Generate STAMP files from raw TPM expression data</span>'
        '</div>'
        '<div style="'
        'background: rgba(102, 126, 234, 0.1);'
        'padding: 18px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(102, 126, 234, 0.2);'
        '">'
        '<strong style="color: #667eea; font-size: 1.1rem;">&#x1F4E4; Upload &amp; Single Analysis</strong><br>'
        '<span style="font-size: 0.9rem; color: #555;">Upload a file and analyze a single tissue</span>'
        '</div>'
        '<div style="'
        'background: rgba(102, 126, 234, 0.1);'
        'padding: 18px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(102, 126, 234, 0.2);'
        '">'
        '<strong style="color: #667eea; font-size: 1.1rem;">&#x1F504; Tissue Comparison</strong><br>'
        '<span style="font-size: 0.9rem; color: #555;">Compare gene switching across two tissues</span>'
        '</div>'
        '<div style="'
        'background: rgba(102, 126, 234, 0.1);'
        'padding: 18px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(102, 126, 234, 0.2);'
        '">'
        '<strong style="color: #667eea; font-size: 1.1rem;">&#x1F9EC; Multi-Tissue Analysis</strong><br>'
        '<span style="font-size: 0.9rem; color: #555;">Simultaneous analysis across multiple tissues</span>'
        '</div>'
        '<div style="'
        'background: rgba(102, 126, 234, 0.1);'
        'padding: 18px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(102, 126, 234, 0.2);'
        '">'
        '<strong style="color: #667eea; font-size: 1.1rem;">&#x1F4C5; Age-Specific Analysis</strong><br>'
        '<span style="font-size: 0.9rem; color: #555;">Analysis for specific age groups</span>'
        '</div>'
        '<div style="'
        'background: rgba(102, 126, 234, 0.1);'
        'padding: 18px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(102, 126, 234, 0.2);'
        '">'
        '<strong style="color: #667eea; font-size: 1.1rem;">&#x1F91D; Gene Sharing Analysis</strong><br>'
        '<span style="font-size: 0.9rem; color: #555;">Discover shared genes across different tissues</span>'
        '</div>'
        '<div style="'
        'background: rgba(102, 126, 234, 0.1);'
        'padding: 18px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(102, 126, 234, 0.2);'
        '">'
        '<strong style="color: #667eea; font-size: 1.1rem;">&#x1F50D; Single Gene Analysis</strong><br>'
        '<span style="font-size: 0.9rem; color: #555;">Detailed analysis of a single gene</span>'
        '</div>'
        '<div style="'
        'background: rgba(102, 126, 234, 0.1);'
        'padding: 18px;'
        'border-radius: 14px;'
        'border: 1px solid rgba(102, 126, 234, 0.2);'
        '">'
        '<strong style="color: #667eea; font-size: 1.1rem;">&#x1F465; Group Comparison</strong><br>'
        '<span style="font-size: 0.9rem; color: #555;">Compare gene groups across tissues</span>'
        '</div>'
        '</div>'
        '</div>',
        unsafe_allow_html=True
    )

    # Important note - prominent and visible
    st.markdown(
        '<div style="'
        'background: linear-gradient(135deg, #e8f5e8 0%, #d4edda 100%);'
        'padding: 25px 30px;'
        'border-radius: 18px;'
        'margin: 25px 0;'
        'border: 2px solid rgba(40, 167, 69, 0.3);'
        'box-shadow: 0 6px 20px rgba(40, 167, 69, 0.15);'
        'text-align: center;'
        '">'
        '<div style="font-size: 1.5rem; margin-bottom: 8px;">&#x1F4A1;</div>'
        '<div style="color: #155724; font-size: 1.15rem; font-weight: bold; margin-bottom: 8px;">'
        'Already generated your STAMP files?'
        '</div>'
        '<div style="color: #155724; font-size: 1.05rem; line-height: 1.6;">'
        'If you have already run the STAMP Generator at least once, your files are ready! '
        'You can skip this step and go directly to any <strong>analysis page</strong> '
        'using the navigation menu on the left.'
        '</div>'
        '</div>',
        unsafe_allow_html=True
    )

    # "Get Started" button
    st.markdown("<br>", unsafe_allow_html=True)
    
    col_left, col_center, col_right = st.columns([1, 2, 1])
    with col_center:
        if st.button(
            "\U0001F680 Get Started \u2014 Go to STAMP Generator",
            use_container_width=True,
            type="primary"
        ):
            st.session_state.current_page = "stamp_generator"
            st.rerun()
