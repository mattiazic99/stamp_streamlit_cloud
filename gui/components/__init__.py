"""
STAMP Gene Switching Explorer - Components Module
===============================================

This module contains reusable components for styling, downloads, and
other utility functions used across all pages of the STAMP application.

Components:
- styling: Custom CSS and visual theme management
- downloads: Comprehensive download handling for all data types

Author: Gene Analysis Team
"""

__version__ = "2.0.0"
__author__ = "Gene Analysis Team"

# Import main components
try:
    from . import styling
    from . import downloads
except ImportError as e:
    import warnings
    warnings.warn(f"Some components could not be imported: {e}")

# Export commonly used functions for easy access
try:
    from .styling import (
        apply_custom_css,
        create_sidebar_style,
        create_status_indicator,
        create_metric_card_html
    )
    
    from .downloads import (
        initialize_download_components,
        create_matplotlib_download,
        create_plotly_download,
        create_csv_download,
        create_gene_list_download,
        create_download_section,
        end_download_section,
        get_export_settings
    )
except ImportError:
    # Graceful degradation if components are not available
    pass

def get_component_status():
    """
    Check the status of all components.
    
    Returns:
        dict: Status of each component
    """
    status = {}
    
    try:
        import components.styling
        status['styling'] = "Available"
    except ImportError as e:
        status['styling'] = f"Error: {str(e)}"
    
    try:
        import components.downloads
        status['downloads'] = "Available"
    except ImportError as e:
        status['downloads'] = f"Error: {str(e)}"
    
    return status
