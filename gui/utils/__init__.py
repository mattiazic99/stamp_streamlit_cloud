"""
STAMP Gene Switching Explorer - Utils Module
==========================================

Utility functions for data parsing, analysis, and plotting.
These are the core computational components of the STAMP application.

Author: Gene Analysis Team
"""

__version__ = "2.0.0"
__author__ = "Gene Analysis Team"

# Import utility modules
try:
    from . import parsing
    from . import analysis
    from . import plots
except ImportError as e:
    import warnings
    warnings.warn(f"Some utility modules could not be imported: {e}")
