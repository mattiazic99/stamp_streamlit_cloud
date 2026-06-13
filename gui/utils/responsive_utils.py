"""
Utility functions for responsive figure sizing and adaptive label formatting.
Used across all analysis modules to ensure consistent, window-adaptive chart rendering.
"""

import matplotlib.pyplot as plt


def get_heatmap_fontsize(n_cells):
    """Return annotation fontsize scaled to the number of cells in the matrix.
    
    Args:
        n_cells: max(n_rows, n_cols) of the heatmap matrix
    Returns:
        int fontsize for heatmap cell annotations
    """
    if n_cells <= 5:
        return 14
    elif n_cells <= 8:
        return 12
    elif n_cells <= 12:
        return 10
    elif n_cells <= 18:
        return 9
    else:
        return 7


def get_figsize(base_w=12, base_h=6, n_items=None, per_item=0.4, min_h=4, max_h=20):
    """Return a reasonable figsize tuple, optionally scaling height to the number of items.
    
    Args:
        base_w: base width in inches
        base_h: base height in inches (used when n_items is None)
        n_items: if set, height = max(min_h, n_items * per_item)
        per_item: inches per item for height scaling
        min_h: minimum height
        max_h: maximum height
    Returns:
        (width, height) tuple
    """
    if n_items is not None:
        h = max(min_h, min(max_h, n_items * per_item))
        return (base_w, h)
    return (base_w, base_h)


def auto_rotate_labels(ax, max_label_len=8, rotation=45):
    """Rotate x-axis tick labels if they are longer than max_label_len characters.
    
    Args:
        ax: matplotlib Axes
        max_label_len: character threshold before rotating
        rotation: degrees to rotate
    """
    labels = [t.get_text() for t in ax.get_xticklabels()]
    if any(len(lbl) > max_label_len for lbl in labels):
        ax.set_xticklabels(labels, rotation=rotation, ha='right')
