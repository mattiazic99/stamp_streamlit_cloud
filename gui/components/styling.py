import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def apply_plot_style():
    """Apply consistent styling to all plots"""
    # Set the style
    plt.style.use('seaborn-v0_8-whitegrid')
    
    # Set color palette
    sns.set_palette("husl")
    
    # Configure matplotlib parameters
    plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'axes.edgecolor': '#333333',
        'axes.linewidth': 1.2,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linewidth': 0.8,
        'font.size': 10,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 16,
        'text.color': 'black',
    })

def get_color_palette(n_colors, palette_name="husl"):
    """Get a color palette with n colors"""
    if palette_name == "husl":
        return sns.color_palette("husl", n_colors)
    elif palette_name == "viridis":
        return sns.color_palette("viridis", n_colors)
    elif palette_name == "Set2":
        return sns.color_palette("Set2", n_colors)
    else:
        return sns.color_palette(palette_name, n_colors)

def format_heatmap(ax, title, xlabel="", ylabel="", cbar_label=""):
    """Format heatmap with consistent styling"""
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    
    # Rotate x-axis labels if they're long
    if len(ax.get_xticklabels()) > 0:
        max_label_length = max([len(label.get_text()) for label in ax.get_xticklabels()])
        if max_label_length > 8:
            plt.xticks(rotation=45, ha='right')
    
    return ax

def format_barplot(ax, title, xlabel="", ylabel="", add_values=True):
    """Format bar plot with consistent styling"""
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    
    # Add value labels on bars if requested
    if add_values:
        for container in ax.containers:
            ax.bar_label(container, fontweight='bold')
    
    ax.grid(True, alpha=0.3, axis='y')
    return ax

def create_styled_figure(figsize=(10, 6)):
    """Create a figure with consistent styling"""
    apply_plot_style()
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax

def save_figure_high_quality(fig, filename, dpi=300):
    """Save figure with high quality settings"""
    fig.savefig(filename, dpi=dpi, bbox_inches='tight', 
                facecolor='white', edgecolor='none')