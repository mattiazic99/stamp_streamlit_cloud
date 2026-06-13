import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import streamlit as st
from typing import List, Dict, Tuple, Optional, Union

# Configure global responsive settings for Matplotlib
plt.rcParams.update({
    'font.size': 14,             # Keep baseline text large so it's readable when scaled down
    'axes.titlesize': 18,        # Large titles
    'axes.labelsize': 14,        # Large axis labels
    'xtick.labelsize': 12,       # Large tick labels
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 20,
    'figure.autolayout': True    # Equivalent to tight_layout globally
})

def plot_pie_common_vs_exclusive(common: int, excl1: int, excl2: int, 
                                group1_name: str, group2_name: str,
                                title: str = None, figsize: Tuple[int, int] = (8, 6)):
    """
    Create an enhanced pie chart comparing shared and exclusive genes between two groups
    
    Args:
        common: Number of shared genes
        excl1: Number of genes exclusive to group 1
        excl2: Number of genes exclusive to group 2
        group1_name: Name of first group
        group2_name: Name of second group
        title: Optional custom title
        figsize: Figure size tuple
    """
    if common + excl1 + excl2 == 0:
        st.warning("No genes to display in pie chart.")
        return
    
    # Data preparation
    sizes = [common, excl1, excl2]
    labels = [f"Shared\n({common})", 
              f"{group1_name}\nExclusive\n({excl1})", 
              f"{group2_name}\nExclusive\n({excl2})"]
    colors = ["#2ecc71", "#3498db", "#e74c3c"]  # Green, blue, red
    explode = (0.05, 0, 0)  # Slightly explode the shared portion
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create pie chart with enhanced styling
    wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', 
                                     startangle=90, colors=colors, explode=explode,
                                     shadow=True, textprops={'fontsize': 10, 'weight': 'bold'})
    
    # Enhance the appearance
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_weight('bold')
    
    # Set title
    if title is None:
        title = f"Gene Distribution: {group1_name} vs {group2_name}"
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Equal aspect ratio ensures that pie is drawn as a circle
    ax.axis('equal')
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def plot_similarity_pie(similarities: List[Tuple[str, float]], ref_tissue: str,
                       title: str = None, figsize: Tuple[int, int] = (8, 6)):
    """
    Create an enhanced pie chart showing similarity distribution relative to a reference tissue
    
    Args:
        similarities: List of (tissue_name, similarity_percentage) tuples
        ref_tissue: Name of reference tissue
        title: Optional custom title
        figsize: Figure size tuple
    """
    if not similarities:
        st.warning("No similarity data to display.")
        return
    
    # Prepare data
    labels = [f"{tissue}\n({sim:.1f}%)" for tissue, sim in similarities]
    sizes = [sim for _, sim in similarities]
    
    # Create color palette
    colors = plt.cm.viridis(np.linspace(0, 1, len(similarities)))
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create pie chart
    wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', 
                                     startangle=90, colors=colors,
                                     textprops={'fontsize': 9, 'weight': 'bold'})
    
    # Enhance appearance
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_weight('bold')
        autotext.set_fontsize(8)
    
    # Set title
    if title is None:
        title = f"Similarity Distribution Relative to '{ref_tissue}'"
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    ax.axis('equal')
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def plot_enhanced_heatmap(matrix: np.ndarray, xticklabels: List[str], yticklabels: List[str],
                         title: str, cmap: str = "coolwarm", annot: bool = True,
                         fmt: str = ".2f", figsize: Tuple[int, int] = (8, 6),
                         mask_upper: bool = False, cbar_label: str = "Value"):
    """
    Create an enhanced heatmap with customizable options
    
    Args:
        matrix: 2D numpy array of values
        xticklabels: Labels for x-axis
        yticklabels: Labels for y-axis  
        title: Plot title
        cmap: Colormap name
        annot: Whether to annotate cells
        fmt: Format string for annotations
        figsize: Figure size tuple
        mask_upper: Whether to mask upper triangle
        cbar_label: Label for colorbar
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create mask if requested
    mask = None
    if mask_upper:
        mask = np.triu(np.ones_like(matrix, dtype=bool), k=1)
    
    # Create heatmap
    heatmap = sns.heatmap(matrix, 
                         annot=annot, 
                         fmt=fmt,
                         xticklabels=xticklabels,
                         yticklabels=yticklabels,
                         cmap=cmap,
                         mask=mask,
                         square=True,
                         linewidths=0.5,
                         cbar_kws={'label': cbar_label},
                         ax=ax)
    
    # Enhance appearance
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    # Rotate labels if they're long
    if any(len(label) > 8 for label in xticklabels):
        plt.xticks(rotation=45, ha='right')
    if any(len(label) > 8 for label in yticklabels):
        plt.yticks(rotation=0)
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def plot_gene_count_comparison(tissues: List[str], counts1: List[int], counts2: List[int] = None,
                              labels: List[str] = None, title: str = "Gene Count Comparison",
                              figsize: Tuple[int, int] = (10, 5), colors: List[str] = None):
    """
    Create an enhanced bar chart comparing gene counts across tissues
    
    Args:
        tissues: List of tissue names
        counts1: Gene counts for first condition/group
        counts2: Optional gene counts for second condition/group
        labels: Labels for the conditions
        title: Plot title
        figsize: Figure size tuple
        colors: Custom colors for bars
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    x = np.arange(len(tissues))
    width = 0.35 if counts2 is not None else 0.6
    
    # Default colors
    if colors is None:
        colors = ['#3498db', '#e74c3c'] if counts2 is not None else ['#3498db']
    
    # Default labels
    if labels is None:
        labels = ['Condition 1', 'Condition 2'] if counts2 is not None else ['Gene Count']
    
    # Create bars
    bars1 = ax.bar(x - width/2 if counts2 is not None else x, counts1, width,
                   label=labels[0], color=colors[0], alpha=0.8, edgecolor='black', linewidth=1.2)
    
    if counts2 is not None:
        bars2 = ax.bar(x + width/2, counts2, width,
                       label=labels[1], color=colors[1], alpha=0.8, edgecolor='black', linewidth=1.2)
    
    # Add value labels on bars
    def add_value_labels(bars):
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(counts1 + (counts2 or []))*0.01,
                   f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    add_value_labels(bars1)
    if counts2 is not None:
        add_value_labels(bars2)
    
    # Enhance appearance
    ax.set_xlabel('Tissue', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(tissues, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    if counts2 is not None:
        ax.legend(frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def plot_age_progression(age_groups: List[str], values: List[float], 
                        tissue_name: str = "Tissue", metric: str = "Gene Count",
                        figsize: Tuple[int, int] = (8, 5), color: str = '#2ecc71'):
    """
    Create an enhanced line plot showing progression across age groups
    
    Args:
        age_groups: List of age group labels
        values: Values for each age group
        tissue_name: Name of tissue being analyzed
        metric: Name of metric being plotted
        figsize: Figure size tuple
        color: Line color
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create line plot with filled area
    ax.plot(age_groups, values, marker='o', linewidth=3, markersize=8, 
           color=color, markerfacecolor='white', markeredgecolor=color, markeredgewidth=2)
    ax.fill_between(age_groups, values, alpha=0.3, color=color)
    
    # Add value labels
    for i, value in enumerate(values):
        ax.text(i, value + max(values)*0.02, f'{value:.1f}', 
               ha='center', va='bottom', fontweight='bold')
    
    # Enhance appearance
    ax.set_title(f"{metric} Progression: {tissue_name}", fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Age Group', fontsize=12, fontweight='bold')
    ax.set_ylabel(metric, fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, max(values) * 1.1 if values else 1)
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def plot_venn_diagram_data(set1: set, set2: set, set1_name: str, set2_name: str,
                          figsize: Tuple[int, int] = (8, 5)):
    """
    Create a visual representation of Venn diagram data using bar charts
    (Since matplotlib-venn is not available, we use bar charts to show overlap)
    
    Args:
        set1: First set of items
        set2: Second set of items
        set1_name: Name of first set
        set2_name: Name of second set
        figsize: Figure size tuple
    """
    # Calculate overlaps
    intersection = len(set1 & set2)
    only_set1 = len(set1 - set2)
    only_set2 = len(set2 - set1)
    
    # Create bar chart representation
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Bar chart of overlaps
    categories = ['Shared', f'Only {set1_name}', f'Only {set2_name}']
    values = [intersection, only_set1, only_set2]
    colors = ['#2ecc71', '#3498db', '#e74c3c']
    
    bars = ax1.bar(categories, values, color=colors, alpha=0.8, edgecolor='black')
    ax1.set_title('Gene Set Overlap', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + max(values)*0.01,
                f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    # Pie chart representation
    if sum(values) > 0:
        ax2.pie(values, labels=categories, autopct='%1.1f%%', colors=colors, startangle=90)
        ax2.set_title('Overlap Proportions', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def plot_similarity_matrix_with_clustering(matrix: np.ndarray, labels: List[str],
                                          linkage_matrix: np.ndarray = None,
                                          title: str = "Similarity Matrix with Clustering",
                                          figsize: Tuple[int, int] = (10, 8)):
    """
    Create a similarity matrix plot with optional dendrogram
    
    Args:
        matrix: Similarity matrix
        labels: Labels for matrix rows/columns
        linkage_matrix: Optional linkage matrix for dendrogram
        title: Plot title
        figsize: Figure size tuple
    """
    if linkage_matrix is not None:
        # Create figure with dendrogram
        fig = plt.figure(figsize=figsize)
        
        # Create grid
        gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[1, 4],
                             hspace=0.05, wspace=0.05)
        
        # Top dendrogram
        ax_dendro_top = fig.add_subplot(gs[0, 1])
        from scipy.cluster.hierarchy import dendrogram
        dendro_top = dendrogram(linkage_matrix, ax=ax_dendro_top, orientation='top',
                               labels=labels, leaf_rotation=90)
        ax_dendro_top.set_xticks([])
        ax_dendro_top.set_yticks([])
        
        # Left dendrogram
        ax_dendro_left = fig.add_subplot(gs[1, 0])
        dendro_left = dendrogram(linkage_matrix, ax=ax_dendro_left, orientation='left',
                                labels=labels)
        ax_dendro_left.set_xticks([])
        ax_dendro_left.set_yticks([])
        
        # Main heatmap
        ax_heatmap = fig.add_subplot(gs[1, 1])
        
        # Reorder matrix according to dendrogram
        order = dendro_top['leaves']
        matrix_ordered = matrix[np.ix_(order, order)]
        labels_ordered = [labels[i] for i in order]
        
    else:
        fig, ax_heatmap = plt.subplots(figsize=figsize)
        matrix_ordered = matrix
        labels_ordered = labels
    
    # Create heatmap
    im = ax_heatmap.imshow(matrix_ordered, cmap='coolwarm', aspect='equal')
    
    # Set ticks and labels
    ax_heatmap.set_xticks(range(len(labels_ordered)))
    ax_heatmap.set_yticks(range(len(labels_ordered)))
    ax_heatmap.set_xticklabels(labels_ordered, rotation=45, ha='right')
    ax_heatmap.set_yticklabels(labels_ordered)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax_heatmap)
    cbar.set_label('Similarity', fontsize=12, fontweight='bold')
    
    # Add value annotations
    for i in range(len(labels_ordered)):
        for j in range(len(labels_ordered)):
            text = ax_heatmap.text(j, i, f'{matrix_ordered[i, j]:.2f}',
                                 ha="center", va="center", color="black", fontweight='bold')
    
    ax_heatmap.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def plot_gene_frequency_distribution(gene_frequencies: Dict[str, int], 
                                   top_n: int = 20, title: str = "Gene Frequency Distribution",
                                   figsize: Tuple[int, int] = (10, 6)):
    """
    Create a bar plot showing gene frequency distribution
    
    Args:
        gene_frequencies: Dictionary mapping gene names to frequencies
        top_n: Number of top genes to show
        title: Plot title
        figsize: Figure size tuple
    """
    if not gene_frequencies:
        st.warning("No gene frequency data to display.")
        return
    
    # Sort genes by frequency
    sorted_genes = sorted(gene_frequencies.items(), key=lambda x: x[1], reverse=True)[:top_n]
    
    genes = [item[0] for item in sorted_genes]
    frequencies = [item[1] for item in sorted_genes]
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create bars with gradient colors
    colors = plt.cm.viridis(np.linspace(0, 1, len(genes)))
    bars = ax.bar(range(len(genes)), frequencies, color=colors, alpha=0.8, edgecolor='black')
    
    # Customize appearance
    ax.set_xlabel('Genes', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + max(frequencies)*0.01,
               f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def plot_correlation_heatmap(correlation_matrix: np.ndarray, labels: List[str],
                           title: str = "Correlation Matrix", figsize: Tuple[int, int] = (8, 6)):
    """
    Create a correlation heatmap with enhanced styling
    
    Args:
        correlation_matrix: Correlation matrix
        labels: Labels for matrix rows/columns
        title: Plot title
        figsize: Figure size tuple
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap with diverging colormap
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool), k=1)
    
    heatmap = sns.heatmap(correlation_matrix,
                         mask=mask,
                         annot=True,
                         fmt='.3f',
                         xticklabels=labels,
                         yticklabels=labels,
                         cmap='RdBu_r',
                         center=0,
                         square=True,
                         linewidths=0.5,
                         cbar_kws={'label': 'Correlation Coefficient'},
                         ax=ax)
    
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    # Rotate labels if needed
    if any(len(label) > 8 for label in labels):
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()

def create_summary_dashboard(data_summary: Dict, figsize: Tuple[int, int] = (14, 9)):
    """
    Create a comprehensive dashboard summarizing the analysis
    
    Args:
        data_summary: Dictionary containing summary statistics
        figsize: Figure size tuple
    """
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    fig.suptitle('STAMP Analysis Summary Dashboard', fontsize=20, fontweight='bold')
    
    # Flatten axes for easier indexing
    axes = axes.flatten()
    
    # Plot 1: Tissue gene counts
    if 'tissue_gene_counts' in data_summary:
        tissues = list(data_summary['tissue_gene_counts'].keys())
        counts = list(data_summary['tissue_gene_counts'].values())
        
        axes[0].bar(tissues, counts, color='skyblue', alpha=0.8)
        axes[0].set_title('Genes per Tissue', fontweight='bold')
        axes[0].set_ylabel('Gene Count')
        axes[0].tick_params(axis='x', rotation=45)
    
    # Plot 2: Age group distribution
    if 'age_gene_counts' in data_summary:
        age_groups = list(data_summary['age_gene_counts'].keys())
        counts = list(data_summary['age_gene_counts'].values())
        
        axes[1].bar(age_groups, counts, color='lightcoral', alpha=0.8)
        axes[1].set_title('Genes per Age Group', fontweight='bold')
        axes[1].set_ylabel('Gene Count')
    
    # Plot 3: Similarity distribution
    if 'similarity_distribution' in data_summary:
        similarities = data_summary['similarity_distribution']
        axes[2].hist(similarities, bins=20, color='lightgreen', alpha=0.8, edgecolor='black')
        axes[2].set_title('Similarity Distribution', fontweight='bold')
        axes[2].set_xlabel('Jaccard Similarity')
        axes[2].set_ylabel('Frequency')
    
    # Plot 4: Top genes
    if 'top_genes' in data_summary:
        genes = list(data_summary['top_genes'].keys())[:10]
        frequencies = list(data_summary['top_genes'].values())[:10]
        
        axes[3].barh(genes, frequencies, color='gold', alpha=0.8)
        axes[3].set_title('Top 10 Most Frequent Genes', fontweight='bold')
        axes[3].set_xlabel('Frequency')
    
    # Plot 5: Dataset statistics
    if 'dataset_stats' in data_summary:
        stats = data_summary['dataset_stats']
        labels = list(stats.keys())
        values = list(stats.values())
        
        axes[4].pie(values, labels=labels, autopct='%1.1f%%', startangle=90)
        axes[4].set_title('Dataset Composition', fontweight='bold')
    
    # Plot 6: Summary metrics
    if 'summary_metrics' in data_summary:
        metrics = data_summary['summary_metrics']
        metric_names = list(metrics.keys())
        metric_values = list(metrics.values())
        
        bars = axes[5].bar(metric_names, metric_values, color='mediumpurple', alpha=0.8)
        axes[5].set_title('Summary Metrics', fontweight='bold')
        axes[5].tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            axes[5].text(bar.get_x() + bar.get_width()/2., height + max(metric_values)*0.01,
                        f'{height:.1f}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)
    plt.close()