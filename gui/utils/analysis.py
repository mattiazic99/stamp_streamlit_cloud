import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom
import pandas as pd

def jaccard(set1, set2):
    """
    Compute Jaccard similarity between two sets
    
    Args:
        set1, set2: Sets to compare
        
    Returns:
        float: Jaccard similarity coefficient (0-1)
    """
    inter = len(set1 & set2)
    union = len(set1 | set2)
    return inter / union if union else 0.0

def compute_jaccard_matrix(data, mode="age"):
    """
    Compute the Jaccard similarity matrix between tissues
    
    Args:
        data: Dictionary with tissue names as keys and list of gene sets as values
        mode: "age" = average across 5 age groups, "life" = union of all age groups
        
    Returns:
        tuple: (similarity_matrix, tissue_names)
    """
    tissues = list(data.keys())
    n = len(tissues)
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if mode == "age":
                # Average similarity across age groups
                sims = []
                for k in range(min(5, len(data[tissues[i]]))):
                    if k < len(data[tissues[j]]):
                        sim = jaccard(data[tissues[i]][k], data[tissues[j]][k])
                        sims.append(sim)
                matrix[i, j] = np.mean(sims) if sims else 0.0
            elif mode == "life":
                # Whole life similarity
                union_i = set.union(*data[tissues[i]]) if data[tissues[i]] else set()
                union_j = set.union(*data[tissues[j]]) if data[tissues[j]] else set()
                matrix[i, j] = jaccard(union_i, union_j)

    return matrix, tissues

def compute_common_genes_matrix(data, idx=None):
    """
    Compute number of genes shared between each pair of tissues
    
    Args:
        data: Dictionary with tissue names as keys and list of gene sets as values
        idx: If given, works on a specific age group (0–4), else uses all data
        
    Returns:
        tuple: (shared_genes_matrix, tissue_names)
    """
    tissues = list(data.keys())
    n = len(tissues)
    matrix = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(n):
            if idx is None:
                # Whole lifespan
                genes_i = set.union(*data[tissues[i]]) if data[tissues[i]] else set()
                genes_j = set.union(*data[tissues[j]]) if data[tissues[j]] else set()
            else:
                # Specific age group
                genes_i = data[tissues[i]][idx] if idx < len(data[tissues[i]]) else set()
                genes_j = data[tissues[j]][idx] if idx < len(data[tissues[j]]) else set()
            
            matrix[i, j] = len(genes_i & genes_j)

    return matrix, tissues

def compute_percent_overlap_matrix(data, idx=None):
    """
    Compute percentage of overlapping genes between tissues
    
    Args:
        data: Dictionary with tissue names as keys and list of gene sets as values
        idx: If given, works on a specific age group (0–4), else uses all data
        
    Returns:
        tuple: (overlap_percentage_matrix, tissue_names)
    """
    tissues = list(data.keys())
    n = len(tissues)
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if idx is None:
                # Whole lifespan
                genes_i = set.union(*data[tissues[i]]) if data[tissues[i]] else set()
                genes_j = set.union(*data[tissues[j]]) if data[tissues[j]] else set()
            else:
                # Specific age group
                genes_i = data[tissues[i]][idx] if idx < len(data[tissues[i]]) else set()
                genes_j = data[tissues[j]][idx] if idx < len(data[tissues[j]]) else set()
            
            inter = len(genes_i & genes_j)
            union = len(genes_i | genes_j)
            matrix[i, j] = round((inter / union * 100), 2) if union else 0.0

    return matrix, tissues

def compute_linkage(matrix, method='average'):
    """
    Perform hierarchical clustering based on the similarity matrix
    
    Args:
        matrix: Similarity matrix
        method: Linkage method ('average', 'complete', 'single', 'ward')
        
    Returns:
        ndarray or None: Linkage matrix for dendrogram plotting
    """
    try:
        # Convert similarity to distance
        dist_matrix = 1 - matrix
        np.fill_diagonal(dist_matrix, 0)
        
        # Ensure the matrix is valid
        if np.any(dist_matrix < 0):
            dist_matrix = np.clip(dist_matrix, 0, None)
        
        # Convert to condensed form
        condensed = squareform(dist_matrix)
        
        # Perform hierarchical clustering
        return linkage(condensed, method=method)
        
    except Exception as e:
        print(f"Error in linkage computation: {e}")
        return None

def compute_overlap_coefficient(set1, set2):
    """
    Compute overlap coefficient (Szymkiewicz–Simpson coefficient)
    
    Args:
        set1, set2: Sets to compare
        
    Returns:
        float: Overlap coefficient (0-1)
    """
    inter = len(set1 & set2)
    min_size = min(len(set1), len(set2))
    return inter / min_size if min_size > 0 else 0.0

def compute_dice_coefficient(set1, set2):
    """
    Compute Dice coefficient (Sørensen–Dice coefficient)
    
    Args:
        set1, set2: Sets to compare
        
    Returns:
        float: Dice coefficient (0-1)
    """
    inter = len(set1 & set2)
    total = len(set1) + len(set2)
    return (2 * inter) / total if total > 0 else 0.0

def compute_statistical_significance(group1_genes, group2_genes, universe_genes):
    """
    Compute statistical significance of gene overlap using hypergeometric test
    
    Args:
        group1_genes: Set of genes in group 1
        group2_genes: Set of genes in group 2
        universe_genes: Set of all possible genes (background)
        
    Returns:
        dict: Statistical test results
    """
    overlap = len(group1_genes & group2_genes)
    group1_size = len(group1_genes)
    group2_size = len(group2_genes)
    universe_size = len(universe_genes)
    
    # Hypergeometric test
    p_value = 1 - hypergeom.cdf(overlap - 1, universe_size, group1_size, group2_size)
    
    # Expected overlap under null hypothesis
    expected_overlap = (group1_size * group2_size) / universe_size if universe_size > 0 else 0
    
    # Fold enrichment
    fold_enrichment = overlap / expected_overlap if expected_overlap > 0 else 0
    
    return {
        'overlap_observed': overlap,
        'overlap_expected': expected_overlap,
        'fold_enrichment': fold_enrichment,
        'p_value': p_value,
        'significant': p_value < 0.05
    }

def compute_diversity_metrics(gene_sets):
    """
    Compute diversity metrics for a collection of gene sets
    
    Args:
        gene_sets: List of gene sets
        
    Returns:
        dict: Diversity metrics
    """
    # Shannon diversity
    total_appearances = {}
    total_sets = len(gene_sets)
    
    for gene_set in gene_sets:
        for gene in gene_set:
            total_appearances[gene] = total_appearances.get(gene, 0) + 1
    
    if not total_appearances:
        return {'shannon_diversity': 0, 'simpson_diversity': 0, 'total_unique_genes': 0}
    
    # Calculate Shannon diversity
    total_genes = sum(total_appearances.values())
    shannon = 0
    simpson = 0
    
    for count in total_appearances.values():
        p = count / total_genes
        if p > 0:
            shannon -= p * np.log(p)
            simpson += p * p
    
    simpson_diversity = 1 - simpson  # Simpson's diversity index
    
    return {
        'shannon_diversity': shannon,
        'simpson_diversity': simpson_diversity,
        'total_unique_genes': len(total_appearances),
        'avg_genes_per_set': np.mean([len(s) for s in gene_sets]),
        'std_genes_per_set': np.std([len(s) for s in gene_sets])
    }

def compute_gene_frequency_analysis(data, min_frequency=1):
    """
    Analyze gene frequency across tissues and age groups
    
    Args:
        data: Dictionary with tissue names as keys and list of gene sets as values
        min_frequency: Minimum frequency threshold for reporting
        
    Returns:
        dict: Gene frequency analysis results
    """
    gene_frequency = {}
    tissue_gene_frequency = {}
    age_gene_frequency = {}
    
    # Count gene appearances
    for tissue, age_sets in data.items():
        tissue_gene_frequency[tissue] = {}
        
        for age_idx, gene_set in enumerate(age_sets):
            for gene in gene_set:
                # Global frequency
                gene_frequency[gene] = gene_frequency.get(gene, 0) + 1
                
                # Tissue-specific frequency
                tissue_gene_frequency[tissue][gene] = tissue_gene_frequency[tissue].get(gene, 0) + 1
                
                # Age-specific frequency
                if age_idx not in age_gene_frequency:
                    age_gene_frequency[age_idx] = {}
                age_gene_frequency[age_idx][gene] = age_gene_frequency[age_idx].get(gene, 0) + 1
    
    # Filter by minimum frequency
    frequent_genes = {gene: freq for gene, freq in gene_frequency.items() if freq >= min_frequency}
    
    # Sort by frequency
    sorted_genes = sorted(frequent_genes.items(), key=lambda x: x[1], reverse=True)
    
    return {
        'gene_frequency': dict(sorted_genes),
        'tissue_frequency': tissue_gene_frequency,
        'age_frequency': age_gene_frequency,
        'total_unique_genes': len(gene_frequency),
        'frequent_genes_count': len(frequent_genes),
        'most_frequent_gene': sorted_genes[0] if sorted_genes else None,
        'frequency_distribution': {
            'mean': np.mean(list(gene_frequency.values())),
            'std': np.std(list(gene_frequency.values())),
            'median': np.median(list(gene_frequency.values())),
            'max': max(gene_frequency.values()) if gene_frequency else 0,
            'min': min(gene_frequency.values()) if gene_frequency else 0
        }
    }

def compute_age_progression_analysis(data):
    """
    Analyze how gene expression changes with age progression
    
    Args:
        data: Dictionary with tissue names as keys and list of gene sets as values
        
    Returns:
        dict: Age progression analysis results
    """
    age_groups = ["30–39", "40–49", "50–59", "60–69", "70–79"]
    
    results = {}
    
    for tissue, age_sets in data.items():
        tissue_results = {
            'age_gene_counts': [],
            'age_unique_genes': [],
            'cumulative_genes': set(),
            'new_genes_per_age': [],
            'lost_genes_per_age': [],
            'persistent_genes': None
        }
        
        previous_genes = set()
        cumulative = set()
        
        for age_idx, gene_set in enumerate(age_sets):
            # Basic counts
            tissue_results['age_gene_counts'].append(len(gene_set))
            tissue_results['age_unique_genes'].append(gene_set)
            
            # Cumulative genes
            cumulative.update(gene_set)
            
            # New and lost genes
            if age_idx > 0:
                new_genes = gene_set - previous_genes
                lost_genes = previous_genes - gene_set
                tissue_results['new_genes_per_age'].append(len(new_genes))
                tissue_results['lost_genes_per_age'].append(len(lost_genes))
            else:
                tissue_results['new_genes_per_age'].append(len(gene_set))
                tissue_results['lost_genes_per_age'].append(0)
            
            previous_genes = gene_set.copy()
        
        # Genes present in all age groups (persistent)
        if age_sets:
            persistent = set.intersection(*[s for s in age_sets if s])
            tissue_results['persistent_genes'] = persistent
        
        tissue_results['cumulative_genes'] = cumulative
        
        # Calculate trends
        counts = tissue_results['age_gene_counts']
        if len(counts) > 1:
            # Linear trend (slope)
            x = np.arange(len(counts))
            trend_slope = np.polyfit(x, counts, 1)[0]
            tissue_results['trend_slope'] = trend_slope
            tissue_results['trend_direction'] = 'increasing' if trend_slope > 0 else 'decreasing' if trend_slope < 0 else 'stable'
        
        results[tissue] = tissue_results
    
    return results

def compute_tissue_similarity_ranking(data, reference_tissue, metric='jaccard'):
    """
    Rank tissues by similarity to a reference tissue
    
    Args:
        data: Dictionary with tissue names as keys and list of gene sets as values
        reference_tissue: Name of reference tissue
        metric: Similarity metric ('jaccard', 'overlap', 'dice')
        
    Returns:
        list: Ranked list of (tissue_name, similarity_score) tuples
    """
    if reference_tissue not in data:
        return []
    
    ref_genes = set.union(*data[reference_tissue]) if data[reference_tissue] else set()
    similarities = []
    
    for tissue, age_sets in data.items():
        if tissue == reference_tissue:
            continue
            
        tissue_genes = set.union(*age_sets) if age_sets else set()
        
        if metric == 'jaccard':
            similarity = jaccard(ref_genes, tissue_genes)
        elif metric == 'overlap':
            similarity = compute_overlap_coefficient(ref_genes, tissue_genes)
        elif metric == 'dice':
            similarity = compute_dice_coefficient(ref_genes, tissue_genes)
        else:
            similarity = jaccard(ref_genes, tissue_genes)  # Default to Jaccard
        
        similarities.append((tissue, similarity))
    
    # Sort by similarity (descending)
    similarities.sort(key=lambda x: x[1], reverse=True)
    
    return similarities

def compute_core_genes(tissue_groups, min_tissues=None):
    """
    Find core genes present in multiple tissues
    
    Args:
        tissue_groups: Dictionary or list of tissue gene sets
        min_tissues: Minimum number of tissues gene must appear in
        
    Returns:
        dict: Core gene analysis results
    """
    if isinstance(tissue_groups, dict):
        tissue_gene_sets = []
        tissue_names = []
        for tissue, age_sets in tissue_groups.items():
            tissue_gene_sets.append(set.union(*age_sets) if age_sets else set())
            tissue_names.append(tissue)
    else:
        tissue_gene_sets = tissue_groups
        tissue_names = [f"Tissue_{i}" for i in range(len(tissue_groups))]
    
    if min_tissues is None:
        min_tissues = len(tissue_gene_sets) // 2  # Default to half
    
    # Count gene appearances across tissues
    gene_counts = {}
    for gene_set in tissue_gene_sets:
        for gene in gene_set:
            gene_counts[gene] = gene_counts.get(gene, 0) + 1
    
    # Filter core genes
    core_genes = {gene: count for gene, count in gene_counts.items() if count >= min_tissues}
    
    # Universal genes (present in ALL tissues)
    universal_genes = set.intersection(*tissue_gene_sets) if tissue_gene_sets else set()
    
    return {
        'core_genes': core_genes,
        'universal_genes': universal_genes,
        'gene_frequency_distribution': gene_counts,
        'total_unique_genes': len(gene_counts),
        'core_genes_count': len(core_genes),
        'universal_genes_count': len(universal_genes),
        'tissue_count': len(tissue_gene_sets)
    }

def compute_enrichment_analysis(target_genes, background_genes, gene_sets, labels=None):
    """
    Perform enrichment analysis of target genes in various gene sets
    
    Args:
        target_genes: Set of genes of interest
        background_genes: Set of background genes (universe)
        gene_sets: List of gene sets to test for enrichment
        labels: Optional labels for gene sets
        
    Returns:
        list: Enrichment results for each gene set
    """
    if labels is None:
        labels = [f"Set_{i}" for i in range(len(gene_sets))]
    
    results = []
    
    for i, gene_set in enumerate(gene_sets):
        label = labels[i] if i < len(labels) else f"Set_{i}"
        
        # Calculate enrichment
        overlap = len(target_genes & gene_set)
        target_size = len(target_genes)
        set_size = len(gene_set)
        background_size = len(background_genes)
        
        # Hypergeometric test
        p_value = 1 - hypergeom.cdf(overlap - 1, background_size, set_size, target_size)
        
        # Expected overlap
        expected = (target_size * set_size) / background_size if background_size > 0 else 0
        
        # Fold enrichment
        fold_enrichment = overlap / expected if expected > 0 else 0
        
        results.append({
            'label': label,
            'overlap': overlap,
            'expected': expected,
            'fold_enrichment': fold_enrichment,
            'p_value': p_value,
            'significant': p_value < 0.05,
            'target_size': target_size,
            'set_size': set_size,
            'enrichment_score': -np.log10(p_value) if p_value > 0 else 0
        })
    
    # Sort by significance
    results.sort(key=lambda x: x['p_value'])
    
    return results