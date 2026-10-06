"""GTEx tissue labels, anatomical groups and plot colours."""

import re
from typing import Dict, List, Tuple

# Full display-name mapping for the GTEx tissues.
TISSUE_MAPPING = {
    # Adipose tissues
    "Adipose - Subcutaneous_sets_mapped": "Adipose - Subcutaneous",
    "Adipose - Visceral (Omentum)_sets_mapped": "Adipose - Visceral",
    
    # Endocrine system
    "Adrenal Gland_sets_mapped": "Adrenal Gland",
    "Pituitary_sets_mapped": "Pituitary",
    "Thyroid_sets_mapped": "Thyroid",
    
    # Cardiovascular system
    "Artery - Aorta_sets_mapped": "Aorta",
    "Artery - Coronary_sets_mapped": "Coronary Artery",
    "Artery - Tibial_sets_mapped": "Tibial Artery",
    "Heart - Atrial Appendage_sets_mapped": "Heart - Atrium",
    "Heart - Left Ventricle_sets_mapped": "Heart - Ventricle",
    
    # Central nervous system
    "Brain - Amygdala_sets_mapped": "Amygdala",
    "Brain - Anterior cingulate cortex (BA24)_sets_mapped": "Anterior Cingulate",
    "Brain - Caudate (basal ganglia)_sets_mapped": "Caudate",
    "Brain - Cerebellar Hemisphere_sets_mapped": "Cerebellum",
    "Brain - Cerebellum_sets_mapped": "Cerebellum",
    "Brain - Cortex_sets_mapped": "Cerebral Cortex",
    "Brain - Frontal Cortex (BA9)_sets_mapped": "Frontal Cortex",
    "Brain - Hippocampus_sets_mapped": "Hippocampus",
    "Brain - Hypothalamus_sets_mapped": "Hypothalamus",
    "Brain - Nucleus accumbens (basal ganglia)_sets_mapped": "Nucleus Accumbens",
    "Brain - Putamen (basal ganglia)_sets_mapped": "Putamen",
    "Brain - Spinal cord (cervical c-1)_sets_mapped": "Spinal Cord",
    "Brain - Substantia nigra_sets_mapped": "Substantia Nigra",
    
    # Peripheral nervous system
    "Nerve - Tibial_sets_mapped": "Tibial Nerve",
    
    # Reproductive system
    "Breast - Mammary Tissue_sets_mapped": "Breast",
    "Cervix - Ectocervix_sets_mapped": "Cervix - Ecto",
    "Cervix - Endocervix_sets_mapped": "Cervix - Endo",
    "Ovary_sets_mapped": "Ovary",
    "Prostate_sets_mapped": "Prostate",
    "Testis_sets_mapped": "Testis",
    "Uterus_sets_mapped": "Uterus",
    "Vagina_sets_mapped": "Vagina",
    "Fallopian Tube_sets_mapped": "Fallopian Tube",
    
    # Digestive system
    "Colon - Sigmoid_sets_mapped": "Colon - Sigmoid",
    "Colon - Transverse_sets_mapped": "Colon - Transverse",
    "Esophagus - Gastroesophageal Junction_sets_mapped": "Esophagus - GE Junction",
    "Esophagus - Mucosa_sets_mapped": "Esophagus - Mucosa",
    "Esophagus - Muscularis_sets_mapped": "Esophagus - Muscle",
    "Liver_sets_mapped": "Liver",
    "Pancreas_sets_mapped": "Pancreas",
    "Small Intestine - Terminal Ileum_sets_mapped": "Small Intestine",
    "Spleen_sets_mapped": "Spleen",
    "Stomach_sets_mapped": "Stomach",
    
    # Urinary system
    "Bladder_sets_mapped": "Bladder",
    "Kidney - Cortex_sets_mapped": "Kidney - Cortex",
    "Kidney - Medulla_sets_mapped": "Kidney - Medulla",
    
    # Respiratory system
    "Lung_sets_mapped": "Lung",
    
    # Musculoskeletal system
    "Muscle - Skeletal_sets_mapped": "Skeletal Muscle",
    
    # Integumentary system
    "Skin - Not Sun Exposed (Suprapubic)_sets_mapped": "Skin - Protected",
    "Skin - Sun Exposed (Lower leg)_sets_mapped": "Skin - Exposed",
    
    # Endocrine glands
    "Minor Salivary Gland_sets_mapped": "Salivary Gland",
    
    # Cell lines
    "Cells - Cultured fibroblasts_sets_mapped": "Fibroblasts",
    "Cells - EBV-transformed lymphocytes_sets_mapped": "Lymphocytes",
    
    # Blood
    "Whole Blood_sets_mapped": "Whole Blood"
}

# Grouping by anatomical system.
TISSUE_CATEGORIES = {
    "Nervous System": [
        "Amygdala", "Anterior Cingulate", "Caudate", "Cerebellum", 
        "Cerebral Cortex", "Frontal Cortex", "Hippocampus", "Hypothalamus", 
        "Nucleus Accumbens", "Putamen", "Spinal Cord", "Substantia Nigra", "Tibial Nerve"
    ],
    "Cardiovascular": [
        "Aorta", "Coronary Artery", "Tibial Artery", "Heart - Atrium", "Heart - Ventricle"
    ],
    "Digestive System": [
        "Colon - Sigmoid", "Colon - Transverse", "Esophagus - GE Junction", 
        "Esophagus - Mucosa", "Esophagus - Muscle", "Liver", "Pancreas", 
        "Small Intestine", "Spleen", "Stomach"
    ],
    "Reproductive": [
        "Breast", "Cervix - Ecto", "Cervix - Endo", "Ovary", "Prostate", 
        "Testis", "Uterus", "Vagina", "Fallopian Tube"
    ],
    "Urinary": [
        "Bladder", "Kidney - Cortex", "Kidney - Medulla"
    ],
    "Endocrine": [
        "Adrenal Gland", "Pituitary", "Thyroid", "Salivary Gland"
    ],
    "Other": [
        "Adipose - Subcutaneous", "Adipose - Visceral", "Lung", "Skeletal Muscle",
        "Skin - Protected", "Skin - Exposed", "Fibroblasts", "Lymphocytes", "Whole Blood"
    ]
}

# One colour per category, so figures stay consistent across pages.
CATEGORY_COLORS = {
    "Nervous System": "#8E44AD",      # Purple
    "Cardiovascular": "#E74C3C",      # Red
    "Digestive System": "#F39C12",    # Orange
    "Reproductive": "#E91E63",        # Pink
    "Urinary": "#3498DB",             # Blue
    "Endocrine": "#9B59B6",           # Light purple
    "Other": "#95A5A6"                # Grey
}

def clean_tissue_name(filename: str) -> str:
    """Return the mapped tissue label, or clean the filename if it is unknown."""
    # Remove common file extensions.
    clean_name = filename.replace('.txt', '').replace('.csv', '')
    
    # Remove pipeline suffixes.
    suffixes_to_remove = ['_sets_mapped', '_sets', '_mapped', '_data', '_switching']
    for suffix in suffixes_to_remove:
        clean_name = clean_name.replace(suffix, '')
    
    # Prefer the explicit tissue mapping.
    if filename in TISSUE_MAPPING:
        return TISSUE_MAPPING[filename]
    
    # Fall back to automatic cleanup for names outside the mapping.
    clean_name = clean_name.replace('_', ' ').replace('-', ' - ')
    
    clean_name = re.sub(r'\s+', ' ', clean_name)  # Collapse whitespace.
    clean_name = clean_name.strip()
    
    # Keep articles and prepositions lowercase within a label.
    words = clean_name.split()
    capitalized_words = []
    
    for i, word in enumerate(words):
        if i == 0 or word.lower() not in ['and', 'or', 'the', 'of', 'in', 'on', 'at', 'to', 'for']:
            capitalized_words.append(word.capitalize())
        else:
            capitalized_words.append(word.lower())
    
    return ' '.join(capitalized_words)

def get_tissue_category(tissue_name: str) -> str:
    """Return the anatomical category, or Other for an unmapped tissue."""
    for category, tissues in TISSUE_CATEGORIES.items():
        if tissue_name in tissues:
            return category
    return "Other"

def get_category_color(tissue_name: str) -> str:
    """Return the hex colour for the tissue category."""
    category = get_tissue_category(tissue_name)
    return CATEGORY_COLORS.get(category, "#95A5A6")

def create_tissue_summary() -> Dict:
    """Return tissue counts, mapped labels and available categories."""
    total_tissues = len(TISSUE_MAPPING)
    categories_count = {cat: len(tissues) for cat, tissues in TISSUE_CATEGORIES.items()}
    
    return {
        "total_tissues": total_tissues,
        "categories": categories_count,
        "mapped_tissues": list(TISSUE_MAPPING.values()),
        "categories_available": list(TISSUE_CATEGORIES.keys())
    }

def batch_clean_tissue_names(filenames: List[str]) -> Dict[str, str]:
    """Return a mapping from filenames to readable tissue labels."""
    return {filename: clean_tissue_name(filename) for filename in filenames}

def get_tissues_by_category() -> Dict[str, List[str]]:
    """Return a shallow copy of the tissue categories."""
    return TISSUE_CATEGORIES.copy()

def suggest_short_names(max_length: int = 15) -> Dict[str, str]:
    """Return abbreviated labels for tissue names longer than max_length."""
    abbreviations = {}
    
    for original, clean in TISSUE_MAPPING.items():
        if len(clean) > max_length:
            words = clean.split()
            if len(words) == 1:
                # Truncate a single long word.
                abbrev = clean[:max_length-3] + "..."
            elif " - " in clean:
                # Keep the main tissue name and abbreviate the subtype.
                parts = clean.split(" - ")
                main = parts[0]
                sub = parts[1][:3] if len(parts[1]) > 3 else parts[1]
                abbrev = f"{main} - {sub}"
            else:
                # Abbreviate each word to its first three letters.
                abbrev = ' '.join([word[:3] for word in words])
            
            if len(abbrev) > max_length:
                abbrev = abbrev[:max_length-3] + "..."
                
            abbreviations[clean] = abbrev
    
    return abbreviations

if __name__ == "__main__":
    # Manual smoke test.
    test_files = [
        "Brain - Anterior cingulate cortex (BA24)_sets_mapped.txt",
        "Skin - Not Sun Exposed (Suprapubic)_sets_mapped.txt",
        "Heart - Left Ventricle_sets_mapped.txt"
    ]
    
    print("🧬 STAMP Tissue Mapping System")
    print("=" * 50)
    
    for file in test_files:
        clean = clean_tissue_name(file)
        category = get_tissue_category(clean)
        color = get_category_color(clean)
        
        print(f"📁 {file}")
        print(f"🏷️  {clean}")
        print(f"📂 {category}")
        print(f"🎨 {color}")
        print("-" * 30)
    
    summary = create_tissue_summary()
    print(f"\n📊 Summary: {summary['total_tissues']} tissues mapped")
    for cat, count in summary['categories'].items():
        print(f"   {cat}: {count} tissues")
