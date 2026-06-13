import os
import numpy as np
import pandas as pd
from tqdm import tqdm


# =========================
# CONFIGURAZIONE
# =========================

VERSION = "v10"
TAU = 0.5
N_PERMUTATIONS = 100
RANDOM_SEED = 42

TPM_PATH = f"data/parquet/{VERSION}/tpm_matrix.parquet"
META_PATH = f"data/parquet/{VERSION}/metadata.parquet"

OUTPUT_DIR = f"output/{VERSION}/validation/permutation_sample_level"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "sample_level_permutation_test_all_tissues_100perm.csv")

AGE_COLUMNS = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]

SELECTED_TISSUES = None

# =========================
# FUNZIONI STAMP
# =========================

def is_switching_binary_vector(values):
    """
    True se il vettore binario ha esattamente una transizione.
    """
    changes = 0

    for i in range(len(values) - 1):
        if values[i] != values[i + 1]:
            changes += 1

    return changes == 1


def count_switching_genes(age_mean_df, tau):
    """
    Conta i geni switching in una matrice gene x age_bracket.
    """
    age_cols = [col for col in AGE_COLUMNS if col in age_mean_df.columns]

    if len(age_cols) < 2:
        return 0

    df = age_mean_df[age_cols].dropna()

    binary_df = (df >= tau).astype(int)

    switching_mask = binary_df.apply(
        lambda row: is_switching_binary_vector(row.values.tolist()),
        axis=1
    )

    return int(switching_mask.sum())


def normalize_minmax_gene_wise(tissue_tpm):
    """
    Normalizzazione min-max per gene sul sottoinsieme di sample del tessuto.
    Input: gene x sample.
    Output: gene x sample normalizzato in [0, 1].
    """
    min_vals = tissue_tpm.min(axis=1)
    max_vals = tissue_tpm.max(axis=1)
    ranges = max_vals - min_vals

    # Evita divisione per zero
    ranges = ranges.replace(0, np.nan)

    norm = tissue_tpm.sub(min_vals, axis=0).div(ranges, axis=0)

    # I geni costanti diventano NaN: li mettiamo a 0
    norm = norm.fillna(0.0)

    return norm


def compute_age_means(norm_tpm, sample_to_age):
    """
    Calcola media per fascia d'età.
    norm_tpm: gene x sample
    sample_to_age: dict sample_id -> age_bracket
    """
    df_t = norm_tpm.T.copy()
    df_t["age_bracket"] = df_t.index.map(sample_to_age)

    df_t = df_t.dropna(subset=["age_bracket"])

    age_mean_df = df_t.groupby("age_bracket").mean().T

    age_cols = [col for col in AGE_COLUMNS if col in age_mean_df.columns]
    age_mean_df = age_mean_df[age_cols]

    return age_mean_df


def empirical_p_value(observed, null_values):
    """
    P-value empirico one-sided:
    P(null >= observed)
    """
    null_values = np.array(null_values)
    return (np.sum(null_values >= observed) + 1) / (len(null_values) + 1)


def process_tissue(tpm, meta, tissue_name, rng):
    """
    Esegue il sample-level permutation test per un tessuto.
    """
    print(f"\n=== Tissue: {tissue_name} ===")

    meta_tissue = meta[meta["tissue"] == tissue_name].copy()

    if meta_tissue.empty:
        print(f"[SKIP] Nessun metadata per tessuto: {tissue_name}")
        return None

    # Sample del tessuto presenti anche nella matrice TPM
    available_samples = [
        sid for sid in meta_tissue["sample_id"].tolist()
        if sid in tpm.columns
    ]

    if len(available_samples) < 10:
        print(f"[SKIP] Troppi pochi sample disponibili: {len(available_samples)}")
        return None

    meta_tissue = meta_tissue[meta_tissue["sample_id"].isin(available_samples)].copy()

    # Consideriamo solo sample con age_bracket valido
    meta_tissue = meta_tissue[meta_tissue["age_bracket"].isin(AGE_COLUMNS)].copy()

    available_samples = meta_tissue["sample_id"].tolist()

    if len(available_samples) < 10:
        print(f"[SKIP] Troppi pochi sample con age valido: {len(available_samples)}")
        return None

    print(f"Sample validi: {len(available_samples)}")
    print("Distribuzione fasce età:")
    print(meta_tissue["age_bracket"].value_counts().sort_index())

    # Estrai TPM del tessuto: gene x sample
    tissue_tpm = tpm[available_samples].copy()

    # Assicura numerico
    tissue_tpm = tissue_tpm.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Normalizzazione per gene dentro tessuto
    norm_tpm = normalize_minmax_gene_wise(tissue_tpm)

    # Osservato reale
    sample_to_age_real = dict(zip(meta_tissue["sample_id"], meta_tissue["age_bracket"]))

    age_mean_real = compute_age_means(norm_tpm, sample_to_age_real)
    observed = count_switching_genes(age_mean_real, TAU)

    print(f"Observed switching genes: {observed}")

    # Permutazioni sample-level
    true_ages = meta_tissue["age_bracket"].values.copy()
    sample_ids = meta_tissue["sample_id"].values.copy()

    null_counts = []

    for _ in tqdm(range(N_PERMUTATIONS), desc=f"Permutazioni {tissue_name}", leave=False):
        permuted_ages = rng.permutation(true_ages)

        sample_to_age_perm = dict(zip(sample_ids, permuted_ages))

        age_mean_perm = compute_age_means(norm_tpm, sample_to_age_perm)
        null_count = count_switching_genes(age_mean_perm, TAU)

        null_counts.append(null_count)

    null_mean = float(np.mean(null_counts))
    null_std = float(np.std(null_counts, ddof=1)) if len(null_counts) > 1 else 0.0

    if null_std > 0:
        z_score = float((observed - null_mean) / null_std)
    else:
        z_score = np.nan

    p_value = empirical_p_value(observed, null_counts)

    print(f"Null mean: {null_mean:.2f}")
    print(f"Null std: {null_std:.2f}")
    print(f"Z-score: {z_score:.2f}")
    print(f"Empirical p-value: {p_value:.4f}")

    return {
        "version": VERSION,
        "tissue": tissue_name,
        "tau": TAU,
        "n_samples": int(len(available_samples)),
        "n_genes": int(tissue_tpm.shape[0]),
        "n_age_brackets": int(age_mean_real.shape[1]),
        "observed_switching": int(observed),
        "null_mean": null_mean,
        "null_std": null_std,
        "z_score": z_score,
        "empirical_p_value": p_value,
        "n_permutations": N_PERMUTATIONS,
        "min_null": int(np.min(null_counts)),
        "max_null": int(np.max(null_counts)),
    }


# =========================
# MAIN
# =========================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    rng = np.random.default_rng(RANDOM_SEED)

    print("Caricamento TPM...")
    tpm = pd.read_parquet(TPM_PATH)

    print("Caricamento metadata...")
    meta = pd.read_parquet(META_PATH)

    print("Preparazione matrice TPM...")
    tpm = tpm.set_index("gene_id")

    # Rimuove description dalle colonne se presente
    if "description" in tpm.columns:
        tpm = tpm.drop(columns=["description"])

    print(f"TPM finale: {tpm.shape[0]} geni x {tpm.shape[1]} sample")
    print(f"Metadata: {meta.shape[0]} sample")
    if SELECTED_TISSUES is None:
        tissues_to_process = sorted(meta["tissue"].dropna().unique().tolist())
    else:
        tissues_to_process = SELECTED_TISSUES

    print(f"Tessuti da analizzare: {len(tissues_to_process)}")
    print(f"Permutazioni: {N_PERMUTATIONS}")

    results = []

    for tissue_name in tissues_to_process:
        result = process_tissue(tpm, meta, tissue_name, rng)

        if result is not None:
            results.append(result)

    df_results = pd.DataFrame(results)
    df_results.to_csv(OUTPUT_FILE, index=False)

    print("\n=== COMPLETATO ===")
    print(f"File salvato in: {OUTPUT_FILE}")

    if not df_results.empty:
        print(df_results)


if __name__ == "__main__":
    main()