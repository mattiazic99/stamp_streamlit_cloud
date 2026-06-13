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

NORMALIZED_DIR = f"output/{VERSION}/normalized"
OUTPUT_DIR = f"output/{VERSION}/validation/permutation_global"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "permutation_summary_simple.csv")

AGE_COLUMNS = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]


# =========================
# FUNZIONI
# =========================

def is_switching_binary_vector(values):
    """
    Ritorna True se il vettore binario ha esattamente una transizione.
    Esempi switching:
    [0, 0, 1, 1, 1, 1]
    [1, 1, 1, 0, 0, 0]

    Esempi non switching:
    [0, 0, 0, 0, 0, 0]
    [1, 1, 1, 1, 1, 1]
    [0, 1, 0, 1, 1, 1]
    """
    changes = 0

    for i in range(len(values) - 1):
        if values[i] != values[i + 1]:
            changes += 1

    return changes == 1


def count_switching_genes(df, tau):
    """
    Conta quanti geni sono switching in una matrice gene x age_bracket.
    """
    binary_df = (df >= tau).astype(int)

    switching_mask = binary_df.apply(
        lambda row: is_switching_binary_vector(row.values.tolist()),
        axis=1
    )

    return int(switching_mask.sum())


def empirical_p_value(observed, null_values):
    """
    P-value empirico one-sided:
    quante volte il nullo produce un numero di switching >= osservato.
    Con correzione +1 per evitare p-value zero.
    """
    null_values = np.array(null_values)

    return (np.sum(null_values >= observed) + 1) / (len(null_values) + 1)


def process_tissue(file_path, tissue_name, rng):
    """
    Esegue il permutation test per un singolo tessuto.
    """
    df = pd.read_parquet(file_path)

    # Teniamo solo le colonne di età presenti
    available_age_cols = [col for col in AGE_COLUMNS if col in df.columns]

    if len(available_age_cols) < 2:
        print(f"[SKIP] {tissue_name}: meno di 2 fasce d'età disponibili")
        return None

    df = df[available_age_cols]

    # Rimuove geni con valori mancanti
    df = df.dropna()

    if df.empty:
        print(f"[SKIP] {tissue_name}: dataframe vuoto dopo dropna")
        return None

    observed = count_switching_genes(df, TAU)

    null_counts = []

    for _ in range(N_PERMUTATIONS):
        permuted_df = df.copy()

        # Permutazione semplice: mischiamo le colonne di età per ogni gene
        permuted_values = np.apply_along_axis(
            lambda x: rng.permutation(x),
            axis=1,
            arr=permuted_df.values
        )

        permuted_df = pd.DataFrame(
            permuted_values,
            index=df.index,
            columns=df.columns
        )

        null_count = count_switching_genes(permuted_df, TAU)
        null_counts.append(null_count)

    null_mean = float(np.mean(null_counts))
    null_std = float(np.std(null_counts, ddof=1)) if len(null_counts) > 1 else 0.0

    if null_std > 0:
        z_score = float((observed - null_mean) / null_std)
    else:
        z_score = np.nan

    p_value = empirical_p_value(observed, null_counts)

    return {
        "version": VERSION,
        "tissue": tissue_name,
        "tau": TAU,
        "n_genes": int(df.shape[0]),
        "n_age_brackets": int(df.shape[1]),
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

    files = sorted([
        f for f in os.listdir(NORMALIZED_DIR)
        if f.endswith(".parquet")
    ])

    print(f"Versione: {VERSION}")
    print(f"Directory input: {NORMALIZED_DIR}")
    print(f"File trovati: {len(files)}")
    print(f"Permutazioni per tessuto: {N_PERMUTATIONS}")
    print(f"Soglia tau: {TAU}")
    print()

    results = []

    for filename in tqdm(files, desc="Tessuti"):
        file_path = os.path.join(NORMALIZED_DIR, filename)
        tissue_name = filename.replace(".parquet", "")

        result = process_tissue(file_path, tissue_name, rng)

        if result is not None:
            results.append(result)

    df_results = pd.DataFrame(results)
    df_results.to_csv(OUTPUT_FILE, index=False)

    print()
    print("Validazione completata.")
    print(f"File salvato in: {OUTPUT_FILE}")

    if not df_results.empty:
        print()
        print("Prime righe:")
        print(df_results.head())


if __name__ == "__main__":
    main()