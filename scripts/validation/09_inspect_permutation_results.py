import os
import pandas as pd

INPUT_FILE = "output/v10/validation/permutation_global/permutation_summary_simple.csv"
OUTPUT_DIR = "output/v10/validation/permutation_global"

df = pd.read_csv(INPUT_FILE)

print("\n=== INFO GENERALI ===")
print(f"Tessuti analizzati: {df.shape[0]}")
print(f"Permutazioni: {df['n_permutations'].iloc[0]}")
print(f"Soglia tau: {df['tau'].iloc[0]}")

print("\n=== P-VALUE ===")
print(df["empirical_p_value"].describe())

print("\n=== TOP 10 tessuti per observed_switching ===")
print(
    df.sort_values("observed_switching", ascending=False)
    [["tissue", "observed_switching", "null_mean", "z_score", "empirical_p_value"]]
    .head(10)
)

print("\n=== TOP 10 tessuti per z_score ===")
print(
    df.sort_values("z_score", ascending=False)
    [["tissue", "observed_switching", "null_mean", "z_score", "empirical_p_value"]]
    .head(10)
)

print("\n=== Tessuti con z_score più basso ===")
print(
    df.sort_values("z_score", ascending=True)
    [["tissue", "observed_switching", "null_mean", "z_score", "empirical_p_value"]]
    .head(10)
)

# Rapporto observed / null_mean
df["observed_over_null"] = df["observed_switching"] / df["null_mean"]

def classify_effect(row):
    if row["z_score"] >= 2:
        return "enriched"
    elif row["z_score"] <= -2:
        return "depleted"
    else:
        return "neutral"

df["effect_type"] = df.apply(classify_effect, axis=1)

print("\n=== TOP 10 observed/null ===")
print(
    df.sort_values("observed_over_null", ascending=False)
    [["tissue", "observed_switching", "null_mean", "observed_over_null", "empirical_p_value"]]
    .head(10)
)

print("\n=== CONTEGGIO EFFECT TYPE ===")
print(df["effect_type"].value_counts())

print("\n=== TESSUTI ENRICHED ===")
print(
    df[df["effect_type"] == "enriched"]
    .sort_values("z_score", ascending=False)
    [["tissue", "observed_switching", "null_mean", "z_score", "empirical_p_value"]]
)

print("\n=== TESSUTI DEPLETED ===")
print(
    df[df["effect_type"] == "depleted"]
    .sort_values("z_score", ascending=True)
    [["tissue", "observed_switching", "null_mean", "z_score", "empirical_p_value"]]
)

output_file = os.path.join(OUTPUT_DIR, "permutation_summary_simple_ranked.csv")
df.sort_values("z_score", ascending=False).to_csv(output_file, index=False)

print(f"\nFile ordinato salvato in: {output_file}")