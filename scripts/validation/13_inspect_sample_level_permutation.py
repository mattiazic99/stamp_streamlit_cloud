import os
import pandas as pd

INPUT_FILE = "output/v10/validation/permutation_sample_level/sample_level_permutation_test_all_tissues_100perm.csv"
OUTPUT_DIR = "output/v10/validation/permutation_sample_level"

df = pd.read_csv(INPUT_FILE)

# Rapporto observed / null
df["observed_over_null"] = df["observed_switching"] / df["null_mean"].replace(0, pd.NA)

# Classificazione principale
def classify(row):
    observed = row["observed_switching"]
    null_mean = row["null_mean"]
    p = row["empirical_p_value"]
    z = row["z_score"]

    if observed > null_mean and p <= 0.05:
        return "significant_enriched"
    elif observed > null_mean and p <= 0.10:
        return "borderline_enriched"
    elif observed > null_mean:
        return "non_significant_enriched"
    elif observed < null_mean:
        return "below_null"
    else:
        return "neutral"

df["validation_class"] = df.apply(classify, axis=1)

# Flag per tessuti potenzialmente instabili
df["low_sample_warning"] = df["n_samples"] < 100
df["missing_age_brackets_warning"] = df["n_age_brackets"] < 6

# Ordina per z-score decrescente
df_ranked = df.sort_values("z_score", ascending=False).copy()

print("\n=== RIEPILOGO CLASSI ===")
print(df["validation_class"].value_counts())

print("\n=== TESSUTI SIGNIFICANT ENRICHED ===")
print(
    df[df["validation_class"] == "significant_enriched"]
    .sort_values("z_score", ascending=False)
    [["tissue", "n_samples", "observed_switching", "null_mean", "z_score", "empirical_p_value", "observed_over_null"]]
)

print("\n=== TESSUTI BORDERLINE ENRICHED ===")
print(
    df[df["validation_class"] == "borderline_enriched"]
    .sort_values("empirical_p_value", ascending=True)
    [["tissue", "n_samples", "observed_switching", "null_mean", "z_score", "empirical_p_value", "observed_over_null"]]
)

print("\n=== TESSUTI BELOW NULL ===")
print(
    df[df["validation_class"] == "below_null"]
    .sort_values("z_score", ascending=True)
    [["tissue", "n_samples", "observed_switching", "null_mean", "z_score", "empirical_p_value"]]
)

print("\n=== TESSUTI CON WARNING CAMPIONI < 100 ===")
print(
    df[df["low_sample_warning"]]
    [["tissue", "n_samples", "n_age_brackets", "observed_switching", "null_mean", "z_score", "empirical_p_value", "validation_class"]]
)

print("\n=== TESSUTI CON FASCE ETA' MANCANTI ===")
print(
    df[df["missing_age_brackets_warning"]]
    [["tissue", "n_samples", "n_age_brackets", "observed_switching", "null_mean", "z_score", "empirical_p_value", "validation_class"]]
)

# Salvataggi
ranked_file = os.path.join(OUTPUT_DIR, "sample_level_permutation_ranked_100perm.csv")
summary_file = os.path.join(OUTPUT_DIR, "sample_level_permutation_class_summary_100perm.csv")

df_ranked.to_csv(ranked_file, index=False)

summary = (
    df.groupby("validation_class")
    .agg(
        n_tissues=("tissue", "count"),
        mean_z_score=("z_score", "mean"),
        mean_p_value=("empirical_p_value", "mean"),
        mean_observed=("observed_switching", "mean"),
        mean_null=("null_mean", "mean"),
    )
    .reset_index()
)

summary.to_csv(summary_file, index=False)

print(f"\nFile ranked salvato in: {ranked_file}")
print(f"File summary salvato in: {summary_file}")