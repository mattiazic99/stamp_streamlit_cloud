import os
import pandas as pd
import matplotlib.pyplot as plt


INPUT_FILE = "output/v10/validation/permutation_sample_level/sample_level_permutation_ranked_100perm.csv"
OUTPUT_DIR = "output/v10/validation/permutation_sample_level"

OUTPUT_PNG = os.path.join(OUTPUT_DIR, "sample_level_permutation_zscore_by_tissue.png")
OUTPUT_PDF = os.path.join(OUTPUT_DIR, "sample_level_permutation_zscore_by_tissue.pdf")

df = pd.read_csv(INPUT_FILE)

# Se manca validation_class, la ricostruiamo
if "validation_class" not in df.columns:
    def classify(row):
        observed = row["observed_switching"]
        null_mean = row["null_mean"]
        p = row["empirical_p_value"]

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

df_sorted = df.sort_values("z_score", ascending=True).copy()

color_map = {
    "significant_enriched": "#1f77b4",
    "borderline_enriched": "#ff7f0e",
    "non_significant_enriched": "#7f7f7f",
    "below_null": "#d62728",
    "neutral": "#bdbdbd",
}

colors = df_sorted["validation_class"].map(color_map)

plt.figure(figsize=(12, 11))

plt.barh(df_sorted["tissue"], df_sorted["z_score"], color=colors)

plt.axvline(0, color="black", linestyle="-", linewidth=0.8)
plt.axvline(2, color="gray", linestyle="--", linewidth=0.8)
plt.axvline(-2, color="gray", linestyle="--", linewidth=0.8)

plt.xlabel("Z-score relative to the sample-level null distribution", fontsize=12)
plt.ylabel("Tissue", fontsize=12)
plt.title("Sample-level permutation validation of STAMP switching gene counts", fontsize=14)

plt.xticks(fontsize=10)
plt.yticks(fontsize=8)

# Legenda manuale
handles = [
    plt.Rectangle((0, 0), 1, 1, color=color_map["significant_enriched"], label="Significant enriched"),
    plt.Rectangle((0, 0), 1, 1, color=color_map["borderline_enriched"], label="Borderline enriched"),
    plt.Rectangle((0, 0), 1, 1, color=color_map["non_significant_enriched"], label="Non-significant enriched"),
    plt.Rectangle((0, 0), 1, 1, color=color_map["below_null"], label="Below null"),
]

plt.legend(handles=handles, loc="lower right", fontsize=9)

plt.tight_layout()

plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
plt.savefig(OUTPUT_PDF, bbox_inches="tight")

plt.close()

print(f"Figura PNG salvata in: {OUTPUT_PNG}")
print(f"Figura PDF salvata in: {OUTPUT_PDF}")

print("\nClassi:")
print(df["validation_class"].value_counts())