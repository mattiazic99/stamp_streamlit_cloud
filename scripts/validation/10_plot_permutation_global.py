import os
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# CONFIGURAZIONE
# =========================

INPUT_FILE = "output/v10/validation/permutation_global/permutation_summary_simple_ranked.csv"
OUTPUT_DIR = "output/v10/validation/permutation_global"

OUTPUT_ZSCORE_PNG = os.path.join(OUTPUT_DIR, "permutation_zscore_by_tissue.png")
OUTPUT_OBS_NULL_PNG = os.path.join(OUTPUT_DIR, "permutation_observed_vs_null.png")


# =========================
# LOAD
# =========================

df = pd.read_csv(INPUT_FILE)

# Se manca effect_type, lo ricostruiamo
if "effect_type" not in df.columns:
    def classify_effect(row):
        if row["z_score"] >= 2:
            return "enriched"
        elif row["z_score"] <= -2:
            return "depleted"
        else:
            return "neutral"

    df["effect_type"] = df.apply(classify_effect, axis=1)

# Ordine per z_score crescente
df_sorted = df.sort_values("z_score", ascending=True).copy()


# =========================
# FIGURA 1: Z-SCORE - PAPER STYLE
# =========================

df_sorted = df.sort_values("z_score", ascending=True).copy()

plt.figure(figsize=(12, 11))

colors = df_sorted["effect_type"].map({
    "enriched": "#1f77b4",
    "depleted": "#d62728",
    "neutral": "#7f7f7f"
})

plt.barh(df_sorted["tissue"], df_sorted["z_score"], color=colors)

plt.axvline(0, color="black", linestyle="-", linewidth=0.8)
plt.axvline(2, color="gray", linestyle="--", linewidth=0.8)
plt.axvline(-2, color="gray", linestyle="--", linewidth=0.8)

plt.xlabel("Z-score relative to the null distribution", fontsize=12)
plt.ylabel("Tissue", fontsize=12)
plt.title("Permutation-based validation of STAMP switching gene counts", fontsize=14)

plt.xticks(fontsize=10)
plt.yticks(fontsize=8)

plt.tight_layout()

output_paper_png = os.path.join(OUTPUT_DIR, "permutation_zscore_by_tissue_paper.png")
output_paper_pdf = os.path.join(OUTPUT_DIR, "permutation_zscore_by_tissue_paper.pdf")

plt.savefig(output_paper_png, dpi=300, bbox_inches="tight")
plt.savefig(output_paper_pdf, bbox_inches="tight")
plt.close()

print(f"Salvata figura paper PNG: {output_paper_png}")
print(f"Salvata figura paper PDF: {output_paper_pdf}")


# =========================
# FIGURA 2: OBSERVED VS NULL
# =========================

df_obs = df.sort_values("observed_switching", ascending=False).copy()

x = range(len(df_obs))

plt.figure(figsize=(16, 7))

plt.plot(x, df_obs["observed_switching"], marker="o", label="Observed switching genes")
plt.plot(x, df_obs["null_mean"], marker="o", label="Mean null switching genes")

plt.xticks(x, df_obs["tissue"], rotation=90)
plt.ylabel("Numero di geni switching")
plt.title("Observed vs null switching gene counts")
plt.legend()
plt.tight_layout()

plt.savefig(OUTPUT_OBS_NULL_PNG, dpi=300)
plt.close()

print(f"Salvata figura observed vs null: {OUTPUT_OBS_NULL_PNG}")


# =========================
# RIEPILOGO
# =========================

print("\n=== RIEPILOGO ===")
print(df["effect_type"].value_counts())

print("\nFigure generate correttamente.")