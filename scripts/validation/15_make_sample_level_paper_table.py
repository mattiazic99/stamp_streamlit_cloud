import os
import pandas as pd


# =========================
# CONFIGURAZIONE
# =========================

INPUT_FILE = "output/v10/validation/permutation_sample_level/sample_level_permutation_ranked_100perm.csv"
OUTPUT_DIR = "output/v10/validation/permutation_sample_level"

OUTPUT_FULL = os.path.join(OUTPUT_DIR, "paper_table_sample_level_permutation_full.csv")
OUTPUT_MAIN = os.path.join(OUTPUT_DIR, "paper_table_sample_level_permutation_main.csv")
OUTPUT_LATEX = os.path.join(OUTPUT_DIR, "paper_table_sample_level_permutation_main.tex")


# =========================
# LOAD
# =========================

df = pd.read_csv(INPUT_FILE)


# =========================
# RICOSTRUZIONE CAMPI UTILI
# =========================

if "observed_over_null" not in df.columns:
    df["observed_over_null"] = df["observed_switching"] / df["null_mean"].replace(0, pd.NA)

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

df["low_sample_warning"] = df["n_samples"] < 100
df["missing_age_brackets_warning"] = df["n_age_brackets"] < 6

def make_warning(row):
    warnings = []

    if row["low_sample_warning"]:
        warnings.append("low sample size")

    if row["missing_age_brackets_warning"]:
        warnings.append("missing age bins")

    if not warnings:
        return ""

    return "; ".join(warnings)

df["warning"] = df.apply(make_warning, axis=1)


# =========================
# TABELLA FULL
# =========================

full_cols = [
    "tissue",
    "n_samples",
    "n_age_brackets",
    "observed_switching",
    "null_mean",
    "null_std",
    "z_score",
    "empirical_p_value",
    "observed_over_null",
    "validation_class",
    "warning",
]

df_full = df[full_cols].copy()

# Arrotondamenti leggibili
df_full["null_mean"] = df_full["null_mean"].round(2)
df_full["null_std"] = df_full["null_std"].round(2)
df_full["z_score"] = df_full["z_score"].round(3)
df_full["empirical_p_value"] = df_full["empirical_p_value"].round(4)
df_full["observed_over_null"] = df_full["observed_over_null"].round(3)

df_full = df_full.sort_values("z_score", ascending=False)

df_full.to_csv(OUTPUT_FULL, index=False)


# =========================
# TABELLA MAIN PAPER
# =========================
# Per la tabella principale teniamo:
# - tutti i significant enriched
# - tutti i borderline
# - tutti i below_null
# I non_significant_enriched possono andare in supplementare.

main_classes = [
    "significant_enriched",
    "borderline_enriched",
    "below_null",
]

df_main = df_full[df_full["validation_class"].isin(main_classes)].copy()

df_main.to_csv(OUTPUT_MAIN, index=False)


# =========================
# EXPORT LATEX
# =========================

latex_df = df_main.copy()

latex_df = latex_df.rename(columns={
    "tissue": "Tissue",
    "n_samples": "N",
    "n_age_brackets": "Age bins",
    "observed_switching": "Observed",
    "null_mean": "Null mean",
    "z_score": "Z-score",
    "empirical_p_value": "Empirical p",
    "validation_class": "Class",
    "warning": "Warning",
})

latex_df = latex_df[
    [
        "Tissue",
        "N",
        "Age bins",
        "Observed",
        "Null mean",
        "Z-score",
        "Empirical p",
        "Class",
        "Warning",
    ]
]

latex_code = latex_df.to_latex(
    index=False,
    escape=True,
    longtable=True,
    caption="Sample-level age-label permutation validation of STAMP switching gene counts across GTEx v10 tissues.",
    label="tab:sample_level_permutation_validation"
)

with open(OUTPUT_LATEX, "w", encoding="utf-8") as f:
    f.write(latex_code)


# =========================
# SUMMARY
# =========================

print("\n=== CLASS SUMMARY ===")
print(df["validation_class"].value_counts())

print("\n=== FILE GENERATI ===")
print(OUTPUT_FULL)
print(OUTPUT_MAIN)
print(OUTPUT_LATEX)

print("\n=== PREVIEW MAIN TABLE ===")
print(latex_df.head(20))