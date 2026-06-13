import pandas as pd

TPM_PATH = "data/parquet/v10/tpm_matrix.parquet"
META_PATH = "data/parquet/v10/metadata.parquet"

print("Caricamento TPM...")
tpm = pd.read_parquet(TPM_PATH)
print("TPM shape:", tpm.shape)
print("TPM columns preview:")
print(tpm.columns[:10].tolist())
print(tpm.head())

print("\nCaricamento metadata...")
meta = pd.read_parquet(META_PATH)
print("Metadata shape:", meta.shape)
print("Metadata columns:")
print(meta.columns.tolist())
print(meta.head())