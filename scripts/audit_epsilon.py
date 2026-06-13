"""Audit del filtro epsilon: quanti geni vengono scartati per ogni tessuto.

Per ogni tessuto, calcola il range (max - min) di ogni gene tra i suoi
sample e suddivide i geni in 4 categorie:

  - all-NaN:        gene completamente mancante per quel tessuto
                    (nessun sample del tessuto ha un valore valido)
  - costanti:       range == 0 (tutti i sample hanno lo stesso valore,
                    tipicamente 0 -> gene "spento" nel tessuto)
  - sotto epsilon:  0 < range <= epsilon
                    (qui sta il rumore che la normalizzazione amplificherebbe)
  - tenuti:         range > epsilon (passano il filtro)

Salva un CSV con i conteggi per tessuto e stampa una tabella riassuntiva
ordinata per % di geni scartati.

Opzionale: con --sweep prova diversi valori di epsilon (0.001, 0.005, 0.01,
0.05, 0.1) e produce una tabella con i conteggi per tutti questi valori,
così vedi quanto è sensibile la scelta del cutoff.

Usage
-----
    python scripts/audit_epsilon.py --version v10
    python scripts/audit_epsilon.py --version v10 --epsilon 0.05
    python scripts/audit_epsilon.py --version v10 --sweep
    python scripts/audit_epsilon.py --version v10 --tissue "Liver"
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from stamp.config import EPSILON, GtexVersion, paths_for
from stamp.io import load_metadata, load_tpm_matrix


SWEEP_EPSILONS = [0.001, 0.005, 0.01, 0.05, 0.1]


# ---------------------------------------------------------------------------
# Core: calcolo del range per ogni gene di un tessuto
# ---------------------------------------------------------------------------

def compute_gene_ranges(
    tpm: pd.DataFrame,
    metadata: pd.DataFrame,
    tissue: str,
) -> pd.Series | None:
    """Restituisce la Serie dei range (max - min) per ogni gene del tessuto.

    NaN nei TPM sono trattati come mancanti (skipna=True). Se nessun sample
    del tessuto e' presente nella matrice, restituisce None.
    """
    samples = metadata.loc[metadata["tissue"] == tissue, "sample_id"].tolist()
    samples = [s for s in samples if s in tpm.columns]
    if not samples:
        return None

    sub = tpm[samples]
    gene_min = sub.min(axis=1, skipna=True)
    gene_max = sub.max(axis=1, skipna=True)
    return gene_max - gene_min


def bucketize(gene_range: pd.Series, epsilon: float) -> dict[str, int]:
    """Conta i geni nelle 4 categorie definite dal filtro epsilon."""
    all_nan = gene_range.isna()
    constant = (gene_range == 0) & ~all_nan
    below_eps = (gene_range > 0) & (gene_range <= epsilon)
    kept = gene_range > epsilon

    return {
        "n_genes_total": int(len(gene_range)),
        "n_all_nan": int(all_nan.sum()),
        "n_constant": int(constant.sum()),
        "n_below_epsilon": int(below_eps.sum()),
        "n_kept": int(kept.sum()),
    }


# ---------------------------------------------------------------------------
# Modalita' 1: report standard per un singolo epsilon
# ---------------------------------------------------------------------------

def run_single_audit(
    tpm: pd.DataFrame,
    metadata: pd.DataFrame,
    tissues: list[str],
    epsilon: float,
) -> pd.DataFrame:
    rows = []
    for tissue in tqdm(tissues, desc="Tessuti", unit="tissue"):
        gene_range = compute_gene_ranges(tpm, metadata, tissue)
        if gene_range is None:
            continue

        counts = bucketize(gene_range, epsilon)
        n_samples = (metadata["tissue"] == tissue).sum()

        discarded = counts["n_all_nan"] + counts["n_constant"] + counts["n_below_epsilon"]
        pct_discarded = 100.0 * discarded / counts["n_genes_total"]

        # Statistiche sui range (utili per scegliere epsilon)
        valid = gene_range.dropna()
        range_stats = {
            "range_p05": float(np.percentile(valid, 5)) if len(valid) else float("nan"),
            "range_p50": float(np.percentile(valid, 50)) if len(valid) else float("nan"),
            "range_p95": float(np.percentile(valid, 95)) if len(valid) else float("nan"),
        }

        rows.append({
            "tissue": tissue,
            "n_samples": int(n_samples),
            **counts,
            "pct_discarded": round(pct_discarded, 2),
            **range_stats,
        })

    return pd.DataFrame(rows).sort_values("pct_discarded", ascending=False)


def print_single_report(df: pd.DataFrame, epsilon: float) -> None:
    print(f"\nReport per tessuto (epsilon={epsilon}, ordinato per % scartata):\n")
    header = (
        f"{'Tessuto':<45s} {'tot':>6s} {'NaN':>6s} "
        f"{'cost':>6s} {'<eps':>6s} {'tenuti':>7s} {'%scart':>8s}"
    )
    print(header)
    print("-" * len(header))
    for _, r in df.iterrows():
        print(
            f"{r['tissue']:<45s} "
            f"{r['n_genes_total']:>6d} "
            f"{r['n_all_nan']:>6d} "
            f"{r['n_constant']:>6d} "
            f"{r['n_below_epsilon']:>6d} "
            f"{r['n_kept']:>7d} "
            f"{r['pct_discarded']:>7.2f}%"
        )

    print("\nMedie su tutti i tessuti:")
    print(f"  Geni totali medi:       {df['n_genes_total'].mean():>8.0f}")
    print(f"  All-NaN medi:           {df['n_all_nan'].mean():>8.0f}")
    print(f"  Costanti medi:          {df['n_constant'].mean():>8.0f}")
    print(f"  Sotto epsilon medi:     {df['n_below_epsilon'].mean():>8.0f}")
    print(f"  Tenuti medi:            {df['n_kept'].mean():>8.0f}")
    print(f"  % scartata media:       {df['pct_discarded'].mean():>7.2f}%")


# ---------------------------------------------------------------------------
# Modalita' 2: sweep su diversi valori di epsilon
# ---------------------------------------------------------------------------

def run_sweep(
    tpm: pd.DataFrame,
    metadata: pd.DataFrame,
    tissues: list[str],
    epsilons: list[float],
) -> pd.DataFrame:
    """Per ogni tessuto, calcola i conteggi per piu' valori di epsilon."""
    rows = []
    for tissue in tqdm(tissues, desc="Tessuti", unit="tissue"):
        gene_range = compute_gene_ranges(tpm, metadata, tissue)
        if gene_range is None:
            continue

        n_total = int(len(gene_range))
        row = {"tissue": tissue, "n_genes_total": n_total}
        for eps in epsilons:
            counts = bucketize(gene_range, eps)
            row[f"kept_eps_{eps}"] = counts["n_kept"]
            row[f"pct_kept_eps_{eps}"] = round(100.0 * counts["n_kept"] / n_total, 2)
        rows.append(row)

    return pd.DataFrame(rows).sort_values(f"pct_kept_eps_{epsilons[-1]}")


def print_sweep_report(df: pd.DataFrame, epsilons: list[float]) -> None:
    print(f"\nSweep epsilon: numero di geni TENUTI per ciascun valore di epsilon\n")
    header_parts = [f"{'Tessuto':<45s}", f"{'tot':>6s}"]
    for eps in epsilons:
        header_parts.append(f"eps={eps:<6}")
    header = " ".join(header_parts)
    print(header)
    print("-" * len(header))

    for _, r in df.iterrows():
        parts = [f"{r['tissue']:<45s}", f"{r['n_genes_total']:>6d}"]
        for eps in epsilons:
            parts.append(f"{r[f'kept_eps_{eps}']:>10d}")
        print(" ".join(parts))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    parser.add_argument("--epsilon", type=float, default=EPSILON,
                        help=f"Valore di epsilon (default: {EPSILON})")
    parser.add_argument("--tissue", default=None,
                        help="Audita solo questo tessuto.")
    parser.add_argument("--sweep", action="store_true",
                        help="Prova piu' valori di epsilon: "
                             f"{SWEEP_EPSILONS}")
    parser.add_argument("--output", type=Path, default=None,
                        help="Path del CSV (default: output/{version}/epsilon_audit.csv)")
    args = parser.parse_args()

    print(f"=== Epsilon audit (GTEx {args.version}) ===")
    print("  Carico metadata...")
    metadata = load_metadata(args.version)

    all_tissues = sorted(metadata["tissue"].dropna().unique().tolist())
    if args.tissue:
        if args.tissue not in all_tissues:
            raise SystemExit(
                f"Tessuto '{args.tissue}' non trovato. "
                f"Disponibili: {all_tissues}"
            )
        tissues = [args.tissue]
    else:
        tissues = all_tissues

    print(f"  Carico matrice TPM...")
    tpm = load_tpm_matrix(args.version)
    print(f"  Shape: {tpm.shape}")

    if args.sweep:
        df = run_sweep(tpm, metadata, tissues, SWEEP_EPSILONS)
        print_sweep_report(df, SWEEP_EPSILONS)
        default_name = "epsilon_audit_sweep.csv"
    else:
        df = run_single_audit(tpm, metadata, tissues, args.epsilon)
        print_single_report(df, args.epsilon)
        default_name = "epsilon_audit.csv"

    out_path = args.output or Path(f"output/{args.version}/{default_name}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"\nCSV salvato: {out_path}")


if __name__ == "__main__":
    main()