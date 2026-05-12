#!/usr/bin/env python3
"""
Multiple-comparisons correction for the Khan parameter-space sweep.
===================================================================

Loads khan_summary.csv (132 configurations) and applies Bonferroni,
Holm-Bonferroni (step-down), and Benjamini-Hochberg (FDR) corrections
to the paired t-test p-values.

Output
------
results/khan_multiple_comparisons.csv
    One row per configuration with original and corrected p-values.
"""

from pathlib import Path

import numpy as np
import pandas as pd

RESULTS_DIR = Path(__file__).resolve().parent.parent.parent / "build" / "results"
INPUT_CSV = RESULTS_DIR / "khan_summary.csv"
OUTPUT_CSV = RESULTS_DIR / "khan_multiple_comparisons.csv"


def bonferroni(pvals: np.ndarray) -> np.ndarray:
    """Bonferroni correction: p_adj = min(p * m, 1)."""
    return np.minimum(pvals * len(pvals), 1.0)


def holm(pvals: np.ndarray) -> np.ndarray:
    """Holm-Bonferroni (step-down) correction."""
    m = len(pvals)
    order = np.argsort(pvals)
    adjusted = np.ones(m)
    for rank, idx in enumerate(order):
        adjusted[idx] = min(pvals[idx] * (m - rank), 1.0)
    # enforce monotonicity
    for rank in range(1, m):
        idx = order[rank]
        prev = order[rank - 1]
        adjusted[idx] = max(adjusted[idx], adjusted[prev])
    return adjusted


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg (FDR) correction."""
    m = len(pvals)
    order = np.argsort(pvals)
    adjusted = np.ones(m)
    for rank, idx in enumerate(order):
        adjusted[idx] = min(pvals[idx] * m / (rank + 1), 1.0)
    # enforce monotonicity (reverse direction)
    for rank in range(m - 2, -1, -1):
        idx = order[rank]
        nxt = order[rank + 1]
        adjusted[idx] = min(adjusted[idx], adjusted[nxt])
    return adjusted


def main() -> None:
    df = pd.read_csv(INPUT_CSV)
    m = len(df)

    pvals = df["p_value_t"].values
    cohen_d = df["cohen_d"].values

    # Apply corrections
    df["bonferroni_p"] = bonferroni(pvals)
    df["holm_p"] = holm(pvals)
    df["bh_p"] = benjamini_hochberg(pvals)

    alpha = 0.05
    df["bonferroni_sig"] = df["bonferroni_p"] < alpha
    df["holm_sig"] = df["holm_p"] < alpha
    df["bh_sig"] = df["bh_p"] < alpha

    # Verdict: combine significance with direction
    df["uncorrected_better"] = (pvals < alpha) & (cohen_d > 0)
    df["uncorrected_worse"] = (pvals < alpha) & (cohen_d < 0)
    df["bonferroni_better"] = df["bonferroni_sig"] & (cohen_d > 0)
    df["bonferroni_worse"] = df["bonferroni_sig"] & (cohen_d < 0)

    df.to_csv(OUTPUT_CSV, index=False)

    # Summary
    n_uncorrected = (pvals < alpha).sum()
    n_better = df["uncorrected_better"].sum()
    n_worse = df["uncorrected_worse"].sum()
    n_bonf = df["bonferroni_sig"].sum()
    n_bonf_better = df["bonferroni_better"].sum()
    n_bonf_worse = df["bonferroni_worse"].sum()
    n_holm = df["holm_sig"].sum()
    n_bh = df["bh_sig"].sum()

    print("=" * 60)
    print("Multiple-Comparisons Correction (Khan Parameter Space)")
    print("=" * 60)
    print(f"Total configurations: {m}")
    print(f"\nSignificant results (alpha = {alpha}):")
    print(f"  Uncorrected:        {n_uncorrected:3d}/{m} "
          f"({n_better} better, {n_worse} worse)")
    print(f"  Bonferroni:         {n_bonf:3d}/{m} "
          f"({n_bonf_better} better, {n_bonf_worse} worse)")
    print(f"  Holm-Bonferroni:    {n_holm:3d}/{m}")
    print(f"  Benjamini-Hochberg: {n_bh:3d}/{m}")
    print(f"\n  Dropped by Bonferroni: {n_uncorrected - n_bonf}")
    print(f"\nOutput: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
