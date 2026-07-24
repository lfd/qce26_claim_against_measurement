#!/usr/bin/env python3
"""
Drift Power Analysis — Intraclass Correlation and Autocorrelation
============================================================================

Quantifies how strongly temporal drift structures a longitudinal ZNE
experiment.

Method
------
1. Intraclass Correlation Coefficient (ICC type-1):
       ICC = sigma_between^2 / (sigma_between^2 + sigma_within^2)
   where *between* is TP-to-TP variance of means, *within* is the pooled
   per-TP variance.  High ICC means drift dominates — repetitions within a
   TP are more alike than repetitions across TPs.

2. Lag-k autocorrelation of the per-TP means, and the Durbin-Watson
   statistic on the linear-detrended residuals.

Inputs
------
  --csv  Path to a raw_data_*.csv file (default: qexa_drift/raw_data_week.csv)
  --outdir  Output directory (default: build/results)

Outputs
-------
  drift_power.csv  — per-TP stats plus ICC summary
  drift_power_summary.txt — plain-text report

Usage
-----
  cd reproduction
  python scripts/power_analysis.py
  python scripts/power_analysis.py --csv data/qexa_drift/raw_data_weekend.csv

Note: expects the qexa_drift CSV schema (timepoint_idx, timestamp, backend, rep,
scale_factor, exp_val, n_shots, ideal). Not compatible with
data/ibm_drift_results/raw_data_ibm_drift.csv, which uses a different schema.
"""

from __future__ import annotations

import argparse
import csv as csv_mod
import math
import os
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
BUILD_DIR = REPO_DIR.parent / "build"
DEFAULT_CSV = REPO_DIR / "data" / "qexa_drift" / "raw_data_week.csv"
DEFAULT_OUTDIR = BUILD_DIR / "results"


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

def load_csv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv_mod.DictReader(f))


def autocorr(x: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Return autocorrelation coefficients rho[1], ..., rho[max_lag]."""
    n = len(x)
    if max_lag is None:
        max_lag = min(n - 1, int(2 * n ** (2 / 3)))
    mu = x.mean()
    var = ((x - mu) ** 2).mean()
    if var < 1e-15:
        return np.zeros(max_lag)
    rhos = []
    for k in range(1, max_lag + 1):
        cov = ((x[:-k] - mu) * (x[k:] - mu)).mean()
        rhos.append(cov / var)
    return np.array(rhos)


def icc_one_way(groups: list[np.ndarray]) -> float:
    """ICC (one-way random effects, type 1) for a list of groups."""
    k = len(groups)
    ns = np.array([len(g) for g in groups])
    n_total = ns.sum()
    grand_mean = sum(g.sum() for g in groups) / n_total

    # Between-group sum of squares
    ss_between = sum(n * (g.mean() - grand_mean) ** 2
                     for n, g in zip(ns, groups))
    # Within-group sum of squares
    ss_within = sum(((g - g.mean()) ** 2).sum() for g in groups)

    df_between = k - 1
    df_within = n_total - k

    ms_between = ss_between / df_between if df_between > 0 else 0.0
    ms_within = ss_within / df_within if df_within > 0 else 1.0

    n0 = (n_total - (ns ** 2).sum() / n_total) / (k - 1)  # harmonic-ish n0
    sigma_between = max((ms_between - ms_within) / n0, 0.0)
    sigma_within = ms_within

    denom = sigma_between + sigma_within
    return sigma_between / denom if denom > 0 else 0.0


def durbin_watson(resid: np.ndarray) -> float:
    """Durbin-Watson statistic (2 = no autocorrelation)."""
    diff = np.diff(resid)
    return float(diff @ diff / (resid @ resid)) if resid @ resid > 0 else 2.0


# ═══════════════════════════════════════════════════════════════════════
# Main analysis
# ═══════════════════════════════════════════════════════════════════════

def analyse(csv_path: Path, outdir: Path) -> None:
    rows = load_csv(csv_path)
    ideal = float(rows[0]["ideal"])

    # Build per-TP groups of E(λ=1) values
    tp_map: dict[int, list[float]] = {}
    scale_factor_col = "scale_factor"

    for r in rows:
        if abs(float(r[scale_factor_col]) - 1.0) > 0.1:
            continue
        tp = int(r["timepoint_idx"])
        tp_map.setdefault(tp, []).append(float(r["exp_val"]))

    tps = sorted(tp_map)
    groups = [np.array(tp_map[tp]) for tp in tps]
    tp_means = np.array([g.mean() for g in groups])
    tp_ses = np.array([g.std() / math.sqrt(len(g)) for g in groups])
    n_reps = np.array([len(g) for g in groups])

    # ── ICC ────────────────────────────────────────────────────────
    icc = icc_one_way(groups)

    # ── Autocorrelation ────────────────────────────────────────────
    rhos = autocorr(tp_means, max_lag=min(20, len(tp_means) // 4))

    # ── Durbin-Watson on mean residuals ───────────────────────────
    trend = np.polyfit(np.arange(len(tp_means)), tp_means, 1)
    resid = tp_means - np.polyval(trend, np.arange(len(tp_means)))
    dw = durbin_watson(resid)

    # ── ZNE Cohen's d per TP ──────────────────────────────────────
    cohen_ds: list[float] = []
    for tp, group in zip(tps, groups):
        tp_rows = [r for r in rows if int(r["timepoint_idx"]) == tp]
        e1 = np.array([float(r["exp_val"]) for r in tp_rows
                       if abs(float(r[scale_factor_col]) - 1.0) < 0.1])
        e3 = np.array([float(r["exp_val"]) for r in tp_rows
                       if abs(float(r[scale_factor_col]) - 3.0) < 0.1])
        e5 = np.array([float(r["exp_val"]) for r in tp_rows
                       if abs(float(r[scale_factor_col]) - 5.0) < 0.1])
        if len(e1) == len(e3) == len(e5) > 0:
            zne = (15 / 8) * e1 - (5 / 4) * e3 + (3 / 8) * e5
            raw_err = np.abs(e1 - ideal)
            mit_err = np.abs(zne - ideal)
            imp = raw_err - mit_err
            std_imp = imp.std()
            d = float(imp.mean() / std_imp) if std_imp > 1e-12 else 0.0
            cohen_ds.append(d)

    # ── Write per-TP CSV ──────────────────────────────────────────
    outdir.mkdir(parents=True, exist_ok=True)
    per_tp_path = outdir / "drift_power.csv"
    with open(per_tp_path, "w", newline="") as f:
        writer = csv_mod.writer(f)
        writer.writerow([
            "timepoint_idx", "n_reps", "mean_e_l1", "se_e_l1",
            "cohen_d_zne", "resid_mean",
        ])
        for i, (tp, g, mean_e, se_e) in enumerate(
                zip(tps, groups, tp_means, tp_ses)):
            d = cohen_ds[i] if i < len(cohen_ds) else float("nan")
            writer.writerow([tp, len(g), f"{mean_e:.6f}", f"{se_e:.6f}",
                             f"{d:.4f}", f"{resid[i]:.6f}"])

    # ── Write summary text ────────────────────────────────────────
    summary_path = outdir / "drift_power_summary.txt"
    N_total = sum(len(g) for g in groups)
    n_tps = len(tps)
    with open(summary_path, "w") as f:
        out = lambda *a, **kw: print(*a, **kw, file=f) or print(*a, **kw)
        out("=" * 60)
        out("Drift Power Analysis")
        out(f"Input:  {csv_path.name}")
        out("=" * 60)
        out(f"\nData summary")
        out(f"  Time points (TPs):          {n_tps}")
        out(f"  Reps per TP (median):       {int(np.median(n_reps))}")
        out(f"  Total observations:         {N_total}")
        out(f"  E(λ=1) range:              [{tp_means.min():.3f}, {tp_means.max():.3f}]")
        out(f"  E(λ=1) mean:               {tp_means.mean():.3f}")
        out()
        out("Temporal dependence")
        out(f"  ICC (one-way, type 1):      {icc:.4f}")
        out(f"  Durbin-Watson:              {dw:.3f}  (2=no autocorr)")
        out(f"  Lag-1 autocorrelation:      {rhos[0]:.4f}")
        if len(rhos) > 1:
            out(f"  Lag-5 autocorrelation:      {rhos[min(4, len(rhos)-1)]:.4f}")
        out()
        out("Observed ZNE effect size")
        out(f"  Observed ZNE d range:       [{min(cohen_ds):.2f}, {max(cohen_ds):.2f}]")
        out(f"  Observed ZNE d mean:        {np.mean(cohen_ds):.2f}")
        out()
        out("Interpretation")
        if icc > 0.5:
            out(f"  HIGH ICC ({icc:.2f}): most variance is between TPs (drift-dominated).")
        elif icc > 0.1:
            out(f"  MODERATE ICC ({icc:.2f}): notable between-TP variance.")
        else:
            out(f"  LOW ICC ({icc:.2f}): within-TP variance dominates.")
        out()
        out(f"Outputs: {per_tp_path.name}, {summary_path.name}")

    print(f"\nWritten: {per_tp_path}")
    print(f"Written: {summary_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV,
                        help="Input raw_data CSV (default: %(default)s)")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                        help="Output directory (default: %(default)s)")
    args = parser.parse_args()

    if not args.csv.exists():
        sys.exit(f"ERROR: CSV not found: {args.csv}")

    analyse(args.csv, args.outdir)


if __name__ == "__main__":
    main()
