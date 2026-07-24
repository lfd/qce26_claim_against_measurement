#!/usr/bin/env python3
"""
Full Statistical Battery — Drift ZNE Study
==========================================

Runs a comprehensive statistical analysis of each time point in a
longitudinal ZNE drift study.  Goes beyond the minimal paired test by
adding:

  - Wilcoxon signed-rank test (non-parametric alternative)
  - Cliff's δ  (non-parametric effect size, distribution-free)
  - Bootstrap 95% CI for Cohen's d and Cliff's δ
  - Durbin-Watson test on the per-TP Cohen's d timeseries
  - Global pooled analysis (treating all reps × TPs as one sample)
  - Summary statistics suitable for the workshop paper Table

Output CSVs
-----------
  drift_stats_per_tp.csv   — one row per time point, full test battery
  drift_stats_summary.csv  — global summary row
  drift_stats_report.txt   — human-readable report

Usage
-----
  cd reproduction
  python scripts/full_statistics.py
  python scripts/full_statistics.py --csv data/qexa_drift/raw_data_weekend.csv

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

# Richardson coefficients for λ = {1, 3, 5}
RICH_C1 = 15 / 8
RICH_C3 = -10 / 8
RICH_C5 = 3 / 8

N_BOOTSTRAP = 2000
RNG_SEED = 42


# ═══════════════════════════════════════════════════════════════════════
# Pure-numpy statistics (no scipy)
# ═══════════════════════════════════════════════════════════════════════

def _t_test_paired_pval(diffs: np.ndarray) -> tuple[float, float]:
    """One-sample t-test for H0: mean(diffs)=0.  Returns (t, p_two_sided)."""
    n = len(diffs)
    mu = diffs.mean()
    se = diffs.std(ddof=1) / math.sqrt(n)
    if se < 1e-15:
        return (float("inf") if mu > 0 else float("-inf"), 0.0)
    t = mu / se
    # Two-sided p-value via the beta-function CDF approximation (accurate for n>3)
    p = 2.0 * _t_cdf_upper(abs(t), df=n - 1)
    return t, p


def _t_cdf_upper(t: float, df: int) -> float:
    """P(T > t) for Student's t distribution using regularised incomplete beta."""
    x = df / (df + t * t)
    # regularised incomplete beta I_x(a, b) where a=df/2, b=0.5
    a, b = df / 2, 0.5
    return 0.5 * _betainc(x, a, b)


def _betainc(x: float, a: float, b: float) -> float:
    """Regularised incomplete beta I_x(a,b) via continued fraction (Lentz)."""
    if x < 0.0 or x > 1.0:
        return float("nan")
    if x == 0.0:
        return 0.0
    if x == 1.0:
        return 1.0
    # Use symmetry if needed for faster convergence
    if x > (a + 1) / (a + b + 2):
        return 1.0 - _betainc(1.0 - x, b, a)
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
    # Lentz continued fraction
    cf = _betacf(x, a, b)
    return min(front * cf, 1.0)


def _betacf(x: float, a: float, b: float, max_iter: int = 200, eps: float = 1e-10) -> float:
    """Continued fraction for regularised incomplete beta."""
    fpmin = 1e-300
    qab = a + b
    qap = a + 1
    qam = a - 1
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _wilcoxon_pval(diffs: np.ndarray) -> tuple[float, float]:
    """Wilcoxon signed-rank test for H0: median=0.  Returns (W, p_two_sided)."""
    diffs = diffs[diffs != 0]
    n = len(diffs)
    if n == 0:
        return 0.0, 1.0
    ranks = np.argsort(np.argsort(np.abs(diffs))) + 1.0
    W_plus = ranks[diffs > 0].sum()
    W = min(W_plus, n * (n + 1) / 2 - W_plus)
    # Normal approximation (valid for n ≥ 10)
    mu_W = n * (n + 1) / 4
    sigma_W = math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z = (W - mu_W) / sigma_W if sigma_W > 0 else 0.0
    p = 2.0 * _norm_cdf(-abs(z))
    return W_plus, p


def _norm_cdf(z: float) -> float:
    """Standard normal CDF via error function."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2)))


def cohen_d(diffs: np.ndarray) -> float:
    s = diffs.std(ddof=1)
    return float(diffs.mean() / s) if s > 1e-15 else 0.0


def cliffs_delta(raw_err: np.ndarray, mit_err: np.ndarray) -> float:
    """Cliff's δ = (N_better - N_worse) / n for paired data."""
    n = len(raw_err)
    diffs = raw_err - mit_err  # positive = mitigation helped
    n_better = (diffs > 0).sum()
    n_worse = (diffs < 0).sum()
    return float(n_better - n_worse) / n


def bootstrap_ci(statfn, data: np.ndarray, n: int = N_BOOTSTRAP,
                 alpha: float = 0.05, rng_seed: int = RNG_SEED):
    """Percentile bootstrap CI for a function of one 1-d array."""
    rng = np.random.default_rng(rng_seed)
    stats = [statfn(rng.choice(data, size=len(data), replace=True))
             for _ in range(n)]
    stats = np.sort(stats)
    lo = int(math.floor(alpha / 2 * n))
    hi = int(math.ceil((1 - alpha / 2) * n)) - 1
    return float(stats[lo]), float(stats[hi])


def durbin_watson(series: np.ndarray) -> float:
    diff = np.diff(series)
    denom = series @ series
    return float(diff @ diff / denom) if denom > 0 else 2.0


# ═══════════════════════════════════════════════════════════════════════
# Core analysis
# ═══════════════════════════════════════════════════════════════════════

def analyse_tp(e1: np.ndarray, e3: np.ndarray, e5: np.ndarray,
               ideal: float) -> dict:
    zne = RICH_C1 * e1 + RICH_C3 * e3 + RICH_C5 * e5
    raw_err = np.abs(e1 - ideal)
    mit_err = np.abs(zne - ideal)
    diffs = raw_err - mit_err  # positive = improved

    t, p_t = _t_test_paired_pval(diffs)
    W, p_w = _wilcoxon_pval(diffs)
    d = cohen_d(diffs)
    delta = cliffs_delta(raw_err, mit_err)
    d_lo, d_hi = bootstrap_ci(cohen_d, diffs)
    delta_lo, delta_hi = bootstrap_ci(
        lambda x: cliffs_delta(np.abs(x[:len(x)//2] - ideal),
                               np.abs(x[len(x)//2:] - ideal)),
        np.concatenate([e1, zne]))

    return dict(
        n=len(e1),
        mean_e1=float(e1.mean()),
        mean_zne=float(zne.mean()),
        mean_raw_err=float(raw_err.mean()),
        mean_mit_err=float(mit_err.mean()),
        t_stat=round(t, 4),
        p_t=round(p_t, 6),
        W_stat=round(W, 2),
        p_wilcoxon=round(p_w, 6),
        cohen_d=round(d, 4),
        cohen_d_ci_lo=round(d_lo, 4),
        cohen_d_ci_hi=round(d_hi, 4),
        cliffs_delta=round(delta, 4),
        cliffs_delta_ci_lo=round(delta_lo, 4),
        cliffs_delta_ci_hi=round(delta_hi, 4),
        significant_t=(p_t < 0.05),
        significant_w=(p_w < 0.05),
    )


def analyse(csv_path: Path, outdir: Path) -> None:
    rows = list(csv_mod.DictReader(open(csv_path, newline="")))
    ideal = float(rows[0]["ideal"])
    sf_col = "scale_factor"

    # Group by TP
    tp_map: dict[int, list] = {}
    for r in rows:
        tp_map.setdefault(int(r["timepoint_idx"]), []).append(r)
    tps = sorted(tp_map)

    per_tp: list[dict] = []
    cohen_ds: list[float] = []
    all_e1 = all_zne = all_raw_err = all_mit_err = None

    for tp in tps:
        tp_rows = tp_map[tp]
        e1 = np.array([float(r["exp_val"]) for r in tp_rows
                       if abs(float(r[sf_col]) - 1.0) < 0.1])
        e3 = np.array([float(r["exp_val"]) for r in tp_rows
                       if abs(float(r[sf_col]) - 3.0) < 0.1])
        e5 = np.array([float(r["exp_val"]) for r in tp_rows
                       if abs(float(r[sf_col]) - 5.0) < 0.1])

        if len(e1) != len(e3) or len(e1) != len(e5) or len(e1) == 0:
            continue

        res = analyse_tp(e1, e3, e5, ideal)
        res["timepoint_idx"] = tp
        per_tp.append(res)
        cohen_ds.append(res["cohen_d"])

        zne = RICH_C1 * e1 + RICH_C3 * e3 + RICH_C5 * e5
        raw_err = np.abs(e1 - ideal)
        mit_err = np.abs(zne - ideal)
        all_e1 = e1 if all_e1 is None else np.concatenate([all_e1, e1])
        all_zne = zne if all_zne is None else np.concatenate([all_zne, zne])
        all_raw_err = (raw_err if all_raw_err is None
                       else np.concatenate([all_raw_err, raw_err]))
        all_mit_err = (mit_err if all_mit_err is None
                       else np.concatenate([all_mit_err, mit_err]))

    # Global pooled analysis
    global_res = {}
    if all_e1 is not None:
        global_res = analyse_tp(all_e1, all_zne - RICH_C1 * all_e1,
                                all_e1 * 0, ideal)
        # redo using pre-computed arrays
        pooled_diff = all_raw_err - all_mit_err
        t_g, p_t_g = _t_test_paired_pval(pooled_diff)
        W_g, p_w_g = _wilcoxon_pval(pooled_diff)
        d_g = cohen_d(pooled_diff)
        delta_g = cliffs_delta(all_raw_err, all_mit_err)
        global_res = dict(
            n=len(all_e1),
            mean_e1=float(all_e1.mean()),
            mean_zne=float(all_zne.mean()),
            mean_raw_err=float(all_raw_err.mean()),
            mean_mit_err=float(all_mit_err.mean()),
            t_stat=round(t_g, 4),
            p_t=round(p_t_g, 6),
            W_stat=round(W_g, 2),
            p_wilcoxon=round(p_w_g, 6),
            cohen_d=round(d_g, 4),
            cliffs_delta=round(delta_g, 4),
        )

    cohen_ds_arr = np.array(cohen_ds)
    dw = durbin_watson(cohen_ds_arr)

    # ── Write per-TP CSV ──────────────────────────────────────────
    outdir.mkdir(parents=True, exist_ok=True)
    fields = ["timepoint_idx", "n", "mean_e1", "mean_zne",
              "mean_raw_err", "mean_mit_err",
              "t_stat", "p_t", "W_stat", "p_wilcoxon",
              "cohen_d", "cohen_d_ci_lo", "cohen_d_ci_hi",
              "cliffs_delta", "cliffs_delta_ci_lo", "cliffs_delta_ci_hi",
              "significant_t", "significant_w"]

    per_tp_path = outdir / "drift_stats_per_tp.csv"
    with open(per_tp_path, "w", newline="") as f:
        writer = csv_mod.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(per_tp)

    # ── Write summary CSV ────────────────────────────────────────
    summary_csv = outdir / "drift_stats_summary.csv"
    with open(summary_csv, "w", newline="") as f:
        writer = csv_mod.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerow(["n_tps", len(per_tp)])
        writer.writerow(["total_reps", global_res.get("n", "n/a")])
        writer.writerow(["cohen_d_min", f"{cohen_ds_arr.min():.4f}"])
        writer.writerow(["cohen_d_max", f"{cohen_ds_arr.max():.4f}"])
        writer.writerow(["cohen_d_mean", f"{cohen_ds_arr.mean():.4f}"])
        writer.writerow(["cohen_d_std", f"{cohen_ds_arr.std():.4f}"])
        writer.writerow(["cohen_d_range_ratio",
                         f"{cohen_ds_arr.max() / max(cohen_ds_arr.min(), 1e-9):.2f}"])
        writer.writerow(["durbin_watson_on_d", f"{dw:.4f}"])
        writer.writerow(["pct_tp_sig_t",
                         f"{sum(r['significant_t'] for r in per_tp)/len(per_tp):.4f}"])
        writer.writerow(["global_pooled_cohen_d",
                         f"{global_res.get('cohen_d', 'n/a')}"])
        writer.writerow(["global_pooled_cliffs_delta",
                         f"{global_res.get('cliffs_delta', 'n/a')}"])

    # ── Human-readable report ────────────────────────────────────
    report_path = outdir / "drift_stats_report.txt"
    with open(report_path, "w") as f:
        def out(*a, **kw):
            print(*a, **kw, file=f)
            print(*a, **kw)

        out("=" * 65)
        out("Full Statistical Battery — Longitudinal ZNE Drift Study")
        out(f"Input:  {csv_path.name}")
        out("=" * 65)
        out(f"\nTime points analysed:       {len(per_tp)}")
        out(f"Reps per TP (first):        {per_tp[0]['n']}")
        out(f"Ideal expectation:          {ideal:.6f}")
        out()
        out("Per-TP ZNE statistics (Richardson {1,3,5})")
        out("-" * 55)
        out(f"  Cohen's d range:          [{cohen_ds_arr.min():.3f}, {cohen_ds_arr.max():.3f}]")
        out(f"  Cohen's d mean ± std:     {cohen_ds_arr.mean():.3f} ± {cohen_ds_arr.std():.3f}")
        out(f"  Range ratio (max/min):    {cohen_ds_arr.max()/max(cohen_ds_arr.min(),1e-9):.2f}×")
        deltas = [r["cliffs_delta"] for r in per_tp]
        out(f"  Cliff's δ range:          [{min(deltas):.3f}, {max(deltas):.3f}]")
        out(f"  TPs significant (t-test): {sum(r['significant_t'] for r in per_tp)}/{len(per_tp)}")
        out(f"  TPs significant (Wilcox): {sum(r['significant_w'] for r in per_tp)}/{len(per_tp)}")
        out()
        out("Temporal structure of Cohen's d")
        out("-" * 55)
        out(f"  Durbin-Watson statistic:  {dw:.3f}  (2=no autocorr)")
        if dw < 1.5:
            out("  → Strong positive autocorrelation (drift)")
        elif dw > 2.5:
            out("  → Negative autocorrelation (oscillation)")
        else:
            out("  → Weak autocorrelation in d timeseries")
        out()
        out("Global pooled analysis (all reps × all TPs)")
        out("-" * 55)
        out(f"  N (total reps):           {global_res.get('n', 'n/a')}")
        out(f"  Pooled Cohen's d:         {global_res.get('cohen_d', 'n/a')}")
        out(f"  Pooled Cliff's δ:         {global_res.get('cliffs_delta', 'n/a')}")
        out(f"  p-value (t-test):         {global_res.get('p_t', 'n/a'):.2e}"
            if isinstance(global_res.get("p_t"), float) else "")
        out()
        out("Note: Pooled analysis inflates N by treating correlated TPs")
        out("as independent.  Per-TP analysis is the correct level.")
        out()
        out(f"Outputs: {per_tp_path.name}, {summary_csv.name}, {report_path.name}")

    print(f"\nWritten: {per_tp_path}")
    print(f"Written: {summary_csv}")
    print(f"Written: {report_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
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
