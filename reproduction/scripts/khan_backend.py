#!/usr/bin/env python3
"""
Khan et al. (2024): Parameter-Space Study
=====================================================================

The comprehensive, final Khan reproduction for the QCE 2026 paper.
Four backends, full parameter sweep, shot-count and repetition-count
sensitivity analysis.

Backends
--------
1. noiseless        — QASM with no noise (Khan's SIM_IBMQ reference)
2. depol_kyoto      — Depolarizing p_2q=0.036 (Khan's documented Kyoto EPLG 3.6%)
3. fake_osaka       — FakeOsaka calibration snapshot (Feb 2024, with readout errors)
4. fake_kyoto       — FakeKyoto calibration snapshot (BROKEN: all ECR errors = 1.0)

Parameter Sweep (one-at-a-time from defaults)
----------------------------------------------
A. folding_strategy    ∈ {from_left, from_right, global}
B. extrapolation_method∈ {richardson, linear, polynomial, exponential}
C. scale_factors       ∈ {[1,3,5], [1,2,3], [1,1.5,2,2.5,3]}
D. transpiler_level    ∈ {1, 3}
E. n_shots             ∈ {1024, 4096, 8192}

Default configuration:
  folding=from_left, extrap=richardson, scales=[1,3,5], transpiler=1,
  n_shots=4096

Shot-Count Sensitivity (post-hoc resampling from QASM probabilities)
---------------------------------------------------------------------
n_shots ∈ {128, 256, 512, 1024, 2048, 4096, 8192}
Applied at the default ZNE config only.

Statistical Power (bootstrap subsampling from N=200)
----------------------------------------------------
n_reps ∈ {5, 10, 20, 30, 50, 100, 200}  x 1000 bootstrap draws each
Applied at the default ZNE config only.

Output CSVs
-----------
khan_summary.csv   — 132 rows (11 configs x 4 backends x 3 TCs)
khan_detail.csv    — 26,400 rows (132 x 200 reps)
khan_shots.csv     — 84 rows (7 shots x 4 backends x 3 TCs)
khan_power.csv     — 84 rows (7 n_reps x 4 backends x 3 TCs)

Usage
-----
    cd reproduction
    python scripts/khan_unified.py
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap
from qiskit_aer import AerSimulator
from scipy import stats as sp_stats

# ── path setup ──
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(REPO_DIR))

from core.circuits import (
    build_khan_trotter,
    calibrate_angles,
    compute_ideal_expectation,
    sample_shot_noise,
)
from core.noise import get_fake_backend, make_noise_model
from core.stats import paired_analysis
from core.zne import extrapolate, fold_circuit, sigma_ci


# ═════════════════════════════════════════════════════════════════════════
#  Constants
# ═════════════════════════════════════════════════════════════════════════

BACKENDS = ["noiseless", "depol_kyoto", "fake_osaka", "fake_kyoto"]

BACKEND_LABELS = {
    "noiseless":   "Noiseless (QASM)",
    "depol_kyoto": "Depol. Kyoto (p₂q=0.036)",
    "fake_osaka":  "FakeOsaka snapshot",
    "fake_kyoto":  "FakeKyoto snapshot",
}

# Khan's documented EPLG values (Table 1 in the paper)
EPLG_PARAMS = {
    "depol_kyoto": {"p_1q": 0.003, "p_2q": 0.036},  # Kyoto EPLG 3.6%
}

TROTTER_STEPS = [1, 3, 5]

# Default ZNE configuration
DEFAULTS = {
    "folding_strategy":     "from_left",
    "extrapolation_method": "richardson",
    "scale_factors":        [1.0, 3.0, 5.0],
    "transpiler_level":     1,
    "n_shots":              4096,
}

# One-at-a-time sweep axes
# Richardson is now default (community standard: Temme 2017, Giurgica-Tiron 2020).
# Transpiler levels reduced to {1,3} — levels 0 and 2 are near-identical to 1.
# Shots added as sweep axis since Khan leaves Q undocumented.
SWEEPS: dict[str, list] = {
    "folding_strategy":     ["from_left", "from_right", "global"],
    "extrapolation_method": ["richardson", "linear", "polynomial", "exponential"],
    "scale_factors":        [[1, 3, 5], [1, 2, 3], [1, 1.5, 2, 2.5, 3]],
    "transpiler_level":     [1, 3],
    "n_shots":              [1024, 4096, 8192],
}

N_SHOTS_DEFAULT = 4096
N_REPS_DEFAULT = 200
SEED = 42

# Shot-count sweep values
SHOTS_SWEEP = [128, 256, 512, 1024, 2048, 4096, 8192]

# Repetition-count sweep values (for power analysis)
REPS_SWEEP = [5, 10, 20, 30, 50, 100, 200]
N_BOOTSTRAP = 1000

BASIS_GATES_CX = ["cx", "id", "rz", "sx", "x"]
BASIS_GATES_ECR = ["ecr", "id", "rz", "sx", "x"]


# ═════════════════════════════════════════════════════════════════════════
#  Backend cache
# ═════════════════════════════════════════════════════════════════════════

_BACKEND_CACHE: dict[str, tuple] = {}


def _get_backend(name: str):
    """Return (noise_model | None, basis_gates, coupling_map_4q | None)."""
    if name in _BACKEND_CACHE:
        return _BACKEND_CACHE[name]

    if name == "noiseless":
        result = (None, BASIS_GATES_CX, None)

    elif name == "depol_kyoto":
        p = EPLG_PARAMS[name]
        nm = make_noise_model("depolarizing", p["p_1q"], p["p_2q"])
        result = (nm, BASIS_GATES_CX, None)

    elif name in ("fake_osaka", "fake_kyoto"):
        backend, nm = get_fake_backend(name)
        full_cm = backend.coupling_map
        edges = {(i, j) for i, j in full_cm.get_edges() if i <= 3 and j <= 3}
        edges |= {(j, i) for i, j in edges}
        cm = CouplingMap(sorted(edges))
        result = (nm, BASIS_GATES_ECR, cm)

    else:
        raise ValueError(f"Unknown backend: {name}")

    _BACKEND_CACHE[name] = result
    return result


# ═════════════════════════════════════════════════════════════════════════
#  QASM batch execution with probability caching
# ═════════════════════════════════════════════════════════════════════════

# Cache: (backend, tc_name, transpiler, folding, scale_factor) → expectation
_QASM_CACHE: dict[tuple, float] = {}


def _qasm_expectation(circuit: QuantumCircuit, noise_model, n_total: int,
                      seed: int) -> float:
    """Run QASM simulation and return ⟨Z⊗…⊗Z⟩."""
    qc = circuit.copy()
    qc.measure_all()
    kw = {"noise_model": noise_model} if noise_model else {}
    sim = AerSimulator(**kw)
    result = sim.run(qc, shots=n_total, seed_simulator=seed).result()
    counts = result.get_counts()
    total = exp_val = 0.0
    for bitstring, count in counts.items():
        parity = (-1) ** bitstring.replace(" ", "").count("1")
        exp_val += parity * count
        total += count
    return exp_val / total


def get_qasm_expectations(
    base_circuit: QuantumCircuit,
    backend_name: str,
    tc_name: str,
    transpiler_level: int,
    folding_strategy: str,
    scale_factors: list[float],
    n_total: int = N_SHOTS_DEFAULT * N_REPS_DEFAULT,
) -> list[float]:
    """Get QASM expectations at each scale factor, with caching.

    Multiple extrapolation methods sharing the same (backend, TC,
    transpiler, folding, scales) will reuse cached QASM results.
    """
    nm, bg, cm = _get_backend(backend_name)
    is_fake = backend_name in ("fake_kyoto", "fake_osaka")

    # Transpile once
    cache_key_base = (backend_name, tc_name, transpiler_level, folding_strategy)
    expectations = []

    for s_idx, lam in enumerate(scale_factors):
        cache_key = cache_key_base + (lam,)
        if cache_key in _QASM_CACHE:
            expectations.append(_QASM_CACHE[cache_key])
            continue

        # Transpile
        if is_fake:
            tc = transpile(base_circuit, basis_gates=bg, coupling_map=cm,
                           optimization_level=transpiler_level,
                           seed_transpiler=SEED)
        else:
            tc = transpile(base_circuit, AerSimulator(),
                           optimization_level=transpiler_level,
                           basis_gates=bg, seed_transpiler=SEED)

        # Fold
        folded = fold_circuit(tc, lam, folding_strategy)

        # Run QASM
        seed = SEED + hash(cache_key) % 100_000
        ev = _qasm_expectation(folded, nm, n_total, seed)

        _QASM_CACHE[cache_key] = ev
        expectations.append(ev)

    return expectations


# ═════════════════════════════════════════════════════════════════════════
#  Resampling: generate n_reps shot-noise samples from QASM probabilities
# ═════════════════════════════════════════════════════════════════════════

def resample_reps(
    qasm_expectations: list[float],
    scale_factors: list[float],
    extrapolation_method: str,
    n_shots: int,
    n_reps: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate (raw, mitigated) arrays from cached QASM expectations."""
    raw = np.zeros(n_reps)
    mit = np.zeros(n_reps)
    for rep in range(n_reps):
        shot_vals = [sample_shot_noise(ev, n_shots, rng) for ev in qasm_expectations]
        raw[rep] = shot_vals[0]  # scale=1 → unmitigated
        mit[rep] = extrapolate(scale_factors, shot_vals, extrapolation_method)
    return raw, mit


# ═════════════════════════════════════════════════════════════════════════
#  Build the one-at-a-time configuration list
# ═════════════════════════════════════════════════════════════════════════

@dataclass
class ZNEConfig:
    """One ZNE parameter configuration."""
    folding_strategy: str
    extrapolation_method: str
    scale_factors: list[float]
    transpiler_level: int
    n_shots: int
    param_axis: str       # which axis is being varied
    param_value: str      # human-readable value of the varied axis
    is_default: bool      # True for the default config


def build_config_list() -> list[ZNEConfig]:
    """Build 11 one-at-a-time configs (1 default + 10 non-default)."""
    configs = []
    seen = set()

    for axis_name, axis_values in SWEEPS.items():
        for val in axis_values:
            # Start from defaults, override this axis
            cfg = dict(DEFAULTS)
            cfg[axis_name] = val

            # Normalise scale_factors to list[float]
            if isinstance(cfg["scale_factors"], list):
                sf_tuple = tuple(float(x) for x in cfg["scale_factors"])
            else:
                sf_tuple = (cfg["scale_factors"],)

            # De-duplicate
            key = (cfg["folding_strategy"], cfg["extrapolation_method"],
                   sf_tuple, cfg["transpiler_level"], cfg["n_shots"])
            if key in seen:
                continue
            seen.add(key)

            is_default = all(
                cfg[k] == DEFAULTS[k] or
                (isinstance(cfg[k], list) and cfg[k] == DEFAULTS[k])
                for k in DEFAULTS
            )

            val_str = (str(val) if not isinstance(val, list)
                       else "[" + ",".join(str(v) for v in val) + "]")

            configs.append(ZNEConfig(
                folding_strategy=cfg["folding_strategy"],
                extrapolation_method=cfg["extrapolation_method"],
                scale_factors=list(float(x) for x in cfg["scale_factors"])
                    if isinstance(cfg["scale_factors"], list) else cfg["scale_factors"],
                transpiler_level=cfg["transpiler_level"],
                n_shots=cfg["n_shots"],
                param_axis="default" if is_default else axis_name,
                param_value="default" if is_default else val_str,
                is_default=is_default,
            ))

    return configs


# ═════════════════════════════════════════════════════════════════════════
#  Phase A: Main parameter sweep
# ═════════════════════════════════════════════════════════════════════════

def run_phase_a(circuits, ideals, configs, n_reps, outdir):
    """Main sweep: 11 configs x 4 backends x 3 TCs."""
    total = len(configs) * len(BACKENDS) * len(TROTTER_STEPS)
    print(f"\n[Phase A] Parameter sweep: {len(configs)} configs x "
          f"{len(BACKENDS)} backends x {len(TROTTER_STEPS)} TCs = {total}")

    summary_rows = []
    detail_rows = []
    idx = 0
    t_start = time.time()

    for backend in BACKENDS:
        print(f"\n  ── {backend} ──")
        for tc_name, circuit in circuits.items():
            ideal = ideals[tc_name]
            for cfg in configs:
                idx += 1
                elapsed = time.time() - t_start
                eta = (elapsed / idx * (total - idx)) if idx > 0 else 0

                tag = f"{tc_name}|{cfg.param_axis}={cfg.param_value}"
                print(f"  [{idx:3d}/{total}] {tag:45s} ", end="", flush=True)

                t0 = time.time()

                # Get QASM expectations (cached)
                expectations = get_qasm_expectations(
                    circuit, backend, tc_name,
                    cfg.transpiler_level, cfg.folding_strategy,
                    cfg.scale_factors,
                )

                # Resample
                rng = np.random.default_rng(SEED + idx * 7)
                raw, mit = resample_reps(
                    expectations, cfg.scale_factors,
                    cfg.extrapolation_method, cfg.n_shots, n_reps, rng,
                )

                # Statistics
                stats = paired_analysis(raw, mit, ideal)
                sci = sigma_ci(cfg.scale_factors)
                dt = time.time() - t0

                # Summary row
                row = {
                    "backend": backend,
                    "circuit": tc_name,
                    "n_steps": int(tc_name[2:]),
                    "param_axis": cfg.param_axis,
                    "param_value": cfg.param_value,
                    "is_default": cfg.is_default,
                    "folding_strategy": cfg.folding_strategy,
                    "extrapolation_method": cfg.extrapolation_method,
                    "scale_factors": str(cfg.scale_factors),
                    "transpiler_level": cfg.transpiler_level,
                    "n_shots": cfg.n_shots,
                    "n_reps": n_reps,
                    "sigma_ci": sci,
                    "ideal": ideal,
                    **stats,
                }
                summary_rows.append(row)

                # Detail rows
                for rep in range(n_reps):
                    detail_rows.append({
                        "backend": backend,
                        "circuit": tc_name,
                        "param_axis": cfg.param_axis,
                        "param_value": cfg.param_value,
                        "rep": rep,
                        "raw_exp": raw[rep],
                        "mit_exp": mit[rep],
                        "ideal": ideal,
                    })

                sig = "***" if stats["significant"] and abs(stats["cohen_d"]) > 0.8 else \
                      "* " if stats["significant"] else "  "
                direction = "✓" if stats["mean_improvement"] > 0 else "✗"
                print(f"d={stats['cohen_d']:+6.2f} {sig}{direction} ({dt:.1f}s) "
                      f"ETA {eta/60:.0f}m")

    # Write
    _write_csv(outdir / "khan_summary.csv", summary_rows)
    _write_csv(outdir / "khan_detail.csv", detail_rows)
    print(f"\n  → summary: {len(summary_rows)} rows, detail: {len(detail_rows)} rows")

    return summary_rows, detail_rows


# ═════════════════════════════════════════════════════════════════════════
#  Phase B: Shot-count sensitivity (post-hoc resampling)
# ═════════════════════════════════════════════════════════════════════════

def run_phase_b(circuits, ideals, outdir):
    """n_shots sweep at default ZNE config, all backends x TCs."""
    total = len(SHOTS_SWEEP) * len(BACKENDS) * len(TROTTER_STEPS)
    print(f"\n[Phase B] Shot-count sweep: {len(SHOTS_SWEEP)} shots x "
          f"{len(BACKENDS)} x {len(TROTTER_STEPS)} = {total}")

    rows = []
    idx = 0

    for backend in BACKENDS:
        for tc_name, circuit in circuits.items():
            ideal = ideals[tc_name]

            # Get QASM expectations at default config (already cached from Phase A)
            expectations = get_qasm_expectations(
                circuit, backend, tc_name,
                DEFAULTS["transpiler_level"],
                DEFAULTS["folding_strategy"],
                DEFAULTS["scale_factors"],
            )

            for n_shots in SHOTS_SWEEP:
                idx += 1
                rng = np.random.default_rng(SEED + 50000 + idx * 13)

                raw, mit = resample_reps(
                    expectations, DEFAULTS["scale_factors"],
                    DEFAULTS["extrapolation_method"],
                    n_shots, N_REPS_DEFAULT, rng,
                )

                stats = paired_analysis(raw, mit, ideal)

                rows.append({
                    "backend": backend,
                    "circuit": tc_name,
                    "n_shots": n_shots,
                    "n_reps": N_REPS_DEFAULT,
                    "ideal": ideal,
                    **stats,
                })

                sig = "*" if stats["significant"] else " "
                print(f"  [{idx:2d}/{total}] {backend:15s} {tc_name} "
                      f"shots={n_shots:5d}  d={stats['cohen_d']:+6.2f}{sig}")

    _write_csv(outdir / "khan_shots.csv", rows)
    print(f"\n  → {len(rows)} rows")
    return rows


# ═════════════════════════════════════════════════════════════════════════
#  Phase C: Statistical power (bootstrap subsampling from N=200)
# ═════════════════════════════════════════════════════════════════════════

def run_phase_c(detail_rows, ideals, outdir):
    """n_reps power analysis via bootstrap from Phase A detail data."""
    total = len(REPS_SWEEP) * len(BACKENDS) * len(TROTTER_STEPS)
    print(f"\n[Phase C] Power analysis: {len(REPS_SWEEP)} n_reps x "
          f"{len(BACKENDS)} x {len(TROTTER_STEPS)} = {total} "
          f"({N_BOOTSTRAP} bootstraps each)")

    # Extract default-config detail data
    default_data: dict[tuple, tuple[np.ndarray, np.ndarray]] = {}
    for backend in BACKENDS:
        for tc_name in [f"TC{s}" for s in TROTTER_STEPS]:
            raws = []
            mits = []
            for r in detail_rows:
                if (r["backend"] == backend and r["circuit"] == tc_name
                        and r["param_axis"] == "default"):
                    raws.append(r["raw_exp"])
                    mits.append(r["mit_exp"])
            if raws:
                default_data[(backend, tc_name)] = (
                    np.array(raws), np.array(mits))

    rows = []
    idx = 0
    rng = np.random.default_rng(SEED + 99999)

    for backend in BACKENDS:
        for tc_name in [f"TC{s}" for s in TROTTER_STEPS]:
            key = (backend, tc_name)
            if key not in default_data:
                continue
            full_raw, full_mit = default_data[key]
            ideal = ideals[tc_name]
            n_full = len(full_raw)

            for n_reps in REPS_SWEEP:
                idx += 1
                n_sig = 0
                d_values = []

                for _ in range(N_BOOTSTRAP):
                    indices = rng.choice(n_full, size=n_reps, replace=True)
                    sub_raw = full_raw[indices]
                    sub_mit = full_mit[indices]

                    raw_err = np.abs(sub_raw - ideal)
                    mit_err = np.abs(sub_mit - ideal)
                    improvement = raw_err - mit_err

                    if n_reps >= 3:
                        _, p = sp_stats.ttest_rel(raw_err, mit_err)
                        std = improvement.std(ddof=1)
                        d = improvement.mean() / std if std > 0 else 0.0
                        if p < 0.05:
                            n_sig += 1
                        d_values.append(d)

                power = n_sig / N_BOOTSTRAP
                rows.append({
                    "backend": backend,
                    "circuit": tc_name,
                    "n_reps": n_reps,
                    "power": power,
                    "median_d": float(np.median(d_values)) if d_values else 0.0,
                    "mean_d": float(np.mean(d_values)) if d_values else 0.0,
                    "d_q25": float(np.percentile(d_values, 25)) if d_values else 0.0,
                    "d_q75": float(np.percentile(d_values, 75)) if d_values else 0.0,
                    "ideal": ideal,
                })

                print(f"  [{idx:2d}/{total}] {backend:15s} {tc_name} "
                      f"n_reps={n_reps:3d}  power={power:.3f}  "
                      f"d̃={rows[-1]['median_d']:+.2f}")

    _write_csv(outdir / "khan_power.csv", rows)
    print(f"\n  → {len(rows)} rows")
    return rows


# ═════════════════════════════════════════════════════════════════════════
#  Summary statistics
# ═════════════════════════════════════════════════════════════════════════

def print_verdict_summary(summary_rows):
    """Print verdict distribution by backend."""
    from collections import defaultdict

    print(f"\n{'=' * 72}")
    print("  VERDICT SUMMARY — by backend")
    print(f"{'=' * 72}")

    by_backend = defaultdict(list)
    for r in summary_rows:
        by_backend[r["backend"]].append(r)

    for backend in BACKENDS:
        rows = by_backend[backend]
        n = len(rows)
        n_better = sum(1 for r in rows
                       if r["mean_improvement"] > 0 and r["significant"])
        n_worse = sum(1 for r in rows
                      if r["mean_improvement"] < 0 and r["significant"])
        n_ns = sum(1 for r in rows if not r["significant"])
        median_d = float(np.median([r["cohen_d"] for r in rows]))

        print(f"\n  {BACKEND_LABELS[backend]:35s} ({n} configs, median d={median_d:+.2f})")
        print(f"    ✓ sig. better : {n_better:3d}  ({100*n_better/n:.0f}%)")
        print(f"    ✗ sig. worse  : {n_worse:3d}  ({100*n_worse/n:.0f}%)")
        print(f"    − not sig.    : {n_ns:3d}  ({100*n_ns/n:.0f}%)")

    # Active parameter analysis
    print(f"\n{'=' * 72}")
    print("  ACTIVE PARAMETER ANALYSIS — verdict flips across backends")
    print(f"{'=' * 72}")

    by_cfg = defaultdict(dict)
    for r in summary_rows:
        key = (r["circuit"], r["param_axis"], r["param_value"])
        verdict = ("better" if r["mean_improvement"] > 0 and r["significant"]
                   else "worse" if r["mean_improvement"] < 0 and r["significant"]
                   else "n.s.")
        by_cfg[key][r["backend"]] = verdict

    n_flips = 0
    for key, verdicts in sorted(by_cfg.items()):
        unique_verdicts = set(verdicts.values())
        if len(unique_verdicts) > 1:
            n_flips += 1
            tc, axis, val = key
            verdict_str = "  ".join(f"{b[:8]}={v}" for b, v in verdicts.items())
            print(f"  FLIP: {tc} {axis}={val:20s}  {verdict_str}")

    print(f"\n  Total configs with verdict flip across backends: "
          f"{n_flips}/{len(by_cfg)}")


# ═════════════════════════════════════════════════════════════════════════
#  I/O helpers
# ═════════════════════════════════════════════════════════════════════════

def _write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  📄 {path.name}")


# ═════════════════════════════════════════════════════════════════════════
#  Main
# ═════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-reps", type=int, default=N_REPS_DEFAULT)
    ap.add_argument("--shots", type=int, default=N_SHOTS_DEFAULT)
    ap.add_argument("--skip-phases", type=str, default="",
                    help="Comma-separated phases to skip: A,B,C")
    ap.add_argument("--outdir", type=str, default=None, required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    skip = set(args.skip_phases.upper().split(","))

    print("=" * 72)
    print("  Khan et al.: Parameter-Space Study")
    print("=" * 72)

    # ── Calibrate angles ─────────────────────────────────────────────
    # Pre-computed: calibrate_angles(target=0.8284, n_grid=80,
    #               trotter_steps=[1,2,3,4,5])
    # yields rx=0.097344, rz=0.133849 → mean(TC1..TC5) = 0.8284
    # TC1=0.981, TC2=0.927, TC3=0.846, TC4=0.746, TC5=0.641
    print("\n[0] Using pre-calibrated angles (mean TC1-TC5 ≈ 0.8284) …")
    t0 = time.time()
    rx_angle, rz_angle = 0.097344, 0.133849
    ideal_tc1 = compute_ideal_expectation(
        build_khan_trotter(4, 1, rx_angle, rz_angle))
    print(f"    rx={rx_angle:.6f}, rz={rz_angle:.6f}")
    print(f"    ideal(TC1) = {ideal_tc1:.6f}  ({time.time()-t0:.1f}s)")

    # ── Load backends ────────────────────────────────────────────────
    print("\n[0] Loading backends …")
    t0 = time.time()
    for b in BACKENDS:
        nm, bg, cm = _get_backend(b)
        print(f"    {b:20s} → {'OK' if nm or b == 'noiseless' else 'MISSING'}")
    print(f"    Done ({time.time()-t0:.1f}s)")

    # ── Build circuits ───────────────────────────────────────────────
    print("\n[0] Building Trotter circuits …")
    circuits: dict[str, QuantumCircuit] = {}
    ideals: dict[str, float] = {}
    for ns in TROTTER_STEPS:
        name = f"TC{ns}"
        qc = build_khan_trotter(4, ns, rx_angle, rz_angle)
        circuits[name] = qc
        ideals[name] = compute_ideal_expectation(qc)
        print(f"    {name}: ideal = {ideals[name]:.6f}")

    # ── Build config list ────────────────────────────────────────────
    configs = build_config_list()
    print(f"\n    {len(configs)} unique ZNE configs "
          f"(1 default + {len(configs)-1} variations)")

    # ── Phase A: Main sweep ──────────────────────────────────────────
    summary_rows = []
    detail_rows = []
    if "A" not in skip:
        summary_rows, detail_rows = run_phase_a(
            circuits, ideals, configs, args.n_reps, outdir)
        print_verdict_summary(summary_rows)
    else:
        print("\n[Phase A] SKIPPED")

    # ── Phase B: Shot-count sweep ────────────────────────────────────
    if "B" not in skip:
        run_phase_b(circuits, ideals, outdir)
    else:
        print("\n[Phase B] SKIPPED")

    # ── Phase C: Power analysis ──────────────────────────────────────
    if "C" not in skip and detail_rows:
        run_phase_c(detail_rows, ideals, outdir)
    else:
        print("\n[Phase C] SKIPPED (need Phase A detail data)")

    # ── Done ─────────────────────────────────────────────────────────
    total_time = time.time() - t0
    print(f"\n{'=' * 72}")
    print(f"  DONE — total wall time: {total_time/60:.1f} min")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
