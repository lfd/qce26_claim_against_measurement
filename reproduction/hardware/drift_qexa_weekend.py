#!/usr/bin/env python3
"""
Experiment: ZNE Weekend Drift — IQM QExa (LRZ)  [48-HOUR CONTINUOUS RUN]
=========================================================================

Extended version of exp_drift_qexa.py for a 48-hour continuous drift study
over the weekend.  Key differences from the 12-hour version:

  1. **Crash-resume**:  On startup, scans existing CSV for already-collected
     time-points and resumes from the next one.  Safe to restart after
     network hiccups, server reboots, or QPU maintenance.

  2. **Connection retry**:  If a QPU submission fails, retries up to 5 times
     with exponential backoff before skipping the time-point.

  3. **Heartbeat logging**:  Every time-point prints a one-line heartbeat
     with elapsed time, ETA, and mean E(λ=1).

  4. **Separate output directory**: ``results/drift_qexa_weekend/`` to avoid
     overwriting the 12-hour data.

Protocol
--------
  Same as the 12-hour study:
  - TC1 Khan Trotter circuit, N=30 reps, 4096 shots, λ ∈ {1,3,5}
  - 30-minute intervals
  - Total: 97 TPs x 30 x 3 x 4096 ≈ 35.8M shots

Usage
-----
  # Production run (deploy in tmux):
  python scripts/exp_drift_qexa_weekend.py

  # Override duration:
  python scripts/exp_drift_qexa_weekend.py --duration-hours 72

  # Quick smoke test (2 reps, 128 shots, 3 TPs):
  python scripts/exp_drift_qexa_weekend.py --test

  # Local simulator test:
  python scripts/exp_drift_qexa_weekend.py --local-test --duration-hours 0.05
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from qiskit import transpile

# ── path setup ──
SCRIPT_DIR = Path(__file__).resolve().parent   # reproduction/hardware/
REPO_DIR = SCRIPT_DIR.parent                    # reproduction/
sys.path.insert(0, str(REPO_DIR))

from core.circuits import build_khan_trotter, compute_ideal_expectation
from core.zne import fold_circuit


# ═══════════════════════════════════════════════════════════════════════
#  Constants
# ═══════════════════════════════════════════════════════════════════════

BACKEND_NAME = "EQE1"
TROTTER_STEPS = 1
N_QUBITS = 4
N_REPS_DEFAULT = 30
N_SHOTS_DEFAULT = 4096
SCALE_FACTORS = [1.0, 3.0, 5.0]
FOLDING_STRATEGY = "from_left"

# Pre-calibrated angles (from khan_backend.py)
RX_ANGLE = 0.097344
RZ_ANGLE = 0.133849

# Defaults
DURATION_HOURS_DEFAULT = 48.0
INTERVAL_MINUTES_DEFAULT = 30.0
EXPERIMENT_NAME_DEFAULT = "drift_qexa_weekend"

# Retry settings
MAX_RETRIES = 5
RETRY_BASE_DELAY_S = 30

# Directories — set at module level, overridden by --experiment-name
RESULTS_DIR = REPO_DIR / "results" / EXPERIMENT_NAME_DEFAULT
LOGS_DIR = REPO_DIR / "logs"


# ═══════════════════════════════════════════════════════════════════════
#  Logging
# ═══════════════════════════════════════════════════════════════════════

def setup_logging(log_dir: Path, experiment_name: str = EXPERIMENT_NAME_DEFAULT) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(experiment_name)
    logger.setLevel(logging.DEBUG)
    if logger.handlers:
        return logger

    fh = logging.FileHandler(log_dir / f"{experiment_name}.log", mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter(
        "%(asctime)s | %(message)s", datefmt="%H:%M:%S",
    ))
    logger.addHandler(ch)

    return logger


# Initial logger — may be reconfigured by main() with --experiment-name
log = setup_logging(LOGS_DIR)


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════

def expectation_from_counts(counts: dict[str, int], n_qubits: int) -> float:
    """Compute ⟨Z⊗…⊗Z⟩ from measurement counts."""
    total = sum(counts.values())
    exp_val = 0.0
    for bitstring, count in counts.items():
        parity = (-1) ** bitstring.replace(" ", "").count("1")
        exp_val += parity * count
    return exp_val / total


def _csv_append(path: Path, row: dict) -> None:
    """Append a single row to a CSV, writing headers if file is new."""
    exists = path.exists() and path.stat().st_size > 0
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not exists:
            w.writeheader()
        w.writerow(row)


def get_completed_timepoints(results_dir: Path) -> set[int]:
    """Scan existing CSV to find already-completed time-point indices."""
    csv_path = results_dir / "raw_data.csv"
    if not csv_path.exists():
        return set()

    completed = set()
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            completed.add(int(row["timepoint_idx"]))

    return completed


def get_qexa_backend(backend_name: str = BACKEND_NAME):
    """Get QExa backend via MQSS adapter with token from env or .env file."""
    token = os.environ.get("MQSS_TOKEN")
    if not token:
        env_file = REPO_DIR / ".env"
        if env_file.exists():
            for line in env_file.read_text().splitlines():
                if line.startswith("MQSS_TOKEN="):
                    token = line.split("=", 1)[1].strip().strip("'\"")
                    break
    if not token:
        raise RuntimeError(
            "MQSS_TOKEN not found. Set MQSS_TOKEN env-var "
            "or add MQSS_TOKEN=... to .env"
        )

    from mqss.qiskit_adapter import MQSSQiskitAdapter
    adapter = MQSSQiskitAdapter(token=token)
    backend = adapter.get_backend(backend_name)
    log.info(f"Connected to {backend_name}")
    return backend


# ═══════════════════════════════════════════════════════════════════════
#  Circuit Preparation
# ═══════════════════════════════════════════════════════════════════════

def prepare_circuits(backend=None) -> dict:
    """Build, transpile, fold circuits.  Return reusable circuit set."""
    qc_base = build_khan_trotter(
        n_qubits=N_QUBITS, n_steps=TROTTER_STEPS,
        rx_angle=RX_ANGLE, rz_angle=RZ_ANGLE,
    )
    ideal = compute_ideal_expectation(qc_base)
    log.info(f"Ideal ⟨ZZZZ⟩ = {ideal:.6f}")

    if backend is not None:
        qc_t = transpile(qc_base, backend=backend, optimization_level=1)
    else:
        qc_t = transpile(qc_base, optimization_level=1)

    log.info(f"Transpiled: {qc_t.size()} gates, "
             f"{qc_t.count_ops().get('cz', 0)} CZ, depth {qc_t.depth()}")

    folded = {}
    for lam in SCALE_FACTORS:
        qc_f = (fold_circuit(qc_t, scale_factor=lam, strategy=FOLDING_STRATEGY)
                if lam > 1 else qc_t.copy())
        qc_f.measure_all()
        folded[lam] = qc_f
        cz = qc_f.count_ops().get("cz", 0)
        log.info(f"  λ={lam}: {qc_f.size()} gates, {cz} CZ")

    backend_name = (backend.name if hasattr(backend, "name")
                    else "local_simulator")

    return {
        "ideal": ideal,
        "base": qc_t,
        "folded": folded,
        "scale_factors": SCALE_FACTORS,
        "backend_name": backend_name,
    }


# ═══════════════════════════════════════════════════════════════════════
#  Data Collection — QExa Hardware (with retry)
# ═══════════════════════════════════════════════════════════════════════

def collect_timepoint_qexa(
    backend,
    circuit_set: dict,
    n_reps: int,
    n_shots: int,
    timepoint_idx: int,
) -> dict:
    """Submit circuits for one time-point with retry logic."""
    folded = circuit_set["folded"]
    scale_factors = circuit_set["scale_factors"]
    ts = datetime.now(timezone.utc).isoformat()

    log.info(f"=== TP {timepoint_idx} @ {ts} | "
             f"{n_reps} reps x {len(scale_factors)} λ x {n_shots} shots ===")

    all_circuits = []
    circuit_map = []
    for rep in range(n_reps):
        for li, lam in enumerate(scale_factors):
            all_circuits.append(folded[lam])
            circuit_map.append((rep, li))

    # Retry loop
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            log.info(f"  Submitting {len(all_circuits)} circuits "
                     f"(attempt {attempt}/{MAX_RETRIES}) …")
            job = backend.run(all_circuits, shots=n_shots, no_modify=True)
            log.info(f"  Job ID: {job.job_id()}")
            result = job.result()
            log.info(f"  Results received.")
            break
        except Exception as e:
            log.warning(f"  Attempt {attempt} failed: {e}")
            if attempt == MAX_RETRIES:
                log.error(f"  ALL {MAX_RETRIES} ATTEMPTS FAILED for TP {timepoint_idx}. "
                          f"Skipping this time-point.")
                raise
            delay = RETRY_BASE_DELAY_S * (2 ** (attempt - 1))
            log.info(f"  Retrying in {delay}s …")
            time.sleep(delay)

    per_rep = {r: {} for r in range(n_reps)}
    for ci, (rep, li) in enumerate(circuit_map):
        counts = result.get_counts(ci)
        lam = scale_factors[li]
        exp_val = expectation_from_counts(counts, N_QUBITS)
        per_rep[rep][lam] = {"exp_val": exp_val, "counts": counts}

    return {
        "timepoint_idx": timepoint_idx,
        "timestamp": ts,
        "backend": circuit_set["backend_name"],
        "n_reps": n_reps,
        "n_shots": n_shots,
        "ideal": float(circuit_set["ideal"]),
        "scale_factors": scale_factors,
        "per_rep": per_rep,
    }


# ═══════════════════════════════════════════════════════════════════════
#  Data Collection — Local Simulator (synthetic drift)
# ═══════════════════════════════════════════════════════════════════════

def collect_timepoint_local(
    circuit_set: dict,
    n_reps: int,
    n_shots: int,
    timepoint_idx: int,
    drift_phase: float = 0.0,
) -> dict:
    """Simulate one time-point with depolarising noise + sinusoidal drift."""
    from qiskit_aer import AerSimulator
    from core.noise import make_noise_model

    scale_factors = circuit_set["scale_factors"]
    ts = datetime.now(timezone.utc).isoformat()

    p_2q = 0.0175 + 0.0125 * np.sin(drift_phase)
    noise_model = make_noise_model(
        "depolarizing", p_1q=p_2q / 10, p_2q=p_2q, n_qubits=N_QUBITS,
    )
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)

    log.info(f"=== TP {timepoint_idx} (local) @ {ts} | "
             f"phase={drift_phase:.2f}  p_2q={p_2q:.4f} ===")

    per_rep = {r: {} for r in range(n_reps)}

    for rep in range(n_reps):
        for lam in scale_factors:
            qc_f = circuit_set["folded"][lam]
            qc_dm = qc_f.remove_final_measurements(inplace=False)
            qc_dm.save_density_matrix()
            res = sim.run(qc_dm, shots=0).result()
            dm = res.data()["density_matrix"]
            probs = np.real(np.diag(dm))
            probs = np.maximum(probs, 0)
            probs /= probs.sum()
            raw_counts = np.random.multinomial(n_shots, probs)
            counts_dict = {
                format(i, f"0{N_QUBITS}b"): int(c)
                for i, c in enumerate(raw_counts) if c > 0
            }
            exp_val = 0.0
            for idx, c in enumerate(raw_counts):
                parity = (-1) ** format(idx, f"0{N_QUBITS}b").count("1")
                exp_val += parity * c
            exp_val /= n_shots

            per_rep[rep][lam] = {"exp_val": exp_val, "counts": counts_dict}

    return {
        "timepoint_idx": timepoint_idx,
        "timestamp": ts,
        "backend": "local_simulator",
        "n_reps": n_reps,
        "n_shots": n_shots,
        "ideal": float(circuit_set["ideal"]),
        "scale_factors": scale_factors,
        "drift_phase": float(drift_phase),
        "p_2q": float(p_2q),
        "per_rep": per_rep,
    }


# ═══════════════════════════════════════════════════════════════════════
#  Persist Results
# ═══════════════════════════════════════════════════════════════════════

def save_timepoint(tp_data: dict, results_dir: Path) -> None:
    """Write granular CSV rows + full JSON snapshot."""
    results_dir.mkdir(parents=True, exist_ok=True)

    # JSON snapshot
    ts_slug = tp_data["timestamp"].replace(":", "-").replace("+", "_")
    json_path = results_dir / f"tp_{tp_data['timepoint_idx']:04d}_{ts_slug}.json"

    serialisable = dict(tp_data)
    serialisable["per_rep"] = {
        str(rep): {str(lam): v for lam, v in lambdas.items()}
        for rep, lambdas in tp_data["per_rep"].items()
    }
    with open(json_path, "w") as f:
        json.dump(serialisable, f, indent=2)
    log.debug(f"  JSON → {json_path.name}")

    # Granular CSV
    csv_path = results_dir / "raw_data.csv"
    for rep, lambdas in tp_data["per_rep"].items():
        for lam, data in lambdas.items():
            row = {
                "timepoint_idx": tp_data["timepoint_idx"],
                "timestamp": tp_data["timestamp"],
                "backend": tp_data["backend"],
                "rep": rep,
                "scale_factor": lam,
                "exp_val": f"{data['exp_val']:.6f}",
                "n_shots": tp_data["n_shots"],
                "ideal": f"{tp_data['ideal']:.6f}",
            }
            if "drift_phase" in tp_data:
                row["drift_phase"] = f"{tp_data['drift_phase']:.4f}"
                row["p_2q"] = f"{tp_data['p_2q']:.6f}"
            _csv_append(csv_path, row)

    lam1_vals = [tp_data["per_rep"][r][1.0]["exp_val"]
                 for r in tp_data["per_rep"]]
    mean_raw = np.mean(lam1_vals)
    log.info(f"  Saved. mean(λ=1) = {mean_raw:.4f}")


# ═══════════════════════════════════════════════════════════════════════
#  Main Run Loop (with crash-resume)
# ═══════════════════════════════════════════════════════════════════════

def run_experiment(args):
    """Longitudinal data collection with crash-resume."""
    global RESULTS_DIR, log

    exp_name = args.experiment_name
    RESULTS_DIR = REPO_DIR / "results" / exp_name

    # Reconfigure logger if name differs from default
    if exp_name != EXPERIMENT_NAME_DEFAULT:
        log = setup_logging(LOGS_DIR, exp_name)

    n_reps = args.n_reps or N_REPS_DEFAULT
    n_shots = args.n_shots or N_SHOTS_DEFAULT
    interval_s = args.interval_minutes * 60
    duration_s = args.duration_hours * 3600
    n_tp = max(1, int(duration_s / interval_s)) + 1

    if args.test:
        n_reps, n_shots = 2, 128
        n_tp = min(n_tp, 3)
        log.info("TEST MODE: 2 reps, 128 shots, max 3 TPs")

    shots_per_tp = n_reps * len(SCALE_FACTORS) * n_shots
    total_shots = n_tp * shots_per_tp

    # Check for crash-resume
    completed_tps = get_completed_timepoints(RESULTS_DIR)
    n_skip = len(completed_tps)

    label = exp_name.upper().replace('_', ' ')
    log.info("╔═════════════════════════════════════════════════════════════╗")
    log.info(f"║  {label}")
    log.info(f"║  Backend: {'LOCAL' if args.local_test else args.backend}")
    log.info(f"║  {n_tp} time-points, "
             f"every {args.interval_minutes:.0f} min, "
             f"{args.duration_hours:.0f}h total                   ║")
    log.info(f"║  {n_reps} reps x {len(SCALE_FACTORS)} λ x {n_shots} shots "
             f"= {shots_per_tp:>7,} shots/tp            ║")
    log.info(f"║  Total: ~{total_shots:>12,.0f} shots "
             f"                              ║")
    if n_skip > 0:
        log.info(f"║  ⚡ RESUMING: {n_skip} TPs already completed, "
                 f"{n_tp - n_skip} remaining         ║")
    log.info("╚═════════════════════════════════════════════════════════════╝")

    # Prepare circuits
    if args.local_test:
        circuit_set = prepare_circuits(backend=None)
        backend = None
    else:
        backend = get_qexa_backend(args.backend)
        circuit_set = prepare_circuits(backend)

    t_start = time.time()
    n_collected = 0
    n_failed = 0

    for tp_idx in range(n_tp):
        # Skip already-collected TPs (crash-resume)
        if tp_idx in completed_tps:
            log.info(f"  TP {tp_idx} already collected, skipping.")
            continue

        elapsed_h = (time.time() - t_start) / 3600
        remaining_tps = n_tp - tp_idx - n_skip
        eta_h = (remaining_tps * interval_s / 3600) if n_collected == 0 else \
                (elapsed_h / max(n_collected, 1)) * remaining_tps

        try:
            if args.local_test:
                drift_phase = 2 * np.pi * tp_idx / max(n_tp - 1, 1)
                tp = collect_timepoint_local(
                    circuit_set, n_reps, n_shots, tp_idx, drift_phase,
                )
            else:
                tp = collect_timepoint_qexa(
                    backend, circuit_set, n_reps, n_shots, tp_idx,
                )

            save_timepoint(tp, RESULTS_DIR)
            n_collected += 1

            # Heartbeat
            lam1_vals = [tp["per_rep"][r][1.0]["exp_val"]
                         for r in tp["per_rep"]]
            log.info(f"  ❤ TP {tp_idx}/{n_tp-1} | "
                     f"elapsed {elapsed_h:.1f}h | "
                     f"ETA {eta_h:.1f}h | "
                     f"E(λ=1)={np.mean(lam1_vals):.4f} | "
                     f"collected={n_collected} failed={n_failed}")

        except Exception as e:
            n_failed += 1
            log.error(f"  ✗ TP {tp_idx} FAILED ({n_failed} total): {e}")
            log.debug(traceback.format_exc())

            # If too many consecutive failures, try reconnecting
            if n_failed >= 3 and not args.local_test:
                log.warning("  3+ failures — attempting to reconnect to QPU …")
                try:
                    backend = get_qexa_backend(args.backend)
                    circuit_set = prepare_circuits(backend)
                    log.info("  Reconnected successfully.")
                except Exception as re_err:
                    log.error(f"  Reconnection failed: {re_err}")

        # Sleep until next time-point
        if tp_idx < n_tp - 1:
            next_t = t_start + (tp_idx + 1 - n_skip) * interval_s
            # If resuming, adjust timing: sleep based on real interval
            if n_skip > 0:
                next_t = time.time() + interval_s
            sleep_s = max(0, next_t - time.time())
            if sleep_s > 0:
                log.info(f"  Sleeping {sleep_s:.0f}s until next TP …")
                time.sleep(sleep_s)

    # Final summary
    log.info("")
    log.info("═" * 60)
    log.info(f"{exp_name.upper()} COLLECTION COMPLETE")
    log.info(f"  Collected: {n_collected}/{n_tp} time-points")
    log.info(f"  Failed:    {n_failed}")
    log.info(f"  Duration:  {(time.time() - t_start)/3600:.1f}h")
    log.info(f"  Data:      {RESULTS_DIR / 'raw_data.csv'}")
    log.info("═" * 60)


# ═══════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description="ZNE Longitudinal Drift — IQM QExa (configurable duration)",
    )
    p.add_argument("--experiment-name", type=str,
                   default=EXPERIMENT_NAME_DEFAULT,
                   help=f"Experiment name — controls results dir and log file "
                        f"(default: {EXPERIMENT_NAME_DEFAULT})")
    p.add_argument("--backend", default=BACKEND_NAME,
                   help=f"QPU backend name (default: {BACKEND_NAME})")
    p.add_argument("--duration-hours", type=float,
                   default=DURATION_HOURS_DEFAULT,
                   help=f"Total run duration in hours (default: {DURATION_HOURS_DEFAULT})")
    p.add_argument("--interval-minutes", type=float,
                   default=INTERVAL_MINUTES_DEFAULT,
                   help=f"Minutes between time-points (default: {INTERVAL_MINUTES_DEFAULT})")
    p.add_argument("--n-reps", type=int, default=None,
                   help=f"Reps per time-point (default: {N_REPS_DEFAULT})")
    p.add_argument("--n-shots", type=int, default=None,
                   help=f"Shots per circuit (default: {N_SHOTS_DEFAULT})")
    p.add_argument("--local-test", action="store_true",
                   help="Local sim with synthetic drift (no QPU)")
    p.add_argument("--test", action="store_true",
                   help="Smoke-test mode (2 reps, 128 shots, 3 TPs)")

    args = p.parse_args()

    log.info(f"Starting {args.experiment_name}: backend={args.backend}, "
             f"duration={args.duration_hours}h, "
             f"interval={args.interval_minutes}min, "
             f"local={args.local_test}")

    run_experiment(args)


if __name__ == "__main__":
    main()
