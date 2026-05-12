#!/usr/bin/env python3
"""
Experiment: ZNE Temporal Drift — IBM Quantum (EU Backend)
=========================================================

Longitudinal ZNE drift study on IBM Quantum hardware, designed as the
IBM counterpart to our 72-hour QExa drift study.  Measures whether ZNE
verdicts are temporally stable when the lambda>1 signal is **preserved**
(unlike QExa, where it is destroyed).

Key design choices
------------------
- **ibm_aachen** (Heron, CZ basis): same architecture as ibm_marrakesh
  → direct comparison with existing cross-vendor data.
- **TC3** (3 Trotter steps, 18 CZ at λ=1): this depth showed d=+0.55
  on Marrakesh, placing it near the verdict boundary.
- **10 reps** per time-point: 89% power at d=1.0 (QExa effects were
  d=3-13, so even moderate IBM drift should be detectable).
- **30 PUBs per TP** (10 reps x 3λ): ~34s QPU per TP.
- The script uses IBM's batch-mode (SamplerV2 per time-point) rather than
session mode, so we only pay for QPU seconds, not wall-clock time.  The
30-minute interval between TPs is pure wall-clock waiting (free).

Crash-resume
------------
On startup, scans existing CSV for completed TPs and skips them.  Safe
to restart after network issues or queue delays.

Output
------
results/drift_ibm/raw_data.csv       — granular: 1 row per (TP x rep x λ)
results/drift_ibm/tp_NNNN_*.json     — full snapshot per time-point
results/drift_ibm/drift_log.csv      — summary: 1 row per TP
results/drift_ibm/calibration.json   — backend calibration snapshots
logs/drift_ibm.log                   — full log

Usage
-----
    cd reproduction
    # Local simulator test (no IBM needed):
    python scripts/exp_drift_ibm.py --local-test --duration-hours 0.05

    # Single time-point on real hardware:
    python scripts/exp_drift_ibm.py --mode single

    # Full longitudinal run (deploy in tmux on server):
    python scripts/exp_drift_ibm.py --mode longitudinal \\
        --duration-hours 8.5 --interval-minutes 30

    # Dry run (transpile only, no submission):
    python scripts/exp_drift_ibm.py --mode single --dry-run
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
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import BasisTranslator
from qiskit.circuit.equivalence_library import SessionEquivalenceLibrary

# ── path setup ──
SCRIPT_DIR = Path(__file__).resolve().parent   # reproduction/hardware/
REPO_DIR = SCRIPT_DIR.parent                    # reproduction/
sys.path.insert(0, str(REPO_DIR))

from core.circuits import (
    build_khan_trotter,
    calibrate_angles,
    compute_ideal_expectation,
    sample_shot_noise,
)
from core.noise import make_noise_model
from core.zne import extrapolate, fold_circuit, sigma_ci
from core.stats import paired_analysis


# ═══════════════════════════════════════════════════════════════════════
#  Constants
# ═══════════════════════════════════════════════════════════════════════

BACKEND_NAME = "ibm_aachen"
BACKEND_ALTERNATIVES = ["ibm_brussel", "ibm_strasbourg"]

TROTTER_STEPS = 3
N_QUBITS = 4
SCALE_FACTORS = [1.0, 3.0, 5.0]
FOLDING_STRATEGY = "from_left"
EXTRAPOLATION_METHOD = "richardson"
TRANSPILER_LEVEL = 1
SEED = 42

# Pre-calibrated angles (from khan_backend.py)
RX_ANGLE = 0.097344
RZ_ANGLE = 0.133849

N_REPS_DEFAULT = 10
N_SHOTS = 4096
N_TIMEPOINTS_DEFAULT = 17  # ~8.5h at 30-min intervals
INTERVAL_MINUTES_DEFAULT = 30.0
DURATION_HOURS_DEFAULT = 8.5

# IBM SamplerV2 limits
MAX_PUBS_PER_JOB = 300

# Retry settings
MAX_RETRIES = 5
RETRY_BASE_DELAY_S = 30

# Directories
RESULTS_DIR = REPO_DIR / "results" / "drift_ibm"
LOGS_DIR = REPO_DIR / "logs"


# ═══════════════════════════════════════════════════════════════════════
#  Logging
# ═══════════════════════════════════════════════════════════════════════

def setup_logging(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("drift_ibm")
    logger.setLevel(logging.DEBUG)
    if logger.handlers:
        return logger

    fh = logging.FileHandler(log_dir / "drift_ibm.log", mode="a")
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


def bitarray_to_counts(bit_array, n_qubits: int) -> dict[str, int]:
    """Convert a SamplerV2 BitArray to a counts dict."""
    if hasattr(bit_array, 'get_counts'):
        return bit_array.get_counts()
    # Fallback for older API
    counts = {}
    for shot_result in bit_array:
        bs = format(int(shot_result), f"0{n_qubits}b")
        counts[bs] = counts.get(bs, 0) + 1
    return counts


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


def get_ibm_backend(backend_name: str):
    """Connect to IBM Quantum and return the backend."""
    from qiskit_ibm_runtime import QiskitRuntimeService

    token = os.environ.get("IBM_QUANTUM_TOKEN")
    if not token:
        env_file = REPO_DIR.parent / ".env"
        if env_file.exists():
            for line in env_file.read_text().splitlines():
                if line.startswith("IBM_QUANTUM_TOKEN="):
                    token = line.split("=", 1)[1].strip().strip("'\"")
                    break
    if token:
        service = QiskitRuntimeService(token=token, channel="ibm_quantum_platform")
    else:
        service = QiskitRuntimeService()

    log.info(f"Connecting to IBM Quantum ({backend_name}) …")
    backend = service.backend(backend_name)
    status = backend.status()
    log.info(f"  Status: operational={status.operational}, "
             f"pending_jobs={status.pending_jobs}")
    return backend, service


def get_backend_calibration_snapshot(backend) -> dict:
    """Extract key calibration metrics from an IBM backend."""
    snapshot = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "backend_name": backend.name,
    }

    try:
        props = backend.properties()
    except Exception:
        props = None

    if props is None:
        snapshot["calibration_date"] = None
        snapshot["gate_errors"] = {}
        snapshot["readout_errors"] = {}
        snapshot["t1_t2"] = {}
        return snapshot

    snapshot["calibration_date"] = str(props.last_update_date)

    # 2Q gate errors for qubits 0-3
    gate_errors = {}
    for q1, q2 in [(0, 1), (1, 2), (2, 3)]:
        for gate_name in ["cz", "ecr", "cx"]:
            try:
                err = props.gate_error(gate_name, [q1, q2])
                gate_errors[f"{gate_name}_{q1}_{q2}"] = err
            except Exception:
                pass
    snapshot["gate_errors"] = gate_errors

    # Readout errors
    readout_errors = {}
    for q in range(N_QUBITS):
        try:
            readout_errors[f"q{q}"] = props.readout_error(q)
        except Exception:
            pass
    snapshot["readout_errors"] = readout_errors

    # T1, T2
    t1_t2 = {}
    for q in range(N_QUBITS):
        try:
            t1_t2[f"T1_q{q}"] = props.t1(q)
            t1_t2[f"T2_q{q}"] = props.t2(q)
        except Exception:
            pass
    snapshot["t1_t2"] = t1_t2

    return snapshot


# ═══════════════════════════════════════════════════════════════════════
#  Circuit Preparation
# ═══════════════════════════════════════════════════════════════════════

def prepare_circuits(backend=None, local_test: bool = False) -> dict:
    """Build, transpile, fold circuits.  Done once, reused across TPs."""
    tc_name = f"TC{TROTTER_STEPS}"
    qc = build_khan_trotter(N_QUBITS, TROTTER_STEPS, RX_ANGLE, RZ_ANGLE)
    ideal = compute_ideal_expectation(qc)
    log.info(f"  {tc_name}: ideal ⟨ZZZZ⟩ = {ideal:.6f}")

    if local_test:
        from qiskit_aer import AerSimulator
        sim = AerSimulator()
        qc_t = transpile(
            qc, sim,
            optimization_level=TRANSPILER_LEVEL,
            basis_gates=["cx", "id", "rz", "sx", "x"],
            seed_transpiler=SEED,
        )
    else:
        qc_t = transpile(
            qc, backend,
            optimization_level=TRANSPILER_LEVEL,
            seed_transpiler=SEED,
            initial_layout=list(range(N_QUBITS)),
        )

    log.info(f"  Transpiled: {qc_t.size()} gates, depth {qc_t.depth()}")

    # Fold at each scale factor
    folded = {}
    for lam in SCALE_FACTORS:
        fc = fold_circuit(qc_t, lam, FOLDING_STRATEGY)
        if not local_test and lam > 1.0:
            # Re-translate folded gates back to backend basis
            basis_pm = PassManager([
                BasisTranslator(SessionEquivalenceLibrary,
                                backend.operation_names)
            ])
            fc = basis_pm.run(fc)
        folded[lam] = fc

        n_2q = sum(v for k, v in fc.count_ops().items()
                   if k in ("cz", "ecr", "cx"))
        log.info(f"  λ={lam}: {fc.size()} gates, {n_2q} 2Q, "
                 f"depth {fc.depth()}")

    return {
        "ideal": ideal,
        "tc_name": tc_name,
        "folded": folded,
        "transpiled_gates": qc_t.size(),
        "transpiled_depth": qc_t.depth(),
    }


# ═══════════════════════════════════════════════════════════════════════
#  Single Time-Point: IBM Hardware
# ═══════════════════════════════════════════════════════════════════════

def collect_timepoint_hw(
    backend,
    circuit_info: dict,
    n_reps: int,
    timepoint_idx: int,
    results_dir: Path,
) -> dict:
    """Execute one TP on IBM hardware via SamplerV2 (batch mode)."""
    from qiskit.circuit import ClassicalRegister
    from qiskit_ibm_runtime import SamplerV2

    timestamp = datetime.now(timezone.utc)
    ts_iso = timestamp.isoformat()
    ts_slug = timestamp.strftime("%Y%m%dT%H%M%SZ")
    ideal = circuit_info["ideal"]
    folded = circuit_info["folded"]
    tc_name = circuit_info["tc_name"]

    log.info(f"═══ TP {timepoint_idx} @ {ts_slug} | "
             f"{n_reps} reps x {len(SCALE_FACTORS)} λ ═══")

    # Prepare circuits with measurements
    folded_with_meas = {}
    for lam, fc in folded.items():
        fc_m = fc.copy()
        if fc_m.num_clbits == 0:
            cr = ClassicalRegister(N_QUBITS, 'meas')
            fc_m.add_register(cr)
            fc_m.measure(list(range(N_QUBITS)), cr)
        folded_with_meas[lam] = fc_m

    # Build PUBs: n_reps x len(SCALE_FACTORS)
    pubs = []
    pub_map = []  # (rep_idx, lambda)
    for rep in range(n_reps):
        for lam in SCALE_FACTORS:
            pubs.append(folded_with_meas[lam])
            pub_map.append((rep, lam))

    # Submit with retry
    job_metadata = []
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            t_submit = time.time()
            sampler = SamplerV2(mode=backend)
            job = sampler.run(pubs, shots=N_SHOTS)
            job_id = job.job_id()
            log.info(f"  Job submitted: {len(pubs)} PUBs, id={job_id} "
                     f"(attempt {attempt})")

            result = job.result()
            elapsed = time.time() - t_submit
            log.info(f"  Job done: {elapsed:.1f}s")

            job_metadata.append({
                "job_id": job_id,
                "n_pubs": len(pubs),
                "elapsed_s": round(elapsed, 1),
                "submit_time": datetime.now(timezone.utc).isoformat(),
            })
            break

        except Exception as e:
            log.warning(f"  Attempt {attempt}/{MAX_RETRIES} failed: {e}")
            if attempt == MAX_RETRIES:
                log.error(f"  ALL ATTEMPTS FAILED for TP {timepoint_idx}")
                raise
            delay = RETRY_BASE_DELAY_S * (2 ** (attempt - 1))
            log.info(f"  Retrying in {delay}s …")
            time.sleep(delay)

    # Extract expectations per (rep, lambda)
    per_lam_exps = {lam: [] for lam in SCALE_FACTORS}
    per_rep_data = {}

    for pub_idx, (rep_idx, lam) in enumerate(pub_map):
        pub_result = result[pub_idx]
        counts = pub_result.data.meas.get_counts()
        exp_val = expectation_from_counts(counts, N_QUBITS)
        per_lam_exps[lam].append(exp_val)

        if rep_idx not in per_rep_data:
            per_rep_data[rep_idx] = {}
        per_rep_data[rep_idx][lam] = {
            "exp_val": exp_val,
            "counts": counts,
        }

    # Compute raw and ZNE-mitigated arrays
    raw_exps = np.array(per_lam_exps[1.0])
    mit_exps = np.zeros(n_reps)
    for rep in range(n_reps):
        shot_vals = [per_lam_exps[lam][rep] for lam in SCALE_FACTORS]
        mit_exps[rep] = extrapolate(SCALE_FACTORS, shot_vals,
                                     EXTRAPOLATION_METHOD)

    # Statistical analysis
    stats = paired_analysis(raw_exps, mit_exps, ideal)

    # Verdict
    if stats["significant"] and stats["cohen_d"] > 0:
        verdict = "sig_better"
    elif stats["significant"] and stats["cohen_d"] < 0:
        verdict = "sig_worse"
    else:
        verdict = "not_sig"

    # Calibration snapshot
    cal_snapshot = get_backend_calibration_snapshot(backend)

    # ── Save raw data JSON ──
    json_path = results_dir / f"tp_{timepoint_idx:04d}_{ts_slug}.json"
    raw_data = {
        "metadata": {
            "timestamp": ts_iso,
            "timepoint_idx": timepoint_idx,
            "backend": backend.name,
            "config": {
                "tc": tc_name,
                "n_reps": n_reps,
                "n_shots": N_SHOTS,
                "scale_factors": SCALE_FACTORS,
                "folding": FOLDING_STRATEGY,
                "extrapolation": EXTRAPOLATION_METHOD,
            },
            "calibration": cal_snapshot,
            "jobs": job_metadata,
        },
        "per_lambda": {str(lam): per_lam_exps[lam]
                       for lam in SCALE_FACTORS},
        "raw_expectations": raw_exps.tolist(),
        "mitigated_expectations": mit_exps.tolist(),
        "statistics": stats,
        "verdict": verdict,
    }
    with open(json_path, "w") as f:
        json.dump(raw_data, f, indent=2, default=str)
    log.debug(f"  JSON → {json_path.name}")

    # ── Granular CSV (1 row per rep x λ) ──
    csv_path = results_dir / "raw_data.csv"
    for rep in range(n_reps):
        for lam in SCALE_FACTORS:
            _csv_append(csv_path, {
                "timepoint_idx": timepoint_idx,
                "timestamp": ts_iso,
                "backend": backend.name,
                "rep": rep,
                "scale_factor": lam,
                "exp_val": f"{per_lam_exps[lam][rep]:.6f}",
                "n_shots": N_SHOTS,
                "ideal": f"{ideal:.6f}",
            })

    # ── Summary CSV (1 row per TP) ──
    summary_row = {
        "timepoint_idx": timepoint_idx,
        "timestamp": ts_iso,
        "backend": backend.name,
        "n_reps": n_reps,
        "n_shots": N_SHOTS,
        "ideal": f"{ideal:.6f}",
        "mean_raw": f"{np.mean(raw_exps):.6f}",
        "std_raw": f"{np.std(raw_exps, ddof=1):.6f}",
        "mean_mit": f"{np.mean(mit_exps):.6f}",
        "std_mit": f"{np.std(mit_exps, ddof=1):.6f}",
        "mean_lam3": f"{np.mean(per_lam_exps[3.0]):.6f}",
        "mean_lam5": f"{np.mean(per_lam_exps[5.0]):.6f}",
        "cohen_d": f"{stats['cohen_d']:.4f}",
        "p_value": f"{stats['p_value_t']:.2e}",
        "significant": stats["significant"],
        "verdict": verdict,
        "sigma_ci": f"{sigma_ci(SCALE_FACTORS):.4f}",
        "calibration_date": cal_snapshot.get("calibration_date"),
        "job_elapsed_s": job_metadata[0]["elapsed_s"] if job_metadata else "",
    }
    _csv_append(results_dir / "drift_log.csv", summary_row)

    # ── Calibration log ──
    cal_path = results_dir / "calibration.json"
    if cal_path.exists():
        with open(cal_path) as f:
            cal_log = json.load(f)
    else:
        cal_log = []
    cal_log.append(cal_snapshot)
    with open(cal_path, "w") as f:
        json.dump(cal_log, f, indent=2, default=str)

    # Print summary
    sig_str = "***" if stats["significant"] and abs(stats["cohen_d"]) > 0.8 \
        else "* " if stats["significant"] else "  "
    log.info(f"  d={stats['cohen_d']:+.4f} {sig_str} "
             f"p={stats['p_value_t']:.2e}  verdict={verdict}")
    log.info(f"  E(λ=1)={np.mean(raw_exps):.4f}  "
             f"E(λ=3)={np.mean(per_lam_exps[3.0]):.4f}  "
             f"E(λ=5)={np.mean(per_lam_exps[5.0]):.4f}")

    return {
        "timepoint_idx": timepoint_idx,
        "timestamp": ts_iso,
        **stats,
        "verdict": verdict,
        "mean_raw": float(np.mean(raw_exps)),
        "mean_lam3": float(np.mean(per_lam_exps[3.0])),
        "mean_lam5": float(np.mean(per_lam_exps[5.0])),
    }


# ═══════════════════════════════════════════════════════════════════════
#  Single Time-Point: Local Simulator
# ═══════════════════════════════════════════════════════════════════════

def collect_timepoint_local(
    circuit_info: dict,
    n_reps: int,
    timepoint_idx: int,
    results_dir: Path,
    noise_drift_factor: float = 1.0,
) -> dict:
    """Simulate one TP with depolarising noise + drift."""
    from qiskit_aer import AerSimulator

    timestamp = datetime.now(timezone.utc)
    ts_iso = timestamp.isoformat()
    ideal = circuit_info["ideal"]
    folded = circuit_info["folded"]

    base_p2q = 0.003  # realistic for Heron
    p2q = base_p2q * noise_drift_factor
    p1q = p2q * 0.1
    nm = make_noise_model("depolarizing", p1q, p2q)

    log.info(f"═══ TP {timepoint_idx} (local) | "
             f"p2q={p2q:.4f}, drift={noise_drift_factor:.2f} ═══")

    sim = AerSimulator(method="density_matrix", noise_model=nm)

    # Get exact noisy expectations
    exact_exps = {}
    for lam in SCALE_FACTORS:
        fc = folded[lam].copy()
        fc.save_density_matrix()
        tc = transpile(fc, sim, optimization_level=0)
        result = sim.run(tc).result()
        dm = np.asarray(result.data()["density_matrix"])
        diag = np.real(np.diag(dm))
        signs = np.array([(-1) ** bin(x).count("1")
                         for x in range(1 << N_QUBITS)], dtype=float)
        exact_exps[lam] = float(signs @ diag)

    # Resample with shot noise
    rng = np.random.default_rng(SEED + timepoint_idx * 1000)
    per_lam_exps = {lam: [] for lam in SCALE_FACTORS}
    for rep in range(n_reps):
        for lam in SCALE_FACTORS:
            val = sample_shot_noise(exact_exps[lam], N_SHOTS, rng)
            per_lam_exps[lam].append(val)

    raw_exps = np.array(per_lam_exps[1.0])
    mit_exps = np.zeros(n_reps)
    for rep in range(n_reps):
        shot_vals = [per_lam_exps[lam][rep] for lam in SCALE_FACTORS]
        mit_exps[rep] = extrapolate(SCALE_FACTORS, shot_vals,
                                     EXTRAPOLATION_METHOD)

    stats = paired_analysis(raw_exps, mit_exps, ideal)

    if stats["significant"] and stats["cohen_d"] > 0:
        verdict = "sig_better"
    elif stats["significant"] and stats["cohen_d"] < 0:
        verdict = "sig_worse"
    else:
        verdict = "not_sig"

    # Save granular CSV
    csv_path = results_dir / "raw_data.csv"
    for rep in range(n_reps):
        for lam in SCALE_FACTORS:
            _csv_append(csv_path, {
                "timepoint_idx": timepoint_idx,
                "timestamp": ts_iso,
                "backend": f"local_depol(p2q={p2q:.4f})",
                "rep": rep,
                "scale_factor": lam,
                "exp_val": f"{per_lam_exps[lam][rep]:.6f}",
                "n_shots": N_SHOTS,
                "ideal": f"{ideal:.6f}",
            })

    # Save summary
    _csv_append(results_dir / "drift_log.csv", {
        "timepoint_idx": timepoint_idx,
        "timestamp": ts_iso,
        "backend": f"local_depol(p2q={p2q:.4f})",
        "n_reps": n_reps,
        "n_shots": N_SHOTS,
        "ideal": f"{ideal:.6f}",
        "mean_raw": f"{np.mean(raw_exps):.6f}",
        "std_raw": f"{np.std(raw_exps, ddof=1):.6f}",
        "mean_mit": f"{np.mean(mit_exps):.6f}",
        "std_mit": f"{np.std(mit_exps, ddof=1):.6f}",
        "mean_lam3": f"{np.mean(per_lam_exps[3.0]):.6f}",
        "mean_lam5": f"{np.mean(per_lam_exps[5.0]):.6f}",
        "cohen_d": f"{stats['cohen_d']:.4f}",
        "p_value": f"{stats['p_value_t']:.2e}",
        "significant": stats["significant"],
        "verdict": verdict,
        "noise_drift_factor": f"{noise_drift_factor:.4f}",
        "p_2q": f"{p2q:.6f}",
    })

    log.info(f"  d={stats['cohen_d']:+.4f}  verdict={verdict}  "
             f"E(λ=1)={np.mean(raw_exps):.4f}  "
             f"E(λ=3)={np.mean(per_lam_exps[3.0]):.4f}")

    return {
        "timepoint_idx": timepoint_idx,
        "timestamp": ts_iso,
        **stats,
        "verdict": verdict,
        "mean_raw": float(np.mean(raw_exps)),
    }


# ═══════════════════════════════════════════════════════════════════════
#  Mode: Single
# ═══════════════════════════════════════════════════════════════════════

def mode_single(args):
    """Run a single time-point (for testing / config check)."""
    results_dir = Path(args.outdir)
    results_dir.mkdir(parents=True, exist_ok=True)

    if args.local_test:
        log.info("▶ LOCAL TEST — single time-point")
        circuit_info = prepare_circuits(local_test=True)
        r = collect_timepoint_local(
            circuit_info, args.n_reps, 0, results_dir,
            noise_drift_factor=1.0,
        )
    elif args.dry_run:
        log.info("▶ DRY RUN — transpile only, no submission")
        backend, service = get_ibm_backend(args.backend)
        circuit_info = prepare_circuits(backend=backend)
        cal = get_backend_calibration_snapshot(backend)
        log.info(f"  Calibration date: {cal.get('calibration_date')}")
        log.info(f"  Gate errors: {cal.get('gate_errors')}")
        log.info(f"  Readout errors: {cal.get('readout_errors')}")
        log.info("  ✓ Dry run complete — circuits ready, no jobs submitted")
        return
    else:
        log.info(f"▶ SINGLE TP on {args.backend}")
        backend, service = get_ibm_backend(args.backend)
        circuit_info = prepare_circuits(backend=backend)
        r = collect_timepoint_hw(
            backend, circuit_info, args.n_reps, 0, results_dir,
        )

    log.info(f"  Result: d={r['cohen_d']:+.4f}  verdict={r['verdict']}")


# ═══════════════════════════════════════════════════════════════════════
#  Mode: Longitudinal
# ═══════════════════════════════════════════════════════════════════════

def mode_longitudinal(args):
    """Run time-points at regular intervals with crash-resume."""
    results_dir = Path(args.outdir)
    results_dir.mkdir(parents=True, exist_ok=True)

    interval_s = args.interval_minutes * 60
    duration_s = args.duration_hours * 3600
    n_tp = max(1, int(duration_s / interval_s)) + 1

    # Cost estimate
    pubs_per_tp = args.n_reps * len(SCALE_FACTORS)
    est_qpu_per_tp = pubs_per_tp * 1.14
    est_total_qpu = n_tp * est_qpu_per_tp

    log.info("╔═════════════════════════════════════════════════════════════╗")
    log.info(f"║  IBM DRIFT STUDY — {args.backend:<20s}               ║")
    log.info(f"║  {n_tp} TPs, every {args.interval_minutes:.0f} min, "
             f"{args.duration_hours:.1f}h total                        ║")
    log.info(f"║  {args.n_reps} reps x {len(SCALE_FACTORS)} λ x "
             f"{N_SHOTS} shots = {pubs_per_tp} PUBs/TP               ║")
    log.info(f"║  Est. QPU: {est_qpu_per_tp:.0f}s/TP, "
             f"{est_total_qpu:.0f}s total ({est_total_qpu/60:.1f} min)       ║")
    log.info("╚═════════════════════════════════════════════════════════════╝")

    # Crash-resume
    completed_tps = get_completed_timepoints(results_dir)
    if completed_tps:
        log.info(f"  ⚡ RESUMING: {len(completed_tps)} TPs already done, "
                 f"{n_tp - len(completed_tps)} remaining")

    # Prepare backend + circuits
    if args.local_test:
        log.info("  Mode: LOCAL SIMULATOR")
        circuit_info = prepare_circuits(local_test=True)
        backend = None
        # Drift factors: sinusoidal to mimic calibration cycle
        drift_factors = 1.0 + 0.5 * np.sin(
            np.linspace(0, 2 * np.pi, n_tp))
    else:
        backend, service = get_ibm_backend(args.backend)
        circuit_info = prepare_circuits(backend=backend)

    results = []
    t_start = time.time()
    n_collected = 0
    n_failed = 0
    consecutive_failures = 0

    for tp_idx in range(n_tp):
        # Skip completed (crash-resume)
        if tp_idx in completed_tps:
            log.info(f"  TP {tp_idx} already done, skipping.")
            continue

        elapsed_h = (time.time() - t_start) / 3600
        remaining = n_tp - tp_idx - len(completed_tps)

        log.info(f"\n{'─' * 60}")
        log.info(f"  TP {tp_idx + 1}/{n_tp} | "
                 f"{elapsed_h:.1f}h elapsed | "
                 f"{remaining} remaining")
        log.info(f"{'─' * 60}")

        try:
            if args.local_test:
                r = collect_timepoint_local(
                    circuit_info, args.n_reps, tp_idx, results_dir,
                    noise_drift_factor=drift_factors[tp_idx],
                )
            else:
                r = collect_timepoint_hw(
                    backend, circuit_info, args.n_reps, tp_idx, results_dir,
                )

            results.append(r)
            n_collected += 1
            consecutive_failures = 0

            # Heartbeat
            log.info(f"  ❤ collected={n_collected} failed={n_failed}")

        except Exception as e:
            n_failed += 1
            consecutive_failures += 1
            log.error(f"  ✗ TP {tp_idx} FAILED "
                      f"({n_failed} total, {consecutive_failures} consecutive)")
            log.debug(traceback.format_exc())
            results.append({"timepoint_idx": tp_idx, "error": str(e)})

            # Re-connect after 3 consecutive failures
            if consecutive_failures >= 3 and not args.local_test:
                log.warning("  ⚡ 3 consecutive failures — reconnecting …")
                try:
                    backend, service = get_ibm_backend(args.backend)
                    circuit_info = prepare_circuits(backend=backend)
                    consecutive_failures = 0
                    log.info("  ✓ Reconnected successfully")
                except Exception as re_e:
                    log.error(f"  Reconnection failed: {re_e}")

        # Sleep until next interval
        if tp_idx < n_tp - 1:
            next_time = t_start + (tp_idx + 1) * interval_s
            wait = max(0, next_time - time.time())
            if wait > 0 and not args.local_test:
                log.info(f"  Sleeping {wait/60:.1f} min until next TP …")
                time.sleep(wait)

    # Final summary
    _print_summary(results)
    return results


# ═══════════════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════════════

def _print_summary(results: list[dict]):
    """Print a summary table of the drift time-series."""
    log.info("\n" + "═" * 72)
    log.info("  IBM DRIFT EXPERIMENT — SUMMARY")
    log.info("═" * 72)
    log.info(f"  {'TP':>3} {'d':>8} {'p':>12} {'Verdict':>12} {'E(λ=1)':>8}")
    log.info("  " + "─" * 50)

    verdicts = []
    ds = []
    for r in results:
        if "error" in r:
            log.info(f"  {r['timepoint_idx']:>3} {'ERROR':>8}")
            continue
        v = r.get("verdict", "?")
        d = r.get("cohen_d", 0)
        p = r.get("p_value_t", 1)
        e1 = r.get("mean_raw", 0)
        verdicts.append(v)
        ds.append(d)
        log.info(f"  {r['timepoint_idx']:>3} {d:>+8.4f} {p:>12.2e} "
                 f"{v:>12} {e1:>8.4f}")

    if verdicts:
        from collections import Counter
        counts = Counter(verdicts)
        log.info(f"\n  Verdict distribution:")
        for v, c in counts.most_common():
            log.info(f"    {v}: {c} ({100*c/len(verdicts):.0f}%)")

        log.info(f"\n  Cohen's d range: [{min(ds):+.3f}, {max(ds):+.3f}]")
        log.info(f"  Cohen's d mean:  {np.mean(ds):+.3f} ± {np.std(ds):.3f}")

        n_unique = len(set(verdicts))
        if n_unique > 1:
            log.info(f"\n  ⚡ VERDICT INSTABILITY: "
                     f"{n_unique} distinct verdicts across "
                     f"{len(verdicts)} TPs")
        else:
            log.info(f"\n  ✓ Verdict stable across {len(verdicts)} TPs")

    log.info("═" * 72)


# ═══════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--mode", choices=["single", "longitudinal"],
        default="single",
        help="Execution mode (default: single)",
    )
    ap.add_argument(
        "--backend", type=str, default=BACKEND_NAME,
        help=f"IBM backend (default: {BACKEND_NAME}). "
             f"Alternatives: {', '.join(BACKEND_ALTERNATIVES)}",
    )
    ap.add_argument(
        "--n-reps", type=int, default=N_REPS_DEFAULT,
        help=f"Repetitions per time-point (default: {N_REPS_DEFAULT})",
    )
    ap.add_argument(
        "--interval-minutes", type=float,
        default=INTERVAL_MINUTES_DEFAULT,
        help=f"Minutes between TPs (default: {INTERVAL_MINUTES_DEFAULT})",
    )
    ap.add_argument(
        "--duration-hours", type=float,
        default=DURATION_HOURS_DEFAULT,
        help=f"Total duration in hours (default: {DURATION_HOURS_DEFAULT})",
    )
    ap.add_argument(
        "--local-test", action="store_true",
        help="Run on local simulator (no IBM connection needed)",
    )
    ap.add_argument(
        "--dry-run", action="store_true",
        help="Transpile circuits + print calibration, no submission",
    )
    ap.add_argument(
        "--outdir", type=str,
        default=str(RESULTS_DIR),
        help=f"Output directory (default: {RESULTS_DIR})",
    )
    args = ap.parse_args()

    log.info(f"exp_drift_ibm.py started | mode={args.mode} | "
             f"backend={args.backend}")

    if args.mode == "single":
        mode_single(args)
    elif args.mode == "longitudinal":
        mode_longitudinal(args)


if __name__ == "__main__":
    main()
