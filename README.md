# Reproduction Package: Claim Against Measurement

This repository contains all scripts, data, and plotting code to reproduce
the experiments and figures presented in the paper:

> **Claim against Measurement:
> Statistical Artefacts in Quantum Error Mitigation Benchmarks** (QCE 2026)

## Quick Start

### Docker (recommended)

Build the container, run the full reproduction pipeline, and compile the paper PDF:

```bash
make repro_docker
```

This is equivalent to running the pipeline inside an isolated Python + R environment
and then compiling the LaTeX source. The final PDF is written to `build/paper/`.

### Local (venv + R)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r reproduction/requirements.txt   # Python: qiskit, numpy, scipy, pandas
Rscript -e "install.packages(c('tidyverse', 'scales', 'patchwork', 'tikzDevice', 'ggnewscale'))"
make repro
```

### flake.nix (Experimental)

```bash
$ nix develop
# make repro
```

### Development shell

```bash
make dev   # drops into the Docker dev container with all dependencies available
```

## Structure

```
claim_against_measurement/
├── Makefile
├── docker/                      # Docker image and compose file
├── paper/                       # LaTeX source, bibliography, precompiled plots
└── reproduction/
    ├── scripts/                 # Python experiment scripts
    │   ├── khan_backend.py                 # Khan et al.: full parameter-space sweep
    │   ├── generate_desdentado_original_boxplot.py
    │   ├── multiple_comparisons.py         # Bonferroni/Holm/BH correction (stdout only)
    │   ├── scan_all_criteria.py            # LLM scan of 81 PDFs for criteria C1-C8
    │   ├── full_statistics.py              # Drift: paired t/Wilcoxon, Cohen's d, Cliff's δ, Durbin-Watson
    │   └── power_analysis.py               # Drift: ICC, effective N, power analysis (Table drift-summary)
    ├── hardware/                # Real-hardware collection scripts (not in Makefile)
    │   ├── drift_qexa_week.py              # QExa longitudinal drift (requires MQSS_TOKEN)
    │   └── drift_ibm.py                    # IBM longitudinal drift (requires IBM_QUANTUM_TOKEN)
    ├── R/                       # R plotting scripts (called by make)
    ├── core/                    # Shared Python library (circuits, noise, stats, zne)
    ├── data/                    # Input datasets (git-tracked)
    │   ├── review_criteria.csv             # Primary: human consensus ratings (C1-C8)
    │   ├── review_criteria_llm.csv         # Reference: LLM automated scan ratings
    │   ├── scan_logs/                      # Per-criterion LLM scan audit trail
    │   ├── desdentado_original/            # Raw data from Desdentado et al.
    │   ├── qexa_drift/                     # QExA drift experiment measurements
    │   └── ibm_drift_results/              # IBM Brussels drift experiment data
    ├── REVIEW_METHODOLOGY.md    # Documentation of the review pipeline
    ├── literature_review.typ    # Typst: per-paper criterion annotations
    └── requirements.txt

build/                           # Generated output — not git-tracked, created by make
├── results/                     # Computed CSV files (input to R plots)
├── plots/                       # Generated TikZ + PDF plots
└── paper/                       # Compiled LaTeX PDF
```

## Make Targets

| Target              | Description                                                |
|---------------------|------------------------------------------------------------|
| `make`              | Compile the paper PDF (assumes plots exist in `build/`)    |
| `make repro`        | Run full pipeline locally: all experiments + plots + paper |
| `make repro_docker` | Same as above, but run inside Docker                       |
| `make dev`          | Start a development environment inside Docker              |
| `make clean`        | Remove `build/`                                            |

## Drift Statistics (`full_statistics.py`, `power_analysis.py`)

`make repro` runs both scripts once per drift session (`first_run`, `day2_full`,
`weekend` — the three sessions in `tab:drift-summary`), writing per-session logs
and CSVs to `build/results/drift_stats/<session>/`:

- `drift_stats_per_tp.csv`, `drift_stats_summary.csv`, `drift_stats_report.txt` —
  paired t-test / Wilcoxon signed-rank, Cohen's d, Cliff's delta, Durbin-Watson
  (`full_statistics.py`).
- `drift_power.csv`, `drift_power_summary.txt` — ICC (one-way), lag-autocorrelation,
  and Durbin-Watson (`power_analysis.py`).

These underpin the ICC, r1, and Cohen's d range figures reported for
each session in `tab:drift-summary`. Both scripts assume the `qexa_drift` CSV
schema (`timepoint_idx, timestamp, backend, rep, scale_factor, exp_val, n_shots, ideal`);
they are **not** compatible with `ibm_drift_results/raw_data_ibm_drift.csv`, which
uses a different, more minimal schema (no `ideal`/`backend`/`n_shots` columns) and
is not part of `tab:drift-summary`. Both scripts can also be run standalone against
any `qexa_drift/raw_data_*.csv` file, e.g.:

```bash
cd reproduction
python scripts/full_statistics.py --csv data/qexa_drift/raw_data_weekend.csv
python scripts/power_analysis.py --csv data/qexa_drift/raw_data_weekend.csv
```


## Backends

The quantum experiment script `khan_backend.py` supports three backends via `--backend`:

| Backend     | Description                                    |
|-------------|------------------------------------------------|
| `simulator` | Qiskit Aer with ideal depolarizing noise model |
| `fake`      | IBM FakeBackend (realistic calibration data)   |
| `ibm`       | Real IBM Quantum hardware via Qiskit Runtime   |

For `ibm`, set the environment variable `IBM_QUANTUM_TOKEN` or pass `--token`.

## Hardware Experiments

The scripts used to collect the real-hardware data included in the paper are
provided under `reproduction/hardware/` for transparency and reproducibility.
They are **not** invoked by the Makefile — the resulting data is already
checked into `reproduction/data/`.

| Script                  | Experiment                                                | Environment         |
|-------------------------|-----------------------------------------------------------|---------------------|
| `drift_qexa_weekend.py` | 48 h+ longitudinal ZNE drift on IQM QExa (LRZ)            | `MQSS_TOKEN`        |
| `drift_ibm.py`          | 12 h longitudinal ZNE drift on IBM Quantum (ibm\_brussel) | `IBM_QUANTUM_TOKEN` |

### QExa week drift (`drift_qexa_weekend.py`)

Requires access to the IQM QExa system at LRZ via the MQSS adapter. Set the
token before running:

```bash
export MQSS_TOKEN="<your-token>"
cd reproduction
python hardware/drift_qexa_weekend.py --duration-hours 48 --backend EQE1

# Local simulator smoke test (no token needed):
python hardware/drift_qexa_weekend.py --local-test --duration-hours 0.1
```

Additionally requires the `mqss` Python package, which is not in the main requirements file
since it's only needed for this specific experiment.

### IBM drift (`drift_ibm.py`)

Requires an IBM Quantum account. Set your API token as an environment
variable or in a `.env` file in the repository root:

```bash
# Option A — environment variable:
export IBM_QUANTUM_TOKEN="<your-token>"

# Option B — .env file in repo root:
echo 'IBM_QUANTUM_TOKEN=<your-token>' > .env
```

```bash
cd reproduction

# Single time-point (dry run — transpile only, no submission):
python hardware/drift_ibm.py --mode single --dry-run

# Full longitudinal run:
python hardware/drift_ibm.py --mode longitudinal --duration-hours 12

# Local simulator test (no token needed):
python hardware/drift_ibm.py --mode single --local-test
```
