# Review Methodology

This document describes the methodology and reproduction steps for the
**systematic literature review** component of the paper *"Claim against
Measurement: Statistical Artefacts in Quantum Error Mitigation
Benchmarks"* (QCE 2026).

For instructions on running the quantum experiments and building the paper
PDF, see the top-level [README](../README.md).

## Prerequisites

Only Python with `pdfplumber` is needed to reproduce the automated
scan. The human consensus ratings are committed directly and require no
additional tooling to use.

| Dependency | Version | Purpose |
|------------|---------|---------|
| Python | 3.11+ | PDF scanning, automated rating |
| `pdfplumber` | (Python) | PDF text extraction |

```bash
python -m venv .venv && source .venv/bin/activate
pip install pdfplumber
```

## Input

| Input | Path | Description |
|-------|------|-------------|
| PDF corpus | `data/review_paper/` | 81 QEM publications (2022–2026) |
| Human consensus ratings | `data/review_criteria.csv` | Primary dataset — committed, no rerun needed |
| LLM reference ratings | `data/review_criteria_llm.csv` | Automated scan output — committed |

The PDF corpus is **not** included in the repository due to copyright. To
rerun the scan, place the 81 PDFs in `reproduction/data/review_paper/`. The
bibliography (`review.bib`) is already committed there. All
downstream analysis uses only the committed CSV files and does not
require the PDFs.

---

## Literature Search

The 81-paper corpus was collected via Google Scholar, arXiv, IEEE Xplore, and
forward/backward citation tracking from the Cai et al. (2023) QEM review. The
following six query strings were used across all databases. Papers were included
if they describe an experimental evaluation of a QEM technique on real or simulated
quantum hardware (2022–2026). We did not exhaust every possible query string but aimed 
to cover a broad range of relevant literature.

**Query 1 — Statistical reporting in QEM**
```
("error mitigation" OR ZNE OR PEC) AND
(variance OR "confidence interval" OR "statistical significance" OR bias)
```

**Query 2 — Negative results and limitations**
```
("quantum error mitigation" OR "zero noise extrapolation" OR "probabilistic error cancellation")
AND (fail* OR limit* OR "does not work" OR "breaks down" OR insufficient OR ineffective)
```

**Query 3 — Temporal drift and calibration instability**
```
("quantum computer" OR "superconducting qubit" OR "trapped ion")
AND ("drift" OR "temporal" OR "time-varying" OR "noise fluctuation"
     OR "calibration drift" OR "instability")
```

**Query 4 — Reproducibility and cross-platform benchmarking**
```
("quantum computing" OR "quantum algorithm")
AND (reproducib* OR replicat* OR "cross-platform" OR "benchmark comparison"
     OR "inter-device")
```

**Query 5 — Uncertainty quantification and statistical inference**
```
("quantum error mitigation" OR "quantum computing")
AND ("bootstrap" OR "confidence interval" OR "uncertainty quantification"
     OR "statistical inference" OR "hypothesis test")
```

**Query 6 — Replication crisis and benchmark variability**
```
("quantum computing" OR "quantum algorithm")
AND ("reproducibility" OR "replication crisis" OR "benchmark variability"
     OR "failed to reproduce")
```

---

## Review Pipeline

```
reproduction/data/review_paper/ (81 PDFs + review.bib)
        │
        ├─ scan_all_criteria.py ──────────────────────────────┐
        │  Regex scan + LLM filter → reference ratings        │
        │  Output: data/review_criteria_llm.csv               │
        └──────────────────────────────────────────────────────┘
        │
        ├─ Human consensus review ────────────────────────────┐
        │  2 authors + 2×15 independent subset reviewers      │
        │  Output: data/review_criteria.csv  (primary)        │
        └──────────────────────────────────────────────────────┘
        │
        ▼
  data/review_criteria.csv (81 papers × 8 criteria)
        │
        ▼
  R/plot_review_compliance.R  →  build/plots/review_compliance.tex
```

---

## Step 1: Automated PDF Scan

Scan all 81 PDFs for textual evidence of the eight statistical-rigour criteria.

```bash
cd reproduction
../.venv/bin/python scripts/scan_all_criteria.py \
    --pdf-dir data/review_paper/ \
    --compare --log --dry-run
```

**What it does:**
- Extracts text from each PDF using `pdfplumber`
- Matches STRONG patterns (→ `yes` candidate) and WEAK patterns (→ `partial`
  candidate) for each of the eight criteria (C1–C8)
- Compares candidates against the committed LLM reference ratings in `data/review_criteria_llm.csv`
- Writes per-criterion match logs to `data/scan_logs/scan_C{1..8}_log.csv`
- Writes a discrepancy report to `data/scan_logs/discrepancy_report.csv`

The `--dry-run` flag prevents writing back to the CSV. Remove it to update
`review_criteria_llm.csv` in place.

**Output files:**

| File | Content |
|------|---------|
| `data/scan_logs/scan_C1_log.csv` … `scan_C8_log.csv` | Every regex match with 150-char context snippet, page number, pattern tag |
| `data/scan_logs/discrepancy_report.csv` | All cases where scan ≠ LLM CSV: paper, criterion, old rating, scan rating, direction |

### Eight Criteria and Pattern Tiers

| ID | Criterion | CSV column | STRONG patterns | WEAK patterns |
|----|-----------|------------|-----------------|---------------|
| C1 | Sample size reported | `sample_size` | 7 (shot counts, N=value) | 5 (budget, repetition) |
| C2 | Variance addressed | `variance` | 10 (error bar, SD, CI, ±) | 5 (uncertainty, fluctuation) |
| C3 | Statistical evidence | `stat_tests` | 23 (p-value, hypothesis test, Bayesian, bootstrap) | 10 (error bar, SD, CI) |
| C4 | Drift control | `drift` | 7 (drift-action, recalibration, interleaved) | 7 (drift mention, calibration) |
| C5 | Overhead quantified | `overhead` | 8 (overhead=value, sampling cost) | 5 (overhead mention, scaling) |
| C6 | Noise model specified | `noise_model` | 11 (depolarizing, Pauli-channel, Pauli-Lindblad) | 6 (noise model mention) |
| C7 | Reproduction package | `reproducibility` | 13 (github URL, zenodo, code available) | 4 (upon request) |
| C8 | Negative results | `neg_results` | 11 (method fails, worse than, degradation) | 8 (limitation, failure, challenge) |

### Rating Scale

- **`yes`:** Criterion fully addressed with specific, quantitative detail.
- **`partial`:** Criterion mentioned but incomplete or only descriptive.
- **`no`:** Not addressed at all.
- **`na`:** Not applicable for the paper's scope.

C3 uses an extended definition: `yes` = inferential statistics (hypothesis
tests, Bayesian posterior analysis, bootstrap for comparison, CI-overlap for
conclusions, scaling analysis); `partial` = descriptive only (error bars,
SDs, CIs shown but not used for inference); `no` = no statistical evidence.

---

## Step 2: Human Consensus Review

The primary dataset is `data/review_criteria.csv`, produced by manual review.
This file is committed and no rerun is required.

Two authors jointly reviewed all 81 papers against the criteria. Two additional
raters each independently rated a (different) subset of 15 papers. Disagreements
were resolved by consensus discussion. The final ratings are committed as
`data/review_criteria.csv`.

The automated LLM review (`data/review_criteria_llm.csv`, generated in Step 1)
achieved **77% raw agreement** with the human consensus across all applicable
paper-criterion pairs. The complete per-paper LLM evidence report is provided as
`literature_review.typ` (compiled PDF in the reproduction package).

**Data files:**

| File | Content |
|------|---------| 
| `data/review_criteria.csv` | Human consensus ratings — **primary dataset** |
| `data/review_criteria_llm.csv` | Automated LLM reference ratings (from Step 1) |

**Verification:**

```bash
# Row count (expect 82 = 81 papers + 1 header)
wc -l data/review_criteria.csv
```

### Expected Final Distribution (Consensus — full scope, n=81)

| Criterion | yes | partial | no | na | Applicable | yes% |
|-----------|-----|---------|----|----|------------|------|
| C1 (Sample Size) | 47 | 10 | 4 | 7 | 61 | 77% |
| C2 (Variance) | 48 | 10 | 4 | 6 | 62 | 77% |
| C3 (Stat. Evidence) | 15 | 25 | 19 | 9 | 59 | 25% |
| C4 (Drift Control) | 12 | 16 | 12 | 28 | 40 | 30% |
| C5 (Overhead) | 43 | 19 | 3 | 3 | 65 | 66% |
| C6 (Noise Model) | 49 | 13 | 1 | 5 | 63 | 78% |
| C7 (Reproducibility) | 37 | 27 | 3 | 1 | 67 | 55% |
| C8 (Neg. Results) | 45 | 23 | 0 | 0 | 68 | 66% |

---

## Prompt for Generating the Scan Script

If you want to recreate `scan_all_criteria.py` from scratch using an AI
assistant, use the following prompt.

```text
Write a Python script that scans a directory of PDF papers for textual
evidence of research-quality criteria. The script should:

1. Accept --pdf-dir (path to PDFs), --csv (path to existing ratings CSV),
   --compare, --log, --dry-run flags
2. For each criterion, define two tiers of regex patterns:
   - STRONG: specific evidence (maps to "yes" candidate)
   - WEAK: vague/indirect mentions (maps to "partial" candidate)
3. Extract text from each PDF page using pdfplumber
4. Match all patterns, recording: paper, criterion, page, pattern tag,
   150-char context snippet before and after each match
5. Compare scan candidates against existing CSV ratings
6. Write per-criterion log CSVs and a discrepancy report

Here are the 8 criteria with their CSV column names:
[paste the criterion table]

Here are example STRONG/WEAK patterns for each criterion:
[paste pattern examples or describe what constitutes evidence]
```
