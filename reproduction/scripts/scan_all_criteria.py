#!/usr/bin/env python3
"""
scan_all_criteria.py — Regex-based scan of QEM papers for all 8 review criteria.

Scans 81 PDF papers for textual evidence of each criterion (C1-C8) in the
eight-criterion statistical rigor framework. Compares scan results against
existing ratings in review_criteria.csv and reports discrepancies.

Criteria:
  C1  Sample size reported           (CSV: sample_size)
  C2  Variance addressed             (CSV: variance)
  C3  Statistical evidence           (CSV: stat_tests)
  C4  Drift control                  (CSV: drift)
  C5  Overhead quantified            (CSV: overhead)
  C6  Noise model specified          (CSV: noise_model)
  C7  Reproduction package           (CSV: reproducibility)
  C8  Negative results reported      (CSV: neg_results)

For each criterion the scan uses two tiers of regex patterns:
  - STRONG patterns  → candidate 'yes'   (clear, specific evidence)
  - WEAK patterns    → candidate 'partial' (vague or indirect mentions)
  - No matches       → candidate 'no'

IMPORTANT: The regex scan provides CANDIDATE classifications only.
Final ratings require human/AI review. See REVIEW_METHODOLOGY.md.

Usage:
    python scan_all_criteria.py --pdf-dir <path> [options]

    # Scan all criteria, compare with CSV, produce logs:
    python scan_all_criteria.py --compare --log --dry-run

    # Scan only C3 and C7:
    python scan_all_criteria.py --criteria C3,C7 --compare --dry-run

    # Update CSV with scan results (careful!):
    python scan_all_criteria.py --update-csv
"""

import argparse
import csv
import json
import os
import re
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path

import pdfplumber

# ═══════════════════════════════════════════════════════════════════════════════
# Regex pattern definitions per criterion
# Each criterion has STRONG patterns (→ 'yes' candidate) and
# WEAK patterns (→ 'partial' candidate).
# ═══════════════════════════════════════════════════════════════════════════════

CRITERIA = {
    # ── C1: Sample Size Reported ──────────────────────────────────────────
    'C1': {
        'csv_column': 'sample_size',
        'label': 'Sample size reported',
        'strong': [
            # Explicit shot counts
            (r'\b\d{3,}\s*shots?\b', 'N-shots'),
            (r'\bshots?\s*=\s*\d{3,}', 'shots=N'),
            (r'\b(?:8192|4096|2048|1024|10000|20000|100000)\b', 'common-shot-count'),
            # Explicit sample/repetition counts
            (r'\b[Nn]\s*=\s*\d{2,}', 'N=value'),
            (r'\b\d+\s*(?:independent\s+)?(?:repetition|run|sample|realization)s?\b',
             'N-repetitions'),
            (r'\b(?:number|num)\s+of\s+(?:shots?|samples?|measurements?|circuits?)\s*'
             r'(?:=|:|\bis\b)\s*\d', 'num-of-shots'),
            (r'\b\d+\s*circuits?\s+(?:per|each|for)', 'circuits-per'),
        ],
        'weak': [
            (r'\bshot\s*(?:budget|count|number)\b', 'shot-budget-mention'),
            (r'\bsample\s*size\b', 'sample-size-mention'),
            (r'\brepetition', 'repetition-mention'),
            (r'\bmeasurement\s*(?:budget|count)', 'measurement-budget'),
            (r'\b(?:execute|run)\w*\s+\d+\s*times?\b', 'run-N-times'),
        ],
    },

    # ── C2: Variance Addressed ────────────────────────────────────────────
    'C2': {
        'csv_column': 'variance',
        'label': 'Variance addressed',
        'strong': [
            (r'\berror\s+bar', 'error-bar'),
            (r'\bstandard\s+deviation\b', 'std-dev'),
            (r'\bstandard\s+error\b', 'std-error'),
            (r'\bconfidence\s+interval', 'CI'),
            (r'±\s*\d', 'plus-minus-value'),
            (r'\bmean\s*±', 'mean-pm'),
            (r'\bbootstrap\b', 'bootstrap'),
            (r'\bjackknife\b', 'jackknife'),
            (r'\bσ\s*=\s*\d', 'sigma-value'),
            (r'\bvariance\b(?!\s*(?:of\s+)?(?:the\s+)?noise)', 'variance-reported'),
        ],
        'weak': [
            (r'\buncertainty\b', 'uncertainty'),
            (r'\bspread\b', 'spread'),
            (r'\bfluctuation', 'fluctuation'),
            (r'\bvariance\b', 'variance-general'),
            (r'\bnoise\s+floor\b', 'noise-floor'),
        ],
    },

    # ── C3: Statistical Evidence ──────────────────────────────────────────
    'C3': {
        'csv_column': 'stat_tests',
        'label': 'Statistical evidence',
        'strong': [
            # Hypothesis tests
            (r'\bp[\s\-]?value', 'p-value'),
            (r'\bp\s*[<>=≤≥]\s*0\.\d', 'p-threshold'),
            (r'hypothesis\s+test', 'hypothesis-test'),
            (r'significance\s+test', 'significance-test'),
            (r'statistically\s+significant', 'statistically-significant'),
            (r'reject\w*\s+(the\s+)?null', 'reject-null'),
            (r'\bt[\s\-]?test\b', 't-test'),
            (r'\bwilcoxon\b', 'Wilcoxon'),
            (r'\bmann[\s\-]?whitney\b', 'Mann-Whitney'),
            (r'\bchi[\s\-]?square\b', 'chi-square'),
            (r'\banova\b', 'ANOVA'),
            (r'\bf[\s\-]?test\b', 'F-test'),
            (r'permutation\s+test', 'permutation-test'),
            (r'kolmogorov[\s\-]?smirnov', 'Kolmogorov-Smirnov'),
            (r'\bkruskal[\s\-]?wallis\b', 'Kruskal-Wallis'),
            # Bootstrap for inference
            (r'bootstrap\w*\s+\w*\s*(?:signif|test|compar|p[\s\-]?value|confidence)',
             'bootstrap-inference'),
            # CI used for comparison
            (r'confidence\s+interval\w*\s+\w*\s*(?:compar|overlap|exclud|contain)',
             'CI-comparison'),
            (r'credible\s+interval', 'credible-interval'),
            (r'\bposterior\s+(?:distribution|probability|mean|median)', 'posterior'),
            (r'\bbayes\s*(?:ian)?\s+factor', 'Bayes-factor'),
            (r'bayesian\s+(?:inference|estimation|analysis|approach|framework|method)',
             'Bayesian-inference'),
            # Effect sizes
            (r"cohen['']?s?\s*d\b", "Cohens-d"),
            (r'\beffect\s+size\b', 'effect-size'),
        ],
        'weak': [
            (r'\berror\s+bar', 'error-bar'),
            (r'\bstandard\s+deviation\b', 'std-dev'),
            (r'\bstandard\s+error\b', 'std-error'),
            (r'\bconfidence\s+interval', 'CI'),
            (r'\buncertainty\b', 'uncertainty'),
            (r'\bvariance\b', 'variance'),
            (r'±', 'plus-minus'),
            (r'\bbootstrap\b', 'bootstrap'),
            (r'\bjackknife\b', 'jackknife'),
            (r'\bσ\b', 'sigma'),
        ],
    },

    # ── C4: Drift Control ─────────────────────────────────────────────────
    'C4': {
        'csv_column': 'drift',
        'label': 'Drift control',
        'strong': [
            (r'\bdrift\s+(?:correct|compensat|mitigat|account|characteriz|monitor)',
             'drift-action'),
            (r'\brecalibrat', 'recalibration'),
            (r'\binterleav\w+\s+(?:measure|circuit|experiment)', 'interleaved'),
            (r'\btemporal\s+(?:stability|variation|fluctuation)\s+\w*\s*(?:measure|monitor|track|account)',
             'temporal-tracking'),
            (r'\bhardware\s+(?:stability|calibration)\s+\w*\s*(?:check|verif|monitor)',
             'hw-stability-check'),
            (r'\bday[\s\-]to[\s\-]day\s+(?:variation|fluctuation|drift)',
             'day-to-day'),
            (r'\btime[\s\-]stamp', 'timestamp'),
        ],
        'weak': [
            (r'\bdrift\b', 'drift-mention'),
            (r'\btemporal\s+(?:fluctuation|variation|instability|noise)',
             'temporal-mention'),
            (r'\bcalibrat', 'calibration-mention'),
            (r'\bT[12]\s+(?:time|relaxation|decay)', 'T1-T2-mention'),
            (r'\bhardware\s+(?:noise|variation|instability)', 'hw-noise-mention'),
            (r'\bcoherence\s+time', 'coherence-time'),
            (r'\bnon[\s\-]?stationar', 'non-stationary'),
        ],
    },

    # ── C5: Overhead Quantified ───────────────────────────────────────────
    'C5': {
        'csv_column': 'overhead',
        'label': 'Overhead quantified',
        'strong': [
            (r'\boverhead\s*(?:=|:|\bis\b|of)\s*\d', 'overhead-value'),
            (r'\b\d+[×x]\s*(?:overhead|more|additional|extra)\s+(?:circuits?|samples?|shots?)',
             'Nx-overhead'),
            (r'\bsampling\s+(?:overhead|cost)\s*(?:=|:|\bis\b|of|scales?)',
             'sampling-cost'),
            (r'\bexponential\s+(?:overhead|cost|scaling)', 'exponential-overhead'),
            (r'\b(?:circuit|gate)\s+overhead', 'circuit-overhead'),
            (r'\bcomputational\s+(?:cost|overhead)\s*(?:=|:|\bis\b|of|scales?)',
             'computational-cost'),
            (r'\b(?:factor|ratio)\s+of\s+\d+', 'factor-of-N'),
            (r'\boverhead\s+(?:grow|scale|increas)', 'overhead-scaling'),
        ],
        'weak': [
            (r'\boverhead\b', 'overhead-mention'),
            (r'\b(?:cost|resource)\s+(?:of|for)\s+(?:mitigation|error)',
             'cost-mention'),
            (r'\bscaling\b', 'scaling-mention'),
            (r'\badditional\s+(?:circuits?|samples?|shots?|measurements?)',
             'additional-resources'),
            (r'\b(?:sample|shot)\s+(?:complexity|requirement)', 'sample-complexity'),
        ],
    },

    # ── C6: Noise Model Specified ─────────────────────────────────────────
    'C6': {
        'csv_column': 'noise_model',
        'label': 'Noise model specified',
        'strong': [
            (r'\bdepolariz\w+\s+(?:channel|noise|model|error)', 'depolarizing'),
            (r'\bpauli[\s\-]?(?:channel|noise|error|twirl)', 'Pauli-channel'),
            (r'\bamplitude\s+damping', 'amplitude-damping'),
            (r'\bphase\s+(?:damping|flip)', 'phase-damping'),
            (r'\bpauli[\s\-]?lindblad', 'Pauli-Lindblad'),
            (r'\bstochastic\s+(?:pauli|noise|error)', 'stochastic-noise'),
            (r'\bnoise\s+model\s*(?:=|:|\bis\b|used|based)', 'noise-model-spec'),
            (r'\btwirl\w*\s+(?:the|noise|into|channel)', 'twirling'),
            (r'\btensor[\s\-]?product\s+noise', 'tensor-product-noise'),
            (r'\blocal\s+(?:noise|depolarizing)', 'local-noise'),
            (r'\berror\s+rate\s*(?:=|:|\bof\b|is)\s*\d', 'error-rate-value'),
        ],
        'weak': [
            (r'\bnoise\s+model\b', 'noise-model-mention'),
            (r'\berror\s+(?:model|channel)', 'error-model-mention'),
            (r'\bnoisy\s+(?:simulation|backend|device)', 'noisy-device-mention'),
            (r'\bfake[\s_]?(?:backend|provider|device)', 'fake-backend'),
            (r'\bnoise\s+(?:level|strength|parameter|rate)', 'noise-parameter'),
            (r'\bgate\s+(?:error|fidelity|infidelity)', 'gate-error'),
        ],
    },

    # ── C7: Reproduction Package ──────────────────────────────────────────
    'C7': {
        'csv_column': 'reproducibility',
        'label': 'Reproduction package',
        'strong': [
            (r'github\.com/\S+', 'github-url'),
            (r'gitlab\.com/\S+', 'gitlab-url'),
            (r'\bzenodo\b', 'zenodo'),
            (r'\bdoi\.org/\S+', 'doi-link'),
            (r'\bcode\s+(?:is\s+)?(?:available|released|provided|open[\s\-]?source)',
             'code-available'),
            (r'\bdata\s+(?:is\s+)?(?:available|released|provided|open[\s\-]?source)',
             'data-available'),
            (r'\bopen[\s\-]?source(?:d)?\s+(?:code|software|implementation|package|library)',
             'open-source'),
            (r'\bpublicly\s+available\b', 'publicly-available'),
            (r'\bsupplementary\s+(?:material|code|data|information)',
             'supplementary'),
            (r'\bmitiq\b', 'mitiq'),
            (r'\bqiskit\b', 'qiskit'),
            (r'\bcirq\b', 'cirq'),
            (r'\bpennylane\b', 'pennylane'),
        ],
        'weak': [
            (r'\bavailable\s+(?:upon|on)\s+request', 'upon-request'),
            (r'\bcode\b.*\brepository\b', 'code-repository-mention'),
            (r'\breproduci', 'reproducibility-mention'),
            (r'\bsoftware\s+(?:package|framework|tool)', 'software-mention'),
        ],
    },

    # ── C8: Negative Results Reported ─────────────────────────────────────
    'C8': {
        'csv_column': 'neg_results',
        'label': 'Negative results reported',
        'strong': [
            (r'\b(?:mitigation|method|technique|approach)\s+(?:fails?|breaks?\s+down|does\s+not\s+(?:work|help|improve))',
             'method-fails'),
            (r'\bworse\s+(?:than|performance|result)', 'worse-than'),
            (r'\bdegradation\s+(?:of|in|at)', 'degradation'),
            (r'\bdiminishing\s+returns?\b', 'diminishing-returns'),
            (r'\b(?:no|without)\s+improvement\b', 'no-improvement'),
            (r'\binfeasible\b', 'infeasible'),
            (r'\bnegative\s+(?:result|finding|outcome)', 'negative-result'),
            (r'\blimitation\s+(?:of|is|include)', 'limitation-stated'),
            (r'\b(?:fails?|failure)\s+(?:for|at|when|in|beyond|above)',
             'failure-condition'),
            (r'\berror\s+(?:increase|amplif|grow)', 'error-increase'),
            (r'\bbreak\w*\s+down\s+(?:at|for|when|beyond)', 'breakdown'),
        ],
        'weak': [
            (r'\blimitation', 'limitation-mention'),
            (r'\bfailure\b', 'failure-mention'),
            (r'\bchalleng', 'challenge-mention'),
            (r'\bdisadvantage', 'disadvantage-mention'),
            (r'\bdrawback', 'drawback-mention'),
            (r'\btrade[\s\-]?off', 'tradeoff-mention'),
            (r'\bcaveat', 'caveat-mention'),
            (r'\bworse\b', 'worse-mention'),
        ],
    },
}

# Context window: characters before/after match to extract
CONTEXT_WINDOW = 150

# CSV column name → criterion key mapping
CSV_TO_CRITERION = {v['csv_column']: k for k, v in CRITERIA.items()}

# ═══════════════════════════════════════════════════════════════════════════════
# PDF text extraction and matching (reused from scan_statistical_evidence.py)
# ═══════════════════════════════════════════════════════════════════════════════


def extract_text_from_pdf(pdf_path: str) -> list[tuple[int, str]]:
    """Extract text from PDF, returning list of (page_number, text) tuples."""
    pages = []
    try:
        with pdfplumber.open(pdf_path) as pdf:
            for i, page in enumerate(pdf.pages):
                text = page.extract_text()
                if text:
                    pages.append((i + 1, text))
    except Exception as e:
        print(f"  WARNING: Could not read {pdf_path}: {e}", file=sys.stderr)
    return pages


def strip_accents(s: str) -> str:
    """Remove diacritical marks from string (e.g. ç→c, é→e)."""
    nfkd = unicodedata.normalize('NFKD', s)
    return ''.join(c for c in nfkd if not unicodedata.combining(c))


def normalize_paper_name(filename: str) -> str:
    """Extract author-year key from PDF filename for matching."""
    stem = Path(filename).stem
    parts = stem.split(' - ')
    if len(parts) >= 2:
        return strip_accents(f"{parts[0].strip()} - {parts[1].strip()}").lower()
    return strip_accents(stem).lower()


def match_pdf_to_csv_row(pdf_name: str, csv_papers: list[dict]) -> dict | None:
    """Match a PDF filename to a CSV row by author name and year.

    Handles disambiguated CSV entries like 'Kim et al. (stabilized)' by
    checking whether the parenthetical keyword appears in the PDF title.
    """
    norm = normalize_paper_name(pdf_name)
    pdf_parts = norm.split(' - ')
    if len(pdf_parts) < 2:
        return None
    pdf_author_part = pdf_parts[0].strip()
    pdf_year = pdf_parts[1].strip()
    # Get full title from raw filename (normalize_paper_name truncates it)
    raw_parts = strip_accents(Path(pdf_name).stem).lower().split(' - ')
    pdf_title = ' - '.join(raw_parts[2:]) if len(raw_parts) > 2 else ''

    pdf_first = re.split(r'\s+et\s+al|,|\s+and\s+', pdf_author_part)[0].strip()

    candidates = []
    for row in csv_papers:
        csv_name = strip_accents(row['paper'].lower())
        csv_year = str(row['year'])
        csv_first = re.split(r'\s+et\s+al|,|\s+and\s+|\s*&\s*|\(', csv_name)[0].strip()

        matched = False
        if pdf_first == csv_first and pdf_year == csv_year:
            matched = True
        elif pdf_first == csv_first and f"({pdf_year})" in csv_name:
            matched = True
        # Fallback: institution prefix before author in PDF filename
        elif pdf_year == csv_year and re.search(r'\b' + re.escape(csv_first) + r'\b', pdf_author_part):
            matched = True

        if matched:
            candidates.append(row)

    if len(candidates) == 1:
        return candidates[0]

    if len(candidates) > 1:
        # Use disambiguating keyword in parentheses to pick the right one
        for row in candidates:
            csv_name = strip_accents(row['paper'].lower())
            m = re.search(r'\(([^)]+)\)', csv_name)
            if m:
                keyword = m.group(1).strip()
                # Skip year-only parentheticals like (2023)
                if re.match(r'^\d{4}$', keyword):
                    continue
                if keyword in pdf_title or keyword in pdf_author_part:
                    return row
        # No keyword match — return first candidate as fallback
        return candidates[0]

    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Core scanning logic
# ═══════════════════════════════════════════════════════════════════════════════


def scan_paper_criterion(pages: list[tuple[int, str]],
                         criterion_key: str) -> dict:
    """
    Scan extracted PDF pages for a single criterion's patterns.

    Returns dict with:
      - 'strong_matches': list of (pattern_name, page, context_snippet)
      - 'weak_matches':   list of (pattern_name, page, context_snippet)
      - 'candidate':      'yes' | 'partial' | 'no'
    """
    crit = CRITERIA[criterion_key]
    strong_matches = []
    weak_matches = []

    for pattern, name in crit['strong']:
        for page_num, page_text in pages:
            for m in re.finditer(pattern, page_text, re.IGNORECASE):
                start = max(0, m.start() - CONTEXT_WINDOW)
                end = min(len(page_text), m.end() + CONTEXT_WINDOW)
                context = page_text[start:end].replace('\n', ' ').strip()
                strong_matches.append((name, page_num, context))

    for pattern, name in crit['weak']:
        for page_num, page_text in pages:
            for m in re.finditer(pattern, page_text, re.IGNORECASE):
                start = max(0, m.start() - CONTEXT_WINDOW)
                end = min(len(page_text), m.end() + CONTEXT_WINDOW)
                context = page_text[start:end].replace('\n', ' ').strip()
                weak_matches.append((name, page_num, context))

    if strong_matches:
        candidate = 'yes'
    elif weak_matches:
        candidate = 'partial'
    else:
        candidate = 'no'

    return {
        'strong_matches': strong_matches,
        'weak_matches': weak_matches,
        'candidate': candidate,
    }


def scan_paper_all(pdf_path: str, criteria_keys: list[str]) -> dict:
    """
    Scan a single PDF for all specified criteria.

    Returns dict mapping criterion_key -> scan result.
    Returns None if PDF is unreadable.
    """
    pages = extract_text_from_pdf(pdf_path)
    if not pages:
        return None

    results = {}
    for key in criteria_keys:
        results[key] = scan_paper_criterion(pages, key)
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Comparison and output
# ═══════════════════════════════════════════════════════════════════════════════


def compare_with_csv(scan_result: dict, csv_row: dict,
                     criteria_keys: list[str]) -> list[dict]:
    """
    Compare scan candidate ratings with existing CSV ratings.
    Returns list of discrepancy dicts.
    """
    discrepancies = []
    for key in criteria_keys:
        col = CRITERIA[key]['csv_column']
        csv_rating = csv_row.get(col, '')
        if csv_rating == 'na':
            continue  # skip n/a (theoretical papers)

        candidate = scan_result[key]['candidate']
        if candidate != csv_rating:
            discrepancies.append({
                'criterion': key,
                'csv_column': col,
                'csv_rating': csv_rating,
                'scan_candidate': candidate,
                'n_strong': len(scan_result[key]['strong_matches']),
                'n_weak': len(scan_result[key]['weak_matches']),
            })
    return discrepancies


def write_criterion_log(log_dir: str, criterion_key: str,
                        all_results: list[dict]):
    """Write detailed match log for a single criterion."""
    log_path = os.path.join(log_dir, f'scan_{criterion_key}_log.csv')
    with open(log_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'pdf_file', 'csv_paper', 'csv_rating', 'scan_candidate',
            'match_tier', 'pattern_name', 'page', 'context_snippet',
        ])
        for r in all_results:
            scan = r.get('scan', {}).get(criterion_key)
            if scan is None:
                writer.writerow([
                    r['pdf_file'], r['csv_paper'], r.get('csv_rating', {}).get(criterion_key, ''),
                    'UNREADABLE', '', '', '', '',
                ])
                continue

            csv_rating = r.get('csv_rating', {}).get(criterion_key, '')
            candidate = scan['candidate']

            if not scan['strong_matches'] and not scan['weak_matches']:
                writer.writerow([
                    r['pdf_file'], r['csv_paper'], csv_rating,
                    candidate, 'none', '', '', '',
                ])
            for name, page, ctx in scan['strong_matches']:
                writer.writerow([
                    r['pdf_file'], r['csv_paper'], csv_rating,
                    candidate, 'strong', name, page, ctx,
                ])
            for name, page, ctx in scan['weak_matches']:
                writer.writerow([
                    r['pdf_file'], r['csv_paper'], csv_rating,
                    candidate, 'weak', name, page, ctx,
                ])
    return log_path


def write_discrepancy_report(report_path: str, all_results: list[dict],
                             criteria_keys: list[str]):
    """Write summary report of all discrepancies between scan and CSV."""
    with open(report_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'paper', 'criterion', 'csv_column', 'csv_rating',
            'scan_candidate', 'direction', 'n_strong', 'n_weak',
        ])
        for r in all_results:
            for d in r.get('discrepancies', []):
                # Direction: upgrade (scan > csv), downgrade (scan < csv), mismatch
                rating_order = {'no': 0, 'partial': 1, 'yes': 2}
                csv_ord = rating_order.get(d['csv_rating'], -1)
                scan_ord = rating_order.get(d['scan_candidate'], -1)
                if scan_ord > csv_ord:
                    direction = 'UPGRADE'
                elif scan_ord < csv_ord:
                    direction = 'DOWNGRADE'
                else:
                    direction = 'MISMATCH'

                writer.writerow([
                    r['csv_paper'], d['criterion'], d['csv_column'],
                    d['csv_rating'], d['scan_candidate'], direction,
                    d['n_strong'], d['n_weak'],
                ])


def print_summary(all_results: list[dict], criteria_keys: list[str]):
    """Print summary statistics."""
    print("\n" + "=" * 78)
    print("SCAN SUMMARY")
    print("=" * 78)

    matched = [r for r in all_results if r['csv_paper'] != 'UNMATCHED']
    unmatched = [r for r in all_results if r['csv_paper'] == 'UNMATCHED']
    unreadable = [r for r in all_results if r.get('scan') is None]

    print(f"\n  PDFs processed:     {len(all_results)}")
    print(f"  Matched to CSV:     {len(matched)}")
    print(f"  Unmatched PDFs:     {len(unmatched)}")
    print(f"  Unreadable PDFs:    {len(unreadable)}")
    if unmatched:
        for r in unmatched:
            print(f"    - {r['pdf_file']}")

    for key in criteria_keys:
        print(f"\n  ── {key}: {CRITERIA[key]['label']} ──")
        candidates = defaultdict(int)
        csv_ratings = defaultdict(int)
        discrepancy_count = 0

        for r in matched:
            scan = r.get('scan', {}).get(key)
            if scan is None:
                continue
            candidates[scan['candidate']] += 1
            csv_col = CRITERIA[key]['csv_column']
            csv_val = r.get('csv_row', {}).get(csv_col, '')
            csv_ratings[csv_val] += 1
            # Count discrepancies (excluding na)
            if csv_val != 'na' and scan['candidate'] != csv_val:
                discrepancy_count += 1

        print(f"    Scan candidates:  yes={candidates['yes']}"
              f"  partial={candidates['partial']}"
              f"  no={candidates['no']}")
        print(f"    CSV ratings:      yes={csv_ratings.get('yes', 0)}"
              f"  partial={csv_ratings.get('partial', 0)}"
              f"  no={csv_ratings.get('no', 0)}"
              f"  na={csv_ratings.get('na', 0)}")
        print(f"    Discrepancies:    {discrepancy_count}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(
        description='Scan QEM papers for all 8 review criteria')
    parser.add_argument(
        '--pdf-dir',
        default=os.path.join(
            os.path.dirname(__file__), '..', 'data', 'review_paper'),
        help='Directory containing PDF papers')
    parser.add_argument(
        '--csv',
        default=os.path.join(
            os.path.dirname(__file__), '..', 'data', 'review_criteria_llm.csv'),
        help='Path to review_criteria_llm.csv (LLM reference ratings)')
    parser.add_argument(
        '--log-dir',
        default=os.path.join(
            os.path.dirname(__file__), '..', 'data', 'scan_logs'),
        help='Directory for per-criterion log files')
    parser.add_argument(
        '--criteria',
        default='all',
        help='Comma-separated list of criteria to scan (e.g. C1,C3,C7). Default: all')
    parser.add_argument(
        '--compare', action='store_true',
        help='Compare scan results with existing CSV ratings')
    parser.add_argument(
        '--log', action='store_true',
        help='Write detailed per-criterion match logs')
    parser.add_argument(
        '--update-csv', action='store_true',
        help='Update CSV with scan candidate ratings (use with caution!)')
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Scan and report only; do not modify any files')
    args = parser.parse_args()

    pdf_dir = Path(args.pdf_dir).resolve()
    csv_path = Path(args.csv).resolve()
    log_dir = Path(args.log_dir).resolve()

    if not pdf_dir.is_dir():
        print(f"ERROR: PDF directory not found: {pdf_dir}", file=sys.stderr)
        sys.exit(1)
    if not csv_path.is_file():
        print(f"ERROR: CSV not found: {csv_path}", file=sys.stderr)
        sys.exit(1)

    # Parse criteria selection
    if args.criteria == 'all':
        criteria_keys = list(CRITERIA.keys())
    else:
        criteria_keys = [c.strip().upper() for c in args.criteria.split(',')]
        for k in criteria_keys:
            if k not in CRITERIA:
                print(f"ERROR: Unknown criterion '{k}'. "
                      f"Valid: {', '.join(CRITERIA.keys())}", file=sys.stderr)
                sys.exit(1)

    # Load CSV
    with open(csv_path, 'r', encoding='utf-8') as f:
        csv_papers = list(csv.DictReader(f))

    # Find PDFs
    pdfs = sorted(pdf_dir.glob('*.pdf'))
    print(f"Found {len(pdfs)} PDFs in {pdf_dir}")
    print(f"CSV has {len(csv_papers)} papers")
    print(f"Scanning criteria: {', '.join(criteria_keys)}")
    print()

    # Scan all PDFs
    all_results = []
    for pdf_path in pdfs:
        print(f"Scanning: {pdf_path.name[:80]}...")
        scan = scan_paper_all(str(pdf_path), criteria_keys)

        # Match to CSV
        matched_row = match_pdf_to_csv_row(pdf_path.name, csv_papers)
        csv_paper = matched_row['paper'] if matched_row else 'UNMATCHED'

        result = {
            'pdf_file': pdf_path.name,
            'csv_paper': csv_paper,
            'csv_row': matched_row or {},
            'scan': scan,
            'csv_rating': {},
            'discrepancies': [],
        }

        # Collect CSV ratings for matched papers
        if matched_row and scan is not None:
            for key in criteria_keys:
                col = CRITERIA[key]['csv_column']
                result['csv_rating'][key] = matched_row.get(col, '')

            # Compare if requested
            if args.compare:
                result['discrepancies'] = compare_with_csv(
                    scan, matched_row, criteria_keys)

        # Print inline status
        if scan is None:
            print(f"  -> {csv_paper}: UNREADABLE")
        else:
            parts = []
            for key in criteria_keys:
                s = scan[key]
                parts.append(f"{key}={s['candidate']}({len(s['strong_matches'])}s/{len(s['weak_matches'])}w)")
            print(f"  -> {csv_paper}: {' '.join(parts)}")
            if result['discrepancies']:
                for d in result['discrepancies']:
                    print(f"     *** DISCREPANCY {d['criterion']}: "
                          f"CSV={d['csv_rating']} vs scan={d['scan_candidate']}")

        all_results.append(result)

    # Print summary
    print_summary(all_results, criteria_keys)

    # Write logs
    if args.log:
        log_dir.mkdir(parents=True, exist_ok=True)
        for key in criteria_keys:
            lp = write_criterion_log(str(log_dir), key, all_results)
            print(f"\n  Log written: {lp}")

    # Write discrepancy report
    if args.compare:
        report_path = log_dir / 'discrepancy_report.csv' if args.log else Path('/tmp/discrepancy_report.csv')
        if args.log:
            log_dir.mkdir(parents=True, exist_ok=True)
        write_discrepancy_report(str(report_path), all_results, criteria_keys)
        print(f"\n  Discrepancy report: {report_path}")

        # Count total discrepancies
        total_disc = sum(len(r['discrepancies']) for r in all_results)
        print(f"  Total discrepancies: {total_disc}")

    # Update CSV
    if args.update_csv and not args.dry_run:
        # Read, update, write CSV
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            rows = list(reader)

        updated = 0
        for r in all_results:
            if r['csv_paper'] == 'UNMATCHED' or r['scan'] is None:
                continue
            for row in rows:
                if row['paper'] == r['csv_paper']:
                    for key in criteria_keys:
                        col = CRITERIA[key]['csv_column']
                        if row[col] == 'na':
                            continue
                        new_val = r['scan'][key]['candidate']
                        if row[col] != new_val:
                            row[col] = new_val
                            updated += 1
                    break

        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\n  CSV updated: {updated} ratings changed")
    elif args.update_csv:
        print("\n  DRY RUN: CSV not modified")


if __name__ == '__main__':
    main()
