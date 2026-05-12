#!/usr/bin/env python3
"""
Generate desdentado_original_boxplot.csv from the published summary logs.

Reads the CSV files in data/desdentado_original/summary_*shots.csv and
reproduces results/desdentado_original_boxplot.csv exactly.

Usage:
    cd reproduction && python scripts/generate_desdentado_original_boxplot.py
"""

from __future__ import annotations

import argparse
import csv
import re
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path


SHOT_LABELS = {
    1024: "Default",
    3502: "Mid-Low",
    5980: "Estimated",
    8458: "Mid-High",
    10936: "High",
}

FIELDNAMES = ["shot_count", "label", "iteration", "mean_err", "success_rate"]
SIX_DECIMALS = Decimal("0.000001")
ONE = Decimal("1")
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
DEFAULT_INPUT_DIR = ROOT_DIR / "data" / "desdentado_original"
DEFAULT_OUTPUT_CSV = ROOT_DIR / "results" / "desdentado_original_boxplot.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rebuild results/desdentado_original_boxplot.csv from summary logs."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing summary_*shots.csv files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_CSV,
        help="Target CSV path to write.",
        required=True,
    )
    return parser.parse_args()


def shot_count_from_name(path: Path) -> int:
    match = re.fullmatch(r"summary_(\d+)shots\.csv", path.name)
    if match is None:
        raise ValueError(f"Unexpected summary filename: {path.name}")
    return int(match.group(1))


def format_decimal(value: Decimal) -> str:
    return format(value.quantize(SIX_DECIMALS, rounding=ROUND_HALF_UP), "f")


def iter_summary_rows(path: Path) -> list[dict[str, str]]:
    shot_count = shot_count_from_name(path)
    if shot_count not in SHOT_LABELS:
        raise ValueError(f"Missing label mapping for shot count {shot_count}")

    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        for row in reader:
            iteration = (row.get("Iteration") or "").strip()
            mean_err_raw = (row.get("Mean_err") or "").strip()
            if not iteration or not mean_err_raw:
                continue

            mean_err = Decimal(mean_err_raw)
            success_rate = ONE - mean_err
            rows.append({
                "shot_count": str(shot_count),
                "label": SHOT_LABELS[shot_count],
                "iteration": iteration,
                "mean_err": format_decimal(mean_err),
                "success_rate": format_decimal(success_rate),
            })

    return rows


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_path = Path(args.output)

    summary_files = sorted(
        input_dir.glob("summary_*shots.csv"),
        key=shot_count_from_name,
    )
    if not summary_files:
        raise FileNotFoundError(f"No summary files found in {input_dir}")

    rows: list[dict[str, str]] = []
    for summary_file in summary_files:
        rows.extend(iter_summary_rows(summary_file))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\r\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output_path}")


if __name__ == "__main__":
    main()
