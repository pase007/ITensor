#!/usr/bin/env python3
"""Merge one-row sigma gap CSV files into a single sorted table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", nargs="?", default="results/raw")
    parser.add_argument("-o", "--output", default="results/sigma_gap_merged.csv")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    rows = []
    for path in sorted(input_dir.glob("sigma_gap_g*_N*.csv")):
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                row["source_file"] = path.name
                rows.append(row)

    rows.sort(key=lambda row: (float(row["g"]), int(row["N"])))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["g", "N", "E0", "E1", "gap", "var0", "var1", "source_file"]
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Merged {len(rows)} rows into {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
