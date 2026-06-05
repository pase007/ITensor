#!/usr/bin/env python3
"""Generate a CSV parameter grid for Slurm array jobs."""

from __future__ import annotations

import argparse
import csv


def parse_values(text: str, cast):
    return [cast(value.strip()) for value in text.split(",") if value.strip()]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--g", required=True, help="Comma-separated g values, e.g. 0.7,0.8,0.9")
    parser.add_argument("--N", required=True, help="Comma-separated N values, e.g. 16,20,24,28")
    parser.add_argument("-o", "--output", default="hpc/gap_params.csv")
    args = parser.parse_args()

    g_values = parse_values(args.g, float)
    n_values = parse_values(args.N, int)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["g", "N"])
        for g in g_values:
            for n in n_values:
                writer.writerow([g, n])

    print(f"Wrote {len(g_values) * len(n_values)} jobs to {args.output}")
    print(f"Use: #SBATCH --array=1-{len(g_values) * len(n_values)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
