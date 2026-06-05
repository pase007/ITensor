#!/bin/bash
#SBATCH --job-name=sigma_gap
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sigma_gap_%A_%a.out
#SBATCH --error=logs/sigma_gap_%A_%a.err

set -euo pipefail

# Submit from the dbqite_itensor directory:
#   mkdir -p logs results/raw
#   sbatch hpc/submit_gap_array.slurm
#
# Keep hpc/gap_params.csv synchronized with the SBATCH --array range.
# Array task 1 reads the first data row after the CSV header.

PARAM_FILE="hpc/gap_params.csv"
OUT_DIR="results/raw"

mkdir -p logs "${OUT_DIR}"

line="$(awk -F',' -v task_id="${SLURM_ARRAY_TASK_ID}" 'NR == task_id + 1 {print $1 "," $2}' "${PARAM_FILE}")"
if [[ -z "${line}" ]]; then
    echo "No parameter row for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
fi

g="${line%,*}"
N="${line#*,}"

echo "Running sigma gap job for g=${g}, N=${N}"
./dbqite_itensor single-gap "${g}" "${N}" "${OUT_DIR}"
