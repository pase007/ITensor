#!/bin/bash
#SBATCH --job-name=sigma_gap
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sigma_gap_%A_%a.out
#SBATCH --error=logs/sigma_gap_%A_%a.err

set -euo pipefail

# Replace these with the exact module names from `module avail` on LEO5.
module purge
module load gcc/...
module load openblas/...

# Submit from the dbqite_itensor directory:
#   mkdir -p logs results/raw
#   sbatch hpc/submit_gap_array.slurm
#
# Keep hpc/gap_params.csv synchronized with the SBATCH --array range.
# Array task 1 reads the first data row after the CSV header.
#
# Edit hpc/gap_params.csv by hand. It must have columns:
#   g,N0,N

PARAM_FILE="hpc/gap_params.csv"
OUT_DIR="results/raw"

mkdir -p logs "${OUT_DIR}"

line="$(awk -F',' -v task_id="${SLURM_ARRAY_TASK_ID}" 'NR == task_id + 1 {print $1 "," $2 "," $3}' "${PARAM_FILE}")"
if [[ -z "${line}" ]]; then
    echo "No parameter row for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
fi

IFS=',' read -r g N0 N <<< "${line}"

echo "Running sigma gap job for g=${g}, N=${N0}..${N}"
srun ./dbqite_itensor gap-range "${g}" "${N0}" "${N}" "${OUT_DIR}"
