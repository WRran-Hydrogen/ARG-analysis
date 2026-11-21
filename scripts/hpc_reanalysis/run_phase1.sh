#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --job-name=sbl_phase1
#SBATCH --output=phase1_%j.out
#
# Submit from scripts/hpc_reanalysis/:   sbatch run_phase1.sh
# Paths come from hpc.env (next to this script or in config/) or exported SBL_HPC_* variables.

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

for env_file in "${SCRIPT_DIR}/../../config/hpc.env" "${SCRIPT_DIR}/hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE (see config/hpc.env.example)}"

mkdir -p "${SBL_HPC_BASE}/01_curated_datasets"
cd "${SBL_HPC_BASE}/01_curated_datasets"
"${SBL_HPC_PYTHON:-python}" -u "${SCRIPT_DIR}/run_phase1.py"
