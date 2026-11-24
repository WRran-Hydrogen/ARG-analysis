#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --job-name=sbl_p4_hgt
#SBATCH --output=p4_hgt_%j.out
#
# HGT inference: pre-registered Fitch + F81 Marginal ML HGT inference,
# Sister-PRE / MAD / MinVar / Midpoint cross-rooting, OAT robustness scoring,
# and CL0013 superfamily LOO rooting validation.
# Submit from scripts/hpc/04_hgt_inference/:   sbatch 02_run_hgt_inference.sh

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"

# Use configurable Python environment (must have ete4, numpy, scipy, pandas)
PY_ETE="${SBL_HPC_PYTHON:-python}"
if ! "${PY_ETE}" -c "import ete4, numpy, scipy, pandas" 2>/dev/null; then
    echo "ERROR: ${PY_ETE} lacks ete4/numpy/scipy/pandas. Set SBL_HPC_PYTHON in config/hpc.env." >&2
    exit 1
fi

mkdir -p "${SBL_HPC_BASE}/04_hgt_inference"
cd "${SCRIPT_DIR}"
echo "[$(date -Is)] Running 01_hgt_inference.py with ${PY_ETE}"
"${PY_ETE}" -u 01_hgt_inference.py
echo "[$(date -Is)] Done"
