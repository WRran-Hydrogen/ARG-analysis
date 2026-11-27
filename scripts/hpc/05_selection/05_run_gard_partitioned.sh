#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --job-name=sbl_p5_gardpart
#SBATCH --output=p5_gardpart_%A_%a.out
#
# Selection analysis, step 5 (sensitivity): MEME and FEL with one tree per GARD segment, for one line of
# gard_partition_tasks.tsv (written by: 05_gard_partitioned_selection.py tasks).
# Submit from scripts/hpc/05_selection/:
#   sbatch --array=1-<n_sets> 05_run_gard_partitioned.sh "$SBL_HPC_BASE/05_selection_pressure/gard_partition_tasks.tsv"

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
TASKS="$1"
HYPHY="${SBL_HPC_HYPHY:-hyphy}"
PY="${SBL_HPC_PYTHON:-python3}"
CPUS="${SLURM_CPUS_PER_TASK:-16}"

line=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n' "${TASKS}")
line="${line%$'\r'}"
[[ -n "${line}" ]] || { echo "No task line $((SLURM_ARRAY_TASK_ID + 1)) in ${TASKS}" >&2; exit 1; }
IFS=$'\t' read -r set_id family set_group track n_seqs n_codons med_ds iss iss_cs iss_ca set_dir <<< "${line}"

"${PY}" "${SCRIPT_DIR}/05_gard_partitioned_selection.py" build "${set_dir}"
cd "${set_dir}"
echo "[$(date -Is)] ${set_id} (${family}, n=${n_seqs}, cod=${n_codons})"
# TOLERATE_NUMERICAL_ERRORS: without it some partitioned fits stop with a ComputeBranchCache LF mismatch
# (differences in the fifth decimal of the log-likelihood), as GARD did. A failure of one method does not
# stop the other; the summary reports a set as 'missing' unless both JSON files are present.
rc=0
for method in meme fel; do
    if [[ ! -s "${method}_gard.json" ]]; then
        "${HYPHY}" CPU="${CPUS}" ENV="TOLERATE_NUMERICAL_ERRORS=1;" "${method}" --alignment gard_partitioned.nex \
            --output "${method}_gard.json" || { echo "${method} failed for ${set_id}" >&2; rm -f "${method}_gard.json"; rc=1; }
    fi
done
echo "[$(date -Is)] ${set_id}: done (rc=${rc})"
exit "${rc}"
