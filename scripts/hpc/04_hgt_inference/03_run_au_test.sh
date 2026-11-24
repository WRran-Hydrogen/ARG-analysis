#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=sbl_p4_au
#SBATCH --output=p4_au_%A_%a.out
#
# AU test of one HGT event (one line of au_test_tasks.tsv, written by 01_hgt_inference.py).
#   1. Constrained ML search: the donor background (all tips below the grandparent except the focal
#      clade) must form a clade that excludes the focal tips, i.e. the focal tips are not nested.
#   2. AU test (RELL, 10,000 replicates) of the unconstrained ML tree against the constrained tree.
# A small p-AU for the constrained tree rejects the "not nested" alternative.
#   usage (from scripts/hpc/04_hgt_inference/): sbatch --array=1-N 03_run_au_test.sh <au_test_tasks.tsv>

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
TASKS="$1"
IQTREE="${SBL_HPC_IQTREE:-iqtree2}"
CPUS="${SLURM_CPUS_PER_TASK:-16}"

line=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n' "${TASKS}")
line="${line%$'\r'}"
[[ -n "${line}" ]] || { echo "No task line $((SLURM_ARRAY_TASK_ID + 1)) in ${TASKS}" >&2; exit 1; }
IFS=$'\t' read -r event_id family aln model ml_tree constraint <<< "${line}"
prefix="$(dirname "${constraint}")/${event_id}"

if [[ -s "${prefix}.au.iqtree" ]] && grep -q "p-AU" "${prefix}.au.iqtree"; then
    echo "Already done: ${prefix}.au.iqtree"; exit 0
fi

echo "[$(date -Is)] ${event_id}: constrained search (${model})"
if [[ ! -s "${prefix}.constr.treefile" ]]; then
    "${IQTREE}" -s "${aln}" -m "${model}" -g "${constraint}" -T "${CPUS}" --seed 20261124 \
        --prefix "${prefix}.constr" -redo
fi

cat "${ml_tree}" "${prefix}.constr.treefile" > "${prefix}.candidates.nwk"
echo "[$(date -Is)] ${event_id}: AU test"
"${IQTREE}" -s "${aln}" -m "${model}" -z "${prefix}.candidates.nwk" -n 0 -zb 10000 -au \
    -T "${CPUS}" --seed 20261124 --prefix "${prefix}.au" -redo
echo "[$(date -Is)] ${event_id}: done"
