#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=sbl_p3_pmsf
#SBATCH --output=p3_pmsf_%A_%a.out
#
# Site-heterogeneous PMSF (LG+C20+F+G) long-branch attraction check.
# Uses the finished Wave-1 primary ML tree as the guide tree (-ft).

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
TASKS="$1"
P3="${SBL_HPC_PHASE3_DIR:-${SBL_HPC_BASE}/03_alignment_and_trees}"
IQTREE="${SBL_HPC_IQTREE:-iqtree2}"

line=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n' "${TASKS}")
line="${line%$'\r'}"
[[ -n "${line}" ]] || { echo "No task line $((SLURM_ARRAY_TASK_ID + 1)) in ${TASKS}" >&2; exit 1; }
IFS=$'\t' read -r name aln guide_tree <<< "${line}"
outdir="${P3}/trees/${name%%.*}"
mkdir -p "${outdir}"
prefix="${outdir}/${name}"

if [[ -s "${prefix}.treefile" && -s "${prefix}.iqtree" ]]; then
    echo "Already done: ${prefix}.treefile"; exit 0
fi

echo "[$(date -Is)] ${name}: ${aln} guide=${guide_tree} model=LG+C20+F+G"
"${IQTREE}" -s "${aln}" -m "LG+C20+F+G" -ft "${guide_tree}" -B 1000 --bnni -alrt 1000 \
    -T "${SLURM_CPUS_PER_TASK:-16}" --seed 20261121 --prefix "${prefix}" -redo
echo "[$(date -Is)] ${name}: done"
