#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --job-name=sbl_p3_cl_c20
#SBATCH --output=p3_cl_c20_%A_%a.out
#
# Structure-guided CL0013 superfamily trees: LG+F+G guide tree -> LG+C20+F+G PMSF + UFBoot 1000 / SH-aLRT 1000
# for each FoldMason+hmmalign CL0013 superfamily alignment (no_PF27310; task file from
# 07_build_structural_cl0013_and_clean_outgroups.py). usage: sbatch --array=1-N 08_run_cl0013_c20.sh <tasks.tsv>

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
TASKS="$1"
P3="${SBL_HPC_PHASE3_DIR:-${SBL_HPC_BASE}/03_alignment_and_trees}"
IQTREE="${SBL_HPC_IQTREE:-iqtree2}"
CPUS="${SLURM_CPUS_PER_TASK:-32}"

line=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n' "${TASKS}")
line="${line%$'\r'}"
[[ -n "${line}" ]] || { echo "No task line $((SLURM_ARRAY_TASK_ID + 1)) in ${TASKS}" >&2; exit 1; }
IFS=$'\t' read -r wave name aln model_spec <<< "${line}"

outdir="${P3}/trees/CL0013"
mkdir -p "${outdir}"
prefix="${outdir}/${name}"

if [[ -s "${prefix}.treefile" && -s "${prefix}.iqtree" ]]; then
    echo "Already done: ${prefix}.treefile"; exit 0
fi

echo "[$(date -Is)] ${name}: step 1/2 guide tree LG+F+G (${CPUS} CPUs)"
"${IQTREE}" -s "${aln}" -m "LG+F+G" -T "${CPUS}" --seed 20261124 --prefix "${prefix}.guide" -redo

echo "[$(date -Is)] ${name}: step 2/2 PMSF LG+C20+F+G + UFBoot 1000 / SH-aLRT 1000"
"${IQTREE}" -s "${aln}" -m "LG+C20+F+G" -ft "${prefix}.guide.treefile" \
    -B 1000 --bnni -alrt 1000 -T "${CPUS}" --seed 20261124 --prefix "${prefix}" -redo
echo "[$(date -Is)] ${name}: done"
