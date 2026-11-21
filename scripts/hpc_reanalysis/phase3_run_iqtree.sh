#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=sbl_p3_tree
#SBATCH --output=p3_tree_%A_%a.out
#
# Phase 3, step 2: one IQ-TREE run per line of a task file (Slurm array; submitted by phase3_submit.sh).
#   usage: sbatch --array=1-N phase3_run_iqtree.sh <tree_tasks_waveK.tsv>
# Wave 1 runs ModelFinder (-m MFP -mset LG,WAG,Q.pfam,VT,JTT). Wave 2 reuses the best-fit model
# of the wave-1 tree named in the task's model_from column (pre-registered compute budget).
# Support: UFBoot 1000 (+ --bnni) and SH-aLRT 1000. Trees stay unrooted; rooting is done in Phase 4.

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
for env_file in "${SCRIPT_DIR}/../../config/hpc.env" "${SCRIPT_DIR}/hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
TASKS="$1"
P3="${SBL_HPC_PHASE3_DIR:-${SBL_HPC_BASE}/03_alignment_and_trees}"
IQTREE="${SBL_HPC_IQTREE:-iqtree2}"

line=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n' "${TASKS}")
IFS=$'\t' read -r wave name aln model_from <<< "${line}"
outdir="${P3}/trees/${name%%.*}"
mkdir -p "${outdir}"
prefix="${outdir}/${name}"

if [[ -s "${prefix}.treefile" && -s "${prefix}.iqtree" ]]; then
    echo "Already done: ${prefix}.treefile"; exit 0
fi

if [[ -z "${model_from}" ]]; then
    model_args=(-m MFP -mset LG,WAG,Q.pfam,VT,JTT)
else
    ref="${P3}/trees/${model_from%%.*}/${model_from}.iqtree"
    model=$(grep -m1 "Best-fit model according to BIC:" "${ref}" | sed 's/.*BIC: *//')
    [[ -n "${model}" ]] || { echo "No best-fit model in ${ref}" >&2; exit 1; }
    model_args=(-m "${model}")
fi

echo "[$(date -Is)] ${name}: ${aln} model=${model_args[*]}"
"${IQTREE}" -s "${aln}" "${model_args[@]}" -B 1000 --bnni -alrt 1000 \
    -T "${SLURM_CPUS_PER_TASK:-16}" --seed 20261121 --prefix "${prefix}" -redo
echo "[$(date -Is)] ${name}: done"
