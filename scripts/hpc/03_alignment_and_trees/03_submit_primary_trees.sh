#!/bin/bash
# Phase 3 driver: alignments -> wave-1 trees (ModelFinder) -> wave-2 trees (fixed models).
# Run on the login node from scripts/hpc/03_alignment_and_trees/:
#   bash 03_submit_primary_trees.sh
# Submits three dependent Slurm jobs; nothing heavy runs on the login node.
# Only run after ../run_curation_and_context.sh has finished.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
P3="${SBL_HPC_PHASE3_DIR:-${SBL_HPC_BASE}/03_alignment_and_trees}"
mkdir -p "${P3}/logs"
cd "${SCRIPT_DIR}"

# Step 1: alignments + task lists (single job)
aln_job=$(sbatch --parsable --partition=compute --cpus-per-task=16 --mem=32G --time=12:00:00 \
    --job-name=sbl_p3_aln --output="${P3}/logs/p3_aln_%j.out" \
    --wrap="set -euo pipefail; cd ${SCRIPT_DIR}; ${SBL_HPC_PYTHON:-python} -u 01_build_alignments.py")
echo "alignments: job ${aln_job}"

# Steps 2-3 are submitted by a small follow-up job, because the array sizes are only known once
# the task lists exist.
sbatch --parsable --dependency=afterok:"${aln_job}" --partition=compute --cpus-per-task=1 --mem=1G \
    --time=00:10:00 --job-name=sbl_p3_dispatch --output="${P3}/logs/p3_dispatch_%j.out" \
    --wrap="set -euo pipefail; cd ${SCRIPT_DIR}
n1=\$((\$(wc -l < ${P3}/tree_tasks_wave1.tsv) - 1)); n2=\$((\$(wc -l < ${P3}/tree_tasks_wave2.tsv) - 1))
w1=\$(sbatch --parsable --array=1-\${n1} --output=${P3}/logs/p3_tree_%A_%a.out 02_run_iqtree.sh ${P3}/tree_tasks_wave1.tsv)
w2=\$(sbatch --parsable --dependency=afterok:\${w1} --array=1-\${n2}%12 --output=${P3}/logs/p3_tree_%A_%a.out 02_run_iqtree.sh ${P3}/tree_tasks_wave2.tsv)
echo \"wave1=\${w1} (\${n1} trees) wave2=\${w2} (\${n2} trees)\""
echo "dispatcher queued; see ${P3}/logs/p3_dispatch_*.out for the tree job ids"
