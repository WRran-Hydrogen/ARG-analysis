#!/bin/bash
# CL0013 superfamily and PMSF driver (run from scripts/hpc/03_alignment_and_trees/):
#   1) Build CL0013 superfamily + LOO alignments and PMSF C20 task list (sbl_p3_cl0013_aln)
#   2) Run PMSF (LG+C20+F+G) array (5 trees)
#   3) Run CL0013 superfamily + LOO IQ-TREE array (4 trees)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
P3="${SBL_HPC_PHASE3_DIR:-${SBL_HPC_BASE}/03_alignment_and_trees}"
mkdir -p "${P3}/logs"
cd "${SCRIPT_DIR}"

aln_job=$(sbatch --parsable --partition=compute --cpus-per-task=16 --mem=32G --time=04:00:00 \
    --job-name=sbl_p3_cl0013_aln --output="${P3}/logs/p3_cl0013_aln_%j.out" \
    --wrap="set -euo pipefail; cd ${SCRIPT_DIR}; ${SBL_HPC_PYTHON:-python} -u 04_build_cl0013_and_pmsf_tasks.py")

pmsf_job=$(sbatch --parsable --dependency=afterok:"${aln_job}" --array=1-5 \
    --output="${P3}/logs/p3_pmsf_%A_%a.out" \
    05_run_pmsf.sh "${P3}/tree_tasks_pmsf.tsv")

cl_job=$(sbatch --parsable --dependency=afterok:"${aln_job}" --array=1-4 \
    --job-name=sbl_p3_cl0013 --output="${P3}/logs/p3_cl0013_%A_%a.out" \
    02_run_iqtree.sh "${P3}/tree_tasks_cl0013.tsv")

echo "cl0013_aln=${aln_job}"
echo "pmsf=${pmsf_job}"
echo "cl0013=${cl_job}"
