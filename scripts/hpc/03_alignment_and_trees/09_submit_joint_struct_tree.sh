#!/bin/bash
# Structure-guided joint PF13354_PF00768 tree (run from scripts/hpc/03_alignment_and_trees/):
#   bash 09_submit_joint_struct_tree.sh [extra sbatch options for the first job]
# 1. 07_build_structural_cl0013_and_clean_outgroups.py --step partC (FoldMason-grid alignment;
#    needs the Part A FoldMason seed MSA)
# 2. 02_run_iqtree.sh wave 1: ModelFinder tree on primary.struct_out
# 3. 02_run_iqtree.sh wave 2: ingroup-only primary.struct_in tree with the wave-1 model

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
P3="${SBL_HPC_PHASE3_DIR:-${SBL_HPC_BASE}/03_alignment_and_trees}"
mkdir -p "${P3}/logs"
cd "${SCRIPT_DIR}"

aln=$(sbatch --parsable "$@" --partition=compute --cpus-per-task=4 --mem=16G --time=04:00:00 \
    --job-name=sbl_p3_joint_aln --output="${P3}/logs/p3_joint_struct_aln_%j.out" \
    --wrap="set -euo pipefail; cd ${SCRIPT_DIR}; ${SBL_HPC_PYTHON:-python} -u 07_build_structural_cl0013_and_clean_outgroups.py --step partC")
w1=$(sbatch --parsable --dependency=afterok:"${aln}" --array=1-1 \
    --output="${P3}/logs/p3_tree_%A_%a.out" 02_run_iqtree.sh "${P3}/tree_tasks_joint_struct_wave1.tsv")
w2=$(sbatch --parsable --dependency=afterok:"${w1}" --array=1-1 \
    --output="${P3}/logs/p3_tree_%A_%a.out" 02_run_iqtree.sh "${P3}/tree_tasks_joint_struct_wave2.tsv")
echo "alignment=${aln} wave1=${w1} wave2=${w2}"
