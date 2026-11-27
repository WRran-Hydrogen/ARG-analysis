#!/bin/bash
# Selection analysis driver (run from scripts/hpc/05_selection/):
#   1) Build 1:1 codon alignments, dS saturation/Xia's test, motif purifying selection tests, & HyPhy task list
#   2) Dispatch HyPhy Slurm array (aBSREL + MEME + FEL + GARD for subclades)
#   3) Summarize HyPhy outputs once array completes

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
P5="${SBL_HPC_BASE}/05_selection_pressure"
mkdir -p "${P5}/logs"
cd "${SCRIPT_DIR}"

PY_ETE="${SBL_HPC_PYTHON:-python}"

prep_job=$(sbatch --parsable --partition=compute --cpus-per-task=16 --mem=32G --time=02:00:00 \
    --job-name=sbl_p5_prep --output="${P5}/logs/p5_prep_%j.out" \
    --wrap="set -euo pipefail; cd ${SCRIPT_DIR}; ${PY_ETE} -u 01_build_codon_and_saturation.py")

disp_job=$(sbatch --parsable --dependency=afterok:"${prep_job}" --partition=compute --cpus-per-task=1 --mem=2G --time=00:10:00 \
    --job-name=sbl_p5_disp --output="${P5}/logs/p5_disp_%j.out" \
    --wrap="set -euo pipefail; cd ${SCRIPT_DIR}
n_tasks=\$((\$(wc -l < ${P5}/hyphy_tasks.tsv) - 1))
hy_job=\$(sbatch --parsable --array=1-\${n_tasks} --output=${P5}/logs/p5_hyphy_%A_%a.out 02_run_hyphy.sh ${P5}/hyphy_tasks.tsv)
sum_job=\$(sbatch --parsable --dependency=afterany:\${hy_job} --partition=compute --cpus-per-task=2 --mem=8G --time=00:30:00 --job-name=sbl_p5_sum --output=${P5}/logs/p5_sum_%j.out --wrap=\"set -euo pipefail; cd ${SCRIPT_DIR}; ${PY_ETE} -u 03_summarize_hyphy.py\")
echo \"hyphy_array=\${hy_job} (\${n_tasks} sets) summary_job=\${sum_job}\"")

echo "prep_job=${prep_job}"
echo "dispatch_job=${disp_job}"
