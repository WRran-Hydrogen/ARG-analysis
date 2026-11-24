#!/bin/bash
# AU-test driver (run from scripts/hpc/04_hgt_inference/ after 01_hgt_inference.py has written
# au_test_tasks.tsv): one 03_run_au_test.sh array task per event, then 04_summarize_au_tests.py.
#   bash 05_submit_au_tests.sh [extra sbatch options, e.g. --dependency=afterok:<jobid>]

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
P4="${SBL_HPC_HGT_DIR:-${SBL_HPC_BASE}/04_hgt_inference}"
mkdir -p "${P4}/logs"
cd "${SCRIPT_DIR}"
PY="${SBL_HPC_PYTHON:-python}"

# The dispatcher reads the task count at run time, so it can wait for 01_hgt_inference.py.
sbatch --parsable "$@" --partition=compute --cpus-per-task=1 --mem=1G --time=00:10:00 \
    --job-name=sbl_p4_au_disp --output="${P4}/logs/p4_au_disp_%j.out" \
    --wrap="set -euo pipefail; cd ${SCRIPT_DIR}
n=\$((\$(wc -l < ${P4}/au_test_tasks.tsv) - 1))
if [[ \${n} -lt 1 ]]; then echo 'no AU tasks'; exit 0; fi
au=\$(sbatch --parsable --array=1-\${n} --output=${P4}/logs/p4_au_%A_%a.out 03_run_au_test.sh ${P4}/au_test_tasks.tsv)
sm=\$(sbatch --parsable --dependency=afterany:\${au} --partition=compute --cpus-per-task=1 --mem=4G --time=00:20:00 --job-name=sbl_p4_au_sum --output=${P4}/logs/p4_au_sum_%j.out --wrap=\"cd ${SCRIPT_DIR}; ${PY} -u 04_summarize_au_tests.py\")
echo \"au_array=\${au} (\${n} events) summary=\${sm}\""
