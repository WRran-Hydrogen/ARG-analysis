#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=sbl_p1p2_rebuild
#SBATCH --output=p1p2_rebuild_%j.out
#
# Curated datasets + genomic context. The three steps must run in this order:
#   1. 01_curation/09_build_curated_datasets.py --stage curate    (CL0013 assignment, own-seed FIMO, CARD/motifs)
#   2. 02_genomic_context/01_build_neighborhoods.py               (strict structural-hallmark tiers)
#   3. 01_curation/09_build_curated_datasets.py --stage finalize  (nr85, tree replicates, outgroups; reads step-2 tiers)
#
# Submit from scripts/hpc/:   sbatch run_curation_and_context.sh
# Paths come from hpc.env (next to this script or in config/) or exported SBL_HPC_* variables.

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

for env_file in "${SCRIPT_DIR}/../../config/hpc.env" "${SCRIPT_DIR}/hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE (see config/hpc.env.example)}"
PY="${SBL_HPC_PYTHON:-python}"

# Inputs needed in addition to the 01_viral_decontamination.py inputs
missing=0
for f in PF00768_seed.faa PF02113_seed.faa card_db/card_true_sbl_refs.faa \
         PF00144_official.hmm PF00905_official.hmm PF13354_official.hmm PF00768_official.hmm PF02113_official.hmm; do
    if [[ ! -s "${SBL_HPC_BASE}/inputs/${f}" ]]; then
        echo "MISSING input: ${SBL_HPC_BASE}/inputs/${f}" >&2
        missing=1
    fi
done
if [[ ${missing} -ne 0 ]]; then
    echo "Upload the missing files (seeds: data/sequences/, HMMs: data/hmm_models/) and resubmit." >&2
    exit 1
fi

cd "${SCRIPT_DIR}"
echo ">>> [1/3] Phase 1 curate  $(date)"
"${PY}" -u 01_curation/09_build_curated_datasets.py --stage curate
echo ">>> [2/3] Phase 2 strict tiers  $(date)"
"${PY}" -u 02_genomic_context/01_build_neighborhoods.py
echo ">>> [3/3] Phase 1 finalize  $(date)"
"${PY}" -u 01_curation/09_build_curated_datasets.py --stage finalize
echo ">>> Done  $(date)"
