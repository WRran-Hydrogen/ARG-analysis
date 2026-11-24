#!/bin/bash
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=sbl_p5_hyphy
#SBATCH --output=p5_hyphy_%A_%a.out
#
# Selection analysis, step 2: run HyPhy aBSREL, MEME, FEL (and GARD for subclades)
# for one dataset line of hyphy_tasks.tsv.

set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
for env_file in "${SCRIPT_DIR}/../../../config/hpc.env" "${SCRIPT_DIR}/../hpc.env"; do
    if [[ -f "${env_file}" ]]; then set -a; source "${env_file}"; set +a; fi
done
: "${SBL_HPC_BASE:?Set SBL_HPC_BASE}"
TASKS="$1"
HYPHY="${SBL_HPC_HYPHY:-hyphy}"
CPUS="${SLURM_CPUS_PER_TASK:-16}"

line=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n' "${TASKS}")
line="${line%$'\r'}"
[[ -n "${line}" ]] || { echo "No task line $((SLURM_ARRAY_TASK_ID + 1)) in ${TASKS}" >&2; exit 1; }
IFS=$'\t' read -r set_id family set_group track n_seqs n_codons med_ds iss iss_cs iss_ca set_dir <<< "${line}"

cd "${set_dir}"
echo "[$(date -Is)] ${set_id} (${family}, ${set_group}, n=${n_seqs}, cod=${n_codons})"

if [[ ! -s "absrel.json" ]]; then
    "${HYPHY}" CPU="${CPUS}" absrel --alignment codon.fna --tree tree.nwk --output absrel.json
fi
if [[ ! -s "meme.json" ]]; then
    "${HYPHY}" CPU="${CPUS}" meme --alignment codon.fna --tree tree.nwk --output meme.json
fi
if [[ ! -s "fel.json" ]]; then
    "${HYPHY}" CPU="${CPUS}" fel --alignment codon.fna --tree tree.nwk --output fel.json
fi
if [[ "${set_group}" == "subclade" && ! -s "gard.json" ]]; then
    timeout 25m "${HYPHY}" CPU="${CPUS}" gard --type nucleotide --model GTR --rv GDD --rate-classes 3 --mode Faster \
        --alignment codon.fna --output gard.json --output-lf gard.bf || true
fi

echo "[$(date -Is)] ${set_id}: done"
