#!/bin/bash
# Reference shell script: HMM search and MEME/FIMO motif filtering (documentation of the commands
# that scripts 01 & 02 run automatically).
# Run inside the staging folder created by scripts 01/02 (.wsl_stage/wsl_m01), where the
# *_seed.faa, *_Viral.faa, *_Bacterial.faa and *.clean.faa inputs are staged.
# Usage: cd .wsl_stage/wsl_m01 && bash ../../scripts/run_hmm_motif_pipeline.sh
# Tool locations: export SBL_MUSCLE / SBL_MEME_BIN_DIR, or set them in config/local.env.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/../config/local.env" ]]; then
    set -a; source "${SCRIPT_DIR}/../config/local.env"; set +a
fi
MUSCLE="${SBL_MUSCLE:-muscle}"
MEME="${SBL_MEME_BIN_DIR:+${SBL_MEME_BIN_DIR}/}meme"
FIMO="${SBL_MEME_BIN_DIR:+${SBL_MEME_BIN_DIR}/}fimo"

FAMILIES=("PF00144" "PF00905" "PF13354")

for FAM in "${FAMILIES[@]}"; do
    echo "=== Processing ${FAM} ==="
    # 1. Align Pfam seed sequences with MUSCLE v5.3
    "${MUSCLE}" -align "${FAM}_seed.faa" -output "${FAM}_seed.aln.faa"

    # 2. Build profile HMM with HMMER v3.4
    hmmbuild -n "${FAM}" "${FAM}.hmm" "${FAM}_seed.aln.faa"

    # 3. Screen viral and bacterial protein sets (rebuilt seed HMM, E-value screen)
    #    Official Pfam gathering thresholds are applied separately with --cut_ga (01_build_and_verify_hmm.py).
    hmmsearch --cpu 4 -E 1e-10 --domE 1e-10 --tblout "${FAM}_Viral.tblout" --domtblout "${FAM}_Viral.domtblout" "${FAM}.hmm" "${FAM}_Viral.faa"
    hmmsearch --cpu 4 -E 1e-10 --domE 1e-10 --tblout "${FAM}_Bacterial.tblout" --domtblout "${FAM}_Bacterial.domtblout" "${FAM}.hmm" "${FAM}_Bacterial.faa"

    # 4. MEME/FIMO motif filtering under identical thresholds for viral and bacterial sets (1e-5 and 1e-10)
    "${MEME}" "${FAM}_seed.clean.faa" -protein -mod zoops -nmotifs 5 -minw 6 -maxw 20 -maxsize 1000000 -oc "meme_r1_${FAM}"
    "${FIMO}" --thresh 1e-5 --oc "fimo_r1_viral_${FAM}_1e5" "meme_r1_${FAM}/meme.txt" "${FAM}_Viral_Initial.clean.faa"
    "${FIMO}" --thresh 1e-5 --oc "fimo_r1_bac_${FAM}_1e5" "meme_r1_${FAM}/meme.txt" "${FAM}_Bacterial_Ref.clean.faa"
    "${FIMO}" --thresh 1e-10 --oc "fimo_r1_viral_${FAM}_1e10" "meme_r1_${FAM}/meme.txt" "${FAM}_Viral_Initial.clean.faa"
    "${FIMO}" --thresh 1e-10 --oc "fimo_r1_bac_${FAM}_1e10" "meme_r1_${FAM}/meme.txt" "${FAM}_Bacterial_Ref.clean.faa"
done
