"""Bacterial reference dataset prophage & viral decontamination audit.

Systematically audits all 2,220 sequences in the bacterial SBL
reference dataset (SARG/CARD + NCBI nr) to identify:
1. Sequences explicitly classified under realm 'Viruses' (e.g., Punavirus/Hendrixvirinae) in SARG metadata,
2. Artificial/synthetic/vector sequences ('other entries'),
3. FASTA header phage/prophage/MGE annotations or Prophage_identify.blastnout matches,
4. NCBI nr proteins sharing >= 90% amino acid identity (>= 80% coverage) with IMG/VR viral SBLs
   (potential unannotated cryptic prophages integrated in bacterial reference genomes).

Generates sensitivity-cleaned bacterial FASTAs for downstream phylogenetic and HGT testing (05_phylogenetic_hgt_inference.py).
"""

from pathlib import Path
import os
import re
import subprocess
import pandas as pd

import sbl_config as cfg
from sbl_config import linux_cmd

RAW_DIR = str(cfg.SEQ_DIR)
CLEAN_OUT_DIR = str(cfg.BAC_CLEAN_DIR)  # CleanStrict / CleanCore FASTAs (output)
META_DIR = cfg.META_DIR
TABLES_DIR = str(cfg.TABLES_DIR)
SUPP_DIR = str(cfg.SUPP_DIR)
_stage, WSL_STAGE_LIN = cfg.stage_dirs("wsl_m01")
WSL_STAGE_WIN = str(_stage)

INFO_XLSX = str(META_DIR / "information.xlsx")
PROPHAGE_BLAST = str(META_DIR / "Prophage_identify.blastnout")
INPHAGE_XLSX = str(META_DIR / "Inphage_gggenome.xlsx")

FAMILIES = ["PF00144", "PF00905", "PF13354"]

PHAGE_KW_RE = re.compile(
    r"\b(bacteriophage|prophage|phage|punavirus|virus|viral|virion|capsid|terminase|integrase)\b",
    re.IGNORECASE,
)
SYNTH_KW_RE = re.compile(
    r"(synthetic construct|artificial sequence|cloning vector|expression vector|pdb\|)",
    re.IGNORECASE,
)


def read_fasta(path):
    records = []
    cur_id, cur_hdr, seq = None, None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, cur_hdr, "".join(seq).replace("*", "")))
                cur_hdr = line[1:].strip()
                cur_id = cur_hdr.split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if cur_id is not None:
            records.append((cur_id, cur_hdr, "".join(seq).replace("*", "")))
    return records


def run_wsl_cmd(cmd):
    proc = subprocess.run(linux_cmd("bash", "-c", cmd), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"WSL command failed:\n{cmd}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc.stdout


def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(SUPP_DIR, exist_ok=True)
    os.makedirs(WSL_STAGE_WIN, exist_ok=True)

    # Load SARG taxonomy table
    sarg_df = pd.read_excel(INFO_XLSX, sheet_name="SARG")
    sarg_df["SARG-id"] = sarg_df["SARG-id"].astype(str).str.strip()
    sarg_map = sarg_df.drop_duplicates("SARG-id").set_index("SARG-id").to_dict("index")

    # Load known prophage bacterial genome accessions
    prophage_bac_accs = set()
    if os.path.exists(PROPHAGE_BLAST):
        with open(PROPHAGE_BLAST, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    # e.g. ref|NZ_CP034181.1|
                    acc = parts[1].replace("ref|", "").strip("|").split(".")[0]
                    prophage_bac_accs.add(acc)
    if os.path.exists(INPHAGE_XLSX):
        inph_bac = pd.read_excel(INPHAGE_XLSX, sheet_name="Bac")
        for acc in inph_bac["seq_id"].dropna().astype(str):
            prophage_bac_accs.add(acc.split(".")[0])

    # Load viral contig checkv metadata so we can see if high-identity matches are to non-prophage viral contigs
    info_all = pd.read_excel(INFO_XLSX, sheet_name="ALL_information")
    info_all["UVIG"] = info_all["UVIG"].astype(str)
    uvig_meta = info_all.set_index("UVIG").to_dict("index")

    # Stage combined Viral_Final and Bacterial_Ref FASTAs for BLASTP
    all_viral_fa = os.path.join(WSL_STAGE_WIN, "all_viral_final.clean.faa")
    all_bac_fa = os.path.join(WSL_STAGE_WIN, "all_bacterial_ref.clean.faa")

    bac_records_by_fam = {}
    with open(all_viral_fa, "w", encoding="utf-8") as fv:
        for fam in FAMILIES:
            for sid, hdr, seq in read_fasta(os.path.join(RAW_DIR, f"{fam}_Viral_Final.faa")):
                fv.write(f">{sid}\n{seq}\n")

    with open(all_bac_fa, "w", encoding="utf-8") as fb:
        for fam in FAMILIES:
            recs = read_fasta(os.path.join(RAW_DIR, f"{fam}_Bacterial.faa"))
            bac_records_by_fam[fam] = recs
            for sid, hdr, seq in recs:
                fb.write(f">{fam}::{sid}\n{seq}\n")

    # Run BLASTP in WSL: bacterial references vs viral final 5,701 SBLs
    cmd_blast = (
        f"cd {WSL_STAGE_LIN} && "
        f"makeblastdb -in all_viral_final.clean.faa -dbtype prot -out viral_sbl_db > /dev/null && "
        f"blastp -query all_bacterial_ref.clean.faa -db viral_sbl_db -num_threads 4 -evalue 1e-10 "
        f"-max_target_seqs 5 -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' "
        f"-out bac_vs_viral.blastp.tsv"
    )
    print("Running BLASTP (Bacterial Ref vs Viral Final 5,701)...")
    run_wsl_cmd(cmd_blast)

    # Parse top viral hit per bacterial sequence
    top_hit = {}
    blast_out = os.path.join(WSL_STAGE_WIN, "bac_vs_viral.blastp.tsv")
    with open(blast_out, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            q_full, sseqid = parts[0], parts[1]
            pident = float(parts[2])
            alen = int(parts[3])
            qlen = int(parts[4])
            slen = int(parts[5])
            qcovs = float(parts[6])
            eval_ = float(parts[7])
            bscore = float(parts[8])
            if q_full not in top_hit or (pident, bscore) > (top_hit[q_full]["pident"], top_hit[q_full]["bitscore"]):
                uvig = sseqid.split("|")[0]
                umeta = uvig_meta.get(uvig, {})
                top_hit[q_full] = {
                    "top_viral_hit": sseqid,
                    "top_viral_uvig": uvig,
                    "pident": pident,
                    "qcovs": qcovs,
                    "aln_cov": round(100.0 * alen / max(qlen, 1), 2),
                    "evalue": eval_,
                    "bitscore": bscore,
                    "viral_contam": umeta.get("Estimated contamination", float("nan")),
                    "viral_miuvig": umeta.get("MIUViG quality", ""),
                    "viral_coords": str(umeta.get("Coordinates ('whole' if the UViG is the entire contig)", "")),
                }

    audit_rows = []
    summary_rows = []

    for fam in FAMILIES:
        recs = bac_records_by_fam[fam]
        clean_a = []
        clean_b = []

        for sid, hdr, seq in recs:
            smeta = sarg_map.get(sid, {})
            realm = str(smeta.get("realm", "Unknown"))
            kingdom = str(smeta.get("kindom", "Unknown"))
            phylum = str(smeta.get("phylum", ""))
            species = str(smeta.get("species", ""))

            is_ncbi_nr = sid.startswith("gi|")
            is_curated_card = ("ARO:" in hdr) or (".gene" in sid)
            flag_sarg_virus = (realm == "Viruses")
            flag_sarg_artificial = (realm == "other entries")
            flag_hdr_phage = bool(PHAGE_KW_RE.search(hdr)) or bool(PHAGE_KW_RE.search(species))
            flag_hdr_synth = bool(SYNTH_KW_RE.search(hdr))

            q_key = f"{fam}::{sid}"
            hit = top_hit.get(q_key, {})
            pident = hit.get("pident", 0.0)
            qcovs = hit.get("qcovs", 0.0)
            # High-identity viral match (>=90% aa identity and >=80% query coverage)
            flag_high_viral_id90 = (pident >= 90.0 and qcovs >= 80.0)
            flag_high_viral_id95 = (pident >= 95.0 and qcovs >= 80.0)

            # Tier A (CleanStrict): Exclude explicit SARG Viruses, Artificial/Synthetic, header phage/synth
            exclude_tier_a = flag_sarg_virus or flag_sarg_artificial or flag_hdr_phage or flag_hdr_synth
            # Tier B (CleanCore): Tier A + exclude any NCBI nr protein with >=90% aa identity to viral SBL
            exclude_tier_b = exclude_tier_a or (is_ncbi_nr and flag_high_viral_id90)

            if not exclude_tier_a:
                clean_a.append((sid, hdr, seq))
            if not exclude_tier_b:
                clean_b.append((sid, hdr, seq))

            audit_rows.append({
                "family": fam,
                "seq_id": sid,
                "source_type": "NCBI_nr" if is_ncbi_nr else ("CARD_SARG_curated" if is_curated_card else "SARG_other"),
                "sarg_realm": realm,
                "sarg_kingdom": kingdom,
                "sarg_phylum": phylum,
                "sarg_species": species,
                "flag_sarg_virus": flag_sarg_virus,
                "flag_sarg_artificial": flag_sarg_artificial,
                "flag_header_phage": flag_hdr_phage,
                "flag_header_synthetic": flag_hdr_synth,
                "top_viral_hit": hit.get("top_viral_hit", ""),
                "top_viral_pident": pident,
                "top_viral_qcovs": qcovs,
                "flag_nr_viral_id90": is_ncbi_nr and flag_high_viral_id90,
                "flag_nr_viral_id95": is_ncbi_nr and flag_high_viral_id95,
                "pass_CleanStrict_TierA": not exclude_tier_a,
                "pass_CleanCore_TierB": not exclude_tier_b,
            })

        # Write cleaned bacterial reference FASTAs
        with open(os.path.join(CLEAN_OUT_DIR, f"{fam}_Bacterial_CleanStrict.faa"), "w", encoding="utf-8") as fa:
            for sid, hdr, seq in clean_a:
                fa.write(f">{hdr}\n{seq}\n")
        with open(os.path.join(CLEAN_OUT_DIR, f"{fam}_Bacterial_CleanCore.faa"), "w", encoding="utf-8") as fb:
            for sid, hdr, seq in clean_b:
                fb.write(f">{hdr}\n{seq}\n")

        fam_df = pd.DataFrame([r for r in audit_rows if r["family"] == fam])
        summary_rows.append({
            "family": fam,
            "total_bacterial_ref": len(fam_df),
            "ncbi_nr_count": int((fam_df["source_type"] == "NCBI_nr").sum()),
            "flag_sarg_virus": int(fam_df["flag_sarg_virus"].sum()),
            "flag_sarg_artificial": int(fam_df["flag_sarg_artificial"].sum()),
            "flag_header_phage_or_synth": int((fam_df["flag_header_phage"] | fam_df["flag_header_synthetic"]).sum()),
            "flag_nr_viral_id90": int(fam_df["flag_nr_viral_id90"].sum()),
            "flag_nr_viral_id95": int(fam_df["flag_nr_viral_id95"].sum()),
            "retained_CleanStrict_TierA": int(fam_df["pass_CleanStrict_TierA"].sum()),
            "retained_CleanCore_TierB": int(fam_df["pass_CleanCore_TierB"].sum()),
        })

    df_audit = pd.DataFrame(audit_rows)
    df_sum = pd.DataFrame(summary_rows)
    # Add pooled total row
    pooled = {
        "family": "All_Pooled",
        "total_bacterial_ref": int(df_sum["total_bacterial_ref"].sum()),
        "ncbi_nr_count": int(df_sum["ncbi_nr_count"].sum()),
        "flag_sarg_virus": int(df_sum["flag_sarg_virus"].sum()),
        "flag_sarg_artificial": int(df_sum["flag_sarg_artificial"].sum()),
        "flag_header_phage_or_synth": int(df_sum["flag_header_phage_or_synth"].sum()),
        "flag_nr_viral_id90": int(df_sum["flag_nr_viral_id90"].sum()),
        "flag_nr_viral_id95": int(df_sum["flag_nr_viral_id95"].sum()),
        "retained_CleanStrict_TierA": int(df_sum["retained_CleanStrict_TierA"].sum()),
        "retained_CleanCore_TierB": int(df_sum["retained_CleanCore_TierB"].sum()),
    }
    df_sum = pd.concat([df_sum, pd.DataFrame([pooled])], ignore_index=True)

    df_audit.to_csv(os.path.join(SUPP_DIR, "PerSequence_Bacterial_Audit.tsv"), sep="\t", index=False)
    df_sum.to_csv(os.path.join(TABLES_DIR, "Table_S03_Bacterial_Reference_Prophage_Audit.csv"), index=False)
    print("\n=== Bacterial Reference Prophage & Viral Decontamination Summary ===")
    print(df_sum.to_string(index=False))


if __name__ == "__main__":
    main()
