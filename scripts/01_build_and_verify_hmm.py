"""Module 01 - Step 1: Official Pfam CL0013 & Seed HMM Verification Across Viral & Bacterial SBLs.

Evaluates all initial viral candidates (n=11,368), final curated viral SBLs (n=5,701),
and bacterial reference sequences (n=2,220) across the three focal SBL families
(PF00144, PF00905, PF13354) using:
1. Official Pfam-A HMM models with curated Gathering Thresholds (--cut_ga) for PF00144.31
   (GA=27.0), PF00905.29 (GA=22.1), and PF13354.14 (GA=22.0).
2. Full Clan CL0013 competitive and auxiliary evolutionary profiling across all 8 CL0013
   members (penicillin-recognizing enzymes: PF00144, PF00905, PF13354, PF00768 [Peptidase_S11],
   PF02113 [Peptidase_S13]; non-penicillin-recognizing: PF04960 [Glutaminase],
   PF06089 [Asparaginase_II], PF27310 [HrcA_C_epsilon]).
3. Rebuilt Pfam seed HMM profiles (MUSCLE v5.3 + HMMER v3.4) for domain coverage and E-value sensitivity.
"""

from pathlib import Path
import os
import shutil
import subprocess
import pandas as pd

import sbl_config as cfg
from sbl_config import linux_cmd

RAW_DIR = str(cfg.SEQ_DIR)
HMM_DIR = str(cfg.HMM_DIR)              # official Pfam / CL0013 HMMs (input)
HMM_OUT_DIR = str(cfg.HMM_REBUILT_DIR)  # rebuilt seed HMMs + MUSCLE MSAs (output)
TABLES_DIR = str(cfg.TABLES_DIR)
SUPP_DIR = str(cfg.SUPP_DIR)
_stage, WSL_STAGE_LIN = cfg.stage_dirs("wsl_m01")
WSL_STAGE_WIN = str(_stage)

FAMILIES = ["PF00144", "PF00905", "PF13354"]
PRE_CL0013 = {"PF00144", "PF00905", "PF13354", "PF00768", "PF02113"}
NON_PRE_CL0013 = {"PF04960", "PF06089", "PF27310"}


def ensure_dirs():
    for d in [RAW_DIR, HMM_DIR, HMM_OUT_DIR, TABLES_DIR, SUPP_DIR, WSL_STAGE_WIN]:
        os.makedirs(d, exist_ok=True)


def stage_inputs():
    """Stage source FASTAs and official CL0013 HMMs into WSL staging folder."""
    for fam in FAMILIES:
        files_to_copy = [
            f"{fam}_seed.faa",
            f"{fam}_Viral.faa",
            f"{fam}_Viral_FIMO_r1.faa",
            f"{fam}_Viral_Final.faa",
            f"{fam}_Bacterial.faa",
        ]
        for fname in files_to_copy:
            src = os.path.join(RAW_DIR, fname)
            if not os.path.exists(src):
                raise FileNotFoundError(f"Missing required input FASTA: {src}")
            shutil.copy2(src, os.path.join(WSL_STAGE_WIN, fname))

    cl_hmm = os.path.join(HMM_DIR, "CL0013_official.hmm")
    if os.path.exists(cl_hmm):
        shutil.copy2(cl_hmm, os.path.join(WSL_STAGE_WIN, "CL0013_official.hmm"))


def run_wsl_cmd(cmd):
    full_cmd = linux_cmd("bash", "-c", cmd)
    proc = subprocess.run(full_cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"WSL command failed:\n{cmd}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc.stdout


def parse_domtblout(path):
    """Parse HMMER3 --domtblout taking the best-scoring domain per sequence ID."""
    best = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            seq_id = parts[0]
            tlen = int(parts[2])
            qlen = int(parts[5])
            seq_evalue = float(parts[6])
            seq_score = float(parts[7])
            dom_ievalue = float(parts[12])
            dom_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            acc = float(parts[21])
            hmm_cov = (hmm_to - hmm_from + 1) / max(qlen, 1)

            rec = {
                "seq_id": seq_id,
                "seq_len": tlen,
                "hmm_len": qlen,
                "seq_evalue": seq_evalue,
                "seq_bitscore": seq_score,
                "dom_ievalue": dom_ievalue,
                "dom_bitscore": dom_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "hmm_coverage": round(hmm_cov, 4),
                "ali_from": ali_from,
                "ali_to": ali_to,
                "mean_acc": acc,
            }
            if seq_id not in best or seq_score > best[seq_id]["seq_bitscore"]:
                best[seq_id] = rec
    return best


def parse_cl0013_multi_hits(path):
    """Parse all official CL0013 hits passing --cut_ga per sequence ID."""
    hits = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            seq_id = parts[0]
            qname = parts[3]
            qacc = parts[4].split(".")[0]
            qlen = int(parts[5])
            seq_ev = float(parts[6])
            seq_bit = float(parts[7])
            dom_bit = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            cov = round((hmm_to - hmm_from + 1) / max(qlen, 1), 4)
            fam_map = hits.setdefault(seq_id, {})
            if qacc not in fam_map or dom_bit > fam_map[qacc]["dom_bit"]:
                fam_map[qacc] = {
                    "qacc": qacc,
                    "qname": qname,
                    "seq_ev": seq_ev,
                    "seq_bit": seq_bit,
                    "dom_bit": dom_bit,
                    "cov": cov,
                }
    return hits


def get_fasta_ids(path):
    ids = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def main():
    ensure_dirs()
    stage_inputs()

    summary_rows = []
    all_seq_rows = []

    for fam in FAMILIES:
        print(f"Processing {fam} (official CL0013 --cut_ga + seed HMM)...")
        bash_script = (
            f"cd {WSL_STAGE_LIN} && "
            f"{cfg.MUSCLE} -align {fam}_seed.faa -output {fam}_seed.aln.faa && "
            f"hmmbuild -n {fam} {fam}.hmm {fam}_seed.aln.faa && "
            f"hmmsearch --cpu 4 -E 10.0 --domE 10.0 --tblout {fam}_Viral.tblout --domtblout {fam}_Viral.domtblout {fam}.hmm {fam}_Viral.faa > /dev/null && "
            f"hmmsearch --cpu 4 -E 10.0 --domE 10.0 --tblout {fam}_Bacterial.tblout --domtblout {fam}_Bacterial.domtblout {fam}.hmm {fam}_Bacterial.faa > /dev/null && "
            f"hmmsearch --cpu 4 --cut_ga --domtblout {fam}_Viral_CL0013_ga.domtblout CL0013_official.hmm {fam}_Viral.faa > /dev/null && "
            f"hmmsearch --cpu 4 --cut_ga --domtblout {fam}_Bacterial_CL0013_ga.domtblout CL0013_official.hmm {fam}_Bacterial.faa > /dev/null"
        )
        run_wsl_cmd(bash_script)

        shutil.copy2(os.path.join(WSL_STAGE_WIN, f"{fam}.hmm"), os.path.join(HMM_OUT_DIR, f"{fam}.hmm"))
        shutil.copy2(os.path.join(WSL_STAGE_WIN, f"{fam}_seed.aln.faa"), os.path.join(HMM_OUT_DIR, f"{fam}_seed.aln.faa"))

        viral_all_ids = get_fasta_ids(os.path.join(RAW_DIR, f"{fam}_Viral.faa"))
        viral_r1_ids = set(get_fasta_ids(os.path.join(RAW_DIR, f"{fam}_Viral_FIMO_r1.faa")))
        viral_final_ids = set(get_fasta_ids(os.path.join(RAW_DIR, f"{fam}_Viral_Final.faa")))
        bac_ids = get_fasta_ids(os.path.join(RAW_DIR, f"{fam}_Bacterial.faa"))

        v_dom = parse_domtblout(os.path.join(WSL_STAGE_WIN, f"{fam}_Viral.domtblout"))
        b_dom = parse_domtblout(os.path.join(WSL_STAGE_WIN, f"{fam}_Bacterial.domtblout"))
        v_cl = parse_cl0013_multi_hits(os.path.join(WSL_STAGE_WIN, f"{fam}_Viral_CL0013_ga.domtblout"))
        b_cl = parse_cl0013_multi_hits(os.path.join(WSL_STAGE_WIN, f"{fam}_Bacterial_CL0013_ga.domtblout"))

        for dom_label, id_list, d_seed, d_cl in [
            ("Viral", viral_all_ids, v_dom, v_cl),
            ("Bacterial", bac_ids, b_dom, b_cl),
        ]:
            for sid in id_list:
                m = d_seed.get(sid, {})
                cl_hits = d_cl.get(sid, {})
                self_cl = cl_hits.get(fam, {})
                best_cl = max(cl_hits.values(), key=lambda x: x["dom_bit"]) if cl_hits else {}
                pass_self_ga = fam in cl_hits
                pass_pre_ga = any(k in PRE_CL0013 for k in cl_hits)
                pass_nonpre_ga = any(k in NON_PRE_CL0013 for k in cl_hits)
                pass_pf00768_ga = "PF00768" in cl_hits
                pass_dual_13354_00768 = ("PF13354" in cl_hits) and ("PF00768" in cl_hits)

                all_seq_rows.append({
                    "family": fam,
                    "domain_source": dom_label,
                    "seq_id": sid,
                    "in_fimo_r1": (sid in viral_r1_ids) if dom_label == "Viral" else False,
                    "in_final_set": (sid in viral_final_ids) if dom_label == "Viral" else True,
                    "pass_official_self_GA": pass_self_ga,
                    "pass_official_CL0013_PRE_GA": pass_pre_ga,
                    "pass_official_nonPRE_CL0013_GA": pass_nonpre_ga,
                    "pass_official_PF00768_GA": pass_pf00768_ga,
                    "pass_dual_PF13354_PF00768_GA": pass_dual_13354_00768,
                    "official_best_CL0013_acc": best_cl.get("qacc", ""),
                    "official_best_CL0013_name": best_cl.get("qname", ""),
                    "official_self_bitscore": self_cl.get("seq_bit", float("nan")),
                    "official_best_bitscore": best_cl.get("seq_bit", float("nan")),
                    "official_best_coverage": best_cl.get("cov", float("nan")),
                    "hmm_matched": sid in d_seed,
                    "seq_evalue": m.get("seq_evalue", float("nan")),
                    "seq_bitscore": m.get("seq_bitscore", float("nan")),
                    "dom_ievalue": m.get("dom_ievalue", float("nan")),
                    "hmm_coverage": m.get("hmm_coverage", float("nan")),
                    "seq_len": m.get("seq_len", float("nan")),
                    "hmm_len": m.get("hmm_len", float("nan")),
                })

        fam_rows = [r for r in all_seq_rows if r["family"] == fam]
        v_init_df = pd.DataFrame([r for r in fam_rows if r["domain_source"] == "Viral"])
        v_final_df = pd.DataFrame([r for r in fam_rows if r["domain_source"] == "Viral" and r["in_final_set"]])
        b_df = pd.DataFrame([r for r in fam_rows if r["domain_source"] == "Bacterial"])

        for grp_name, sub_df in [
            ("Viral_Initial", v_init_df),
            ("Viral_Final_5701", v_final_df),
            ("Bacterial_Ref", b_df),
        ]:
            n_tot = len(sub_df)
            n_self_ga = int(sub_df["pass_official_self_GA"].sum())
            n_pre_ga = int(sub_df["pass_official_CL0013_PRE_GA"].sum())
            n_nonpre_ga = int(sub_df["pass_official_nonPRE_CL0013_GA"].sum())
            n_pf00768 = int(sub_df["pass_official_PF00768_GA"].sum())
            n_dual = int(sub_df["pass_dual_PF13354_PF00768_GA"].sum())
            summary_rows.append({
                "family": fam,
                "dataset": grp_name,
                "n_sequences": n_tot,
                "n_official_self_GA": n_self_ga,
                "pct_official_self_GA": round(100.0 * n_self_ga / max(n_tot, 1), 2),
                "n_official_CL0013_PRE_GA": n_pre_ga,
                "pct_official_CL0013_PRE_GA": round(100.0 * n_pre_ga / max(n_tot, 1), 2),
                "n_official_PF00768_PBP_GA": n_pf00768,
                "n_dual_PF13354_PF00768_GA": n_dual,
                "n_official_nonPRE_CL0013_GA": n_nonpre_ga,
                "median_official_bitscore": round(float(sub_df["official_best_bitscore"].median()), 2),
                "n_hmm_E1e10": int((sub_df["seq_evalue"] <= 1e-10).sum()),
                "n_hmm_E1e15": int((sub_df["seq_evalue"] <= 1e-15).sum()),
                "n_hmm_E1e20": int((sub_df["seq_evalue"] <= 1e-20).sum()),
                "median_bitscore": round(float(sub_df["seq_bitscore"].median()), 2),
                "median_hmm_coverage": round(float(sub_df["hmm_coverage"].median()), 4),
            })

    df_sum = pd.DataFrame(summary_rows)
    df_all = pd.DataFrame(all_seq_rows)
    df_sum.to_csv(os.path.join(TABLES_DIR, "Table_S01_HMM_Profile_Summary.csv"), index=False)
    df_all.to_csv(os.path.join(SUPP_DIR, "PerSequence_HMM_Scores.tsv"), sep="\t", index=False)
    print("\n=== Official Pfam CL0013 & Seed HMM Verification Summary ===")
    print(df_sum.to_string(index=False))


if __name__ == "__main__":
    main()
