"""Official Pfam --cut_ga + unified MEME/FIMO threshold sensitivity.

Evaluates viral and bacterial SBL candidates under strictly
identical MEME/FIMO thresholds (p <= 1e-5 and p <= 1e-10 across Round 1 and Round 2)
combined with official Pfam --cut_ga gathering thresholds (without imposing an artificial
hard-coded SXXK regex pre-filter, relying instead on motif-PWM filtering to eliminate
non-catalytic false positives), alongside Clan CL0013 penicillin-recognizing enzyme retention.
"""

from pathlib import Path
import os
import re
import subprocess
import pandas as pd

import sbl_config as cfg
from sbl_config import linux_cmd

RAW_DIR = str(cfg.SEQ_DIR)
TABLES_DIR = str(cfg.TABLES_DIR)
SUPP_DIR = str(cfg.SUPP_DIR)
_stage, WSL_STAGE_LIN = cfg.stage_dirs("wsl_m01")
WSL_STAGE_WIN = str(_stage)
MEME, FIMO = cfg.MEME, cfg.FIMO

FAMILIES = ["PF00144", "PF00905", "PF13354"]


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


def parse_fimo_best_pvals(tsv_path):
    """Return dict: seq_id -> (min_pvalue, min_qvalue, max_score, n_motif_hits)."""
    stats = {}
    with open(tsv_path, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline()
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            motif_id = parts[0]
            seq_id = parts[2]
            try:
                score = float(parts[6])
                pval = float(parts[7])
                qval = float(parts[8]) if parts[8] else 1.0
            except ValueError:
                continue
            if seq_id not in stats:
                stats[seq_id] = {
                    "min_p": pval,
                    "min_q": qval,
                    "max_score": score,
                    "motifs": {motif_id},
                }
            else:
                st = stats[seq_id]
                if pval < st["min_p"]:
                    st["min_p"] = pval
                if qval < st["min_q"]:
                    st["min_q"] = qval
                if score > st["max_score"]:
                    st["max_score"] = score
                st["motifs"].add(motif_id)
    return stats


def has_canonical_sxxk(seq):
    return bool(re.search(r"S[A-Z]{2}K", seq))


def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(SUPP_DIR, exist_ok=True)
    os.makedirs(WSL_STAGE_WIN, exist_ok=True)
    hmm_df = pd.read_csv(os.path.join(SUPP_DIR, "PerSequence_HMM_Scores.tsv"), sep="\t")
    summary_rows = []
    per_seq_rows = []

    for fam in FAMILIES:
        r1_v_tsv = os.path.join(WSL_STAGE_WIN, f"fimo_r1_viral_{fam}", "fimo.tsv")
        r1_b_tsv = os.path.join(WSL_STAGE_WIN, f"fimo_r1_bac_{fam}", "fimo.tsv")
        r2_v_tsv = os.path.join(WSL_STAGE_WIN, f"fimo_r2_viral_{fam}", "fimo.tsv")
        r2_b_tsv = os.path.join(WSL_STAGE_WIN, f"fimo_r2_bac_{fam}", "fimo.tsv")

        v_recs = read_fasta(os.path.join(RAW_DIR, f"{fam}_Viral.faa"))
        b_recs = read_fasta(os.path.join(RAW_DIR, f"{fam}_Bacterial.faa"))

        if os.path.exists(r1_v_tsv) and os.path.exists(r1_b_tsv):
            print(f"Reusing cached MEME/FIMO Round 1 results for {fam}...")
        else:
            print(f"Running MEME/FIMO Round 1 for {fam}...")
            seed_recs = read_fasta(os.path.join(RAW_DIR, f"{fam}_seed.faa"))
            seed_clean = os.path.join(WSL_STAGE_WIN, f"{fam}_seed.clean.faa")
            with open(seed_clean, "w", encoding="utf-8") as f:
                for sid, hdr, seq in seed_recs:
                    f.write(f">{sid}\n{seq}\n")

            for kind, recs in [("Viral_Initial", v_recs), ("Bacterial_Ref", b_recs)]:
                out_fa = os.path.join(WSL_STAGE_WIN, f"{fam}_{kind}.clean.faa")
                with open(out_fa, "w", encoding="utf-8") as f:
                    for sid, hdr, seq in recs:
                        f.write(f">{sid}\n{seq}\n")

            cmd_r1 = (
                f"cd {WSL_STAGE_LIN} && "
                f"{MEME} {fam}_seed.clean.faa -protein -mod zoops -nmotifs 5 -minw 6 -maxw 20 -maxsize 1000000 -oc meme_r1_{fam} && "
                f"{FIMO} --thresh 1e-4 --oc fimo_r1_viral_{fam} meme_r1_{fam}/meme.txt {fam}_Viral_Initial.clean.faa && "
                f"{FIMO} --thresh 1e-4 --oc fimo_r1_bac_{fam} meme_r1_{fam}/meme.txt {fam}_Bacterial_Ref.clean.faa"
            )
            run_wsl_cmd(cmd_r1)

        r1_viral = parse_fimo_best_pvals(r1_v_tsv)
        r1_bac = parse_fimo_best_pvals(r1_b_tsv)

        if os.path.exists(r2_v_tsv) and os.path.exists(r2_b_tsv):
            print(f"Reusing cached MEME/FIMO Round 2 results for {fam}...")
        else:
            print(f"Running MEME/FIMO Round 2 for {fam}...")
            r1_pass_v = [r for r in v_recs if r1_viral.get(r[0], {}).get("min_p", 1.0) <= 1e-5]
            r1_pass_b = [r for r in b_recs if r1_bac.get(r[0], {}).get("min_p", 1.0) <= 1e-5]

            step_v = max(1, len(r1_pass_v) // 80)
            step_b = max(1, len(r1_pass_b) // 40)
            r2_train = r1_pass_v[::step_v][:80] + r1_pass_b[::step_b][:40]
            r2_train_path = os.path.join(WSL_STAGE_WIN, f"{fam}_r2_train.clean.faa")
            with open(r2_train_path, "w", encoding="utf-8") as f:
                seen = set()
                for sid, hdr, seq in r2_train:
                    safe_id = sid.replace("|", "_")[:50]
                    if safe_id not in seen and len(seq) >= 50:
                        seen.add(safe_id)
                        f.write(f">{safe_id}\n{seq}\n")

            cmd_r2 = (
                f"cd {WSL_STAGE_LIN} && "
                f"{MEME} {fam}_r2_train.clean.faa -protein -mod zoops -nmotifs 5 -minw 6 -maxw 20 -maxsize 1000000 -oc meme_r2_{fam} && "
                f"{FIMO} --thresh 1e-4 --oc fimo_r2_viral_{fam} meme_r2_{fam}/meme.txt {fam}_Viral_Initial.clean.faa && "
                f"{FIMO} --thresh 1e-4 --oc fimo_r2_bac_{fam} meme_r2_{fam}/meme.txt {fam}_Bacterial_Ref.clean.faa"
            )
            run_wsl_cmd(cmd_r2)

        r2_viral = parse_fimo_best_pvals(r2_v_tsv)
        r2_bac = parse_fimo_best_pvals(r2_b_tsv)

        final_v_ids = {r[0] for r in read_fasta(os.path.join(RAW_DIR, f"{fam}_Viral_Final.faa"))}
        r1_v_ids = {r[0] for r in read_fasta(os.path.join(RAW_DIR, f"{fam}_Viral_FIMO_r1.faa"))}

        hmm_sub = hmm_df[hmm_df["family"] == fam].drop_duplicates(subset=["domain_source", "seq_id"]).set_index(["domain_source", "seq_id"])

        for dom_src, recs, d_r1, d_r2 in [
            ("Viral_Initial", v_recs, r1_viral, r2_viral),
            ("Viral_Final_5701", [r for r in v_recs if r[0] in final_v_ids], r1_viral, r2_viral),
            ("Bacterial_Ref", b_recs, r1_bac, r2_bac),
        ]:
            hmm_key = "Viral" if dom_src.startswith("Viral") else "Bacterial"
            n_tot = len(recs)
            c_r1_1e5 = 0
            c_r1_1e10 = 0
            c_r2_1e5 = 0
            c_r2_1e10 = 0
            c_sxxk = 0
            c_self_ga_r2_1e5 = 0
            c_self_ga_r2_1e10 = 0
            c_cl_ga_r2_1e5 = 0
            c_cl_ga_r2_1e10 = 0

            for sid, hdr, seq in recs:
                p1 = d_r1.get(sid, {}).get("min_p", 1.0)
                p2 = d_r2.get(sid, {}).get("min_p", 1.0)
                m2 = len(d_r2.get(sid, {}).get("motifs", set()))
                sxxk = has_canonical_sxxk(seq)
                if (hmm_key, sid) in hmm_sub.index:
                    hrow = hmm_sub.loc[(hmm_key, sid)]
                    hmm_e = float(hrow["seq_evalue"])
                    pass_self_ga = bool(hrow["pass_official_self_GA"])
                    pass_cl_ga = bool(hrow["pass_official_CL0013_PRE_GA"])
                else:
                    hmm_e = 10.0
                    pass_self_ga = False
                    pass_cl_ga = False

                if p1 <= 1e-5:
                    c_r1_1e5 += 1
                if p1 <= 1e-10:
                    c_r1_1e10 += 1
                if p2 <= 1e-5:
                    c_r2_1e5 += 1
                if p2 <= 1e-10:
                    c_r2_1e10 += 1
                if sxxk:
                    c_sxxk += 1
                if pass_self_ga and p2 <= 1e-5:
                    c_self_ga_r2_1e5 += 1
                if pass_self_ga and p2 <= 1e-10:
                    c_self_ga_r2_1e10 += 1
                if pass_cl_ga and p2 <= 1e-5:
                    c_cl_ga_r2_1e5 += 1
                if pass_cl_ga and p2 <= 1e-10:
                    c_cl_ga_r2_1e10 += 1

                if dom_src in ("Viral_Initial", "Bacterial_Ref"):
                    per_seq_rows.append({
                        "family": fam,
                        "domain_source": hmm_key,
                        "seq_id": sid,
                        "in_viral_r1_orig": sid in r1_v_ids if hmm_key == "Viral" else False,
                        "in_viral_final_5701": sid in final_v_ids if hmm_key == "Viral" else True,
                        "pass_official_self_GA": pass_self_ga,
                        "pass_official_CL0013_PRE_GA": pass_cl_ga,
                        "hmm_evalue": hmm_e,
                        "fimo_r1_min_p": p1,
                        "fimo_r2_min_p": p2,
                        "fimo_r2_n_motifs": m2,
                        "has_SXXK_descriptive": sxxk,
                        "pass_unified_1e5_r1_or_r2": (p1 <= 1e-5) or (p2 <= 1e-5),
                        "pass_unified_1e10_r1_or_r2": (p1 <= 1e-10) or (p2 <= 1e-10),
                    })

            summary_rows.append({
                "family": fam,
                "dataset": dom_src,
                "total_n": n_tot,
                "fimo_r1_p1e5": c_r1_1e5,
                "fimo_r1_p1e5_pct": round(100 * c_r1_1e5 / max(n_tot, 1), 2),
                "fimo_r1_p1e10": c_r1_1e10,
                "fimo_r1_p1e10_pct": round(100 * c_r1_1e10 / max(n_tot, 1), 2),
                "fimo_r2_p1e5": c_r2_1e5,
                "fimo_r2_p1e5_pct": round(100 * c_r2_1e5 / max(n_tot, 1), 2),
                "fimo_r2_p1e10": c_r2_1e10,
                "fimo_r2_p1e10_pct": round(100 * c_r2_1e10 / max(n_tot, 1), 2),
                "has_SXXK_n": c_sxxk,
                "has_SXXK_pct": round(100 * c_sxxk / max(n_tot, 1), 2),
                "official_self_GA_and_r2_p1e5": c_self_ga_r2_1e5,
                "official_self_GA_and_r2_p1e10": c_self_ga_r2_1e10,
                "official_CL0013_GA_and_r2_p1e5": c_cl_ga_r2_1e5,
                "official_CL0013_GA_and_r2_p1e10": c_cl_ga_r2_1e10,
            })

    df_sum = pd.DataFrame(summary_rows)
    df_seq = pd.DataFrame(per_seq_rows)
    df_sum.to_csv(os.path.join(TABLES_DIR, "Table_S02_Filtering_Threshold_Sensitivity.csv"), index=False)
    df_seq.to_csv(os.path.join(SUPP_DIR, "PerSequence_FIMO_Sensitivity.tsv"), sep="\t", index=False)
    print("\n=== Official Pfam --cut_ga + Unified MEME/FIMO Threshold Sensitivity Summary ===")
    print(df_sum.to_string(index=False))


if __name__ == "__main__":
    main()
