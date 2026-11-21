"""Script 08: Quantitative 3D Structural Localization of Positively Selected vs. Purifying Sites (Comment 6).

Addresses Reviewer Comment 6 (structural basis for assigning positively selected residues to
flexible loops or solvent-accessible surfaces vs. buried catalytic center):
1. Computes residue-level 3D biophysical metrics directly from PDB atomic coordinates across
   experimentally determined bacterial reference crystal structures (1BLS for PF00144,
   1K25 for PF00905, 1BUE for PF13354) and 40 representative ESMFold-predicted viral SBL
   structures per family (123 3D structures total):
   - Relative Solvent Accessibility (RSA = Shrake-Rupley SASA / Tien et al. 2013 max ASA),
   - Backbone Secondary Structure Class (Loop/Coil vs. alpha-Helix / beta-Sheet via phi/psi dihedrals),
   - 3D Euclidean Distance (Angstroms) from each residue's C-alpha atom to the catalytic Serine.
2. Maps 3D biophysical metrics onto all codon alignment sites in PerSite_Selection_Metrics.tsv
   and tests whether MEME positively selected sites (p < 0.05) exhibit significantly higher
   RSA, greater spatial distance from the catalytic Serine, and enrichment in flexible loops/coil
   relative to catalytic motif residues and buried core sites.
"""

import os
import shutil
import subprocess
from pathlib import Path

import sbl_config as cfg
from sbl_config import linux_cmd

DATA_SEQ = cfg.SEQ_DIR
DATA_SEL = cfg.SELECTION_DIR
TABLES_DIR = cfg.TABLES_DIR
SUPP_DIR = cfg.SUPP_DIR
FIG_DIR = cfg.FIG_DIR

WSL_STAGE_WIN, WSL_STAGE_LIN = cfg.stage_dirs("wsl_m05")
WSL_PY = cfg.WSL_PY

WSL_WORKER = r'''
import os
import math
import glob
import numpy as np
import pandas as pd
from scipy import stats
from Bio import SeqIO, Align
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.SASA import ShrakeRupley
from Bio.Data.IUPACData import protein_letters_3to1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STAGE = os.path.dirname(os.path.abspath(__file__))
REPO_LIN = os.path.abspath(os.path.join(STAGE, "..", ".."))
PDB_LIN = f"{REPO_LIN}/data/structures/ProteinStructure"
DOCK_LIN = f"{REPO_LIN}/data/structures/ref_pdbs"

FAMILIES = ["PF00144", "PF00905", "PF13354"]

MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

REF_PDBS = {
    "PF00144": f"{DOCK_LIN}/1bls_P99.pdb",
    "PF00905": f"{DOCK_LIN}/1k25_pneumoniae_PF00905.pdb",
    "PF13354": f"{DOCK_LIN}/1bue_ENTEROBACTER_PF13354.pdb",
}

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def trans_codons(nt):
    return "".join(
        "-" if ("-" in nt[i:i+3] or "N" in nt[i:i+3].upper()) else CODON_TABLE.get(nt[i:i+3].upper(), "X")
        for i in range(0, len(nt) - 2, 3)
    )

def res3to1(resname):
    return protein_letters_3to1.get(resname.strip().capitalize(), "X")

def classify_ss_phi_psi(phi_rad, psi_rad):
    if phi_rad is None or psi_rad is None:
        return "Loop_Coil"
    phi = math.degrees(phi_rad)
    psi = math.degrees(psi_rad)
    if -110.0 <= phi <= -25.0 and -75.0 <= psi <= -5.0:
        return "Helix"
    if (-180.0 <= phi <= -60.0 and (psi >= 80.0 or psi <= -140.0)) or (phi >= 150.0 and psi >= 100.0):
        return "Sheet"
    return "Loop_Coil"

def extract_pdb_metrics(pdb_path):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    model = struct[0]
    best_chain = None
    best_res = []
    for ch in model:
        res = [r for r in ch.get_residues() if r.id[0] == " " and "CA" in r]
        if len(res) > len(best_res):
            best_res = res
            best_chain = ch
    if not best_res or len(best_res) < 50:
        return None

    sr = ShrakeRupley()
    sr.compute(best_chain, level="R")

    ppb = PPBuilder()
    phi_psi_map = {}
    for pp in ppb.build_peptides(best_chain):
        for res_obj, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            phi_psi_map[res_obj.id] = (phi, psi)

    rows = []
    for r in best_res:
        aa = res3to1(r.resname)
        sasa = float(getattr(r, "sasa", 0.0))
        rsa = min(1.0, sasa / MAX_ASA.get(aa, 180.0))
        phi, psi = phi_psi_map.get(r.id, (None, None))
        ss = classify_ss_phi_psi(phi, psi)
        rows.append({
            "aa": aa,
            "rsa": rsa,
            "ss": ss,
            "is_loop": 1.0 if ss == "Loop_Coil" else 0.0,
            "ca": r["CA"].get_coord(),
        })
    return rows

def find_catalytic_ser_in_pdb(pdb_rows):
    """Locate catalytic Serine in SXXK motif inside the 3D PDB sequence."""
    import re
    seq = "".join(r["aa"] for r in pdb_rows)
    # Match canonical S[TVIAFGM][STFVIKL]K or any S..K
    for pat in [r"S[TVIAFGM][STFVIKLA]K", r"S[A-Z]{2}K"]:
        m = re.search(pat, seq)
        if m:
            return m.start()
    return None

def map_pdb_to_msa(pdb_rows, msa_seq_dict):
    pdb_seq = "".join(r["aa"] for r in pdb_rows).replace("X", "A")
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    pdb_kmers = {pdb_seq[i:i+5] for i in range(len(pdb_seq) - 4)}
    cands = []
    for sid, m_seq in msa_seq_dict.items():
        ungap = m_seq.replace("-", "").replace("X", "A")
        u_kmers = {ungap[i:i+5] for i in range(len(ungap) - 4)}
        cands.append((len(pdb_kmers & u_kmers), sid, m_seq, ungap))
    cands.sort(key=lambda x: -x[0])

    best_aln, best_mseq, best_sc = None, None, -1e9
    for _, sid, m_seq, ungap in cands[:3]:
        alns = aligner.align(ungap, pdb_seq)
        if alns.score > best_sc:
            best_sc = alns.score
            best_aln = alns[0]
            best_mseq = m_seq
    if best_aln is None:
        return None

    ungap_to_col = [i + 1 for i, ch in enumerate(best_mseq) if ch != "-"]
    msa_to_pdb = {}
    for (u0, u1), (p0, p1) in zip(best_aln.aligned[0], best_aln.aligned[1]):
        for u_i, p_i in zip(range(u0, u1), range(p0, p1)):
            if u_i < len(ungap_to_col) and p_i < len(pdb_rows):
                msa_to_pdb[ungap_to_col[u_i]] = p_i

    cat_idx = find_catalytic_ser_in_pdb(pdb_rows)
    if cat_idx is None:
        return None
    cat_ca = pdb_rows[cat_idx]["ca"]

    mapped = {}
    for col, p_i in msa_to_pdb.items():
        pr = pdb_rows[p_i]
        mapped[col] = {
            "rsa": pr["rsa"],
            "is_loop": pr["is_loop"],
            "dist_cat_ser_A": float(np.linalg.norm(pr["ca"] - cat_ca)),
        }
    return mapped

def main():
    site_df = pd.read_csv(f"{STAGE}/PerSite_Selection_Metrics.tsv", sep="\t")
    uvig_to_fam = {}
    with open(f"{STAGE}/uvig_to_family.tsv") as f:
        for l in f:
            uv, fm = l.strip().split("\t")
            uvig_to_fam[uv] = fm

    all_pdbs = sorted(glob.glob(f"{PDB_LIN}/*.pdb"))
    vir_pdbs_by_fam = {f: [] for f in FAMILIES}
    for p in all_pdbs:
        bn = os.path.basename(p)
        if bn.startswith("IMGVR"):
            uv = bn.split("__")[0]
            fm = uvig_to_fam.get(uv)
            if fm in vir_pdbs_by_fam:
                vir_pdbs_by_fam[fm].append(p)

    site_records = []
    summary_rows = []

    for fam in FAMILIES:
        fl = fam.lower()
        sub_sites = site_df[site_df["family"] == fam].copy()
        msa_recs = {
            r.id: trans_codons(str(r.seq))
            for r in SeqIO.parse(f"{STAGE}/{fl}.codon_trim.fna", "fasta")
        }

        v_list = vir_pdbs_by_fam[fam]
        step = max(1, len(v_list) // 40)
        selected_pdbs = [("Ref_Crystal", REF_PDBS[fam])] + [("Viral_ESMFold", p) for p in v_list[::step][:40]]

        site_rsa = {s: [] for s in sub_sites["Site"]}
        site_loop = {s: [] for s in sub_sites["Site"]}
        site_dist = {s: [] for s in sub_sites["Site"]}
        ref_rsa, ref_loop, ref_dist = {}, {}, {}

        n_mapped = 0
        for kind, pdb_p in selected_pdbs:
            if not os.path.exists(pdb_p):
                continue
            prows = extract_pdb_metrics(pdb_p)
            if not prows:
                continue
            mapped = map_pdb_to_msa(prows, msa_recs)
            if not mapped or len(mapped) < 80:
                continue
            n_mapped += 1
            for col, vals in mapped.items():
                if col in site_rsa:
                    site_rsa[col].append(vals["rsa"])
                    site_loop[col].append(vals["is_loop"])
                    site_dist[col].append(vals["dist_cat_ser_A"])
                    if kind == "Ref_Crystal":
                        ref_rsa[col] = vals["rsa"]
                        ref_loop[col] = vals["is_loop"]
                        ref_dist[col] = vals["dist_cat_ser_A"]

        print(f"{fam}: mapped {n_mapped} 3D PDB structures onto {len(sub_sites)} codon sites.", flush=True)

        for _, r in sub_sites.iterrows():
            s = int(r["Site"])
            mean_rsa = float(np.mean(site_rsa[s])) if site_rsa[s] else np.nan
            mean_loop = float(np.mean(site_loop[s])) if site_loop[s] else np.nan
            mean_dist = float(np.mean(site_dist[s])) if site_dist[s] else np.nan
            is_motif = bool(r["in_catalytic_motif_core"])
            is_pos = bool(r["p-value"] < 0.05)

            if is_motif:
                cat = "1_Catalytic_Motif"
            elif is_pos:
                cat = "2_MEME_Positive_Selection"
            else:
                cat = "3_Background_NonSelected"

            site_records.append({
                "family": fam,
                "Site": s,
                "site_category": cat,
                "in_catalytic_motif": is_motif,
                "meme_p_value": float(r["p-value"]),
                "meme_pos_p05": is_pos,
                "codon_omega_ng": float(r["codon_omega_ng"]),
                "mean_rsa": mean_rsa,
                "is_solvent_exposed_rsa20": bool(mean_rsa >= 0.20) if pd.notna(mean_rsa) else False,
                "mean_loop_coil_fraction": mean_loop,
                "is_predom_loop_coil": bool(mean_loop >= 0.50) if pd.notna(mean_loop) else False,
                "mean_dist_to_cat_Ser_A": mean_dist,
                "ref_crystal_rsa": ref_rsa.get(s, np.nan),
                "ref_crystal_is_loop": ref_loop.get(s, np.nan),
                "ref_crystal_dist_A": ref_dist.get(s, np.nan),
                "n_structs_mapped": len(site_rsa[s]),
            })

    df_3d = pd.DataFrame(site_records).dropna(subset=["mean_rsa", "mean_dist_to_cat_Ser_A"])
    df_3d.to_csv(f"{STAGE}/PerSite_Selection_and_3D_Biophysics.tsv", sep="\t", index=False)

    # Also compare Active-Center Proximal (< 10 A from catalytic Ser) vs Distal Surface (RSA >= 0.20 & >= 10 A)
    for fam_lbl, sub in [
        ("PF00144", df_3d[df_3d["family"] == "PF00144"]),
        ("PF00905", df_3d[df_3d["family"] == "PF00905"]),
        ("PF13354", df_3d[df_3d["family"] == "PF13354"]),
        ("All_Pooled", df_3d),
    ]:
        mot = sub[sub["site_category"] == "1_Catalytic_Motif"]
        pos = sub[sub["site_category"] == "2_MEME_Positive_Selection"]
        bkg = sub[sub["site_category"] == "3_Background_NonSelected"]

        _, p_rsa_pos_mot = stats.mannwhitneyu(pos["mean_rsa"], mot["mean_rsa"], alternative="greater") if (len(pos) and len(mot)) else (0, np.nan)
        _, p_dist_pos_mot = stats.mannwhitneyu(pos["mean_dist_to_cat_Ser_A"], mot["mean_dist_to_cat_Ser_A"], alternative="greater") if (len(pos) and len(mot)) else (0, np.nan)
        _, p_loop_pos_mot = stats.mannwhitneyu(pos["mean_loop_coil_fraction"], mot["mean_loop_coil_fraction"], alternative="greater") if (len(pos) and len(mot)) else (0, np.nan)
        _, p_rsa_pos_bkg = stats.mannwhitneyu(pos["mean_rsa"], bkg["mean_rsa"], alternative="greater") if (len(pos) and len(bkg)) else (0, np.nan)
        _, p_loop_pos_bkg = stats.mannwhitneyu(pos["mean_loop_coil_fraction"], bkg["mean_loop_coil_fraction"], alternative="greater") if (len(pos) and len(bkg)) else (0, np.nan)

        e_pos = int(pos["is_solvent_exposed_rsa20"].sum())
        e_mot = int(mot["is_solvent_exposed_rsa20"].sum())
        _, p_exp_fisher = stats.fisher_exact([[e_pos, len(pos) - e_pos], [e_mot, len(mot) - e_mot]], alternative="greater") if (len(pos) and len(mot)) else (0, np.nan)

        # Also check 3D active-pocket (<12 A from cat Ser) vs distal outer shell (>=12 A) MEME positive rate
        prox = sub[sub["mean_dist_to_cat_Ser_A"] < 12.0]
        dist_sh = sub[sub["mean_dist_to_cat_Ser_A"] >= 12.0]
        pos_prox = int(prox["meme_pos_p05"].sum())
        pos_dist = int(dist_sh["meme_pos_p05"].sum())
        or_shell, p_shell = stats.fisher_exact([[pos_dist, len(dist_sh) - pos_dist], [pos_prox, len(prox) - pos_prox]], alternative="greater") if (len(prox) and len(dist_sh)) else (np.nan, np.nan)

        for grp_name, gdf in [("Catalytic_Motifs", mot), ("MEME_Positive_Selection", pos), ("Background_NonSelected", bkg)]:
            summary_rows.append({
                "family": fam_lbl,
                "site_group": grp_name,
                "n_sites": len(gdf),
                "median_RSA": round(float(gdf["mean_rsa"].median()), 4),
                "iqr_RSA": f"{gdf['mean_rsa'].quantile(0.25):.3f}-{gdf['mean_rsa'].quantile(0.75):.3f}",
                "pct_solvent_exposed_RSA20": round(100 * float(gdf["is_solvent_exposed_rsa20"].mean()), 2) if len(gdf) else np.nan,
                "median_loop_coil_fraction": round(float(gdf["mean_loop_coil_fraction"].median()), 4),
                "pct_predom_loop_coil": round(100 * float(gdf["is_predom_loop_coil"].mean()), 2) if len(gdf) else np.nan,
                "median_dist_to_cat_Ser_A": round(float(gdf["mean_dist_to_cat_Ser_A"].median()), 2),
                "iqr_dist_A": f"{gdf['mean_dist_to_cat_Ser_A'].quantile(0.25):.1f}-{gdf['mean_dist_to_cat_Ser_A'].quantile(0.75):.1f}",
                "MW_p_RSA_Pos_gt_Motif": p_rsa_pos_mot if grp_name == "MEME_Positive_Selection" else np.nan,
                "MW_p_Dist_Pos_gt_Motif": p_dist_pos_mot if grp_name == "MEME_Positive_Selection" else np.nan,
                "MW_p_Loop_Pos_gt_Motif": p_loop_pos_mot if grp_name == "MEME_Positive_Selection" else np.nan,
                "MW_p_RSA_Pos_gt_Bkg": p_rsa_pos_bkg if grp_name == "MEME_Positive_Selection" else np.nan,
                "MW_p_Loop_Pos_gt_Bkg": p_loop_pos_bkg if grp_name == "MEME_Positive_Selection" else np.nan,
                "Fisher_p_Exposed_Pos_gt_Motif": p_exp_fisher if grp_name == "MEME_Positive_Selection" else np.nan,
                "Distal_gt12A_vs_Proximal_Pos_OR": round(float(or_shell), 3) if grp_name == "MEME_Positive_Selection" and pd.notna(or_shell) else np.nan,
                "Distal_gt12A_vs_Proximal_Fisher_p": float(p_shell) if grp_name == "MEME_Positive_Selection" and pd.notna(p_shell) else np.nan,
            })

    df_sum = pd.DataFrame(summary_rows)
    df_sum.to_csv(f"{STAGE}/Table_S11_Structural_Assessment_Selected_Sites.csv", index=False)
    print("\n=== Quantitative 3D Structural Assessment Summary ===")
    print(df_sum.to_string(index=False))

    # Plot Figure S04
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), dpi=300)
    cats = ["1_Catalytic_Motif", "2_MEME_Positive_Selection", "3_Background_NonSelected"]
    lbls = ["Catalytic\nMotifs", "MEME Pos.\nSelected", "Background\nNon-Selected"]
    cols = ["#E64A19", "#1B9E77", "#90A4AE"]

    ax = axes[0]
    bp = ax.boxplot([df_3d[df_3d["site_category"] == c]["mean_rsa"].values for c in cats],
                    tick_labels=lbls, patch_artist=True, showfliers=False)
    for pt, c in zip(bp["boxes"], cols):
        pt.set_facecolor(c)
        pt.set_alpha(0.85)
    ax.axhline(0.20, color="#333333", linestyle="--", linewidth=0.9, label="Solvent-Exposed Cutoff (RSA=0.20)")
    ax.set_ylabel("Relative Solvent Accessibility (RSA)", fontweight="bold")
    ax.set_title("(A) Solvent Exposure of Selected Sites", fontweight="bold", loc="left")
    ax.legend(fontsize=8)

    ax = axes[1]
    bp = ax.boxplot([df_3d[df_3d["site_category"] == c]["mean_dist_to_cat_Ser_A"].values for c in cats],
                    tick_labels=lbls, patch_artist=True, showfliers=False)
    for pt, c in zip(bp["boxes"], cols):
        pt.set_facecolor(c)
        pt.set_alpha(0.85)
    ax.set_ylabel("3D Distance to Catalytic Serine Ca (Angstrom)", fontweight="bold")
    ax.set_title("(B) Spatial Distance from Active Center", fontweight="bold", loc="left")

    ax = axes[2]
    exp_pcts = [100 * df_3d[df_3d["site_category"] == c]["is_solvent_exposed_rsa20"].mean() for c in cats]
    loop_pcts = [100 * df_3d[df_3d["site_category"] == c]["is_predom_loop_coil"].mean() for c in cats]
    x = np.arange(3)
    w = 0.35
    ax.bar(x - w/2, exp_pcts, w, label="Solvent-Exposed (RSA >= 0.20)", color="#2B5C8F")
    ax.bar(x + w/2, loop_pcts, w, label="Flexible Loop/Coil Structure", color="#FFB74D")
    ax.set_xticks(x)
    ax.set_xticklabels(lbls, fontweight="bold")
    ax.set_ylabel("Proportion of Sites (%)", fontweight="bold")
    ax.set_title("(C) Surface & Loop Enrichment", fontweight="bold", loc="left")
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{STAGE}/Fig_S04_Structural_RSA_Distance_Selection.pdf", bbox_inches="tight")
    plt.savefig(f"{STAGE}/Fig_S04_Structural_RSA_Distance_Selection.png", bbox_inches="tight")
    plt.close()

main()
'''

def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(SUPP_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(WSL_STAGE_WIN, exist_ok=True)
    cfg.require(cfg.PDB_DIR, cfg.PDB_HINT)

    shutil.copy2(SUPP_DIR / "PerSite_Selection_Metrics.tsv", WSL_STAGE_WIN / "PerSite_Selection_Metrics.tsv")
    uv_map = []
    for fam in ["PF00144", "PF00905", "PF13354"]:
        fl = fam.lower()
        shutil.copy2(DATA_SEL / f"{fl}.codon_trim.fna", WSL_STAGE_WIN / f"{fl}.codon_trim.fna")
        with open(DATA_SEQ / f"{fam}_Viral_Final.faa", encoding="utf-8", errors="ignore") as f:
            for l in f:
                if l.startswith(">"):
                    uv = l[1:].strip().split("|")[0]
                    uv_map.append(f"{uv}\t{fam}\n")
    (WSL_STAGE_WIN / "uvig_to_family.tsv").write_text("".join(sorted(set(uv_map))), encoding="utf-8")

    worker_path = WSL_STAGE_WIN / "run_m05_rsa.py"
    worker_path.write_text(WSL_WORKER, encoding="utf-8")

    proc = subprocess.run(linux_cmd(WSL_PY, f"{WSL_STAGE_LIN}/run_m05_rsa.py"), capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(f"Script 08 failed:\n{proc.stderr}")

    shutil.copy2(WSL_STAGE_WIN / "Table_S11_Structural_Assessment_Selected_Sites.csv",
                 TABLES_DIR / "Table_S11_Structural_Assessment_Selected_Sites.csv")
    shutil.copy2(WSL_STAGE_WIN / "PerSite_Selection_and_3D_Biophysics.tsv",
                 SUPP_DIR / "PerSite_Selection_and_3D_Biophysics.tsv")
    shutil.copy2(WSL_STAGE_WIN / "Fig_S04_Structural_RSA_Distance_Selection.pdf",
                 FIG_DIR / "Fig_S04_Structural_RSA_Distance_Selection.pdf")
    shutil.copy2(WSL_STAGE_WIN / "Fig_S04_Structural_RSA_Distance_Selection.png",
                 FIG_DIR / "Fig_S04_Structural_RSA_Distance_Selection.png")


if __name__ == "__main__":
    main()
