"""PF00144 MEME diagnostics and statistical tests of purifying selection at catalytic motifs.

1. MEME diagnostics:
   - Quantifies threshold sensitivity (p < 0.10 [41.87%] vs p < 0.05 [29.76%] vs p < 0.01 vs BH-FDR q < 0.05),
   - Quantifies tree size/power (1,382 vs 916 vs 252 tested branches; Newick tree total branch length and
     occupancy-weighted effective site branch length),
   - Quantifies sequence divergence (pairwise AA p-distance), alignment quality (column Shannon
     entropy H_aa, occupancy, and SXXK motif alignment shift vs hmmalign domain-anchored MSAs),
   - Tests intragenic homoplasy/recombination via the Bruen et al. (2006) Pairwise Homoplasy Index (PHI).
2. Purifying selection & depletion tests:
   - Performs Fisher's exact, exact binomial, and 10,000-iteration permutation tests evaluating depletion
     of MEME positively selected sites within catalytic motifs (Motifs 1-3) relative to background codons.
   - Tests codon-level purifying selection at catalytic motifs from domain-anchored
     (hmmalign back-translated) Nei-Gojobori codon dN/dS (omega_domain), amino-acid Shannon entropy
     (H_aa_domain), and MEME purifying-regime parameters (omega^- = beta^1 / alpha, p^1).
"""

import os
import shutil
import subprocess
from pathlib import Path

import sbl_config as cfg
from sbl_config import linux_cmd

DATA_SEL = cfg.SELECTION_DIR
DATA_HMM = cfg.HMM_REBUILT_DIR  # seed HMMs rebuilt by 01_build_and_verify_hmm.py (used for hmmalign)
TABLES_DIR = cfg.TABLES_DIR
SUPP_DIR = cfg.SUPP_DIR
FIG_DIR = cfg.FIG_DIR

SEL_SRC = cfg.SELECTION_DIR
SP_SRC = cfg.SELECTION_DIR

WSL_STAGE_WIN, WSL_STAGE_LIN = cfg.stage_dirs("wsl_m04_fast")
WSL_PY = cfg.WSL_PY

WSL_WORKER = r'''
import os
import json
import math
import re
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
from Bio import SeqIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STAGE = os.path.dirname(os.path.abspath(__file__))
FAMILIES = ["PF00144", "PF00905", "PF13354"]

# Core catalytic motif positions in codon_trim.fna (1-based)
CORE_MOTIFS = {
    "PF00144": list(range(61, 65)) + list(range(120, 125)) + list(range(218, 221)),
    "PF00905": list(range(64, 68)) + list(range(143, 147)) + list(range(231, 234)),
    "PF13354": list(range(75, 79)) + list(range(134, 138)) + list(range(239, 242)),
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

def build_codon_ns_table():
    """Precompute Nei-Gojobori (1986) nonsynonymous/synonymous site counts per sense codon."""
    tbl = {}
    for c, aa in CODON_TABLE.items():
        if aa == "*":
            continue
        s_cnt, n_cnt = 0.0, 0.0
        for pos in range(3):
            for b in ("A", "C", "G", "T"):
                if b == c[pos]:
                    continue
                c2 = c[:pos] + b + c[pos+1:]
                aa2 = CODON_TABLE[c2]
                if aa2 == "*":
                    continue
                if aa2 == aa:
                    s_cnt += 1.0 / 3.0
                else:
                    n_cnt += 1.0 / 3.0
        tbl[c] = (n_cnt, s_cnt)
    return tbl

CODON_NS = build_codon_ns_table()

def nei_gojobori_column_omega(col_codons, min_seqs=10):
    """Frequency-weighted Nei-Gojobori pairwise dN/dS (omega) for one codon alignment column."""
    from collections import Counter
    valid = [c for c in col_codons if c in CODON_NS]
    if len(valid) < min_seqs:
        return np.nan
    cnt = Counter(valid)
    items = list(cnt.items())
    tot_nd, tot_sd, tot_en, tot_es = 0.0, 0.0, 0.0, 0.0
    for i in range(len(items)):
        c1, w1 = items[i]
        n1, s1 = CODON_NS[c1]
        tot_en += w1 * (w1 - 1) * n1
        tot_es += w1 * (w1 - 1) * s1
        for k in range(i + 1, len(items)):
            c2, w2 = items[k]
            n2, s2 = CODON_NS[c2]
            w = 2.0 * w1 * w2
            tot_en += w * 0.5 * (n1 + n2)
            tot_es += w * 0.5 * (s1 + s2)
            diffs = sum(1 for p in range(3) if c1[p] != c2[p])
            if CODON_TABLE[c1] == CODON_TABLE[c2]:
                tot_sd += w * diffs
            else:
                tot_nd += w * (1.0 + (diffs - 1) * 0.75)
                tot_sd += w * ((diffs - 1) * 0.25)
    dn = tot_nd / max(tot_en, 1e-6)
    ds = tot_sd / max(tot_es, 1e-6)
    if ds < 1e-4:
        return np.nan
    return float(min(dn / ds, 10.0))

def shannon_entropy(chars):
    valid = [c for c in chars if c not in ("-", "X", "*")]
    if not valid:
        return 0.0
    n = len(valid)
    counts = {}
    for c in valid:
        counts[c] = counts.get(c, 0) + 1
    return -sum((v / n) * math.log2(v / n) for v in counts.values())

def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0.0, 1.0)
    return out

def compute_phi_homoplasy(nt_seqs, n_perm=500, seed=42):
    """Bruen et al. (2006) Pairwise Homoplasy Index (PHI) permutation test."""
    rng = np.random.default_rng(seed)
    idx = np.linspace(0, len(nt_seqs) - 1, min(80, len(nt_seqs)), dtype=int)
    seqs = [nt_seqs[i].upper() for i in idx]
    L = len(seqs[0])
    inf_cols = []
    for col_i in range(L):
        col = [s[col_i] for s in seqs]
        counts = {}
        for ch in col:
            if ch in ("A", "C", "G", "T"):
                counts[ch] = counts.get(ch, 0) + 1
        if sum(1 for v in counts.values() if v >= 2) >= 2:
            inf_cols.append(col)
    m = len(inf_cols)
    if m < 10:
        return np.nan, np.nan, np.nan, m
    if m > 220:
        sub_idx = np.linspace(0, m - 1, 220, dtype=int)
        inf_cols = [inf_cols[k] for k in sub_idx]
        m = len(inf_cols)

    top2_masks = []
    for col in inf_cols:
        counts = {}
        for i, ch in enumerate(col):
            if ch in ("A", "C", "G", "T"):
                counts.setdefault(ch, set()).add(i)
        ordered = sorted(counts.values(), key=lambda s: -len(s))
        top2_masks.append((ordered[0], ordered[1]))

    full_mat = np.zeros((m, m), dtype=float)
    for j1 in range(m):
        a0, a1 = top2_masks[j1]
        for j2 in range(j1 + 1, m):
            b0, b1 = top2_masks[j2]
            gam = int(bool(a0 & b0)) + int(bool(a0 & b1)) + int(bool(a1 & b0)) + int(bool(a1 & b1))
            if gam == 4:
                full_mat[j1, j2] = 1.0
                full_mat[j2, j1] = 1.0

    w = max(5, min(50, m // 3))
    def phi_of_order(order):
        sub = full_mat[np.ix_(order, order)]
        s = 0.0
        c = 0
        for d in range(1, w + 1):
            diag = np.diagonal(sub, offset=d)
            s += diag.sum()
            c += len(diag)
        return s / max(c, 1)

    obs_phi = phi_of_order(np.arange(m))
    base_ord = np.arange(m)
    null_phis = np.array([phi_of_order(rng.permutation(base_ord)) for _ in range(n_perm)])
    p_phi = float((np.sum(null_phis <= obs_phi) + 1) / (n_perm + 1))
    return float(obs_phi), float(np.mean(null_phis)), p_phi, m

def hmmalign_sxxk_lock_pct(fam):
    """Run hmmalign on full protein set to measure motif lock percentage vs unanchored MSA."""
    subprocess.run(
        f"hmmalign --trim --outformat afa {STAGE}/{fam}.hmm {STAGE}/{fam}_rename.faa > {STAGE}/{fam}_hmm.afa",
        shell=True, check=True
    )
    recs = list(SeqIO.parse(f"{STAGE}/{fam}_hmm.afa", "fasta"))
    ref = str(recs[0].seq)
    mcols = [i for i, ch in enumerate(ref) if not ch.islower() and ch != "."]
    cnt = {}
    for r in recs:
        s = str(r.seq)
        aa = "".join(s[i].upper() if s[i] != "." else "-" for i in mcols)
        for m in re.finditer(r"S[A-Z]{2}K", aa):
            cnt[m.start() + 1] = cnt.get(m.start() + 1, 0) + 1
    top = max(cnt.values()) if cnt else 0
    return round(100.0 * top / max(len(recs), 1), 2), len(recs)

def compute_domain_anchored_metrics(fam, recs):
    """Back-translate hmmalign match columns of codon_trim.fna and compute site dN/dS and entropy."""
    from collections import Counter
    fl = fam.lower()
    ungap_codons = {}
    faa_path = f"{STAGE}/{fl}_codon_ungap.faa"
    afa_path = f"{STAGE}/{fl}_codon_hmm.afa"
    with open(faa_path, "w") as out:
        for r in recs:
            nt = str(r.seq).upper()
            cods = [
                nt[i:i+3] for i in range(0, len(nt) - 2, 3)
                if "-" not in nt[i:i+3] and "N" not in nt[i:i+3] and CODON_TABLE.get(nt[i:i+3], "*") != "*"
            ]
            aa = "".join(CODON_TABLE[c] for c in cods)
            ungap_codons[r.id] = cods
            out.write(f">{r.id}\n{aa}\n")

    subprocess.run(f"hmmalign --outformat afa {STAGE}/{fam}.hmm {faa_path} > {afa_path}", shell=True, check=True)
    arecs = list(SeqIO.parse(afa_path, "fasta"))
    ref = str(arecs[0].seq)
    mcols = [i for i, ch in enumerate(ref) if not ch.islower() and ch != "."]
    mset = set(mcols)
    L = len(mcols)

    cod_mat, aa_mat = [], []
    for r in arecs:
        cods = ungap_codons[r.id]
        idx = 0
        row_c, row_a = [], []
        for i, ch in enumerate(str(r.seq)):
            if ch in ("-", "."):
                if i in mset:
                    row_c.append("---")
                    row_a.append("-")
            else:
                c = cods[idx] if idx < len(cods) else "---"
                idx += 1
                if i in mset:
                    row_c.append(c)
                    row_a.append(ch.upper())
        cod_mat.append(row_c)
        aa_mat.append("".join(row_a))

    core_pos = set()
    ext_pos = set()
    for pat, plen in [(r"S[A-Z]{2}K", 4), (r"[SY][A-Z]N", 3), (r"[KH][TS]G", 3)]:
        cnt = Counter()
        for s in aa_mat:
            for m in re.finditer(pat, s):
                cnt[m.start()] += 1
        top_col = cnt.most_common(1)[0][0]
        for k in range(top_col, top_col + plen):
            core_pos.add(k)
        for k in range(max(0, top_col - 1), min(L, top_col + plen + 1)):
            ext_pos.add(k)

    rows = []
    for j in range(L):
        col_c = [row[j] for row in cod_mat]
        col_a = [row[j] for row in aa_mat]
        rows.append({
            "family": fam,
            "hmm_col": j + 1,
            "in_catalytic_motif_core": j in core_pos,
            "in_catalytic_pocket_ext": j in ext_pos,
            "domain_codon_omega": nei_gojobori_column_omega(col_c, min_seqs=10),
            "domain_entropy_aa": shannon_entropy(col_a),
            "domain_occupancy": sum(1 for c in col_a if c not in ("-", "X")) / max(len(col_a), 1),
        })
    return pd.DataFrame(rows)

def main():
    diag_rows = []
    proof_rows = []
    site_dfs = []
    dom_dfs = []

    for fam in FAMILIES:
        fl = fam.lower()
        recs = list(SeqIO.parse(f"{STAGE}/{fl}.codon_trim.fna", "fasta"))
        nt_seqs = [str(r.seq).upper() for r in recs]
        aa_seqs = [trans_codons(s) for s in nt_seqs]
        n_taxa = len(aa_seqs)
        n_codons = len(aa_seqs[0])

        entropies = [shannon_entropy([a[j] for a in aa_seqs]) for j in range(n_codons)]
        occs = [sum(1 for a in aa_seqs if a[j] not in ("-", "X")) / n_taxa for j in range(n_codons)]
        codon_omegas = [
            nei_gojobori_column_omega([nt[3*j:3*j+3] for nt in nt_seqs], min_seqs=10)
            for j in range(n_codons)
        ]

        # Tree total branch length from Newick file (since MEME.bf v4.1 leaves JSON Total branch length at 0)
        with open(f"{STAGE}/{fl}.nwk") as f:
            nwk_txt = f.read()
        tree_total_bl = float(sum(float(x) for x in re.findall(r":([0-9.eE+-]+)", nwk_txt)))

        # SXXK alignment concentration in published codon_trim.fna
        sxxk_cols = {}
        for a in aa_seqs:
            for m in re.finditer(r"S[A-Z]{2}K", a):
                sxxk_cols[m.start() + 1] = sxxk_cols.get(m.start() + 1, 0) + 1
        top_sxxk_pub = round(100.0 * (max(sxxk_cols.values()) if sxxk_cols else 0) / n_taxa, 2)
        hmm_lock_pct, n_hmm_seqs = hmmalign_sxxk_lock_pct(fam)

        dom_dfs.append(compute_domain_anchored_metrics(fam, recs))

        # Pairwise AA p-distance
        rng = np.random.default_rng(42)
        samp = [aa_seqs[i] for i in rng.choice(n_taxa, min(80, n_taxa), replace=False)]
        pdists = []
        for i in range(len(samp)):
            for k in range(i + 1, len(samp)):
                comp = sum(1 for x, y in zip(samp[i], samp[k]) if x != "-" and y != "-")
                diff = sum(1 for x, y in zip(samp[i], samp[k]) if x != "-" and y != "-" and x != y)
                if comp > 30:
                    pdists.append(diff / comp)

        obs_phi, null_phi, p_phi, n_inf = compute_phi_homoplasy(nt_seqs, n_perm=500)

        # Load HyPhy MEME JSON
        with open(f"{STAGE}/{fl}_meme.json") as f:
            mm = json.load(f)
        headers = [x[0] for x in mm["MLE"]["headers"]]
        mdf = pd.DataFrame(mm["MLE"]["content"]["0"], columns=headers)
        mdf["Site"] = np.arange(1, len(mdf) + 1)
        mdf["family"] = fam
        mdf["q_value"] = bh_fdr(mdf["p-value"].values)
        mdf["entropy_aa"] = entropies
        mdf["occupancy"] = occs
        mdf["tree_total_bl"] = tree_total_bl
        mdf["site_effective_branch_length"] = tree_total_bl * np.asarray(occs)
        mdf["codon_omega_ng"] = codon_omegas
        mdf["meme_omega_minus"] = np.where(
            mdf["&alpha;"] > 1e-3,
            np.clip(mdf["&beta;<sup>1</sup>"] / mdf["&alpha;"], 0.0, 10.0),
            np.nan,
        )
        mdf["n_tested_branches"] = len(mm["branch attributes"]["0"]) if "branch attributes" in mm else (2 * n_taxa - 3)

        core_set = set(CORE_MOTIFS[fam])
        ext_set = set()
        for s in core_set:
            ext_set.update([s - 1, s, s + 1])
        mdf["in_catalytic_motif_core"] = mdf["Site"].isin(core_set)
        mdf["in_catalytic_pocket_ext"] = mdf["Site"].isin(ext_set)
        site_dfs.append(mdf)

        rho_tbl, p_tbl = stats.spearmanr(mdf["site_effective_branch_length"], mdf["LRT"])
        rho_ent, p_ent = stats.spearmanr(mdf["entropy_aa"], mdf["LRT"])

        q25_ent = float(mdf["entropy_aa"].quantile(0.25))
        q75_ent = float(mdf["entropy_aa"].quantile(0.75))
        low_ent = mdf[(mdf["entropy_aa"] <= q25_ent) & (mdf["occupancy"] >= 0.90)]
        if len(low_ent) < 5:
            low_ent = mdf[mdf["entropy_aa"] <= q25_ent]
        high_ent = mdf[mdf["entropy_aa"] >= q75_ent]

        n_p10 = int((mdf["p-value"] < 0.10).sum())
        n_p05 = int((mdf["p-value"] < 0.05).sum())
        n_p01 = int((mdf["p-value"] < 0.01).sum())
        n_q05 = int((mdf["q_value"] < 0.05).sum())

        diag_rows.append({
            "family": fam,
            "n_taxa": n_taxa,
            "n_branches_tested": int(mdf["n_tested_branches"].iloc[0]),
            "tree_total_branch_length": round(tree_total_bl, 2),
            "n_codons": n_codons,
            "median_site_effective_branch_length": round(float(mdf["site_effective_branch_length"].median()), 2),
            "mean_pairwise_aa_pdist": round(float(np.mean(pdists)), 4),
            "mean_col_entropy_bits": round(float(np.mean(entropies)), 3),
            "mean_col_occupancy": round(float(np.mean(occs)), 4),
            "SXXK_col_lock_codon_trim_pct": top_sxxk_pub,
            "SXXK_col_lock_hmmalign_pct": hmm_lock_pct,
            "phi_obs": round(obs_phi, 4) if pd.notna(obs_phi) else np.nan,
            "phi_null_mean": round(null_phi, 4) if pd.notna(null_phi) else np.nan,
            "phi_p_value": round(p_phi, 4) if pd.notna(p_phi) else np.nan,
            "n_sig_p10": n_p10,
            "pct_sig_p10": round(100 * n_p10 / n_codons, 2),
            "n_sig_p05": n_p05,
            "pct_sig_p05": round(100 * n_p05 / n_codons, 2),
            "n_sig_p01": n_p01,
            "pct_sig_p01": round(100 * n_p01 / n_codons, 2),
            "n_sig_FDR_q05": n_q05,
            "pct_sig_FDR_q05": round(100 * n_q05 / n_codons, 2),
            "spearman_LRT_vs_EffBranchLen": round(float(rho_tbl), 3),
            "spearman_LRT_vs_EffBranchLen_p": float(p_tbl),
            "spearman_LRT_vs_Entropy": round(float(rho_ent), 3),
            "spearman_LRT_vs_Entropy_p": float(p_ent),
            "pct_p05_low_entropy_q25_cols": round(100 * float((low_ent["p-value"] < 0.05).mean()), 2) if len(low_ent) else np.nan,
            "pct_p05_high_entropy_q75_cols": round(100 * float((high_ent["p-value"] < 0.05).mean()), 2) if len(high_ent) else np.nan,
        })

    df_all = pd.concat(site_dfs, ignore_index=True)
    df_dom = pd.concat(dom_dfs, ignore_index=True)
    df_all.to_csv(f"{STAGE}/PerSite_Selection_Metrics.tsv", sep="\t", index=False)

    # Depletion (codon_trim/MEME) + domain-anchored purifying selection tests
    for region_lbl, col_mask in [
        ("Core_Catalytic_Motifs", "in_catalytic_motif_core"),
        ("Extended_Catalytic_Pocket_pm1", "in_catalytic_pocket_ext"),
    ]:
        for fam_lbl in FAMILIES + ["All_Pooled"]:
            sub = df_all if fam_lbl == "All_Pooled" else df_all[df_all["family"] == fam_lbl]
            dsub = df_dom if fam_lbl == "All_Pooled" else df_dom[df_dom["family"] == fam_lbl]

            in_m = sub[sub[col_mask]]
            out_m = sub[~sub[col_mask]]

            pos05_in = int((in_m["p-value"] < 0.05).sum())
            pos05_out = int((out_m["p-value"] < 0.05).sum())
            pos10_in = int((in_m["p-value"] < 0.10).sum())
            pos10_out = int((out_m["p-value"] < 0.10).sum())

            or05, p_fish05 = stats.fisher_exact([[pos05_in, len(in_m) - pos05_in], [pos05_out, len(out_m) - pos05_out]], alternative="less")
            or10, p_fish10 = stats.fisher_exact([[pos10_in, len(in_m) - pos10_in], [pos10_out, len(out_m) - pos10_out]], alternative="less")

            bg_rate05 = float((sub["p-value"] < 0.05).mean())
            bg_rate10 = float((sub["p-value"] < 0.10).mean())
            p_binom05 = float(stats.binomtest(pos05_in, len(in_m), bg_rate05, alternative="less").pvalue)
            p_binom10 = float(stats.binomtest(pos10_in, len(in_m), bg_rate10, alternative="less").pvalue)

            rng = np.random.default_rng(42)
            arr05 = (sub["p-value"] < 0.05).values.astype(int)
            arr10 = (sub["p-value"] < 0.10).values.astype(int)
            p_perm05 = float((np.sum([rng.choice(arr05, len(in_m), replace=False).sum() <= pos05_in for _ in range(10000)]) + 1) / 10001)
            p_perm10 = float((np.sum([rng.choice(arr10, len(in_m), replace=False).sum() <= pos10_in for _ in range(10000)]) + 1) / 10001)

            in_wm = in_m["meme_omega_minus"].dropna()
            out_wm = out_m["meme_omega_minus"].dropna()
            _, p_p1_mw = stats.mannwhitneyu(in_m["p<sup>1</sup>"], out_m["p<sup>1</sup>"], alternative="greater")
            _, p_wm_mw = stats.mannwhitneyu(in_wm, out_wm, alternative="less")

            # Direct domain-anchored codon omega (dN/dS) & entropy at homologous HMM match columns
            din = dsub[dsub[col_mask]]
            dout = dsub[~dsub[col_mask]]
            din_w = din["domain_codon_omega"].dropna()
            dout_w = dout["domain_codon_omega"].dropna()
            _, p_dom_w = stats.mannwhitneyu(din_w, dout_w, alternative="less") if (len(din_w) and len(dout_w)) else (0, np.nan)
            _, p_dom_h = stats.mannwhitneyu(din["domain_entropy_aa"], dout["domain_entropy_aa"], alternative="less")

            proof_rows.append({
                "region_definition": region_lbl,
                "family": fam_lbl,
                "n_motif_sites_codon_trim": len(in_m),
                "n_nonmotif_sites_codon_trim": len(out_m),
                "MEME_p05_in_motif": f"{pos05_in}/{len(in_m)} ({100*pos05_in/len(in_m):.2f}%)",
                "MEME_p05_outside": f"{pos05_out}/{len(out_m)} ({100*pos05_out/len(out_m):.2f}%)",
                "Fisher_OR_p05": round(float(or05), 4),
                "Fisher_p_depletion_p05": float(p_fish05),
                "Binom_p_depletion_p05": p_binom05,
                "Perm10k_p_depletion_p05": p_perm05,
                "MEME_p10_in_motif": f"{pos10_in}/{len(in_m)} ({100*pos10_in/len(in_m):.2f}%)",
                "MEME_p10_outside": f"{pos10_out}/{len(out_m)} ({100*pos10_out/len(out_m):.2f}%)",
                "Fisher_p_depletion_p10": float(p_fish10),
                "Perm10k_p_depletion_p10": p_perm10,
                "n_domain_motif_sites": len(din_w),
                "median_domain_omega_motif": round(float(din_w.median()), 4),
                "iqr_domain_omega_motif": f"{din_w.quantile(0.25):.4f}-{din_w.quantile(0.75):.4f}",
                "median_domain_omega_outside": round(float(dout_w.median()), 4),
                "domain_omega_MannWhitney_p": float(p_dom_w),
                "median_domain_Haa_motif": round(float(din["domain_entropy_aa"].median()), 3),
                "median_domain_Haa_outside": round(float(dout["domain_entropy_aa"].median()), 3),
                "domain_Haa_MannWhitney_p": float(p_dom_h),
                "median_MEME_p1_motif": round(float(in_m["p<sup>1</sup>"].median()), 4),
                "median_MEME_p1_outside": round(float(out_m["p<sup>1</sup>"].median()), 4),
                "MEME_p1_MannWhitney_p": float(p_p1_mw),
                "median_MEME_omega_minus_motif": round(float(in_wm.median()), 4),
                "median_MEME_omega_minus_outside": round(float(out_wm.median()), 4),
                "MEME_omega_minus_MannWhitney_p": float(p_wm_mw),
            })

    df_diag = pd.DataFrame(diag_rows)
    df_proof = pd.DataFrame(proof_rows)
    df_diag.to_csv(f"{STAGE}/Table_S09_MEME_Divergence_Alignment_Recombination_Diagnostics.csv", index=False)
    df_proof.to_csv(f"{STAGE}/Table_S10_Catalytic_Motif_Purifying_Selection_Proof.csv", index=False)

    print("=== MEME Divergence, Alignment & Recombination Diagnostics ===")
    print(df_diag.to_string(index=False))
    print("\n=== Catalytic Motif Depletion & Purifying Selection Tests ===")
    print(df_proof[df_proof["region_definition"] == "Core_Catalytic_Motifs"].to_string(index=False))

    # Plot Figure S03: MEME diagnostics & purifying selection tests
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), dpi=300)

    # Panel A: Threshold & Alignment Stratification of MEME Positive Selection Rates
    ax = axes[0]
    x = np.arange(3)
    w = 0.20
    ax.bar(x - 1.5*w, df_diag["pct_sig_p10"], w, label="p < 0.10 (Fig S12)", color="#E64A19", alpha=0.85)
    ax.bar(x - 0.5*w, df_diag["pct_sig_p05"], w, label="p < 0.05", color="#FFB74D", alpha=0.9)
    ax.bar(x + 0.5*w, df_diag["pct_sig_FDR_q05"], w, label="BH-FDR q < 0.05", color="#386CB0", alpha=0.9)
    ax.bar(x + 1.5*w, df_diag["pct_p05_low_entropy_q25_cols"], w, label="Conserved Q25 Cols (p<0.05)", color="#1B9E77", alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(FAMILIES, fontweight="bold")
    ax.set_ylabel("Positively Selected Sites (%)", fontweight="bold")
    ax.set_title("(A) MEME Sensitivity to Threshold & Divergence", fontweight="bold", loc="left")
    ax.legend(fontsize=7.5)

    # Panel B: Depletion of Positive Selection Inside Core Motifs vs Outside
    ax = axes[1]
    core_p = df_proof[df_proof["region_definition"] == "Core_Catalytic_Motifs"].set_index("family").loc[FAMILIES + ["All_Pooled"]]
    in_r = [float(s.split("(")[1].rstrip("%)")) for s in core_p["MEME_p10_in_motif"]]
    out_r = [float(s.split("(")[1].rstrip("%)")) for s in core_p["MEME_p10_outside"]]
    x2 = np.arange(4)
    w2 = 0.35
    ax.bar(x2 - w2/2, in_r, w2, label="Inside Core Catalytic Motifs", color="#1B9E77")
    ax.bar(x2 + w2/2, out_r, w2, label="Outside Catalytic Motifs", color="#D95F02")
    ax.set_xticks(x2)
    ax.set_xticklabels(["PF00144", "PF00905", "PF13354", "Pooled"], fontweight="bold")
    ax.set_ylabel("MEME Positive Sites (p < 0.10, %)", fontweight="bold")
    ax.set_title("(B) Positive Selection Depletion in Motifs", fontweight="bold", loc="left")
    ax.legend(fontsize=8)

    # Panel C: Domain-Anchored Codon dN/dS Inside Motifs vs Outside
    ax = axes[2]
    b_data = []
    b_lbls = []
    b_cols = []
    for fam in FAMILIES:
        sub = df_dom[df_dom["family"] == fam]
        b_data.extend([
            np.clip(sub[sub["in_catalytic_motif_core"]]["domain_codon_omega"].dropna().values, 0, 2.5),
            np.clip(sub[~sub["in_catalytic_motif_core"]]["domain_codon_omega"].dropna().values, 0, 2.5),
        ])
        b_lbls.extend([f"{fam}\nMotifs", f"{fam}\nOutside"])
        b_cols.extend(["#1B9E77", "#B0BEC5"])
    bp = ax.boxplot(b_data, tick_labels=b_lbls, patch_artist=True, showfliers=False)
    for pt, col in zip(bp["boxes"], b_cols):
        pt.set_facecolor(col)
        pt.set_alpha(0.85)
    ax.axhline(1.0, color="#D32F2F", linestyle="--", linewidth=0.9, label="Neutrality (omega = 1.0)")
    ax.set_ylabel("Domain-Anchored Codon omega (dN/dS)", fontweight="bold")
    ax.set_title("(C) Purifying Selection at Motifs (hmmalign)", fontweight="bold", loc="left")
    ax.tick_params(axis="x", labelsize=8)
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{STAGE}/Fig_S03_MEME_Diagnostics_and_Motif_Purifying.pdf", bbox_inches="tight")
    plt.savefig(f"{STAGE}/Fig_S03_MEME_Diagnostics_and_Motif_Purifying.png", bbox_inches="tight")
    plt.close()

main()
'''

def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(SUPP_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(WSL_STAGE_WIN, exist_ok=True)

    for fam in ["PF00144", "PF00905", "PF13354"]:
        fl = fam.lower()
        shutil.copy2(DATA_SEL / f"{fl}.codon_trim.fna", WSL_STAGE_WIN / f"{fl}.codon_trim.fna")
        shutil.copy2(DATA_SEL / f"{fl}.nwk", WSL_STAGE_WIN / f"{fl}.nwk")
        shutil.copy2(SEL_SRC / f"{fl}_meme.json", WSL_STAGE_WIN / f"{fl}_meme.json")
        shutil.copy2(DATA_HMM / f"{fam}.hmm", WSL_STAGE_WIN / f"{fam}.hmm")
        shutil.copy2(SP_SRC / f"{fam}_rename.faa", WSL_STAGE_WIN / f"{fam}_rename.faa")

    worker_file = WSL_STAGE_WIN / "run_m04_fast.py"
    worker_file.write_text(WSL_WORKER, encoding="utf-8")

    proc = subprocess.run(linux_cmd(WSL_PY, f"{WSL_STAGE_LIN}/run_m04_fast.py"), capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(f"07_selection_meme_and_catalytic_motifs.py failed:\n{proc.stderr}")

    shutil.copy2(WSL_STAGE_WIN / "Table_S09_MEME_Divergence_Alignment_Recombination_Diagnostics.csv",
                 TABLES_DIR / "Table_S09_MEME_Divergence_Alignment_Recombination_Diagnostics.csv")
    shutil.copy2(WSL_STAGE_WIN / "Table_S10_Catalytic_Motif_Purifying_Selection_Proof.csv",
                 TABLES_DIR / "Table_S10_Catalytic_Motif_Purifying_Selection_Proof.csv")
    shutil.copy2(WSL_STAGE_WIN / "PerSite_Selection_Metrics.tsv",
                 SUPP_DIR / "PerSite_Selection_Metrics.tsv")
    shutil.copy2(WSL_STAGE_WIN / "Fig_S03_MEME_Diagnostics_and_Motif_Purifying.pdf",
                 FIG_DIR / "Fig_S03_MEME_Diagnostics_and_Motif_Purifying.pdf")
    shutil.copy2(WSL_STAGE_WIN / "Fig_S03_MEME_Diagnostics_and_Motif_Purifying.png",
                 FIG_DIR / "Fig_S03_MEME_Diagnostics_and_Motif_Purifying.png")


if __name__ == "__main__":
    main()
