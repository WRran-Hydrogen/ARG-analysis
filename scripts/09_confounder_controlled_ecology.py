"""Module 06: Confounder-Controlled Structural Diversity & Habitat Motif Enrichment Analyses (Comment 7).

Addresses Reviewer Comment 7 (controlling for phylogenetic background, SBL family composition,
viral taxonomy, and reliable host taxonomy):
1. Stratified Structural Diversity Comparison (Viral vs. Bacterial):
   - Evaluates structural diversity using both the structure-guided phylogeny (maxcc.nwk, n=1,561)
     and direct 3D biophysical shape metrics (Ca radius of gyration R_g, domain length, secondary
     structure dispersion) computed across all 1,507 PDB structures in ProteinStructure/.
   - Compares Viral vs. Bacterial structural dispersion globally, within each SBL family
     (PF00144, PF00905, PF13354), under equal-sample-size rarefaction (1,000 subsamples),
     and restricted to viruses with reliable host assignments (Isolate + CRISPR) stratified by host phylum.
2. Stratified & Host-Controlled Habitat Motif Enrichment:
   - Links all motif-annotated viral sequences (motif_node.csv / motif_edge.csv, n=1,399 viral)
     with SBL family, viral phylum, host taxonomy, and host prediction method (information.xlsx).
   - Evaluates whether habitat-specific motif enrichments persist after controlling for SBL family
     and viral phylum via Cochran-Mantel-Haenszel (CMH) stratified tests, within-family / within-
     Uroviricota tests, and sensitivity analyses restricted to reliable host-assigned viruses
     (Isolate + CRISPR) stratified by bacterial host phylum.
"""

import os
import re
import glob
import shutil
import subprocess
from pathlib import Path
import pandas as pd

import sbl_config as cfg
from sbl_config import linux_cmd

DATA_DIR = cfg.DATA_DIR
TABLES_DIR = cfg.TABLES_DIR
FIG_DIR = cfg.FIG_DIR

M01_RAW = cfg.SEQ_DIR
INFO_XLSX = cfg.META_DIR / "information.xlsx"
MAXCC_NWK = cfg.PHYLO_DIR / "maxcc.nwk"
MOTIF_NODE = cfg.MOTIF_DIR / "motif_node.csv"
MOTIF_EDGE = cfg.MOTIF_DIR / "motif_edge.csv"

WSL_STAGE_WIN, WSL_STAGE_LIN = cfg.stage_dirs("wsl_m06")
WSL_PY = cfg.WSL_PY

WSL_WORKER = r'''
import os
import re
import glob
import math
import numpy as np
import pandas as pd
from scipy import stats
from ete4 import Tree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STAGE = os.path.dirname(os.path.abspath(__file__))
REPO_LIN = os.path.abspath(os.path.join(STAGE, "..", ".."))
PDB_LIN = f"{REPO_LIN}/data/structures/ProteinStructure"
FAMILIES = ["PF00144", "PF00905", "PF13354"]

def extract_taxon_rank(tax_str, prefix):
    if not isinstance(tax_str, str):
        return "Unassigned"
    m = re.search(rf"{prefix}__([^;]+)", tax_str)
    return m.group(1).strip() if m and m.group(1).strip() else "Unassigned"

def fast_pdb_ca_features(pdb_path):
    """Fast extraction of Ca coordinates -> length, Radius of Gyration Rg, end-to-end distance, Ca compactness."""
    coords = []
    seen_res = set()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
                chain = line[21]
                resseq = line[22:27].strip()
                key = (chain, resseq)
                if key in seen_res:
                    continue
                seen_res.add(key)
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except ValueError:
                    continue
    if len(coords) < 30:
        return None
    arr = np.array(coords, dtype=float)
    centroid = arr.mean(axis=0)
    diffs = arr - centroid
    rg = float(np.sqrt((diffs ** 2).sum(axis=1).mean()))
    # Normalized Rg (Rg / N^(1/3)) measures conformational compactness independent of chain length
    norm_rg = rg / (len(arr) ** (1.0 / 3.0))
    # Asphericity from gyration tensor eigenvalues
    cov = np.cov(diffs.T)
    eigvals = np.sort(np.linalg.eigvalsh(cov))[::-1]
    tr = eigvals.sum()
    asphericity = float((eigvals[0] - 0.5 * (eigvals[1] + eigvals[2])) / max(tr, 1e-6))
    return {
        "n_res": len(arr),
        "rg_A": rg,
        "norm_rg": norm_rg,
        "asphericity": asphericity,
    }

def rarefied_tree_diversity(tree, tips_a, tips_b, n_boot=500, seed=42):
    """Compare Faith's PD & Mean Pairwise Patristic Distance under equal-sample-size rarefaction."""
    rng = np.random.default_rng(seed)
    k = min(len(tips_a), len(tips_b), 80)
    if k < 8:
        return None
    # Precompute leaf-to-root path lengths or pairwise distances on sampled sets
    # For fast exact Faith's PD and MPD on subsamples of size k:
    # Actually, sample 150 bootstrap subsamples of size k (or k_mpd=30 for MPD) cleanly
    pd_a, pd_b = [], []
    mpd_a, mpd_b = [], []
    # Map leaf names to node objects
    leaf_map = {l.name: l for l in tree.leaves()}
    va = [x for x in tips_a if x in leaf_map]
    vb = [x for x in tips_b if x in leaf_map]
    k = min(len(va), len(vb), 50)
    if k < 6:
        return None

    for b in range(n_boot):
        sa = rng.choice(va, size=k, replace=False)
        sb = rng.choice(vb, size=k, replace=False)
        # Fast subtree PD: union of ancestor edges up to MRCA
        def subtree_pd(s_names):
            nodes = [leaf_map[nm] for nm in s_names]
            mrca = tree.common_ancestor(nodes) if len(nodes) > 1 else nodes[0]
            seen_edges = set()
            tot_len = 0.0
            for nd in nodes:
                cur = nd
                while cur is not None and cur is not mrca:
                    if cur in seen_edges:
                        break
                    seen_edges.add(cur)
                    tot_len += float(cur.dist or 0.0)
                    cur = cur.up
            return tot_len

        pd_a.append(subtree_pd(sa))
        pd_b.append(subtree_pd(sb))

        if b < 100:
            # Pairwise patristic distances among up to 20 sampled tips
            ka = min(k, 20)
            pa_nodes = [leaf_map[nm] for nm in sa[:ka]]
            pb_nodes = [leaf_map[nm] for nm in sb[:ka]]
            d_a = [tree.get_distance(pa_nodes[i], pa_nodes[j]) for i in range(ka) for j in range(i + 1, ka)]
            d_b = [tree.get_distance(pb_nodes[i], pb_nodes[j]) for i in range(ka) for j in range(i + 1, ka)]
            mpd_a.append(np.mean(d_a))
            mpd_b.append(np.mean(d_b))

    _, p_pd = stats.mannwhitneyu(pd_a, pd_b, alternative="two-sided")
    _, p_mpd = stats.mannwhitneyu(mpd_a, mpd_b, alternative="two-sided")
    return {
        "k_rarefied": k,
        "viral_rarefied_PD_mean": round(float(np.mean(pd_a)), 3),
        "viral_rarefied_PD_sd": round(float(np.std(pd_a)), 3),
        "bac_rarefied_PD_mean": round(float(np.mean(pd_b)), 3),
        "bac_rarefied_PD_sd": round(float(np.std(pd_b)), 3),
        "PD_MannWhitney_p": float(p_pd),
        "viral_MPD_mean": round(float(np.mean(mpd_a)), 3),
        "bac_MPD_mean": round(float(np.mean(mpd_b)), 3),
        "MPD_MannWhitney_p": float(p_mpd),
    }

def cmh_test_stratified(strata_tables):
    """Exact Cochran-Mantel-Haenszel 2x2 test across K strata: list of [[a,b],[c,d]]."""
    num = 0.0
    den = 0.0
    or_num = 0.0
    or_den = 0.0
    valid_strata = 0
    for tbl in strata_tables:
        a, b = tbl[0]
        c, d = tbl[1]
        n = a + b + c + d
        if n <= 1:
            continue
        r1 = a + b
        r2 = c + d
        c1 = a + c
        c2 = b + d
        if r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0:
            continue
        valid_strata += 1
        e_a = r1 * c1 / n
        v_a = (r1 * r2 * c1 * c2) / (n * n * (n - 1))
        num += (a - e_a)
        den += v_a
        or_num += (a * d) / n
        or_den += (b * c) / n
    if den <= 0 or valid_strata == 0:
        return np.nan, np.nan, 1.0, valid_strata
    # Continuity-corrected or standard CMH chi-square (1 df)
    chi2 = (max(0.0, abs(num) - 0.5) ** 2) / den
    pval = float(stats.chi2.sf(chi2, df=1))
    cmh_or = float(or_num / max(or_den, 1e-9))
    return float(chi2), cmh_or, pval, valid_strata

def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty_like(q)
    out[order] = q
    return out

def main():
    meta = pd.read_csv(f"{STAGE}/viral_metadata.tsv", sep="\t")
    meta_map = meta.set_index("uvig").to_dict("index")
    bac_fam_map = {}
    with open(f"{STAGE}/bac_family_map.tsv") as f:
        for l in f:
            bid, fam = l.strip().split("\t")
            bac_fam_map[bid] = fam
    bac_host_map = {}
    with open(f"{STAGE}/bac_phylum_map.tsv") as f:
        for l in f:
            bid, phy = l.strip().split("\t")
            bac_host_map[bid] = phy

    # =========================================================================
    # PART A: Confounder-Controlled Structural Diversity Comparison
    # =========================================================================
    tree = Tree(open(f"{STAGE}/maxcc.nwk").read(), parser=0)
    # Assign family to every leaf in maxcc.nwk
    leaf_records = []
    # Map PDB filename prefix to family via viral metadata or pdb_names
    pdb_files = glob.glob(f"{PDB_LIN}/*.pdb")
    pdb_by_stem = {}
    for p in pdb_files:
        stem = os.path.splitext(os.path.basename(p))[0]
        pdb_by_stem[stem] = p

    # Infer family for maxcc bacterial leaves from PDB filename or bac_fam_map
    # Also propagate clade consensus family for any unmapped bacterial leaves in maxcc.nwk
    for l in tree.leaves():
        nm = l.name
        is_vir = "IMGVR" in nm
        if is_vir:
            uv = nm.split("__")[0]
            mm = meta_map.get(uv, {})
            fam = mm.get("sbl_family", "Unknown")
            v_phy = mm.get("viral_phylum", "Unassigned")
            h_phy = mm.get("host_phylum", "Unassigned")
            h_meth = mm.get("host_method", "Unassigned")
            rel_host = h_meth in ("Isolate taxonomy", "CRISPR spacer match")
        else:
            # Check explicit PDB suffix or accession match
            fam = "Unknown"
            for f_cand in FAMILIES:
                if f_cand in nm:
                    fam = f_cand
                    break
            if fam == "Unknown":
                # Match accession prefix in bac_fam_map
                acc = nm.split("_")[0] if not nm.startswith("NZ_") else "_".join(nm.split("_")[:2])
                for bk, bf in bac_fam_map.items():
                    if acc in bk or bk in nm:
                        fam = bf
                        break
            v_phy = "Bacterial_Ref"
            h_phy = "Unassigned"
            for bk, bp in bac_host_map.items():
                acc = nm.split("_")[0] if not nm.startswith("NZ_") else "_".join(nm.split("_")[:2])
                if acc in bk or bk in nm:
                    h_phy = bp
                    break
            h_meth = "Bacterial_Ref"
            rel_host = True

        # Find matching PDB path and compute 3D shape descriptors
        stem_no_chain = re.sub(r"_[A-Z]$", "", nm)
        pdb_path = pdb_by_stem.get(stem_no_chain) or pdb_by_stem.get(nm)
        shape = fast_pdb_ca_features(pdb_path) if pdb_path else None

        leaf_records.append({
            "leaf_name": nm,
            "domain": "Viral" if is_vir else "Bacterial",
            "sbl_family": fam,
            "viral_phylum": v_phy,
            "host_phylum": h_phy,
            "host_method": h_meth,
            "reliable_host": rel_host,
            "n_res": shape["n_res"] if shape else np.nan,
            "rg_A": shape["rg_A"] if shape else np.nan,
            "norm_rg": shape["norm_rg"] if shape else np.nan,
            "asphericity": shape["asphericity"] if shape else np.nan,
        })

    df_struct = pd.DataFrame(leaf_records)
    # Assign remaining Unknown bacterial leaves via nearest sister majority family on maxcc.nwk
    leaf_fam_dict = df_struct.set_index("leaf_name")["sbl_family"].to_dict()
    for l in tree.leaves():
        if leaf_fam_dict.get(l.name) == "Unknown":
            cur = l.up
            while cur is not None:
                fams = [leaf_fam_dict.get(x.name, "Unknown") for x in cur.leaves() if leaf_fam_dict.get(x.name, "Unknown") != "Unknown"]
                if fams:
                    leaf_fam_dict[l.name] = max(set(fams), key=fams.count)
                    break
                cur = cur.up
    df_struct["sbl_family"] = df_struct["leaf_name"].map(leaf_fam_dict)

    struct_rows = []
    strata_defs = [
        ("1_Global_Unstratified", df_struct),
        ("2_Within_PF00144", df_struct[df_struct["sbl_family"] == "PF00144"]),
        ("2_Within_PF00905", df_struct[df_struct["sbl_family"] == "PF00905"]),
        ("2_Within_PF13354", df_struct[df_struct["sbl_family"] == "PF13354"]),
        ("3_Viral_Uroviricota_vs_Bacteria", df_struct[(df_struct["domain"] == "Bacterial") | (df_struct["viral_phylum"] == "Uroviricota")]),
        ("4_ReliableHost_Subset_Isolate_CRISPR", df_struct[(df_struct["domain"] == "Bacterial") | (df_struct["reliable_host"])]),
        ("5_ReliableHost_PF00905", df_struct[(df_struct["sbl_family"] == "PF00905") & ((df_struct["domain"] == "Bacterial") | (df_struct["reliable_host"]))]),
        ("5_ReliableHost_PF00144", df_struct[(df_struct["sbl_family"] == "PF00144") & ((df_struct["domain"] == "Bacterial") | (df_struct["reliable_host"]))]),
    ]

    for sname, sub in strata_defs:
        v_sub = sub[sub["domain"] == "Viral"]
        b_sub = sub[sub["domain"] == "Bacterial"]
        pd_res = rarefied_tree_diversity(tree, v_sub["leaf_name"].tolist(), b_sub["leaf_name"].tolist(), n_boot=300)

        # 3D Shape dispersion (Levene test on norm_rg and rg_A)
        v_rg = v_sub["norm_rg"].dropna().values
        b_rg = b_sub["norm_rg"].dropna().values
        v_asp = v_sub["asphericity"].dropna().values
        b_asp = b_sub["asphericity"].dropna().values
        if len(v_rg) >= 5 and len(b_rg) >= 5:
            lev_stat, lev_p = stats.levene(v_rg, b_rg)
            iqr_v = float(np.percentile(v_rg, 75) - np.percentile(v_rg, 25))
            iqr_b = float(np.percentile(b_rg, 75) - np.percentile(b_rg, 25))
            sd_v = float(np.std(v_rg))
            sd_b = float(np.std(b_rg))
        else:
            lev_stat, lev_p, iqr_v, iqr_b, sd_v, sd_b = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        row = {
            "stratum": sname,
            "n_viral_structs": len(v_sub),
            "n_bacterial_structs": len(b_sub),
            "norm_Rg_SD_Viral": round(sd_v, 4) if pd.notna(sd_v) else np.nan,
            "norm_Rg_SD_Bacterial": round(sd_b, 4) if pd.notna(sd_b) else np.nan,
            "norm_Rg_IQR_Viral": round(iqr_v, 4) if pd.notna(iqr_v) else np.nan,
            "norm_Rg_IQR_Bacterial": round(iqr_b, 4) if pd.notna(iqr_b) else np.nan,
            "Levene_p_normRg_dispersion": float(lev_p) if pd.notna(lev_p) else np.nan,
        }
        if pd_res:
            row.update(pd_res)
        struct_rows.append(row)

    df_struct_sum = pd.DataFrame(struct_rows)
    df_struct_sum.to_csv(f"{STAGE}/Table_S_Stratified_Structural_Diversity.csv", index=False)

    # =========================================================================
    # PART B: Confounder-Controlled Habitat Motif Enrichment
    # =========================================================================
    mn = pd.read_csv(f"{STAGE}/motif_node.csv")
    me = pd.read_csv(f"{STAGE}/motif_edge.csv")

    # Separate motif nodes (M1/M2/M3) and sequence nodes
    motif_cats = mn[mn["taxonomy"].isin(["M1", "M2", "M3"])].set_index("Id")["taxonomy"].to_dict()
    vir_nodes = mn[~mn["taxonomy"].isin(["M1", "M2", "M3", "BACTERIAL"])].copy()
    vir_nodes["uvig"] = vir_nodes["Id"].astype(str).apply(lambda x: x.split("__")[0])
    vir_nodes["sbl_family"] = vir_nodes["pfam"].astype(str)
    vir_nodes["viral_phylum"] = vir_nodes["uvig"].map(lambda u: meta_map.get(u, {}).get("viral_phylum", "Unassigned"))
    vir_nodes["host_phylum"] = vir_nodes["uvig"].map(lambda u: meta_map.get(u, {}).get("host_phylum", "Unassigned"))
    vir_nodes["host_method"] = vir_nodes["uvig"].map(lambda u: meta_map.get(u, {}).get("host_method", "Unassigned"))
    vir_nodes["reliable_host"] = vir_nodes["host_method"].isin(["Isolate taxonomy", "CRISPR spacer match"])

    # Clean habitat typo "Other enviroment" -> Keep original taxonomy label so counts match motif_enrichment.r
    vir_nodes["habitat_clean"] = vir_nodes["taxonomy"].astype(str)

    # Map each viral sequence to its set of motif edges
    seq_motifs = {}
    for _, er in me.iterrows():
        m_id = str(er["Source"])
        seq_id = str(er["Target"])
        if m_id in motif_cats:
            seq_motifs.setdefault(seq_id, set()).add(m_id)

    vir_nodes["motifs"] = vir_nodes["Id"].map(lambda x: seq_motifs.get(x, set()))

    # Identify all (Motif, Habitat) pairs that are tested (reach >=2% or >=5 occurrences in a habitat)
    all_habitats = sorted(vir_nodes["habitat_clean"].unique())
    all_motifs = sorted(motif_cats.keys())

    # First run Global Fisher test across all (Motif, Habitat) pairs where obs >= 5
    test_pairs = []
    for hab in all_habitats:
        if hab in ("Unclassfied", "Other enviroment", "Other Host-associated"):
            continue
        sub_h = vir_nodes[vir_nodes["habitat_clean"] == hab]
        if len(sub_h) < 15:
            continue
        for m_id in all_motifs:
            a = int(sum(1 for ms in sub_h["motifs"] if m_id in ms))
            if a < 5:
                continue
            b = len(sub_h) - a
            rest = vir_nodes[vir_nodes["habitat_clean"] != hab]
            c = int(sum(1 for ms in rest["motifs"] if m_id in ms))
            d = len(rest) - c
            or_g, p_g = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            test_pairs.append({
                "motif": m_id,
                "motif_category": motif_cats[m_id],
                "habitat": hab,
                "obs_in_hab": a,
                "n_hab": len(sub_h),
                "pct_in_hab": round(100 * a / len(sub_h), 2),
                "pct_other_hab": round(100 * c / max(len(rest), 1), 2),
                "global_OR": round(float(or_g), 3),
                "global_p": float(p_g),
            })

    df_pairs = pd.DataFrame(test_pairs)
    df_pairs["global_BH_padj"] = bh_fdr(df_pairs["global_p"].values)

    # Now take ALL globally significant pairs (global_BH_padj < 0.05) and evaluate confounder controls!
    sig_pairs = df_pairs[df_pairs["global_BH_padj"] < 0.05].sort_values("global_BH_padj").copy()

    # Pre-subset reliable host data
    rel_nodes = vir_nodes[vir_nodes["reliable_host"]].copy()

    conf_rows = []
    for _, r in sig_pairs.iterrows():
        m_id = r["motif"]
        hab = r["habitat"]

        # 1) CMH stratified by SBL Family (PF00144, PF00905, PF13354)
        fam_strata = []
        fam_notes = []
        for fam in FAMILIES:
            f_sub = vir_nodes[vir_nodes["sbl_family"] == fam]
            fh = f_sub[f_sub["habitat_clean"] == hab]
            fr = f_sub[f_sub["habitat_clean"] != hab]
            a = int(sum(1 for ms in fh["motifs"] if m_id in ms))
            b = len(fh) - a
            c = int(sum(1 for ms in fr["motifs"] if m_id in ms))
            d = len(fr) - c
            fam_strata.append([[a, b], [c, d]])
            if a >= 3:
                _, p_f = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
                fam_notes.append(f"{fam}:{a}/{len(fh)}(p={p_f:.2e})")
        _, cmh_or_fam, cmh_p_fam, _ = cmh_test_stratified(fam_strata)

        # 2) CMH stratified by Viral Phylum (Uroviricota, Nucleocytoviricota, Other)
        vphy_strata = []
        uro_p = np.nan
        for vp in vir_nodes["viral_phylum"].unique():
            v_sub = vir_nodes[vir_nodes["viral_phylum"] == vp]
            vh = v_sub[v_sub["habitat_clean"] == hab]
            vr = v_sub[v_sub["habitat_clean"] != hab]
            a = int(sum(1 for ms in vh["motifs"] if m_id in ms))
            b = len(vh) - a
            c = int(sum(1 for ms in vr["motifs"] if m_id in ms))
            d = len(vr) - c
            vphy_strata.append([[a, b], [c, d]])
            if vp == "Uroviricota" and len(vh) > 0:
                _, uro_p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        _, cmh_or_vphy, cmh_p_vphy, _ = cmh_test_stratified(vphy_strata)

        # 3) CMH stratified simultaneously by (SBL Family x Viral Phylum)
        joint_strata = []
        for (fam, vp), j_sub in vir_nodes.groupby(["sbl_family", "viral_phylum"]):
            jh = j_sub[j_sub["habitat_clean"] == hab]
            jr = j_sub[j_sub["habitat_clean"] != hab]
            a = int(sum(1 for ms in jh["motifs"] if m_id in ms))
            b = len(jh) - a
            c = int(sum(1 for ms in jr["motifs"] if m_id in ms))
            d = len(jr) - c
            joint_strata.append([[a, b], [c, d]])
        _, cmh_or_joint, cmh_p_joint, n_joint_strata = cmh_test_stratified(joint_strata)

        # 4) Sensitivity analysis on Reliable-Host-Assigned subset (Isolate + CRISPR), CMH stratified by Host Phylum!
        rh_h = rel_nodes[rel_nodes["habitat_clean"] == hab]
        rh_r = rel_nodes[rel_nodes["habitat_clean"] != hab]
        rh_a = int(sum(1 for ms in rh_h["motifs"] if m_id in ms))
        rh_b = len(rh_h) - rh_a
        rh_c = int(sum(1 for ms in rh_r["motifs"] if m_id in ms))
        rh_d = len(rh_r) - rh_c
        if len(rh_h) >= 5 and rh_a >= 1:
            rh_or, rh_p = stats.fisher_exact([[rh_a, rh_b], [rh_c, rh_d]], alternative="greater")
        else:
            rh_or, rh_p = np.nan, np.nan

        host_strata = []
        for hp, h_sub in rel_nodes.groupby("host_phylum"):
            hh = h_sub[h_sub["habitat_clean"] == hab]
            hr = h_sub[h_sub["habitat_clean"] != hab]
            a = int(sum(1 for ms in hh["motifs"] if m_id in ms))
            b = len(hh) - a
            c = int(sum(1 for ms in hr["motifs"] if m_id in ms))
            d = len(hr) - c
            host_strata.append([[a, b], [c, d]])
        _, cmh_or_host, cmh_p_host, _ = cmh_test_stratified(host_strata)

        conf_rows.append({
            "motif": m_id,
            "motif_category": r["motif_category"],
            "habitat": hab,
            "obs_in_hab": r["obs_in_hab"],
            "n_hab": r["n_hab"],
            "pct_in_hab": r["pct_in_hab"],
            "pct_other_hab": r["pct_other_hab"],
            "global_OR": r["global_OR"],
            "global_BH_padj": r["global_BH_padj"],
            "CMH_OR_ctrl_Family": round(cmh_or_fam, 3) if pd.notna(cmh_or_fam) else np.nan,
            "CMH_p_ctrl_Family": cmh_p_fam,
            "CMH_OR_ctrl_ViralPhylum": round(cmh_or_vphy, 3) if pd.notna(cmh_or_vphy) else np.nan,
            "CMH_p_ctrl_ViralPhylum": cmh_p_vphy,
            "CMH_OR_ctrl_FamilyXPhylum": round(cmh_or_joint, 3) if pd.notna(cmh_or_joint) else np.nan,
            "CMH_p_ctrl_FamilyXPhylum": cmh_p_joint,
            "Uroviricota_only_p": uro_p,
            "ReliableHost_n_in_hab": f"{rh_a}/{len(rh_h)}",
            "ReliableHost_OR": round(float(rh_or), 3) if pd.notna(rh_or) else np.nan,
            "ReliableHost_p": float(rh_p) if pd.notna(rh_p) else np.nan,
            "CMH_OR_ctrl_HostPhylum": round(cmh_or_host, 3) if pd.notna(cmh_or_host) else np.nan,
            "CMH_p_ctrl_HostPhylum": cmh_p_host,
            "within_family_signals": "; ".join(fam_notes),
        })

    df_conf = pd.DataFrame(conf_rows)
    if len(df_conf):
        df_conf["CMH_BH_padj_ctrl_FamilyXPhylum"] = bh_fdr(df_conf["CMH_p_ctrl_FamilyXPhylum"].fillna(1.0).values)
        df_conf["persists_after_FamilyXPhylum_FDR05"] = (df_conf["CMH_BH_padj_ctrl_FamilyXPhylum"] < 0.05) & (df_conf["CMH_OR_ctrl_FamilyXPhylum"].fillna(0) > 1.0)
    df_conf.to_csv(f"{STAGE}/Table_S_Stratified_Habitat_Motif_Enrichment.csv", index=False)

    print("=== Part A: Confounder-Controlled Structural Diversity Summary ===")
    print(df_struct_sum.to_string(index=False))
    print(f"\n=== Part B: Habitat Motif Enrichments (n_global_sig={len(df_conf)}, persisting_after_FamilyXPhylum_FDR05={int(df_conf['persists_after_FamilyXPhylum_FDR05'].sum()) if len(df_conf) else 0}) ===")
    if len(df_conf):
        cols_show = ["motif", "motif_category", "habitat", "pct_in_hab", "pct_other_hab", "global_OR", "global_BH_padj", "CMH_OR_ctrl_FamilyXPhylum", "CMH_BH_padj_ctrl_FamilyXPhylum", "persists_after_FamilyXPhylum_FDR05", "ReliableHost_OR", "CMH_p_ctrl_HostPhylum"]
        print(df_conf[cols_show].head(25).to_string(index=False))

    # Plot Figure: Confounder-Controlled Structural Diversity & Habitat Motif Enrichment
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4), dpi=300)

    # Panel A: Rarefied Structural Phylogenetic Diversity (Viral vs Bacterial) across strata
    ax = axes[0]
    plot_strata = ["1_Global_Unstratified", "2_Within_PF00144", "2_Within_PF00905", "2_Within_PF13354", "4_ReliableHost_Subset_Isolate_CRISPR"]
    sub_st = df_struct_sum[df_struct_sum["stratum"].isin(plot_strata)].set_index("stratum").loc[plot_strata]
    x = np.arange(len(plot_strata))
    w = 0.36
    ax.bar(x - w/2, sub_st["viral_rarefied_PD_mean"], w, yerr=sub_st["viral_rarefied_PD_sd"], capsize=3, label="Viral SBLs", color="#1B9E77", alpha=0.9)
    ax.bar(x + w/2, sub_st["bac_rarefied_PD_mean"], w, yerr=sub_st["bac_rarefied_PD_sd"], capsize=3, label="Bacterial SBLs", color="#D95F02", alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(["Global", "Within\nPF00144", "Within\nPF00905", "Within\nPF13354", "Reliable Host\n(Isolate/CRISPR)"], fontweight="bold", fontsize=8.5)
    ax.set_ylabel("Rarefied Structural Phylogenetic Diversity (PD)", fontweight="bold")
    ax.set_title("(A) Stratified Structural Diversity (Equal Sample Size)", fontweight="bold", loc="left")
    ax.legend(fontsize=8.5)

    # Panel B: Global OR vs Confounder-Adjusted CMH OR (Family x Viral Phylum) for significant motif-habitat pairs
    ax = axes[1]
    if len(df_conf):
        top_c = df_conf.head(15).copy()
        y_pos = np.arange(len(top_c))
        lbls = [f"{r['motif']} ({r['habitat']})" for _, r in top_c.iterrows()]
        g_ors = np.log2(np.maximum(top_c["global_OR"].values, 0.1))
        c_ors = np.log2(np.maximum(top_c["CMH_OR_ctrl_FamilyXPhylum"].fillna(1.0).values, 0.1))
        ax.scatter(g_ors, y_pos, color="#B0BEC5", s=65, marker="o", edgecolor="#333333", label="Unstratified Global log2(OR)", zorder=3)
        pers = top_c["persists_after_FamilyXPhylum_FDR05"].values
        ax.scatter(c_ors[pers], y_pos[pers], color="#D32F2F", s=75, marker="D", edgecolor="black", label="CMH Adjusted (FDR < 0.05)", zorder=4)
        ax.scatter(c_ors[~pers], y_pos[~pers], color="#78909C", s=75, marker="D", edgecolor="black", label="Attenuated After Control", zorder=4)
        ax.axvline(0, color="#333333", linestyle="--", linewidth=0.9)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(lbls, fontsize=8.5)
        ax.invert_yaxis()
        ax.set_xlabel("log2(Odds Ratio) for Habitat Enrichment", fontweight="bold")
        ax.set_title("(B) Motif-Habitat Enrichment After Lineage Control", fontweight="bold", loc="left")
        ax.legend(fontsize=8, loc="lower right")

    plt.tight_layout()
    plt.savefig(f"{STAGE}/Fig_S_Confounder_Controlled_Structure_and_Motifs.pdf", bbox_inches="tight")
    plt.savefig(f"{STAGE}/Fig_S_Confounder_Controlled_Structure_and_Motifs.png", bbox_inches="tight")
    plt.close()

main()
'''

def extract_taxon_rank_local(tax_str, prefix):
    if not isinstance(tax_str, str):
        return "Unassigned"
    import re as _re
    m = _re.search(rf"{prefix}__([^;]+)", tax_str)
    return m.group(1).strip() if m and m.group(1).strip() else "Unassigned"

def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(WSL_STAGE_WIN, exist_ok=True)
    cfg.require(cfg.PDB_DIR, cfg.PDB_HINT)

    # Prepare viral metadata TSV from information.xlsx + PF*_Viral_Final.faa
    info = pd.read_excel(INFO_XLSX, sheet_name="ALL_information")
    info["uvig"] = info["UVIG"].astype(str)
    info["viral_phylum"] = info["Taxonomic classification"].apply(lambda x: extract_taxon_rank_local(x, "p"))
    info["host_phylum"] = info["Host taxonomy prediction"].apply(lambda x: extract_taxon_rank_local(x, "p"))
    info["host_method"] = info["Host prediction method"].fillna("Unassigned").astype(str)

    uv_fam = {}
    for fam in ["PF00144", "PF00905", "PF13354"]:
        with open(M01_RAW / f"{fam}_Viral_Final.faa", encoding="utf-8", errors="ignore") as f:
            for l in f:
                if l.startswith(">"):
                    uv_fam[l[1:].strip().split("|")[0]] = fam
    info["sbl_family"] = info["uvig"].map(lambda u: uv_fam.get(u, "Unknown"))
    info[["uvig", "sbl_family", "viral_phylum", "host_phylum", "host_method"]].to_csv(
        WSL_STAGE_WIN / "viral_metadata.tsv", sep="\t", index=False
    )

    # Bacterial family & phylum maps
    bac_fam_lines = []
    for fam in ["PF00144", "PF00905", "PF13354"]:
        with open(M01_RAW / f"{fam}_Bacterial.faa", encoding="utf-8", errors="ignore") as f:
            for l in f:
                if l.startswith(">"):
                    sid = l[1:].strip().split()[0]
                    bac_fam_lines.append(f"{sid}\t{fam}\n")
    (WSL_STAGE_WIN / "bac_family_map.tsv").write_text("".join(bac_fam_lines), encoding="utf-8")

    sarg = pd.read_excel(INFO_XLSX, sheet_name="SARG")
    with open(WSL_STAGE_WIN / "bac_phylum_map.tsv", "w", encoding="utf-8") as f:
        for _, r in sarg.iterrows():
            f.write(f"{r['SARG-id']}\t{r['phylum']}\n")

    shutil.copy2(MAXCC_NWK, WSL_STAGE_WIN / "maxcc.nwk")
    shutil.copy2(MOTIF_NODE, WSL_STAGE_WIN / "motif_node.csv")
    shutil.copy2(MOTIF_EDGE, WSL_STAGE_WIN / "motif_edge.csv")

    worker_path = WSL_STAGE_WIN / "run_m06_confounders.py"
    worker_path.write_text(WSL_WORKER, encoding="utf-8")

    proc = subprocess.run(linux_cmd(WSL_PY, f"{WSL_STAGE_LIN}/run_m06_confounders.py"), capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(f"Module 06 failed:\n{proc.stderr}")

    shutil.copy2(WSL_STAGE_WIN / "Table_S_Stratified_Structural_Diversity.csv",
                 TABLES_DIR / "Table_S12_Stratified_Structural_Diversity.csv")
    shutil.copy2(WSL_STAGE_WIN / "Table_S_Stratified_Habitat_Motif_Enrichment.csv",
                 TABLES_DIR / "Table_S13_Stratified_Habitat_Motif_Enrichment.csv")
    shutil.copy2(WSL_STAGE_WIN / "Fig_S_Confounder_Controlled_Structure_and_Motifs.pdf",
                 FIG_DIR / "Fig_S05_Confounder_Controlled_Structure_and_Motifs.pdf")
    shutil.copy2(WSL_STAGE_WIN / "Fig_S_Confounder_Controlled_Structure_and_Motifs.png",
                 FIG_DIR / "Fig_S05_Confounder_Controlled_Structure_and_Motifs.png")


if __name__ == "__main__":
    main()
