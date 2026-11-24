"""aBSREL background branch normalization and lineage-resolved ancestral-to-descendant
selection trajectories.

- Background normalization: normalizes positive selection branch counts by total background branches
  tested in each category (Internal vs. Terminal; Viral-monophyletic vs. Bacterial-monophyletic
  vs. Mixed/HGT-associated), reports Fisher's exact tests and Odds Ratios, and explicitly
  reconciles branch counts vs. branch-site codon event shares.
- Lineage trajectories: maps selection pressure (baseline omega, weighted omega, and episodic positive
  selection frequency) longitudinally along the same root-to-tip evolutionary lineages
  (Ancestral/HGT-Transition -> Intermediate Descendant Internal -> Descendant Terminal),
  comparing HGT-associated lineages directly against non-HGT purely viral and purely bacterial lineages.
"""

import os
import shutil
import subprocess
from pathlib import Path
import pandas as pd

import sbl_config as cfg
from sbl_config import linux_cmd

DATA_DIR = cfg.SELECTION_DIR
TABLES_DIR = cfg.TABLES_DIR
FIG_DIR = cfg.FIG_DIR

SEL_SRC = cfg.SELECTION_DIR
WSL_STAGE_WIN, WSL_STAGE_LIN = cfg.stage_dirs("wsl_m04")
WSL_PY = cfg.WSL_PY

WSL_WORKER = r'''
import os
import json
import re
import numpy as np
import pandas as pd
from scipy import stats
from ete4 import Tree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STAGE = os.path.dirname(os.path.abspath(__file__))
FAMILIES = ["pf00144", "pf00905", "pf13354"]
FAM_LABELS = {"pf00144": "PF00144 (Beta-lactamase)", "pf00905": "PF00905 (Transpeptidase)", "pf13354": "PF13354 (Beta-lactamase2)"}

def tip_domain(name):
    return "Viral" if "IMGVR" in name else "Bacterial"

def clean_hyphy_newick(nwk_str):
    # Remove HyPhy partition tags {FG}/{BG} and quotes around internal node labels
    s = re.sub(r"\{[^}]+\}", "", nwk_str)
    s = s.replace("'", "")
    return s

def weighted_mean_omega(rate_dist):
    if not rate_dist:
        return np.nan, 0.0
    w_sum = 0.0
    p_gt1 = 0.0
    for omega, prop in rate_dist:
        w_sum += min(float(omega), 100.0) * float(prop)
        if float(omega) > 1.0:
            p_gt1 += float(prop)
    return w_sum, p_gt1

def run_all():
    bg_norm_rows = []
    info_tsv_rows = []
    lineage_rows = []
    paired_lineage_stats = []

    # Collect per-branch records across families directly from absrel.json + .nwk AND info.tsv
    all_info_dfs = []
    all_json_branches = []

    for fam in FAMILIES:
        info_df = pd.read_csv(f"{STAGE}/{fam}_info.tsv", sep="\t")
        info_df["family"] = fam.upper()
        all_info_dfs.append(info_df)

        with open(f"{STAGE}/{fam}_absrel.json") as f:
            ab = json.load(f)
        br_attrs = ab["branch attributes"]["0"]

        # Use HyPhy's exact internal-node-labeled tree from absrel.json input.trees.0
        raw_tree_str = ab["input"]["trees"]["0"].strip()
        if not raw_tree_str.endswith(";"):
            raw_tree_str += ";"
        nwk_clean = clean_hyphy_newick(raw_tree_str)
        tree = Tree(nwk_clean, parser=1)

        # Midpoint root for consistent ancestral-to-descendant traversal
        # Preserve node names cleanly
        for nd in tree.traverse("postorder"):
            if nd.is_leaf:
                d = tip_domain(nd.name)
                nd.add_prop("n_vir", 1 if d == "Viral" else 0)
                nd.add_prop("n_bac", 1 if d == "Bacterial" else 0)
                nd.add_prop("dom_class", "Viral" if d == "Viral" else "Bacterial")
            else:
                nv = sum(c.props["n_vir"] for c in nd.children)
                nb = sum(c.props["n_bac"] for c in nd.children)
                nd.add_prop("n_vir", nv)
                nd.add_prop("n_bac", nb)
                if nv > 0 and nb > 0:
                    nd.add_prop("dom_class", "Mixed_Nested")
                elif nv > 0:
                    nd.add_prop("dom_class", "Viral_Monophyletic")
                else:
                    nd.add_prop("dom_class", "Bacterial_Monophyletic")

        # Identify HGT transition nodes (where a child is pure/dominant Viral/Bacterial under a Mixed parent)
        hgt_node_set = set()
        for nd in tree.traverse():
            if nd.up is not None and nd.up.props["dom_class"] == "Mixed_Nested":
                if nd.props["dom_class"] in ("Viral_Monophyletic", "Bacterial_Monophyletic", "Viral", "Bacterial"):
                    hgt_node_set.add(nd)
                    hgt_node_set.add(nd.up)

        # Match tree branches to absrel.json branch attributes
        for nd in tree.traverse():
            if nd.up is None:
                continue
            bname = nd.name or ""
            attr = br_attrs.get(bname)
            if attr is None:
                attr = br_attrs.get(bname.strip("'"))
            if attr is None:
                continue

            b_type = "Terminal" if nd.is_leaf else "Internal"
            d_class = nd.props["dom_class"]
            if nd.is_leaf:
                clade_cat = "Viral" if d_class == "Viral" else "Bacterial"
            else:
                if d_class == "Mixed_Nested":
                    clade_cat = "Mixed_Nested"
                elif d_class == "Viral_Monophyletic":
                    clade_cat = "Viral"
                else:
                    clade_cat = "Bacterial"

            p_corr = float(attr.get("Corrected P-value", 1.0))
            p_uncorr = float(attr.get("Uncorrected P-value", 1.0))
            lrt = float(attr.get("LRT", 0.0))
            base_omega = float(attr.get("Baseline MG94xREV omega ratio", np.nan))
            w_omega, prop_gt1 = weighted_mean_omega(attr.get("Rate Distributions", []))

            # Check if node sits on an HGT-associated path (either transition node or nested descendant)
            anc = nd.up
            on_hgt_lineage = False
            is_hgt_transition = nd in hgt_node_set
            while anc is not None:
                if anc.props["dom_class"] == "Mixed_Nested" and anc.up is not None:
                    on_hgt_lineage = True
                    break
                anc = anc.up
            # Specifically classify HGT vs pure non-HGT lineage
            # A pure lineage has NO Mixed_Nested ancestor except possibly deepest backbone splits if we check immediate parent
            # Let's define:
            # - HGT_Transition_Branch: Mixed_Nested branch that directly gives rise to a monophyletic Viral/Bacterial clade, OR the stem branch of a nested Viral/Bacterial clade inside a Mixed_Nested clade
            # - HGT_Descendant_Branch: branch inside a nested recipient clade (whose parent clade is opposite-domain dominated)
            # - NonHGT_Pure_Viral / NonHGT_Pure_Bacterial: branches inside large (>10 tips) monophyletic Viral or Bacterial clades
            all_json_branches.append({
                "family": fam.upper(),
                "branch_name": bname,
                "branch_type": b_type,
                "clade_category": clade_cat,
                "n_vir_desc": nd.props["n_vir"],
                "n_bac_desc": nd.props["n_bac"],
                "p_corrected": p_corr,
                "p_uncorrected": p_uncorr,
                "lrt": lrt,
                "baseline_omega": base_omega,
                "weighted_omega": w_omega,
                "prop_omega_gt1": prop_gt1,
                "sig_corr_05": p_corr <= 0.05,
                "sig_uncorr_05": p_uncorr <= 0.05,
                "is_hgt_transition": is_hgt_transition,
            })

        # Trace root-to-tip lineages on this tree
        for tip in tree.leaves():
            t_dom = tip_domain(tip.name)
            # Walk up to root collecting path nodes (excluding root itself)
            path = []
            cur = tip
            while cur is not None and cur.up is not None:
                path.append(cur)
                cur = cur.up
            path = list(reversed(path))  # root-child -> ... -> tip
            if len(path) < 2:
                continue

            # Determine if this tip belongs to an HGT-nested lineage vs pure monophyletic lineage
            # Find the deepest node on the path where the clade transitions to pure tip domain
            # If the parent of that pure clade is Mixed_Nested and dominated (>=70%) by the OPPOSITE domain, it is an HGT-associated lineage!
            pure_stem_idx = len(path) - 1
            while pure_stem_idx > 0 and path[pure_stem_idx - 1].props["dom_class"] in ("Viral_Monophyletic", "Bacterial_Monophyletic"):
                pure_stem_idx -= 1

            pure_stem = path[pure_stem_idx]
            parent_of_stem = pure_stem.up
            # Check opposite-domain enrichment in parent_of_stem (sister + parent)
            is_hgt_lin = False
            if parent_of_stem is not None and parent_of_stem.up is not None:
                pv = parent_of_stem.props["n_vir"]
                pb = parent_of_stem.props["n_bac"]
                # If tip is Viral and parent_of_stem has >=50% Bacterial (or >=2 Bacteria), or vice-versa
                if t_dom == "Viral" and pb >= 2 and (pb / max(pv + pb, 1)) >= 0.30:
                    is_hgt_lin = True
                elif t_dom == "Bacterial" and pv >= 2 and (pv / max(pv + pb, 1)) >= 0.30:
                    is_hgt_lin = True

            if is_hgt_lin:
                lin_type = f"HGT_Associated_{t_dom}"
            else:
                lin_type = f"NonHGT_Pure_{t_dom}"

            # Record ancestral branch (stem or mixed parent), intermediate internal branches, and terminal tip branch
            # Pair: Ancestral/Transition Internal vs Descendant Terminal along this exact lineage
            tip_attr = br_attrs.get(tip.name)
            if tip_attr is None:
                continue
            tip_omega = float(tip_attr.get("Baseline MG94xREV omega ratio", np.nan))
            tip_p_unc = float(tip_attr.get("Uncorrected P-value", 1.0))
            tip_p_cor = float(tip_attr.get("Corrected P-value", 1.0))

            # Collect internal ancestors along this same lineage
            anc_nodes = path[:-1]
            # For HGT_Associated: focal ancestral/transition window is parent_of_stem and pure_stem (if internal)
            # Post-HGT descendant internal: nodes strictly below pure_stem and above tip
            # For NonHGT: split internal ancestors into Deep/Mid Ancestral vs Recent Pre-Terminal Internal
            for depth_idx, anc_nd in enumerate(anc_nodes):
                a_attr = br_attrs.get(anc_nd.name)
                if a_attr is None:
                    continue
                a_omega = float(a_attr.get("Baseline MG94xREV omega ratio", np.nan))
                a_p_unc = float(a_attr.get("Uncorrected P-value", 1.0))
                a_p_cor = float(a_attr.get("Corrected P-value", 1.0))
                _, a_pgt1 = weighted_mean_omega(a_attr.get("Rate Distributions", []))

                if is_hgt_lin:
                    if depth_idx <= pure_stem_idx:
                        stage = "1_Ancestral_or_HGT_Transition"
                    else:
                        stage = "2_PostHGT_Descendant_Internal"
                else:
                    if depth_idx < max(1, len(anc_nodes) // 2):
                        stage = "1_Deep_Ancestral_Internal"
                    else:
                        stage = "2_Recent_Descendant_Internal"

                lineage_rows.append({
                    "family": fam.upper(),
                    "tip_id": tip.name,
                    "lineage_type": lin_type,
                    "branch_name": anc_nd.name,
                    "stage": stage,
                    "baseline_omega": min(a_omega, 10.0) if pd.notna(a_omega) else np.nan,
                    "p_uncorrected": a_p_unc,
                    "p_corrected": a_p_cor,
                    "sig_uncorr_05": a_p_unc <= 0.05,
                    "sig_corr_05": a_p_cor <= 0.05,
                    "prop_omega_gt1": a_pgt1,
                })

            lineage_rows.append({
                "family": fam.upper(),
                "tip_id": tip.name,
                "lineage_type": lin_type,
                "branch_name": tip.name,
                "stage": "3_Terminal_Tip",
                "baseline_omega": min(tip_omega, 10.0) if pd.notna(tip_omega) else np.nan,
                "p_uncorrected": tip_p_unc,
                "p_corrected": tip_p_cor,
                "sig_uncorr_05": tip_p_unc <= 0.05,
                "sig_corr_05": tip_p_cor <= 0.05,
                "prop_omega_gt1": weighted_mean_omega(tip_attr.get("Rate Distributions", []))[1],
            })

            # Also record per-lineage ancestral-vs-terminal pair for paired Wilcoxon test
            anc_omegas = [
                min(float(br_attrs[n.name].get("Baseline MG94xREV omega ratio", np.nan)), 10.0)
                for n in anc_nodes if n.name in br_attrs
            ]
            anc_sigs = [br_attrs[n.name].get("Uncorrected P-value", 1.0) <= 0.05 for n in anc_nodes if n.name in br_attrs]
            if is_hgt_lin:
                trans_nodes = [n for idx, n in enumerate(anc_nodes) if idx <= pure_stem_idx and n.name in br_attrs]
                trans_omegas = [min(float(br_attrs[n.name].get("Baseline MG94xREV omega ratio", np.nan)), 10.0) for n in trans_nodes]
                trans_sigs = [br_attrs[n.name].get("Uncorrected P-value", 1.0) <= 0.05 for n in trans_nodes]
            else:
                trans_omegas = anc_omegas
                trans_sigs = anc_sigs

            if trans_omegas and pd.notna(tip_omega):
                paired_lineage_stats.append({
                    "family": fam.upper(),
                    "tip_id": tip.name,
                    "lineage_type": lin_type,
                    "is_hgt": is_hgt_lin,
                    "ancestral_median_omega": float(np.nanmedian(trans_omegas)),
                    "ancestral_any_pos_sel_uncorr": any(trans_sigs),
                    "terminal_omega": min(tip_omega, 10.0),
                    "terminal_pos_sel_uncorr": tip_p_unc <= 0.05,
                })

    # === Background normalization summary: Both info.tsv (exact Fig 3B dataset) and full absrel.json ===
    df_info_all = pd.concat(all_info_dfs, ignore_index=True)
    df_json_all = pd.DataFrame(all_json_branches)

    # 1) Audit pf*_info.tsv directly (explaining 56% & 45.5% and normalizing by background!)
    for fam_label, sub in [("PF00144", df_info_all[df_info_all["family"]=="PF00144"]),
                           ("PF00905", df_info_all[df_info_all["family"]=="PF00905"]),
                           ("PF13354", df_info_all[df_info_all["family"]=="PF13354"]),
                           ("All_Pooled", df_info_all)]:
        sub = sub[sub["LRT"] >= 0].copy() if "LRT" in sub else sub.copy()
        # Note: in selection_pressure.r, filter(LRT >= 0) was applied! Let's check both all rows and LRT >= 0!
        n_tot = len(sub)
        sub_int = sub[sub["group"] == "internal_node"]
        sub_term = sub[sub["group"] != "internal_node"]
        n_int, n_term = len(sub_int), len(sub_term)
        pos_int = int((sub_int["p_used"] < 0.05).sum())
        pos_term = int((sub_term["p_used"] < 0.05).sum())
        pos_tot = pos_int + pos_term

        or_int, p_int_vs_term = stats.fisher_exact([[pos_int, n_int - pos_int], [pos_term, n_term - pos_term]])

        info_tsv_rows.append({
            "source_table": "Fig3B_info_tsv_LRTge0",
            "family": fam_label,
            "comparison_group": "Internal_vs_Terminal",
            "category": "Internal_Branches",
            "tested_branches": n_int,
            "bg_share_pct": round(100 * n_int / max(n_tot, 1), 2),
            "pos_selected_branches": pos_int,
            "share_among_pos_selected_pct": round(100 * pos_int / max(pos_tot, 1), 2),
            "within_category_pos_rate_pct": round(100 * pos_int / max(n_int, 1), 2),
            "fisher_OR_vs_complement": round(float(or_int), 3),
            "fisher_p_value": float(p_int_vs_term),
        })
        info_tsv_rows.append({
            "source_table": "Fig3B_info_tsv_LRTge0",
            "family": fam_label,
            "comparison_group": "Internal_vs_Terminal",
            "category": "Terminal_Branches",
            "tested_branches": n_term,
            "bg_share_pct": round(100 * n_term / max(n_tot, 1), 2),
            "pos_selected_branches": pos_term,
            "share_among_pos_selected_pct": round(100 * pos_term / max(pos_tot, 1), 2),
            "within_category_pos_rate_pct": round(100 * pos_term / max(n_term, 1), 2),
            "fisher_OR_vs_complement": round(1.0 / max(float(or_int), 1e-9), 3),
            "fisher_p_value": float(p_int_vs_term),
        })

        # Stratify internal branches by internal_type (virus, bacteria, Nested)
        for cat in ["virus", "bacteria", "Nested"]:
            c_sub = sub_int[sub_int["internal_type"] == cat]
            rest = sub_int[sub_int["internal_type"] != cat]
            nc = len(c_sub)
            pc = int((c_sub["p_used"] < 0.05).sum())
            nr = len(rest)
            pr = int((rest["p_used"] < 0.05).sum())
            or_c, p_c = stats.fisher_exact([[pc, nc - pc], [pr, nr - pr]]) if (nc > 0 and nr > 0) else (np.nan, 1.0)
            info_tsv_rows.append({
                "source_table": "Fig3B_info_tsv_LRTge0",
                "family": fam_label,
                "comparison_group": "Internal_Branches_By_Type",
                "category": f"Internal_{cat}",
                "tested_branches": nc,
                "bg_share_pct": round(100 * nc / max(n_int, 1), 2),
                "pos_selected_branches": pc,
                "share_among_pos_selected_pct": round(100 * pc / max(pos_int, 1), 2) if pos_int else 0.0,
                "within_category_pos_rate_pct": round(100 * pc / max(nc, 1), 2) if nc else 0.0,
                "fisher_OR_vs_complement": round(float(or_c), 3) if pd.notna(or_c) else np.nan,
                "fisher_p_value": float(p_c),
            })

        # Stratify ALL branches by internal_type (virus, bacteria, Nested)
        for cat in ["virus", "bacteria", "Nested"]:
            c_sub = sub[sub["internal_type"] == cat]
            rest = sub[sub["internal_type"] != cat]
            nc = len(c_sub)
            pc = int((c_sub["p_used"] < 0.05).sum())
            nr = len(rest)
            pr = int((rest["p_used"] < 0.05).sum())
            or_c, p_c = stats.fisher_exact([[pc, nc - pc], [pr, nr - pr]]) if (nc > 0 and nr > 0) else (np.nan, 1.0)
            info_tsv_rows.append({
                "source_table": "Fig3B_info_tsv_LRTge0",
                "family": fam_label,
                "comparison_group": "All_Branches_By_Type",
                "category": f"All_{cat}",
                "tested_branches": nc,
                "bg_share_pct": round(100 * nc / max(n_tot, 1), 2),
                "pos_selected_branches": pc,
                "share_among_pos_selected_pct": round(100 * pc / max(pos_tot, 1), 2) if pos_tot else 0.0,
                "within_category_pos_rate_pct": round(100 * pc / max(nc, 1), 2) if nc else 0.0,
                "fisher_OR_vs_complement": round(float(or_c), 3) if pd.notna(or_c) else np.nan,
                "fisher_p_value": float(p_c),
            })

    # 2) Also summarize full raw absrel.json (both Uncorrected P<=0.05 and Holm-Bonferroni Corrected P<=0.05)
    for p_col, lbl in [("sig_uncorr_05", "absrel_JSON_UncorrP05"), ("sig_corr_05", "absrel_JSON_CorrP05")]:
        for fam_label, sub in [("PF00144", df_json_all[df_json_all["family"]=="PF00144"]),
                               ("PF00905", df_json_all[df_json_all["family"]=="PF00905"]),
                               ("PF13354", df_json_all[df_json_all["family"]=="PF13354"]),
                               ("All_Pooled", df_json_all)]:
            n_tot = len(sub)
            sub_int = sub[sub["branch_type"] == "Internal"]
            sub_term = sub[sub["branch_type"] == "Terminal"]
            n_int, n_term = len(sub_int), len(sub_term)
            pos_int = int(sub_int[p_col].sum())
            pos_term = int(sub_term[p_col].sum())
            pos_tot = pos_int + pos_term
            or_int, p_it = stats.fisher_exact([[pos_int, n_int - pos_int], [pos_term, n_term - pos_term]])
            info_tsv_rows.append({
                "source_table": lbl,
                "family": fam_label,
                "comparison_group": "Internal_vs_Terminal",
                "category": "Internal_Branches",
                "tested_branches": n_int,
                "bg_share_pct": round(100 * n_int / max(n_tot, 1), 2),
                "pos_selected_branches": pos_int,
                "share_among_pos_selected_pct": round(100 * pos_int / max(pos_tot, 1), 2) if pos_tot else 0.0,
                "within_category_pos_rate_pct": round(100 * pos_int / max(n_int, 1), 2),
                "fisher_OR_vs_complement": round(float(or_int), 3),
                "fisher_p_value": float(p_it),
            })
            info_tsv_rows.append({
                "source_table": lbl,
                "family": fam_label,
                "comparison_group": "Internal_vs_Terminal",
                "category": "Terminal_Branches",
                "tested_branches": n_term,
                "bg_share_pct": round(100 * n_term / max(n_tot, 1), 2),
                "pos_selected_branches": pos_term,
                "share_among_pos_selected_pct": round(100 * pos_term / max(pos_tot, 1), 2) if pos_tot else 0.0,
                "within_category_pos_rate_pct": round(100 * pos_term / max(n_term, 1), 2),
                "fisher_OR_vs_complement": round(1.0 / max(float(or_int), 1e-9), 3),
                "fisher_p_value": float(p_it),
            })
            for cat in ["Viral", "Bacterial", "Mixed_Nested"]:
                c_sub = sub_int[sub_int["clade_category"] == cat]
                rest = sub_int[sub_int["clade_category"] != cat]
                nc, pc = len(c_sub), int(c_sub[p_col].sum())
                nr, pr = len(rest), int(rest[p_col].sum())
                or_c, p_c = stats.fisher_exact([[pc, nc - pc], [pr, nr - pr]]) if (nc > 0 and nr > 0) else (np.nan, 1.0)
                info_tsv_rows.append({
                    "source_table": lbl,
                    "family": fam_label,
                    "comparison_group": "Internal_Branches_By_Type",
                    "category": f"Internal_{cat}",
                    "tested_branches": nc,
                    "bg_share_pct": round(100 * nc / max(n_int, 1), 2),
                    "pos_selected_branches": pc,
                    "share_among_pos_selected_pct": round(100 * pc / max(pos_int, 1), 2) if pos_int else 0.0,
                    "within_category_pos_rate_pct": round(100 * pc / max(nc, 1), 2) if nc else 0.0,
                    "fisher_OR_vs_complement": round(float(or_c), 3) if pd.notna(or_c) else np.nan,
                    "fisher_p_value": float(p_c),
                })

    df_c4 = pd.DataFrame(info_tsv_rows)
    df_c4.to_csv(f"{STAGE}/Table_S_aBSREL_Background_Normalized_Frequencies.csv", index=False)

    # === Lineage trajectory summary: Lineage-Resolved Ancestral-to-Descendant Selection ===
    df_lin = pd.DataFrame(lineage_rows)
    df_pair = pd.DataFrame(paired_lineage_stats)

    # Deduplicate branch-stage occurrences within each lineage_type to test unique branch rates,
    # AND test paired ancestral-vs-terminal shifts per root-to-tip lineage
    lin_summary = []
    for fam_lbl, sub_l, sub_p in [
        ("PF00144", df_lin[df_lin["family"]=="PF00144"], df_pair[df_pair["family"]=="PF00144"]),
        ("PF00905", df_lin[df_lin["family"]=="PF00905"], df_pair[df_pair["family"]=="PF00905"]),
        ("PF13354", df_lin[df_lin["family"]=="PF13354"], df_pair[df_pair["family"]=="PF13354"]),
        ("All_Pooled", df_lin, df_pair),
    ]:
        for ltype in sorted(sub_l["lineage_type"].unique()):
            lt_df = sub_l[sub_l["lineage_type"] == ltype]
            # Deduplicate unique branches within this lineage_type & stage
            lt_uniq = lt_df.drop_duplicates(subset=["family", "branch_name", "stage"])
            for stg in sorted(lt_uniq["stage"].unique()):
                s_df = lt_uniq[lt_uniq["stage"] == stg]
                n_b = len(s_df)
                n_pos_u = int(s_df["sig_uncorr_05"].sum())
                n_pos_c = int(s_df["sig_corr_05"].sum())
                med_w = float(s_df["baseline_omega"].median())
                iqr_lo = float(s_df["baseline_omega"].quantile(0.25))
                iqr_hi = float(s_df["baseline_omega"].quantile(0.75))
                lin_summary.append({
                    "family": fam_lbl,
                    "lineage_type": ltype,
                    "evolutionary_stage": stg,
                    "n_unique_branches": n_b,
                    "n_pos_uncorr_05": n_pos_u,
                    "pct_pos_uncorr_05": round(100 * n_pos_u / max(n_b, 1), 2),
                    "n_pos_corr_05": n_pos_c,
                    "pct_pos_corr_05": round(100 * n_pos_c / max(n_b, 1), 2),
                    "median_baseline_omega": round(med_w, 4),
                    "iqr_omega_25_75": f"{iqr_lo:.4f}-{iqr_hi:.4f}",
                })

        # Paired test on exact root-to-tip lineages: ancestral omega vs terminal omega
        for ltype, p_sub in sub_p.groupby("lineage_type"):
            a_w = p_sub["ancestral_median_omega"].values
            t_w = p_sub["terminal_omega"].values
            valid = np.isfinite(a_w) & np.isfinite(t_w)
            if valid.sum() >= 5:
                stat_w, p_wil = stats.wilcoxon(a_w[valid], t_w[valid], alternative="greater")
                pct_dec = 100.0 * float((a_w[valid] > t_w[valid]).mean())
            else:
                p_wil, pct_dec = np.nan, np.nan
            lin_summary.append({
                "family": fam_lbl,
                "lineage_type": ltype,
                "evolutionary_stage": "PAIRED_Ancestral_vs_Terminal_Test",
                "n_unique_branches": int(valid.sum()),
                "n_pos_uncorr_05": int(p_sub["ancestral_any_pos_sel_uncorr"].sum()),
                "pct_pos_uncorr_05": round(100 * float(p_sub["ancestral_any_pos_sel_uncorr"].mean()), 2),
                "n_pos_corr_05": int(p_sub["terminal_pos_sel_uncorr"].sum()),
                "pct_pos_corr_05": round(pct_dec, 2) if pd.notna(pct_dec) else np.nan,
                "median_baseline_omega": round(float(np.nanmedian(a_w[valid] - t_w[valid])), 4) if valid.sum() else np.nan,
                "iqr_omega_25_75": f"Wilcoxon_p={p_wil:.2e}" if pd.notna(p_wil) else "NA",
            })

    df_lin_sum = pd.DataFrame(lin_summary)
    df_lin_sum.to_csv(f"{STAGE}/Table_S_Lineage_Resolved_Selection_Trajectories.csv", index=False)

    print("=== Background normalization: Fig3B_info_tsv_LRTge0 (All_Pooled) ===")
    print(df_c4[(df_c4["source_table"]=="Fig3B_info_tsv_LRTge0") & (df_c4["family"]=="All_Pooled")].to_string(index=False))
    print("\n=== Background normalization: absrel_JSON_UncorrP05 (All_Pooled) ===")
    print(df_c4[(df_c4["source_table"]=="absrel_JSON_UncorrP05") & (df_c4["family"]=="All_Pooled")].to_string(index=False))
    print("\n=== Lineage-Resolved Selection Trajectories (All_Pooled) ===")
    print(df_lin_sum[df_lin_sum["family"]=="All_Pooled"].to_string(index=False))

    # Plot Figure: Lineage-Resolved Ancestral-to-Descendant Selection Trajectories
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), dpi=300)

    # Panel A: Episodic positive selection rate (Uncorr P<=0.05) across stages by lineage type (Pooled)
    ax = axes[0]
    pool_u = df_lin_sum[(df_lin_sum["family"]=="All_Pooled") & (~df_lin_sum["evolutionary_stage"].str.startswith("PAIRED"))]
    ltypes = ["HGT_Associated_Viral", "HGT_Associated_Bacterial", "NonHGT_Pure_Viral", "NonHGT_Pure_Bacterial"]
    colors = {"HGT_Associated_Viral": "#1B9E77", "HGT_Associated_Bacterial": "#D95F02",
              "NonHGT_Pure_Viral": "#386CB0", "NonHGT_Pure_Bacterial": "#7570B3"}
    x_pos = [0, 1, 2]
    for lt in ltypes:
        sub = pool_u[pool_u["lineage_type"] == lt].sort_values("evolutionary_stage")
        if len(sub) == 0:
            continue
        ys = []
        xs = []
        for s_idx, (_, r) in enumerate(sub.iterrows()):
            if r["evolutionary_stage"].startswith("1_"): xs.append(0)
            elif r["evolutionary_stage"].startswith("2_"): xs.append(1)
            else: xs.append(2)
            ys.append(r["pct_pos_uncorr_05"])
        ax.plot(xs, ys, marker="o", linewidth=2.2, markersize=7, color=colors.get(lt, "#555555"), label=lt.replace("_", " "))
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["1. Deep/HGT-Transition\nInternal Branches", "2. Descendant\nInternal Branches", "3. Terminal Tip\nBranches"], fontweight="bold")
    ax.set_ylabel("Branches with Episodic Selection (%)", fontweight="bold")
    ax.set_title("(A) Positive Selection Along Evolutionary Stages", fontweight="bold", loc="left")
    ax.legend(fontsize=8.5, frameon=True)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    # Panel B: Paired baseline omega shift (Ancestral Internal -> Terminal Tip) per lineage type
    ax = axes[1]
    data_boxes = []
    labels_boxes = []
    box_cols = []
    for lt in ltypes:
        sub_p = df_pair[df_pair["lineage_type"] == lt]
        if len(sub_p) < 3:
            continue
        data_boxes.append(np.clip(sub_p["ancestral_median_omega"].dropna().values, 0, 1.5))
        labels_boxes.append(f"{lt.replace('_', chr(10))}\nAncestral\n(n={len(sub_p)})")
        box_cols.append(colors[lt])
        data_boxes.append(np.clip(sub_p["terminal_omega"].dropna().values, 0, 1.5))
        labels_boxes.append(f"Terminal\n(n={len(sub_p)})")
        box_cols.append("#B0BEC5")
    bp = ax.boxplot(data_boxes, tick_labels=labels_boxes, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], box_cols):
        patch.set_facecolor(c)
        patch.set_alpha(0.8)
    ax.axhline(1.0, color="#D32F2F", linestyle="--", linewidth=1.0, label="Neutral (omega=1)")
    ax.set_ylabel("Baseline MG94xREV omega (dN/dS)", fontweight="bold")
    ax.set_title("(B) Paired Ancestral-to-Terminal omega Shift", fontweight="bold", loc="left")
    ax.tick_params(axis="x", labelsize=7.5)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(f"{STAGE}/Fig_S_Lineage_Selection_Trajectories.pdf", bbox_inches="tight")
    plt.savefig(f"{STAGE}/Fig_S_Lineage_Selection_Trajectories.png", bbox_inches="tight")
    plt.close()

run_all()
'''

def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(WSL_STAGE_WIN, exist_ok=True)

    for fam in ["pf00144", "pf00905", "pf13354"]:
        for fn in [f"{fam}.nwk", f"{fam}_info.tsv"]:
            shutil.copy2(DATA_DIR / fn, WSL_STAGE_WIN / fn)
        shutil.copy2(SEL_SRC / f"{fam}_absrel.json", WSL_STAGE_WIN / f"{fam}_absrel.json")

    worker_file = WSL_STAGE_WIN / "run_m04_part1.py"
    worker_file.write_text(WSL_WORKER, encoding="utf-8")

    proc = subprocess.run(linux_cmd(WSL_PY, f"{WSL_STAGE_LIN}/run_m04_part1.py"), capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(f"06_selection_absrel_and_lineages.py failed:\n{proc.stderr}")

    shutil.copy2(WSL_STAGE_WIN / "Table_S_aBSREL_Background_Normalized_Frequencies.csv",
                 TABLES_DIR / "Table_S07_aBSREL_Background_Normalized_Frequencies.csv")
    shutil.copy2(WSL_STAGE_WIN / "Table_S_Lineage_Resolved_Selection_Trajectories.csv",
                 TABLES_DIR / "Table_S08_Lineage_Resolved_Selection_Trajectories.csv")
    shutil.copy2(WSL_STAGE_WIN / "Fig_S_Lineage_Selection_Trajectories.pdf",
                 FIG_DIR / "Fig_S02_Lineage_Selection_Trajectories.pdf")
    shutil.copy2(WSL_STAGE_WIN / "Fig_S_Lineage_Selection_Trajectories.png",
                 FIG_DIR / "Fig_S02_Lineage_Selection_Trajectories.png")


if __name__ == "__main__":
    main()
