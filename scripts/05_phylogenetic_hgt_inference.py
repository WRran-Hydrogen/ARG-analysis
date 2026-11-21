"""Module 03: Deterministic, Reproducible Phylogenetic HGT Detection & Sensitivity Analysis.

Addresses Reviewer Comments 1, 2, and 10:
1. Implements an explicit, deterministic phylogenetic HGT inference framework modeled after
   Irwin et al. (Nature Microbiology 2022) and Dmitrijeva et al. (Nature Ecology & Evolution 2024):
   - Roots each family tree deterministically on the primary inter-domain backbone split
     maximizing Viral-vs-Bacterial clade separation (preventing midpoint-rooting artifacts
     where an entire ancient kingdom backbone is mislabeled as an HGT clade).
   - Collapses unresolved zero-length polytomy nodes (dist <= 1e-6 & UFBoot <= 1%) so
     unresolved sister tips are not artificially inflated into separate HGT events.
   - Combines Fitch maximum-parsimony ancestral domain-state reconstruction with explicit
     double-outgroup topological nesting criteria: focal clade C of recipient domain D1
     (>=80% purity, <25% of total tree tips, non-root-child) nested within sister S1 and
     cousin/outgroup S2 dominated by donor domain D2 (donor >=1 in both S1 and S2,
     combined outgroup donor fraction >=75%).
   - Stratifies HGT events by UFBoot support (>=70%, >=80%, >=90%) on the bounding nesting
     nodes (parent P and grandparent G) and focal monophyly node.
   - Explicitly separates Independent HGT Transition Events from Descendant Nested Branches & Tips.
2. Performs Comment 10 Sensitivity Analyses by pruning:
   - Bacterial Tier A contaminants (61 SARG-misclassified Viruses, 153 artificial/synthetic, header phage),
   - Bacterial Tier B putative cryptic prophages (NCBI nr hits >=90% aa identity to viral SBLs),
   - Viral Tier 3 borderline fragments (retaining Module 02 Tier 1 + Tier 2 viral loci).
"""

from pathlib import Path
import os
import shutil
import subprocess
import pandas as pd

import sbl_config as cfg
from sbl_config import linux_cmd

TABLES_DIR = str(cfg.TABLES_DIR)
SUPP_DIR = str(cfg.SUPP_DIR)
DATA_DIR = str(cfg.PHYLO_DIR)
_stage, WSL_STAGE_LIN = cfg.stage_dirs("wsl_m03")
WSL_STAGE_WIN = str(_stage)

NW_DIR = DATA_DIR
BAC_AUDIT_TSV = str(cfg.SUPP_DIR / "PerSequence_Bacterial_Audit.tsv")
VIR_AUDIT_TSV = str(cfg.SUPP_DIR / "Supplementary_Data_Viral_Genomic_Neighborhoods.tsv")

WSL_PY = cfg.WSL_PY

WSL_WORKER_CODE = r'''
import os
import pandas as pd
from ete4 import Tree

STAGE_DIR = os.path.dirname(os.path.abspath(__file__))
FAMILIES = ["PF00144", "PF00905", "PF13354"]

def tip_domain(name):
    return "Viral" if "IMGVR" in name else "Bacterial"

def find_primary_interdomain_root(tree):
    """Find the internal node whose bipartition maximizes Youden's J (|TPR_viral - FPR_bacterial|)."""
    all_leaves = list(tree.leaves())
    tot_v = max(1, sum(1 for l in all_leaves if tip_domain(l.name) == "Viral"))
    tot_b = max(1, len(all_leaves) - tot_v)
    tot = len(all_leaves)

    for nd in tree.traverse("postorder"):
        if nd.is_leaf:
            d = tip_domain(nd.name)
            nd.add_prop("_v", 1 if d == "Viral" else 0)
            nd.add_prop("_b", 1 if d == "Bacterial" else 0)
        else:
            nd.add_prop("_v", sum(c.props["_v"] for c in nd.children))
            nd.add_prop("_b", sum(c.props["_b"] for c in nd.children))

    best_node = None
    best_score = -1.0
    for nd in tree.traverse():
        if nd.up is None:
            continue
        v_in = nd.props["_v"]
        b_in = nd.props["_b"]
        n_in = v_in + b_in
        n_out = tot - n_in
        if n_in < 10 or n_out < 10:
            continue
        # Youden's J separation between Viral and Bacterial domains across the split
        score = abs((v_in / tot_v) - (b_in / tot_b))
        if score > best_score:
            best_score = score
            best_node = nd
    return best_node

def collapse_unresolved_polytomies(tree, eps=1e-6):
    """Collapse internal branches that have ~0 branch length and no bootstrap support (<=1%)."""
    to_delete = []
    for node in tree.traverse("postorder"):
        if not node.is_leaf and node.up is not None:
            ub = float(node.props.get("ufboot", 0.0))
            dist = float(node.dist) if node.dist is not None else 0.0
            if dist <= eps and ub <= 1.0:
                to_delete.append(node)
    for node in to_delete:
        node.delete(prevent_nondicotomic=False)

def fitch_parsimony_transitions(tree):
    """Run two-pass Fitch parsimony for binary state {Viral, Bacterial}."""
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            node.add_prop("fitch_set", {tip_domain(node.name)})
        else:
            child_sets = [c.props["fitch_set"] for c in node.children]
            inter = set.intersection(*child_sets)
            if inter:
                node.add_prop("fitch_set", inter)
            else:
                node.add_prop("fitch_set", set.union(*child_sets))

    root = tree
    # At root: assign each root child independently if root fitch_set is ambiguous,
    # otherwise assign root_state
    if len(root.props["fitch_set"]) == 1:
        root.add_prop("asr_state", next(iter(root.props["fitch_set"])))
    else:
        # Ambiguous root between primary Viral and Bacterial sides: mark root as "Root_Split"
        root.add_prop("asr_state", "Root_Split")

    transitions = []
    for node in tree.traverse("preorder"):
        if node is root:
            continue
        p_state = node.up.props["asr_state"]
        fset = node.props["fitch_set"]
        if p_state == "Root_Split":
            # Primary children of root establish their own clade backbone domain
            if len(fset) == 1:
                state = next(iter(fset))
            else:
                state = "Viral" if node.props["n_vir"] >= node.props["n_bac"] else "Bacterial"
            node.add_prop("asr_state", state)
            continue
        if p_state in fset:
            state = p_state
        else:
            # Break tie toward majority tip domain inside this clade
            if len(fset) == 1:
                state = next(iter(fset))
            else:
                state = "Viral" if node.props["n_vir"] >= node.props["n_bac"] else "Bacterial"
        node.add_prop("asr_state", state)
        if state != p_state:
            transitions.append((node, p_state, state))
    return transitions

def detect_nested_hgt(tree, fam, dataset_name):
    """Identify formal HGT events combining Fitch ASR transitions with double-outgroup nesting."""
    tot_tree_tips = len(list(tree.leaves()))
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            d = tip_domain(node.name)
            node.add_prop("n_vir", 1 if d == "Viral" else 0)
            node.add_prop("n_bac", 1 if d == "Bacterial" else 0)
        else:
            node.add_prop("n_vir", sum(c.props["n_vir"] for c in node.children))
            node.add_prop("n_bac", sum(c.props["n_bac"] for c in node.children))

    transitions = fitch_parsimony_transitions(tree)
    events = []
    ev_idx = 0

    for node, from_dom, to_dom in transitions:
        n_v = node.props["n_vir"]
        n_b = node.props["n_bac"]
        n_tot = n_v + n_b
        n_recip = n_v if to_dom == "Viral" else n_b
        purity = n_recip / max(n_tot, 1)

        parent = node.up
        grandparent = parent.up if parent is not None else None
        is_root_child = (parent is tree)

        # Sister S1 (other children of parent)
        s1_v, s1_b = 0, 0
        if parent is not None:
            for sib in parent.children:
                if sib is not node:
                    s1_v += sib.props["n_vir"]
                    s1_b += sib.props["n_bac"]
        s1_tot = s1_v + s1_b
        s1_donor = s1_v if from_dom == "Viral" else s1_b
        s1_donor_frac = s1_donor / max(s1_tot, 1)

        # Second outgroup S2 (other children of grandparent)
        s2_v, s2_b = 0, 0
        if grandparent is not None:
            for aunt in grandparent.children:
                if aunt is not parent:
                    s2_v += aunt.props["n_vir"]
                    s2_b += aunt.props["n_bac"]
        s2_tot = s2_v + s2_b
        s2_donor = s2_v if from_dom == "Viral" else s2_b
        s2_donor_frac = s2_donor / max(s2_tot, 1)

        bg_tot = s1_tot + s2_tot
        bg_donor = s1_donor + s2_donor
        bg_donor_frac = bg_donor / max(bg_tot, 1)

        tot_recip_in_tree = tree.props["n_vir"] if to_dom == "Viral" else tree.props["n_bac"]
        # Strict Irwin et al. / Dmitrijeva et al. nesting criteria:
        # 1) Not a direct root-split backbone child, focal clade < 25% of tree tips,
        #    and contains < 50% of all recipient-domain sequences in the tree
        # 2) Focal recipient purity >= 80%
        # 3) Sister S1 and Second-Outgroup S2 both contain donor sequences and combined outgroup >= 75% donor
        is_nested_topology = (
            (not is_root_child)
            and (n_tot < 0.25 * tot_tree_tips)
            and (n_recip < 0.50 * max(tot_recip_in_tree, 1))
            and (purity >= 0.80)
            and (s1_donor >= 1 and s2_donor >= 1)
            and (bg_donor_frac >= 0.75)
        )

        sup_focal = float(node.props.get("ufboot", 100.0)) if not node.is_leaf else 100.0
        sup_parent = float(parent.props.get("ufboot", 0.0)) if (parent is not None and not parent.is_leaf) else 0.0
        sup_grand = float(grandparent.props.get("ufboot", 0.0)) if (grandparent is not None and not grandparent.is_leaf) else 0.0
        max_nesting_sup = max(sup_parent, sup_grand)
        min_nesting_sup = min(sup_parent, sup_grand) if grandparent is not None else sup_parent

        n_desc_branches = sum(1 for _ in node.traverse())
        leaves_list = [l.name for l in node.leaves() if tip_domain(l.name) == to_dom]

        ev_idx += 1
        events.append({
            "dataset": dataset_name,
            "family": fam,
            "event_id": f"{fam}_{dataset_name}_HGT{ev_idx:03d}",
            "direction": f"{from_dom}_to_{to_dom}",
            "donor_domain": from_dom,
            "recipient_domain": to_dom,
            "n_focal_tips_total": n_tot,
            "n_focal_recipient_tips": n_recip,
            "focal_recipient_purity": round(purity, 4),
            "n_descendant_branches": n_desc_branches,
            "s1_donor_n": s1_donor,
            "s1_total_n": s1_tot,
            "s1_donor_frac": round(s1_donor_frac, 4),
            "s2_donor_n": s2_donor,
            "s2_total_n": s2_tot,
            "s2_donor_frac": round(s2_donor_frac, 4),
            "combined_outgroup_donor_frac": round(bg_donor_frac, 4),
            "passes_double_outgroup_nesting": is_nested_topology,
            "focal_node_ufboot": sup_focal,
            "parent_node_ufboot": sup_parent,
            "grandparent_node_ufboot": sup_grand,
            "max_bounding_ufboot": max_nesting_sup,
            "min_bounding_ufboot": min_nesting_sup,
            "focal_is_singleton_tip": node.is_leaf,
            "representative_tips": ";".join(leaves_list[:5]),
            "all_recipient_tips": ";".join(leaves_list),
        })
    return pd.DataFrame(events)

def run_analysis():
    bac_audit = pd.read_csv(f"{STAGE_DIR}/PerSequence_Bacterial_Audit.tsv", sep="\t")
    vir_audit = pd.read_csv(f"{STAGE_DIR}/Supplementary_Data_Viral_Genomic_Neighborhoods.tsv", sep="\t")

    bac_tierA_keep = set(bac_audit[bac_audit["pass_CleanStrict_TierA"]]["seq_id"].astype(str))
    bac_tierB_keep = set(bac_audit[bac_audit["pass_CleanCore_TierB"]]["seq_id"].astype(str))
    vir_tier12_uvigs = set(vir_audit[vir_audit["viral_origin_tier"].isin(["Tier1_HighConfidence_Viral", "Tier2_Supported_Viral"])]["uvig"].astype(str))

    all_events = []
    summary_rows = []
    sens_rows = []

    for fam in FAMILIES:
        raw_nw = open(f"{STAGE_DIR}/{fam}.nwk").read()
        scenarios = [
            ("Baseline_Published_Tree", None, None),
            ("Sensitivity_A_BacCleanStrict", bac_tierA_keep, None),
            ("Sensitivity_B_BacCleanCore_NoNRProphage", bac_tierB_keep, None),
            ("Sensitivity_C_BacCleanCore_And_ViralTier12", bac_tierB_keep, vir_tier12_uvigs),
        ]

        for scen_name, bac_keep, vir_keep in scenarios:
            t = Tree(raw_nw, parser=0)
            for nd in t.traverse():
                if not nd.is_leaf:
                    nd.add_prop("ufboot", float(nd.support) if nd.support is not None else 100.0)
                else:
                    nd.add_prop("ufboot", 100.0)
                if "support" in nd.props:
                    nd.del_prop("support")
            if "dist" in t.props:
                t.del_prop("dist")

            if bac_keep is not None or vir_keep is not None:
                keep_leaves = []
                for l in t.leaves():
                    if tip_domain(l.name) == "Bacterial":
                        if bac_keep is None or l.name in bac_keep:
                            keep_leaves.append(l.name)
                    else:
                        uvig = l.name.split("|")[0]
                        if vir_keep is None or uvig in vir_keep:
                            keep_leaves.append(l.name)
                t.prune(keep_leaves, preserve_branch_length=True)

            # Root deterministically on the primary inter-domain backbone split
            backbone_node = find_primary_interdomain_root(t)
            if backbone_node is not None:
                t.set_outgroup(backbone_node)

            collapse_unresolved_polytomies(t)

            n_leaves = len(list(t.leaves()))
            n_vir_tips = sum(1 for l in t.leaves() if tip_domain(l.name) == "Viral")
            n_bac_tips = n_leaves - n_vir_tips

            ev_df = detect_nested_hgt(t, fam, scen_name)
            all_events.append(ev_df)

            if len(ev_df) > 0:
                nested = ev_df[ev_df["passes_double_outgroup_nesting"]]
                nested_b2v = nested[nested["direction"] == "Bacterial_to_Viral"]
                nested_v2b = nested[nested["direction"] == "Viral_to_Bacterial"]
                sup70 = nested[nested["max_bounding_ufboot"] >= 70]
                sup80 = nested[nested["max_bounding_ufboot"] >= 80]
                sup90 = nested[nested["max_bounding_ufboot"] >= 90]
                strict_p70 = nested[nested["parent_node_ufboot"] >= 70]
            else:
                nested = nested_b2v = nested_v2b = sup70 = sup80 = sup90 = strict_p70 = ev_df

            row = {
                "family": fam,
                "scenario": scen_name,
                "total_tips": n_leaves,
                "viral_tips": n_vir_tips,
                "bacterial_tips": n_bac_tips,
                "all_fitch_transitions": len(ev_df),
                "nested_hgt_events_total": len(nested),
                "nested_B_to_V_events": len(nested_b2v),
                "nested_V_to_B_events": len(nested_v2b),
                "nested_recipient_tips_total": int(nested["n_focal_recipient_tips"].sum()) if len(nested) else 0,
                "nested_B_to_V_viral_tips": int(nested_b2v["n_focal_recipient_tips"].sum()) if len(nested_b2v) else 0,
                "nested_V_to_B_bac_tips": int(nested_v2b["n_focal_recipient_tips"].sum()) if len(nested_v2b) else 0,
                "nested_descendant_branches_total": int(nested["n_descendant_branches"].sum()) if len(nested) else 0,
                "hgt_events_bounding_UFBoot_ge70": len(sup70),
                "hgt_events_bounding_UFBoot_ge80": len(sup80),
                "hgt_events_bounding_UFBoot_ge90": len(sup90),
                "hgt_events_parent_UFBoot_ge70": len(strict_p70),
                "tips_in_UFBoot_ge70_events": int(sup70["n_focal_recipient_tips"].sum()) if len(sup70) else 0,
            }
            if scen_name == "Baseline_Published_Tree":
                summary_rows.append(row)
            sens_rows.append(row)

    df_all_ev = pd.concat(all_events, ignore_index=True)
    df_sum = pd.DataFrame(summary_rows)
    df_sens = pd.DataFrame(sens_rows)

    num_cols = [c for c in df_sum.columns if c not in ("family", "scenario")]
    pooled = {"family": "All_Pooled", "scenario": "Baseline_Published_Tree"}
    for c in num_cols:
        pooled[c] = int(df_sum[c].sum())
    df_sum = pd.concat([df_sum, pd.DataFrame([pooled])], ignore_index=True)

    pooled_sens = []
    for scen in df_sens["scenario"].unique():
        sub = df_sens[df_sens["scenario"] == scen]
        p = {"family": "All_Pooled", "scenario": scen}
        for c in num_cols:
            p[c] = int(sub[c].sum())
        pooled_sens.append(p)
    df_sens = pd.concat([df_sens, pd.DataFrame(pooled_sens)], ignore_index=True)

    df_all_ev.to_csv(f"{STAGE_DIR}/Supplementary_Data_HGT_Events_Catalog.tsv", sep="\t", index=False)
    df_sum.to_csv(f"{STAGE_DIR}/Table_S_HGT_Summary_By_Criteria.csv", index=False)
    df_sens.to_csv(f"{STAGE_DIR}/Table_S_HGT_Prophage_Pruning_Sensitivity.csv", index=False)

    print("=== Baseline Rigorous HGT Detection Summary (Backbone-Rooted, Polytomy-Collapsed) ===")
    print(df_sum.to_string(index=False))
    print("\n=== Comment 10 Sensitivity (Pooled Across Families) ===")
    print(df_sens[df_sens["family"] == "All_Pooled"].to_string(index=False))
    print("\n=== Comment 10 Sensitivity (By Family) ===")
    print(df_sens[df_sens["family"] != "All_Pooled"][["family", "scenario", "bacterial_tips", "nested_hgt_events_total", "nested_B_to_V_events", "nested_V_to_B_events", "hgt_events_bounding_UFBoot_ge70", "nested_recipient_tips_total"]].to_string(index=False))

run_analysis()
'''

def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(SUPP_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(WSL_STAGE_WIN, exist_ok=True)

    shutil.copy2(BAC_AUDIT_TSV, os.path.join(WSL_STAGE_WIN, "PerSequence_Bacterial_Audit.tsv"))
    shutil.copy2(VIR_AUDIT_TSV, os.path.join(WSL_STAGE_WIN, "Supplementary_Data_Viral_Genomic_Neighborhoods.tsv"))
    for fam in ["PF00144", "PF00905", "PF13354"]:
        shutil.copy2(os.path.join(NW_DIR, f"{fam}.nwk"), os.path.join(WSL_STAGE_WIN, f"{fam}.nwk"))

    worker_path = os.path.join(WSL_STAGE_WIN, "run_hgt_detection.py")
    with open(worker_path, "w", encoding="utf-8") as f:
        f.write(WSL_WORKER_CODE)

    proc = subprocess.run(linux_cmd(WSL_PY, f"{WSL_STAGE_LIN}/run_hgt_detection.py"), capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(f"HGT detection failed:\n{proc.stderr}")

    shutil.copy2(os.path.join(WSL_STAGE_WIN, "Supplementary_Data_HGT_Events_Catalog.tsv"),
                 os.path.join(SUPP_DIR, "Supplementary_Data_HGT_Events_Catalog.tsv"))
    shutil.copy2(os.path.join(WSL_STAGE_WIN, "Table_S_HGT_Summary_By_Criteria.csv"),
                 os.path.join(TABLES_DIR, "Table_S05_HGT_Summary_By_Criteria.csv"))
    shutil.copy2(os.path.join(WSL_STAGE_WIN, "Table_S_HGT_Prophage_Pruning_Sensitivity.csv"),
                 os.path.join(TABLES_DIR, "Table_S06_HGT_Prophage_Pruning_Sensitivity.csv"))


if __name__ == "__main__":
    main()
