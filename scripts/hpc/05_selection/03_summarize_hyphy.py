#!/usr/bin/env python3
"""Selection analysis, step 3: summarize HyPhy aBSREL, MEME, FEL & GARD outputs.

Parses completed HyPhy JSON outputs across hyphy_tasks.tsv and writes:
  1. PF00144_2x2_Alignment_And_SampleSize_Benchmark_2026.tsv
  2. aBSREL_Background_Normalized_Summary_2026.tsv
  3. aBSREL_Lineage_Trajectories_2026.tsv
  4. Subclade_GARD_And_HyPhy_Summary_2026.tsv
  5. Updates PerSite_Domain_Selection_Metrics_2026.tsv & Catalytic_Motif_Purifying_Proof_2026.tsv
  6. Stage4_Selection_Summary_Report.txt
"""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from ete4 import Tree
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

OUT_DIR = hpc.BASE / "05_selection_pressure"
HYPHY_DIR = OUT_DIR / "hyphy_sets"
TASKS_TSV = OUT_DIR / "hyphy_tasks.tsv"
SAT_TSV = OUT_DIR / "dS_Saturation_And_Xia_Test_2026.tsv"
PERSITE_TSV = OUT_DIR / "PerSite_Domain_Selection_Metrics_2026.tsv"
PROOF_TSV = OUT_DIR / "Catalytic_Motif_Purifying_Proof_2026.tsv"


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def parse_meme(meme_json: Path) -> tuple[pd.DataFrame, dict]:
    data = json.loads(meme_json.read_text(encoding="utf-8"))
    headers = [h[0] for h in data["MLE"]["headers"]]
    content = data["MLE"]["content"]["0"]
    df = pd.DataFrame(content, columns=headers)
    df["hyphy_site_1based"] = np.arange(1, len(df) + 1)
    pvals = df["p-value"].astype(float).values
    df["bh_q_value"] = bh_fdr(pvals)
    bl_tot = sum(float(v.get("total branch length", 0.0)) for v in data.get("branch attributes", {}).values() if isinstance(v, dict))
    # In HyPhy JSON, total tree length is also under fits -> Global MG94xREV -> tree length
    fits = data.get("fits", {})
    tree_len = np.nan
    for mk in ("Global MG94xREV", "Baseline MG94xREV", "Standard MG94"):
        if mk in fits and "tree length" in fits[mk]:
            tree_len = float(fits[mk]["tree length"])
            break
    if math.isnan(tree_len):
        br0 = data.get("branch attributes", {}).get("0", {})
        tree_len = sum(float(b.get("Global MG94xREV", 0.0)) for b in br0.values() if isinstance(b, dict))

    summary = {
        "meme_tree_length": round(float(tree_len), 3) if not math.isnan(tree_len) else np.nan,
        "meme_sites_p10": int((pvals < 0.10).sum()),
        "meme_sites_p05": int((pvals < 0.05).sum()),
        "meme_sites_p01": int((pvals < 0.01).sum()),
        "meme_sites_q10": int((df["bh_q_value"].values < 0.10).sum()),
        "meme_pct_p05": round(100.0 * float((pvals < 0.05).mean()), 2) if len(pvals) else 0.0,
    }
    return df, summary


def parse_fel(fel_json: Path) -> tuple[pd.DataFrame, dict]:
    data = json.loads(fel_json.read_text(encoding="utf-8"))
    headers = [h[0] for h in data["MLE"]["headers"]]
    content = data["MLE"]["content"]["0"]
    df = pd.DataFrame(content, columns=headers)
    df["hyphy_site_1based"] = np.arange(1, len(df) + 1)
    pvals = df["p-value"].astype(float).values
    alpha = df["alpha"].astype(float).values
    beta = df["beta"].astype(float).values
    df["bh_q_value"] = bh_fdr(pvals)
    pos_05 = int(((pvals < 0.05) & (beta > alpha)).sum())
    pur_05 = int(((pvals < 0.05) & (alpha > beta)).sum())
    summary = {
        "fel_pos_p05": pos_05,
        "fel_purifying_p05": pur_05,
        "fel_purifying_pct": round(100.0 * pur_05 / max(len(df), 1), 2),
    }
    return df, summary


def parse_gard(gard_json: Path) -> dict:
    """Breakpoints of the best GARD model (segments in breakpointData - 1) and its AICc gain over one tree.

    gard.status (written by 02_run_hyphy.sh) says whether the search finished; GARD rewrites gard.json during
    the search, so a crashed or timed-out run leaves a partial file whose breakpoint count is a lower bound.
    """
    status_p = gard_json.with_name("gard.status")
    status = status_p.read_text(encoding="utf-8").strip() if status_p.is_file() else "no_status"
    empty = {"gard_status": status, "gard_breakpoints": np.nan, "gard_delta_cAIC": np.nan}
    if not gard_json.is_file():
        return empty
    try:
        data = json.loads(gard_json.read_text(encoding="utf-8"))
        segs = data.get("breakpointData", {})
        n_bp = max(len(segs) - 1, 0) if isinstance(segs, (list, dict)) else 0
        d_aic = float(data["singleTreeAICc"]) - float(data["bestModelAICc"])
        return {"gard_status": status, "gard_breakpoints": int(n_bp), "gard_delta_cAIC": round(d_aic, 2)}
    except Exception:
        return empty


def weighted_omega(rate_dist) -> tuple[float, float]:
    if not rate_dist:
        return np.nan, 0.0
    w_sum = p_gt1 = 0.0
    for om, prop in rate_dist:
        w_sum += min(float(om), 100.0) * float(prop)
        if float(om) > 1.0:
            p_gt1 += float(prop)
    return w_sum, p_gt1


def parse_absrel(
    absrel_json: Path, id_map_tsv: Path, family: str, set_id: str
) -> tuple[list[dict], list[dict], dict]:
    data = json.loads(absrel_json.read_text(encoding="utf-8"))
    id_dom = {}
    with open(id_map_tsv, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            id_dom[r["hyphy_id"]] = r["domain"]

    br_attrs = data["branch attributes"]["0"]
    raw_tree = data["input"]["trees"]["0"].strip()
    if not raw_tree.endswith(";"):
        raw_tree += ";"
    nwk_clean = re.sub(r"\{[^}]+\}", "", raw_tree).replace("'", "")
    tree = Tree(nwk_clean, parser=1)

    for nd in tree.traverse("postorder"):
        if nd.is_leaf:
            d = id_dom.get(nd.name, "Bacterial")
            nd.add_prop("n_vir", 1 if d == "Viral" else 0)
            nd.add_prop("n_bac", 0 if d == "Viral" else 1)
            nd.add_prop("dom_class", d)
        else:
            nv = sum(c.props["n_vir"] for c in nd.children)
            nb = sum(c.props["n_bac"] for c in nd.children)
            nd.add_prop("n_vir", nv)
            nd.add_prop("n_bac", nb)
            if nv > 0 and nb > 0:
                nd.add_prop("dom_class", "Mixed_Nested")
            elif nv > 0:
                nd.add_prop("dom_class", "Viral")
            else:
                nd.add_prop("dom_class", "Bacterial")

    branch_rows = []
    for nd in tree.traverse():
        if nd.up is None:
            continue
        bname = nd.name or ""
        attr = br_attrs.get(bname) or br_attrs.get(bname.strip("'"))
        if attr is None:
            continue
        p_corr = float(attr.get("Corrected P-value", 1.0))
        p_unc = float(attr.get("Uncorrected P-value", 1.0))
        base_om = float(attr.get("Baseline MG94xREV omega ratio", np.nan))
        w_om, p_gt1 = weighted_omega(attr.get("Rate Distributions", []))
        branch_rows.append({
            "set_id": set_id,
            "family": family,
            "branch_name": bname,
            "branch_type": "Terminal" if nd.is_leaf else "Internal",
            "clade_category": nd.props["dom_class"],
            "p_corrected": p_corr,
            "p_uncorrected": p_unc,
            "sig_corr_05": p_corr <= 0.05,
            "sig_uncorr_05": p_unc <= 0.05,
            "baseline_omega": min(base_om, 10.0) if not math.isnan(base_om) else np.nan,
            "weighted_omega": w_om,
            "prop_omega_gt1": p_gt1,
        })

    # Root-to-tip longitudinal trajectories
    traj_rows = []
    for tip in tree.leaves():
        t_dom = id_dom.get(tip.name, "Bacterial")
        path = []
        cur = tip
        while cur is not None and cur.up is not None:
            path.append(cur)
            cur = cur.up
        path = list(reversed(path))
        if len(path) < 2:
            continue
        pure_stem_idx = len(path) - 1
        while pure_stem_idx > 0 and path[pure_stem_idx - 1].props["dom_class"] in ("Viral", "Bacterial"):
            pure_stem_idx -= 1
        parent_of_stem = path[pure_stem_idx].up
        is_hgt_lin = False
        if parent_of_stem is not None and parent_of_stem.up is not None:
            pv, pb = parent_of_stem.props["n_vir"], parent_of_stem.props["n_bac"]
            if t_dom == "Viral" and pb >= 2 and pb / max(pv + pb, 1) >= 0.30:
                is_hgt_lin = True
            elif t_dom == "Bacterial" and pv >= 2 and pv / max(pv + pb, 1) >= 0.30:
                is_hgt_lin = True
        lin_type = f"HGT_Associated_{t_dom}" if is_hgt_lin else f"NonHGT_Pure_{t_dom}"

        for idx, nd in enumerate(path):
            attr = br_attrs.get(nd.name or "")
            if attr is None:
                continue
            if nd.is_leaf:
                stage = "3_Terminal_Tip"
            elif is_hgt_lin:
                stage = "1_Ancestral_or_HGT_Transition" if idx <= pure_stem_idx else "2_PostHGT_Descendant_Internal"
            else:
                stage = "1_Deep_Ancestral_Internal" if idx < max(1, (len(path) - 1) // 2) else "2_Recent_Descendant_Internal"
            bom = float(attr.get("Baseline MG94xREV omega ratio", np.nan))
            traj_rows.append({
                "set_id": set_id,
                "family": family,
                "tip_id": tip.name,
                "lineage_type": lin_type,
                "stage": stage,
                "baseline_omega": min(bom, 10.0) if not math.isnan(bom) else np.nan,
                "p_uncorrected": float(attr.get("Uncorrected P-value", 1.0)),
                "p_corrected": float(attr.get("Corrected P-value", 1.0)),
                "sig_uncorr_05": float(attr.get("Uncorrected P-value", 1.0)) <= 0.05,
                "sig_corr_05": float(attr.get("Corrected P-value", 1.0)) <= 0.05,
                "prop_omega_gt1": weighted_omega(attr.get("Rate Distributions", []))[1],
            })

    summary = {
        "absrel_tested_branches": len(branch_rows),
        "absrel_sig_corr_05": sum(1 for r in branch_rows if r["sig_corr_05"]),
        "absrel_sig_uncorr_05": sum(1 for r in branch_rows if r["sig_uncorr_05"]),
        "absrel_pct_corr_05": round(100.0 * sum(1 for r in branch_rows if r["sig_corr_05"]) / max(len(branch_rows), 1), 2),
    }
    return branch_rows, traj_rows, summary


def _nonempty_json(p: Path) -> bool:
    return p.is_file() and p.stat().st_size > 100


def main() -> int:
    if not TASKS_TSV.is_file():
        raise FileNotFoundError(f"Missing {TASKS_TSV}")
    tasks_df = pd.read_csv(TASKS_TSV, sep="\t")
    sat_df = pd.read_csv(SAT_TSV, sep="\t")
    sat_map = sat_df.set_index("dataset").to_dict("index")

    bench_rows = []
    subclade_rows = []
    all_branches = []
    all_trajs = []
    site_meme_by_fam = defaultdict(list)
    site_fel_by_fam = defaultdict(list)

    for _, t in tasks_df.iterrows():
        sid = t["set_id"]
        fam = t["family"]
        sgrp = t["set_group"]
        sdir = Path(t["set_dir"])
        ab_p, me_p, fe_p, ga_p = sdir / "absrel.json", sdir / "meme.json", sdir / "fel.json", sdir / "gard.json"
        if not (_nonempty_json(me_p) and _nonempty_json(fe_p)):
            print(f"[SKIP incomplete MEME/FEL] {sid}")
            continue

        col_map = pd.read_csv(sdir / "col_map.tsv", sep="\t")
        df_meme, meme_sum = parse_meme(me_p)
        df_fel, fel_sum = parse_fel(fe_p)
        gard_sum = parse_gard(ga_p) if (sgrp == "subclade" and _nonempty_json(ga_p)) else {}
        if _nonempty_json(ab_p):
            br_rows, tr_rows, ab_sum = parse_absrel(ab_p, sdir / "id_map.tsv", fam, sid)
        else:
            print(f"[PARTIAL: MEME+FEL ready, aBSREL pending] {sid}")
            br_rows, tr_rows = [], []
            ab_sum = {
                "absrel_tested_branches": np.nan,
                "absrel_sig_corr_05": np.nan,
                "absrel_sig_uncorr_05": np.nan,
                "absrel_pct_corr_05": np.nan,
            }

        all_branches.extend(br_rows)
        all_trajs.extend(tr_rows)

        # Map site p-values back to original HMM match column
        df_m = df_meme.merge(col_map, on="hyphy_site_1based", how="inner")
        df_m["set_id"] = sid
        df_f = df_fel.merge(col_map, on="hyphy_site_1based", how="inner")
        df_f["set_id"] = sid
        if sgrp in ("benchmark_N90", "subclade"):
            site_meme_by_fam[fam].append(df_m[["orig_col_1based", "p-value", "bh_q_value"]])
            site_fel_by_fam[fam].append(df_f[["orig_col_1based", "alpha", "beta", "p-value"]])

        sinfo = sat_map.get(sid, {})
        row_combined = {
            "set_id": sid,
            "family": fam,
            "set_group": sgrp,
            "track": t["track"],
            "n_seqs": int(t["n_seqs"]),
            "n_codons": int(t["n_codons"]),
            "median_dS": sinfo.get("median_dS", np.nan),
            "pos3_Iss": sinfo.get("pos3_Iss", np.nan),
            "pos3_Iss_c_sym": sinfo.get("pos3_Iss_c_sym", np.nan),
            "pos3_unsaturated_sym": sinfo.get("pos3_unsaturated_sym", False),
            "PHI_p": sinfo.get("PHI_p", np.nan),
            **ab_sum,
            **meme_sum,
            **fel_sum,
            **gard_sum,
        }
        if sgrp in ("benchmark_N90", "benchmark_2x2"):
            bench_rows.append(row_combined)
        else:
            subclade_rows.append(row_combined)

    df_bench = pd.DataFrame(bench_rows)
    df_sub = pd.DataFrame(subclade_rows)
    df_br = pd.DataFrame(all_branches)
    df_tr = pd.DataFrame(all_trajs)

    df_bench.to_csv(OUT_DIR / "PF00144_2x2_Alignment_And_SampleSize_Benchmark_2026.tsv", sep="\t", index=False)
    df_sub.to_csv(OUT_DIR / "Subclade_GARD_And_HyPhy_Summary_2026.tsv", sep="\t", index=False)

    # Background-normalized aBSREL branch summary (Internal vs Terminal x Viral/Bacterial/Mixed_Nested)
    bg_rows = []
    if len(df_br) > 0:
        for (fam, sgrp), sub in df_br.merge(tasks_df[["set_id", "set_group"]], on="set_id").groupby(["family", "set_group"]):
            for btype in ("Internal", "Terminal", "All"):
                sb = sub if btype == "All" else sub[sub["branch_type"] == btype]
                for ccat in ("Viral", "Bacterial", "Mixed_Nested", "All"):
                    sc = sb if ccat == "All" else sb[sb["clade_category"] == ccat]
                    if len(sc) == 0:
                        continue
                    bg_rows.append({
                        "family": fam,
                        "set_group": sgrp,
                        "branch_type": btype,
                        "clade_category": ccat,
                        "n_tested_branches": len(sc),
                        "n_sig_corr_05": int(sc["sig_corr_05"].sum()),
                        "pct_sig_corr_05": round(100.0 * float(sc["sig_corr_05"].mean()), 2),
                        "n_sig_uncorr_05": int(sc["sig_uncorr_05"].sum()),
                        "pct_sig_uncorr_05": round(100.0 * float(sc["sig_uncorr_05"].mean()), 2),
                        "median_baseline_omega": round(float(sc["baseline_omega"].median()), 4),
                    })
    df_bg = pd.DataFrame(bg_rows)
    df_bg.to_csv(OUT_DIR / "aBSREL_Background_Normalized_Summary_2026.tsv", sep="\t", index=False)

    # Lineage Trajectories Summary
    traj_sum_rows = []
    if len(df_tr) > 0:
        for (fam, ltype, stage), sub in df_tr.groupby(["family", "lineage_type", "stage"]):
            traj_sum_rows.append({
                "family": fam,
                "lineage_type": ltype,
                "stage": stage,
                "n_branches": len(sub),
                "median_baseline_omega": round(float(sub["baseline_omega"].median()), 4),
                "iqr_baseline_omega": f"{sub['baseline_omega'].quantile(0.25):.3f}-{sub['baseline_omega'].quantile(0.75):.3f}",
                "pct_sig_uncorr_05": round(100.0 * float(sub["sig_uncorr_05"].mean()), 2),
                "pct_sig_corr_05": round(100.0 * float(sub["sig_corr_05"].mean()), 2),
                "mean_prop_omega_gt1": round(float(sub["prop_omega_gt1"].mean()), 4),
            })
    df_traj_sum = pd.DataFrame(traj_sum_rows)
    df_traj_sum.to_csv(OUT_DIR / "aBSREL_Lineage_Trajectories_2026.tsv", sep="\t", index=False)

    # Enrich PerSite_Domain_Selection_Metrics_2026.tsv with subclade/benchmark MEME & FEL p-values
    if PERSITE_TSV.is_file():
        df_ps = pd.read_csv(PERSITE_TSV, sep="\t")
        meme_agg_list = []
        fel_agg_list = []
        for fam, dfs in site_meme_by_fam.items():
            dm = pd.concat(dfs, ignore_index=True).groupby("orig_col_1based", as_index=False).agg(
                meme_min_p=("p-value", "min"),
                meme_median_p=("p-value", "median"),
                meme_n_sets_p05=("p-value", lambda s: int((s < 0.05).sum())),
            )
            dm["family"] = fam
            dm["hmm_col"] = dm["orig_col_1based"]
            meme_agg_list.append(dm[["family", "hmm_col", "meme_min_p", "meme_median_p", "meme_n_sets_p05"]])
        for fam, dfs in site_fel_by_fam.items():
            dfel = pd.concat(dfs, ignore_index=True).groupby("orig_col_1based", as_index=False).agg(
                fel_min_p=("p-value", "min"),
                fel_mean_alpha=("alpha", "mean"),
                fel_mean_beta=("beta", "mean"),
            )
            dfel["family"] = fam
            dfel["hmm_col"] = dfel["orig_col_1based"]
            fel_agg_list.append(dfel[["family", "hmm_col", "fel_min_p", "fel_mean_alpha", "fel_mean_beta"]])
        # Drop old columns if re-running
        for col in ("meme_min_p", "meme_median_p", "meme_n_sets_p05", "fel_min_p", "fel_mean_alpha", "fel_mean_beta"):
            if col in df_ps.columns:
                df_ps = df_ps.drop(columns=[col])
        if meme_agg_list:
            df_ps = df_ps.merge(pd.concat(meme_agg_list, ignore_index=True), on=["family", "hmm_col"], how="left")
        if fel_agg_list:
            df_ps = df_ps.merge(pd.concat(fel_agg_list, ignore_index=True), on=["family", "hmm_col"], how="left")
        df_ps.to_csv(PERSITE_TSV, sep="\t", index=False)

        # Update Catalytic_Motif_Purifying_Proof_2026.tsv with MEME depletion check
        if PROOF_TSV.is_file() and "meme_min_p" in df_ps.columns:
            df_pr = pd.read_csv(PROOF_TSV, sep="\t")
            dep_rows = {}
            for fam_k in list(df_pr["family"].unique()):
                sub_p = df_ps if fam_k == "All_Pooled" else df_ps[df_ps["family"] == fam_k]
                sub_p = sub_p.dropna(subset=["meme_min_p"])
                mot_p = sub_p[sub_p["in_catalytic_motif_core"]]
                bkg_p = sub_p[~sub_p["in_catalytic_motif_core"]]
                m_pos = int((mot_p["meme_min_p"] < 0.05).sum())
                b_pos = int((bkg_p["meme_min_p"] < 0.05).sum())
                _, p_dep = stats.fisher_exact(
                    [[m_pos, len(mot_p) - m_pos], [b_pos, len(bkg_p) - b_pos]],
                    alternative="less",
                ) if (len(mot_p) and len(bkg_p)) else (np.nan, np.nan)
                dep_rows[fam_k] = (
                    round(100.0 * m_pos / max(len(mot_p), 1), 2),
                    round(100.0 * b_pos / max(len(bkg_p), 1), 2),
                    float(p_dep),
                )
            df_pr["motif_meme_pos_pct"] = df_pr["family"].map(lambda k: dep_rows.get(k, (np.nan, np.nan, np.nan))[0])
            df_pr["bg_meme_pos_pct"] = df_pr["family"].map(lambda k: dep_rows.get(k, (np.nan, np.nan, np.nan))[1])
            df_pr["Fisher_p_MEME_motif_depletion"] = df_pr["family"].map(lambda k: dep_rows.get(k, (np.nan, np.nan, np.nan))[2])
            df_pr.to_csv(PROOF_TSV, sep="\t", index=False)

    lines = [
        "=== Stage 4A: Selection Pressure & Saturation Summary ===",
        "",
        "--- 1. Benchmarks (N=90 across 4 families + PF00144 2x2 factorial) ---",
        df_bench.to_string(index=False) if len(df_bench) else "(pending)",
        "",
        "--- 2. Low-Divergence Subclades (GARD + HyPhy) ---",
        df_sub.to_string(index=False) if len(df_sub) else "(pending)",
        "",
        "--- 3. aBSREL Background-Normalized Summary (subclades) ---",
        df_bg[df_bg["set_group"] == "subclade"].to_string(index=False) if len(df_bg) else "(pending)",
    ]
    rep = "\n".join(lines) + "\n"
    (OUT_DIR / "Stage4_Selection_Summary_Report.txt").write_text(rep, encoding="utf-8")
    print("\n" + rep)
    return 0


if __name__ == "__main__":
    sys.exit(main())
