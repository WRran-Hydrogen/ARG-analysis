#!/usr/bin/env python3
"""Pre-registered phylogenetic HGT inference & OAT robustness pipeline.

Implements Section 4.1 of docs/METHODS_AND_PIPELINE.md on HPC Phase 3 IQ-TREE outputs:
  1. Rooting across three unbiased rooting families:
       - Primary: Sister-PRE biological outgroup (_out.treefile)
       - MAD (Minimal Ancestor Deviation; Tria et al. 2017) on ingroup-only (_in.treefile)
       - MinVar (Minimum Variance; Mai et al. 2017) & Midpoint on ingroup-only (_in.treefile)
     Plus Youden's J recorded strictly as an appendix comparator (excluded from cross-validation).
  2. Unresolved near-zero polytomy collapsing (dist <= 1e-6 & ufboot <= 1.0).
  3. Concordant Ancestral State Reconstruction:
       - Two-pass Fitch maximum parsimony
       - Exact F81 continuous-time Markov Marginal ML ASR (Felsenstein pruning + ML rate mu,
         equivalent to PastML F81; also invokes pastml package if installed for audit).
  4. Irwin et al. (2022) double-outgroup topological nesting + strict Tier-1 gating:
       - Tier-1 HGT: double-outgroup nested + Fitch/ML concordant + min(P, G) UFBoot >= 95 &
         SH-aLRT >= 80 + singleton viral tips must pass Phase 2 Strict Tier1_HighConfidence_Viral
         + not Prophage_Borne.
       - Prophage_Borne: V->B where 100% of bacterial recipients are in_bacterial_prophage=True
         (reported separately, excluded from primary V->B counts).
       - Functional Stratification: Putative_Hydrolytic_BetaLactamase vs. PBP_or_PRE_Homolog,
         plus expanded viral_cluster_size reporting.
  5. OAT Robustness Recurrence across all completed sensitivity trees (Jaccard >= 0.70):
       - Root_Invariant = Sister_PRE + MAD + (MinVar or Midpoint)
       - Robust_Core_HGT = OAT recurrence >= 80%
  6. CL0013 Superfamily + LOO outgroup rooting validation (when CL0013 trees exist).
"""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from ete4 import Tree
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

csv.field_size_limit(10_000_000)

CURATED_DIR = hpc.BASE / "01_final_curated"
P3_DIR = hpc.PHASE3_DIR
ALN_DIR = P3_DIR / "alignments"
TREES_DIR = P3_DIR / "trees"
OUT_DIR = hpc.BASE / "04_hgt_inference"

FAMILIES = ("PF00144", "PF00905", "PF13354", "PF00768", "PF13354_PF00768")
CL0013_TREES = (
    ("CL0013.struct_c20.no_PF27310", "outgroup_struct_no_PF27310.txt"),
    ("CL0013.struct_c20.only_PF04960", "outgroup_struct_only_PF04960.txt"),
    ("CL0013.struct_c20.only_PF06089", "outgroup_struct_only_PF06089.txt"),
    ("CL0013.all3_out.mafft_out", "outgroup_taxa_all3.txt"),
    ("CL0013.loo_no_PF04960.mafft_out", "outgroup_taxa_loo_no_PF04960.txt"),
    ("CL0013.loo_no_PF06089.mafft_out", "outgroup_taxa_loo_no_PF06089.txt"),
    ("CL0013.loo_no_PF27310.mafft_out", "outgroup_taxa_loo_no_PF27310.txt"),
)


# ---------------------------------------------------------------------------
# Metadata & Tree I/O Helpers
# ---------------------------------------------------------------------------
def as_bool(v) -> bool:
    return str(v).strip().lower() in ("true", "t", "1", "yes")


def load_catalog_metadata() -> tuple[dict[str, dict], dict[str, dict]]:
    v_meta: dict[str, dict] = {}
    with open(CURATED_DIR / "Final_Viral_SBL_Catalog.tsv", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            v_meta[row["uid"]] = {
                "origin": "Viral",
                "functional_tier": row.get("functional_tier", "PBP_or_PRE_Homolog"),
                "viral_origin_tier": row.get("viral_origin_tier", ""),
                "viral_cluster_size": int(float(row.get("viral_cluster_size") or 1)),
                "in_bacterial_prophage": False,
                "winning_acc": row.get("winning_acc", ""),
                "card_pident": float(row.get("card_pident") or 0.0),
                "card_cov": float(row.get("card_cov") or 0.0),
                "card_gene_family": row.get("card_gene_family", "NONE"),
            }

    b_meta: dict[str, dict] = {}
    with open(CURATED_DIR / "Final_Bacterial_TripleSource_Catalog.tsv", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            b_meta[row["uid"]] = {
                "origin": "Bacterial",
                "functional_tier": row.get("functional_tier", "PBP_or_PRE_Homolog"),
                "viral_origin_tier": "",
                "viral_cluster_size": 1,
                "in_bacterial_prophage": as_bool(row.get("in_bacterial_prophage", "False")),
                "winning_acc": row.get("winning_acc", ""),
                "card_pident": float(row.get("card_pident") or 0.0),
                "card_cov": float(row.get("card_cov") or 0.0),
                "card_gene_family": row.get("card_gene_family", "NONE"),
            }
    return v_meta, b_meta


def tip_domain(uid: str, v_meta: dict[str, dict]) -> str:
    if uid in v_meta or uid.startswith("IMGVR_") or uid.startswith("MOSAiC"):
        return "Viral"
    return "Bacterial"


def parse_iqtree_support(label: str | None) -> tuple[float, float]:
    """Parse IQ-TREE 'SH-aLRT/UFBoot' node label; empty/unparseable labels default to (0.0, 0.0)."""
    if not label:
        return (0.0, 0.0)
    s = str(label).strip()
    if "/" in s:
        a, b = s.split("/", 1)
        try:
            return (float(a), float(b))
        except ValueError:
            return (0.0, 0.0)
    try:
        val = float(s)
        return (val, val)
    except ValueError:
        return (0.0, 0.0)


def _has_valid_iqtree_label(label: str | None) -> bool:
    if not label:
        return False
    s = str(label).strip()
    if not s:
        return False
    if "/" in s:
        a, b = s.split("/", 1)
        try:
            float(a)
            float(b)
            return True
        except ValueError:
            return False
    try:
        float(s)
        return True
    except ValueError:
        return False


def verify_tree_bootstrap_complete(tpath: Path, nw_str: str) -> bool:
    """Reject unfinished/unbootstrapped IQ-TREE files instead of inventing support."""
    iq_path = Path(str(tpath)[:-len(".treefile")] + ".iqtree") if str(tpath).endswith(".treefile") else None
    if iq_path is not None and iq_path.is_file():
        tail = iq_path.read_text(encoding="utf-8", errors="ignore")[-4096:]
        if "Analysis results written to" not in tail:
            print(f"  [SKIP incomplete .iqtree] {tpath.name}", flush=True)
            return False
    t = Tree(nw_str, parser=1)
    internal_nonroot = [
        nd for nd in t.traverse()
        if not nd.is_leaf and nd.up is not None and float(nd.dist or 0.0) > 1e-6
    ]
    if len(internal_nonroot) >= 10:
        n_labeled = sum(1 for nd in internal_nonroot if _has_valid_iqtree_label(nd.name))
        if n_labeled / len(internal_nonroot) < 0.80:
            print(
                f"  [SKIP unbootstrapped treefile] {tpath.name} "
                f"(labeled={n_labeled}/{len(internal_nonroot)})",
                flush=True,
            )
            return False
    return True


def read_outgroup_ids(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def extract_unrooted_split_supports(tree: Tree, leaf_subset: set[str] | None = None) -> tuple[dict[int, tuple[float, float]], dict[str, int], int]:
    """Map every unrooted bipartition of `leaf_subset` (bitmask) to its IQ-TREE (sh_alrt, ufboot)."""
    target_leaves = sorted(l.name for l in tree.leaves() if leaf_subset is None or l.name in leaf_subset)
    bit_map = {name: i for i, name in enumerate(target_leaves)}
    all_mask = (1 << len(target_leaves)) - 1

    split_sup: dict[int, tuple[float, float]] = {}
    for nd in tree.traverse("postorder"):
        if nd.is_leaf:
            m = (1 << bit_map[nd.name]) if nd.name in bit_map else 0
        else:
            m = 0
            for c in nd.children:
                m |= c.props["_bmask"]
        nd.add_prop("_bmask", m)
        if not nd.is_leaf and nd.up is not None and 0 < m < all_mask:
            key = min(m, all_mask ^ m)
            sh, uf = parse_iqtree_support(nd.name)
            # If pruning collapsed an adjacent outgroup-only side into the same ingroup split,
            # retain the highest-supported representative split
            if key not in split_sup or uf > split_sup[key][1]:
                split_sup[key] = (sh, uf)
    return split_sup, bit_map, all_mask


def restore_split_supports(tree: Tree, split_sup: dict[int, tuple[float, float]], bit_map: dict[str, int], all_mask: int) -> None:
    """Re-attach canonical split supports (sh_alrt, ufboot) onto rooted tree nodes."""
    for nd in tree.traverse("postorder"):
        if nd.is_leaf:
            m = 1 << bit_map[nd.name]
            nd.add_prop("_bmask", m)
            nd.add_prop("sh_alrt", 100.0)
            nd.add_prop("ufboot", 100.0)
        else:
            m = 0
            for c in nd.children:
                m |= c.props["_bmask"]
            nd.add_prop("_bmask", m)
            if nd.up is None:
                nd.add_prop("sh_alrt", 0.0)
                nd.add_prop("ufboot", 0.0)
            else:
                key = min(m, all_mask ^ m)
                sh, uf = split_sup.get(key, (0.0, 0.0))
                nd.add_prop("sh_alrt", sh)
                nd.add_prop("ufboot", uf)


def write_rooted_treefile(tree: Tree, out_path: Path) -> None:
    """Write rooted Newick with SH-aLRT/UFBoot internal labels (before polytomy collapsing)."""
    t_copy = tree.copy()
    for nd in t_copy.traverse():
        if nd.is_leaf:
            continue
        if nd.up is None:
            nd.name = ""
        else:
            sh = float(nd.props.get("sh_alrt", 0.0))
            uf = float(nd.props.get("ufboot", 0.0))
            nd.name = f"{sh:.1f}/{uf:.1f}"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(t_copy.write(parser=1) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Rooting Algorithms: Sister-PRE, MAD, MinVar, Midpoint, Youden's J
# ---------------------------------------------------------------------------
def root_sister_pre(nw_str: str, og_ids: set[str]) -> tuple[Tree, dict]:
    t = Tree(nw_str, parser=1)
    all_leaves = {l.name for l in t.leaves()}
    og_present = og_ids & all_leaves
    in_leaves = all_leaves - og_present
    if not og_present:
        raise RuntimeError("No outgroup taxa found in tree")

    split_sup, bit_map, all_mask = extract_unrooted_split_supports(t, leaf_subset=in_leaves)

    # Identify best split separating outgroup from ingroup
    n_og = len(og_present)
    n_tot = len(all_leaves)
    for nd in t.traverse("postorder"):
        if nd.is_leaf:
            nd.add_prop("_og", 1 if nd.name in og_present else 0)
            nd.add_prop("_sz", 1)
        else:
            nd.add_prop("_og", sum(c.props["_og"] for c in nd.children))
            nd.add_prop("_sz", sum(c.props["_sz"] for c in nd.children))

    best_node = None
    best_f1 = -1.0
    best_on_og_side = True
    for nd in t.traverse():
        if nd.up is None:
            continue
        k_og = nd.props["_og"]
        k_sz = nd.props["_sz"]
        # Check side A (nd)
        if k_og > 0:
            f1_a = 2.0 * k_og / (k_sz + n_og)
            if f1_a > best_f1:
                best_f1 = f1_a
                best_node = nd
                best_on_og_side = True
        # Check complement side B (everything outside nd)
        comp_og = n_og - k_og
        comp_sz = n_tot - k_sz
        if comp_og > 0:
            f1_b = 2.0 * comp_og / (comp_sz + n_og)
            if f1_b > best_f1:
                best_f1 = f1_b
                best_node = nd
                best_on_og_side = False

    if best_node is not None:
        t.set_outgroup(best_node)
    # Remove outgroup taxa while preserving ingroup branch lengths
    t.prune(sorted(in_leaves), preserve_branch_length=True)
    restore_split_supports(t, split_sup, bit_map, all_mask)
    return t, {"root_f1_outgroup": round(best_f1, 4), "outgroup_n": len(og_present)}


def _compute_patristic_and_dfs(tree: Tree) -> tuple[list[Tree], np.ndarray, np.ndarray, list[tuple[Tree, int, int, float, float]]]:
    """Assign DFS order 0..N-1 to leaves, build NxN patristic matrix D, and return node spans."""
    leaves: list[Tree] = []
    node_span: dict[int, tuple[int, int]] = {}

    def assign_dfs(nd: Tree) -> tuple[int, int]:
        if nd.is_leaf:
            idx = len(leaves)
            leaves.append(nd)
            node_span[id(nd)] = (idx, idx + 1)
            return (idx, idx + 1)
        lo = len(leaves)
        for c in nd.children:
            assign_dfs(c)
        hi = len(leaves)
        node_span[id(nd)] = (lo, hi)
        return (lo, hi)

    assign_dfs(tree)
    n = len(leaves)

    # Root-to-node depth relative to current root
    depth_map: dict[int, float] = {id(tree): 0.0}
    for nd in tree.traverse("preorder"):
        if nd.up is not None:
            depth_map[id(nd)] = depth_map[id(nd.up)] + max(0.0, float(nd.dist or 0.0))

    x_leaf = np.array([depth_map[id(l)] for l in leaves], dtype=np.float64)

    # Compute NxN patristic distance matrix in O(N^2) via sibling subtree Cartesian updates
    D = np.zeros((n, n), dtype=np.float64)
    for nd in tree.traverse("postorder"):
        if nd.is_leaf or len(nd.children) < 2:
            continue
        du = depth_map[id(nd)]
        spans = [node_span[id(c)] for c in nd.children]
        for a_i in range(len(spans)):
            s1, t1 = spans[a_i]
            if s1 >= t1:
                continue
            for b_i in range(a_i + 1, len(spans)):
                s2, t2 = spans[b_i]
                if s2 >= t2:
                    continue
                blk = (x_leaf[s1:t1, None] - du) + (x_leaf[None, s2:t2] - du)
                D[s1:t1, s2:t2] = blk
                D[s2:t2, s1:t1] = blk.T

    edges = []
    for nd in tree.traverse("preorder"):
        if nd.up is not None:
            s, t = node_span[id(nd)]
            L = max(0.0, float(nd.dist or 0.0))
            edges.append((nd, s, t, L, depth_map[id(nd)]))
    return leaves, x_leaf, D, edges


def _prefix_2d(M: np.ndarray) -> np.ndarray:
    n = M.shape[0]
    P = np.zeros((n + 1, n + 1), dtype=np.float64)
    P[1:, 1:] = np.cumsum(np.cumsum(M, axis=0), axis=1)
    return P


def _rect_sum(P: np.ndarray, s: int, t: int) -> float:
    return float(P[t, t] - P[s, t] - P[t, s] + P[s, s])


def root_mad(nw_str: str) -> tuple[Tree, dict]:
    """Exact vectorized Minimal Ancestor Deviation (MAD; Tria et al. 2017) rooting."""
    t = Tree(nw_str, parser=1)
    split_sup, bit_map, all_mask = extract_unrooted_split_supports(t)
    leaves, x_leaf, D, edges = _compute_patristic_and_dfs(t)
    n = len(leaves)

    # Zero-safe inverse distance matrices
    with np.errstate(divide="ignore", invalid="ignore"):
        D_inv = np.where(D > 1e-12, 1.0 / D, 0.0)
        D_pos = np.where(D > 1e-12, 1.0, 0.0)
    D_inv2 = D_inv * D_inv

    # Precompute 1D row sums and 2D prefix sums for O(1) per-edge MAD evaluation
    matrices = [
        D_pos,
        D_inv,
        x_leaf[:, None] * D_inv,
        D_inv2,
        x_leaf[:, None] * D_inv2,
        (x_leaf[:, None] ** 2) * D_inv2,
    ]
    row_csum = [np.concatenate([[0.0], np.cumsum(M.sum(axis=1))]) for M in matrices]
    pref2d = [_prefix_2d(M) for M in matrices]

    def cross_sum(k: int, s: int, t_end: int) -> float:
        return float(row_csum[k][t_end] - row_csum[k][s]) - _rect_sum(pref2d[k], s, t_end)

    ads: list[tuple[float, Tree, float]] = []
    for nd, s, t_end, L, y_v in edges:
        n_a = t_end - s
        n_b = n - n_a
        if n_a == 0 or n_b == 0:
            continue
        cnt = cross_sum(0, s, t_end)
        if cnt <= 0.0:
            continue
        s_inv = cross_sum(1, s, t_end)
        s_xinv = cross_sum(2, s, t_end) - y_v * s_inv
        s_inv2 = cross_sum(3, s, t_end)
        s_xinv2 = cross_sum(4, s, t_end) - y_v * s_inv2
        s_x2inv2 = cross_sum(5, s, t_end) - 2.0 * y_v * cross_sum(4, s, t_end) + (y_v * y_v) * s_inv2

        if L > 1e-12 and s_inv2 > 0.0:
            num = 2.0 * L * (s_inv - 2.0 * s_xinv2)
            den = 4.0 * L * L * s_inv2
            alpha = min(1.0, max(0.0, num / den))
        else:
            alpha = 0.5
            num = 0.0
            den = 0.0

        rss = cnt - 4.0 * s_xinv + 4.0 * s_x2inv2 - 2.0 * alpha * num + (alpha * alpha) * den
        ad = math.sqrt(max(0.0, rss) / cnt)
        ads.append((ad, nd, alpha))

    ads.sort(key=lambda z: z[0])
    best_ad, best_nd, best_alpha = ads[0]
    second_ad = ads[1][0] if len(ads) > 1 else best_ad
    rai = best_ad / max(second_ad, 1e-12)

    orig_L = max(0.0, float(best_nd.dist or 0.0))
    t.set_outgroup(best_nd)
    for c in t.children:
        if c is best_nd:
            c.dist = best_alpha * orig_L
        else:
            c.dist = (1.0 - best_alpha) * orig_L
    restore_split_supports(t, split_sup, bit_map, all_mask)
    return t, {"mad_score": round(best_ad, 5), "mad_rai": round(rai, 4)}


def root_minvar_or_midpoint(nw_str: str, method: str) -> tuple[Tree, dict]:
    """Exact O(N) FastRoot MinVar (Mai et al. 2017) or Midpoint rooting."""
    t = Tree(nw_str, parser=1)
    split_sup, bit_map, all_mask = extract_unrooted_split_supports(t)
    n = len(list(t.leaves()))

    # Bottom-up postorder DP
    for nd in t.traverse("postorder"):
        if nd.is_leaf:
            nd.add_prop("_n_in", 1)
            nd.add_prop("_sum_in", 0.0)
            nd.add_prop("_sq_in", 0.0)
            nd.add_prop("_h_in", 0.0)
        else:
            n_in = 0
            sum_in = 0.0
            sq_in = 0.0
            h_in = 0.0
            for c in nd.children:
                L = max(0.0, float(c.dist or 0.0))
                cn = c.props["_n_in"]
                cs = c.props["_sum_in"]
                cq = c.props["_sq_in"]
                ch = c.props["_h_in"]
                n_in += cn
                sum_in += cs + cn * L
                sq_in += cq + 2.0 * L * cs + cn * L * L
                h_in = max(h_in, ch + L)
            nd.add_prop("_n_in", n_in)
            nd.add_prop("_sum_in", sum_in)
            nd.add_prop("_sq_in", sq_in)
            nd.add_prop("_h_in", h_in)

    # Top-down preorder DP
    t.add_prop("_sum_out", 0.0)
    t.add_prop("_sq_out", 0.0)
    t.add_prop("_h_out", 0.0)

    best_var = float("inf")
    best_mp_span = float("inf")
    best_nd = None
    best_alpha = 0.5

    for p in t.traverse("preorder"):
        if p.is_leaf:
            continue
        p_sum_all = p.props["_sum_out"] + p.props["_sum_in"]
        p_sq_all = p.props["_sq_out"] + p.props["_sq_in"]
        # Top-2 heights among siblings + up_height for O(N) Midpoint
        cand_h = [(-1, p.props["_h_out"])]
        for c in p.children:
            cand_h.append((id(c), c.props["_h_in"] + max(0.0, float(c.dist or 0.0))))
        cand_h.sort(key=lambda x: x[1], reverse=True)

        for c in p.children:
            L = max(0.0, float(c.dist or 0.0))
            n_a = c.props["_n_in"]
            n_b = n - n_a
            if n_a == 0 or n_b == 0:
                continue
            s_a = c.props["_sum_in"]
            q_a = c.props["_sq_in"]
            s_p_other = p_sum_all - (s_a + n_a * L)
            q_p_other = p_sq_all - (q_a + 2.0 * L * s_a + n_a * L * L)

            c.add_prop("_sum_out", s_p_other + n_b * L)
            c.add_prop("_sq_out", q_p_other + 2.0 * L * s_p_other + n_b * L * L)

            h_a = c.props["_h_in"]
            h_p_other = cand_h[0][1] if cand_h[0][0] != id(c) else cand_h[1][1]
            c.add_prop("_h_out", h_p_other + L)

            if method == "MinVar":
                # Variance(alpha) for root at distance alpha*L from c, (1-alpha)*L from p:
                # mean = [s_a + s_p_other + n_b*L + (n_a - n_b)*alpha*L] / n
                # E[Z^2] = [q_a + q_p_other + 2*L*s_p_other + n_b*L^2
                #          + 2*L*(s_a - s_p_other - n_b*L)*alpha + n*L^2*alpha^2] / n
                c0 = s_a + s_p_other + n_b * L
                c1 = (n_a - n_b) * L
                d0 = q_a + q_p_other + 2.0 * L * s_p_other + n_b * L * L
                d1 = 2.0 * L * (s_a - s_p_other - n_b * L)
                d2 = n * L * L

                a_q = (d2 / n) - (c1 / n) ** 2
                b_q = (d1 / n) - 2.0 * (c0 * c1) / (n * n)
                if a_q > 1e-18:
                    alpha = min(1.0, max(0.0, -b_q / (2.0 * a_q)))
                else:
                    alpha = 0.5
                var_val = max(0.0, (d0 + d1 * alpha + d2 * alpha * alpha) / n - ((c0 + c1 * alpha) / n) ** 2)
                if var_val < best_var:
                    best_var = var_val
                    best_nd = c
                    best_alpha = alpha
            else:  # Midpoint
                # Root at alpha*L from c balances h_a + alpha*L == h_p_other + (1-alpha)*L
                if L > 1e-12:
                    alpha = (h_p_other - h_a + L) / (2.0 * L)
                else:
                    alpha = 0.5
                if 0.0 <= alpha <= 1.0:
                    max_dist = max(h_a + alpha * L, h_p_other + (1.0 - alpha) * L)
                    if max_dist < best_mp_span:
                        best_mp_span = max_dist
                        best_nd = c
                        best_alpha = alpha

    if best_nd is None:
        # Fallback if Midpoint lands on an internal zero-length boundary
        best_nd = max((c for c in t.traverse() if c.up is not None), key=lambda x: x.props["_h_in"])
        best_alpha = 0.5

    orig_L = max(0.0, float(best_nd.dist or 0.0))
    t.set_outgroup(best_nd)
    for c in t.children:
        if c is best_nd:
            c.dist = best_alpha * orig_L
        else:
            c.dist = (1.0 - best_alpha) * orig_L
    restore_split_supports(t, split_sup, bit_map, all_mask)
    metric = best_var if method == "MinVar" else best_mp_span
    return t, {"root_metric": round(metric, 5)}


def root_youden_j(nw_str: str, v_meta: dict[str, dict]) -> tuple[Tree, dict]:
    """Appendix comparator ONLY (excluded from primary cross-validation due to circularity)."""
    t = Tree(nw_str, parser=1)
    split_sup, bit_map, all_mask = extract_unrooted_split_supports(t)
    leaves = list(t.leaves())
    tot_v = max(1, sum(1 for l in leaves if tip_domain(l.name, v_meta) == "Viral"))
    tot_b = max(1, len(leaves) - tot_v)
    tot = len(leaves)

    for nd in t.traverse("postorder"):
        if nd.is_leaf:
            is_v = tip_domain(nd.name, v_meta) == "Viral"
            nd.add_prop("_v", 1 if is_v else 0)
            nd.add_prop("_b", 0 if is_v else 1)
        else:
            nd.add_prop("_v", sum(c.props["_v"] for c in nd.children))
            nd.add_prop("_b", sum(c.props["_b"] for c in nd.children))

    best_node = None
    best_j = -1.0
    for nd in t.traverse():
        if nd.up is None:
            continue
        n_in = nd.props["_v"] + nd.props["_b"]
        if n_in < 10 or (tot - n_in) < 10:
            continue
        score = abs((nd.props["_v"] / tot_v) - (nd.props["_b"] / tot_b))
        if score > best_j:
            best_j = score
            best_node = nd
    if best_node is not None:
        t.set_outgroup(best_node)
    restore_split_supports(t, split_sup, bit_map, all_mask)
    return t, {"youden_j": round(best_j, 4)}


def extract_superfamily_crown_splits() -> dict[str, dict]:
    """Derive basal crown-group bipartition anchors (S1, S2) per family from CL0013 superfamily tree."""
    meta_files = [
        ALN_DIR / "CL0013" / "cl0013_taxa_metadata.tsv",
        ALN_DIR / "CL0013" / "cl0013_struct_taxa_metadata.tsv",
    ]
    taxa_fam: dict[str, str] = {}
    for mp in meta_files:
        if mp.is_file():
            with open(mp, encoding="utf-8") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    taxa_fam[r["uid"]] = r["cl0013_family"]
    if not taxa_fam:
        return {}

    candidates = [
        ("CL0013.struct_c20.no_PF27310", "outgroup_struct_no_PF27310.txt"),
        ("CL0013.loo_no_PF27310.mafft_out", "outgroup_taxa_loo_no_PF27310.txt"),
    ]
    chosen = None
    for tname, og_fn in candidates:
        tp = TREES_DIR / "CL0013" / f"{tname}.treefile"
        ogp = ALN_DIR / "CL0013" / og_fn
        if tp.is_file() and ogp.is_file():
            nw = tp.read_text(encoding="utf-8")
            if verify_tree_bootstrap_complete(tp, nw):
                chosen = (tname, nw, read_outgroup_ids(ogp))
                break
    if chosen is None:
        return {}

    tname, nw, og_ids = chosen
    t_sup, _ = root_sister_pre(nw, og_ids)
    fam_tips: dict[str, set[str]] = defaultdict(set)
    for lf in t_sup.leaves():
        f = taxa_fam.get(lf.name)
        if f:
            fam_tips[f].add(lf.name)

    crown_map: dict[str, dict] = {}
    for fam in ("PF00144", "PF00905", "PF13354", "PF00768"):
        n_f = len(fam_tips[fam])
        if n_f < 4:
            continue
        best_f1 = -1.0
        best_nd = None
        for nd in t_sup.traverse("postorder"):
            if nd.is_leaf:
                nd.add_prop("_cf", 1 if taxa_fam.get(nd.name) == fam else 0)
                nd.add_prop("_csz", 1)
            else:
                nd.add_prop("_cf", sum(c.props["_cf"] for c in nd.children))
                nd.add_prop("_csz", sum(c.props["_csz"] for c in nd.children))
                fc = nd.props["_cf"]
                fsz = nd.props["_csz"]
                f1 = 2.0 * fc / (fsz + n_f)
                if f1 > best_f1:
                    best_f1 = f1
                    best_nd = nd
        if best_nd is None:
            continue

        u = best_nd
        while not u.is_leaf and len(u.children) >= 2:
            ranked_ch = sorted(u.children, key=lambda c: c.props["_cf"], reverse=True)
            dom_c = ranked_ch[0]
            rest_cf = sum(c.props["_cf"] for c in ranked_ch[1:])
            if rest_cf >= 2 or dom_c.props["_cf"] < max(4, int(0.80 * u.props["_cf"])):
                break
            u = dom_c

        ranked_ch = sorted(u.children, key=lambda c: c.props["_cf"], reverse=True)
        if len(ranked_ch) < 2:
            continue
        s1 = {lf.name for lf in ranked_ch[0].leaves() if taxa_fam.get(lf.name) == fam}
        s2 = {
            lf.name
            for ch in ranked_ch[1:]
            for lf in ch.leaves()
            if taxa_fam.get(lf.name) == fam
        }
        if s1 and s2:
            crown_map[fam] = {
                "s1": s1,
                "s2": s2,
                "cl0013_source": tname,
                "crown_clade_f1": round(best_f1, 4),
                "crown_ufboot": float(u.props.get("ufboot", 0.0)),
            }

    if fam_tips["PF13354"] and fam_tips["PF00768"]:
        crown_map["PF13354_PF00768"] = {
            "s1": set(fam_tips["PF13354"]),
            "s2": set(fam_tips["PF00768"]),
            "cl0013_source": tname,
            "crown_clade_f1": 1.0,
            "crown_ufboot": 100.0,
        }
    return crown_map


def root_crown_projection(nw_str: str, crown_info: dict) -> tuple[Tree, dict]:
    """Two-stage rooting: project CL0013 superfamily crown split (S1 vs S2) onto unrooted ingroup tree."""
    t = Tree(nw_str, parser=1)
    leaf_names = {l.name for l in t.leaves()}
    s1 = crown_info["s1"] & leaf_names
    s2 = crown_info["s2"] & leaf_names
    if not s1 or not s2:
        raise RuntimeError(f"Insufficient crown anchor overlap for {crown_info.get('cl0013_source')}: {len(s1)}/{len(s2)}")

    split_sup, bit_map, all_mask = extract_unrooted_split_supports(t)
    for nd in t.traverse("postorder"):
        if nd.is_leaf:
            nd.add_prop("_n1", 1 if nd.name in s1 else 0)
            nd.add_prop("_n2", 1 if nd.name in s2 else 0)
        else:
            nd.add_prop("_n1", sum(c.props["_n1"] for c in nd.children))
            nd.add_prop("_n2", sum(c.props["_n2"] for c in nd.children))

    leaves, x_leaf, D, edges = _compute_patristic_and_dfs(t)
    n = len(leaves)
    with np.errstate(divide="ignore", invalid="ignore"):
        D_inv = np.where(D > 1e-12, 1.0 / D, 0.0)
        D_pos = np.where(D > 1e-12, 1.0, 0.0)
    D_inv2 = D_inv * D_inv
    matrices = [
        D_pos,
        D_inv,
        x_leaf[:, None] * D_inv,
        D_inv2,
        x_leaf[:, None] * D_inv2,
        (x_leaf[:, None] ** 2) * D_inv2,
    ]
    row_csum = [np.concatenate([[0.0], np.cumsum(M.sum(axis=1))]) for M in matrices]
    pref2d = [_prefix_2d(M) for M in matrices]

    def cross_sum(k: int, s: int, t_end: int) -> float:
        return float(row_csum[k][t_end] - row_csum[k][s]) - _rect_sum(pref2d[k], s, t_end)

    tot1, tot2 = len(s1), len(s2)
    tot_s = tot1 + tot2
    cands: list[tuple[float, float, Tree, float]] = []

    for nd, s_idx, t_end, L, y_v in edges:
        n_a = t_end - s_idx
        n_b = n - n_a
        if n_a == 0 or n_b == 0:
            continue
        x1 = nd.props["_n1"]
        x2 = nd.props["_n2"]
        xs = x1 + x2
        y1 = tot1 - x1
        y2 = tot2 - x2
        ys = tot_s - xs

        f1_a1 = (2.0 * x1 / (xs + tot1)) if (xs + tot1) > 0 else 0.0
        f1_a2 = (2.0 * y2 / (ys + tot2)) if (ys + tot2) > 0 else 0.0
        f1_b1 = (2.0 * y1 / (ys + tot1)) if (ys + tot1) > 0 else 0.0
        f1_b2 = (2.0 * x2 / (xs + tot2)) if (xs + tot2) > 0 else 0.0
        macro_f1 = max(0.5 * (f1_a1 + f1_a2), 0.5 * (f1_b1 + f1_b2))

        cnt = cross_sum(0, s_idx, t_end)
        if cnt <= 0.0:
            ad = 0.0
            alpha = 0.5
        else:
            s_inv = cross_sum(1, s_idx, t_end)
            s_xinv = cross_sum(2, s_idx, t_end) - y_v * s_inv
            s_inv2 = cross_sum(3, s_idx, t_end)
            s_xinv2 = cross_sum(4, s_idx, t_end) - y_v * s_inv2
            s_x2inv2 = cross_sum(5, s_idx, t_end) - 2.0 * y_v * cross_sum(4, s_idx, t_end) + (y_v * y_v) * s_inv2
            if L > 1e-12 and s_inv2 > 0.0:
                num = 2.0 * L * (s_inv - 2.0 * s_xinv2)
                den = 4.0 * L * L * s_inv2
                alpha = min(1.0, max(0.0, num / den))
            else:
                alpha = 0.5
                num = 0.0
                den = 0.0
            rss = cnt - 4.0 * s_xinv + 4.0 * s_x2inv2 - 2.0 * alpha * num + (alpha * alpha) * den
            ad = math.sqrt(max(0.0, rss) / cnt)
        cands.append((macro_f1, -ad, nd, alpha))

    best_macro = max(z[0] for z in cands)
    best_ad_neg, best_nd, best_alpha = max((z[1], z[2], z[3]) for z in cands if z[0] >= best_macro - 1e-6)
    best_ad = -best_ad_neg

    orig_L = max(0.0, float(best_nd.dist or 0.0))
    t.set_outgroup(best_nd)
    for c in t.children:
        if c is best_nd:
            c.dist = best_alpha * orig_L
        else:
            c.dist = (1.0 - best_alpha) * orig_L
    restore_split_supports(t, split_sup, bit_map, all_mask)
    return t, {"root_f1_outgroup": round(best_macro, 4), "outgroup_n": tot_s, "mad_score": round(best_ad, 5)}


# ---------------------------------------------------------------------------
# Polytomy Collapsing, Fitch Parsimony & F81 Marginal ML ASR (PastML Model)
# ---------------------------------------------------------------------------
def collapse_unresolved_polytomies(tree: Tree, eps: float = 1e-6) -> int:
    to_del = []
    for nd in tree.traverse("postorder"):
        if not nd.is_leaf and nd.up is not None:
            ub = float(nd.props.get("ufboot", 0.0))
            dist = float(nd.dist) if nd.dist is not None else 0.0
            if dist <= eps and ub <= 1.0:
                to_del.append(nd)
    for nd in to_del:
        nd.delete(prevent_nondicotomic=False)
    return len(to_del)


def run_fitch_asr(tree: Tree, v_meta: dict[str, dict]) -> None:
    """Two-pass Fitch parsimony for binary state {Bacterial, Viral}."""
    for nd in tree.traverse("postorder"):
        if nd.is_leaf:
            d = tip_domain(nd.name, v_meta)
            nd.add_prop("fitch_set", {d})
            nd.add_prop("n_vir", 1 if d == "Viral" else 0)
            nd.add_prop("n_bac", 0 if d == "Viral" else 1)
        else:
            nd.add_prop("n_vir", sum(c.props["n_vir"] for c in nd.children))
            nd.add_prop("n_bac", sum(c.props["n_bac"] for c in nd.children))
            csets = [c.props["fitch_set"] for c in nd.children]
            inter = set.intersection(*csets)
            nd.add_prop("fitch_set", inter if inter else set.union(*csets))

    rset = tree.props["fitch_set"]
    if len(rset) == 1:
        tree.add_prop("fitch_state", next(iter(rset)))
    else:
        tree.add_prop("fitch_state", "Root_Split")

    for nd in tree.traverse("preorder"):
        if nd is tree:
            continue
        p_state = nd.up.props["fitch_state"]
        fset = nd.props["fitch_set"]
        if p_state == "Root_Split":
            state = next(iter(fset)) if len(fset) == 1 else (
                "Viral" if nd.props["n_vir"] >= nd.props["n_bac"] else "Bacterial"
            )
        elif p_state in fset:
            state = p_state
        elif len(fset) == 1:
            state = next(iter(fset))
        else:
            state = "Viral" if nd.props["n_vir"] >= nd.props["n_bac"] else "Bacterial"
        nd.add_prop("fitch_state", state)


def run_f81_marginal_ml_asr(tree: Tree, v_meta: dict[str, dict]) -> tuple[float, float]:
    """Exact Felsenstein pruning + Marginal Posterior ASR under binary F81 CTMC (PastML model)."""
    post_nodes = list(tree.traverse("postorder"))
    pre_nodes = list(tree.traverse("preorder"))

    n_v = sum(1 for l in tree.leaves() if tip_domain(l.name, v_meta) == "Viral")
    n_b = sum(1 for l in tree.leaves() if tip_domain(l.name, v_meta) != "Viral")
    tot = max(1, n_v + n_b)
    pi1 = min(0.999, max(0.001, n_v / tot))  # Viral
    pi0 = 1.0 - pi1                          # Bacterial
    pi = np.array([pi0, pi1], dtype=np.float64)

    # Index nodes and pre-extract leaf states & branch lengths
    nid = {id(nd): i for i, nd in enumerate(post_nodes)}
    n_nodes = len(post_nodes)
    children_idx = [[nid[id(c)] for c in nd.children] for nd in post_nodes]
    dists = np.array([max(0.0, float(nd.dist or 0.0)) for nd in post_nodes], dtype=np.float64)
    is_leaf = np.array([nd.is_leaf for nd in post_nodes], dtype=bool)
    leaf_state = np.array(
        [1 if (nd.is_leaf and tip_domain(nd.name, v_meta) == "Viral") else 0 for nd in post_nodes],
        dtype=np.int32,
    )

    def eval_upward(mu: float) -> tuple[float, np.ndarray, np.ndarray]:
        L = np.zeros((n_nodes, 2), dtype=np.float64)
        C = np.ones((n_nodes, 2), dtype=np.float64)
        exp_t = np.exp(-mu * dists)
        logL = 0.0

        for idx in range(n_nodes):
            if is_leaf[idx]:
                L[idx, leaf_state[idx]] = 1.0
            else:
                acc = np.ones(2, dtype=np.float64)
                for c_idx in children_idx[idx]:
                    e = exp_t[c_idx]
                    # P[i, j] = pi[j] + (I[i,j] - pi[j]) * e
                    # C_c[i] = sum_j P[i, j] * L_c[j] = e * L_c[i] + (1 - e) * (pi @ L_c)
                    avg_lc = pi0 * L[c_idx, 0] + pi1 * L[c_idx, 1]
                    c_msg = e * L[c_idx] + (1.0 - e) * avg_lc
                    C[c_idx] = c_msg
                    acc *= c_msg
                mx = acc.max()
                if mx > 0.0:
                    L[idx] = acc / mx
                    logL += math.log(mx)
                else:
                    L[idx] = np.array([0.5, 0.5])
        root_idx = n_nodes - 1
        root_lik = pi0 * L[root_idx, 0] + pi1 * L[root_idx, 1]
        logL += math.log(max(root_lik, 1e-300))
        return logL, L, C

    res = minimize_scalar(
        lambda log_mu: -eval_upward(math.exp(log_mu))[0],
        bounds=(math.log(1e-4), math.log(50.0)),
        method="bounded",
    )
    mu_hat = math.exp(float(res.x))
    best_logL, L, C = eval_upward(mu_hat)

    # Downward pass (division-free sibling products)
    D_msg = np.zeros((n_nodes, 2), dtype=np.float64)
    M = np.zeros((n_nodes, 2), dtype=np.float64)
    root_idx = n_nodes - 1
    D_msg[root_idx] = pi
    M[root_idx] = D_msg[root_idx] * L[root_idx]
    s = M[root_idx].sum()
    if s > 0.0:
        M[root_idx] /= s

    exp_t = np.exp(-mu_hat * dists)
    for nd in pre_nodes:
        u_idx = nid[id(nd)]
        ch = children_idx[u_idx]
        if not ch:
            continue
        k = len(ch)
        pref = np.ones((k, 2), dtype=np.float64)
        suff = np.ones((k, 2), dtype=np.float64)
        for j in range(1, k):
            pref[j] = pref[j - 1] * C[ch[j - 1]]
        for j in range(k - 2, -1, -1):
            suff[j] = suff[j + 1] * C[ch[j + 1]]

        for pos, c_idx in enumerate(ch):
            E_u_c = D_msg[u_idx] * pref[pos] * suff[pos]
            e = exp_t[c_idx]
            # D_c[j] = sum_i E[i] * P[i, j] = e * E[j] + (1 - e) * pi[j] * sum(E)
            sum_e = E_u_c[0] + E_u_c[1]
            d_c = e * E_u_c + (1.0 - e) * pi * sum_e
            mx_d = d_c.max()
            if mx_d > 0.0:
                d_c /= mx_d
            D_msg[c_idx] = d_c
            m_c = d_c * L[c_idx]
            ms = m_c.sum()
            if ms > 0.0:
                m_c /= ms
            else:
                m_c = np.array([0.5, 0.5])
            M[c_idx] = m_c

    for idx, nd in enumerate(post_nodes):
        p_bac, p_vir = float(M[idx, 0]), float(M[idx, 1])
        ml_state = "Viral" if p_vir >= 0.50 else "Bacterial"
        nd.add_prop("ml_state", ml_state)
        nd.add_prop("ml_p_vir", p_vir)
        nd.add_prop("ml_p_bac", p_bac)

    return mu_hat, best_logL


# ---------------------------------------------------------------------------
# Irwin et al. / Pre-Locked HGT Detection & Classification
# ---------------------------------------------------------------------------
def detect_hgt_events(
    tree: Tree,
    family: str,
    tree_label: str,
    rooting_method: str,
    v_meta: dict[str, dict],
    b_meta: dict[str, dict],
) -> tuple[list[dict], dict]:
    n_collapsed = collapse_unresolved_polytomies(tree)
    run_fitch_asr(tree, v_meta)
    mu_hat, ml_logL = run_f81_marginal_ml_asr(tree, v_meta)

    tot_tips = len(list(tree.leaves()))
    tot_vir = tree.props["n_vir"]
    tot_bac = tree.props["n_bac"]

    events: list[dict] = []
    ev_idx = 0

    for nd in tree.traverse("preorder"):
        if nd is tree or nd.up is None:
            continue
        p = nd.up
        p_fitch = p.props["fitch_state"]
        c_fitch = nd.props["fitch_state"]
        p_ml = p.props["ml_state"]
        c_ml = nd.props["ml_state"]

        # Candidate transition if Fitch OR Marginal ML shows a domain change (skip Root_Split backbone)
        fitch_trans = (p_fitch != "Root_Split") and (c_fitch != p_fitch)
        ml_trans = (c_ml != p_ml)
        if not (fitch_trans or ml_trans):
            continue

        if fitch_trans:
            from_dom, to_dom = p_fitch, c_fitch
        else:
            from_dom, to_dom = p_ml, c_ml

        fitch_ml_concordant = (p_fitch == from_dom == p_ml) and (c_fitch == to_dom == c_ml)
        ml_parent_prob = float(p.props["ml_p_vir"] if from_dom == "Viral" else p.props["ml_p_bac"])
        ml_focal_prob = float(nd.props["ml_p_vir"] if to_dom == "Viral" else nd.props["ml_p_bac"])

        n_v = nd.props["n_vir"]
        n_b = nd.props["n_bac"]
        n_tot = n_v + n_b
        n_recip = n_v if to_dom == "Viral" else n_b
        purity = n_recip / max(n_tot, 1)

        gp = p.up
        is_root_child = (p is tree)

        s1_v = s1_b = 0
        s1_vir_tips: set[str] = set()
        for sib in p.children:
            if sib is not nd:
                s1_v += sib.props["n_vir"]
                s1_b += sib.props["n_bac"]
                for lf in sib.leaves():
                    if tip_domain(lf.name, v_meta) == "Viral":
                        s1_vir_tips.add(lf.name)
        s1_tot = s1_v + s1_b
        s1_donor = s1_v if from_dom == "Viral" else s1_b

        s2_v = s2_b = 0
        s2_vir_tips: set[str] = set()
        if gp is not None:
            for aunt in gp.children:
                if aunt is not p:
                    s2_v += aunt.props["n_vir"]
                    s2_b += aunt.props["n_bac"]
                    for lf in aunt.leaves():
                        if tip_domain(lf.name, v_meta) == "Viral":
                            s2_vir_tips.add(lf.name)
        s2_tot = s2_v + s2_b
        s2_donor = s2_v if from_dom == "Viral" else s2_b

        bg_tot = s1_tot + s2_tot
        bg_donor = s1_donor + s2_donor
        bg_donor_frac = bg_donor / max(bg_tot, 1)
        tot_recip_tree = tot_vir if to_dom == "Viral" else tot_bac

        passes_nesting = (
            (not is_root_child)
            and (n_tot < 0.25 * tot_tips)
            and (n_recip < 0.50 * max(tot_recip_tree, 1))
            and (purity >= 0.80)
            and (s1_donor >= 1 and s2_donor >= 1)
            and (bg_donor_frac >= 0.75)
        )

        f_sh = float(nd.props.get("sh_alrt", 100.0 if nd.is_leaf else 0.0))
        f_uf = float(nd.props.get("ufboot", 100.0 if nd.is_leaf else 0.0))
        p_sh = float(p.props.get("sh_alrt", 0.0)) if not p.is_leaf else 0.0
        p_uf = float(p.props.get("ufboot", 0.0)) if not p.is_leaf else 0.0
        g_sh = float(gp.props.get("sh_alrt", 0.0)) if (gp is not None and not gp.is_leaf) else 0.0
        g_uf = float(gp.props.get("ufboot", 0.0)) if (gp is not None and not gp.is_leaf) else 0.0

        min_pg_uf = min(p_uf, g_uf) if gp is not None else p_uf
        min_pg_sh = min(p_sh, g_sh) if gp is not None else p_sh
        max_pg_uf = max(p_uf, g_uf)

        recip_leaves = [lf.name for lf in nd.leaves() if tip_domain(lf.name, v_meta) == to_dom]
        is_singleton = (n_tot == 1) or (n_recip == 1)

        # Singleton Viral Tier-1 gating & Prophage_Borne & Functional Stratification
        n_hydro_sbl = 0
        n_pbp_pre = 0
        n_card_40_70 = 0
        expanded_v_cluster = 0
        n_v_tier1 = 0
        n_v_tier2 = 0
        n_v_tier3 = 0
        n_v_mosaic = 0
        n_bac_prophage = 0
        card_fams: list[str] = []

        for uid in recip_leaves:
            if to_dom == "Viral":
                info = v_meta.get(uid, {})
                ft = info.get("functional_tier", "PBP_or_PRE_Homolog")
                vt = info.get("viral_origin_tier", "")
                csz = int(info.get("viral_cluster_size", 1))
                expanded_v_cluster += csz
                if vt == "Tier1_HighConfidence_Viral":
                    n_v_tier1 += 1
                elif vt == "Tier2_Supported_Viral":
                    n_v_tier2 += 1
                elif "MOSAiC" in vt:
                    n_v_mosaic += 1
                else:
                    n_v_tier3 += 1
            else:
                info = b_meta.get(uid, {})
                ft = info.get("functional_tier", "PBP_or_PRE_Homolog")
                if info.get("in_bacterial_prophage", False):
                    n_bac_prophage += 1

            if ft == "Putative_Hydrolytic_BetaLactamase":
                n_hydro_sbl += 1
            else:
                n_pbp_pre += 1
            if ft == "Putative_Hydrolytic_BetaLactamase" and info.get("card_pident", 0.0) >= 40.0 and info.get("card_cov", 0.0) >= 70.0:
                n_card_40_70 += 1
                cf = info.get("card_gene_family", "NONE")
                if cf and cf != "NONE":
                    card_fams.append(cf.split(";")[0].strip())

        if to_dom == "Viral" and is_singleton:
            singleton_viral_pass = (n_v_tier1 == n_recip and n_recip > 0)
        else:
            singleton_viral_pass = True

        prophage_frac = (n_bac_prophage / max(n_recip, 1)) if to_dom == "Bacterial" else 0.0
        is_prophage_borne = (to_dom == "Bacterial") and (n_recip > 0) and (n_bac_prophage == n_recip)

        is_tier1 = (
            passes_nesting
            and fitch_ml_concordant
            and (min_pg_uf >= 95.0)
            and (min_pg_sh >= 80.0)
            and singleton_viral_pass
            and (not is_prophage_borne)
        )
        is_tier2 = (
            passes_nesting
            and fitch_ml_concordant
            and (max_pg_uf >= 80.0)
            and singleton_viral_pass
            and (not is_prophage_borne)
            and (not is_tier1)
        )

        if is_prophage_borne and passes_nesting and fitch_ml_concordant:
            ev_class = "Prophage_Borne"
        elif is_tier1:
            ev_class = "Tier1_HGT"
        elif is_tier2:
            ev_class = "Tier2_Supported_HGT"
        elif passes_nesting and fitch_ml_concordant:
            ev_class = "Nested_LowSupport_Or_Filtered"
        else:
            ev_class = "Unnested_Or_Discordant"

        func_category = (
            "Putative_Hydrolytic_BetaLactamase"
            if n_hydro_sbl >= n_pbp_pre and n_hydro_sbl > 0
            else "PBP_or_PRE_Homolog"
        )

        ev_idx += 1
        events.append({
            "family": family,
            "tree_label": tree_label,
            "rooting_method": rooting_method,
            "event_id": f"{family}__{tree_label}__{rooting_method}__HGT{ev_idx:03d}",
            "direction": f"{from_dom}_to_{to_dom}",
            "donor_domain": from_dom,
            "recipient_domain": to_dom,
            "event_class": ev_class,
            "is_tier1_hgt": is_tier1,
            "is_tier2_hgt": is_tier2,
            "is_prophage_borne": is_prophage_borne,
            "passes_double_outgroup_nesting": passes_nesting,
            "fitch_transition": fitch_trans,
            "ml_transition": ml_trans,
            "fitch_ml_concordant": fitch_ml_concordant,
            "ml_parent_prob": round(ml_parent_prob, 4),
            "ml_focal_prob": round(ml_focal_prob, 4),
            "n_focal_tips_total": n_tot,
            "n_focal_recipient_tips": n_recip,
            "focal_recipient_purity": round(purity, 4),
            "focal_is_singleton_tip": is_singleton,
            "singleton_viral_tier1_pass": singleton_viral_pass,
            "s1_donor_n": s1_donor,
            "s1_total_n": s1_tot,
            "s2_donor_n": s2_donor,
            "s2_total_n": s2_tot,
            "combined_outgroup_donor_frac": round(bg_donor_frac, 4),
            "focal_ufboot": f_uf,
            "focal_shalrt": f_sh,
            "parent_ufboot": p_uf,
            "parent_shalrt": p_sh,
            "grandparent_ufboot": g_uf,
            "grandparent_shalrt": g_sh,
            "min_PG_ufboot": min_pg_uf,
            "min_PG_shalrt": min_pg_sh,
            "max_PG_ufboot": max_pg_uf,
            "functional_category": func_category,
            "n_hydrolytic_sbl_tips": n_hydro_sbl,
            "n_pbp_pre_tips": n_pbp_pre,
            "n_card_40_70_tips": n_card_40_70,
            "card_families": ";".join(sorted(set(card_fams))[:5]),
            "expanded_viral_cluster_seqs": expanded_v_cluster if to_dom == "Viral" else 0,
            "n_viral_tier1": n_v_tier1,
            "n_viral_tier2": n_v_tier2,
            "n_viral_tier3": n_v_tier3,
            "n_viral_mosaic": n_v_mosaic,
            "prophage_recipient_fraction": round(prophage_frac, 4),
            "representative_tips": ";".join(recip_leaves[:5]),
            "all_recipient_tips": ";".join(recip_leaves),
            "_recip_set": set(recip_leaves),
            "_donor_vir_sister_set": s1_vir_tips | s2_vir_tips,
        })

    stats = {
        "total_tips": tot_tips,
        "viral_tips": tot_vir,
        "bacterial_tips": tot_bac,
        "collapsed_polytomies": n_collapsed,
        "f81_mu_hat": round(mu_hat, 4),
        "f81_logL": round(ml_logL, 2),
    }
    return events, stats


# ---------------------------------------------------------------------------
# OAT Matching & CL0013 Validation
# ---------------------------------------------------------------------------
def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def match_event_in_comparator(
    prim_ev: dict,
    prim_tips: set[str],
    comp_events: list[dict],
    comp_tips: set[str],
    allow_donor_sister_fallback: bool = False,
) -> tuple[bool, str]:
    p_dir = prim_ev["direction"]
    shared_tips = prim_tips & comp_tips
    p_recip = prim_ev["_recip_set"] & shared_tips
    best_match = None

    for cev in comp_events:
        if not (cev["passes_double_outgroup_nesting"] and cev["fitch_ml_concordant"]):
            continue
        if cev["direction"] != p_dir:
            continue

        c_recip = cev["_recip_set"] & shared_tips
        if len(p_recip) > 0 and len(c_recip) > 0:
            if jaccard(p_recip, c_recip) >= 0.70:
                best_match = cev
                break
        elif allow_donor_sister_fallback and p_dir == "Viral_to_Bacterial":
            p_sis = prim_ev["_donor_vir_sister_set"] & shared_tips
            c_sis = cev["_donor_vir_sister_set"] & shared_tips
            if len(p_sis) > 0 and jaccard(p_sis, c_sis) >= 0.50:
                best_match = cev
                break

    if best_match is None:
        return (False, "Unmatched")
    return (True, best_match["event_class"])


def validate_cl0013_superfamily() -> list[dict]:
    meta_files = [
        ALN_DIR / "CL0013" / "cl0013_taxa_metadata.tsv",
        ALN_DIR / "CL0013" / "cl0013_struct_taxa_metadata.tsv",
    ]
    taxa_fam: dict[str, str] = {}
    taxa_role: dict[str, str] = {}
    for mp in meta_files:
        if mp.is_file():
            with open(mp, encoding="utf-8") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    taxa_fam[r["uid"]] = r["cl0013_family"]
                    taxa_role[r["uid"]] = r["role"]
    if not taxa_fam:
        return []

    rows: list[dict] = []
    for tree_name, og_fn in CL0013_TREES:
        tpath = TREES_DIR / "CL0013" / f"{tree_name}.treefile"
        og_path = ALN_DIR / "CL0013" / og_fn
        if not tpath.is_file() or not og_path.is_file():
            continue
        nw_str = tpath.read_text(encoding="utf-8")
        if not verify_tree_bootstrap_complete(tpath, nw_str):
            continue
        og_ids = read_outgroup_ids(og_path)
        t = Tree(nw_str, parser=1)
        all_leaves = {l.name for l in t.leaves()}
        og_present = og_ids & all_leaves
        in_leaves = all_leaves - og_present
        split_sup, bit_map, all_mask = extract_unrooted_split_supports(t, leaf_subset=in_leaves)

        # Find best non-PRE outgroup split and record F1
        n_og = len(og_present)
        n_tot = len(all_leaves)
        for nd in t.traverse("postorder"):
            if nd.is_leaf:
                nd.add_prop("_og", 1 if nd.name in og_present else 0)
                nd.add_prop("_sz", 1)
            else:
                nd.add_prop("_og", sum(c.props["_og"] for c in nd.children))
                nd.add_prop("_sz", sum(c.props["_sz"] for c in nd.children))

        best_nd = None
        best_f1 = -1.0
        for nd in t.traverse():
            if nd.up is None:
                continue
            k_og = nd.props["_og"]
            k_sz = nd.props["_sz"]
            if k_og > 0:
                f1 = 2.0 * k_og / (k_sz + n_og)
                if f1 > best_f1:
                    best_f1 = f1
                    best_nd = nd
            comp_og = n_og - k_og
            comp_sz = n_tot - k_sz
            if comp_og > 0:
                f1 = 2.0 * comp_og / (comp_sz + n_og)
                if f1 > best_f1:
                    best_f1 = f1
                    best_nd = nd
        if best_nd is not None:
            t.set_outgroup(best_nd)
        t.prune(sorted(in_leaves), preserve_branch_length=True)
        restore_split_supports(t, split_sup, bit_map, all_mask)

        # Evaluate best clade F1 for each of the 5 PRE families inside the rooted PRE ingroup
        pre_fams = ("PF00144", "PF00905", "PF13354", "PF00768", "PF02113")
        fam_stats: dict[str, tuple[int, float, float, float]] = {}
        fam_counts = defaultdict(int)
        for l in t.leaves():
            fam_counts[taxa_fam.get(l.name, "Unknown")] += 1

        for fam in pre_fams:
            n_f = fam_counts.get(fam, 0)
            if n_f == 0:
                fam_stats[fam] = (0, 0.0, 0.0, 0.0)
                continue
            best_cl_f1 = 0.0
            best_sh = 0.0
            best_uf = 0.0
            for nd in t.traverse("postorder"):
                if nd.is_leaf:
                    nd.add_prop("_fc", 1 if taxa_fam.get(nd.name) == fam else 0)
                    nd.add_prop("_fsz", 1)
                else:
                    nd.add_prop("_fc", sum(c.props["_fc"] for c in nd.children))
                    nd.add_prop("_fsz", sum(c.props["_fsz"] for c in nd.children))
                fc = nd.props["_fc"]
                fsz = nd.props["_fsz"]
                f1 = 2.0 * fc / (fsz + n_f)
                if f1 > best_cl_f1:
                    best_cl_f1 = f1
                    best_sh = float(nd.props.get("sh_alrt", 100.0 if nd.is_leaf else 0.0))
                    best_uf = float(nd.props.get("ufboot", 100.0 if nd.is_leaf else 0.0))
            fam_stats[fam] = (n_f, round(best_cl_f1, 4), best_sh, best_uf)

        row = {
            "cl0013_tree": tree_name,
            "n_ingroup_pre": len(in_leaves),
            "n_outgroup_nonpre": len(og_present),
            "nonpre_outgroup_split_F1": round(best_f1, 4),
        }
        for fam in pre_fams:
            nf, bf1, bsh, buf = fam_stats[fam]
            row[f"{fam}_n"] = nf
            row[f"{fam}_clade_F1"] = bf1
            row[f"{fam}_SH_aLRT"] = bsh
            row[f"{fam}_UFBoot"] = buf
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Main Pipeline
# ---------------------------------------------------------------------------
def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    v_meta, b_meta = load_catalog_metadata()
    crown_map = extract_superfamily_crown_splits()
    if crown_map:
        print("[CL0013 Crown Anchors] Projected from superfamily backbone:", flush=True)
        for cfam, cinfo in crown_map.items():
            print(
                f"  {cfam:15s}: S1={len(cinfo['s1'])}, S2={len(cinfo['s2'])} "
                f"(source={cinfo['cl0013_source']}, F1={cinfo['crown_clade_f1']}, UFBoot={cinfo['crown_ufboot']:.0f})",
                flush=True,
            )

    all_events_rows: list[dict] = []
    tree_matrix_rows: list[dict] = []
    primary_catalog_rows: list[dict] = []

    # Cache events & leaf sets per (family, comp_label) for OAT matching
    events_by_fam: dict[str, dict[str, tuple[list[dict], set[str]]]] = {}

    for fam in FAMILIES:
        primary_track = "mafft" if fam == "PF13354_PF00768" else "hmm"
        og_path = ALN_DIR / fam / "outgroup_taxa.txt"
        og_ids = read_outgroup_ids(og_path) if og_path.is_file() else set()
        events_by_fam[fam] = {}

        # Define all candidate trees for this family: (oat_key, rooting_method, tree_label, tree_path, og_id_set)
        tree_jobs: list[tuple[str, str, str, Path, set[str]]] = []
        prim_out_path = TREES_DIR / fam / f"{fam}.primary.{primary_track}_out.treefile"
        prim_in_path = TREES_DIR / fam / f"{fam}.primary.{primary_track}_in.treefile"

        if not prim_in_path.is_file() and not prim_out_path.is_file():
            print(f"[SKIP] Missing primary trees for {fam}: {prim_in_path} / {prim_out_path}")
            continue

        if prim_in_path.is_file() and fam in crown_map:
            tree_jobs.append(("primary_out", "Crown_Projected", "Primary_Crown_Projected", prim_in_path, set()))
            if prim_out_path.is_file():
                tree_jobs.append(("legacy_sarg24_out", "Sister_PRE", "OAT_Legacy_SARG24_Outgroup", prim_out_path, og_ids))
        elif prim_out_path.is_file():
            tree_jobs.append(("primary_out", "Sister_PRE", "Primary_Sister_PRE", prim_out_path, og_ids))

        clean_og_p = ALN_DIR / fam / "outgroup_clean_gtdb12.txt"
        clean_tp = TREES_DIR / fam / f"{fam}.clean_prof_out.{primary_track}_out.treefile"
        if clean_og_p.is_file() and clean_tp.is_file():
            tree_jobs.append(("clean_prof_out", "Sister_PRE", "OAT_Clean_GTDB12_Profile", clean_tp, read_outgroup_ids(clean_og_p)))

        if prim_in_path.is_file():
            tree_jobs.append(("primary_in_MAD", "MAD", "Primary_Ingroup_MAD", prim_in_path, set()))
            tree_jobs.append(("primary_in_MinVar", "MinVar", "Primary_Ingroup_MinVar", prim_in_path, set()))
            tree_jobs.append(("primary_in_Midpoint", "Midpoint", "Primary_Ingroup_Midpoint", prim_in_path, set()))
            tree_jobs.append(("primary_in_YoudenJ", "Youden_J_Appendix", "Appendix_Youden_J", prim_in_path, set()))

        if fam != "PF13354_PF00768":
            mafft_out = TREES_DIR / fam / f"{fam}.primary.mafft_out.treefile"
            if mafft_out.is_file():
                tree_jobs.append(("mafft_out", "Sister_PRE", "OAT_MAFFT_ClipKIT", mafft_out, og_ids))

        for sname in ("rep2", "rep3", "rep4", "rep5", "phylo", "nopro", "cov60", "cov70"):
            tp = TREES_DIR / fam / f"{fam}.{sname}.{primary_track}_out.treefile"
            if tp.is_file():
                tree_jobs.append((sname, "Sister_PRE", f"OAT_{sname}", tp, og_ids))

        pmsf_tp = TREES_DIR / fam / f"{fam}.pmsf_c20.{primary_track}_out.treefile"
        if pmsf_tp.is_file():
            tree_jobs.append(("pmsf_c20", "Sister_PRE", "OAT_PMSF_C20", pmsf_tp, og_ids))

        for key, rmethod, label, tpath, job_og_ids in tree_jobs:
            nw_str = tpath.read_text(encoding="utf-8")
            if not verify_tree_bootstrap_complete(tpath, nw_str):
                continue

            if rmethod == "Crown_Projected":
                tr, rmeta = root_crown_projection(nw_str, crown_map[fam])
            elif rmethod == "Sister_PRE":
                tr, rmeta = root_sister_pre(nw_str, job_og_ids)
            elif rmethod == "MAD":
                tr, rmeta = root_mad(nw_str)
            elif rmethod in ("MinVar", "Midpoint"):
                tr, rmeta = root_minvar_or_midpoint(nw_str, rmethod)
            else:
                tr, rmeta = root_youden_j(nw_str, v_meta)

            if key == "primary_out":
                write_rooted_treefile(tr, TREES_DIR / fam / f"{fam}.primary.rooted.treefile")

            ev_list, stats = detect_hgt_events(tr, fam, label, rmethod, v_meta, b_meta)
            tip_set = {l.name for l in tr.leaves()}
            events_by_fam[fam][key] = (ev_list, tip_set)

            for ev in ev_list:
                row = {k: v for k, v in ev.items() if not k.startswith("_")}
                all_events_rows.append(row)

            # Summary row for OAT Sensitivity Matrix
            nested_all = [e for e in ev_list if e["passes_double_outgroup_nesting"] and e["fitch_ml_concordant"]]
            t1_all = [e for e in nested_all if e["is_tier1_hgt"]]
            t2_all = [e for e in nested_all if e["is_tier2_hgt"]]
            pb_all = [e for e in nested_all if e["is_prophage_borne"]]
            t1_b2v = [e for e in t1_all if e["direction"] == "Bacterial_to_Viral"]
            t1_v2b = [e for e in t1_all if e["direction"] == "Viral_to_Bacterial"]
            t1_sbl = [e for e in t1_all if e["functional_category"] == "Putative_Hydrolytic_BetaLactamase"]

            mat_row = {
                "family": fam,
                "oat_key": key,
                "tree_label": label,
                "rooting_method": rmethod,
                **stats,
                **rmeta,
                "nested_fitch_ml_events": len(nested_all),
                "tier1_hgt_events": len(t1_all),
                "tier1_B_to_V_events": len(t1_b2v),
                "tier1_V_to_B_events": len(t1_v2b),
                "tier1_hydrolytic_sbl_events": len(t1_sbl),
                "tier2_hgt_events": len(t2_all),
                "prophage_borne_events": len(pb_all),
                "tier1_B_to_V_viral_nr85_tips": sum(e["n_focal_recipient_tips"] for e in t1_b2v),
                "tier1_B_to_V_expanded_seqs": sum(e["expanded_viral_cluster_seqs"] for e in t1_b2v),
                "tier1_V_to_B_bac_tips": sum(e["n_focal_recipient_tips"] for e in t1_v2b),
            }
            tree_matrix_rows.append(mat_row)
            print(
                f"  [{fam:15s} | {key:20s}] tips={stats['total_tips']:4d} "
                f"nested={len(nested_all):3d} Tier1={len(t1_all):2d} "
                f"(B->V:{len(t1_b2v):2d}, V->B:{len(t1_v2b):2d}, SBL:{len(t1_sbl):2d}) "
                f"Tier2={len(t2_all):2d} Prophage={len(pb_all):2d}",
                flush=True,
            )

    # Cross-tree OAT Robustness Scoring on Primary Nested Concordant Events
    # Also include Joint PF13354_PF00768 as an OAT comparator for PF13354 & PF00768
    oat_order = [
        "primary_in_MAD", "primary_in_MinVar", "primary_in_Midpoint", "clean_prof_out",
        "mafft_out", "rep2", "rep3", "rep4", "rep5", "phylo", "nopro", "cov60", "cov70", "pmsf_c20",
    ]
    resample_keys = {"rep2", "rep3", "rep4", "rep5", "phylo"}

    for fam in FAMILIES:
        if fam not in events_by_fam or "primary_out" not in events_by_fam[fam]:
            continue
        prim_evs, prim_tips = events_by_fam[fam]["primary_out"]
        prim_nested = [e for e in prim_evs if e["passes_double_outgroup_nesting"] and e["fitch_ml_concordant"]]

        for pev in prim_nested:
            rec = {k: v for k, v in pev.items() if not k.startswith("_")}
            # Cross-rooting checks
            mad_ok, mad_cls = (False, "NA")
            minvar_ok, minvar_cls = (False, "NA")
            mid_ok, mid_cls = (False, "NA")
            youd_ok, youd_cls = (False, "NA")

            if "primary_in_MAD" in events_by_fam[fam]:
                cevs, ctips = events_by_fam[fam]["primary_in_MAD"]
                mad_ok, mad_cls = match_event_in_comparator(pev, prim_tips, cevs, ctips)
            if "primary_in_MinVar" in events_by_fam[fam]:
                cevs, ctips = events_by_fam[fam]["primary_in_MinVar"]
                minvar_ok, minvar_cls = match_event_in_comparator(pev, prim_tips, cevs, ctips)
            if "primary_in_Midpoint" in events_by_fam[fam]:
                cevs, ctips = events_by_fam[fam]["primary_in_Midpoint"]
                mid_ok, mid_cls = match_event_in_comparator(pev, prim_tips, cevs, ctips)
            if "primary_in_YoudenJ" in events_by_fam[fam]:
                cevs, ctips = events_by_fam[fam]["primary_in_YoudenJ"]
                youd_ok, youd_cls = match_event_in_comparator(pev, prim_tips, cevs, ctips)
            leg_ok = False
            if "legacy_sarg24_out" in events_by_fam[fam]:
                cevs, ctips = events_by_fam[fam]["legacy_sarg24_out"]
                leg_ok, _ = match_event_in_comparator(pev, prim_tips, cevs, ctips)

            root_invariant = mad_ok and (minvar_ok or mid_ok)
            rec["root_MAD_match"] = mad_ok
            rec["root_MinVar_match"] = minvar_ok
            rec["root_Midpoint_match"] = mid_ok
            rec["root_Invariant"] = root_invariant
            rec["appendix_YoudenJ_match"] = youd_ok
            rec["legacy_sarg24_match"] = leg_ok

            # OAT recurrence across completed comparators (excluding Midpoint to avoid double-counting clock family,
            # while including MAD + MinVar + MSA + sampling + prophage + cov + PMSF + joint tree)
            oat_eval_keys = [
                k for k in oat_order
                if k != "primary_in_Midpoint" and k in events_by_fam[fam]
            ]
            n_avail = 0
            n_matched = 0
            matched_labels: list[str] = []
            for k in oat_eval_keys:
                # Skip nopro for strictly Prophage_Borne events where all recipient tips were deleted
                if k == "nopro" and pev["is_prophage_borne"]:
                    continue
                cevs, ctips = events_by_fam[fam][k]
                # Skip cov60/cov70 if recipient tip wasn't included in cov70
                if k in ("cov60", "cov70") and len(pev["_recip_set"] & ctips) == 0:
                    continue
                n_avail += 1
                ok, _ = match_event_in_comparator(
                    pev, prim_tips, cevs, ctips, allow_donor_sister_fallback=(k in resample_keys)
                )
                rec[f"oat_{k}"] = ok
                if ok:
                    n_matched += 1
                    matched_labels.append(k)

            # Joint PF13354_PF00768 comparator for PF13354 and PF00768
            if fam in ("PF13354", "PF00768") and "primary_out" in events_by_fam.get("PF13354_PF00768", {}):
                j_evs, j_tips = events_by_fam["PF13354_PF00768"]["primary_out"]
                if len(pev["_recip_set"] & j_tips) > 0:
                    n_avail += 1
                    ok, _ = match_event_in_comparator(pev, prim_tips, j_evs, j_tips, allow_donor_sister_fallback=True)
                    rec["oat_joint_tree"] = ok
                    if ok:
                        n_matched += 1
                        matched_labels.append("joint_tree")

            oat_rate = n_matched / max(n_avail, 1)
            rec["oat_n_matched"] = n_matched
            rec["oat_n_evaluated"] = n_avail
            rec["oat_recurrence_rate"] = round(oat_rate, 4)
            rec["robust_core_hgt"] = (oat_rate >= 0.80) and (pev["is_tier1_hgt"] or pev["is_tier2_hgt"])
            rec["matched_oat_trees"] = ";".join(matched_labels)
            primary_catalog_rows.append(rec)

    # Write outputs
    df_prim = pd.DataFrame(primary_catalog_rows)
    df_all = pd.DataFrame(all_events_rows)
    df_mat = pd.DataFrame(tree_matrix_rows)

    df_prim.to_csv(OUT_DIR / "Primary_HGT_Events_Catalog_2026.tsv", sep="\t", index=False)
    df_all.to_csv(OUT_DIR / "All_Trees_HGT_Events_Long.tsv", sep="\t", index=False)
    df_mat.to_csv(OUT_DIR / "OAT_Sensitivity_Summary_Matrix_2026.tsv", sep="\t", index=False)

    # Family x Functional Tier summary table
    sum_rows = []
    if len(df_prim) > 0:
        for fam in FAMILIES:
            sub_f = df_prim[df_prim["family"] == fam]
            for d in ("Bacterial_to_Viral", "Viral_to_Bacterial"):
                for fc in ("Putative_Hydrolytic_BetaLactamase", "PBP_or_PRE_Homolog"):
                    s = sub_f[(sub_f["direction"] == d) & (sub_f["functional_category"] == fc)]
                    t1 = s[s["is_tier1_hgt"]]
                    t2 = s[s["is_tier2_hgt"]]
                    pb = s[s["is_prophage_borne"]]
                    rc = s[s["robust_core_hgt"]]
                    ri = t1[t1["root_Invariant"]]
                    sum_rows.append({
                        "family": fam,
                        "direction": d,
                        "functional_category": fc,
                        "all_nested_concordant": len(s),
                        "tier1_events": len(t1),
                        "tier1_root_invariant": len(ri),
                        "tier1_robust_core": len(t1[t1["robust_core_hgt"]]),
                        "tier1_recip_nr85_tips": int(t1["n_focal_recipient_tips"].sum()),
                        "tier1_expanded_viral_seqs": int(t1["expanded_viral_cluster_seqs"].sum()),
                        "tier2_events": len(t2),
                        "robust_core_total_t1_t2": len(rc),
                        "prophage_borne_events": len(pb),
                    })
    df_sum = pd.DataFrame(sum_rows)
    df_sum.to_csv(OUT_DIR / "HGT_Summary_By_Family_And_Function_2026.tsv", sep="\t", index=False)

    # CL0013 Superfamily Validation
    cl_rows = validate_cl0013_superfamily()
    if cl_rows:
        df_cl = pd.DataFrame(cl_rows)
        df_cl.to_csv(OUT_DIR / "CL0013_Superfamily_Rooting_Validation.tsv", sep="\t", index=False)
    else:
        df_cl = pd.DataFrame()

    # Human-readable report
    lines = [
        "=== Stage 3 (Phase 4): Pre-Locked Primary HGT Inference & OAT Robustness Summary ===",
        f"Trees analyzed across all families/perturbations: {len(df_mat)}",
        "",
        "--- Primary Crown-Projected Nested Concordant HGT Summary (by Family x Direction x Function) ---",
        df_sum[df_sum["all_nested_concordant"] > 0].to_string(index=False) if len(df_sum) else "(none)",
        "",
        "--- Primary Trees Overview ---",
        df_mat[df_mat["oat_key"] == "primary_out"][[
            "family", "total_tips", "viral_tips", "bacterial_tips",
            "nested_fitch_ml_events", "tier1_hgt_events", "tier1_B_to_V_events",
            "tier1_V_to_B_events", "tier1_hydrolytic_sbl_events",
            "tier1_B_to_V_viral_nr85_tips", "tier1_B_to_V_expanded_seqs",
            "tier2_hgt_events", "prophage_borne_events",
        ]].to_string(index=False),
    ]
    if len(df_cl) > 0:
        lines.extend(["", "--- CL0013 Superfamily & LOO Rooting Validation ---", df_cl.to_string(index=False)])
    else:
        lines.extend(["", "--- CL0013 Superfamily & LOO Rooting Validation ---", "(pending CL0013 IQ-TREE runs)"])

    rep_txt = "\n".join(lines) + "\n"
    (OUT_DIR / "Stage3_HGT_Summary_Report.txt").write_text(rep_txt, encoding="utf-8")
    print("\n" + rep_txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
