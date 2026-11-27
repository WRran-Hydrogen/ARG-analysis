#!/usr/bin/env python3
"""Audit MOSAiC CL0013/SBL candidates via overlapping-envelope Pfam-A (--cut_ga)."""

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

WORK = hpc.BASE / "01_curated_datasets" / "mosaic_audit"
CAND_FAA = WORK / "mosaic_cl0013_candidates.faa"
CL13_DOMTBL = WORK / "mosaic_cl0013.domtbl"
PFAM_ALL_DOMTBL = WORK / "mosaic_pfam_all.domtbl"

PFAM_CLANS_TSV = hpc.PFAM_CLANS_TSV
VIR_INFO = hpc.RAW_INPUT_DIR / "vir.info"
ORIG_TSV = hpc.RAW_INPUT_DIR / "sbl_mosaic.tsv"
ORIG_FAA = hpc.RAW_INPUT_DIR / "mosaic_virus_pf00144.faa"

FOCAL_FAMILIES = ("PF00144", "PF00905", "PF13354")
AUX_PRE_FAMILIES = ("PF00768", "PF02113")
NON_PRE_CL0013 = ("PF04960", "PF06089", "PF27310")
OV_THRESHOLD_AA = 30
DOMAIN_PAD_AA = 10


def read_fasta(path):
    recs = {}
    cur_id, cur_hdr, buf = None, None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if cur_id is not None:
                    recs[cur_id] = (cur_hdr, "".join(buf).rstrip("*"))
                cur_hdr = s[1:].strip()
                cur_id = cur_hdr.split()[0]
                buf = []
            else:
                buf.append(s)
        if cur_id is not None:
            recs[cur_id] = (cur_hdr, "".join(buf).rstrip("*"))
    return recs


def write_fasta(path, items):
    with open(path, "w", encoding="utf-8") as f:
        for hdr, seq in items:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def interval_overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def load_pfam_clans():
    acc_to_clan = {}
    with open(PFAM_CLANS_TSV, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                acc_to_clan[p[0]] = p[1] if p[1] else "NO_CLAN"
    return acc_to_clan


def parse_hmmsearch_domtbl(path):
    hits = defaultdict(list)
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 22:
                continue
            hits[p[0]].append({
                "pfam_name": p[3],
                "pfam_acc": p[4].split(".")[0],
                "tlen": int(p[2]),
                "dom_score": float(p[13]),
                "hmm_from": int(p[15]),
                "hmm_to": int(p[16]),
                "env_from": int(p[19]),
                "env_to": int(p[20]),
            })
    return hits


def resolve_overlapping_domains(hit_list, acc_to_clan, ov_thresh=OV_THRESHOLD_AA):
    n = len(hit_list)
    if n == 0:
        return []
    adj = [[] for _ in range(n)]
    for i in range(n):
        a = hit_list[i]
        for j in range(i + 1, n):
            b = hit_list[j]
            if interval_overlap(a["env_from"], a["env_to"], b["env_from"], b["env_to"]) >= ov_thresh:
                adj[i].append(j)
                adj[j].append(i)
    seen = [False] * n
    winners = []
    for i in range(n):
        if seen[i]:
            continue
        stack = [i]
        seen[i] = True
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
        comp_hits = [hit_list[idx] for idx in comp]
        best = max(comp_hits, key=lambda h: (h["dom_score"], h["env_to"] - h["env_from"]))
        winners.append({
            **best,
            "clan": acc_to_clan.get(best["pfam_acc"], "NO_CLAN"),
            "comp_candidates": ",".join(sorted({h["pfam_acc"] for h in comp_hits})),
        })
    winners.sort(key=lambda h: h["env_from"])
    return winners


def classify_winners(winners):
    if not winners:
        return "drop_below_GA_all_Pfam", "NONE", None
    cl13 = [w for w in winners if w["clan"] == "CL0013"]
    if not cl13:
        top_w = max(winners, key=lambda w: w["dom_score"])
        return "drop_non_CL0013_only", top_w["pfam_acc"], top_w
    best_cl = max(cl13, key=lambda w: w["dom_score"])
    acc = best_cl["pfam_acc"]
    if acc in FOCAL_FAMILIES:
        return f"focal_{acc}", acc, best_cl
    if acc in AUX_PRE_FAMILIES:
        return f"aux_{acc}", acc, best_cl
    if acc in NON_PRE_CL0013:
        return f"drop_non_PRE_CL0013_{acc}", acc, best_cl
    return f"drop_other_CL0013_{acc}", acc, best_cl


def main():
    cands = read_fasta(CAND_FAA)
    acc_to_clan = load_pfam_clans()
    cl13_hits = parse_hmmsearch_domtbl(CL13_DOMTBL)
    pfam_all_hits = parse_hmmsearch_domtbl(PFAM_ALL_DOMTBL)

    orig_tsv_hits = {}
    with open(ORIG_TSV, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                p = line.split()
                orig_tsv_hits[p[0]] = f"{p[3]}({p[5]})"

    orig_fa_ids = set(read_fasta(ORIG_FAA).keys())

    vir_map = {}
    target_set = set(cands.keys())
    with open(VIR_INFO, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            p = line.strip().split(",")
            if len(p) >= 2 and p[1] in target_set:
                vir_map[p[1]] = p[0]

    audit_rows = []
    cat_counts = Counter()
    tax_counts = Counter()
    full_exports = defaultdict(list)
    dom_exports = defaultdict(list)

    for sid in sorted(cands):
        hdr, prot = cands[sid]
        vbin = vir_map.get(sid, "NOT_IN_vir.info")
        # e.g. vir_metabat_out_v1.38885_caudoviricetes or vir_unbinned_unknown
        tax_tag = vbin.rsplit("_", 1)[-1] if "_" in vbin else "unknown"
        bin_type = "unbinned" if "unbinned" in vbin else ("binned_vMAG" if vbin != "NOT_IN_vir.info" else "unknown")

        m_contig = re.match(r"^(.+)_\d+_\d+$", sid)
        contig_id = m_contig.group(1) if m_contig else sid
        partial_flag = "partial" in hdr.lower()
        partial_detail = ""
        mp = re.search(r"(partial=[^\s;]+)", hdr)
        if mp:
            partial_detail = mp.group(1)
        gene_caller = "GeneMark" if "GeneMark" in hdr else ("Prodigal" if "Prodigal" in hdr else "unknown")

        winners = resolve_overlapping_domains(pfam_all_hits.get(sid, []), acc_to_clan, OV_THRESHOLD_AA)
        cat, win_acc, best_dom = classify_winners(winners)
        cat_counts[cat] += 1
        tax_counts[(cat, tax_tag)] += 1

        arch_str = "+".join(f"{w['pfam_acc']}:{w['pfam_name']}({w['env_from']}-{w['env_to']})" for w in winners)
        cl13_ga_str = ",".join(sorted({h["pfam_acc"] for h in cl13_hits.get(sid, [])}))
        dom_start = max(1, best_dom["env_from"] - DOMAIN_PAD_AA) if best_dom else 1
        dom_end = min(len(prot), best_dom["env_to"] + DOMAIN_PAD_AA) if best_dom else len(prot)
        dom_seq = prot[dom_start - 1 : dom_end]
        has_sxxk = bool(re.search(r"S[A-Z]{2}K", dom_seq))

        audit_rows.append({
            "seq_id": sid,
            "inferred_contig": contig_id,
            "vir_bin": vbin,
            "bin_type": bin_type,
            "vir_tax_suffix": tax_tag,
            "gene_caller": gene_caller,
            "is_partial_orf": partial_flag,
            "partial_detail": partial_detail,
            "seq_len": len(prot),
            "in_orig_sbl_mosaic_tsv": sid in orig_tsv_hits,
            "orig_sbl_mosaic_hit": orig_tsv_hits.get(sid, ""),
            "in_orig_mosaic_pf00144_faa": sid in orig_fa_ids,
            "official_cl0013_ga_hits": cl13_ga_str,
            "pfam_category": cat,
            "winning_acc": win_acc,
            "win_dom_score": best_dom["dom_score"] if best_dom else 0.0,
            "env_from": best_dom["env_from"] if best_dom else 0,
            "env_to": best_dom["env_to"] if best_dom else 0,
            "dom_slice_start": dom_start,
            "dom_slice_end": dom_end,
            "has_SXXK_in_dom": has_sxxk,
            "domain_architecture": arch_str,
            "full_header": hdr,
        })

        if cat.startswith("focal_") or cat.startswith("aux_"):
            tag = f"{sid} fam={win_acc} source=MOSAiC bin={vbin} dom={dom_start}..{dom_end} score={best_dom['dom_score']:.1f}"
            full_exports[cat].append((tag, prot))
            dom_exports[cat].append((tag, dom_seq))

    with open(WORK / "MOSAiC_SBL_Audit.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(audit_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(audit_rows)

    for fam in FOCAL_FAMILIES:
        cat = f"focal_{fam}"
        write_fasta(WORK / f"mosaic_{fam}_full.faa", full_exports.get(cat, []))
        write_fasta(WORK / f"mosaic_{fam}_dom.faa", dom_exports.get(cat, []))
    for acc in AUX_PRE_FAMILIES:
        cat = f"aux_{acc}"
        write_fasta(WORK / f"mosaic_CL0013_aux_{acc}_full.faa", full_exports.get(cat, []))
        write_fasta(WORK / f"mosaic_CL0013_aux_{acc}_dom.faa", dom_exports.get(cat, []))

    lines = [
        "=== MOSAiC Viral Proteome CL0013 & Pfam-A (--cut_ga) Audit Summary ===",
        f"Total union candidates checked: {len(cands)}",
        f"In original sbl_mosaic.tsv: {len(orig_tsv_hits)} | In mosaic_virus_pf00144.faa: {len(orig_fa_ids)}",
        "\nPfam-A (--cut_ga) Overlapping-Envelope Competition Results:",
    ]
    for k, v in sorted(cat_counts.items()):
        lines.append(f"  {k}: {v}")
    lines.append("\nBreakdown by (category, vir.info taxonomy suffix):")
    for (k, t), v in sorted(tax_counts.items()):
        lines.append(f"  ({k}, {t}): {v}")
    summary = "\n".join(lines) + "\n"
    (WORK / "MOSAiC_Audit_Summary.txt").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
