#!/usr/bin/env python3
"""Build Final Phase 1 + Phase 1B Curated Viral & Triple-Source Bacterial SBL Datasets.

Rules implemented (user-approved):
1. Viral Decontamination: Scheme B surviving IMG/VR (4,695) + Euk-purged MOSAiC (25).
2. Bacterial Pool (Triple-Source):
   - GTDB/RefSeq 3,000 representative genomes (geNomad prophage-excised, native 1:1 verified CDS)
   - 259 dereplicated marine bacterial MAGs (geNomad prophage-excised, complete ORF, 1:1 verified CDS)
   - CARD/SARG clinical gold-standard bacterial references (CleanCore non-NR subset, plasmid/MGE tagged)
3. Focal-Family Priority Pfam-A (--cut_ga) Overlapping-Envelope Resolution:
   - Any CL0013-winning cluster containing a --cut_ga hit to {PF00144, PF00905, PF13354} is
     retained in the focal family (logging dual_PF13354_PF00768_GA flag).
   - Non-focal PRE CL0013 {PF00768, PF02113} -> aux pool.
   - Non-PRE CL0013 {PF04960, PF06089, PF27310} -> outgroup pool.
4. Single-Round Official-Pfam-Seed MEME/FIMO 5.5.9 (p <= 1e-5):
   - Uses official seed MEME PWMs identically across viral & bacterial catalytic domain cores.
5. Taxonomically Balanced MMseqs2 Deduplication (--min-seq-id 0.85 -c 0.80) for Bacterial Pool:
   - Exports both full nr85 bacterial sets and balanced tree-ready bacterial sets (target <=1,500)
     preserving 100% of CARD/SARG clinical 85%-reps + all GTDB families + all bacterial phyla.
"""

import csv
import glob
import gzip
import os
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import hpc_config as hpc

BASE = hpc.BASE
OUT = BASE / "01_final_curated"
WORK = OUT / "work"
INPUTS = BASE / "inputs"

P1_DIR = BASE / "01_curated_datasets"
PRODIGAL_DIR = hpc.PRODIGAL_DIR
VIRAL_FAA = hpc.VIRAL_PFAM_DIR / "sbl_proteins.faa"
VIRAL_DOMTBL = hpc.VIRAL_PFAM_DIR / "sbl_proteins_pfam.domtbl.gz"
STEP1_DIR = hpc.STEP1_DIR
PFAM_CLANS = hpc.PFAM_CLANS_TSV

GTDB_DIR = BASE / "01b_bacterial_ref_db"
MAR_DIR = GTDB_DIR / "marine_mags"

FIMO_BIN = hpc.FIMO_BIN
MMSEQS_BIN = hpc.MMSEQS_BIN

FOCAL = ("PF00144", "PF00905", "PF13354")
AUX = ("PF00768", "PF02113")
OUTGROUP = ("PF04960", "PF06089", "PF27310")
OV_THRESH = 30
DOM_PAD = 10

CODON_STD = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
ALT_START = {"GTG", "TTG", "ATT", "ATC", "CTG"}
ASM_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


def open_text(path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8", errors="ignore") if p.endswith(".gz") else open(p, "rt", encoding="utf-8", errors="ignore")


def read_fasta(path):
    recs = []
    hdr, buf = None, []
    with open_text(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if hdr is not None:
                    recs.append((hdr.split()[0], hdr, "".join(buf)))
                hdr = s[1:].strip()
                buf = []
            else:
                buf.append(s)
        if hdr is not None:
            recs.append((hdr.split()[0], hdr, "".join(buf)))
    return recs


def write_fasta(path, items):
    with open(path, "w", encoding="utf-8") as f:
        for hdr, seq in items:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def translate_dna(dna, alt_start=True, tga_w=False):
    dna = dna.upper()
    n = len(dna) - (len(dna) % 3)
    out = []
    for i in range(0, n, 3):
        cod = dna[i : i + 3]
        if i == 0 and alt_start and cod in ALT_START:
            out.append("M")
        elif tga_w and cod == "TGA":
            out.append("W")
        else:
            out.append(CODON_STD.get(cod, "X"))
    return "".join(out)


def match_protein_to_dna(prot_seq, dna_seq):
    target = prot_seq.rstrip("*")
    dna = dna_seq.upper()
    cands = [dna]
    if len(dna) >= 3 and len(dna) % 3 == 0 and dna[-3:] in ("TAA", "TAG", "TGA"):
        cands.append(dna[:-3])
    for cand in cands:
        for tga_w in (False, True):
            aa = translate_dna(cand, alt_start=True, tga_w=tga_w).rstrip("*")
            if aa == target and len(cand) >= 3 * len(target):
                return cand[: 3 * len(target)]
            aa_lit = translate_dna(cand, alt_start=False, tga_w=tga_w).rstrip("*")
            if aa_lit == target and len(cand) >= 3 * len(target):
                return cand[: 3 * len(target)]
            if len(target) >= 30:
                idx = aa_lit.find(target)
                if idx >= 0:
                    return cand[idx * 3 : (idx + len(target)) * 3]
                if target.startswith("M") and len(target) > 30:
                    idx2 = aa_lit.find(target[1:])
                    if idx2 >= 1:
                        s_nt = cand[(idx2 - 1) * 3 : (idx2 - 1 + len(target)) * 3]
                        if len(s_nt) == 3 * len(target) and s_nt[:3] in (ALT_START | {"ATG"}):
                            return s_nt
    return None


def norm_tax(val, prefix):
    v = (val or "").strip()
    if not v or v in ("NA", "NONE", prefix):
        return f"{prefix}unclassified"
    return v if v.startswith(prefix) else f"{prefix}{v}"


def interval_overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def load_pfam_clans():
    acc_to_clan = {}
    with open_text(PFAM_CLANS) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                acc_to_clan[p[0]] = p[1] if p[1] else "NO_CLAN"
    return acc_to_clan


def parse_hmmsearch_domtbl(path, filter_ids=None):
    hits = defaultdict(list)
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 22:
                continue
            tid = p[0]
            if filter_ids is not None and tid not in filter_ids:
                continue
            hits[tid].append({
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


def resolve_focal_priority(hit_list, acc_to_clan):
    """Overlapping-envelope Pfam competition with Focal-Family Priority inside CL0013 clusters."""
    n = len(hit_list)
    if n == 0:
        return "drop_no_pfam_ga", "NONE", None, False, "NONE", ""
    adj = [[] for _ in range(n)]
    for i in range(n):
        a = hit_list[i]
        for j in range(i + 1, n):
            b = hit_list[j]
            if interval_overlap(a["env_from"], a["env_to"], b["env_from"], b["env_to"]) >= OV_THRESH:
                adj[i].append(j)
                adj[j].append(i)

    seen = [False] * n
    clusters = []
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
        top_hit = max(comp_hits, key=lambda h: (h["dom_score"], h["env_to"] - h["env_from"]))
        clusters.append((top_hit, comp_hits))

    clusters.sort(key=lambda x: x[0]["env_from"])
    arch_str = "+".join(
        f"{top['pfam_acc']}:{top['pfam_name']}({top['env_from']}-{top['env_to']})"
        for top, _ in clusters
    )

    cl13_clusters = [(top, ch) for top, ch in clusters if acc_to_clan.get(top["pfam_acc"], "NO_CLAN") == "CL0013"]
    if not cl13_clusters:
        top_non = max(clusters, key=lambda x: x[0]["dom_score"])[0]
        return "drop_non_CL0013_only", top_non["pfam_acc"], top_non, False, top_non["pfam_acc"], arch_str

    # Evaluate across CL0013-winning clusters using Focal-Family Priority
    best_top_cl13 = max((t for t, _ in cl13_clusters), key=lambda h: h["dom_score"])
    focal_cands = []
    aux_cands = []
    out_cands = []
    all_cl_accs = set()
    for top, ch in cl13_clusters:
        for h in ch:
            acc = h["pfam_acc"]
            if acc_to_clan.get(acc, "NO_CLAN") != "CL0013":
                continue
            all_cl_accs.add(acc)
            if acc in FOCAL:
                focal_cands.append(h)
            elif acc in AUX:
                aux_cands.append(h)
            else:
                out_cands.append(h)

    dual_flag = ("PF13354" in all_cl_accs) and ("PF00768" in all_cl_accs)
    if focal_cands:
        chosen = max(focal_cands, key=lambda h: (h["dom_score"], h["env_to"] - h["env_from"]))
        return f"focal_{chosen['pfam_acc']}", chosen["pfam_acc"], chosen, dual_flag, best_top_cl13["pfam_acc"], arch_str
    if aux_cands:
        chosen = max(aux_cands, key=lambda h: (h["dom_score"], h["env_to"] - h["env_from"]))
        return f"aux_{chosen['pfam_acc']}", chosen["pfam_acc"], chosen, dual_flag, best_top_cl13["pfam_acc"], arch_str
    chosen = max(out_cands, key=lambda h: (h["dom_score"], h["env_to"] - h["env_from"]))
    return f"outgroup_{chosen['pfam_acc']}", chosen["pfam_acc"], chosen, dual_flag, best_top_cl13["pfam_acc"], arch_str


def parse_fimo_tsv(path, safe_map):
    stats = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        _ = f.readline()
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            sid = safe_map.get(p[2], p[2])
            try:
                score = float(p[6])
                pval = float(p[7])
            except ValueError:
                continue
            if sid not in stats:
                stats[sid] = {"min_p": pval, "max_score": score, "motifs": {p[0]}}
            else:
                if pval < stats[sid]["min_p"]:
                    stats[sid]["min_p"] = pval
                if score > stats[sid]["max_score"]:
                    stats[sid]["max_score"] = score
                stats[sid]["motifs"].add(p[0])
    return stats


def run_cmd(cmd):
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")
    return p.stdout


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    acc_to_clan = load_pfam_clans()
    summary = ["=== Final Phase 1 + Phase 1B Curated Dataset Assembly Summary ==="]

    # ------------------------------------------------------------------
    # PART 1: VIRAL POOL (IMG/VR Scheme B + MOSAiC Euk-Purged)
    # ------------------------------------------------------------------
    print("[1/6] Loading IMG/VR Scheme B survivors & verifying viral CDSs...")
    imgvr_keep_meta = {}
    with open(P1_DIR / "Phase1_Viral_Decontamination_Audit.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["schemeB_status"].startswith("KEEP"):
                imgvr_keep_meta[r["sbl_protein_id"]] = r

    # Reconstruct verified viral CDSs from Prodigal orf/*.fna + ARG_cds.fna
    prodigal_by_uvig = defaultdict(list)
    target_uvigs = {r["uvig"] for r in imgvr_keep_meta.values()}
    for fna_p in sorted(glob.glob(str(PRODIGAL_DIR / "orf" / "*.fna"))):
        for oid, _, dna in read_fasta(fna_p):
            uv = oid.rsplit("_", 1)[0]
            if uv in target_uvigs:
                dna_u = dna.upper()
                prodigal_by_uvig[uv].append((
                    dna_u,
                    translate_dna(dna_u, True, False).rstrip("*"),
                    translate_dna(dna_u, False, False).rstrip("*"),
                ))

    arg_cds = {sid: dna.upper() for sid, _, dna in read_fasta(INPUTS / "ARG_cds.fna") if sid in imgvr_keep_meta}
    viral_prots = {uid: seq.rstrip("*") for uid, _, seq in read_fasta(VIRAL_FAA) if uid in imgvr_keep_meta}

    viral_cds_map = {}
    for uid, prot in viral_prots.items():
        uv = imgvr_keep_meta[uid]["uvig"]
        cslice = None
        for dna_u, aa_std, aa_lit in prodigal_by_uvig.get(uv, []):
            if aa_std == prot or aa_lit == prot:
                cslice = match_protein_to_dna(prot, dna_u)
                if cslice:
                    break
        if cslice is None:
            for dna_u, aa_std, aa_lit in prodigal_by_uvig.get(uv, []):
                if (len(prot) >= 30 and (prot in aa_lit or prot[1:] in aa_lit)) or (len(aa_lit) >= 30 and aa_lit in prot):
                    cslice = match_protein_to_dna(prot, dna_u)
                    if cslice:
                        break
        if cslice is None and uid in arg_cds:
            cslice = match_protein_to_dna(prot, arg_cds[uid])
        viral_cds_map[uid] = cslice

    v_hits = parse_hmmsearch_domtbl(VIRAL_DOMTBL, filter_ids=set(imgvr_keep_meta))

    # Load clean MOSAiC (25 seqs)
    mosaic_meta = {}
    with open(P1_DIR / "mosaic_audit" / "MOSAiC_SBL_Clean_Audit.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["final_mosaic_status"] == "KEEP_MOSAiC_CLEAN":
                mosaic_meta[r["seq_id"]] = r
    mosaic_prots = {sid: seq.rstrip("*") for sid, _, seq in read_fasta(P1_DIR / "mosaic_audit" / "mosaic_cl0013_candidates.faa") if sid in mosaic_meta}
    mosaic_hits = parse_hmmsearch_domtbl(P1_DIR / "mosaic_audit" / "mosaic_pfam_all.domtbl", filter_ids=set(mosaic_meta))

    viral_records = []
    v_cat_counts = Counter()
    for src, meta_d, prot_d, hit_d in [
        ("IMG_VR", imgvr_keep_meta, viral_prots, v_hits),
        ("MOSAiC", mosaic_meta, mosaic_prots, mosaic_hits),
    ]:
        for uid, prot in prot_d.items():
            cat, win_acc, dom, dual, top_cl, arch = resolve_focal_priority(hit_d.get(uid, []), acc_to_clan)
            v_cat_counts[(src, cat)] += 1
            if cat.startswith("drop_"):
                continue
            ds = max(1, dom["env_from"] - DOM_PAD)
            de = min(len(prot), dom["env_to"] + DOM_PAD)
            dom_seq = prot[ds - 1 : de]
            cds = viral_cds_map.get(uid) if src == "IMG_VR" else None
            dom_cds = cds[(ds - 1) * 3 : de * 3] if (cds and len(cds) >= 3 * de) else None
            # Build clean tree-safe uid
            safe_uid = f"MOSAIC_{uid}" if src == "MOSAiC" else uid
            safe_uid = re.sub(r"[\s:(),;\[\]'\"]+", "_", safe_uid)
            viral_records.append({
                "uid": safe_uid,
                "orig_uid": uid,
                "source": src,
                "pfam_category": cat,
                "winning_acc": win_acc,
                "top_overall_cl0013_acc": top_cl,
                "dual_PF13354_PF00768_GA": dual,
                "dom_score": dom["dom_score"],
                "seq_len": len(prot),
                "env_from": dom["env_from"],
                "env_to": dom["env_to"],
                "dom_start": ds,
                "dom_end": de,
                "prot": prot,
                "dom_seq": dom_seq,
                "cds": cds,
                "dom_cds": dom_cds,
                "cds_verified": cds is not None,
                "has_SXXK_in_dom": bool(re.search(r"S[A-Z]{2}K", dom_seq)),
                "domain_architecture": arch,
                "vir_sub_info": meta_d[uid].get("schemeB_status", "") if src == "IMG_VR" else meta_d[uid].get("vir_bin", ""),
            })

    print("Viral Focal-Priority breakdown:", dict(v_cat_counts))
    summary.append("Viral Pfam-A (--cut_ga) Focal-Priority Breakdown:")
    for k, v in sorted(v_cat_counts.items()):
        summary.append(f"  {k}: {v}")

    # ------------------------------------------------------------------
    # PART 2: TRIPLE-SOURCE BACTERIAL POOL
    # ------------------------------------------------------------------
    print("[2/6] Loading Triple-Source Bacterial Pool (GTDB_3000 + Marine_MAGs_259 + CARD_SARG_Clinical)...")
    # 2a. Merged geNomad across 8 GTDB shards
    pv_ivs = defaultdict(list)
    full_vir = set()
    plasmids = set()
    for i in range(8):
        cdir = GTDB_DIR / "genomad_out" / f"chunk_{i}"
        with open(cdir / f"chunk_{i}_find_proviruses" / f"chunk_{i}_provirus.tsv", "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                pv_ivs[r["source_seq"]].append((int(r["start"]), int(r["end"])))
        with open(cdir / f"chunk_{i}_summary" / f"chunk_{i}_virus_summary.tsv", "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                p = r["seq_name"].split("|")
                base = f"{p[0]}|{p[1]}"
                topo = r.get("topology", "")
                coords = r.get("coordinates", "NA")
                if topo == "Provirus" or (coords and coords != "NA" and "-" in coords):
                    for seg in coords.split(","):
                        m = re.match(r"^(\d+)-(\d+)$", seg.strip())
                        if m:
                            pv_ivs[base].append((int(m.group(1)), int(m.group(2))))
                else:
                    full_vir.add(base)
        with open(cdir / f"chunk_{i}_summary" / f"chunk_{i}_plasmid_summary.tsv", "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                p = r["seq_name"].split("|")
                plasmids.add(f"{p[0]}|{p[1]}")

    bac_records = []
    b_counts = Counter()

    # Source 1: GTDB_RefSeq_3000
    gtdb_hits_tsv = {}
    with open(GTDB_DIR / "sbl_contigs" / "bac_gtdb_sbl_hits.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gtdb_hits_tsv[r["uid"]] = r
    gtdb_prots = {uid: seq.rstrip("*") for uid, _, seq in read_fasta(GTDB_DIR / "sbl_contigs" / "bac_gtdb_sbl_proteins.faa")}
    gtdb_cds = {uid: seq.upper() for uid, _, seq in read_fasta(GTDB_DIR / "sbl_contigs" / "bac_gtdb_sbl_cds_verified.fna")}
    gtdb_domtbl = parse_hmmsearch_domtbl(GTDB_DIR / "sbl_contigs" / "bac_gtdb_pfam_all.domtbl", filter_ids=set(gtdb_hits_tsv))

    for uid, r in gtdb_hits_tsv.items():
        # Verify ALL loci of this protein aren't solely inside prophage; keep if >=1 chromosomal/plasmid non-prophage locus exists
        loci_strs = [x for x in r["all_verified_loci"].split(";") if ":" in x and x != "NA"]
        if not loci_strs:
            key = f"{r['acc']}|{r['contig_acc']}"
            st, en = int(r["cds_start"]), int(r["cds_end"])
            in_pv = (key in full_vir) or (key in pv_ivs and st > 0 and any(interval_overlap(st, en, ps, pe) > 0 for ps, pe in pv_ivs[key]))
            if in_pv:
                b_counts[("GTDB_3000", "DROP_prophage")] += 1
                continue
            replicon = "plasmid" if key in plasmids else "chromosome"
        else:
            clean_rep = []
            for loc in loci_strs:
                contig_c, span = loc.split(":", 1)
                m = re.match(r"^(\d+)\.\.(\d+)", span)
                st, en = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
                key = f"{r['acc']}|{contig_c}"
                in_pv = (key in full_vir) or (key in pv_ivs and st > 0 and any(interval_overlap(st, en, ps, pe) > 0 for ps, pe in pv_ivs[key]))
                if not in_pv:
                    clean_rep.append("plasmid" if key in plasmids else "chromosome")
            if not clean_rep:
                b_counts[("GTDB_3000", "DROP_prophage")] += 1
                continue
            replicon = "chromosome" if "chromosome" in clean_rep else "plasmid"

        prot = gtdb_prots[uid]
        cat, win_acc, dom, dual, top_cl, arch = resolve_focal_priority(gtdb_domtbl.get(uid, []), acc_to_clan)
        b_counts[("GTDB_3000", cat)] += 1
        if cat.startswith("drop_"):
            continue
        ds = max(1, dom["env_from"] - DOM_PAD)
        de = min(len(prot), dom["env_to"] + DOM_PAD)
        dom_seq = prot[ds - 1 : de]
        cds = gtdb_cds.get(uid)
        dom_cds = cds[(ds - 1) * 3 : de * 3] if (cds and len(cds) >= 3 * de) else None
        safe_uid = re.sub(r"[\s:(),;\[\]'\"]+", "_", f"GTDB|{uid}")
        bac_records.append({
            "uid": safe_uid,
            "orig_uid": uid,
            "source_db": "GTDB_RefSeq_3000",
            "source_priority": 2,
            "asm_rank": ASM_RANK.get(r["assembly_level"], 4),
            "pfam_category": cat,
            "winning_acc": win_acc,
            "top_overall_cl0013_acc": top_cl,
            "dual_PF13354_PF00768_GA": dual,
            "dom_score": dom["dom_score"],
            "seq_len": len(prot),
            "env_from": dom["env_from"],
            "env_to": dom["env_to"],
            "dom_start": ds,
            "dom_end": de,
            "prot": prot,
            "dom_seq": dom_seq,
            "cds": cds,
            "dom_cds": dom_cds,
            "cds_verified": cds is not None,
            "replicon": replicon,
            "phylum": norm_tax(r["gtdb_phylum"], "p__"),
            "family": norm_tax(r["gtdb_family"], "f__"),
            "genus": norm_tax(r["gtdb_genus"], "g__"),
            "species": r["gtdb_species"],
            "is_clinical_mobile_bla": False,
            "has_SXXK_in_dom": bool(re.search(r"S[A-Z]{2}K", dom_seq)),
            "domain_architecture": arch,
        })

    # Source 2: Marine_MAGs_259
    mar_rows = {}
    with open(MAR_DIR / "Marine_MAGs_SBL_Audit.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["final_status"].startswith("KEEP_"):
                mar_rows[r["uid"]] = r
    mar_prots = {uid: seq.rstrip("*") for uid, _, seq in read_fasta(MAR_DIR / "marine_mags_sbl_proteins.faa")}
    mar_cds = {uid: seq.upper() for uid, _, seq in read_fasta(MAR_DIR / "marine_mags_sbl_cds_verified.fna")}
    mar_domtbl = parse_hmmsearch_domtbl(MAR_DIR / "marine_mags_pfam_all.domtbl", filter_ids=set(mar_rows))

    for uid, r in mar_rows.items():
        prot = mar_prots[uid]
        cat, win_acc, dom, dual, top_cl, arch = resolve_focal_priority(mar_domtbl.get(uid, []), acc_to_clan)
        b_counts[("MarineMAG_259", cat)] += 1
        if cat.startswith("drop_"):
            continue
        ds = max(1, dom["env_from"] - DOM_PAD)
        de = min(len(prot), dom["env_to"] + DOM_PAD)
        dom_seq = prot[ds - 1 : de]
        cds = mar_cds.get(uid)
        dom_cds = cds[(ds - 1) * 3 : de * 3] if (cds and len(cds) >= 3 * de) else None
        safe_uid = re.sub(r"[\s:(),;\[\]'\"]+", "_", f"MMAG|{uid}")
        bac_records.append({
            "uid": safe_uid,
            "orig_uid": uid,
            "source_db": "Marine_MAGs_259",
            "source_priority": 3,
            "asm_rank": 5,
            "pfam_category": cat,
            "winning_acc": win_acc,
            "top_overall_cl0013_acc": top_cl,
            "dual_PF13354_PF00768_GA": dual,
            "dom_score": dom["dom_score"],
            "seq_len": len(prot),
            "env_from": dom["env_from"],
            "env_to": dom["env_to"],
            "dom_start": ds,
            "dom_end": de,
            "prot": prot,
            "dom_seq": dom_seq,
            "cds": cds,
            "dom_cds": dom_cds,
            "cds_verified": cds is not None,
            "replicon": r["genomad_status"],
            "phylum": norm_tax(r["gtdb_phylum"], "p__"),
            "family": norm_tax(r["gtdb_family"], "f__"),
            "genus": norm_tax(r["gtdb_genus"], "g__"),
            "species": r["gtdb_species"],
            "is_clinical_mobile_bla": False,
            "has_SXXK_in_dom": bool(re.search(r"S[A-Z]{2}K", dom_seq)),
            "domain_architecture": arch,
        })

    # Source 3: CARD_SARG_Clinical (1,219 CleanCore non-NR subset)
    sarg_seqs, sarg_hdrs = {}, {}
    for fam in FOCAL:
        for sid, hdr, seq in read_fasta(INPUTS / f"{fam}_SARG.faa"):
            sarg_seqs[(fam, sid)] = seq.rstrip("*")
            sarg_hdrs[(fam, sid)] = hdr
    idmap_bac = {}
    with open(STEP1_DIR / "step1_idmap.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["origin"] != "virus":
                idmap_bac[(r["assigned_family"], r["orig_id"])] = r["safe_id"]

    card_rows = []
    with open(INPUTS / "PerSequence_Bacterial_Audit.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["pass_CleanCore_TierB"] == "True" and r["source_type"] in ("CARD_SARG_curated", "SARG_other"):
                card_rows.append(r)
    card_safe_ids = {idmap_bac[(r["family"], r["seq_id"])] for r in card_rows}
    step1_domtbl = parse_hmmsearch_domtbl(STEP1_DIR / "step1_pfam.domtbl.gz", filter_ids=card_safe_ids)

    seen_card = set()
    for r in card_rows:
        ofam, sid = r["family"], r["seq_id"]
        safe_id = idmap_bac[(ofam, sid)]
        prot = sarg_seqs[(ofam, sid)]
        hdr = sarg_hdrs[(ofam, sid)]
        cat, win_acc, dom, dual, top_cl, arch = resolve_focal_priority(step1_domtbl.get(safe_id, []), acc_to_clan)
        b_counts[("CARD_SARG", cat)] += 1
        if cat.startswith("drop_"):
            continue
        ds = max(1, dom["env_from"] - DOM_PAD)
        de = min(len(prot), dom["env_to"] + DOM_PAD)
        dom_seq = prot[ds - 1 : de]
        is_pl = bool(re.search(r"plasmid|transposon|integron|conjugat|mobile", hdr, re.I))
        is_clin = bool(re.search(r"\b(TEM|SHV|CTX-M|OXA|KPC|PER|VEB|GES|CMY|DHA|FOX|ACC|MIR|ACT|SME|IMI|NMC|BEL|TLA|CARB|PSE|ROB|BRO)[-_0-9]", hdr, re.I))
        # Deduplicate exact identical sequence+family if SARG files overlapped
        raw_tag = f"SARG|{sid}" if sid not in seen_card else f"SARG|{ofam}|{sid}"
        seen_card.add(sid)
        safe_uid = re.sub(r"[\s:(),;\[\]'\"]+", "_", raw_tag)
        bac_records.append({
            "uid": safe_uid,
            "orig_uid": sid,
            "source_db": "CARD_SARG_Clinical",
            "source_priority": 1 if r["source_type"] == "CARD_SARG_curated" else 1.5,
            "asm_rank": 0,
            "pfam_category": cat,
            "winning_acc": win_acc,
            "top_overall_cl0013_acc": top_cl,
            "dual_PF13354_PF00768_GA": dual,
            "dom_score": dom["dom_score"],
            "seq_len": len(prot),
            "env_from": dom["env_from"],
            "env_to": dom["env_to"],
            "dom_start": ds,
            "dom_end": de,
            "prot": prot,
            "dom_seq": dom_seq,
            "cds": None,
            "dom_cds": None,
            "cds_verified": False,
            "replicon": "plasmid_or_MGE" if (is_pl or is_clin) else "chromosome_or_unspecified",
            "phylum": norm_tax(r["sarg_phylum"], "p__"),
            "family": "f__CARD_SARG_clinical",
            "genus": "g__CARD_SARG",
            "species": r["sarg_species"],
            "is_clinical_mobile_bla": is_clin or is_pl,
            "has_SXXK_in_dom": bool(re.search(r"S[A-Z]{2}K", dom_seq)),
            "domain_architecture": arch,
        })

    print("Bacterial Triple-Source breakdown:", dict(b_counts))
    summary.append("\nBacterial Triple-Source Pfam-A (--cut_ga) Breakdown:")
    for k, v in sorted(b_counts.items()):
        summary.append(f"  {k}: {v}")

    # ------------------------------------------------------------------
    # PART 3: SINGLE-ROUND OFFICIAL-PFAM-SEED FIMO (p <= 1e-5)
    # ------------------------------------------------------------------
    print("[3/6] Running unified single-round official-Pfam-seed FIMO (p <= 1e-5)...")
    fimo_stats = {}
    for fam in FOCAL:
        meme_txt = P1_DIR / "meme_fimo_work" / fam / "r1_meme" / "meme.txt"
        v_cands = [r for r in viral_records if r["pfam_category"] == f"focal_{fam}"]
        b_cands = [r for r in bac_records if r["pfam_category"] == f"focal_{fam}"]

        for kind, cands in [("viral", v_cands), ("bacterial", b_cands)]:
            safe_map = {}
            fa_items = []
            for idx, r in enumerate(cands):
                sid = f"F_{idx:06d}"
                safe_map[sid] = r["uid"]
                fa_items.append((sid, r["dom_seq"]))
            fa_p = WORK / f"{fam}_{kind}_dom_for_fimo.faa"
            fdir = WORK / f"fimo_{fam}_{kind}"
            write_fasta(fa_p, fa_items)
            run_cmd(f"{FIMO_BIN} --thresh 1e-4 --oc {fdir} {meme_txt} {fa_p}")
            st_map = parse_fimo_tsv(fdir / "fimo.tsv", safe_map)
            for r in cands:
                st = st_map.get(r["uid"], {"min_p": 1.0, "max_score": 0.0, "motifs": set()})
                r["fimo_r1_min_p"] = st["min_p"]
                r["fimo_r1_n_motifs"] = len(st["motifs"])
                r["pass_fimo_1e5"] = st["min_p"] <= 1e-5
                r["pass_fimo_1e10"] = st["min_p"] <= 1e-10
                fimo_stats[(fam, kind, "total")] = fimo_stats.get((fam, kind, "total"), 0) + 1
                if r["pass_fimo_1e5"]:
                    fimo_stats[(fam, kind, "p1e5")] = fimo_stats.get((fam, kind, "p1e5"), 0) + 1
                if r["pass_fimo_1e10"]:
                    fimo_stats[(fam, kind, "p1e10")] = fimo_stats.get((fam, kind, "p1e10"), 0) + 1

    print("Single-round official-seed FIMO (p <= 1e-5) summary:")
    summary.append("\nSingle-Round Official-Seed FIMO (p <= 1e-5) Summary:")
    for fam in FOCAL:
        for kind in ("viral", "bacterial"):
            tot = fimo_stats.get((fam, kind, "total"), 0)
            p5 = fimo_stats.get((fam, kind, "p1e5"), 0)
            p10 = fimo_stats.get((fam, kind, "p1e10"), 0)
            line = f"  {fam} {kind:9s}: {p5:,} / {tot:,} ({100*p5/max(tot,1):.2f}%) pass p<=1e-5 (p<=1e-10: {p10:,})"
            print(line)
            summary.append(line)

    # ------------------------------------------------------------------
    # PART 4: EXPORT VIRAL FOCAL FASTAs & RUN MMseqs2 BALANCED BACTERIAL DEDUP
    # ------------------------------------------------------------------
    print("[4/6] Exporting final viral pools & running MMseqs2 (--min-seq-id 0.85 -c 0.80) on Triple-Source bacterial pools...")
    summary.append("\nFinal Exported Focal Datasets (Viral vs Bacterial nr85 & Tree-Balanced):")

    bac_keep_uids_nr85 = set()
    bac_keep_uids_balanced = set()
    bac_cluster_id = {}

    for fam in FOCAL:
        v_pass = [r for r in viral_records if r["pfam_category"] == f"focal_{fam}" and r["pass_fimo_1e5"]]
        write_fasta(OUT / f"{fam}_viral_clean.faa", [(r["uid"], r["prot"]) for r in v_pass])
        write_fasta(OUT / f"{fam}_viral_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_pass])
        v_cds = [(r["uid"], r["cds"]) for r in v_pass if r["cds"] is not None]
        v_dcds = [(r["uid"], r["dom_cds"]) for r in v_pass if r["dom_cds"] is not None]
        write_fasta(OUT / f"{fam}_viral_cds.fna", v_cds)
        write_fasta(OUT / f"{fam}_viral_dom_cds.fna", v_dcds)

        b_pass = [r for r in bac_records if r["pfam_category"] == f"focal_{fam}" and r["pass_fimo_1e5"]]
        b_by_uid = {r["uid"]: r for r in b_pass}
        b_dom_fa = WORK / f"{fam}_bac_all_fimo1e5_dom.faa"
        write_fasta(b_dom_fa, [(r["uid"], r["dom_seq"]) for r in b_pass])

        clu_pref = WORK / f"mmseqs85_{fam}"
        tmp_d = WORK / f"tmp_{fam}"
        tmp_d.mkdir(exist_ok=True)
        run_cmd(f"{MMSEQS_BIN} easy-cluster {b_dom_fa} {clu_pref} {tmp_d} --min-seq-id 0.85 -c 0.80 --cov-mode 0 --threads 16 >/dev/null")

        clusters = defaultdict(list)
        with open(f"{clu_pref}_cluster.tsv", "r", encoding="utf-8") as f:
            for line in f:
                rep, mem = line.rstrip("\n").split("\t")
                clusters[rep].append(mem)

        # Deterministic best representative per 85% cluster: CARD/SARG > GTDB (CDS verified, Complete>Chr) > MarineMAG (CDS verified)
        nr85_reps = []
        hyphy_cds_reps = []
        for idx_c, (rep_k, mems) in enumerate(sorted(clusters.items())):
            m_recs = [b_by_uid[m] for m in mems if m in b_by_uid]
            best_r = min(
                m_recs,
                key=lambda x: (
                    x["source_priority"],
                    0 if x["cds_verified"] else 1,
                    x["asm_rank"],
                    -x["dom_score"],
                    x["uid"],
                ),
            )
            # Inherit exact verified CDS/dom_cds if an identical GTDB/MarineMAG sequence sits in the same cluster
            if not best_r["cds_verified"]:
                exact_prot_donors = [m for m in m_recs if m["cds_verified"] and m["prot"] == best_r["prot"] and m["dom_cds"]]
                if exact_prot_donors:
                    d0 = min(exact_prot_donors, key=lambda x: (x["source_priority"], x["asm_rank"], x["uid"]))
                    best_r["cds"] = d0["cds"]
                    best_r["dom_cds"] = d0["dom_cds"]
                    best_r["cds_verified"] = True
                else:
                    exact_dom_donors = [m for m in m_recs if m["cds_verified"] and m["dom_seq"] == best_r["dom_seq"] and m["dom_cds"]]
                    if exact_dom_donors:
                        d0 = min(exact_dom_donors, key=lambda x: (x["source_priority"], x["asm_rank"], x["uid"]))
                        best_r["dom_cds"] = d0["dom_cds"]
            # Best CDS-verified member in cluster for Phase 5 HyPhy codon analyses
            cds_mems = [m for m in m_recs if m["cds_verified"] and m["dom_cds"]]
            if cds_mems:
                best_cds = min(cds_mems, key=lambda x: (x["asm_rank"], x["source_priority"], -x["dom_score"], x["uid"]))
                hyphy_cds_reps.append(best_cds)
            nr85_reps.append(best_r)
            bac_keep_uids_nr85.add(best_r["uid"])
            for m in mems:
                bac_cluster_id[m] = f"{fam}_c85_{idx_c:05d}"

        write_fasta(OUT / f"{fam}_bacterial_nr85_clean.faa", [(r["uid"], r["prot"]) for r in nr85_reps])
        write_fasta(OUT / f"{fam}_bacterial_nr85_dom.faa", [(r["uid"], r["dom_seq"]) for r in nr85_reps])
        write_fasta(OUT / f"{fam}_bacterial_nr85_cds.fna", [(r["uid"], r["cds"]) for r in nr85_reps if r["cds"]])
        write_fasta(OUT / f"{fam}_bacterial_nr85_dom_cds.fna", [(r["uid"], r["dom_cds"]) for r in nr85_reps if r["dom_cds"]])
        write_fasta(OUT / f"{fam}_bacterial_nr85_hyphy_dom.faa", [(r["uid"], r["dom_seq"]) for r in hyphy_cds_reps])
        write_fasta(OUT / f"{fam}_bacterial_nr85_hyphy_dom_cds.fna", [(r["uid"], r["dom_cds"]) for r in hyphy_cds_reps])
        write_fasta(OUT / f"{fam}_bacterial_nr85_hyphy_cds.fna", [(r["uid"], r["cds"]) for r in hyphy_cds_reps if r["cds"]])

        # Balanced tree representative set (target 1,500 if len(nr85_reps) > 1600)
        TREE_CAP = 1500
        if len(nr85_reps) <= 1600:
            bal_reps = nr85_reps
        else:
            sel = []
            sel_uids = set()
            # 1. Keep 100% CARD/SARG clinical 85%-cluster reps
            card_reps = [r for r in nr85_reps if r["source_db"] == "CARD_SARG_Clinical"]
            card_reps.sort(key=lambda x: (-x["dom_score"], x["uid"]))
            for r in card_reps:
                sel.append(r)
                sel_uids.add(r["uid"])
            # 2. Keep #1 best GTDB rep per unique GTDB family, round-robin across phyla so all phyla stay represented
            gtdb_by_fam = defaultdict(list)
            for r in nr85_reps:
                if r["source_db"] == "GTDB_RefSeq_3000" and r["uid"] not in sel_uids:
                    gtdb_by_fam[r["family"]].append(r)
            fam_by_phy = defaultdict(list)
            for fk, lst in gtdb_by_fam.items():
                lst.sort(key=lambda x: (x["asm_rank"], -x["dom_score"], x["uid"]))
                fam_by_phy[lst[0]["phylum"]].append(lst[0])
            for phy in fam_by_phy:
                fam_by_phy[phy].sort(key=lambda x: (x["asm_rank"], -x["dom_score"], x["uid"]))
            f_idx = defaultdict(int)
            while len(sel) < TREE_CAP:
                prog = 0
                for phy in sorted(fam_by_phy):
                    if len(sel) >= TREE_CAP:
                        break
                    lst = fam_by_phy[phy]
                    if f_idx[phy] < len(lst):
                        r = lst[f_idx[phy]]
                        f_idx[phy] += 1
                        sel.append(r)
                        sel_uids.add(r["uid"])
                        prog += 1
                if prog == 0:
                    break
            # 3. Keep up to 150 diverse Marine MAG 85%-reps spread across phyla
            mar_by_phy = defaultdict(list)
            for r in nr85_reps:
                if r["source_db"] == "Marine_MAGs_259" and r["uid"] not in sel_uids:
                    mar_by_phy[r["phylum"]].append(r)
            for phy in mar_by_phy:
                mar_by_phy[phy].sort(key=lambda x: (-x["dom_score"], x["uid"]))
            mar_added = 0
            mar_idx = defaultdict(int)
            while mar_added < 150 and len(sel) < TREE_CAP:
                prog = 0
                for phy in sorted(mar_by_phy):
                    if mar_added >= 150 or len(sel) >= TREE_CAP:
                        break
                    lst = mar_by_phy[phy]
                    if mar_idx[phy] < len(lst):
                        r = lst[mar_idx[phy]]
                        mar_idx[phy] += 1
                        sel.append(r)
                        sel_uids.add(r["uid"])
                        mar_added += 1
                        prog += 1
                if prog == 0:
                    break
            # 4. Fill remaining slots round-robin across bacterial phyla from remaining nr85_reps
            rem_by_phy = defaultdict(list)
            for r in nr85_reps:
                if r["uid"] not in sel_uids:
                    rem_by_phy[r["phylum"]].append(r)
            for phy in rem_by_phy:
                rem_by_phy[phy].sort(key=lambda x: (x["source_priority"], x["asm_rank"], -x["dom_score"], x["uid"]))
            r_idx = defaultdict(int)
            while len(sel) < TREE_CAP:
                prog = 0
                for phy in sorted(rem_by_phy):
                    if len(sel) >= TREE_CAP:
                        break
                    lst = rem_by_phy[phy]
                    if r_idx[phy] < len(lst):
                        r = lst[r_idx[phy]]
                        r_idx[phy] += 1
                        sel.append(r)
                        sel_uids.add(r["uid"])
                        prog += 1
                if prog == 0:
                    break
            bal_reps = sel

        for r in bal_reps:
            bac_keep_uids_balanced.add(r["uid"])

        write_fasta(OUT / f"{fam}_bacterial_tree_clean.faa", [(r["uid"], r["prot"]) for r in bal_reps])
        write_fasta(OUT / f"{fam}_bacterial_tree_dom.faa", [(r["uid"], r["dom_seq"]) for r in bal_reps])
        write_fasta(OUT / f"{fam}_bacterial_tree_cds.fna", [(r["uid"], r["cds"]) for r in bal_reps if r["cds"]])
        write_fasta(OUT / f"{fam}_bacterial_tree_dom_cds.fna", [(r["uid"], r["dom_cds"]) for r in bal_reps if r["dom_cds"]])
        write_fasta(OUT / f"{fam}_combined_tree_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_pass] + [(r["uid"], r["dom_seq"]) for r in bal_reps])

        v_src_c = Counter(r["source"] for r in v_pass)
        b_src_nr = Counter(r["source_db"] for r in nr85_reps)
        b_src_tr = Counter(r["source_db"] for r in bal_reps)
        n_dual = sum(1 for r in v_pass if r["dual_PF13354_PF00768_GA"])
        n_phy = len({r["phylum"] for r in bal_reps})
        sline = (
            f"  {fam}: Viral={len(v_pass):,} ({dict(v_src_c)}, CDS={len(v_cds):,}/{v_src_c.get('IMG_VR',0)}, dual_PF00768={n_dual}) | "
            f"Bacterial raw_fimo={len(b_pass):,} -> nr85={len(nr85_reps):,} ({dict(b_src_nr)}) -> "
            f"tree_balanced={len(bal_reps):,} ({dict(b_src_tr)}, phyla={n_phy}, cds_verified={sum(1 for r in bal_reps if r['cds'])})"
        )
        print(sline)
        summary.append(sline)

    # ------------------------------------------------------------------
    # PART 5: AUXILIARY CL0013 (PF00768/PF02113) & OUTGROUP (PF04960) EXPORTS
    # ------------------------------------------------------------------
    print("[5/6] Exporting Auxiliary CL0013 (PF00768/PF02113) & PF04960 Biological Outgroups...")
    aux_v = [r for r in viral_records if r["pfam_category"].startswith("aux_")]
    aux_b = [r for r in bac_records if r["pfam_category"].startswith("aux_")]
    out_b = [r for r in bac_records if r["pfam_category"] == "outgroup_PF04960"]

    write_fasta(OUT / "CL0013_aux_viral_clean.faa", [(r["uid"], r["prot"]) for r in aux_v])
    write_fasta(OUT / "CL0013_aux_viral_dom.faa", [(r["uid"], r["dom_seq"]) for r in aux_v])

    # Pick taxonomically balanced bacterial aux (MMseqs2 60%) + 20 PF04960 complete-genome outgroups
    def cluster_reps(recs, tag, id_thresh=0.60, max_keep=None):
        if not recs:
            return []
        by_u = {r["uid"]: r for r in recs}
        fa_p = WORK / f"{tag}_for_clu.faa"
        write_fasta(fa_p, [(r["uid"], r["dom_seq"]) for r in recs])
        cpref = WORK / f"mmseqs_{tag}"
        td = WORK / f"tmp_{tag}"
        td.mkdir(exist_ok=True)
        run_cmd(f"{MMSEQS_BIN} easy-cluster {fa_p} {cpref} {td} --min-seq-id {id_thresh} -c 0.80 --cov-mode 0 --threads 16 >/dev/null")
        clus = defaultdict(list)
        with open(f"{cpref}_cluster.tsv", "r", encoding="utf-8") as f:
            for line in f:
                rep, mem = line.rstrip("\n").split("\t")
                clus[rep].append(mem)
        reps = []
        for rep_k, mems in sorted(clus.items()):
            best = min((by_u[m] for m in mems), key=lambda x: (x["asm_rank"], -x["dom_score"], x["uid"]))
            reps.append(best)
        if max_keep and len(reps) > max_keep:
            # Round-robin across phyla
            by_phy = defaultdict(list)
            for r in reps:
                by_phy[r["phylum"]].append(r)
            for p in by_phy:
                by_phy[p].sort(key=lambda x: (x["asm_rank"], -x["dom_score"], x["uid"]))
            sel = []
            idx = defaultdict(int)
            while len(sel) < max_keep:
                prog = 0
                for p in sorted(by_phy):
                    if len(sel) >= max_keep:
                        break
                    if idx[p] < len(by_phy[p]):
                        sel.append(by_phy[p][idx[p]])
                        idx[p] += 1
                        prog += 1
                if prog == 0:
                    break
            reps = sel
        return reps

    aux_b_reps = cluster_reps(aux_b, "bac_aux60", id_thresh=0.60, max_keep=300)
    out_b_reps = cluster_reps(out_b, "bac_pf04960_out", id_thresh=0.50, max_keep=24)

    write_fasta(OUT / "CL0013_aux_bacterial_rep_dom.faa", [(r["uid"], r["dom_seq"]) for r in aux_b_reps])
    write_fasta(OUT / "CL0013_aux_bacterial_rep_clean.faa", [(r["uid"], r["prot"]) for r in aux_b_reps])
    write_fasta(OUT / "CL0013_outgroup_PF04960_dom.faa", [(r["uid"], r["dom_seq"]) for r in out_b_reps])
    write_fasta(OUT / "CL0013_outgroup_PF04960_clean.faa", [(r["uid"], r["prot"]) for r in out_b_reps])

    summary.append(f"\nAuxiliary & Biological Outgroup Exports:")
    summary.append(f"  Viral aux_PF00768: {len(aux_v)} (IMG/VR + MOSAiC)")
    summary.append(f"  Bacterial aux (PF00768+PF02113) raw={len(aux_b):,} -> diverse 60%-reps={len(aux_b_reps)}")
    summary.append(f"  Bacterial non-PRE outgroup PF04960:Glutaminase raw={len(out_b):,} -> diverse 50%-reps={len(out_b_reps)}")

    # ------------------------------------------------------------------
    # PART 6: WRITE UNIFIED AUDIT CATALOGS & SUMMARY REPORT
    # ------------------------------------------------------------------
    print("[6/6] Writing Final_Viral_SBL_Catalog.tsv, Final_Bacterial_TripleSource_Catalog.tsv & Summary...")
    v_cols = [
        "uid", "orig_uid", "source", "pfam_category", "winning_acc", "top_overall_cl0013_acc",
        "dual_PF13354_PF00768_GA", "dom_score", "seq_len", "env_from", "env_to", "dom_start", "dom_end",
        "fimo_r1_min_p", "fimo_r1_n_motifs", "pass_fimo_1e5", "pass_fimo_1e10",
        "cds_verified", "has_SXXK_in_dom", "vir_sub_info", "domain_architecture",
    ]
    with open(OUT / "Final_Viral_SBL_Catalog.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=v_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in viral_records:
            row = dict(r)
            row.setdefault("fimo_r1_min_p", "NA")
            row.setdefault("fimo_r1_n_motifs", "NA")
            row.setdefault("pass_fimo_1e5", "AUX_EXEMPT")
            row.setdefault("pass_fimo_1e10", "AUX_EXEMPT")
            w.writerow(row)

    b_cols = [
        "uid", "orig_uid", "source_db", "pfam_category", "winning_acc", "top_overall_cl0013_acc",
        "dual_PF13354_PF00768_GA", "dom_score", "seq_len", "env_from", "env_to", "dom_start", "dom_end",
        "fimo_r1_min_p", "fimo_r1_n_motifs", "pass_fimo_1e5", "pass_fimo_1e10",
        "mmseqs85_cluster_id", "in_nr85_rep", "in_tree_balanced",
        "cds_verified", "replicon", "is_clinical_mobile_bla", "phylum", "family", "genus", "species",
        "has_SXXK_in_dom", "domain_architecture",
    ]
    with open(OUT / "Final_Bacterial_TripleSource_Catalog.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=b_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in bac_records:
            row = dict(r)
            row.setdefault("fimo_r1_min_p", "NA")
            row.setdefault("fimo_r1_n_motifs", "NA")
            row.setdefault("pass_fimo_1e5", "NON_FOCAL")
            row.setdefault("pass_fimo_1e10", "NON_FOCAL")
            row["mmseqs85_cluster_id"] = bac_cluster_id.get(r["uid"], "NA")
            row["in_nr85_rep"] = r["uid"] in bac_keep_uids_nr85
            row["in_tree_balanced"] = r["uid"] in bac_keep_uids_balanced
            w.writerow(row)

    rep = "\n".join(summary) + "\n"
    (OUT / "Final_Phase1_Dataset_Summary.txt").write_text(rep, encoding="utf-8")
    print("[Done] Final Phase 1 + Phase 1B assembly complete!")


if __name__ == "__main__":
    main()
