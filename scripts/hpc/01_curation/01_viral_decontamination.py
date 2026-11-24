#!/usr/bin/env python3
"""Phase 1 Curation & Joint Decontamination Pipeline for Viral/Bacterial SBLs.

Implements:
A. Scheme B viral decontamination (geNomad v1.7 + CheckV v1.5 joint coordinate filter):
   - Drop geNomad plasmid contigs
   - Drop SBL loci outside geNomad provirus boundaries (host flank)
   - Drop SBL loci outside CheckV viral regions / inside host flanks
   - Drop non-virus geNomad contigs dominated by chromosome/plasmid score
   - Rescue non-virus geNomad virus-score-max fragments ONLY IF CheckV viral_genes >= 1 and host_genes == 0
B. Overlapping-envelope full Pfam-A (--cut_ga) domain competition (ov >= 30 aa):
   - Preserves genuine non-overlapping multi-domain fusions (e.g., PF00912+PF00905, PF38217+PF00144)
   - Routes focal CL0013 winners to {PF00144, PF00905, PF13354} and auxiliary PRE CL0013 winners to {PF00768, PF02113}
C. Unified two-round MEME/FIMO 5.5.9 motif filtering (p <= 1e-5; no hard SXXK regex pre-filter)
D. Bacterial CleanCore plasmid/MGE tagging (retaining authentic bacterial plasmid SBLs)
E. Exports clean full proteins, catalytic domain cores (env_from-10..env_to+10),
   1:1 translation-verified viral CDS nucleotides & domain CDS slices, and auxiliary CL0013 sequences.
"""

import csv
import glob
import gzip
import os
import re
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

BASE_DIR = hpc.BASE
INPUT_DIR = BASE_DIR / "inputs"
CHECKV_DIR = BASE_DIR / "00_checkv" / "out"
OUT_DIR = BASE_DIR / "01_curated_datasets"
WORK_DIR = OUT_DIR / "meme_fimo_work"

GENOMAD_DIR = hpc.GENOMAD_DIR
PRODIGAL_DIR = hpc.PRODIGAL_DIR
VIRAL_PFAM_DOMTBL = hpc.VIRAL_PFAM_DIR / "sbl_proteins_pfam.domtbl.gz"
VIRAL_FAA = hpc.VIRAL_PFAM_DIR / "sbl_proteins.faa"
STEP1_DIR = hpc.STEP1_DIR
PFAM_CLANS_TSV = hpc.PFAM_CLANS_TSV

MEME_BIN = hpc.MEME_BIN
FIMO_BIN = hpc.FIMO_BIN

FOCAL_FAMILIES = ("PF00144", "PF00905", "PF13354")
AUX_PRE_FAMILIES = ("PF00768", "PF02113")
NON_PRE_CL0013 = ("PF04960", "PF06089", "PF27310")
OV_THRESHOLD_AA = 30
DOMAIN_PAD_AA = 10

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
ALT_START_CODONS = {"GTG", "TTG", "ATT", "ATC", "CTG"}


def open_text(path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8", errors="ignore") if p.endswith(".gz") else open(p, "rt", encoding="utf-8", errors="ignore")


def read_fasta(path):
    records = []
    cur_hdr = None
    buf = []
    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_hdr is not None:
                    records.append((cur_hdr.split()[0], cur_hdr, "".join(buf)))
                cur_hdr = line[1:].strip()
                buf = []
            else:
                buf.append(line)
        if cur_hdr is not None:
            records.append((cur_hdr.split()[0], cur_hdr, "".join(buf)))
    return records


def write_fasta(path, items):
    """items: list of (header, seq)."""
    with open(path, "w", encoding="utf-8") as f:
        for hdr, seq in items:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def read_tsv(path):
    with open_text(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def revcomp(dna):
    return dna.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def translate_dna(dna, alt_start=True, tga_w=False):
    dna = dna.upper()
    n = len(dna) - (len(dna) % 3)
    aa = []
    for i in range(0, n, 3):
        cod = dna[i : i + 3]
        if i == 0 and alt_start and cod in ALT_START_CODONS:
            aa.append("M")
        elif tga_w and cod == "TGA":
            aa.append("W")
        else:
            aa.append(CODON_STD.get(cod, "X"))
    return "".join(aa)


def match_protein_to_dna(prot_seq, dna_seq):
    """Return exact coding DNA slice (len == 3*len(prot_seq)) if translation matches, else None."""
    target = prot_seq.rstrip("*")
    dna = dna_seq.upper()
    # Strip trailing stop codon if present
    candidates = [dna]
    if len(dna) >= 3 and len(dna) % 3 == 0 and dna[-3:] in ("TAA", "TAG", "TGA"):
        candidates.append(dna[:-3])

    for cand in candidates:
        for tga_w in (False, True):
            aa = translate_dna(cand, alt_start=True, tga_w=tga_w).rstrip("*")
            if aa == target and len(cand) >= 3 * len(target):
                return cand[: 3 * len(target)]
            # Also try exact start codon literal (without alt_start -> M)
            aa_lit = translate_dna(cand, alt_start=False, tga_w=tga_w).rstrip("*")
            if aa_lit == target and len(cand) >= 3 * len(target):
                return cand[: 3 * len(target)]
            # Substring match (e.g. different upstream start codon in Prodigal)
            if len(target) >= 30:
                idx = aa_lit.find(target)
                if idx >= 0:
                    return cand[idx * 3 : (idx + len(target)) * 3]
                # Allow N-terminal M mismatch due to GTG/TTG start codon
                if target.startswith("M") and len(target) > 30:
                    idx2 = aa_lit.find(target[1:])
                    if idx2 >= 1:
                        s_nt = cand[(idx2 - 1) * 3 : (idx2 - 1 + len(target)) * 3]
                        if len(s_nt) == 3 * len(target) and s_nt[:3] in (ALT_START_CODONS | {"ATG"}):
                            return s_nt
    return None


def interval_overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def parse_checkv_regions(region_types_str, region_coords_str):
    """Parse CheckV comma-separated region types and 1-based inclusive start-end intervals."""
    viral_iv = []
    host_iv = []
    if not region_types_str or region_types_str == "NA" or not region_coords_str or region_coords_str == "NA":
        return viral_iv, host_iv
    rtypes = [x.strip().lower() for x in region_types_str.split(",")]
    rcoords = [x.strip() for x in region_coords_str.split(",")]
    for rt, rc in zip(rtypes, rcoords):
        m = re.match(r"^(\d+)-(\d+)$", rc)
        if not m:
            continue
        st, en = int(m.group(1)), int(m.group(2))
        if rt == "viral":
            viral_iv.append((st, en))
        elif rt == "host":
            host_iv.append((st, en))
    return viral_iv, host_iv


def load_pfam_clans():
    acc_to_clan = {}
    acc_to_name = {}
    with open_text(PFAM_CLANS_TSV) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                acc, clan, _, short_id = parts[0], parts[1], parts[2], parts[3]
                acc_to_clan[acc] = clan if clan else "NO_CLAN"
                acc_to_name[acc] = short_id
    return acc_to_clan, acc_to_name


def parse_hmmsearch_domtbl(path):
    """Return dict: target_id -> list of domain hit dicts."""
    hits = defaultdict(list)
    with open_text(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 22:
                continue
            target_id = p[0]
            tlen = int(p[2])
            pfam_name = p[3]
            pfam_acc = p[4].split(".")[0]
            dom_score = float(p[13])
            hmm_from = int(p[15])
            hmm_to = int(p[16])
            env_from = int(p[19])
            env_to = int(p[20])
            hits[target_id].append({
                "pfam_acc": pfam_acc,
                "pfam_name": pfam_name,
                "tlen": tlen,
                "dom_score": dom_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "env_from": env_from,
                "env_to": env_to,
            })
    return hits


def resolve_overlapping_domains(hit_list, acc_to_clan, ov_thresh=OV_THRESHOLD_AA):
    """Cluster hits with envelope overlap >= ov_thresh aa and pick highest dom_score per cluster."""
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
        # Winner is highest dom_score in connected overlapping component
        comp_hits = [hit_list[idx] for idx in comp]
        best = max(comp_hits, key=lambda h: (h["dom_score"], h["env_to"] - h["env_from"]))
        comp_accs = sorted({h["pfam_acc"] for h in comp_hits})
        winners.append({
            **best,
            "clan": acc_to_clan.get(best["pfam_acc"], "NO_CLAN"),
            "comp_candidates": ",".join(comp_accs),
            "n_comp_hits": len(comp_hits),
        })
    winners.sort(key=lambda h: h["env_from"])
    return winners


def classify_protein_by_winners(winners):
    """Classify protein from non-overlapping winning domains."""
    if not winners:
        return "drop_no_pfam_ga", "NONE", None
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


def parse_fimo_best(tsv_path, id_map):
    """Map FIMO safe IDs back to original IDs and record best p-value, q-value, score, motifs."""
    out = {}
    with open(tsv_path, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline()
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            motif_id = p[0]
            safe_id = p[2]
            orig_id = id_map.get(safe_id, safe_id)
            try:
                score = float(p[6])
                pval = float(p[7])
                qval = float(p[8]) if len(p) > 8 and p[8] else 1.0
            except ValueError:
                continue
            if orig_id not in out:
                out[orig_id] = {"min_p": pval, "min_q": qval, "max_score": score, "motifs": {motif_id}}
            else:
                st = out[orig_id]
                if pval < st["min_p"]:
                    st["min_p"] = pval
                if qval < st["min_q"]:
                    st["min_q"] = qval
                if score > st["max_score"]:
                    st["max_score"] = score
                st["motifs"].add(motif_id)
    return out


def run_cmd(cmd):
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {cmd}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc.stdout


def run_meme_fimo_family(fam, candidates):
    """Run unified 2-round MEME/FIMO motif filter on focal family candidates.

    candidates: list of dicts with keys:
      'uid', 'origin' ('viral'|'bacterial'), 'dom_seq', 'dom_score'
    Uses domain-extracted core sequences for both MEME Round 2 refinement and FIMO scanning
    so passenger fusion domains cannot produce false motif matches.
    """
    fwork = WORK_DIR / fam
    fwork.mkdir(parents=True, exist_ok=True)
    seed_faa = INPUT_DIR / f"{fam}_seed.faa"

    # Write safe-ID FASTA splits for viral & bacterial domain cores
    safe_to_uid = {}
    v_items, b_items = [], []
    for idx, c in enumerate(candidates):
        sid = f"S_{idx:05d}"
        safe_to_uid[sid] = c["uid"]
        if c["origin"] == "viral":
            v_items.append((sid, c["dom_seq"]))
        else:
            b_items.append((sid, c["dom_seq"]))

    v_fa = fwork / f"{fam}_viral_dom_safe.faa"
    b_fa = fwork / f"{fam}_bac_dom_safe.faa"
    write_fasta(v_fa, v_items)
    write_fasta(b_fa, b_items)

    # Round 1: MEME on official Pfam seed -> FIMO on viral & bacterial identically
    r1_meme = fwork / "r1_meme"
    r1_v_fimo = fwork / "r1_fimo_viral"
    r1_b_fimo = fwork / "r1_fimo_bac"
    run_cmd(f"{MEME_BIN} {seed_faa} -protein -nmotifs 5 -minw 6 -maxw 20 -mod zoops -maxsize 1000000 -nostatus -oc {r1_meme}")
    if v_items:
        run_cmd(f"{FIMO_BIN} --thresh 1e-4 --oc {r1_v_fimo} {r1_meme}/meme.txt {v_fa}")
    if b_items:
        run_cmd(f"{FIMO_BIN} --thresh 1e-4 --oc {r1_b_fimo} {r1_meme}/meme.txt {b_fa}")

    r1_stats = {}
    if v_items:
        r1_stats.update(parse_fimo_best(r1_v_fimo / "fimo.tsv", safe_to_uid))
    if b_items:
        r1_stats.update(parse_fimo_best(r1_b_fimo / "fimo.tsv", safe_to_uid))

    # Build balanced Round-2 training set: up to 125 non-redundant viral + 125 bacterial R1 passers (p<=1e-5)
    def select_top_nr(origin_str, cap=125):
        pool = [c for c in candidates if c["origin"] == origin_str and r1_stats.get(c["uid"], {}).get("min_p", 1.0) <= 1e-5 and 40 <= len(c["dom_seq"]) < 1000]
        pool.sort(key=lambda c: (-c["dom_score"], r1_stats[c["uid"]]["min_p"], c["uid"]))
        chosen = []
        seen_seq = set()
        for c in pool:
            if c["dom_seq"] not in seen_seq:
                seen_seq.add(c["dom_seq"])
                chosen.append(c)
                if len(chosen) >= cap:
                    break
        return chosen

    r2_v = select_top_nr("viral", 125)
    r2_b = select_top_nr("bacterial", 125)
    r2_items = []
    seen_all = set()
    for idx, c in enumerate(r2_v + r2_b):
        if c["dom_seq"] not in seen_all:
            seen_all.add(c["dom_seq"])
            r2_items.append((f"T_{idx:04d}", c["dom_seq"]))
    r2_train_fa = fwork / f"{fam}_r2_train.faa"
    write_fasta(r2_train_fa, r2_items)

    # Round 2: MEME on balanced decontaminated SBL set -> FIMO identically
    r2_meme = fwork / "r2_meme"
    r2_v_fimo = fwork / "r2_fimo_viral"
    r2_b_fimo = fwork / "r2_fimo_bac"
    run_cmd(f"{MEME_BIN} {r2_train_fa} -protein -nmotifs 5 -minw 6 -maxw 20 -mod zoops -maxsize 1000000 -nostatus -oc {r2_meme}")
    if v_items:
        run_cmd(f"{FIMO_BIN} --thresh 1e-4 --oc {r2_v_fimo} {r2_meme}/meme.txt {v_fa}")
    if b_items:
        run_cmd(f"{FIMO_BIN} --thresh 1e-4 --oc {r2_b_fimo} {r2_meme}/meme.txt {b_fa}")

    r2_stats = {}
    if v_items:
        r2_stats.update(parse_fimo_best(r2_v_fimo / "fimo.tsv", safe_to_uid))
    if b_items:
        r2_stats.update(parse_fimo_best(r2_b_fimo / "fimo.tsv", safe_to_uid))

    return r1_stats, r2_stats, len(r2_v), len(r2_b)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    summary_lines = ["=== SBL Phase 1 Re-analysis Curation Summary ==="]

    # ------------------------------------------------------------------
    # STAGE 1: Load geNomad v1.7, CheckV v1.5, Contigs, Prodigal ORFs & ARG_cds
    # ------------------------------------------------------------------
    print("[Stage 1] Loading geNomad v1.7, CheckV v1.5, contigs, and ORF CDSs...")
    vs_rows = read_tsv(GENOMAD_DIR / "sbl_contigs_virus_summary.tsv.gz")
    pv_rows = read_tsv(GENOMAD_DIR / "sbl_contigs_provirus.tsv.gz")
    pl_rows = read_tsv(GENOMAD_DIR / "genomad_out/sbl_contigs_summary/sbl_contigs_plasmid_summary.tsv")
    agg_rows = read_tsv(GENOMAD_DIR / "genomad_out/sbl_contigs_aggregated_classification/sbl_contigs_aggregated_classification.tsv")

    gn_virus_contigs = {r["seq_name"].split("|")[0] for r in vs_rows}
    gn_plasmid_contigs = {r["seq_name"].split("|")[0] for r in pl_rows}
    gn_pv_intervals = defaultdict(list)
    for r in pv_rows:
        gn_pv_intervals[r["source_seq"]].append((int(r["start"]), int(r["end"])))
    gn_agg = {
        r["seq_name"]: (float(r["chromosome_score"]), float(r["plasmid_score"]), float(r["virus_score"]))
        for r in agg_rows
    }

    cv_qual = {r["contig_id"]: r for r in read_tsv(CHECKV_DIR / "quality_summary.tsv")}
    cv_contam = {r["contig_id"]: r for r in read_tsv(CHECKV_DIR / "contamination.tsv")}
    cv_regions = {}
    for cid, r in cv_contam.items():
        if r.get("provirus", "No") == "Yes":
            cv_regions[cid] = parse_checkv_regions(r.get("region_types", ""), r.get("region_coords_bp", ""))

    # Load contig DNA sequences (for fallback coordinate lookup)
    contig_seqs = {}
    for cid, _, dna in read_fasta(PRODIGAL_DIR / "sbl_contigs.fna"):
        contig_seqs[cid] = dna.upper()

    # Load all Prodigal ORFs (.faa headers for coords + .fna for nucleotide sequence)
    prodigal_coords = {}
    with open(PRODIGAL_DIR / "all_orfs.faa", "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                p = line[1:].strip().split(" # ")
                if len(p) >= 4:
                    oid = p[0].strip()
                    prodigal_coords[oid] = (int(p[1]), int(p[2]), int(p[3]))

    prodigal_by_uvig = defaultdict(list)
    for fna_path in sorted(glob.glob(str(PRODIGAL_DIR / "orf" / "*.fna"))):
        for oid, _, dna in read_fasta(fna_path):
            uvig = oid.rsplit("_", 1)[0]
            st, en, strand = prodigal_coords.get(oid, (0, 0, 1))
            dna_u = dna.upper()
            aa_std = translate_dna(dna_u, alt_start=True, tga_w=False).rstrip("*")
            aa_lit = translate_dna(dna_u, alt_start=False, tga_w=False).rstrip("*")
            prodigal_by_uvig[uvig].append((oid, st, en, strand, dna_u, aa_std, aa_lit))

    # Load ARG_cds.fna as verified fallback for non-standard GVMAG / partial ORFs
    arg_cds_by_id = {}
    for sid, _, dna in read_fasta(INPUT_DIR / "ARG_cds.fna"):
        arg_cds_by_id[sid] = dna.upper()

    # Read 5,701 viral proteins
    viral_sbl_raw = read_fasta(VIRAL_FAA)
    print(f"Loaded {len(viral_sbl_raw)} viral SBL proteins across {len(contig_seqs)} contigs.")

    # Map each viral protein to (st, en, strand, exact_cds_nt, cds_source)
    viral_loci = {}
    cds_source_counts = Counter()
    for uid, hdr, prot in viral_sbl_raw:
        prot_clean = prot.rstrip("*")
        uvig = uid.split("|")[0]
        orig_fam = "UNK"
        for tok in hdr.split():
            if tok.startswith("fam="):
                orig_fam = tok.split("=", 1)[1]
        matched_cds = None
        matched_coord = None
        cds_src = "unmatched"

        # 1. Exact/substring match in Prodigal ORFs on same contig
        for oid, st, en, strand, dna_u, aa_std, aa_lit in prodigal_by_uvig.get(uvig, []):
            if aa_std == prot_clean or aa_lit == prot_clean:
                cslice = match_protein_to_dna(prot_clean, dna_u)
                if cslice is not None:
                    matched_cds = cslice
                    matched_coord = (st, en, strand)
                    cds_src = "prodigal_exact"
                    break
        if matched_cds is None:
            for oid, st, en, strand, dna_u, aa_std, aa_lit in prodigal_by_uvig.get(uvig, []):
                if (len(prot_clean) >= 30 and (prot_clean in aa_lit or prot_clean[1:] in aa_lit)) or (
                    len(aa_lit) >= 30 and aa_lit in prot_clean
                ):
                    cslice = match_protein_to_dna(prot_clean, dna_u)
                    if cslice is not None:
                        matched_cds = cslice
                        matched_coord = (st, en, strand)
                        cds_src = "prodigal_substr"
                        break
                    elif matched_coord is None:
                        matched_coord = (st, en, strand)

        # 2. Fallback to ARG_cds.fna with strict 1:1 translation check + exact contig position
        if matched_cds is None and uid in arg_cds_by_id:
            cand_dna = arg_cds_by_id[uid]
            cslice = match_protein_to_dna(prot_clean, cand_dna)
            if cslice is not None:
                matched_cds = cslice
                cds_src = "arg_cds_verified"
                # Locate exact coordinates on contig_seqs[uvig]
                if matched_coord is None and uvig in contig_seqs:
                    cdna = contig_seqs[uvig]
                    idx_f = cdna.find(cslice[: min(60, len(cslice))])
                    if idx_f >= 0:
                        matched_coord = (idx_f + 1, idx_f + len(cslice), 1)
                    else:
                        rc = revcomp(cslice)
                        idx_r = cdna.find(rc[: min(60, len(rc))])
                        if idx_r >= 0:
                            matched_coord = (idx_r + 1, idx_r + len(cslice), -1)

        cds_source_counts[cds_src] += 1
        viral_loci[uid] = {
            "uid": uid,
            "uvig": uvig,
            "orig_fam": orig_fam,
            "prot": prot_clean,
            "coord": matched_coord,
            "cds": matched_cds,
            "cds_src": cds_src,
        }

    print("Viral ORF-to-CDS mapping breakdown:", dict(cds_source_counts))
    summary_lines.append(f"Total input viral SBL proteins: {len(viral_loci)}")
    summary_lines.append(f"Viral CDS mapping: {dict(cds_source_counts)}")

    # ------------------------------------------------------------------
    # STAGE 2: Scheme B Joint Viral Decontamination (geNomad v1.7 + CheckV v1.5)
    # ------------------------------------------------------------------
    print("[Stage 2] Running Scheme B viral decontamination (geNomad v1.7 + CheckV v1.5)...")
    decontam_counts = Counter()
    kept_viral = {}
    viral_decontam_audit = []

    for uid, info in viral_loci.items():
        uvig = info["uvig"]
        coord = info["coord"]
        in_gn_virus = uvig in gn_virus_contigs
        in_gn_plasmid = uvig in gn_plasmid_contigs
        cs, ps, vs_ = gn_agg.get(uvig, (0.0, 0.0, 0.0))
        qrow = cv_qual.get(uvig, {})
        cv_vg = int(float(qrow.get("viral_genes", 0) or 0))
        cv_hg = int(float(qrow.get("host_genes", 0) or 0))
        cv_contam_pct = float(qrow.get("contamination", 0.0) or 0.0)
        cv_pv = qrow.get("provirus", "No")

        status = "KEEP_VIRAL"
        gn_pv_ov_bp = 0
        cv_v_ov_bp = 0
        cv_h_ov_bp = 0

        if in_gn_plasmid:
            status = "DROP_genomad_plasmid"
        else:
            # Check geNomad provirus coordinate overlap
            if uvig in gn_pv_intervals:
                if coord is not None:
                    st, en, _ = coord
                    gn_pv_ov_bp = max(interval_overlap(st, en, ps0, pe0) for ps0, pe0 in gn_pv_intervals[uvig])
                    if gn_pv_ov_bp == 0:
                        status = "DROP_genomad_host_flank"
            # Check CheckV v1.5 provirus viral vs host coordinate overlap
            if status == "KEEP_VIRAL" and uvig in cv_regions and coord is not None:
                st, en, _ = coord
                v_ivs, h_ivs = cv_regions[uvig]
                cv_v_ov_bp = sum(interval_overlap(st, en, vs0, ve0) for vs0, ve0 in v_ivs)
                cv_h_ov_bp = sum(interval_overlap(st, en, hs0, he0) for hs0, he0 in h_ivs)
                if cv_v_ov_bp == 0 or cv_h_ov_bp > cv_v_ov_bp:
                    status = "DROP_checkv_host_flank"

            # Check contigs NOT called virus by geNomad v1.7
            if status == "KEEP_VIRAL" and not in_gn_virus:
                if vs_ >= cs and vs_ >= ps:
                    # Scheme B rescue rule: CheckV viral_genes >= 1 and host_genes == 0
                    if cv_vg >= 1 and cv_hg == 0:
                        status = "KEEP_Tier3_CheckV_Rescued"
                    else:
                        status = "DROP_nonviral_unverified_fragment"
                else:
                    status = "DROP_nonviral_chrom_or_plasmid"

        decontam_counts[status] += 1
        if status.startswith("KEEP"):
            kept_viral[uid] = info

        viral_decontam_audit.append({
            "sbl_protein_id": uid,
            "uvig": uvig,
            "orig_fam": info["orig_fam"],
            "cds_match": info["cds_src"],
            "orf_start": coord[0] if coord else "NA",
            "orf_end": coord[1] if coord else "NA",
            "genomad_in_virus_summary": in_gn_virus,
            "genomad_in_plasmid_summary": in_gn_plasmid,
            "genomad_chrom_score": cs,
            "genomad_plasmid_score": ps,
            "genomad_virus_score": vs_,
            "genomad_pv_overlap_bp": gn_pv_ov_bp,
            "checkv_provirus": cv_pv,
            "checkv_viral_genes": cv_vg,
            "checkv_host_genes": cv_hg,
            "checkv_contamination_pct": cv_contam_pct,
            "checkv_viral_overlap_bp": cv_v_ov_bp,
            "checkv_host_overlap_bp": cv_h_ov_bp,
            "schemeB_status": status,
        })

    print("Scheme B Viral Decontamination Breakdown:", dict(decontam_counts))
    summary_lines.append("Viral Decontamination (Scheme B):")
    for k, v in sorted(decontam_counts.items()):
        summary_lines.append(f"  {k}: {v}")
    summary_lines.append(f"Total viral proteins surviving Scheme B: {len(kept_viral)}")

    with open(OUT_DIR / "Phase1_Viral_Decontamination_Audit.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(viral_decontam_audit[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(viral_decontam_audit)

    # ------------------------------------------------------------------
    # STAGE 3: Load Bacterial CleanCore (1,775) + SARG Plasmid Metadata
    # ------------------------------------------------------------------
    print("[Stage 3] Loading Bacterial CleanCore & Plasmid/MGE tags...")
    bac_audit = read_tsv(INPUT_DIR / "PerSequence_Bacterial_Audit.tsv")
    core_rows = [r for r in bac_audit if r["pass_CleanCore_TierB"] == "True"]

    # Map full headers and protein sequences from PF*_SARG.faa
    sarg_seqs = {}
    sarg_hdrs = {}
    for fam in FOCAL_FAMILIES:
        for sid, hdr, seq in read_fasta(INPUT_DIR / f"{fam}_SARG.faa"):
            sarg_seqs[(fam, sid)] = seq.rstrip("*")
            sarg_hdrs[(fam, sid)] = hdr

    # Map (orig_family, seq_id) -> safe_id in step1_pfam.domtbl.gz
    idmap_lookup = {}
    with open(STEP1_DIR / "step1_idmap.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["origin"] != "virus":
                idmap_lookup[(r["assigned_family"], r["orig_id"])] = r["safe_id"]

    # ------------------------------------------------------------------
    # STAGE 4: Overlapping-Envelope Full Pfam-A (--cut_ga) Competition
    # ------------------------------------------------------------------
    print("[Stage 4] Running overlapping-envelope Pfam-A (--cut_ga) competition...")
    acc_to_clan, acc_to_name = load_pfam_clans()
    viral_hits = parse_hmmsearch_domtbl(VIRAL_PFAM_DOMTBL)
    step1_hits = parse_hmmsearch_domtbl(STEP1_DIR / "step1_pfam.domtbl.gz")

    pfam_class_counts = Counter()
    arch_rows = []
    focal_pool = defaultdict(list)
    aux_pool = defaultdict(list)

    # Process surviving viral proteins
    for uid, info in kept_viral.items():
        winners = resolve_overlapping_domains(viral_hits.get(uid, []), acc_to_clan, OV_THRESHOLD_AA)
        cat, win_acc, best_dom = classify_protein_by_winners(winners)
        pfam_class_counts[("viral", cat)] += 1
        arch_str = "+".join(f"{w['pfam_acc']}:{w['pfam_name']}({w['env_from']}-{w['env_to']})" for w in winners)
        cl13_cands = ";".join(w["comp_candidates"] for w in winners if w["clan"] == "CL0013")
        prot = info["prot"]

        dom_start = max(1, best_dom["env_from"] - DOMAIN_PAD_AA) if best_dom else 1
        dom_end = min(len(prot), best_dom["env_to"] + DOMAIN_PAD_AA) if best_dom else len(prot)
        dom_seq = prot[dom_start - 1 : dom_end]

        dom_cds = None
        if info["cds"] is not None and len(info["cds"]) >= 3 * dom_end:
            dom_cds = info["cds"][(dom_start - 1) * 3 : dom_end * 3]

        arch_rows.append({
            "uid": uid,
            "origin": "viral",
            "orig_family": info["orig_fam"],
            "pfam_category": cat,
            "winning_acc": win_acc,
            "seq_len": len(prot),
            "n_winning_domains": len(winners),
            "domain_architecture": arch_str,
            "cl0013_competitors": cl13_cands,
            "win_dom_score": best_dom["dom_score"] if best_dom else 0.0,
            "env_from": best_dom["env_from"] if best_dom else 0,
            "env_to": best_dom["env_to"] if best_dom else 0,
            "dom_slice_start": dom_start,
            "dom_slice_end": dom_end,
            "hmm_from": best_dom["hmm_from"] if best_dom else 0,
            "hmm_to": best_dom["hmm_to"] if best_dom else 0,
        })

        rec = {
            "uid": uid,
            "origin": "viral",
            "orig_family": info["orig_fam"],
            "assigned_family": win_acc,
            "prot": prot,
            "dom_seq": dom_seq,
            "dom_start": dom_start,
            "dom_end": dom_end,
            "dom_score": best_dom["dom_score"] if best_dom else 0.0,
            "env_from": best_dom["env_from"] if best_dom else 0,
            "env_to": best_dom["env_to"] if best_dom else 0,
            "cds": info["cds"],
            "dom_cds": dom_cds,
            "cds_src": info["cds_src"],
        }
        if cat.startswith("focal_"):
            focal_pool[win_acc].append(rec)
        elif cat.startswith("aux_"):
            aux_pool[win_acc].append(rec)

    # Process bacterial CleanCore proteins
    bac_meta_rows = []
    seen_bac_uids = set()
    for r in core_rows:
        ofam, sid = r["family"], r["seq_id"]
        prot = sarg_seqs[(ofam, sid)]
        hdr = sarg_hdrs[(ofam, sid)]
        safe_id = idmap_lookup[(ofam, sid)]
        # Disambiguate if same seq_id appeared across SARG files
        uid = sid if sid not in seen_bac_uids else f"{ofam}|{sid}"
        seen_bac_uids.add(uid)

        winners = resolve_overlapping_domains(step1_hits.get(safe_id, []), acc_to_clan, OV_THRESHOLD_AA)
        cat, win_acc, best_dom = classify_protein_by_winners(winners)
        pfam_class_counts[("bacterial", cat)] += 1
        arch_str = "+".join(f"{w['pfam_acc']}:{w['pfam_name']}({w['env_from']}-{w['env_to']})" for w in winners)
        cl13_cands = ";".join(w["comp_candidates"] for w in winners if w["clan"] == "CL0013")

        dom_start = max(1, best_dom["env_from"] - DOMAIN_PAD_AA) if best_dom else 1
        dom_end = min(len(prot), best_dom["env_to"] + DOMAIN_PAD_AA) if best_dom else len(prot)
        dom_seq = prot[dom_start - 1 : dom_end]

        is_pl = bool(re.search(r"plasmid", hdr, re.I))
        is_mge = bool(re.search(r"plasmid|transposon|integron|conjugat|mobile|insertion sequence", hdr, re.I))
        is_clin_bla = bool(re.search(r"\b(TEM|SHV|CTX-M|OXA|KPC|PER|VEB|GES|CMY|DHA|FOX|ACC|MIR|ACT|SME|IMI|NMC|BEL|TLA|CARB|PSE|ROB|BRO)[-_0-9]", hdr, re.I))

        bac_meta_rows.append({
            "uid": uid,
            "seq_id": sid,
            "orig_family": ofam,
            "pfam_category": cat,
            "reassigned_family": win_acc,
            "source_type": r["source_type"],
            "sarg_phylum": r["sarg_phylum"],
            "sarg_species": r["sarg_species"],
            "is_plasmid_header": is_pl,
            "is_mge_header": is_mge,
            "is_clinical_mobile_bla": is_clin_bla,
            "win_dom_score": best_dom["dom_score"] if best_dom else 0.0,
            "domain_architecture": arch_str,
            "full_header": hdr,
        })

        arch_rows.append({
            "uid": uid,
            "origin": "bacterial",
            "orig_family": ofam,
            "pfam_category": cat,
            "winning_acc": win_acc,
            "seq_len": len(prot),
            "n_winning_domains": len(winners),
            "domain_architecture": arch_str,
            "cl0013_competitors": cl13_cands,
            "win_dom_score": best_dom["dom_score"] if best_dom else 0.0,
            "env_from": best_dom["env_from"] if best_dom else 0,
            "env_to": best_dom["env_to"] if best_dom else 0,
            "dom_slice_start": dom_start,
            "dom_slice_end": dom_end,
            "hmm_from": best_dom["hmm_from"] if best_dom else 0,
            "hmm_to": best_dom["hmm_to"] if best_dom else 0,
        })

        rec = {
            "uid": uid,
            "origin": "bacterial",
            "orig_family": ofam,
            "assigned_family": win_acc,
            "prot": prot,
            "dom_seq": dom_seq,
            "dom_start": dom_start,
            "dom_end": dom_end,
            "dom_score": best_dom["dom_score"] if best_dom else 0.0,
            "env_from": best_dom["env_from"] if best_dom else 0,
            "env_to": best_dom["env_to"] if best_dom else 0,
            "cds": None,
            "dom_cds": None,
            "cds_src": "deferred_ncbi_ipg",
        }
        if cat.startswith("focal_"):
            focal_pool[win_acc].append(rec)
        elif cat.startswith("aux_"):
            aux_pool[win_acc].append(rec)

    print("Pfam-A (--cut_ga) envelope competition summary:")
    summary_lines.append("\nPfam-A (--cut_ga) Overlapping-Envelope Competition:")
    for (orig, cat), cnt in sorted(pfam_class_counts.items()):
        print(f"  {orig:10s} {cat}: {cnt}")
        summary_lines.append(f"  {orig:10s} {cat}: {cnt}")

    with open(OUT_DIR / "Phase1_MultiDomain_Architecture_Catalog.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(arch_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(arch_rows)

    with open(OUT_DIR / "Phase1_Bacterial_Metadata.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(bac_meta_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(bac_meta_rows)

    # ------------------------------------------------------------------
    # STAGE 5: Unified Two-Round MEME/FIMO 5.5.9 Filtering (Parallel across 3 families)
    # ------------------------------------------------------------------
    print("[Stage 5] Running unified 2-round MEME/FIMO 5.5.9 across PF00144, PF00905, PF13354...")
    motif_results = {}
    with ThreadPoolExecutor(max_workers=3) as ex:
        futs = {fam: ex.submit(run_meme_fimo_family, fam, focal_pool[fam]) for fam in FOCAL_FAMILIES}
        for fam, fut in futs.items():
            motif_results[fam] = fut.result()

    motif_audit_rows = []
    summary_lines.append("\nUnified Two-Round MEME/FIMO (p <= 1e-5, no SXXK regex pre-filter) & Final Clean Counts:")

    for fam in FOCAL_FAMILIES:
        r1_stats, r2_stats, n_r2v, n_r2b = motif_results[fam]
        v_clean_full, b_clean_full = [], []
        v_clean_dom, b_clean_dom = [], []
        v_clean_cds, v_clean_dom_cds = [], []
        cnts = Counter()

        for c in focal_pool[fam]:
            uid = c["uid"]
            org = c["origin"]
            p1 = r1_stats.get(uid, {}).get("min_p", 1.0)
            p2 = r2_stats.get(uid, {}).get("min_p", 1.0)
            n_m2 = len(r2_stats.get(uid, {}).get("motifs", set()))
            has_sxxk = bool(re.search(r"S[A-Z]{2}K", c["dom_seq"]))
            pass_r1 = p1 <= 1e-5
            pass_r2 = p2 <= 1e-5
            pass_both_1e5 = pass_r1 and pass_r2
            pass_either_1e5 = pass_r1 or pass_r2
            pass_both_1e10 = (p1 <= 1e-10) and (p2 <= 1e-10)

            cnts[(org, "post_pfam")] += 1
            if pass_r1:
                cnts[(org, "pass_r1_1e5")] += 1
            if pass_r2:
                cnts[(org, "pass_r2_1e5")] += 1
            if pass_both_1e5:
                cnts[(org, "pass_both_1e5")] += 1
            if pass_either_1e5:
                cnts[(org, "pass_either_1e5")] += 1
            if pass_both_1e10:
                cnts[(org, "pass_both_1e10")] += 1
            if pass_both_1e5 and has_sxxk:
                cnts[(org, "pass_both_1e5_with_SXXK")] += 1

            motif_audit_rows.append({
                "family": fam,
                "origin": org,
                "uid": uid,
                "orig_family": c["orig_family"],
                "dom_score": c["dom_score"],
                "dom_len": len(c["dom_seq"]),
                "r1_min_p": p1,
                "r2_min_p": p2,
                "r2_n_motifs": n_m2,
                "has_SXXK_in_dom": has_sxxk,
                "pass_r1_1e5": pass_r1,
                "pass_r2_1e5": pass_r2,
                "pass_both_1e5": pass_both_1e5,
                "pass_either_1e5": pass_either_1e5,
                "pass_both_1e10": pass_both_1e10,
                "has_verified_cds": c["cds"] is not None,
            })

            if pass_both_1e5:
                hdr_meta = f"{uid} fam={fam} origin={org} dom={c['dom_start']}..{c['dom_end']} score={c['dom_score']:.1f}"
                if org == "viral":
                    v_clean_full.append((hdr_meta, c["prot"]))
                    v_clean_dom.append((hdr_meta, c["dom_seq"]))
                    if c["cds"] is not None:
                        v_clean_cds.append((hdr_meta, c["cds"]))
                        cnts[("viral", "final_verified_cds")] += 1
                    if c["dom_cds"] is not None:
                        v_clean_dom_cds.append((hdr_meta, c["dom_cds"]))
                        cnts[("viral", "final_verified_dom_cds")] += 1
                else:
                    b_clean_full.append((hdr_meta, c["prot"]))
                    b_clean_dom.append((hdr_meta, c["dom_seq"]))

        # Export FASTAs for this focal family
        write_fasta(OUT_DIR / f"{fam}_viral_clean.faa", v_clean_full)
        write_fasta(OUT_DIR / f"{fam}_bacterial_clean.faa", b_clean_full)
        write_fasta(OUT_DIR / f"{fam}_combined_clean.faa", v_clean_full + b_clean_full)
        write_fasta(OUT_DIR / f"{fam}_viral_dom.faa", v_clean_dom)
        write_fasta(OUT_DIR / f"{fam}_bacterial_dom.faa", b_clean_dom)
        write_fasta(OUT_DIR / f"{fam}_combined_dom.faa", v_clean_dom + b_clean_dom)
        write_fasta(OUT_DIR / f"{fam}_viral_cds.fna", v_clean_cds)
        write_fasta(OUT_DIR / f"{fam}_viral_dom_cds.fna", v_clean_dom_cds)

        line = (
            f"  {fam} (R2 train: {n_r2v}V + {n_r2b}B): "
            f"Viral {cnts[('viral','post_pfam')]} -> both_1e5={cnts[('viral','pass_both_1e5')]} "
            f"(either_1e5={cnts[('viral','pass_either_1e5')]}, both_1e10={cnts[('viral','pass_both_1e10')]}, "
            f"SXXK={cnts[('viral','pass_both_1e5_with_SXXK')]}, CDS={cnts[('viral','final_verified_cds')]}/{cnts[('viral','pass_both_1e5')]}) | "
            f"Bacterial {cnts[('bacterial','post_pfam')]} -> both_1e5={cnts[('bacterial','pass_both_1e5')]} "
            f"(either_1e5={cnts[('bacterial','pass_either_1e5')]}, both_1e10={cnts[('bacterial','pass_both_1e10')]})"
        )
        print(line)
        summary_lines.append(line)

    with open(OUT_DIR / "Phase1_Motif_Filter_Audit.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(motif_audit_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(motif_audit_rows)

    # ------------------------------------------------------------------
    # STAGE 6: Auxiliary Clan CL0013 Exports (PF00768, PF02113)
    # ------------------------------------------------------------------
    all_aux_v, all_aux_b = [], []
    summary_lines.append("\nAuxiliary Clan CL0013 Penicillin-Recognizing Enzyme Exports:")
    for acc in AUX_PRE_FAMILIES:
        recs = aux_pool.get(acc, [])
        v_full = [(f"{r['uid']} fam={acc} orig={r['orig_family']} origin=viral dom={r['dom_start']}..{r['dom_end']}", r["prot"]) for r in recs if r["origin"] == "viral"]
        v_dom = [(f"{r['uid']} fam={acc} orig={r['orig_family']} origin=viral dom={r['dom_start']}..{r['dom_end']}", r["dom_seq"]) for r in recs if r["origin"] == "viral"]
        b_full = [(f"{r['uid']} fam={acc} orig={r['orig_family']} origin=bacterial dom={r['dom_start']}..{r['dom_end']}", r["prot"]) for r in recs if r["origin"] == "bacterial"]
        b_dom = [(f"{r['uid']} fam={acc} orig={r['orig_family']} origin=bacterial dom={r['dom_start']}..{r['dom_end']}", r["dom_seq"]) for r in recs if r["origin"] == "bacterial"]
        write_fasta(OUT_DIR / f"CL0013_aux_{acc}_viral.faa", v_full)
        write_fasta(OUT_DIR / f"CL0013_aux_{acc}_viral_dom.faa", v_dom)
        if b_full:
            write_fasta(OUT_DIR / f"CL0013_aux_{acc}_bacterial.faa", b_full)
            write_fasta(OUT_DIR / f"CL0013_aux_{acc}_bacterial_dom.faa", b_dom)
        all_aux_v.extend(v_full)
        all_aux_b.extend(b_full)
        sline = f"  {acc}: Viral={len(v_full)}, Bacterial={len(b_full)}"
        print(sline)
        summary_lines.append(sline)

    write_fasta(OUT_DIR / "CL0013_aux_all_viral.faa", all_aux_v)
    if all_aux_b:
        write_fasta(OUT_DIR / "CL0013_aux_all_bacterial.faa", all_aux_b)

    summary_lines.append(
        "\nNote: Bacterial CDS nucleotide FASTAs were intentionally deferred from Phase 1 because legacy "
        "local PF*_cds.fna files collapsed alleles to 87-170 representative BLAST hits (2.9-5.0% exact translation match). "
        "Before Phase 5 HyPhy, exact 1:1 translation-verified bacterial CDSs will be retrieved via NCBI Entrez/IPG for the CleanCore/MMseqs2 representatives."
    )

    with open(OUT_DIR / "Phase1_Summary_Report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")
    print("[Done] Phase 1 curation completed cleanly.")


if __name__ == "__main__":
    main()
