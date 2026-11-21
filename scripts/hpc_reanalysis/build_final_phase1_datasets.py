#!/usr/bin/env python3
"""Build Revised Phase 1 + Phase 1B Curated Viral & Triple-Source Bacterial SBL/PRE Datasets.

Revision 2026 Fixes Implemented:
1. Strict Pfam Clan CL0013 Highest-Bitscore Assignment (top_overall_cl0013_acc):
   No artificial focal-family override inside PF13354/PF00768 dual hits; assigns strictly
   to the CL0013 family with highest domain bitscore while logging dual_PF13354_PF00768_GA.
2. Strict Own-Family Official-Seed MEME/FIMO 5.5.9 (p <= 1e-5):
   Each sequence in {PF00144, PF00905, PF13354, PF00768, PF02113} is filtered ONLY by FIMO
   using its own assigned family's official Pfam seed MEME model across domain cores (env +/- 10 aa).
3. Retain & Tag Bacterial Prophage Sequences (in_bacterial_prophage):
   Retains loci inside geNomad provirus intervals with in_bacterial_prophage=True /
   replicon=prophage (included in Primary, pruned in the *_nopro sensitivity set). Whole contigs
   called viral inside a GTDB assembly or Marine MAG (likely free phage genomes) are tagged
   replicon=viral_contig / is_viral_contig_in_genome=True, kept in the catalog, and excluded
   from all bacterial tree/HyPhy sets. Phase 4 must classify V->B events whose bacterial
   recipients are all prophage-borne separately (docs/METHODS_AND_PIPELINE.md, Section 4.1).
4. HMM Coverage Filter (hmm_cov >= 0.65 Primary; >=0.60 and >=0.70 tracked):
   Eliminates fragment sequences likely to distort phylogenetic placement.
5. Independent Domain-Slice MMseqs2 nr85 Clustering (--min-seq-id 0.85 -c 0.80 --cov-mode 0):
   Clusters viral domain slices separately from bacteria (prioritizing Phase 2 Strict Tier 1),
   records viral_cluster_size (for reporting only, never leaf-weighting), and exports joint
   PF13354_PF00768 datasets plus 6 bacterial tree replicates (rep1..rep5 + phylo_rep) and
   all 3 non-PRE CL0013 biological outgroups (PF04960, PF06089, PF27310).
6. SBL vs. PBP/PRE Functional Stratification Layer:
   Profiles every sequence against CARD beta-lactamases with PBP-type entries removed
   (mecA/mecC, PBP1a/2b/2x/3 variants; see work/card_excluded_PBP_entries.txt) and
   hmmalign (A2M, match-state verified against HMM LENG) catalytic motifs (SXXK,
   SDN/YAN/YGN/SXN, KTG, Class A Omega-loop E166-x-x-L-N) to separate
   Putative_Hydrolytic_BetaLactamase from PBP_or_PRE_Homolog.
7. Own-seed MEME models for PF00768/PF02113 are built from INPUTS/<fam>_seed.faa when absent,
   and viral-vs-bacterial FIMO pass rates are compared per family (bias check in summary).

Run order (run_phase1_phase2_rebuild.sh):
   --stage curate  ->  build_phase2_neighborhoods.py  ->  --stage finalize
"""

import argparse
import csv
import glob
import gzip
import os
import random
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
P2_TSV = BASE / "02_genomic_context" / "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"
PRODIGAL_DIR = hpc.PRODIGAL_DIR
VIRAL_FAA = hpc.VIRAL_PFAM_DIR / "sbl_proteins.faa"
VIRAL_DOMTBL = hpc.VIRAL_PFAM_DIR / "sbl_proteins_pfam.domtbl.gz"
STEP1_DIR = hpc.STEP1_DIR
PFAM_CLANS = hpc.PFAM_CLANS_TSV
CARD_FAA = INPUTS / "card_db" / "card_true_sbl_refs.faa"

GTDB_DIR = BASE / "01b_bacterial_ref_db"
MAR_DIR = GTDB_DIR / "marine_mags"

FIMO_BIN = hpc.FIMO_BIN
MEME_BIN = hpc.MEME_BIN
MMSEQS_BIN = hpc.MMSEQS_BIN
HMMSEARCH_BIN = hpc.HMMSEARCH_BIN
HMMALIGN_BIN = os.path.join(os.path.dirname(HMMSEARCH_BIN), "hmmalign") if os.path.dirname(HMMSEARCH_BIN) else "hmmalign"
DIAMOND_BIN = "diamond"

FOCAL = ("PF00144", "PF00905", "PF13354")
AUX_PRE = ("PF00768", "PF02113")
PRE_FAMS = ("PF00144", "PF00905", "PF13354", "PF00768", "PF02113")
EXPORT_FAMS = ("PF00144", "PF00905", "PF13354", "PF00768")
OUTGROUP = ("PF04960", "PF06089", "PF27310")
# PF13354 and PF00768 form a dual-hit continuum and may not be reciprocally monophyletic,
# so neither may serve as the rooting outgroup for the other (or for their joint tree).
SISTER_CONTINUUM = {"PF13354", "PF00768"}

# CARD references describing penicillin-binding proteins (mecA/mecC PBP2a, PBP1a/2b/2x/3
# resistance variants) are NOT hydrolytic beta-lactamases and are removed before DIAMOND.
CARD_PBP_RE = re.compile(r"(penicillin[-_ ]?binding|\bPBP\d*[a-z]?\b|\bpbp\d*[a-z]?\b|\bmec[ACR]\d*\b|\bmecA\b|\bmecC\b|ftsI)", re.I)
# Class A Omega-loop signature E166-x-x-L169-N170 (TEM numbering; L may be I/V/M in some enzymes)
OMEGA_E166_RE = re.compile(r"E[A-Z]{2}[LIVM]N")
OMEGA_WINDOW_PAD = 3   # tolerate +/-3 match states around the calibrated Omega window
# Flag own-seed FIMO filters that remove viral sequences much more often than bacterial ones
FIMO_BIAS_WARN_PP = 10.0

OV_THRESH = 30
DOM_PAD = 10
PRIMARY_COV = 0.65

CL0013_HMM_LENS = {
    "PF00144": 329,
    "PF00905": 300,
    "PF13354": 256,
    "PF00768": 239,
    "PF02113": 445,
    "PF04960": 286,
    "PF06089": 322,
    "PF27310": 188,
}

# Calibrated 1-based HMM match-state windows from CARD Class A/C/D prototypes
MOTIF_WINDOWS = {
    "PF00144": {"m1": (49, 60), "m2": (139, 150), "omega_e": (175, 186), "m3": (271, 291)},
    "PF13354": {"m1": (37, 48), "m2": (98, 109), "omega_e": (129, 145), "m3": (199, 211)},
    "PF00905": {"m1": (39, 50), "m2": (123, 134), "omega_e": None, "m3": (233, 244)},
}

TIER_PRIORITY = {
    "Tier1_HighConfidence_Viral": 0,
    "MOSAiC_Proteome_Verified_vMAG": 1,
    "Tier2_Supported_Viral": 2,
    "Tier3_Fragment_or_Borderline": 3,
}

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
STOP_CODONS = {"TAA", "TAG", "TGA"}
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


def revcomp(dna):
    return dna.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


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


def load_hmm_lengths():
    lens = dict(CL0013_HMM_LENS)
    for acc in list(lens):
        p = INPUTS / f"{acc}_official.hmm"
        if p.exists():
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith("LENG "):
                        lens[acc] = int(line.split()[1])
                        break
    return lens


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
                "qlen": int(p[5]),
                "dom_score": float(p[13]),
                "hmm_from": int(p[15]),
                "hmm_to": int(p[16]),
                "env_from": int(p[19]),
                "env_to": int(p[20]),
            })
    return hits


def resolve_strict_cl0013(hit_list, acc_to_clan, hmm_lens):
    """Overlapping-envelope Pfam competition strictly choosing the highest CL0013 bitscore."""
    n = len(hit_list)
    if n == 0:
        return "drop_no_pfam_ga", "NONE", None, False, 0.0, 0.0, 0.0, ""
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
        return "drop_non_CL0013_only", top_non["pfam_acc"], top_non, False, 0.0, 0.0, 0.0, arch_str

    all_cl_accs = set()
    pf13354_sc = 0.0
    pf00768_sc = 0.0
    for _, ch in cl13_clusters:
        for h in ch:
            acc = h["pfam_acc"]
            if acc_to_clan.get(acc, "NO_CLAN") == "CL0013":
                all_cl_accs.add(acc)
                if acc == "PF13354" and h["dom_score"] > pf13354_sc:
                    pf13354_sc = h["dom_score"]
                elif acc == "PF00768" and h["dom_score"] > pf00768_sc:
                    pf00768_sc = h["dom_score"]

    dual_flag = ("PF13354" in all_cl_accs) and ("PF00768" in all_cl_accs)
    chosen = max((top for top, _ in cl13_clusters), key=lambda h: (h["dom_score"], h["env_to"] - h["env_from"]))
    acc = chosen["pfam_acc"]
    hlen = hmm_lens.get(acc, chosen["qlen"])
    hmm_cov = round(min(1.0, (chosen["hmm_to"] - chosen["hmm_from"] + 1) / max(hlen, 1)), 4)

    if acc in FOCAL:
        cat = f"focal_{acc}"
    elif acc in AUX_PRE:
        cat = f"aux_{acc}"
    elif acc in OUTGROUP:
        cat = f"outgroup_{acc}"
    else:
        cat = f"other_CL0013_{acc}"
    return cat, acc, chosen, dual_flag, pf13354_sc, pf00768_sc, hmm_cov, arch_str


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


def ensure_meme_model(fam):
    """Return the round-1 MEME model for `fam`, building it from the official Pfam seed if absent.

    run_phase1.py only built models for the three focal families; PF00768 and PF02113 are built
    here with identical MEME parameters from INPUTS/<fam>_seed.faa (unaligned official seed).
    """
    meme_dir = P1_DIR / "meme_fimo_work" / fam / "r1_meme"
    meme_txt = meme_dir / "meme.txt"
    if meme_txt.exists():
        return meme_txt
    seed_faa = INPUTS / f"{fam}_seed.faa"
    if not seed_faa.exists():
        raise FileNotFoundError(
            f"No MEME model for {fam} and no seed to build one: upload data/sequences/{fam}_seed.faa to {seed_faa}"
        )
    n_seed = sum(1 for _ in read_fasta(seed_faa))
    print(f"  Building MEME model for {fam} from {n_seed} official seed sequences...")
    if n_seed < 20:
        print(f"  [WARN] {fam} seed has only {n_seed} sequences; motifs may be over-specific (see bias check).")
    meme_dir.parent.mkdir(parents=True, exist_ok=True)
    run_cmd(f"{MEME_BIN} {seed_faa} -protein -nmotifs 5 -minw 6 -maxw 20 -mod zoops -maxsize 1000000 -nostatus -oc {meme_dir}")
    return meme_txt


def extract_domain_slice_and_build_record(uid, orig_uid, prot, dom, dual, pf13354_sc, pf00768_sc, hmm_cov, arch, extra):
    ds = max(1, dom["env_from"] - DOM_PAD)
    de = min(len(prot), dom["env_to"] + DOM_PAD)
    dom_seq = prot[ds - 1 : de]
    cds = extra.pop("cds", None)
    dom_cds = cds[(ds - 1) * 3 : de * 3] if (cds and len(cds) >= 3 * de) else None
    return {
        "uid": uid,
        "orig_uid": orig_uid,
        "winning_acc": dom["pfam_acc"],
        "dual_PF13354_PF00768_GA": dual,
        "pf13354_dom_score": pf13354_sc,
        "pf00768_dom_score": pf00768_sc,
        "dom_score": dom["dom_score"],
        "hmm_cov": hmm_cov,
        "pass_cov60": hmm_cov >= 0.60,
        "pass_cov65": hmm_cov >= 0.65,
        "pass_cov70": hmm_cov >= 0.70,
        "seq_len": len(prot),
        "env_from": dom["env_from"],
        "env_to": dom["env_to"],
        "dom_start": ds,
        "dom_end": de,
        "dom_len": len(dom_seq),
        "prot": prot,
        "dom_seq": dom_seq,
        "cds": cds,
        "dom_cds": dom_cds,
        "cds_verified": cds is not None,
        "has_SXXK_in_dom": bool(re.search(r"S[A-Z]{2}K", dom_seq)),
        "domain_architecture": arch,
        **extra,
    }


def split_a2m_match_states(a2m_row):
    """Split one A2M row into match-state characters and the insertions that follow each state.

    Returns (states, inserts): states[k-1] is the residue (or '-') at HMM match state k;
    inserts[k-1] lists inserted residues between states k and k+1. Leading insertions
    (before state 1) are discarded.
    """
    states, inserts = [], []
    for ch in a2m_row:
        if ch == ".":
            continue
        if ch.isupper() or ch == "-":
            states.append(ch)
            inserts.append([])
        elif ch.islower() and inserts:
            inserts[-1].append(ch.upper())
    return states, inserts


def annotate_card_and_motifs(all_records):
    """Profile all curated PRE sequences against official CARD beta-lactamases & hmmalign catalytic motifs."""
    print("[CARD+Motif] Filtering PBP-type entries out of the CARD reference set...")
    card_kept, card_pbp = [], []
    for sid, hdr, seq in read_fasta(CARD_FAA):
        (card_pbp if CARD_PBP_RE.search(hdr) else card_kept).append((hdr, seq))
    card_clean = WORK / "card_true_sbl_refs.noPBP.faa"
    write_fasta(card_clean, card_kept)
    (WORK / "card_excluded_PBP_entries.txt").write_text("\n".join(h for h, _ in card_pbp) + "\n", encoding="utf-8")
    print(f"  CARD refs kept: {len(card_kept):,} | PBP-type removed: {len(card_pbp):,} (listed in card_excluded_PBP_entries.txt)")

    print("[CARD+Motif] Running DIAMOND blastp vs PBP-free CARD beta-lactamases...")
    dmnd_db = WORK / "card_true_sbl.dmnd"
    run_cmd(f"{DIAMOND_BIN} makedb --in {card_clean} -d {dmnd_db} >/dev/null 2>&1")

    q_map = {}
    q_fa = WORK / "all_for_card_dom.faa"
    items = []
    for idx, r in enumerate(all_records):
        sid = f"Q_{idx:06d}"
        q_map[sid] = r
        items.append((sid, r["dom_seq"]))
    write_fasta(q_fa, items)

    blast_out = WORK / "all_vs_card.tsv"
    run_cmd(
        f"{DIAMOND_BIN} blastp --query {q_fa} --db {dmnd_db} --very-sensitive -e 1e-5 -k 1 --threads 16 "
        f"--outfmt 6 qseqid stitle pident qcovhsp scovhsp bitscore -o {blast_out} >/dev/null 2>&1"
    )
    for r in all_records:
        r.update({
            "best_card_hit": "NONE",
            "card_gene_family": "NONE",
            "card_pident": 0.0,
            "card_cov": 0.0,
        })
    with open(blast_out, "r", encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            r = q_map[p[0]]
            stitle = p[1]
            pident = float(p[2])
            qcov = max(float(p[3]), float(p[4]))
            parts = stitle.split("|")
            gf = parts[1].replace("_", " ") if len(parts) >= 2 else stitle
            nm = parts[2].replace("_", " ") if len(parts) >= 3 else stitle
            r["best_card_hit"] = f"{parts[0]}|{nm}"
            r["card_gene_family"] = gf
            r["card_pident"] = round(pident, 2)
            r["card_cov"] = round(qcov, 2)

    print("[CARD+Motif] Running domain-anchored hmmalign catalytic triad + Omega-E166 extraction...")
    hmm_lens = load_hmm_lengths()
    by_fam = defaultdict(list)
    for r in all_records:
        by_fam[r["winning_acc"]].append(r)

    for fam, recs in by_fam.items():
        hmm_path = INPUTS / f"{fam}_official.hmm"
        win_cfg = MOTIF_WINDOWS.get(fam)
        if not hmm_path.exists() or win_cfg is None:
            for r in recs:
                m1 = bool(re.search(r"S[A-Z]{2}K", r["dom_seq"]))
                m2 = "SXN_PBP" if re.search(r"S[A-Z]N", r["dom_seq"]) else ("YAN" if re.search(r"Y[AS]N", r["dom_seq"]) else "NONE")
                m3 = bool(re.search(r"[KH][TS]G", r["dom_seq"]))
                r.update({"motif1_SXXK_hmm": m1, "motif2_class": m2, "motif3_KTG_hmm": m3, "has_omega_E166": False})
            continue

        fa_p = WORK / f"{fam}_for_motif_align.faa"
        sid_to_rec = {}
        fa_items = []
        for idx, r in enumerate(recs):
            sid = f"M_{idx:06d}"
            sid_to_rec[sid] = r
            fa_items.append((sid, r["dom_seq"]))
        write_fasta(fa_p, fa_items)
        # A2M marks match states unambiguously (uppercase residue or '-'), insertions as lowercase
        # and insert-column padding as '.', so match-state coordinates can be read per sequence.
        a2m_p = WORK / f"{fam}_for_motif_align.a2m"
        run_cmd(f"{HMMALIGN_BIN} --outformat A2M {hmm_path} {fa_p} > {a2m_p}")
        aln_seqs = {sid: seq for sid, _, seq in read_fasta(a2m_p)}
        hmm_len = hmm_lens[fam]

        for sid, row in aln_seqs.items():
            r = sid_to_rec[sid]
            states, inserts = split_a2m_match_states(row)
            if len(states) != hmm_len:
                raise RuntimeError(
                    f"{fam}: {sid} has {len(states)} match states but the HMM has LENG={hmm_len}; "
                    "motif windows would be misplaced (check the hmmalign A2M output)."
                )
            mseq = "".join(states)
            s1, e1 = win_cfg["m1"]
            seg1 = mseq[s1 - 1 : e1].replace("-", "")
            has_m1 = bool(re.search(r"S[A-Z]{2}K", seg1))

            s2, e2 = win_cfg["m2"]
            seg2 = mseq[s2 - 1 : e2].replace("-", "")
            if re.search(r"[SY]DN", seg2):
                m2_cls = "SDN_ClassA"
            elif re.search(r"Y[AS]N", seg2):
                m2_cls = "YAN_ClassC"
            elif re.search(r"YGN", seg2):
                m2_cls = "YGN_ClassD"
            elif re.search(r"[SY][A-Z]N", seg2):
                m2_cls = "SXN_PBP"
            else:
                m2_cls = "NONE"

            s3, e3 = win_cfg["m3"]
            seg3 = mseq[s3 - 1 : e3].replace("-", "")
            has_m3 = bool(re.search(r"[KH][TS]G", seg3))

            # Omega-loop: require the E166-x-x-L-N signature inside the calibrated window (+/- pad),
            # not merely any glutamate (which occurs by chance in about half of random windows).
            has_e166 = False
            if win_cfg["omega_e"] is not None:
                w_lo = max(1, win_cfg["omega_e"][0] - OMEGA_WINDOW_PAD)
                w_hi = min(hmm_len, win_cfg["omega_e"][1] + OMEGA_WINDOW_PAD)
                raw_om = "".join(
                    (states[k - 1] if states[k - 1] != "-" else "") + "".join(inserts[k - 1])
                    for k in range(w_lo, w_hi + 1)
                )
                has_e166 = bool(OMEGA_E166_RE.search(raw_om))

            r["motif1_SXXK_hmm"] = has_m1
            r["motif2_class"] = m2_cls
            r["motif3_KTG_hmm"] = has_m3
            r["has_omega_E166"] = has_e166

    # Assign final functional stratification label
    for r in all_records:
        card_ok = (r["card_pident"] >= 40.0) and (r["card_cov"] >= 70.0)
        fam = r["winning_acc"]
        is_cls_a = (fam in ("PF00144", "PF13354")) and r["motif1_SXXK_hmm"] and (r["motif2_class"] == "SDN_ClassA") and r["has_omega_E166"] and r["motif3_KTG_hmm"]
        is_cls_c = (fam in ("PF00144", "PF13354")) and r["motif1_SXXK_hmm"] and (r["motif2_class"] == "YAN_ClassC") and r["motif3_KTG_hmm"]
        is_cls_d = (fam == "PF00905") and card_ok and r["motif1_SXXK_hmm"] and (r["motif2_class"] == "YGN_ClassD")
        if card_ok or is_cls_a or is_cls_c or is_cls_d:
            r["functional_tier"] = "Putative_Hydrolytic_BetaLactamase"
            ev = []
            if card_ok:
                ev.append("CARD_Match")
            if is_cls_a:
                ev.append("ClassA_SXXK+SDN+E166+KTG")
            if is_cls_c:
                ev.append("ClassC_SXXK+YAN+KTG")
            if is_cls_d:
                ev.append("ClassD_OXA_YGN")
            r["functional_evidence"] = ";".join(ev)
        else:
            r["functional_tier"] = "PBP_or_PRE_Homolog"
            r["functional_evidence"] = f"PBP_PRE_{fam}"


def mmseqs_cluster_records(recs, tag, id_thresh=0.85):
    """Run MMseqs2 clustering on domain slices and return clusters dict {rep_uid: [rec, ...]}."""
    if not recs:
        return {}
    by_u = {r["uid"]: r for r in recs}
    fa_p = WORK / f"{tag}_in.faa"
    write_fasta(fa_p, [(r["uid"], r["dom_seq"]) for r in recs])
    cpref = WORK / f"mmseqs_{tag}"
    td = WORK / f"tmp_{tag}"
    td.mkdir(exist_ok=True)
    run_cmd(f"{MMSEQS_BIN} easy-cluster {fa_p} {cpref} {td} --min-seq-id {id_thresh} -c 0.80 --cov-mode 0 --threads 16 >/dev/null")
    clusters = defaultdict(list)
    with open(f"{cpref}_cluster.tsv", "r", encoding="utf-8") as f:
        for line in f:
            rep, mem = line.rstrip("\n").split("\t")
            if mem in by_u:
                clusters[rep].append(by_u[mem])
    return clusters


def pick_balanced_bacterial_sample(nr85_reps, cap=1500, seed=None, prefer_phylo=False):
    """Select <=cap taxonomy-balanced bacterial reps; seed=None gives deterministic primary (rep1)."""
    if len(nr85_reps) <= cap:
        return list(nr85_reps)

    rng = random.Random(seed) if seed is not None else None
    sel = []
    sel_uids = set()

    # Always preserve 100% of CARD/SARG clinical reps + 100% of bacterial prophage reps (needed for sensitivity)
    mandatory = [r for r in nr85_reps if r["source_db"] == "CARD_SARG_Clinical" or r["in_bacterial_prophage"]]
    mandatory.sort(key=lambda x: (-x["dom_score"], x["uid"]))
    for r in mandatory:
        sel.append(r)
        sel_uids.add(r["uid"])

    def sort_or_shuffle(lst):
        if rng is not None:
            out = list(lst)
            rng.shuffle(out)
            return out
        return sorted(lst, key=lambda x: (x["source_priority"], x["asm_rank"], -x["dom_score"], x["uid"]))

    # Stage 1: 1 GTDB rep per unique GTDB family, round-robin across bacterial phyla
    gtdb_by_fam = defaultdict(list)
    for r in nr85_reps:
        if r["source_db"] == "GTDB_RefSeq_3000" and r["uid"] not in sel_uids:
            gtdb_by_fam[r["family"]].append(r)
    fam_by_phy = defaultdict(list)
    for _, lst in gtdb_by_fam.items():
        ordered = sort_or_shuffle(lst)
        fam_by_phy[ordered[0]["phylum"]].append(ordered[0])

    # Stage 2: Marine MAG reps (up to 150) & Stage 3: remaining nr85_reps across phyla
    mar_by_phy = defaultdict(list)
    rem_by_phy = defaultdict(list)
    for r in nr85_reps:
        if r["uid"] not in sel_uids:
            if r["source_db"] == "Marine_MAGs_259":
                mar_by_phy[r["phylum"]].append(r)
            else:
                rem_by_phy[r["phylum"]].append(r)

    def round_robin_fill(bucket_map, max_add=None):
        added = 0
        phyla = sorted(bucket_map)
        if rng is not None:
            rng.shuffle(phyla)
        ordered_map = {p: sort_or_shuffle(bucket_map[p]) for p in phyla}
        idx_map = defaultdict(int)
        while len(sel) < cap and (max_add is None or added < max_add):
            prog = 0
            for p in phyla:
                if len(sel) >= cap or (max_add is not None and added >= max_add):
                    break
                lst = ordered_map[p]
                while idx_map[p] < len(lst) and lst[idx_map[p]]["uid"] in sel_uids:
                    idx_map[p] += 1
                if idx_map[p] < len(lst):
                    r = lst[idx_map[p]]
                    idx_map[p] += 1
                    sel.append(r)
                    sel_uids.add(r["uid"])
                    added += 1
                    prog += 1
            if prog == 0:
                break

    round_robin_fill(fam_by_phy)
    round_robin_fill(mar_by_phy, max_add=150)
    round_robin_fill(rem_by_phy)
    return sel


def build_viral_and_bacterial_pools():
    """Stage A: Curate IMG/VR + MOSAiC viruses and Triple-Source bacteria under strict CL0013 bitscore & own-seed FIMO."""
    acc_to_clan = load_pfam_clans()
    hmm_lens = load_hmm_lengths()

    # 1. VIRAL POOL
    print("[1/5] Loading IMG/VR Scheme B survivors & verifying 1:1 viral CDSs...")
    imgvr_keep_meta = {}
    with open(P1_DIR / "Phase1_Viral_Decontamination_Audit.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["schemeB_status"].startswith("KEEP"):
                imgvr_keep_meta[r["sbl_protein_id"]] = r

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
            cat, win_acc, dom, dual, p13_sc, p07_sc, hcov, arch = resolve_strict_cl0013(hit_d.get(uid, []), acc_to_clan, hmm_lens)
            v_cat_counts[(src, cat)] += 1
            if cat.startswith("drop_"):
                continue
            safe_uid = f"MOSAIC_{uid}" if src == "MOSAiC" else uid
            safe_uid = re.sub(r"[\s:(),;\[\]'\"]+", "_", safe_uid)
            rec = extract_domain_slice_and_build_record(
                safe_uid, uid, prot, dom, dual, p13_sc, p07_sc, hcov, arch,
                {
                    "origin": "Viral",
                    "source": src,
                    "pfam_category": cat,
                    "cds": viral_cds_map.get(uid) if src == "IMG_VR" else None,
                    "vir_sub_info": meta_d[uid].get("schemeB_status", "") if src == "IMG_VR" else meta_d[uid].get("vir_bin", ""),
                    "contig_context_available": src == "IMG_VR",
                },
            )
            viral_records.append(rec)

    print("Viral Strict CL0013 breakdown:", dict(v_cat_counts))

    # 2. TRIPLE-SOURCE BACTERIAL POOL (RETAINING & TAGGING PROPHAGES)
    print("[2/5] Loading Triple-Source Bacterial Pool (prophages retained with explicit tag)...")
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

    def classify_bac_locus(key, st, en):
        """viral_contig (whole contig called viral) | prophage (inside a provirus interval) | plasmid | chromosome."""
        if key in full_vir:
            return "viral_contig"
        if key in pv_ivs and st > 0 and any(interval_overlap(st, en, ps, pe) > 0 for ps, pe in pv_ivs[key]):
            return "prophage"
        return "plasmid" if key in plasmids else "chromosome"

    for uid, r in gtdb_hits_tsv.items():
        loci_strs = [x for x in r["all_verified_loci"].split(";") if ":" in x and x != "NA"]
        if not loci_strs:
            loc_classes = [classify_bac_locus(f"{r['acc']}|{r['contig_acc']}", int(r["cds_start"]), int(r["cds_end"]))]
        else:
            loc_classes = []
            for loc in loci_strs:
                contig_c, span = loc.split(":", 1)
                m = re.match(r"^(\d+)\.\.(\d+)", span)
                st, en = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
                loc_classes.append(classify_bac_locus(f"{r['acc']}|{contig_c}", st, en))
        # A protein counts as bacterial-genomic if any verified copy lies on a chromosome/plasmid.
        if "chromosome" in loc_classes:
            replicon = "chromosome"
        elif "plasmid" in loc_classes:
            replicon = "plasmid"
        elif "prophage" in loc_classes:
            replicon = "prophage"
        else:
            replicon = "viral_contig"
        in_pv_only = replicon == "prophage"

        prot = gtdb_prots[uid]
        cat, win_acc, dom, dual, p13_sc, p07_sc, hcov, arch = resolve_strict_cl0013(gtdb_domtbl.get(uid, []), acc_to_clan, hmm_lens)
        b_counts[("GTDB_3000", cat, replicon if replicon in ("prophage", "viral_contig") else "non_prophage")] += 1
        if cat.startswith("drop_"):
            continue
        safe_uid = re.sub(r"[\s:(),;\[\]'\"]+", "_", f"GTDB|{uid}")
        bac_records.append(extract_domain_slice_and_build_record(
            safe_uid, uid, prot, dom, dual, p13_sc, p07_sc, hcov, arch,
            {
                "origin": "Bacterial",
                "source_db": "GTDB_RefSeq_3000",
                "source_priority": 2,
                "asm_rank": ASM_RANK.get(r["assembly_level"], 4),
                "pfam_category": cat,
                "cds": gtdb_cds.get(uid),
                "replicon": replicon,
                "in_bacterial_prophage": in_pv_only,
                "is_viral_contig_in_genome": replicon == "viral_contig",
                "phylum": norm_tax(r["gtdb_phylum"], "p__"),
                "family": norm_tax(r["gtdb_family"], "f__"),
                "genus": norm_tax(r["gtdb_genus"], "g__"),
                "species": r["gtdb_species"],
                "is_clinical_mobile_bla": False,
            },
        ))

    # Source 2: Marine_MAGs_259 (keep KEEP_* AND non-partial DROP_genomad_full_viral_contig as tagged prophage)
    mar_rows = {}
    mar_pv_recover = {}
    with open(MAR_DIR / "Marine_MAGs_SBL_Audit.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["final_status"].startswith("KEEP_"):
                mar_rows[r["uid"]] = (r, "non_prophage")
            elif r["final_status"] in ("DROP_genomad_full_viral_contig", "DROP_genomad_prophage_interval") and r["is_partial_orf"] == "False" and r["cds_verified"] == "True":
                # A whole viral contig binned into a MAG is most likely a free phage genome, not an
                # integrated prophage: keep it tagged as viral_contig and out of the primary bacterial set.
                mar_rows[r["uid"]] = (r, "viral_contig" if r["final_status"] == "DROP_genomad_full_viral_contig" else "prophage")
                mar_pv_recover[r["uid"]] = r

    mar_prots = {uid: seq.rstrip("*") for uid, _, seq in read_fasta(MAR_DIR / "marine_mags_raw_candidates.faa") if uid in mar_rows}
    mar_cds = {uid: seq.upper() for uid, _, seq in read_fasta(MAR_DIR / "marine_mags_sbl_cds_verified.fna")}
    # Recover exact CDS from contig DNA for recovered complete prophage ORFs
    for uid, r in mar_pv_recover.items():
        fna_p = hpc.MARINE_MAGS_DIR / "sbl_fna" / "contigs_containing_pfam_from_data_descriptor_mags" / f"{r['mag_id']}.sbl_contigs.fna"
        if fna_p.exists():
            cdna_map = {cid: dna for cid, _, dna in read_fasta(fna_p)}
            cdna = cdna_map.get(r["contig_id"])
            if cdna:
                st, en = int(r["cds_start"]), int(r["cds_end"])
                seg = cdna[st - 1 : en]
                for dna_c in (seg, revcomp(seg)):
                    m_c = match_protein_to_dna(mar_prots[uid], dna_c)
                    if m_c:
                        mar_cds[uid] = m_c
                        break

    mar_domtbl = parse_hmmsearch_domtbl(MAR_DIR / "marine_mags_pfam_all.domtbl", filter_ids=set(mar_rows))
    for uid, (r, loc_class) in mar_rows.items():
        prot = mar_prots[uid]
        is_pv = loc_class == "prophage"
        cat, win_acc, dom, dual, p13_sc, p07_sc, hcov, arch = resolve_strict_cl0013(mar_domtbl.get(uid, []), acc_to_clan, hmm_lens)
        b_counts[("MarineMAG_259", cat, loc_class)] += 1
        if cat.startswith("drop_"):
            continue
        safe_uid = re.sub(r"[\s:(),;\[\]'\"]+", "_", f"MMAG|{uid}")
        bac_records.append(extract_domain_slice_and_build_record(
            safe_uid, uid, prot, dom, dual, p13_sc, p07_sc, hcov, arch,
            {
                "origin": "Bacterial",
                "source_db": "Marine_MAGs_259",
                "source_priority": 3,
                "asm_rank": 5,
                "pfam_category": cat,
                "cds": mar_cds.get(uid),
                "replicon": loc_class if loc_class != "non_prophage" else r["genomad_status"],
                "in_bacterial_prophage": is_pv,
                "is_viral_contig_in_genome": loc_class == "viral_contig",
                "phylum": norm_tax(r["gtdb_phylum"], "p__"),
                "family": norm_tax(r["gtdb_family"], "f__"),
                "genus": norm_tax(r["gtdb_genus"], "g__"),
                "species": r["gtdb_species"],
                "is_clinical_mobile_bla": False,
            },
        ))

    # Source 3: CARD_SARG_Clinical
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
        cat, win_acc, dom, dual, p13_sc, p07_sc, hcov, arch = resolve_strict_cl0013(step1_domtbl.get(safe_id, []), acc_to_clan, hmm_lens)
        b_counts[("CARD_SARG", cat, "non_prophage")] += 1
        if cat.startswith("drop_"):
            continue
        is_pl = bool(re.search(r"plasmid|transposon|integron|conjugat|mobile", hdr, re.I))
        is_clin = bool(re.search(r"\b(TEM|SHV|CTX-M|OXA|KPC|PER|VEB|GES|CMY|DHA|FOX|ACC|MIR|ACT|SME|IMI|NMC|BEL|TLA|CARB|PSE|ROB|BRO)[-_0-9]", hdr, re.I))
        raw_tag = f"SARG|{sid}" if sid not in seen_card else f"SARG|{ofam}|{sid}"
        seen_card.add(sid)
        safe_uid = re.sub(r"[\s:(),;\[\]'\"]+", "_", raw_tag)
        bac_records.append(extract_domain_slice_and_build_record(
            safe_uid, sid, prot, dom, dual, p13_sc, p07_sc, hcov, arch,
            {
                "origin": "Bacterial",
                "source_db": "CARD_SARG_Clinical",
                "source_priority": 1 if r["source_type"] == "CARD_SARG_curated" else 1.5,
                "asm_rank": 0,
                "pfam_category": cat,
                "cds": None,
                "replicon": "plasmid_or_MGE" if (is_pl or is_clin) else "chromosome_or_unspecified",
                "in_bacterial_prophage": False,
                "is_viral_contig_in_genome": False,
                "phylum": norm_tax(r["sarg_phylum"], "p__"),
                "family": "f__CARD_SARG_clinical",
                "genus": "g__CARD_SARG",
                "species": r["sarg_species"],
                "is_clinical_mobile_bla": is_clin or is_pl,
            },
        ))

    print("Bacterial Triple-Source Strict CL0013 breakdown:", dict(b_counts))

    # 3. OWN-FAMILY OFFICIAL-PFAM-SEED FIMO (p <= 1e-5) FOR ALL 5 PRE FAMILIES
    print("[3/5] Running symmetric own-family official-seed FIMO (p <= 1e-5) for all 5 PRE families...")
    fimo_summary = []
    fimo_pass_rate = {}
    for fam in PRE_FAMS:
        meme_txt = ensure_meme_model(fam)
        for kind, recs_all in [("viral", viral_records), ("bacterial", bac_records)]:
            cands = [r for r in recs_all if r["winning_acc"] == fam]
            if not cands:
                continue
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
            p5, p10, c65 = 0, 0, 0
            for r in cands:
                st = st_map.get(r["uid"], {"min_p": 1.0, "max_score": 0.0, "motifs": set()})
                r["fimo_r1_min_p"] = st["min_p"]
                r["fimo_r1_n_motifs"] = len(st["motifs"])
                r["pass_fimo_1e5"] = st["min_p"] <= 1e-5
                r["pass_fimo_1e10"] = st["min_p"] <= 1e-10
                if r["pass_fimo_1e5"]:
                    p5 += 1
                    if r["pass_cov65"]:
                        c65 += 1
                if r["pass_fimo_1e10"]:
                    p10 += 1
            fimo_pass_rate[(fam, kind)] = 100 * p5 / max(len(cands), 1)
            fimo_summary.append(
                f"  {fam} {kind:9s}: raw_GA={len(cands):,} -> FIMO<=1e-5={p5:,} ({100*p5/max(len(cands),1):.2f}%) -> &cov>=0.65={c65:,} (p<=1e-10: {p10:,})"
            )

    # Asymmetry check: an own-seed MEME model trained on few, closely related bacterial seed
    # sequences can reject divergent viral homologues preferentially.
    for fam in PRE_FAMS:
        v_rate, b_rate = fimo_pass_rate.get((fam, "viral")), fimo_pass_rate.get((fam, "bacterial"))
        if v_rate is None or b_rate is None:
            continue
        gap = b_rate - v_rate
        flag = "  <-- WARNING: viral pass rate much lower; consider a larger/more diverse MEME training set" if gap > FIMO_BIAS_WARN_PP else ""
        fimo_summary.append(f"  [bias check] {fam}: bacterial {b_rate:.2f}% - viral {v_rate:.2f}% = {gap:+.2f} pp{flag}")

    for line in fimo_summary:
        print(line)

    # Set default FIMO columns for non-PRE outgroups
    for r in viral_records + bac_records:
        if r["winning_acc"] not in PRE_FAMS:
            r["fimo_r1_min_p"] = 1.0
            r["fimo_r1_n_motifs"] = 0
            r["pass_fimo_1e5"] = "OUTGROUP_EXEMPT"
            r["pass_fimo_1e10"] = "OUTGROUP_EXEMPT"

    # 4. SBL vs PBP/PRE CARD + HMMALIGN CATALYTIC FUNCTIONAL STRATIFICATION
    print("[4/5] Running CARD homology + hmmalign catalytic stratification across PRE pools...")
    pre_all = [r for r in (viral_records + bac_records) if r["winning_acc"] in PRE_FAMS]
    annotate_card_and_motifs(pre_all)
    for r in viral_records + bac_records:
        if r["winning_acc"] not in PRE_FAMS:
            r.update({
                "best_card_hit": "OUTGROUP", "card_gene_family": "OUTGROUP",
                "card_pident": 0.0, "card_cov": 0.0,
                "motif1_SXXK_hmm": False, "motif2_class": "NONE",
                "motif3_KTG_hmm": False, "has_omega_E166": False,
                "functional_tier": "CL0013_NonPRE_Outgroup",
                "functional_evidence": f"Outgroup_{r['winning_acc']}",
            })

    return viral_records, bac_records, fimo_summary, v_cat_counts, b_counts


CKPT_PATH = WORK / "_curated_pools_ckpt.pkl"
P2_REQUIRED_COLUMN = "flanked_strict_structural_hallmark_any"   # only written by the revised Phase 2


def load_current_phase2_tiers(viral_records):
    """Load Phase 2 tiers, refusing stale output from the pre-revision Phase 2 run.

    Required order: --stage curate -> build_phase2_neighborhoods.py -> --stage finalize.
    """
    order = ("Run order: build_final_phase1_datasets.py --stage curate -> "
             "build_phase2_neighborhoods.py -> build_final_phase1_datasets.py --stage finalize "
             "(see run_phase1_phase2_rebuild.sh).")
    if not P2_TSV.exists():
        raise FileNotFoundError(f"Phase 2 table not found: {P2_TSV}\n{order}")
    with open(P2_TSV, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    cols = set(rows[0].keys()) if rows else set()
    if P2_REQUIRED_COLUMN not in cols:
        raise RuntimeError(
            f"{P2_TSV} was produced by the OLD Phase 2 (no '{P2_REQUIRED_COLUMN}' column; "
            f"tiers still count integrases as viral evidence).\n{order}"
        )
    if CKPT_PATH.exists() and P2_TSV.stat().st_mtime < CKPT_PATH.stat().st_mtime:
        raise RuntimeError(f"{P2_TSV} is older than the curation checkpoint; re-run Phase 2 first.\n{order}")
    tiers = {r["uid"]: r["viral_origin_tier"] for r in rows}
    imgvr = [r["uid"] for r in viral_records if r["source"] == "IMG_VR"]
    covered = sum(1 for u in imgvr if u in tiers)
    if imgvr and covered / len(imgvr) < 0.95:
        raise RuntimeError(
            f"Phase 2 covers only {covered:,}/{len(imgvr):,} IMG/VR loci of this curation (uid mismatch?).\n{order}"
        )
    print(f"  Phase 2 tiers loaded for {covered:,}/{len(imgvr):,} IMG/VR loci.")
    return tiers


def finalize_and_export(viral_records, bac_records, fimo_summary, v_cat_counts, b_counts):
    """Stage C: Cluster viral & bacterial pools independently (nr85), integrate Phase 2 Strict Tier 1, and export tree/HyPhy inputs."""
    p2_tier_map = load_current_phase2_tiers(viral_records)

    for r in viral_records:
        r["viral_origin_tier"] = p2_tier_map.get(r["uid"], "MOSAiC_Proteome_Verified_vMAG" if r["source"] == "MOSAiC" else "UNSCORED")

    summary = [
        "=== Revised Final Phase 1 + Phase 1B Curated Dataset Assembly (2026) ===",
        "Rules: Strict CL0013 Highest Bitscore | Own-Seed FIMO (p<=1e-5) | hmm_cov>=0.65 |",
        "       Independent Viral & Bacterial Domain-Slice nr85 | Tagged Bacterial Prophages |",
        "       Strict Structural Hallmark Phase 2 Tier 1 | SBL vs PBP/PRE Functional Stratification.",
        "\nViral Strict CL0013 Breakdown:",
    ] + [f"  {k}: {v}" for k, v in sorted(v_cat_counts.items())] + [
        "\nBacterial Triple-Source Strict CL0013 Breakdown:",
    ] + [f"  {k}: {v}" for k, v in sorted(b_counts.items())] + [
        "\nOwn-Family Official-Seed FIMO (p <= 1e-5) & Coverage >= 0.65 Filter:",
    ] + fimo_summary

    print("[5/5] Clustering independent viral & bacterial nr85 domain slices & building multi-replicate tree inputs...")
    summary.append("\nIndependent Domain-Slice nr85 Clustering & Functional Tier Breakdown (cov >= 0.65 & FIMO <= 1e-5):")

    viral_cluster_info = {}
    bac_cluster_info = {}
    bac_rep_sets = defaultdict(set)  # uid -> set of tags ('nr85', 'rep1'..'rep5', 'phylo')
    rep1_by_fam = {}
    viral_nr85_by_fam = {}

    family_groups = [(fam, (fam,)) for fam in EXPORT_FAMS] + [("PF13354_PF00768", ("PF13354", "PF00768"))]

    for tag, accs in family_groups:
        # 1. Primary curated pools (own-seed FIMO <= 1e-5 & hmm_cov >= 0.65)
        v_cur = [r for r in viral_records if r["winning_acc"] in accs and r["pass_fimo_1e5"] is True and r["pass_cov65"]]
        # Whole viral contigs found inside bacterial assemblies/MAGs are kept in the catalog but not used as bacteria.
        b_cur = [r for r in bac_records if r["winning_acc"] in accs and r["pass_fimo_1e5"] is True and r["pass_cov65"]
                 and not r["is_viral_contig_in_genome"]]

        # Export unclustered full primary viral sets (for ecology / abundance)
        write_fasta(OUT / f"{tag}_viral_clean.faa", [(r["uid"], r["prot"]) for r in v_cur])
        write_fasta(OUT / f"{tag}_viral_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_cur])
        write_fasta(OUT / f"{tag}_viral_cds.fna", [(r["uid"], r["cds"]) for r in v_cur if r["cds"]])
        write_fasta(OUT / f"{tag}_viral_dom_cds.fna", [(r["uid"], r["dom_cds"]) for r in v_cur if r["dom_cds"]])

        # 2. Independent VIRAL domain-slice MMseqs2 nr85 (NEVER mixed with bacteria)
        v_clus = mmseqs_cluster_records(v_cur, f"{tag}_viral_nr85", id_thresh=0.85)
        v_nr85_reps = []
        for idx_c, (_, mems) in enumerate(sorted(v_clus.items())):
            best_v = min(
                mems,
                key=lambda x: (
                    TIER_PRIORITY.get(x["viral_origin_tier"], 9),
                    0 if x["cds_verified"] else 1,
                    -x["dom_score"],
                    x["uid"],
                ),
            )
            v_nr85_reps.append((best_v, len(mems)))
            if tag != "PF13354_PF00768":
                cid_str = f"{tag}_v85_{idx_c:05d}"
                for m in mems:
                    viral_cluster_info[m["uid"]] = (cid_str, best_v["uid"], len(mems))

        v_nr85_only = [x[0] for x in v_nr85_reps]
        viral_nr85_by_fam[tag] = v_nr85_only
        write_fasta(OUT / f"{tag}_viral_nr85_clean.faa", [(r["uid"], r["prot"]) for r in v_nr85_only])
        write_fasta(OUT / f"{tag}_viral_nr85_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_nr85_only])
        write_fasta(OUT / f"{tag}_viral_nr85_cds.fna", [(r["uid"], r["cds"]) for r in v_nr85_only if r["cds"]])
        write_fasta(OUT / f"{tag}_viral_nr85_dom_cds.fna", [(r["uid"], r["dom_cds"]) for r in v_nr85_only if r["dom_cds"]])

        # 3. Independent BACTERIAL domain-slice MMseqs2 nr85
        b_clus = mmseqs_cluster_records(b_cur, f"{tag}_bac_nr85", id_thresh=0.85)
        b_nr85_reps = []
        hyphy_cds_reps = []
        for idx_c, (_, mems) in enumerate(sorted(b_clus.items())):
            best_b = min(
                mems,
                key=lambda x: (
                    x["source_priority"],
                    0 if x["cds_verified"] else 1,
                    x["asm_rank"],
                    -x["dom_score"],
                    x["uid"],
                ),
            )
            # Only copy CDS if identical sequence exists in cluster
            if not best_b["cds_verified"]:
                ex_p = [m for m in mems if m["cds_verified"] and m["prot"] == best_b["prot"] and m["dom_cds"]]
                if ex_p:
                    d0 = min(ex_p, key=lambda x: (x["source_priority"], x["asm_rank"], x["uid"]))
                    best_b["cds"] = d0["cds"]
                    best_b["dom_cds"] = d0["dom_cds"]
                    best_b["cds_verified"] = True
                else:
                    ex_d = [m for m in mems if m["cds_verified"] and m["dom_seq"] == best_b["dom_seq"] and m["dom_cds"]]
                    if ex_d:
                        d0 = min(ex_d, key=lambda x: (x["source_priority"], x["asm_rank"], x["uid"]))
                        best_b["dom_cds"] = d0["dom_cds"]
            # If any non-prophage chromosome/plasmid copy shares this exact nr85 cluster, cluster isn't exclusively prophage
            any_non_pv = any(not m["in_bacterial_prophage"] for m in mems)
            if not any_non_pv:
                best_b["in_bacterial_prophage"] = True

            cds_mems = [m for m in mems if m["cds_verified"] and m["dom_cds"]]
            if cds_mems:
                hyphy_cds_reps.append(min(cds_mems, key=lambda x: (x["asm_rank"], x["source_priority"], -x["dom_score"], x["uid"])))

            b_nr85_reps.append(best_b)
            if tag != "PF13354_PF00768":
                bac_rep_sets[best_b["uid"]].add("nr85")
                cid_str = f"{tag}_b85_{idx_c:05d}"
                for m in mems:
                    bac_cluster_info[m["uid"]] = (cid_str, len(mems))

        write_fasta(OUT / f"{tag}_bacterial_nr85_dom.faa", [(r["uid"], r["dom_seq"]) for r in b_nr85_reps])
        write_fasta(OUT / f"{tag}_bacterial_nr85_dom_cds.fna", [(r["uid"], r["dom_cds"]) for r in b_nr85_reps if r["dom_cds"]])
        write_fasta(OUT / f"{tag}_bacterial_nr85_hyphy_dom.faa", [(r["uid"], r["dom_seq"]) for r in hyphy_cds_reps])
        write_fasta(OUT / f"{tag}_bacterial_nr85_hyphy_dom_cds.fna", [(r["uid"], r["dom_cds"]) for r in hyphy_cds_reps])

        # 4. 6 Bacterial Downsampled Tree Replicates (rep1 Primary, rep2..rep5 Random Stratified, phylo_rep)
        rep1 = pick_balanced_bacterial_sample(b_nr85_reps, cap=1500, seed=None)
        rep1_by_fam[tag] = rep1
        for r in rep1:
            if tag != "PF13354_PF00768":
                bac_rep_sets[r["uid"]].add("rep1")

        write_fasta(OUT / f"{tag}_bacterial_tree_clean.faa", [(r["uid"], r["prot"]) for r in rep1])
        write_fasta(OUT / f"{tag}_bacterial_tree_dom.faa", [(r["uid"], r["dom_seq"]) for r in rep1])
        write_fasta(OUT / f"{tag}_bacterial_tree_cds.fna", [(r["uid"], r["cds"]) for r in rep1 if r["cds"]])
        write_fasta(OUT / f"{tag}_bacterial_tree_dom_cds.fna", [(r["uid"], r["dom_cds"]) for r in rep1 if r["dom_cds"]])
        write_fasta(OUT / f"{tag}_combined_tree_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_nr85_only] + [(r["uid"], r["dom_seq"]) for r in rep1])

        # Prophage-excluded sensitivity version of rep1
        rep1_nopro = [r for r in rep1 if not r["in_bacterial_prophage"]]
        write_fasta(OUT / f"{tag}_combined_tree_nopro_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_nr85_only] + [(r["uid"], r["dom_seq"]) for r in rep1_nopro])

        # rep2..rep5
        for rep_i, seed in enumerate((202602, 202603, 202604, 202605), start=2):
            rep_k = pick_balanced_bacterial_sample(b_nr85_reps, cap=1500, seed=seed)
            write_fasta(OUT / f"{tag}_combined_tree_rep{rep_i}_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_nr85_only] + [(r["uid"], r["dom_seq"]) for r in rep_k])

        # phylo_rep: deep-homology nr65 pre-clustering before 1,500 cap
        b_clus65 = mmseqs_cluster_records(b_nr85_reps, f"{tag}_bac_nr65", id_thresh=0.65)
        b_nr65_reps = [
            min(mems, key=lambda x: (x["source_priority"], 0 if x["cds_verified"] else 1, x["asm_rank"], -x["dom_score"], x["uid"]))
            for _, mems in sorted(b_clus65.items())
        ]
        phylo_rep = pick_balanced_bacterial_sample(b_nr65_reps, cap=1500, seed=None)
        write_fasta(OUT / f"{tag}_combined_tree_phylo_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_nr85_only] + [(r["uid"], r["dom_seq"]) for r in phylo_rep])

        # Coverage sensitivity (cov >= 0.70 and cov >= 0.60 on focal families)
        if tag in FOCAL:
            v_c70 = [r for r in v_nr85_only if r["pass_cov70"]]
            b_c70 = [r for r in rep1 if r["pass_cov70"]]
            write_fasta(OUT / f"{tag}_combined_tree_cov70_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_c70] + [(r["uid"], r["dom_seq"]) for r in b_c70])
            v_c60_all = [r for r in viral_records if r["winning_acc"] in accs and r["pass_fimo_1e5"] is True and r["pass_cov60"]]
            v_c60_clus = mmseqs_cluster_records(v_c60_all, f"{tag}_viral_cov60_nr85", id_thresh=0.85)
            v_c60_reps = [
                min(mems, key=lambda x: (TIER_PRIORITY.get(x["viral_origin_tier"], 9), 0 if x["cds_verified"] else 1, -x["dom_score"], x["uid"]))
                for _, mems in sorted(v_c60_clus.items())
            ]
            # Apply the relaxed 0.60 threshold symmetrically to bacteria as well
            b_c60_all = [r for r in bac_records if r["winning_acc"] in accs and r["pass_fimo_1e5"] is True and r["pass_cov60"]
                         and not r["is_viral_contig_in_genome"]]
            b_c60_clus = mmseqs_cluster_records(b_c60_all, f"{tag}_bac_cov60_nr85", id_thresh=0.85)
            b_c60_reps = [
                min(mems, key=lambda x: (x["source_priority"], 0 if x["cds_verified"] else 1, x["asm_rank"], -x["dom_score"], x["uid"]))
                for _, mems in sorted(b_c60_clus.items())
            ]
            rep1_c60 = pick_balanced_bacterial_sample(b_c60_reps, cap=1500, seed=None)
            write_fasta(OUT / f"{tag}_combined_tree_cov60_dom.faa", [(r["uid"], r["dom_seq"]) for r in v_c60_reps] + [(r["uid"], r["dom_seq"]) for r in rep1_c60])

        v_sbl_cnt = sum(1 for r in v_cur if r["functional_tier"] == "Putative_Hydrolytic_BetaLactamase")
        b_sbl_cnt = sum(1 for r in b_nr85_reps if r["functional_tier"] == "Putative_Hydrolytic_BetaLactamase")
        v85_sbl_cnt = sum(1 for r in v_nr85_only if r["functional_tier"] == "Putative_Hydrolytic_BetaLactamase")
        n_pv_rep1 = sum(1 for r in rep1 if r["in_bacterial_prophage"])
        sline = (
            f"  {tag:14s}: Viral cov65+FIMO={len(v_cur):,} (SBL={v_sbl_cnt:,}) -> Viral nr85={len(v_nr85_only):,} (SBL={v85_sbl_cnt:,}, MOSAiC={sum(1 for r in v_nr85_only if r['source']=='MOSAiC')}) | "
            f"Bac cov65+FIMO={len(b_cur):,} -> Bac nr85={len(b_nr85_reps):,} (SBL={b_sbl_cnt:,}) -> Tree rep1={len(rep1):,} (prophage={n_pv_rep1}, phylo_rep={len(phylo_rep):,})"
        )
        print(sline)
        summary.append(sline)

    # 5. EXPORT SISTER-PRE OUTGROUPS & ALL 3 NON-PRE CL0013 OUTGROUPS (PF04960, PF06089, PF27310)
    print("Exporting Sister-PRE outgroups & all non-PRE CL0013 biological outgroups...")
    # Collect diverse bacterial 60%-identity reps for Sister-PRE outgroups
    sister_pool = {}
    for fam in EXPORT_FAMS:
        fam_clus60 = mmseqs_cluster_records(rep1_by_fam[fam], f"{fam}_sister60", id_thresh=0.60)
        reps60 = [min(mems, key=lambda x: (x["asm_rank"], -x["dom_score"], x["uid"])) for _, mems in sorted(fam_clus60.items())]
        sister_pool[fam] = pick_balanced_bacterial_sample(reps60, cap=20, seed=None)

    for tag, accs in family_groups:
        # Trees touching the PF13354/PF00768 continuum are rooted only with PF00144 + PF00905.
        excluded = set(accs) | (SISTER_CONTINUUM if SISTER_CONTINUUM & set(accs) else set())
        sister_fams = [sfam for sfam in EXPORT_FAMS if sfam not in excluded]
        take_n = 24 // max(len(sister_fams), 1)
        sisters = []
        for sfam in sister_fams:
            sisters.extend(sister_pool[sfam][:take_n])
        write_fasta(OUT / f"{tag}_sister_pre_outgroup_dom.faa", [(r["uid"], r["dom_seq"]) for r in sisters])
        summary.append(f"  Sister-PRE outgroup for {tag}: {', '.join(sister_fams)} ({len(sisters)} seqs)")

    # Auxiliary PF02113 & ALL 3 non-PRE CL0013 outgroups (PF04960, PF06089, PF27310)
    summary.append("\nCL0013 Auxiliary (PF02113) & Non-PRE Biological Outgroups (PF04960, PF06089, PF27310):")
    aux_2113_v = [r for r in viral_records if r["winning_acc"] == "PF02113" and r["pass_fimo_1e5"] is True and r["pass_cov65"]]
    aux_2113_b = [r for r in bac_records if r["winning_acc"] == "PF02113" and r["pass_fimo_1e5"] is True and r["pass_cov65"]
                  and not r["is_viral_contig_in_genome"]]
    clus_2113 = mmseqs_cluster_records(aux_2113_b, "bac_PF02113_60", id_thresh=0.60)
    reps_2113 = pick_balanced_bacterial_sample(
        [min(mems, key=lambda x: (x["asm_rank"], -x["dom_score"], x["uid"])) for _, mems in sorted(clus_2113.items())],
        cap=100, seed=None,
    )
    write_fasta(OUT / "CL0013_aux_PF02113_viral_dom.faa", [(r["uid"], r["dom_seq"]) for r in aux_2113_v])
    write_fasta(OUT / "CL0013_aux_PF02113_bacterial_rep_dom.faa", [(r["uid"], r["dom_seq"]) for r in reps_2113])
    summary.append(f"  PF02113 (Peptidase_S13): Viral={len(aux_2113_v)}, Bacterial FIMO+cov65={len(aux_2113_b):,} -> 60%-rep={len(reps_2113)}")

    for out_fam in OUTGROUP:
        out_b = [r for r in bac_records if r["winning_acc"] == out_fam and r["pass_cov65"] and not r["is_viral_contig_in_genome"]]
        clus_o = mmseqs_cluster_records(out_b, f"out_{out_fam}_50", id_thresh=0.50)
        reps_o = pick_balanced_bacterial_sample(
            [min(mems, key=lambda x: (x["asm_rank"], -x["dom_score"], x["uid"])) for _, mems in sorted(clus_o.items())],
            cap=24, seed=None,
        )
        write_fasta(OUT / f"CL0013_outgroup_{out_fam}_dom.faa", [(r["uid"], r["dom_seq"]) for r in reps_o])
        write_fasta(OUT / f"CL0013_outgroup_{out_fam}_clean.faa", [(r["uid"], r["prot"]) for r in reps_o])
        summary.append(f"  {out_fam}: raw cov>=0.65={len(out_b):,} -> balanced 50%-rep={len(reps_o)}")

    # 6. WRITE UPDATED CATALOGS & SUMMARY
    v_cols = [
        "uid", "orig_uid", "source", "pfam_category", "winning_acc",
        "dual_PF13354_PF00768_GA", "pf13354_dom_score", "pf00768_dom_score",
        "dom_score", "hmm_cov", "pass_cov60", "pass_cov65", "pass_cov70",
        "seq_len", "env_from", "env_to", "dom_start", "dom_end", "dom_len",
        "fimo_r1_min_p", "fimo_r1_n_motifs", "pass_fimo_1e5", "pass_fimo_1e10",
        "viral_origin_tier", "mmseqs85_cluster_id", "viral_nr85_rep_uid", "viral_cluster_size", "in_viral_nr85_rep",
        "cds_verified", "has_SXXK_in_dom", "motif1_SXXK_hmm", "motif2_class", "motif3_KTG_hmm", "has_omega_E166",
        "best_card_hit", "card_gene_family", "card_pident", "card_cov",
        "functional_tier", "functional_evidence", "vir_sub_info", "domain_architecture",
    ]
    with open(OUT / "Final_Viral_SBL_Catalog.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=v_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in viral_records:
            row = dict(r)
            cid, rep_u, csz = viral_cluster_info.get(r["uid"], ("NA", "NA", 0))
            row["mmseqs85_cluster_id"] = cid
            row["viral_nr85_rep_uid"] = rep_u
            row["viral_cluster_size"] = csz
            row["in_viral_nr85_rep"] = (r["uid"] == rep_u)
            w.writerow(row)

    b_cols = [
        "uid", "orig_uid", "source_db", "pfam_category", "winning_acc",
        "dual_PF13354_PF00768_GA", "pf13354_dom_score", "pf00768_dom_score",
        "dom_score", "hmm_cov", "pass_cov60", "pass_cov65", "pass_cov70",
        "seq_len", "env_from", "env_to", "dom_start", "dom_end", "dom_len",
        "fimo_r1_min_p", "fimo_r1_n_motifs", "pass_fimo_1e5", "pass_fimo_1e10",
        "mmseqs85_cluster_id", "bac_cluster_size", "in_nr85_rep", "in_tree_balanced",
        "cds_verified", "replicon", "in_bacterial_prophage", "is_viral_contig_in_genome", "is_clinical_mobile_bla",
        "phylum", "family", "genus", "species",
        "has_SXXK_in_dom", "motif1_SXXK_hmm", "motif2_class", "motif3_KTG_hmm", "has_omega_E166",
        "best_card_hit", "card_gene_family", "card_pident", "card_cov",
        "functional_tier", "functional_evidence", "domain_architecture",
    ]
    with open(OUT / "Final_Bacterial_TripleSource_Catalog.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=b_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in bac_records:
            row = dict(r)
            cid, csz = bac_cluster_info.get(r["uid"], ("NA", 0))
            row["mmseqs85_cluster_id"] = cid
            row["bac_cluster_size"] = csz
            rtags = bac_rep_sets.get(r["uid"], set())
            row["in_nr85_rep"] = "nr85" in rtags
            row["in_tree_balanced"] = "rep1" in rtags
            w.writerow(row)

    rep_txt = "\n".join(summary) + "\n"
    (OUT / "Final_Phase1_Dataset_Summary.txt").write_text(rep_txt, encoding="utf-8")
    print("[Done] Revised Final Phase 1 + Phase 1B assembly complete!")


def write_prelim_viral_catalog(viral_records):
    v_cols = [
        "uid", "orig_uid", "source", "pfam_category", "winning_acc",
        "dual_PF13354_PF00768_GA", "pf13354_dom_score", "pf00768_dom_score",
        "dom_score", "hmm_cov", "pass_cov60", "pass_cov65", "pass_cov70",
        "seq_len", "env_from", "env_to", "dom_start", "dom_end", "dom_len",
        "fimo_r1_min_p", "fimo_r1_n_motifs", "pass_fimo_1e5", "pass_fimo_1e10",
        "cds_verified", "has_SXXK_in_dom", "motif1_SXXK_hmm", "motif2_class", "motif3_KTG_hmm", "has_omega_E166",
        "best_card_hit", "card_gene_family", "card_pident", "card_cov",
        "functional_tier", "functional_evidence", "vir_sub_info", "domain_architecture",
    ]
    with open(OUT / "Final_Viral_SBL_Catalog.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=v_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in viral_records:
            w.writerow(r)


def main():
    parser = argparse.ArgumentParser(
        description="Revised Phase 1 + Phase 1B Dataset Assembly. Run --stage curate, then "
                    "build_phase2_neighborhoods.py, then --stage finalize (run_phase1_phase2_rebuild.sh)."
    )
    parser.add_argument("--stage", choices=["curate", "finalize"], required=True)
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    import pickle
    if args.stage == "curate":
        v_rec, b_rec, f_sum, vc, bc = build_viral_and_bacterial_pools()
        write_prelim_viral_catalog(v_rec)
        with open(CKPT_PATH, "wb") as f:
            pickle.dump((v_rec, b_rec, f_sum, vc, bc), f)
        print(f"Saved curation checkpoint ({len(v_rec)} viral, {len(b_rec)} bacterial) & preliminary Final_Viral_SBL_Catalog.tsv.")
        print("Next: build_phase2_neighborhoods.py, then --stage finalize.")
    else:
        with open(CKPT_PATH, "rb") as f:
            v_rec, b_rec, f_sum, vc, bc = pickle.load(f)
        finalize_and_export(v_rec, b_rec, f_sum, vc, bc)


if __name__ == "__main__":
    main()
