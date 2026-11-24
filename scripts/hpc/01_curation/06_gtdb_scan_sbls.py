#!/usr/bin/env python3
"""Scan downloaded GTDB/RefSeq representative bacterial proteomes with official CL0013 & Pfam-A (--cut_ga),
apply overlapping-envelope (ov >= 30 aa) competition, verify 1:1 native bacterial CDS translation,
and generate genomic .fna download list for geNomad prophage/plasmid excision.
"""

import csv
import gzip
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

BDIR = hpc.BASE / "01b_bacterial_ref_db"
META_TSV = BDIR / "metadata" / "selected_3000_genomes.tsv"
FAA_DIR = BDIR / "downloads" / "faa"
CDS_DIR = BDIR / "downloads" / "cds"
FNA_DIR = BDIR / "downloads" / "fna"
OUT_DIR = BDIR / "sbl_contigs"

CL0013_HMM = hpc.BASE / "inputs" / "CL0013_official.hmm"
PFAM_ALL_HMM = hpc.PFAM_ALL_HMM
PFAM_CLANS_TSV = hpc.PFAM_CLANS_TSV
HMMSEARCH_BIN = hpc.HMMSEARCH_BIN

FOCAL_FAMILIES = ("PF00144", "PF00905", "PF13354")
AUX_PRE_FAMILIES = ("PF00768", "PF02113")
OUTGROUP_CL0013 = ("PF04960", "PF06089", "PF27310")
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


def verify_cds_translation(prot_seq, dna_seq):
    """Return exact coding DNA slice (len == 3*len(prot)) if translation matches Table 11/4."""
    target = prot_seq.rstrip("*")
    dna = dna_seq.upper()
    cands = [dna]
    if len(dna) >= 3 and len(dna) % 3 == 0 and dna[-3:] in ("TAA", "TAG", "TGA"):
        cands.append(dna[:-3])
    for cand in cands:
        if len(cand) < 3 * len(target):
            continue
        for tga_w in (False, True):
            for alt_st in (True, False):
                aa = translate_dna(cand, alt_start=alt_st, tga_w=tga_w).rstrip("*")
                if aa == target and len(cand) >= 3 * len(target):
                    return cand[: 3 * len(target)]
                # Allow selenocysteine U -> X mismatch
                if len(aa) == len(target) and "U" in target:
                    aa2 = "".join(t if t == "U" else a for a, t in zip(aa, target))
                    if aa2 == target:
                        return cand[: 3 * len(target)]
    return None


def open_gz_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore")


def iter_fasta(fobj):
    hdr, buf = None, []
    for line in fobj:
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            if hdr is not None:
                yield hdr.split()[0], hdr, "".join(buf)
            hdr = s[1:].strip()
            buf = []
        else:
            buf.append(s)
    if hdr is not None:
        yield hdr.split()[0], hdr, "".join(buf)


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


def classify_cl0013(winners):
    if not winners:
        return "drop_no_pfam_ga", "NONE", None
    cl13 = [w for w in winners if w["clan"] == "CL0013"]
    if not cl13:
        top_w = max(winners, key=lambda w: w["dom_score"])
        return "drop_non_CL0013_only", top_w["pfam_acc"], top_w
    best = max(cl13, key=lambda w: w["dom_score"])
    acc = best["pfam_acc"]
    if acc in FOCAL_FAMILIES:
        return f"focal_{acc}", acc, best
    if acc in AUX_PRE_FAMILIES:
        return f"aux_{acc}", acc, best
    if acc in OUTGROUP_CL0013:
        return f"outgroup_{acc}", acc, best
    return f"other_CL0013_{acc}", acc, best


def run_cmd(cmd):
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd}\nSTDERR:\n{p.stderr}")


def parse_cds_genome(cds_gz_path, target_prot_ids):
    """Map protein_id -> list of (contig_acc, start, end, strand, dna_seq)."""
    out = defaultdict(list)
    if not cds_gz_path.exists():
        return out
    try:
        with open_gz_text(cds_gz_path) as f:
            for tok, hdr, dna in iter_fasta(f):
                mp = re.search(r"\[protein_id=([^\]]+)\]", hdr)
                if not mp:
                    continue
                pid = mp.group(1).strip()
                if pid not in target_prot_ids:
                    continue
                # tok looks like lcl|NZ_CP012345.1_cds_WP_012345678.1_123
                mc = re.match(r"^lcl\|(.+?)_cds_", tok)
                contig_acc = mc.group(1) if mc else tok
                ml = re.search(r"\[location=([^\]]+)\]", hdr)
                if ml:
                    loc = ml.group(1)
                    nums = [int(x) for x in re.findall(r"\d+", loc)]
                    st, en = (min(nums), max(nums)) if nums else (0, 0)
                    strand = -1 if "complement" in loc else 1
                else:
                    st, en, strand = 0, 0, 1
                out[pid].append((contig_acc, st, en, strand, dna.upper()))
    except Exception as e:
        print(f"WARN parsing {cds_gz_path.name}: {e}")
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(META_TSV, "r", encoding="utf-8") as f:
        meta = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}
    print(f"[1/6] Loaded metadata for {len(meta)} GTDB/RefSeq genomes.")

    # Concatenate all proteomes with accession|prot_id headers
    all_faa = BDIR / "downloads" / "all_3000_proteins.faa"
    valid_genomes = 0
    bad_genomes = []
    n_prots = 0
    with open(all_faa, "w", encoding="utf-8") as out_f:
        for acc in sorted(meta):
            p = FAA_DIR / f"{acc}_protein.faa.gz"
            if not p.exists() or p.stat().st_size < 100:
                bad_genomes.append((acc, "missing_or_empty"))
                continue
            try:
                cnt_g = 0
                with open_gz_text(p) as gf:
                    for pid, hdr, seq in iter_fasta(gf):
                        desc = hdr.split(None, 1)[1] if " " in hdr else ""
                        out_f.write(f">{acc}|{pid} {desc}\n{seq.rstrip('*')}\n")
                        cnt_g += 1
                        n_prots += 1
                if cnt_g > 0:
                    valid_genomes += 1
                else:
                    bad_genomes.append((acc, "zero_proteins"))
            except Exception as e:
                bad_genomes.append((acc, str(e)))
    print(f"[2/6] Combined {n_prots:,} proteins across {valid_genomes} valid genomes (bad/missing={len(bad_genomes)}).")

    # Run CL0013 official HMMER search
    cl13_dom = OUT_DIR / "bac_gtdb_cl0013.domtbl"
    print("[3/6] Running hmmsearch CL0013_official.hmm (--cut_ga) on bacterial proteome...")
    run_cmd(f"{HMMSEARCH_BIN} --cpu 32 --cut_ga --domtblout {cl13_dom} {CL0013_HMM} {all_faa} >/dev/null")

    cand_ids = set()
    with open(cl13_dom, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                cand_ids.add(line.split()[0])
    print(f"  CL0013 GA candidate proteins: {len(cand_ids)}")

    # Extract candidate FASTA & store sequences
    cand_seqs = {}
    cand_faa = OUT_DIR / "bac_cl0013_candidates.faa"
    with open(all_faa, "r", encoding="utf-8") as inf, open(cand_faa, "w", encoding="utf-8") as outf:
        for uid, hdr, seq in iter_fasta(inf):
            if uid in cand_ids:
                cand_seqs[uid] = (hdr, seq)
                outf.write(f">{hdr}\n{seq}\n")

    # Run full Pfam-A (--cut_ga) on candidates
    pfam_all_dom = OUT_DIR / "bac_gtdb_pfam_all.domtbl"
    print(f"[4/6] Running hmmsearch full Pfam-A.hmm (--cut_ga) on {len(cand_seqs)} candidates...")
    run_cmd(f"{HMMSEARCH_BIN} --cpu 32 --cut_ga --domtblout {pfam_all_dom} {PFAM_ALL_HMM} {cand_faa} >/dev/null")

    # Overlapping-envelope competition
    acc_to_clan = load_pfam_clans()
    pfam_hits = parse_hmmsearch_domtbl(pfam_all_dom)
    cat_counts = Counter()
    kept_by_genome = defaultdict(dict)  # acc -> prot_id -> rec
    all_records = []

    for uid in sorted(cand_seqs):
        hdr, prot = cand_seqs[uid]
        acc, prot_id = uid.split("|", 1)
        winners = resolve_overlapping_domains(pfam_hits.get(uid, []), acc_to_clan, OV_THRESHOLD_AA)
        cat, win_acc, best_dom = classify_cl0013(winners)
        cat_counts[cat] += 1
        if cat.startswith("drop_"):
            continue
        dom_start = max(1, best_dom["env_from"] - DOMAIN_PAD_AA)
        dom_end = min(len(prot), best_dom["env_to"] + DOMAIN_PAD_AA)
        dom_seq = prot[dom_start - 1 : dom_end]
        arch_str = "+".join(f"{w['pfam_acc']}:{w['pfam_name']}({w['env_from']}-{w['env_to']})" for w in winners)
        rec = {
            "uid": uid,
            "acc": acc,
            "prot_id": prot_id,
            "category": cat,
            "winning_acc": win_acc,
            "dom_score": best_dom["dom_score"],
            "env_from": best_dom["env_from"],
            "env_to": best_dom["env_to"],
            "dom_start": dom_start,
            "dom_end": dom_end,
            "prot": prot,
            "dom_seq": dom_seq,
            "has_sxxk": bool(re.search(r"S[A-Z]{2}K", dom_seq)),
            "arch": arch_str,
        }
        kept_by_genome[acc][prot_id] = rec
        all_records.append(rec)

    print("  Overlapping-envelope Pfam competition:", dict(cat_counts))

    # Extract 1:1 verified CDS + coordinates from _cds_from_genomic.fna.gz
    print(f"[5/6] Verifying 1:1 CDS translation & genomic coordinates across {len(kept_by_genome)} hit genomes...")
    cds_ok = 0
    cds_fail = 0
    tsv_rows = []
    prot_fasta = []
    dom_fasta = []
    cds_fasta = []
    dom_cds_fasta = []
    fna_url_lines = []

    for acc in sorted(kept_by_genome):
        mrow = meta[acc]
        fna_url = mrow["fna_url"]
        fna_url_lines.append(f"{fna_url}\n  dir={FNA_DIR}\n  out={acc}_genomic.fna.gz\n")
        target_pids = set(kept_by_genome[acc].keys())
        cds_map = parse_cds_genome(CDS_DIR / f"{acc}_cds_from_genomic.fna.gz", target_pids)

        for prot_id, rec in kept_by_genome[acc].items():
            prot = rec["prot"]
            loci = cds_map.get(prot_id, [])
            verified_loci = []
            for contig_acc, st, en, strand, dna in loci:
                cslice = verify_cds_translation(prot, dna)
                if cslice is not None:
                    verified_loci.append((contig_acc, st, en, strand, cslice))
            if verified_loci:
                cds_ok += 1
                contig_acc, st, en, strand, cslice = verified_loci[0]
                all_locs_str = ";".join(f"{c}:{a}..{b}({'+' if s > 0 else '-'})" for c, a, b, s, _ in verified_loci)
                cds_verified = True
            else:
                cds_fail += 1
                contig_acc, st, en, strand, cslice = ("NA", 0, 0, 0, None)
                all_locs_str = "NA"
                cds_verified = False

            hdr_tag = (
                f"{rec['uid']} fam={rec['winning_acc']} cat={rec['category']} origin=bacterial_GTDB "
                f"phylum={mrow['phylum']} contig={contig_acc} loc={st}..{en} dom={rec['dom_start']}..{rec['dom_end']} score={rec['dom_score']:.1f}"
            )
            prot_fasta.append((hdr_tag, prot))
            dom_fasta.append((hdr_tag, rec["dom_seq"]))
            if cslice is not None:
                cds_fasta.append((hdr_tag, cslice))
                if len(cslice) >= 3 * rec["dom_end"]:
                    dom_cds_fasta.append((hdr_tag, cslice[(rec["dom_start"] - 1) * 3 : rec["dom_end"] * 3]))

            tsv_rows.append({
                "uid": rec["uid"],
                "acc": acc,
                "prot_id": prot_id,
                "pfam_category": rec["category"],
                "winning_acc": rec["winning_acc"],
                "dom_score": rec["dom_score"],
                "seq_len": len(prot),
                "env_from": rec["env_from"],
                "env_to": rec["env_to"],
                "dom_slice_start": rec["dom_start"],
                "dom_slice_end": rec["dom_end"],
                "has_SXXK_in_dom": rec["has_sxxk"],
                "domain_architecture": rec["arch"],
                "contig_acc": contig_acc,
                "cds_start": st,
                "cds_end": en,
                "cds_strand": strand,
                "n_genome_loci": len(verified_loci),
                "all_verified_loci": all_locs_str,
                "cds_verified": cds_verified,
                "gtdb_phylum": mrow["phylum"],
                "gtdb_class": mrow["class"],
                "gtdb_order": mrow["order"],
                "gtdb_family": mrow["family"],
                "gtdb_genus": mrow["genus"],
                "gtdb_species": mrow["species"],
                "assembly_level": mrow["assembly_level"],
            })

    with open(OUT_DIR / "bac_gtdb_sbl_hits.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(tsv_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(tsv_rows)
    write_fasta(OUT_DIR / "bac_gtdb_sbl_proteins.faa", prot_fasta)
    write_fasta(OUT_DIR / "bac_gtdb_sbl_dom.faa", dom_fasta)
    write_fasta(OUT_DIR / "bac_gtdb_sbl_cds_verified.fna", cds_fasta)
    write_fasta(OUT_DIR / "bac_gtdb_sbl_dom_cds_verified.fna", dom_cds_fasta)
    (BDIR / "metadata" / "urls_sbl_genomes_fna.txt").write_text("".join(fna_url_lines), encoding="utf-8")

    summary = [
        "=== GTDB/RefSeq 3,000 Bacterial Representative Scan Summary ===",
        f"Valid proteomes scanned: {valid_genomes} / {len(meta)} (bad/missing: {len(bad_genomes)})",
        f"Total bacterial proteins scanned: {n_prots:,}",
        f"CL0013 --cut_ga hits: {len(cand_seqs)}",
        f"Genomes carrying >=1 winning CL0013 hit: {len(kept_by_genome)}",
        f"1:1 Translation-verified bacterial CDSs: {cds_ok} / {cds_ok + cds_fail} ({100*cds_ok/max(cds_ok+cds_fail,1):.2f}%)",
        "\nCategory breakdown:",
    ] + [f"  {k}: {v}" for k, v in sorted(cat_counts.items())]
    rep = "\n".join(summary) + "\n"
    (OUT_DIR / "bac_gtdb_scan_summary.txt").write_text(rep, encoding="utf-8")
    print("[6/6] " + rep)


if __name__ == "__main__":
    main()
