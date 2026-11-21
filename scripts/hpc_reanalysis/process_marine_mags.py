#!/usr/bin/env python3
"""Audit, prophage-excise (geNomad coordinates), plasmid-tag, and 1:1 CDS-verify
the 259 dereplicated near-complete marine bacterial MAGs in $SBL_HPC_MARINE_MAGS_DIR.
"""

import csv
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import hpc_config as hpc

SFY = hpc.MARINE_MAGS_DIR
BEST_CSV = SFY / "best_genomes_quals.csv"
TAX_TSV = SFY / "data_descriptor_mags_cc.tsv"
FAA_DIR = SFY / "sbl_faa" / "pfam_hits_from_data_descriptor_mags"
FNA_DIR = SFY / "sbl_fna" / "contigs_containing_pfam_from_data_descriptor_mags"
PV_TSV = SFY / "plasmid_search_genomad" / "sbl_contigs_concat_find_proviruses" / "sbl_contigs_concat_provirus.tsv"
VS_TSV = SFY / "plasmid_search_genomad" / "sbl_contigs_concat_summary" / "sbl_contigs_concat_virus_summary.tsv"
PL_TSV = SFY / "plasmid_search_genomad" / "sbl_contigs_concat_summary" / "sbl_contigs_concat_plasmid_summary.tsv"

OUT_DIR = hpc.BASE / "01b_bacterial_ref_db" / "marine_mags"
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
START_CODONS = {"ATG", "GTG", "TTG", "ATT", "ATC", "CTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    recs = []
    hdr, buf = None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
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


def translate_dna(dna, alt_start=True):
    dna = dna.upper()
    n = len(dna) - (len(dna) % 3)
    aa = []
    for i in range(0, n, 3):
        cod = dna[i : i + 3]
        if i == 0 and alt_start and cod in START_CODONS:
            aa.append("M")
        else:
            aa.append(CODON_STD.get(cod, "X"))
    return "".join(aa)


def match_contig_slice_to_prot(contig_dna, st, en, prot_seq):
    """Extract [st-1:en] on + or - strand and verify 100% translation match."""
    target = prot_seq.rstrip("*")
    seg = contig_dna[st - 1 : en].upper()
    contig_len = len(contig_dna)
    is_edge = (st <= 3) or (en >= contig_len - 2)

    for strand, dna in ((1, seg), (-1, revcomp(seg))):
        if len(dna) < 3 * len(target):
            continue
        for cand in (dna, dna[:-3] if len(dna) >= 3 and dna[-3:] in STOP_CODONS else dna):
            if len(cand) < 3 * len(target):
                continue
            for alt_st in (True, False):
                aa = translate_dna(cand, alt_start=alt_st).rstrip("*")
                if aa == target:
                    has_start = dna[:3] in START_CODONS
                    has_stop = dna[-3:] in STOP_CODONS
                    return cand[: 3 * len(target)], strand, has_start, has_stop, is_edge
    return None, 0, False, False, is_edge


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


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 2a. Load 259 MAGs & taxonomy
    mags = {}
    with open(BEST_CSV, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            mid = r["genome"].replace(".fna", "").strip()
            mags[mid] = {"completeness": r["completeness"], "contamination": r["contamination"], "tax": []}

    with open(TAX_TSV, "r", encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] in mags:
                tax_parts = p[1].split(";")
                while len(tax_parts) < 7:
                    tax_parts.append("")
                mags[p[0]]["tax"] = tax_parts[:7]
    print(f"[1/5] Loaded {len(mags)} dereplicated marine bacterial MAGs.")

    # 2b. Collect deduplicated candidate proteins across 259 MAGs
    raw_cands = {}
    missing_faa = 0
    for mid in sorted(mags):
        p = FAA_DIR / f"{mid}.sbl_pfams.faa"
        if not p.exists():
            missing_faa += 1
            continue
        for gid, hdr, seq in read_fasta(p):
            uid = f"{mid}|{gid}"
            if uid not in raw_cands:
                raw_cands[uid] = (mid, gid, seq.rstrip("*"))

    raw_faa = OUT_DIR / "marine_mags_raw_candidates.faa"
    write_fasta(raw_faa, [(uid, v[2]) for uid, v in sorted(raw_cands.items())])
    print(f"[2/5] Collected {len(raw_cands)} unique proteins across {len(mags) - missing_faa} MAGs (missing .faa={missing_faa}).")

    # 2c. Full Pfam-A (--cut_ga) + overlapping-envelope competition
    pfam_dom = OUT_DIR / "marine_mags_pfam_all.domtbl"
    print("[3/5] Running hmmsearch full Pfam-A.hmm (--cut_ga) on marine MAG candidates...")
    run_cmd(f"{HMMSEARCH_BIN} --cpu 16 --cut_ga --domtblout {pfam_dom} {PFAM_ALL_HMM} {raw_faa} >/dev/null")

    acc_to_clan = load_pfam_clans()
    pfam_hits = parse_hmmsearch_domtbl(pfam_dom)

    # 2d. Load geNomad prophage, full-virus, and plasmid annotations
    pv_intervals = defaultdict(list)
    with open(PV_TSV, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            pv_intervals[r["source_seq"]].append((int(r["start"]), int(r["end"])))

    full_viral_contigs = set()
    with open(VS_TSV, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            base = r["seq_name"].split("|")[0]
            topo = r.get("topology", "")
            coords = r.get("coordinates", "NA")
            if topo == "Provirus" or (coords and coords != "NA" and "-" in coords):
                for seg in coords.split(","):
                    m = re.match(r"^(\d+)-(\d+)$", seg.strip())
                    if m:
                        pv_intervals[base].append((int(m.group(1)), int(m.group(2))))
            else:
                full_viral_contigs.add(base)

    plasmid_contigs = set()
    with open(PL_TSV, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            plasmid_contigs.add(r["seq_name"].split("|")[0])

    # Load contig DNA sequences lazily per MAG when verifying CL0013 winners
    contig_cache = {}
    def get_mag_contigs(mid):
        if mid not in contig_cache:
            cpath = FNA_DIR / f"{mid}.sbl_contigs.fna"
            contig_cache[mid] = {cid: dna for cid, _, dna in read_fasta(cpath)} if cpath.exists() else {}
        return contig_cache[mid]

    # 2e. Evaluate Pfam competition + geNomad excision + 1:1 CDS verification + partial-ORF filter
    audit_rows = []
    status_counts = Counter()
    cds_ok = 0
    cds_fail = 0
    prot_out, dom_out, cds_out, dom_cds_out = [], [], [], []

    for uid in sorted(raw_cands):
        mid, gid, prot = raw_cands[uid]
        winners = resolve_overlapping_domains(pfam_hits.get(uid, []), acc_to_clan, OV_THRESHOLD_AA)
        cat, win_acc, best_dom = classify_cl0013(winners)
        arch_str = "+".join(f"{w['pfam_acc']}:{w['pfam_name']}({w['env_from']}-{w['env_to']})" for w in winners)

        m_gid = re.match(r"^(.+)_(\d+)_(\d+)$", gid)
        contig_id = m_gid.group(1) if m_gid else gid
        st = int(m_gid.group(2)) if m_gid else 0
        en = int(m_gid.group(3)) if m_gid else 0
        concat_contig = f"{mid}_{contig_id}"

        dom_start = max(1, best_dom["env_from"] - DOMAIN_PAD_AA) if best_dom else 1
        dom_end = min(len(prot), best_dom["env_to"] + DOMAIN_PAD_AA) if best_dom else len(prot)
        dom_seq = prot[dom_start - 1 : dom_end]
        has_sxxk = bool(re.search(r"S[A-Z]{2}K", dom_seq))

        # Check geNomad prophage/virus vs plasmid/chromosome
        pv_ov = 0
        if concat_contig in full_viral_contigs:
            gn_flag = "DROP_genomad_full_viral_contig"
        elif concat_contig in pv_intervals and st > 0:
            pv_ov = max(interval_overlap(st, en, ps, pe) for ps, pe in pv_intervals[concat_contig])
            gn_flag = "DROP_genomad_prophage_interval" if pv_ov > 0 else (
                "plasmid" if concat_contig in plasmid_contigs else "chromosome"
            )
        else:
            gn_flag = "plasmid" if concat_contig in plasmid_contigs else "chromosome"

        # Verify 1:1 CDS from contig .fna + detect partial/incomplete ORFs
        cslice = None
        strand = 0
        has_start = False
        has_stop = False
        is_edge = False
        cds_verified = False
        if not cat.startswith("drop_") and st > 0:
            cdna = get_mag_contigs(mid).get(contig_id)
            if cdna is not None:
                cslice, strand, has_start, has_stop, is_edge = match_contig_slice_to_prot(cdna, st, en, prot)
                if cslice is not None:
                    cds_verified = True
                    cds_ok += 1
                else:
                    cds_fail += 1
            else:
                cds_fail += 1

        is_partial_orf = is_edge or (cds_verified and not (has_start and has_stop))

        if cat.startswith("drop_"):
            final_status = cat
        elif gn_flag.startswith("DROP_"):
            final_status = gn_flag
        elif not cds_verified:
            final_status = "DROP_unverified_cds"
        elif is_partial_orf:
            final_status = "DROP_partial_orf"
        else:
            final_status = f"KEEP_{cat}_{gn_flag}"

        status_counts[final_status] += 1
        tax = mags[mid]["tax"] or [""] * 7
        audit_rows.append({
            "uid": uid,
            "mag_id": mid,
            "gene_id": gid,
            "contig_id": contig_id,
            "concat_contig": concat_contig,
            "cds_start": st,
            "cds_end": en,
            "cds_strand": strand,
            "seq_len": len(prot),
            "pfam_category": cat,
            "winning_acc": win_acc,
            "win_dom_score": best_dom["dom_score"] if best_dom else 0.0,
            "env_from": best_dom["env_from"] if best_dom else 0,
            "env_to": best_dom["env_to"] if best_dom else 0,
            "dom_slice_start": dom_start,
            "dom_slice_end": dom_end,
            "has_SXXK_in_dom": has_sxxk,
            "domain_architecture": arch_str,
            "genomad_status": gn_flag,
            "genomad_pv_overlap_bp": pv_ov,
            "cds_verified": cds_verified,
            "has_start_codon": has_start,
            "has_stop_codon": has_stop,
            "is_contig_edge": is_edge,
            "is_partial_orf": is_partial_orf,
            "final_status": final_status,
            "gtdb_phylum": tax[1],
            "gtdb_class": tax[2],
            "gtdb_order": tax[3],
            "gtdb_family": tax[4],
            "gtdb_genus": tax[5],
            "gtdb_species": tax[6],
            "mag_completeness": mags[mid]["completeness"],
            "mag_contamination": mags[mid]["contamination"],
        })

        if final_status.startswith("KEEP_"):
            hdr_tag = (
                f"{uid} fam={win_acc} cat={cat} origin=bacterial_MarineMAG replicon={gn_flag} "
                f"phylum={tax[1]} dom={dom_start}..{dom_end} score={best_dom['dom_score']:.1f}"
            )
            prot_out.append((hdr_tag, prot))
            dom_out.append((hdr_tag, dom_seq))
            cds_out.append((hdr_tag, cslice))
            if len(cslice) >= 3 * dom_end:
                dom_cds_out.append((hdr_tag, cslice[(dom_start - 1) * 3 : dom_end * 3]))

    with open(OUT_DIR / "Marine_MAGs_SBL_Audit.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(audit_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(audit_rows)

    write_fasta(OUT_DIR / "marine_mags_sbl_proteins.faa", prot_out)
    write_fasta(OUT_DIR / "marine_mags_sbl_dom.faa", dom_out)
    write_fasta(OUT_DIR / "marine_mags_sbl_cds_verified.fna", cds_out)
    write_fasta(OUT_DIR / "marine_mags_sbl_dom_cds_verified.fna", dom_cds_out)

    lines = [
        "=== 259 Dereplicated Marine Bacterial MAGs Audit Summary ===",
        f"Total raw candidate proteins across {len(mags)} MAGs: {len(raw_cands)}",
        f"CDS translation verification among CL0013 GA non-drop hits: {cds_ok} / {cds_ok + cds_fail} ({100*cds_ok/max(cds_ok+cds_fail,1):.2f}%)",
        f"Total clean non-prophage, complete-ORF bacterial marine MAG SBLs kept: {len(prot_out)}",
        "\nDetailed status breakdown:",
    ] + [f"  {k}: {v}" for k, v in sorted(status_counts.items())]
    rep = "\n".join(lines) + "\n"
    (OUT_DIR / "Marine_MAGs_Summary.txt").write_text(rep, encoding="utf-8")
    print("[5/5] " + rep)


if __name__ == "__main__":
    main()
