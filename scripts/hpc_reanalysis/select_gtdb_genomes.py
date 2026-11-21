#!/usr/bin/env python3
"""Select 3,000 phylogenetically balanced, high-quality GTDB/RefSeq bacterial species representatives."""

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

import hpc_config as hpc

BDIR = hpc.BASE / "01b_bacterial_ref_db"
GTDB_GZ = BDIR / "metadata" / "bac120_metadata.tsv.gz"
REFSEQ_TXT = BDIR / "metadata" / "refseq_assembly_summary.txt"
OUT_META = BDIR / "metadata" / "selected_3000_genomes.tsv"
OUT_URLS = BDIR / "metadata" / "urls_faa_cds.txt"
OUT_SUM = BDIR / "metadata" / "selection_summary.txt"

FAA_DIR = BDIR / "downloads" / "faa"
CDS_DIR = BDIR / "downloads" / "cds"
TARGET_N = 3000

ASM_RANK = {
    "Complete Genome": 0,
    "Chromosome": 1,
    "Scaffold": 2,
    "Contig": 3,
}

PRIORITY_PHYLA = [
    "p__Pseudomonadota",
    "p__Bacillota",
    "p__Bacillota_A",
    "p__Actinomycetota",
    "p__Bacteroidota",
    "p__Cyanobacteriota",
    "p__Planctomycetota",
    "p__Verrucomicrobiota",
    "p__Chloroflexota",
    "p__Acidobacteriota",
    "p__Desulfobacterota",
    "p__Campylobacterota",
    "p__Spirochaetota",
    "p__Myxococcota",
    "p__Deinococcota",
    "p__Fusobacteriota",
]


def safe_float(val, default=-1.0):
    if not val or val.lower() in ("na", "none", "nan"):
        return None
    try:
        return float(val)
    except ValueError:
        return None


def parse_taxonomy(tax_str):
    parts = [x.strip() for x in tax_str.split(";")]
    while len(parts) < 7:
        parts.append("")
    return tuple(parts[:7])


def main():
    print("[1/4] Parsing NCBI RefSeq bacteria assembly_summary.txt...")
    refseq = {}
    with open(REFSEQ_TXT, "r", encoding="utf-8", errors="ignore") as f:
        header = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#assembly_accession"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header):
                continue
            acc = p[0]
            ftp_path = p[header.index("ftp_path")].strip()
            if not (ftp_path.startswith("https://") or ftp_path.startswith("ftp://")):
                continue
            if ftp_path.startswith("ftp://"):
                ftp_path = "https://" + ftp_path[6:]
            refseq[acc] = {
                "assembly_level": p[header.index("assembly_level")].strip(),
                "refseq_category": p[header.index("refseq_category")].strip(),
                "organism_name": p[header.index("organism_name")].strip(),
                "version_status": p[header.index("version_status")].strip(),
                "ftp_path": ftp_path.rstrip("/"),
            }
    print(f"  Loaded {len(refseq):,} RefSeq bacterial assemblies with valid URLs.")

    print("[2/4] Filtering GTDB bac120_metadata.tsv.gz (gtdb_representative=='t', comp>=95, contam<=2)...")
    candidates = []
    stats = Counter()
    with gzip.open(GTDB_GZ, "rt", encoding="utf-8", errors="ignore") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx_acc = hdr.index("accession")
        idx_rep = hdr.index("gtdb_representative")
        idx_tax = hdr.index("gtdb_taxonomy")
        idx_c2 = hdr.index("checkm2_completeness") if "checkm2_completeness" in hdr else -1
        idx_k2 = hdr.index("checkm2_contamination") if "checkm2_contamination" in hdr else -1
        idx_c1 = hdr.index("checkm_completeness")
        idx_k1 = hdr.index("checkm_contamination")

        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= idx_tax:
                continue
            if p[idx_rep] != "t":
                continue
            stats["gtdb_reps"] += 1
            raw_acc = p[idx_acc]
            if not raw_acc.startswith("RS_GCF_"):
                continue
            acc = raw_acc[3:]  # strip RS_
            stats["rs_gcf_reps"] += 1
            rinfo = refseq.get(acc)
            if rinfo is None:
                stats["missing_in_refseq_summary"] += 1
                continue

            comp = safe_float(p[idx_c2]) if idx_c2 >= 0 else None
            contam = safe_float(p[idx_k2]) if idx_k2 >= 0 else None
            if comp is None or contam is None:
                comp = safe_float(p[idx_c1])
                contam = safe_float(p[idx_k1])
            if comp is None or contam is None:
                stats["missing_checkm"] += 1
                continue
            if comp < 95.0 or contam > 2.0:
                stats["failed_qc"] += 1
                continue

            dom, phy, cls, ord_, fam, gen, sp = parse_taxonomy(p[idx_tax])
            if dom != "d__Bacteria":
                continue
            stats["passed_qc"] += 1

            asm_lvl = rinfo["assembly_level"]
            base_url = rinfo["ftp_path"]
            asm_name = base_url.split("/")[-1]
            candidates.append({
                "accession": acc,
                "domain": dom,
                "phylum": phy or "p__unclassified",
                "class": cls or "c__unclassified",
                "order": ord_ or "o__unclassified",
                "family": fam or f"f__unclassified({ord_})",
                "genus": gen or f"g__unclassified({fam})",
                "species": sp or acc,
                "organism_name": rinfo["organism_name"],
                "assembly_level": asm_lvl,
                "refseq_category": rinfo["refseq_category"],
                "completeness": comp,
                "contamination": contam,
                "sort_key": (ASM_RANK.get(asm_lvl, 4), -(comp - 5.0 * contam), acc),
                "faa_url": f"{base_url}/{asm_name}_protein.faa.gz",
                "cds_url": f"{base_url}/{asm_name}_cds_from_genomic.fna.gz",
                "fna_url": f"{base_url}/{asm_name}_genomic.fna.gz",
            })

    print("  Filter stats:", dict(stats))
    if len(candidates) < 1500:
        raise RuntimeError(f"Too few high-quality RefSeq GTDB representatives: {len(candidates)} < 1500")

    print("[3/4] Running 2-pass phylogenetically balanced selection (target N=3,000)...")
    # Pass 1: #1 best-ranked genome per GTDB (phylum, family)
    by_family = defaultdict(list)
    for c in candidates:
        by_family[(c["phylum"], c["family"])].append(c)
    for k in by_family:
        by_family[k].sort(key=lambda x: x["sort_key"])

    selected = []
    selected_acc = set()
    used_genera = set()

    pass1_pool = [lst[0] for lst in by_family.values()]
    pass1_pool.sort(key=lambda x: x["sort_key"])
    for c in pass1_pool:
        if len(selected) >= TARGET_N:
            break
        c["select_pass"] = "Pass1_FamilyBackbone"
        selected.append(c)
        selected_acc.add(c["accession"])
        used_genera.add((c["phylum"], c["genus"]))

    n_pass1 = len(selected)

    # Pass 2: 1 best-ranked genome per remaining unselected GTDB genus, round-robin across phyla
    by_genus = defaultdict(list)
    for c in candidates:
        if c["accession"] in selected_acc:
            continue
        gk = (c["phylum"], c["genus"])
        if gk in used_genera:
            continue
        by_genus[gk].append(c)

    genus_queues = defaultdict(list)  # phylum -> sorted list of best genus reps
    for gk, lst in by_genus.items():
        lst.sort(key=lambda x: x["sort_key"])
        genus_queues[gk[0]].append(lst[0])
    for phy in genus_queues:
        genus_queues[phy].sort(key=lambda x: x["sort_key"])

    # Phylum round-robin order: PRIORITY_PHYLA first (with weight 2 per cycle for major phage-host phyla), then remaining phyla
    other_phyla = sorted(p for p in genus_queues if p not in PRIORITY_PHYLA)
    cycle_order = []
    for p in PRIORITY_PHYLA:
        if p in genus_queues:
            cycle_order.extend([p, p])
    for p in other_phyla:
        cycle_order.append(p)

    q_idx = defaultdict(int)
    while len(selected) < TARGET_N:
        added_this_round = 0
        for phy in cycle_order:
            if len(selected) >= TARGET_N:
                break
            lst = genus_queues.get(phy, [])
            i = q_idx[phy]
            while i < len(lst) and lst[i]["accession"] in selected_acc:
                i += 1
            if i < len(lst):
                c = lst[i]
                q_idx[phy] = i + 1
                c["select_pass"] = "Pass2_GenusHostFill"
                selected.append(c)
                selected_acc.add(c["accession"])
                used_genera.add((c["phylum"], c["genus"]))
                added_this_round += 1
        if added_this_round == 0:
            break

    print(f"[4/4] Selected {len(selected)} genomes (Pass1={n_pass1}, Pass2={len(selected) - n_pass1}). Writing outputs...")
    cols = [
        "accession", "select_pass", "phylum", "class", "order", "family", "genus", "species",
        "organism_name", "assembly_level", "refseq_category", "completeness", "contamination",
        "faa_url", "cds_url", "fna_url",
    ]
    with open(OUT_META, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(selected)

    with open(OUT_URLS, "w", encoding="utf-8") as f:
        for c in selected:
            acc = c["accession"]
            f.write(f"{c['faa_url']}\n  dir={FAA_DIR}\n  out={acc}_protein.faa.gz\n")
            f.write(f"{c['cds_url']}\n  dir={CDS_DIR}\n  out={acc}_cds_from_genomic.fna.gz\n")

    asm_counts = Counter(c["assembly_level"] for c in selected)
    phy_counts = Counter(c["phylum"] for c in selected)
    u_phy = len({c["phylum"] for c in selected})
    u_cls = len({c["class"] for c in selected})
    u_ord = len({c["order"] for c in selected})
    u_fam = len({c["family"] for c in selected})
    u_gen = len({c["genus"] for c in selected})

    sum_lines = [
        "=== GTDB/RefSeq Balanced 3,000 Bacterial Genomes Selection Summary ===",
        f"Total GTDB bacterial species reps in RS_GCF: {stats['rs_gcf_reps']:,}",
        f"Passed QC (comp >= 95%, contam <= 2%, valid RefSeq URL): {stats['passed_qc']:,}",
        f"Total selected genomes: {len(selected)} (Pass1 FamilyBackbone={n_pass1}, Pass2 GenusFill={len(selected)-n_pass1})",
        f"Taxonomic breadth: {u_phy} phyla | {u_cls} classes | {u_ord} orders | {u_fam} families | {u_gen} genera",
        f"Assembly levels: {dict(asm_counts)}",
        "\nTop 25 Bacterial Phyla Represented:",
    ]
    for phy, cnt in phy_counts.most_common(25):
        sum_lines.append(f"  {phy:30s}: {cnt}")
    rep = "\n".join(sum_lines) + "\n"
    OUT_SUM.write_text(rep, encoding="utf-8")
    print(rep)


if __name__ == "__main__":
    main()
