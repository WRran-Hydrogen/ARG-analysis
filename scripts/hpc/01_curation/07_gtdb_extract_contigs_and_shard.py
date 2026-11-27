#!/usr/bin/env python3
"""Extract SBL-bearing bacterial contigs from GTDB/RefSeq _genomic.fna.gz files,
rename headers uniquely as >{acc}|{contig_acc}, and pack into 8 balanced bp shards for geNomad.
"""

import csv
import gzip
from collections import defaultdict
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

BDIR = hpc.BASE / "01b_bacterial_ref_db"
HITS_TSV = BDIR / "sbl_contigs" / "bac_gtdb_sbl_hits.tsv"
FNA_DIR = BDIR / "downloads" / "fna"
SHARD_DIR = BDIR / "sbl_contigs" / "shards"
NUM_SHARDS = 8


def main():
    SHARD_DIR.mkdir(parents=True, exist_ok=True)
    targets_by_acc = defaultdict(set)
    with open(HITS_TSV, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            acc = r["acc"]
            cid = r["contig_acc"]
            if cid and cid != "NA":
                targets_by_acc[acc].add(cid)
            # Also include any additional verified paralog contigs
            for loc in r.get("all_verified_loci", "").split(";"):
                if ":" in loc and loc != "NA":
                    c2 = loc.split(":", 1)[0]
                    if c2:
                        targets_by_acc[acc].add(c2)

    total_targets = sum(len(v) for v in targets_by_acc.values())
    print(f"[1/3] Targeting {total_targets:,} unique (genome, contig) entries across {len(targets_by_acc):,} genomes.")

    # Extract contigs directly to an index of (uniq_hdr, seq) sorted by descending length
    extracted = []
    missing_fna = 0
    missing_contigs = 0
    for acc in sorted(targets_by_acc):
        want = targets_by_acc[acc]
        fpath = FNA_DIR / f"{acc}_genomic.fna.gz"
        if not fpath.exists() or fpath.stat().st_size < 100:
            missing_fna += 1
            missing_contigs += len(want)
            continue
        found_here = set()
        cur_cid = None
        keep = False
        buf = []
        try:
            with gzip.open(fpath, "rt", encoding="utf-8", errors="ignore") as gf:
                for line in gf:
                    s = line.strip()
                    if not s:
                        continue
                    if s.startswith(">"):
                        if keep and cur_cid is not None:
                            seq = "".join(buf)
                            extracted.append((f"{acc}|{cur_cid}", seq))
                            found_here.add(cur_cid)
                        cur_cid = s[1:].split()[0]
                        keep = cur_cid in want
                        buf = []
                    elif keep:
                        buf.append(s)
                if keep and cur_cid is not None:
                    seq = "".join(buf)
                    extracted.append((f"{acc}|{cur_cid}", seq))
                    found_here.add(cur_cid)
        except Exception as e:
            print(f"WARN reading {fpath.name}: {e}")
        missing_contigs += len(want - found_here)

    print(f"[2/3] Extracted {len(extracted):,} contigs (missing_fna_files={missing_fna}, missing_contigs={missing_contigs}).")
    if len(extracted) == 0:
        raise RuntimeError("No contigs extracted; aborting sharding.")

    # Greedy longest-first bin packing across NUM_SHARDS shards
    extracted.sort(key=lambda x: (-len(x[1]), x[0]))
    shard_bp = [0] * NUM_SHARDS
    shard_cnt = [0] * NUM_SHARDS
    shard_files = [open(SHARD_DIR / f"chunk_{i}.fna", "w", encoding="utf-8") for i in range(NUM_SHARDS)]
    manifest_rows = []

    try:
        for hdr, seq in extracted:
            idx = min(range(NUM_SHARDS), key=lambda i: (shard_bp[i], shard_cnt[i], i))
            fobj = shard_files[idx]
            fobj.write(f">{hdr}\n")
            for k in range(0, len(seq), 80):
                fobj.write(seq[k : k + 80] + "\n")
            shard_bp[idx] += len(seq)
            shard_cnt[idx] += 1
            manifest_rows.append((idx, hdr, len(seq)))
    finally:
        for fobj in shard_files:
            fobj.close()

    with open(SHARD_DIR / "shard_manifest.tsv", "w", newline="", encoding="utf-8") as mf:
        w = csv.writer(mf, delimiter="\t")
        w.writerow(["shard_id", "seq_id", "length_bp"])
        w.writerows(manifest_rows)

    print("[3/3] Balanced shard statistics:")
    for i in range(NUM_SHARDS):
        print(f"  chunk_{i}.fna: {shard_cnt[i]:,} contigs | {shard_bp[i]:,} bp ({shard_bp[i]/1e6:.2f} Mb)")


if __name__ == "__main__":
    main()
