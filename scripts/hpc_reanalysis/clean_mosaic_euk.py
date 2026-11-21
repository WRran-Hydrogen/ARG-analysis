#!/usr/bin/env python3
"""Purge MOSAiC viral SBL candidates whose contigs overlap euk_info.tsv eukaryotic bins."""

import csv
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import hpc_config as hpc

WORK = hpc.BASE / "01_curated_datasets" / "mosaic_audit"
AUDIT_TSV = WORK / "MOSAiC_SBL_Audit.tsv"
CAND_FAA = WORK / "mosaic_cl0013_candidates.faa"
EUK_INFO = hpc.RAW_INPUT_DIR / "euk_info.tsv"


def read_fasta(path):
    recs = {}
    cur_id, buf = None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if cur_id is not None:
                    recs[cur_id] = "".join(buf).rstrip("*")
                cur_id = s[1:].strip().split()[0]
                buf = []
            else:
                buf.append(s)
        if cur_id is not None:
            recs[cur_id] = "".join(buf).rstrip("*")
    return recs


def write_fasta(path, items):
    with open(path, "w", encoding="utf-8") as f:
        for hdr, seq in items:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def main():
    with open(AUDIT_TSV, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    keep_rows = [r for r in rows if r["pfam_category"].startswith(("focal_", "aux_"))]
    contigs = sorted({r["inferred_contig"] for r in keep_rows})
    seqs = read_fasta(CAND_FAA)

    euk_hits = defaultdict(set)
    with tempfile.NamedTemporaryFile("w", delete=False) as tf:
        for c in contigs:
            tf.write(c + "_\n")
        pat_file = tf.name

    try:
        proc = subprocess.run(["grep", "-F", "-f", pat_file, str(EUK_INFO)], capture_output=True, text=True)
        if proc.returncode not in (0, 1):
            raise RuntimeError(f"grep failed ({proc.returncode}): {proc.stderr}")
        for line in proc.stdout.splitlines():
            p = line.split("\t")
            if len(p) < 2:
                continue
            gid, genome = p[0], p[1]
            for c in contigs:
                if gid.startswith(c + "_"):
                    euk_hits[c].add(genome)
    finally:
        Path(pat_file).unlink(missing_ok=True)

    clean_rows = []
    dropped = []
    surv_counts = Counter()
    surv_by_cat_bin = Counter()
    full_out = defaultdict(list)
    dom_out = defaultdict(list)

    for r in keep_rows:
        c = r["inferred_contig"]
        genomes = sorted(euk_hits.get(c, set()))
        in_euk = len(genomes) > 0
        r_out = dict(r)
        r_out["in_euk_info"] = str(in_euk)
        r_out["euk_genomes"] = ";".join(genomes) if genomes else "NONE"
        r_out["final_mosaic_status"] = "DROP_euk_info_overlap" if in_euk else "KEEP_MOSAiC_CLEAN"
        clean_rows.append(r_out)

        cat = r["pfam_category"]
        if in_euk:
            dropped.append((r["seq_id"], cat, r["vir_bin"], ";".join(genomes)))
            surv_counts[("DROP_euk", cat)] += 1
        else:
            surv_counts[("KEEP", cat)] += 1
            surv_by_cat_bin[(cat, r["bin_type"], r["vir_tax_suffix"])] += 1
            prot = seqs[r["seq_id"]]
            st, en = int(r["dom_slice_start"]), int(r["dom_slice_end"])
            dom_seq = prot[st - 1 : en]
            hdr = f"{r['seq_id']} fam={r['winning_acc']} source=MOSAiC bin={r['vir_bin']} dom={st}..{en} score={float(r['win_dom_score']):.1f}"
            full_out[cat].append((hdr, prot))
            dom_out[cat].append((hdr, dom_seq))

    with open(WORK / "MOSAiC_SBL_Clean_Audit.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(clean_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(clean_rows)

    write_fasta(WORK / "mosaic_clean_PF00144_full.faa", full_out.get("focal_PF00144", []))
    write_fasta(WORK / "mosaic_clean_PF00144_dom.faa", dom_out.get("focal_PF00144", []))
    write_fasta(WORK / "mosaic_clean_PF00905_full.faa", full_out.get("focal_PF00905", []))
    write_fasta(WORK / "mosaic_clean_PF00905_dom.faa", dom_out.get("focal_PF00905", []))
    write_fasta(WORK / "mosaic_clean_aux_PF00768_full.faa", full_out.get("aux_PF00768", []))
    write_fasta(WORK / "mosaic_clean_aux_PF00768_dom.faa", dom_out.get("aux_PF00768", []))

    lines = [
        "=== MOSAiC Eukaryotic Overlap Purge Summary ===",
        f"Input passing CL0013 focal/aux candidates: {len(keep_rows)}",
        f"Purged (in_euk_info == True): {len(dropped)}",
    ]
    for sid, cat, vbin, eg in dropped:
        lines.append(f"  DROPPED {sid} | {cat} | vir_bin={vbin} | euk={eg}")
    lines.append(f"\nSurviving Clean MOSAiC Viral SBLs: {len(keep_rows) - len(dropped)}")
    for (st_k, cat), cnt in sorted(surv_counts.items()):
        lines.append(f"  {st_k:8s} {cat}: {cnt}")
    lines.append("\nSurviving Clean Breakdown (category, bin_type, vir_tax_suffix):")
    for k, cnt in sorted(surv_by_cat_bin.items()):
        lines.append(f"  {k}: {cnt}")

    summary = "\n".join(lines) + "\n"
    (WORK / "MOSAiC_Clean_Summary.txt").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
