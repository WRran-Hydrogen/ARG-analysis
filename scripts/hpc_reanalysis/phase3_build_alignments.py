#!/usr/bin/env python3
"""Phase 3, step 1: build the alignments and the IQ-TREE task list.

Inputs (Stage 1 outputs in $SBL_HPC_BASE/01_final_curated/):
  {tag}_combined_tree*_dom.faa   viral nr85 + bacterial tree set (primary = rep1, plus
                                 rep2..rep5, phylo, nopro, cov60, cov70 sensitivity sets)
  {tag}_sister_pre_outgroup_dom.faa   Sister-PRE outgroup (PF13354/PF00768 trees: PF00144+PF00905 only)

Alignment tracks (pre-registered in docs/METHODS_AND_PIPELINE.md, Section 4.1):
  hmm    hmmalign against the family's official Pfam HMM; insert (lowercase) states removed,
         so every row has exactly LENG match columns. Primary track for single families.
  mafft  MAFFT --auto + ClipKIT gap trimming (-m gappy -g 0.5). Sensitivity track; primary for the joint
         PF13354_PF00768 tree, which has no single family HMM.
Every "_out" alignment includes the Sister-PRE outgroup; "_in" (ingroup-only) alignments are the
same alignment with the outgroup rows and then all-gap columns removed (no re-alignment), for
outgroup-free rooting (MAD / MinVar) in Phase 4.

Outputs ($SBL_HPC_PHASE3_DIR):
  alignments/{tag}/{set}.{track}_{out|in}.fa, alignments/{tag}/outgroup_taxa.txt
  alignment_qc.tsv   rows, columns, mean occupancy, SXXK / KTG column lock-in rates
  tree_tasks.tsv     wave 1 (ModelFinder) and wave 2 (model fixed from the matching wave-1 tree)
"""

import csv
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import hpc_config as hpc

IN_DIR = hpc.BASE / "01_final_curated"
INPUTS = hpc.BASE / "inputs"
OUT_DIR = hpc.PHASE3_DIR
ALN_DIR = OUT_DIR / "alignments"

SINGLE_FAMS = ("PF00144", "PF00905", "PF13354", "PF00768")
JOINT = "PF13354_PF00768"
SENS_SETS = ("rep2", "rep3", "rep4", "rep5", "phylo", "nopro", "cov60", "cov70")
THREADS = 16
CLIPKIT_GAP_THRESHOLD = 0.5


def run(cmd, stdout_path=None):
    if stdout_path is None:
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    else:
        with open(stdout_path, "w", encoding="utf-8") as fh:
            p = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd}\n{p.stderr}")


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path, recs):
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def hmm_leng(hmm_path):
    with open(hmm_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("LENG "):
                return int(line.split()[1])
    raise ValueError(f"No LENG in {hmm_path}")


def hmmalign_match_states(fasta_in, fam, fasta_out):
    """hmmalign (A2M) -> keep match states only (uppercase residues and '-')."""
    hmm = INPUTS / f"{fam}_official.hmm"
    a2m = fasta_out.with_suffix(".a2m")
    run(f"{hpc.HMMALIGN_BIN} --outformat A2M {hmm} {fasta_in}", stdout_path=a2m)
    leng = hmm_leng(hmm)
    out = []
    for name, row in read_fasta(a2m):
        m = "".join(c for c in row if c.isupper() or c == "-")
        if len(m) != leng:
            raise RuntimeError(f"{fasta_out.name}: {name} has {len(m)} match states, HMM LENG={leng}")
        out.append((name, m))
    write_fasta(fasta_out, out)
    a2m.unlink()
    return out


def mafft_clipkit(fasta_in, fasta_out):
    """MAFFT --auto, then ClipKIT gap-only trimming (drop columns with >= CLIPKIT_GAP_THRESHOLD gaps).

    kpic-smart-gap was rejected: on these insertion-rich domain slices its data-driven threshold
    kept >70% gap columns (occupancy 0.16-0.35), and kpic modes also remove singleton sites,
    which shortens terminal branches that MAD/MinVar rooting relies on. The untrimmed MAFFT
    alignment is kept next to the trimmed one.
    """
    raw = fasta_out.with_suffix(".mafft_raw.fa")
    run(f"{hpc.MAFFT_BIN} --auto --anysymbol --thread {THREADS} {fasta_in}", stdout_path=raw)
    run(f"{hpc.CLIPKIT_BIN} {raw} -m gappy -g {CLIPKIT_GAP_THRESHOLD} -o {fasta_out}")
    return read_fasta(fasta_out)


def drop_rows_and_empty_columns(recs, drop):
    kept = [(n, s) for n, s in recs if n not in drop]
    if not kept:
        return kept
    keep_cols = [i for i in range(len(kept[0][1])) if any(s[i] not in "-.Xx?" for _, s in kept)]
    return [(n, "".join(s[i] for i in keep_cols)) for n, s in kept]


def motif_lock_rate(recs, pattern):
    """Fraction of rows with an occurrence of `pattern` starting at the modal alignment column.

    All (overlapping) occurrences per row are considered, so a chance match elsewhere in the
    sequence does not hide a correctly aligned catalytic motif.
    """
    rx = re.compile(f"(?=({pattern}))")
    row_cols = []
    for _, row in recs:
        pos = [i for i, c in enumerate(row) if c not in "-."]
        ungapped = "".join(row[i] for i in pos).upper()
        row_cols.append({pos[m.start()] for m in rx.finditer(ungapped)})
    tally = Counter(c for cols in row_cols for c in cols)
    if not tally:
        return 0.0
    modal = tally.most_common(1)[0][0]
    return round(sum(modal in cols for cols in row_cols) / len(recs), 4)


def qc_row(tag, set_name, track, recs, outgroup):
    n_cols = len(recs[0][1]) if recs else 0
    occ = sum(sum(c not in "-." for c in s) for _, s in recs) / max(len(recs) * n_cols, 1)
    ingroup = [(n, s) for n, s in recs if n not in outgroup]
    return {
        "tag": tag, "set": set_name, "track": track, "n_rows": len(recs), "n_outgroup": len(recs) - len(ingroup),
        "n_columns": n_cols, "mean_occupancy": round(occ, 4),
        "SXXK_lock_in_rate": motif_lock_rate(ingroup, r"S..K"),
        "KTG_lock_in_rate": motif_lock_rate(ingroup, r"[KH][TS]G"),
    }


def main():
    ALN_DIR.mkdir(parents=True, exist_ok=True)
    qc, tasks = [], []

    def add_task(wave, name, aln, model_from=""):
        tasks.append({"wave": wave, "name": name, "alignment": str(aln), "model_from": model_from})

    for tag in SINGLE_FAMS + (JOINT,):
        tdir = ALN_DIR / tag
        tdir.mkdir(exist_ok=True)
        outgroup_recs = read_fasta(IN_DIR / f"{tag}_sister_pre_outgroup_dom.faa")
        outgroup = {n for n, _ in outgroup_recs}
        (tdir / "outgroup_taxa.txt").write_text("\n".join(sorted(outgroup)) + "\n", encoding="utf-8")

        sets = {"primary": IN_DIR / f"{tag}_combined_tree_dom.faa"}
        for s in SENS_SETS:
            p = IN_DIR / f"{tag}_combined_tree_{s}_dom.faa"
            if p.exists():
                sets[s] = p
        primary_track = "hmm" if tag != JOINT else "mafft"

        for set_name, src in sets.items():
            seqs = read_fasta(src) + outgroup_recs
            dup = [n for n, c in Counter(n for n, _ in seqs).items() if c > 1]
            if dup:
                raise RuntimeError(f"{tag}/{set_name}: duplicate taxa between ingroup and outgroup: {dup[:3]}")
            unaligned = tdir / f"{set_name}.unaligned.faa"
            write_fasta(unaligned, seqs)
            # primary set: both tracks; sensitivity sets: primary track only
            tracks = ("hmm", "mafft") if (set_name == "primary" and tag != JOINT) else (primary_track,)
            for track in tracks:
                out_fa = tdir / f"{set_name}.{track}_out.fa"
                recs = hmmalign_match_states(unaligned, tag, out_fa) if track == "hmm" else mafft_clipkit(unaligned, out_fa)
                qc.append(qc_row(tag, set_name, f"{track}_out", recs, outgroup))
                if set_name == "primary":
                    in_fa = tdir / f"{set_name}.{track}_in.fa"
                    in_recs = drop_rows_and_empty_columns(recs, outgroup)
                    write_fasta(in_fa, in_recs)
                    qc.append(qc_row(tag, set_name, f"{track}_in", in_recs, outgroup))
            unaligned.unlink()

        ref = f"{tag}.primary.{primary_track}_out"
        add_task(1, ref, tdir / f"primary.{primary_track}_out.fa")
        if tag != JOINT:
            add_task(1, f"{tag}.primary.mafft_out", tdir / "primary.mafft_out.fa")
        add_task(2, f"{tag}.primary.{primary_track}_in", tdir / f"primary.{primary_track}_in.fa", model_from=ref)
        for set_name in sets:
            if set_name != "primary":
                add_task(2, f"{tag}.{set_name}.{primary_track}_out", tdir / f"{set_name}.{primary_track}_out.fa", model_from=ref)
        print(f"{tag}: {len(sets)} sets aligned, outgroup={len(outgroup)}", flush=True)

    with open(OUT_DIR / "alignment_qc.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(qc[0]), delimiter="\t")
        w.writeheader()
        w.writerows(qc)
    for wave in (1, 2):
        with open(OUT_DIR / f"tree_tasks_wave{wave}.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["wave", "name", "alignment", "model_from"], delimiter="\t")
            w.writeheader()
            w.writerows(t for t in tasks if t["wave"] == wave)
    n1 = sum(t["wave"] == 1 for t in tasks)
    print(f"Wrote alignment_qc.tsv, tree_tasks_wave1.tsv ({n1} trees), tree_tasks_wave2.tsv ({len(tasks) - n1} trees)")


if __name__ == "__main__":
    sys.exit(main())
