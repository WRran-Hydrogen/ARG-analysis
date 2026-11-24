#!/usr/bin/env python3
"""Build CL0013 superfamily + LOO alignments and PMSF C20 task list.

Part A: Builds balanced CL0013 superfamily alignments (5 PRE families: PF00144, PF00905,
PF13354, PF00768, PF02113; plus 3 non-PRE biological outgroup families: PF04960, PF06089,
PF27310) and derives Leave-One-Outgroup-Family-Out (LOO) alignments for outgroup-rooting
validation. Writes $SBL_HPC_PHASE3_DIR/tree_tasks_cl0013.tsv.

Part B: Writes $SBL_HPC_PHASE3_DIR/tree_tasks_pmsf.tsv for the 5 primary trees under
the site-heterogeneous PMSF (LG+C20+F+G) model using finished Wave-1 ML guide trees.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

csv.field_size_limit(10_000_000)

IN_DIR = hpc.BASE / "01_final_curated"
OUT_DIR = hpc.PHASE3_DIR
CL_ALN_DIR = OUT_DIR / "alignments" / "CL0013"

PRE_FOUR = ("PF00144", "PF00905", "PF13354", "PF00768")
OUTGROUP_FAMS = ("PF04960", "PF06089", "PF27310")
THREADS = 16
CLIPKIT_GAP_THRESHOLD = 0.5


def run(cmd: str, stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    else:
        with open(stdout_path, "w", encoding="utf-8") as fh:
            p = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd}\n{p.stderr}")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, str]] = []
    name: str | None = None
    buf: list[str] = []
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


def write_fasta(path: Path, recs: list[tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def evenly_spaced(items: list, k: int) -> list:
    if k <= 0 or not items:
        return []
    if len(items) <= k:
        return list(items)
    if k == 1:
        return [items[0]]
    idxs = sorted({round(i * (len(items) - 1) / (k - 1)) for i in range(k)})
    return [items[i] for i in idxs]


def drop_rows_and_empty_columns(
    recs: list[tuple[str, str]], drop: set[str]
) -> list[tuple[str, str]]:
    kept = [(n, s) for n, s in recs if n not in drop]
    if not kept:
        return kept
    keep_cols = [
        i for i in range(len(kept[0][1]))
        if any(s[i] not in "-.Xx?" for _, s in kept)
    ]
    return [(n, "".join(s[i] for i in keep_cols)) for n, s in kept]


def load_catalogs() -> tuple[dict[str, dict], dict[str, dict]]:
    v_meta: dict[str, dict] = {}
    with open(IN_DIR / "Final_Viral_SBL_Catalog.tsv", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            v_meta[row["uid"]] = row
    b_meta: dict[str, dict] = {}
    with open(IN_DIR / "Final_Bacterial_TripleSource_Catalog.tsv", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            b_meta[row["uid"]] = row
    return v_meta, b_meta


def build_cl0013_alignments(v_meta: dict[str, dict], b_meta: dict[str, dict]) -> list[dict]:
    CL_ALN_DIR.mkdir(parents=True, exist_ok=True)

    selected_seqs: dict[str, str] = {}
    meta_rows: dict[str, dict] = {}

    def register(uid: str, seq: str, fam: str, origin: str, role: str, f_tier: str) -> None:
        if uid not in selected_seqs:
            selected_seqs[uid] = seq
            meta_rows[uid] = {
                "uid": uid,
                "cl0013_family": fam,
                "origin": origin,
                "role": role,
                "functional_tier": f_tier,
            }

    # 1. Include all Sister-PRE outgroup taxa (so every primary outgroup taxon is anchored in CL0013)
    sister_uids_by_fam: dict[str, list[tuple[str, str]]] = {f: [] for f in PRE_FOUR + ("PF02113",)}
    for fam in PRE_FOUR:
        for uid, seq in read_fasta(IN_DIR / f"{fam}_sister_pre_outgroup_dom.faa"):
            info = b_meta.get(uid, {})
            true_fam = info.get("winning_acc", "")
            if true_fam in sister_uids_by_fam:
                sister_uids_by_fam[true_fam].append((uid, seq))
            register(
                uid,
                seq,
                true_fam or "Sister_PRE",
                "Bacterial",
                "PRE_ingroup",
                info.get("functional_tier", "PBP_or_PRE_Homolog"),
            )

    # 2. Balanced bacterial & viral reps from PF00144, PF00905, PF13354, PF00768
    for fam in PRE_FOUR:
        phylo_recs = read_fasta(IN_DIR / f"{fam}_combined_tree_phylo_dom.faa")
        v_nr85_recs = read_fasta(IN_DIR / f"{fam}_viral_nr85_dom.faa")
        v_nr85_map = dict(v_nr85_recs)

        existing_bac = [
            uid for uid, row in meta_rows.items()
            if row["cl0013_family"] == fam and row["origin"] == "Bacterial"
        ]
        bac_cands = [
            (uid, seq) for uid, seq in phylo_recs
            if uid in b_meta and uid not in v_meta and uid not in v_nr85_map and uid not in selected_seqs
        ]
        need_bac = max(0, 60 - len(existing_bac))
        for uid, seq in evenly_spaced(bac_cands, need_bac):
            info = b_meta.get(uid, {})
            register(
                uid,
                seq,
                fam,
                "Bacterial",
                "PRE_ingroup",
                info.get("functional_tier", "PBP_or_PRE_Homolog"),
            )

        # Up to 30 viral nr85 reps per family: prioritize hydrolytic SBLs + Tier 1 viral
        v_sbl_tier1: list[tuple[str, str, int, int]] = []
        v_other: list[tuple[str, str]] = []
        for uid, seq in v_nr85_recs:
            info = v_meta.get(uid, {})
            is_sbl = info.get("functional_tier") == "Putative_Hydrolytic_BetaLactamase"
            is_t1 = info.get("viral_origin_tier") == "Tier1_HighConfidence_Viral"
            if is_sbl or is_t1:
                v_sbl_tier1.append((uid, seq, 0 if is_sbl else 1, 0 if is_t1 else 1))
            else:
                v_other.append((uid, seq))
        v_sbl_tier1.sort(key=lambda x: (x[2], x[3], x[0]))
        if len(v_sbl_tier1) >= 30:
            # Keep all hydrolytic SBLs (up to 30) and evenly sample Tier 1 remainder
            sbl_only = [(u, s) for u, s, r0, _ in v_sbl_tier1 if r0 == 0]
            t1_only = [(u, s) for u, s, r0, _ in v_sbl_tier1 if r0 == 1]
            v_picked = sbl_only[:30] + evenly_spaced(t1_only, max(0, 30 - len(sbl_only)))
        else:
            v_picked = [(u, s) for u, s, _, _ in v_sbl_tier1] + evenly_spaced(
                v_other, max(0, 30 - len(v_sbl_tier1))
            )
        for uid, seq in v_picked:
            info = v_meta.get(uid, {})
            register(
                uid,
                seq,
                fam,
                "Viral",
                "PRE_ingroup",
                info.get("functional_tier", "PBP_or_PRE_Homolog"),
            )

    # 3. Auxiliary PF02113 (first 40 bacterial representatives; deduplicated against Sister-PRE)
    existing_2113 = sum(1 for r in meta_rows.values() if r["cl0013_family"] == "PF02113")
    pf02113_recs = read_fasta(IN_DIR / "CL0013_aux_PF02113_bacterial_rep_dom.faa")
    added_2113 = 0
    for uid, seq in pf02113_recs:
        if existing_2113 + added_2113 >= 40:
            break
        if uid in selected_seqs:
            continue
        info = b_meta.get(uid, {})
        register(
            uid,
            seq,
            "PF02113",
            "Bacterial",
            "PRE_ingroup",
            info.get("functional_tier", "PBP_or_PRE_Homolog"),
        )
        added_2113 += 1

    # 4. Non-PRE biological outgroup families (PF04960, PF06089, PF27310)
    outgroup_by_fam: dict[str, list[str]] = {}
    for og_fam in OUTGROUP_FAMS:
        og_recs = read_fasta(IN_DIR / f"CL0013_outgroup_{og_fam}_dom.faa")
        outgroup_by_fam[og_fam] = []
        for uid, seq in og_recs:
            if uid in selected_seqs:
                raise RuntimeError(f"Outgroup {og_fam} sequence {uid} collides with ingroup")
            info = b_meta.get(uid, {})
            register(
                uid,
                seq,
                og_fam,
                "Bacterial",
                "NonPRE_outgroup",
                info.get("functional_tier", "NonPRE_Outgroup"),
            )
            outgroup_by_fam[og_fam].append(uid)

    # Write taxa metadata table
    meta_path = CL_ALN_DIR / "cl0013_taxa_metadata.tsv"
    with open(meta_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["uid", "cl0013_family", "origin", "role", "functional_tier"],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        for uid in selected_seqs:
            w.writerow(meta_rows[uid])

    fam_counts = Counter(r["cl0013_family"] for r in meta_rows.values())
    print(f"[CL0013] Selected {len(selected_seqs)} unique taxa: {dict(sorted(fam_counts.items()))}")

    # Align once with MAFFT --auto + ClipKIT -m gappy -g 0.5
    unaligned_fa = CL_ALN_DIR / "all3_out.unaligned.faa"
    raw_fa = CL_ALN_DIR / "all3_out.mafft_raw.fa"
    all3_fa = CL_ALN_DIR / "all3_out.mafft_out.fa"
    write_fasta(unaligned_fa, list(selected_seqs.items()))
    print(f"[CL0013] Running MAFFT --auto + ClipKIT on {len(selected_seqs)} sequences...")
    run(f"{hpc.MAFFT_BIN} --auto --anysymbol --thread {THREADS} {unaligned_fa}", stdout_path=raw_fa)
    run(f"{hpc.CLIPKIT_BIN} {raw_fa} -m gappy -g {CLIPKIT_GAP_THRESHOLD} -o {all3_fa}")
    unaligned_fa.unlink()

    all3_recs = read_fasta(all3_fa)
    all3_og = [u for f in OUTGROUP_FAMS for u in outgroup_by_fam[f]]
    (CL_ALN_DIR / "outgroup_taxa_all3.txt").write_text("\n".join(sorted(all3_og)) + "\n", encoding="utf-8")
    (CL_ALN_DIR / "outgroup_taxa.txt").write_text("\n".join(sorted(all3_og)) + "\n", encoding="utf-8")

    tasks = [
        {"wave": 1, "name": "CL0013.all3_out.mafft_out", "alignment": str(all3_fa), "model_from": ""},
    ]

    # Derive the 3 Leave-One-Outgroup-Family-Out (LOO) alignments
    for drop_fam in OUTGROUP_FAMS:
        drop_set = set(outgroup_by_fam[drop_fam])
        kept_og = [u for f in OUTGROUP_FAMS if f != drop_fam for u in outgroup_by_fam[f]]
        loo_fa = CL_ALN_DIR / f"loo_no_{drop_fam}.mafft_out.fa"
        loo_recs = drop_rows_and_empty_columns(all3_recs, drop_set)
        write_fasta(loo_fa, loo_recs)
        (CL_ALN_DIR / f"outgroup_taxa_loo_no_{drop_fam}.txt").write_text(
            "\n".join(sorted(kept_og)) + "\n", encoding="utf-8"
        )
        print(
            f"[CL0013 LOO] Drop {drop_fam} -> {len(loo_recs)} rows, {len(loo_recs[0][1])} cols, "
            f"remaining outgroup={len(kept_og)}"
        )
        tasks.append(
            {
                "wave": 1,
                "name": f"CL0013.loo_no_{drop_fam}.mafft_out",
                "alignment": str(loo_fa),
                "model_from": "",
            }
        )

    cl_tasks_path = OUT_DIR / "tree_tasks_cl0013.tsv"
    with open(cl_tasks_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["wave", "name", "alignment", "model_from"],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(tasks)
    print(f"[CL0013] Wrote {cl_tasks_path} ({len(tasks)} trees)")
    return tasks


def write_pmsf_tasks() -> None:
    specs = [
        ("PF00144.pmsf_c20.hmm_out", "PF00144", "primary.hmm_out.fa", "PF00144.primary.hmm_out.treefile"),
        ("PF00905.pmsf_c20.hmm_out", "PF00905", "primary.hmm_out.fa", "PF00905.primary.hmm_out.treefile"),
        ("PF13354.pmsf_c20.hmm_out", "PF13354", "primary.hmm_out.fa", "PF13354.primary.hmm_out.treefile"),
        ("PF00768.pmsf_c20.hmm_out", "PF00768", "primary.hmm_out.fa", "PF00768.primary.hmm_out.treefile"),
        ("PF13354_PF00768.pmsf_c20.mafft_out", "PF13354_PF00768", "primary.mafft_out.fa", "PF13354_PF00768.primary.mafft_out.treefile"),
    ]
    rows = []
    for name, tag, aln_fn, tree_fn in specs:
        aln = OUT_DIR / "alignments" / tag / aln_fn
        guide = OUT_DIR / "trees" / tag / tree_fn
        if not aln.is_file():
            raise FileNotFoundError(f"Missing PMSF alignment: {aln}")
        if not guide.is_file():
            raise FileNotFoundError(f"Missing PMSF guide tree: {guide}")
        rows.append({"name": name, "alignment": str(aln), "guide_tree": str(guide)})

    pmsf_path = OUT_DIR / "tree_tasks_pmsf.tsv"
    with open(pmsf_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["name", "alignment", "guide_tree"],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)
    print(f"[PMSF] Wrote {pmsf_path} ({len(rows)} trees)")


def main() -> int:
    v_meta, b_meta = load_catalogs()
    build_cl0013_alignments(v_meta, b_meta)
    write_pmsf_tasks()
    return 0


if __name__ == "__main__":
    sys.exit(main())
