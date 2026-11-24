#!/usr/bin/env python3
"""Structure-guided (FoldMason 3D + profile) CL0013 superfamily alignment & clean GTDB sister outgroups.

Two complementary rooting checks that do not require re-running the 50 sensitivity trees:
  Part A (structure-guided CL0013 superfamily, n = 448):
    1. Extract & HMM-slice Chain A crystal domains for 14 InterPro-verified PDB/mmCIF anchors
       across the 5 PRE families (PF00144, PF00905, PF13354, PF00768, PF02113) and the
       2 full-length non-PRE outgroup families (PF04960 [~301 aa], PF06089 [~336 aa]),
       dropping the short (~200 aa) PF27310 family that degraded PF00905 monophyly.
    2. Align the 14 sliced crystal structures in 3D via FoldMason (`foldmason easy-msa`).
    3. Sample 90 PRE reps per focal family strictly from `{fam}_combined_tree_dom.faa`
       (prioritizing `phylo_dom` diversity + all hydrolytic viral SBLs & Tier 1 loci ->
       100% tip overlap with `{fam}.primary.{track}_in.treefile`) + 40 PF02113 + 24 PF04960
       + 24 PF06089 (448 sequences total).
    4. Run per-family `hmmalign --outformat A2M` anchored on each family's primary crystal
       structure, project HMM match columns onto the 3D FoldMason master grid, and trim
       gaps with ClipKIT (-m gappy -g 0.5). Derive LOO alignments (`only_PF04960`,
       `only_PF06089`) and write `tree_tasks_cl0013_struct.tsv`.

  Part B (Clean Slow-Evolving GTDB Chromosomal Single-Family Outgroups):
    1. From MAD-rooted single-family ingroup trees (`PF00144` & `PF00905`), select 12
       taxonomically diverse, short-branch, chromosomal GTDB sequences (`hmm_cov >= 0.85`,
       zero SARG/plasmid/prophage sequences).
    2. Align each 12-sequence sister outgroup with its own family HMM and merge onto
       `alignments/{fam}/primary.{track}_in.fa` via `mafft --localpair --maxiterate 1000
       --addprofile` (`primary.clean_prof_out.fa`), preserving ingroup HMM columns.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import math
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from ete4 import Tree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

csv.field_size_limit(10_000_000)

IN_DIR = hpc.BASE / "01_final_curated"
INPUTS = hpc.BASE / "inputs"
PDB_RAW_DIR = INPUTS / "cl0013_ref_pdbs"
PDB_SLICED_DIR = PDB_RAW_DIR / "sliced"
OUT_DIR = hpc.PHASE3_DIR
ALN_DIR = OUT_DIR / "alignments"
CL_ALN_DIR = ALN_DIR / "CL0013"
TREES_DIR = OUT_DIR / "trees"

PRE_FOUR = ("PF00144", "PF00905", "PF13354", "PF00768")
OUTGROUP_TWO = ("PF04960", "PF06089")
ALL_SEVEN = PRE_FOUR + ("PF02113",) + OUTGROUP_TWO
THREADS = 16
CLIPKIT_GAP_THRESHOLD = 0.5

# First item in each list is the primary HMM->3D bridge anchor for that family
CRYSTAL_ANCHORS: dict[str, list[tuple[str, str]]] = {
    "PF00144": [("PF00144_1BLS", "A"), ("PF00144_1CEF", "A")],
    "PF00905": [("PF00905_1K25", "A"), ("PF00905_1E3U", "A")],
    "PF13354": [("PF13354_1BUE", "A"), ("PF13354_1ALQ", "A")],
    "PF00768": [("PF00768_1ES2", "A"), ("PF00768_1ES4", "A")],
    "PF02113": [("PF02113_1W79", "A"), ("PF02113_1W5D", "A")],
    "PF04960": [("PF04960_1U60", "A"), ("PF04960_1MKI", "A")],
    "PF06089": [("PF06089_7OS5", "A"), ("PF06089_7OS3", "A")],
}

SISTER_OUTGROUP_MAP: dict[str, str] = {
    "PF00144": "PF00905",
    "PF00905": "PF00144",
    "PF13354": "PF00905",
    "PF00768": "PF00905",
    "PF13354_PF00768": "PF00905",
}

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
    "VAL": "V", "SEC": "C", "PYL": "K",
}


def run(cmd: str, stdout_path: Path | None = None) -> str:
    if stdout_path is None:
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    else:
        with open(stdout_path, "w", encoding="utf-8", newline="\n") as fh:
            p = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd}\n{p.stderr}")
    return p.stdout


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
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def hmm_leng(hmm_path: Path) -> int:
    with open(hmm_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("LENG "):
                return int(line.split()[1])
    raise ValueError(f"No LENG in {hmm_path}")


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


def motif_lock_rate(recs: list[tuple[str, str]], pattern: str) -> float:
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


# ---------------------------------------------------------------------------
# Structure Parsing & HMM Domain Slicing (PDB + mmCIF)
# ---------------------------------------------------------------------------
def _tokenize_cif_line(line: str) -> list[str]:
    tokens: list[str] = []
    i, n = 0, len(line)
    while i < n:
        while i < n and line[i].isspace():
            i += 1
        if i >= n:
            break
        q = line[i]
        if q in ("'", '"'):
            i += 1
            start = i
            while i < n and not (line[i] == q and (i + 1 == n or line[i + 1].isspace())):
                i += 1
            tokens.append(line[start:i])
            if i < n and line[i] == q:
                i += 1
        else:
            start = i
            while i < n and not line[i].isspace():
                i += 1
            tokens.append(line[start:i])
    return tokens


def _parse_pdb_residues(pdb_path: Path, pref_chain: str) -> list[tuple[str, str, list[dict]]]:
    by_chain: dict[str, dict[tuple[str, str], dict]] = defaultdict(dict)
    order: dict[str, list[tuple[str, str]]] = defaultdict(list)

    with open(pdb_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            rec = line[:6].strip()
            resname = line[17:20].strip()
            if rec != "ATOM" and not (rec == "HETATM" and resname == "MSE"):
                continue
            if resname not in AA3TO1:
                continue
            alt = line[16] if len(line) > 16 else " "
            if alt not in (" ", "A", ""):
                continue
            chain = (line[21] if len(line) > 21 else "A").strip() or "A"
            resseq = line[22:26].strip()
            icode = line[26] if len(line) > 26 else " "
            key = (resseq, icode)
            if key not in by_chain[chain]:
                by_chain[chain][key] = {"aa": AA3TO1[resname], "resname": "MET" if resname == "MSE" else resname, "atoms": []}
                order[chain].append(key)
            aname = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occ = float(line[54:60]) if len(line) >= 60 and line[54:60].strip() else 1.0
            bfac = float(line[60:66]) if len(line) >= 66 and line[60:66].strip() else 20.0
            elem = (line[76:78].strip() if len(line) >= 78 else aname[0]) or aname[0]
            by_chain[chain][key]["atoms"].append({
                "aname": aname, "resname": "MET" if resname == "MSE" else resname,
                "resseq": int(re.sub(r"[^0-9\-]", "", resseq) or 1), "icode": icode if icode != " " else " ",
                "x": x, "y": y, "z": z, "occ": occ, "bfac": bfac, "elem": elem,
            })

    def chain_res(ch: str) -> list[tuple[str, str, list[dict]]]:
        out = []
        for k in order.get(ch, []):
            entry = by_chain[ch][k]
            if any(a["aname"] == "CA" for a in entry["atoms"]):
                out.append((entry["aa"], entry["resname"], entry["atoms"]))
        return out

    res = chain_res(pref_chain)
    if len(res) >= 100:
        return res
    best_ch = max(order.keys(), key=lambda c: len(chain_res(c)))
    return chain_res(best_ch)


def _parse_cif_residues(cif_path: Path, pref_chain: str) -> list[tuple[str, str, list[dict]]]:
    cols: list[str] = []
    by_chain: dict[str, dict[tuple[str, str], dict]] = defaultdict(dict)
    order: dict[str, list[tuple[str, str]]] = defaultdict(list)

    in_atom_loop = False
    reading_headers = False
    first_model: str | None = None

    with open(cif_path, encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith("#"):
                if in_atom_loop and not reading_headers and by_chain:
                    break
                continue
            if s == "loop_":
                if in_atom_loop and by_chain:
                    break
                in_atom_loop = False
                reading_headers = True
                cols = []
                continue
            if reading_headers and s.startswith("_atom_site."):
                in_atom_loop = True
                cols.append(s.split(".", 1)[1])
                continue
            elif reading_headers and s.startswith("_"):
                in_atom_loop = False
                reading_headers = False
                cols = []
                continue
            elif in_atom_loop:
                reading_headers = False
                if s.startswith("_") or s == "loop_":
                    break
                toks = _tokenize_cif_line(s)
                if len(toks) < len(cols):
                    continue
                row = dict(zip(cols, toks))
                mdl = row.get("pdbx_PDB_model_num", "1")
                if first_model is None:
                    first_model = mdl
                elif mdl != first_model:
                    break
                grp = row.get("group_PDB", "ATOM")
                resname = row.get("label_comp_id") or row.get("auth_comp_id", "")
                if grp != "ATOM" and not (grp == "HETATM" and resname == "MSE"):
                    continue
                if resname not in AA3TO1:
                    continue
                alt = row.get("label_alt_id", ".")
                if alt not in (".", "?", "", "A"):
                    continue
                chain = row.get("auth_asym_id") or row.get("label_asym_id", "A")
                seq_id = row.get("auth_seq_id") or row.get("label_seq_id", "1")
                ins = row.get("pdbx_PDB_ins_code", ".")
                icode = " " if ins in (".", "?", "") else ins[0]
                key = (seq_id, icode)
                if key not in by_chain[chain]:
                    by_chain[chain][key] = {"aa": AA3TO1[resname], "resname": "MET" if resname == "MSE" else resname, "atoms": []}
                    order[chain].append(key)
                aname = (row.get("label_atom_id") or row.get("auth_atom_id", "CA")).strip()
                resseq_int = int(re.sub(r"[^0-9\-]", "", seq_id) or 1)
                by_chain[chain][key]["atoms"].append({
                    "aname": aname,
                    "resname": "MET" if resname == "MSE" else resname,
                    "resseq": resseq_int,
                    "icode": icode,
                    "x": float(row["Cartn_x"]),
                    "y": float(row["Cartn_y"]),
                    "z": float(row["Cartn_z"]),
                    "occ": float(row["occupancy"]) if row.get("occupancy", ".") not in (".", "?") else 1.0,
                    "bfac": float(row["B_iso_or_equiv"]) if row.get("B_iso_or_equiv", ".") not in (".", "?") else 20.0,
                    "elem": (row.get("type_symbol", "C") or aname[0]).strip(),
                })

    def chain_res(ch: str) -> list[tuple[str, str, list[dict]]]:
        out = []
        for k in order.get(ch, []):
            entry = by_chain[ch][k]
            if any(a["aname"] == "CA" for a in entry["atoms"]):
                out.append((entry["aa"], entry["resname"], entry["atoms"]))
        return out

    res = chain_res(pref_chain)
    if len(res) >= 100:
        return res
    best_ch = max(order.keys(), key=lambda c: len(chain_res(c)))
    return chain_res(best_ch)


def extract_and_slice_crystal(fam: str, stem: str, chain: str) -> tuple[Path, str]:
    pdb_path = PDB_RAW_DIR / f"{stem}.pdb"
    cif_path = PDB_RAW_DIR / f"{stem}.cif"
    if pdb_path.is_file():
        residues = _parse_pdb_residues(pdb_path, chain)
    elif cif_path.is_file():
        residues = _parse_cif_residues(cif_path, chain)
    else:
        raise FileNotFoundError(f"Missing structure for {stem} in {PDB_RAW_DIR}")

    full_seq = "".join(aa for aa, _, _ in residues)
    if len(full_seq) < 120:
        raise RuntimeError(f"{stem}: extracted chain {chain} too short ({len(full_seq)} aa)")

    PDB_SLICED_DIR.mkdir(parents=True, exist_ok=True)
    tmp_fa = PDB_SLICED_DIR / f"_tmp_{stem}.fa"
    tmp_dom = PDB_SLICED_DIR / f"_tmp_{stem}.domtbl"
    write_fasta(tmp_fa, [(stem, full_seq)])
    hmm = INPUTS / f"{fam}_official.hmm"
    run(f"{hpc.HMMSEARCH_BIN} --domtblout {tmp_dom} --noali {hmm} {tmp_fa}")

    best_score = -1.0
    env_from, env_to = 1, len(full_seq)
    with open(tmp_dom, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            score = float(p[13])
            if score > best_score:
                best_score = score
                env_from = int(p[19])
                env_to = int(p[20])
    tmp_fa.unlink(missing_ok=True)
    tmp_dom.unlink(missing_ok=True)
    if best_score < 0:
        raise RuntimeError(f"{stem}: no {fam} HMM hit in extracted chain sequence")

    lo = max(0, env_from - 1 - 5)
    hi = min(len(residues), env_to + 5)
    sliced_res = residues[lo:hi]
    sliced_seq = "".join(aa for aa, _, _ in sliced_res)

    out_pdb = PDB_SLICED_DIR / f"{stem}_A_dom.pdb"
    serial = 1
    with open(out_pdb, "w", encoding="utf-8", newline="\n") as fh:
        for res_idx, (_, _, atoms) in enumerate(sliced_res, start=1):
            for a in atoms:
                aname = a["aname"]
                if len(aname) < 4 and not aname[0].isdigit() and len(a["elem"]) == 1:
                    aname_fmt = f" {aname:<3s}"
                else:
                    aname_fmt = f"{aname:<4s}"
                fh.write(
                    f"ATOM  {serial:5d} {aname_fmt:4s} {a['resname']:>3s} A{res_idx:4d}    "
                    f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{a['occ']:6.2f}{a['bfac']:6.2f}"
                    f"          {a['elem']:>2s}  \n"
                )
                serial += 1
        fh.write("END\n")

    print(f"  [Crystal Slice] {stem} -> {out_pdb.name} ({len(sliced_seq)} aa, hmm_score={best_score:.1f})")
    return out_pdb, sliced_seq


# ---------------------------------------------------------------------------
# Part A: FoldMason 3D + Profile-Anchored CL0013 Superfamily (no_PF27310)
# ---------------------------------------------------------------------------
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


def run_hmmalign_raw(fam: str, recs: list[tuple[str, str]], work_dir: Path) -> tuple[int, dict[str, str]]:
    in_fa = work_dir / f"_{fam}_in.fa"
    a2m_path = work_dir / f"_{fam}.a2m"
    hmm = INPUTS / f"{fam}_official.hmm"
    leng = hmm_leng(hmm)
    write_fasta(in_fa, recs)
    run(f"{hpc.HMMALIGN_BIN} --outformat A2M {hmm} {in_fa}", stdout_path=a2m_path)
    out = dict(read_fasta(a2m_path))
    in_fa.unlink(missing_ok=True)
    a2m_path.unlink(missing_ok=True)
    return leng, out


def map_anchor_seq_to_fm_cols(sliced_seq: str, fm_row: str) -> dict[int, int]:
    fm_pos = [i for i, c in enumerate(fm_row) if c not in "-."]
    fm_ungap = "".join(fm_row[i] for i in fm_pos).upper()
    if fm_ungap == sliced_seq:
        return {r: fm_pos[r] for r in range(len(sliced_seq))}
    sm = difflib.SequenceMatcher(None, sliced_seq, fm_ungap, autojunk=False)
    mapping: dict[int, int] = {}
    for a0, b0, size in sm.get_matching_blocks():
        for k in range(size):
            mapping[a0 + k] = fm_pos[b0 + k]
    if len(mapping) < 0.85 * len(sliced_seq):
        raise RuntimeError(f"Low FoldMason-to-crystal identity overlap: {len(mapping)}/{len(sliced_seq)}")
    return mapping


def build_cl0013_struct_alignments(v_meta: dict[str, dict], b_meta: dict[str, dict]) -> None:
    CL_ALN_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Slice 14 crystal structures to Chain A HMM domain cores & align with FoldMason
    sliced_pdbs: list[Path] = []
    anchor_seqs: dict[str, str] = {}
    print("[Part A.1] Extracting & HMM-slicing 14 crystal anchor structures...")
    for fam in ALL_SEVEN:
        for stem, chain in CRYSTAL_ANCHORS[fam]:
            spdb, sseq = extract_and_slice_crystal(fam, stem, chain)
            sliced_pdbs.append(spdb)
            anchor_seqs[stem] = sseq

    fm_prefix = CL_ALN_DIR / "cl0013_foldmason_seed"
    fm_tmp = CL_ALN_DIR / "tmp_foldmason"
    if fm_tmp.exists():
        shutil.rmtree(fm_tmp)
    pdb_args = " ".join(str(p) for p in sliced_pdbs)
    run(
        f"{hpc.FOLDMASON_BIN} easy-msa {pdb_args} {fm_prefix} {fm_tmp} "
        f"--threads {THREADS} --filter-msa 0 --recompute-scores 1 --refine-iters 2"
    )
    fm_aa_fa = CL_ALN_DIR / "cl0013_foldmason_seed_aa.fa"
    if not fm_aa_fa.is_file():
        cands = list(CL_ALN_DIR.glob("cl0013_foldmason_seed*aa*.fa*"))
        if not cands:
            raise FileNotFoundError(f"FoldMason output FASTA not found at {fm_prefix}*")
        fm_aa_fa = cands[0]

    fm_recs = {}
    for hdr, seq in read_fasta(fm_aa_fa):
        key = re.sub(r"_A_dom(\.pdb)?$", "", hdr)
        fm_recs[key] = seq
    if len(fm_recs) != len(sliced_pdbs):
        raise RuntimeError(f"Expected {len(sliced_pdbs)} anchors in FoldMason MSA, got {len(fm_recs)}: {sorted(fm_recs)}")

    # 2. Sample 448 CL0013 sequences (400 PRE with 100% primary-tree overlap + 48 full-length non-PRE)
    selected_seqs: dict[str, dict[str, str]] = {}  # uid -> seq + meta
    seqs_by_fam: dict[str, list[tuple[str, str]]] = defaultdict(list)

    def register(uid: str, seq: str, fam: str, origin: str, role: str, f_tier: str) -> None:
        if uid not in selected_seqs:
            selected_seqs[uid] = {
                "uid": uid, "seq": seq, "cl0013_family": fam,
                "origin": origin, "role": role, "functional_tier": f_tier,
            }
            seqs_by_fam[fam].append((uid, seq))

    for fam in PRE_FOUR:
        prim_recs = read_fasta(IN_DIR / f"{fam}_combined_tree_dom.faa")
        phylo_uids = {u for u, _ in read_fasta(IN_DIR / f"{fam}_combined_tree_phylo_dom.faa")}
        v_nr85_recs = read_fasta(IN_DIR / f"{fam}_viral_nr85_dom.faa")
        v_nr85_set = {u for u, _ in v_nr85_recs}

        # Bacterial 60: strictly from primary tree (so 100% overlap for crown projection),
        # preferring sequences also in the phylo 65%-diversity set
        prim_bac_phylo = [
            (u, s) for u, s in prim_recs
            if u in b_meta and u not in v_meta and u not in v_nr85_set and u in phylo_uids and u not in selected_seqs
        ]
        prim_bac_other = [
            (u, s) for u, s in prim_recs
            if u in b_meta and u not in v_meta and u not in v_nr85_set and u not in phylo_uids and u not in selected_seqs
        ]
        picked_bac = evenly_spaced(prim_bac_phylo, 60)
        if len(picked_bac) < 60:
            picked_bac += evenly_spaced(prim_bac_other, 60 - len(picked_bac))
        for uid, seq in picked_bac:
            info = b_meta.get(uid, {})
            register(uid, seq, fam, "Bacterial", "PRE_ingroup", info.get("functional_tier", "PBP_or_PRE_Homolog"))

        # Viral 30: all verified hydrolytic SBLs first, then Tier1, then evenly spaced remainder
        v_pri = []
        v_other = []
        for uid, seq in v_nr85_recs:
            if uid in selected_seqs:
                continue
            info = v_meta.get(uid, {})
            is_sbl = info.get("functional_tier") == "Putative_Hydrolytic_BetaLactamase"
            is_t1 = info.get("viral_origin_tier") == "Tier1_HighConfidence_Viral"
            if is_sbl or is_t1:
                v_pri.append((uid, seq, 0 if is_sbl else 1, 0 if is_t1 else 1))
            else:
                v_other.append((uid, seq))
        v_pri.sort(key=lambda x: (x[2], x[3], x[0]))
        sbl_only = [(u, s) for u, s, r0, _ in v_pri if r0 == 0]
        t1_only = [(u, s) for u, s, r0, _ in v_pri if r0 == 1]
        v_picked = sbl_only[:30]
        if len(v_picked) < 30:
            v_picked += evenly_spaced(t1_only, 30 - len(v_picked))
        if len(v_picked) < 30:
            v_picked += evenly_spaced(v_other, 30 - len(v_picked))
        for uid, seq in v_picked:
            info = v_meta.get(uid, {})
            register(uid, seq, fam, "Viral", "PRE_ingroup", info.get("functional_tier", "PBP_or_PRE_Homolog"))

    # Auxiliary PF02113 (40 reps)
    for uid, seq in read_fasta(IN_DIR / "CL0013_aux_PF02113_bacterial_rep_dom.faa"):
        if len(seqs_by_fam["PF02113"]) >= 40:
            break
        info = b_meta.get(uid, {})
        register(uid, seq, "PF02113", "Bacterial", "PRE_ingroup", info.get("functional_tier", "PBP_or_PRE_Homolog"))

    # Full-length non-PRE biological outgroups (PF04960 + PF06089; PF27310 excluded)
    outgroup_by_fam: dict[str, list[str]] = {}
    for og_fam in OUTGROUP_TWO:
        og_recs = read_fasta(IN_DIR / f"CL0013_outgroup_{og_fam}_dom.faa")
        outgroup_by_fam[og_fam] = []
        for uid, seq in og_recs:
            info = b_meta.get(uid, {})
            register(uid, seq, og_fam, "Bacterial", "NonPRE_outgroup", info.get("functional_tier", "CL0013_NonPRE_Outgroup"))
            outgroup_by_fam[og_fam].append(uid)

    fam_cnts = {f: len(seqs_by_fam[f]) for f in ALL_SEVEN}
    print(f"[Part A.2] Selected {len(selected_seqs)} CL0013 taxa (no PF27310): {fam_cnts}")

    # Preserve baseline cl0013_taxa_metadata.tsv entries (like PF27310) while adding/updating all 448 taxa
    meta_path = CL_ALN_DIR / "cl0013_taxa_metadata.tsv"
    merged_meta: dict[str, dict[str, str]] = {}
    if meta_path.is_file():
        with open(meta_path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                merged_meta[r["uid"]] = r
    for uid, item in selected_seqs.items():
        merged_meta[uid] = {k: item[k] for k in ("uid", "cl0013_family", "origin", "role", "functional_tier")}

    with open(meta_path, "w", newline="\n", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["uid", "cl0013_family", "origin", "role", "functional_tier"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in merged_meta.values():
            w.writerow(r)

    struct_meta_path = CL_ALN_DIR / "cl0013_struct_taxa_metadata.tsv"
    with open(struct_meta_path, "w", newline="\n", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["uid", "cl0013_family", "origin", "role", "functional_tier"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        for uid, item in selected_seqs.items():
            w.writerow({k: item[k] for k in ("uid", "cl0013_family", "origin", "role", "functional_tier")})

    # 3. Per-family hmmalign projected onto FoldMason 3D master coordinate grid
    all_slots: set[tuple[int, int, str, int]] = set()
    # fam -> {uid: {slot_tuple: aa_char}}
    fam_slot_rows: dict[str, dict[str, dict[tuple[int, int, str, int], str]]] = {}

    for fam in ALL_SEVEN:
        primary_stem = CRYSTAL_ANCHORS[fam][0][0]
        anchor_seq = anchor_seqs[primary_stem]
        res_to_fm = map_anchor_seq_to_fm_cols(anchor_seq, fm_recs[primary_stem])
        f_recs = [(primary_stem, anchor_seq)] + seqs_by_fam[fam]
        leng, a2m_map = run_hmmalign_raw(fam, f_recs, CL_ALN_DIR)

        # Build monotonic mapping HMM match column m (0..LENG-1) -> master slot tuple
        anchor_a2m = a2m_map[primary_stem]
        col_to_slot: list[tuple[int, int, str, int]] = []
        r_idx = 0
        m_idx = 0
        last_fm = -1
        ins_offset = 0
        for ch in anchor_a2m:
            if ch in (".",):
                continue
            if ch == "-":
                # Match deletion in anchor
                ins_offset += 1
                col_to_slot.append((last_fm, 1, fam, ins_offset))
                m_idx += 1
            elif ch.isupper():
                fm_c = res_to_fm.get(r_idx)
                if fm_c is not None and fm_c > last_fm:
                    col_to_slot.append((fm_c, 0, "", 0))
                    last_fm = fm_c
                    ins_offset = 0
                else:
                    ins_offset += 1
                    col_to_slot.append((last_fm, 1, fam, ins_offset))
                r_idx += 1
                m_idx += 1
            else:
                # lowercase insertion state
                r_idx += 1

        if m_idx != leng or len(col_to_slot) != leng:
            raise RuntimeError(f"{fam}: match state count {m_idx}/{len(col_to_slot)} != LENG {leng}")

        for sl in col_to_slot:
            all_slots.add(sl)

        fam_rows: dict[str, dict[tuple[int, int, str, int], str]] = {}
        for uid, _ in seqs_by_fam[fam]:
            m_str = "".join(c for c in a2m_map[uid] if c.isupper() or c == "-")
            if len(m_str) != leng:
                raise RuntimeError(f"{fam}/{uid}: match len {len(m_str)} != {leng}")
            fam_rows[uid] = {col_to_slot[i]: m_str[i] for i in range(leng)}
        fam_slot_rows[fam] = fam_rows

    sorted_slots = sorted(all_slots)
    raw_master_fa = CL_ALN_DIR / "struct_c20.no_PF27310.raw.fa"
    out_master_fa = CL_ALN_DIR / "struct_c20.no_PF27310.fa"
    assembled: list[tuple[str, str]] = []
    for fam in ALL_SEVEN:
        for uid, _ in seqs_by_fam[fam]:
            smap = fam_slot_rows[fam][uid]
            assembled.append((uid, "".join(smap.get(sl, "-") for sl in sorted_slots)))

    write_fasta(raw_master_fa, assembled)
    run(f"{hpc.CLIPKIT_BIN} {raw_master_fa} -m gappy -g {CLIPKIT_GAP_THRESHOLD} -o {out_master_fa}")
    trimmed = read_fasta(out_master_fa)
    pre_trimmed = [(u, s) for u, s in trimmed if selected_seqs[u]["cl0013_family"] not in OUTGROUP_TWO]
    n_cols = len(trimmed[0][1])
    occ = sum(sum(c not in "-." for c in s) for _, s in trimmed) / max(len(trimmed) * n_cols, 1)
    sxxk_rate = motif_lock_rate(pre_trimmed, r"S..K")
    print(
        f"[Part A.3] CL0013 FoldMason+hmmalign MSA: raw_cols={len(sorted_slots)} -> ClipKIT_cols={n_cols}, "
        f"occ={occ:.4f}, PRE_SXXK_lock={sxxk_rate:.4f}"
    )

    # Derive outgroup files & LOO alignments
    og_no27 = sorted(u for f in OUTGROUP_TWO for u in outgroup_by_fam[f])
    (CL_ALN_DIR / "outgroup_struct_no_PF27310.txt").write_text("\n".join(og_no27) + "\n", encoding="utf-8")

    tasks = [
        {"wave": 1, "name": "CL0013.struct_c20.no_PF27310", "alignment": str(out_master_fa), "model_from": "LG+C20+F+G"},
    ]
    for keep_og in OUTGROUP_TWO:
        drop_og = [f for f in OUTGROUP_TWO if f != keep_og][0]
        drop_set = set(outgroup_by_fam[drop_og])
        kept_og = sorted(outgroup_by_fam[keep_og])
        loo_fa = CL_ALN_DIR / f"struct_c20.only_{keep_og}.fa"
        loo_recs = drop_rows_and_empty_columns(trimmed, drop_set)
        write_fasta(loo_fa, loo_recs)
        (CL_ALN_DIR / f"outgroup_struct_only_{keep_og}.txt").write_text("\n".join(kept_og) + "\n", encoding="utf-8")
        tasks.append({
            "wave": 1, "name": f"CL0013.struct_c20.only_{keep_og}",
            "alignment": str(loo_fa), "model_from": "LG+C20+F+G",
        })
        print(f"  [LOO only_{keep_og}] rows={len(loo_recs)}, cols={len(loo_recs[0][1])}, og={len(kept_og)}")

    tasks_path = OUT_DIR / "tree_tasks_cl0013_struct.tsv"
    with open(tasks_path, "w", newline="\n", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["wave", "name", "alignment", "model_from"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(tasks)
    print(f"[Part A Done] Wrote {tasks_path} ({len(tasks)} trees)")


# ---------------------------------------------------------------------------
# Part B: Clean Chromosomal Slow-Evolving GTDB Single-Family Profile Outgroup
# ---------------------------------------------------------------------------
def _root_to_tip_distances_from_mad(tree_path: Path) -> dict[str, float]:
    """Quickly approximate evolutionary rate via FastRoot/MinVar centroid depth (O(N))."""
    t = Tree(tree_path.read_text(encoding="utf-8"), parser=1)
    n = len(list(t.leaves()))
    for nd in t.traverse("postorder"):
        if nd.is_leaf:
            nd.add_prop("_n_in", 1)
            nd.add_prop("_sum_in", 0.0)
            nd.add_prop("_sq_in", 0.0)
        else:
            n_in = 0
            s_in = 0.0
            q_in = 0.0
            for c in nd.children:
                L = max(0.0, float(c.dist or 0.0))
                cn, cs, cq = c.props["_n_in"], c.props["_sum_in"], c.props["_sq_in"]
                n_in += cn
                s_in += cs + cn * L
                q_in += cq + 2.0 * L * cs + cn * L * L
            nd.add_prop("_n_in", n_in)
            nd.add_prop("_sum_in", s_in)
            nd.add_prop("_sq_in", q_in)

    t.add_prop("_sum_out", 0.0)
    t.add_prop("_sq_out", 0.0)
    best_var = float("inf")
    best_nd = None
    best_alpha = 0.5
    for p in t.traverse("preorder"):
        if p.is_leaf:
            continue
        ps = p.props["_sum_out"] + p.props["_sum_in"]
        pq = p.props["_sq_out"] + p.props["_sq_in"]
        for c in p.children:
            L = max(0.0, float(c.dist or 0.0))
            na = c.props["_n_in"]
            nb = n - na
            if na == 0 or nb == 0:
                continue
            sa, qa = c.props["_sum_in"], c.props["_sq_in"]
            spo = ps - (sa + na * L)
            qpo = pq - (qa + 2.0 * L * sa + na * L * L)
            c.add_prop("_sum_out", spo + nb * L)
            c.add_prop("_sq_out", qpo + 2.0 * L * spo + nb * L * L)
            c0 = sa + spo + nb * L
            c1 = (na - nb) * L
            d0 = qa + qpo + 2.0 * L * spo + nb * L * L
            d1 = 2.0 * L * (sa - spo - nb * L)
            d2 = n * L * L
            aq = (d2 / n) - (c1 / n) ** 2
            bq = (d1 / n) - 2.0 * (c0 * c1) / (n * n)
            alpha = min(1.0, max(0.0, -bq / (2.0 * aq))) if aq > 1e-18 else 0.5
            v = max(0.0, (d0 + d1 * alpha + d2 * alpha * alpha) / n - ((c0 + c1 * alpha) / n) ** 2)
            if v < best_var:
                best_var = v
                best_nd = c
                best_alpha = alpha

    if best_nd is not None:
        orig_L = max(0.0, float(best_nd.dist or 0.0))
        t.set_outgroup(best_nd)
        for c in t.children:
            c.dist = best_alpha * orig_L if c is best_nd else (1.0 - best_alpha) * orig_L

    dists: dict[str, float] = {id(t): 0.0}
    out: dict[str, float] = {}
    for nd in t.traverse("preorder"):
        if nd.up is not None:
            dists[id(nd)] = dists[id(nd.up)] + max(0.0, float(nd.dist or 0.0))
        if nd.is_leaf:
            out[nd.name] = dists[id(nd)]
    return out


def build_clean_gtdb_profile_outgroups(b_meta: dict[str, dict]) -> None:
    print("[Part B] Building Clean Slow-Evolving Chromosomal GTDB Single-Family Outgroups via mafft --addprofile...")
    # Cache root-to-tip distances and domain sequences for PF00905 & PF00144 sister donors
    sis_rates: dict[str, dict[str, float]] = {}
    sis_seqs: dict[str, dict[str, str]] = {}
    for sis in ("PF00144", "PF00905"):
        sis_in_tree = TREES_DIR / sis / f"{sis}.primary.hmm_in.treefile"
        sis_rates[sis] = _root_to_tip_distances_from_mad(sis_in_tree)
        sis_seqs[sis] = dict(read_fasta(IN_DIR / f"{sis}_combined_tree_dom.faa"))

    tasks: list[dict[str, str]] = []
    for fam, sis in SISTER_OUTGROUP_MAP.items():
        track = "mafft" if fam == "PF13354_PF00768" else "hmm"
        fdir = ALN_DIR / fam
        in_fa = fdir / f"primary.{track}_in.fa"
        if not in_fa.is_file():
            raise FileNotFoundError(f"Missing ingroup alignment: {in_fa}")
        in_uids = {u for u, _ in read_fasta(in_fa)}

        # Rank clean chromosomal GTDB reps by ascending branch length
        cands: list[tuple[float, str, str]] = []
        for uid, rdist in sis_rates[sis].items():
            if not uid.startswith("GTDB|") or uid in in_uids or uid not in sis_seqs[sis]:
                continue
            info = b_meta.get(uid, {})
            if info.get("replicon", "") != "chromosome":
                continue
            if str(info.get("in_bacterial_prophage", "False")).lower() in ("true", "1", "yes"):
                continue
            if str(info.get("is_viral_contig_in_genome", "False")).lower() in ("true", "1", "yes"):
                continue
            if float(info.get("hmm_cov") or 0.0) < 0.85:
                continue
            cands.append((rdist, uid, info.get("family", "")))
        cands.sort(key=lambda x: (x[0], x[1]))

        picked: list[tuple[str, str]] = []
        seen_tax_fams: set[str] = set()
        for rdist, uid, tfam in cands:
            if tfam and tfam in seen_tax_fams:
                continue
            picked.append((uid, sis_seqs[sis][uid]))
            if tfam:
                seen_tax_fams.add(tfam)
            if len(picked) >= 12:
                break
        if len(picked) < 12:
            used = {u for u, _ in picked}
            for rdist, uid, _ in cands:
                if uid not in used:
                    picked.append((uid, sis_seqs[sis][uid]))
                    used.add(uid)
                if len(picked) >= 12:
                    break

        if len(picked) < 8:
            raise RuntimeError(f"{fam}: only found {len(picked)} clean GTDB sister candidates from {sis}")

        # Align the 12 clean outgroup sequences with the sister family's official HMM
        og_raw_fa = fdir / "_clean_og12.unaligned.fa"
        og_a2m = fdir / "_clean_og12.a2m"
        og_hmm_fa = fdir / "clean_og12.hmm.fa"
        write_fasta(og_raw_fa, picked)
        hmm = INPUTS / f"{sis}_official.hmm"
        run(f"{hpc.HMMALIGN_BIN} --outformat A2M {hmm} {og_raw_fa}", stdout_path=og_a2m)
        og_recs = [(n, "".join(c for c in s if c.isupper() or c == "-")) for n, s in read_fasta(og_a2m)]
        write_fasta(og_hmm_fa, og_recs)
        og_raw_fa.unlink(missing_ok=True)
        og_a2m.unlink(missing_ok=True)

        prof_out_fa = fdir / f"primary.clean_prof_out.{track}_out.fa"
        run(
            f"{hpc.MAFFT_BIN} --localpair --maxiterate 1000 --thread {THREADS} "
            f"--addprofile {og_hmm_fa} {in_fa}",
            stdout_path=prof_out_fa,
        )
        og_txt = fdir / "outgroup_clean_gtdb12.txt"
        og_txt.write_text("\n".join(sorted(u for u, _ in picked)) + "\n", encoding="utf-8")

        tasks.append({
            "wave": "2",
            "name": f"{fam}.clean_prof_out.{track}_out",
            "alignment": str(prof_out_fa),
            "model_from": f"{fam}.primary.{track}_out",
        })
        med_dist = cands[min(len(cands) - 1, 5)][0]
        print(
            f"  [{fam:15s} <- {sis}] {len(picked)} clean chromosomal GTDB reps "
            f"({len(seen_tax_fams)} tax families, median_root_dist={med_dist:.3f}) -> {prof_out_fa.name}"
        )

    tp = OUT_DIR / "tree_tasks_clean_prof.tsv"
    with open(tp, "w", newline="\n", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["wave", "name", "alignment", "model_from"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(tasks)
    print(f"[Part B Done] Wrote {tp} ({len(tasks)} trees)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", choices=("all", "partA", "partB"), default="all")
    args = ap.parse_args()
    v_meta, b_meta = load_catalogs()
    if args.step in ("all", "partA"):
        build_cl0013_struct_alignments(v_meta, b_meta)
    if args.step in ("all", "partB"):
        build_clean_gtdb_profile_outgroups(b_meta)
    return 0


if __name__ == "__main__":
    sys.exit(main())
