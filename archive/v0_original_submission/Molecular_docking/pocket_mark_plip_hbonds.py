#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import glob
import json
import re
from pathlib import Path
import pandas as pd
import xml.etree.ElementTree as ET

AA3_TO_AA1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "SEC":"U","PYL":"O","ASX":"B","GLX":"Z","XLE":"J","UNK":"X"
}

def aa3_to_aa1(x: str) -> str:
    if x is None:
        return "X"
    x = str(x).strip().upper()
    return AA3_TO_AA1.get(x, "X")

def parse_pockets(pocket_str: str):
    parts = []
    for chunk in pocket_str.replace(" ", ",").split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts.append(int(chunk))
    if not parts:
        raise argparse.ArgumentTypeError("No valid pocket given, e.g. -p 1,2")
    return sorted(set(parts))

def iter_files(items):
    files = []
    for it in items:
        it = os.path.expanduser(it)
        m = glob.glob(it)
        if m:
            files.extend(m)
        elif os.path.exists(it):
            files.append(it)
    # de-duplicate, keep order
    seen, uniq = set(), []
    for f in files:
        if f not in seen:
            uniq.append(f); seen.add(f)
    return uniq

def normalize_pocket_csv(df: pd.DataFrame) -> pd.DataFrame:
    # CSV column names may have leading spaces: strip, then lower-case
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={c: c.strip().lower() for c in df.columns})

    # required columns of the pocket-prediction CSV
    need = {
        "chain": ["chain"],
        "residue_label": ["residue_label", "residue_number", "resid", "resnum", "resi"],
        "residue_name": ["residue_name", "resname", "aa3"],
        "pocket": ["pocket"],
    }
    found = {}
    for std, cands in need.items():
        for c in cands:
            if c in df.columns:
                found[std] = c
                break
        if std not in found:
            raise ValueError(f"Missing required column {std}; columns = {list(df.columns)}")

    opt = {
        "score": ["score"],
        "zscore": ["zscore", "z_score"],
        "probability": ["probability", "prob"],
    }
    for std, cands in opt.items():
        for c in cands:
            if c in df.columns:
                found[std] = c
                break

    keep = list(found.values())
    df2 = df[keep].copy()
    inv = {v: k for k, v in found.items()}
    df2 = df2.rename(columns=inv)

    # type conversion
    df2["pocket"] = pd.to_numeric(df2["pocket"], errors="coerce")
    df2["residue_label"] = pd.to_numeric(df2["residue_label"], errors="coerce")
    df2 = df2.dropna(subset=["pocket", "residue_label", "chain", "residue_name"])
    df2["pocket"] = df2["pocket"].astype(int)
    df2["residue_label"] = df2["residue_label"].astype(int)
    df2["chain"] = df2["chain"].astype(str).str.strip()
    df2["residue_name"] = df2["residue_name"].astype(str).str.strip().str.upper()
    return df2

# ------------------------ PLIP parsing ------------------------

def _collect_from_json(obj, hits):
    """
    Recursively search a JSON report for hydrogen_bonds and collect protein residues as (chain, resnr).
    Accepts the common key names resnr/reschain or residue_number/chain.
    """
    if isinstance(obj, dict):
        for k, v in obj.items():
            lk = str(k).lower()
            if "hydrogen" in lk and "bond" in lk:
                # this level may itself be a list or dict
                _collect_from_json(v, hits)
            else:
                _collect_from_json(v, hits)
    elif isinstance(obj, list):
        for it in obj:
            _collect_from_json(it, hits)
    else:
        return

def parse_plip_json(path: Path):
    data = json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    hits = []
    # first collect every candidate interaction dict
    def walk(o):
        if isinstance(o, dict):
            # if this dict looks like an hbond record, read its residue
            keys = {str(x).lower() for x in o.keys()}
            if ("resnr" in keys and "reschain" in keys) or ("residue_number" in keys and "chain" in keys):
                chain = o.get("reschain", o.get("chain", None))
                resnr = o.get("resnr", o.get("residue_number", None))
                if chain is not None and resnr is not None:
                    try:
                        hits.append((str(chain).strip(), int(resnr)))
                    except Exception:
                        pass
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for it in o:
                walk(it)
    walk(data)
    return hits

def parse_plip_xml(path: Path):
    tree = ET.parse(path)
    root = tree.getroot()
    hits = []

    # find every node whose tag/name contains hydrogen_bond
    for elem in root.iter():
        tag = elem.tag.lower()
        if "hydrogen" in tag and "bond" in tag:
            # elem may be a <hydrogen_bond> record or a container
            # read reschain/resnr from child nodes
            sub = {child.tag.lower(): (child.text.strip() if child.text else "") for child in list(elem)}
            # common PLIP XML field names
            chain = sub.get("reschain") or sub.get("chain") or sub.get("protchain") or sub.get("protein_chain")
            resnr = sub.get("resnr") or sub.get("resnum") or sub.get("residue_number") or sub.get("protein_resnr")
            if chain and resnr:
                try:
                    hits.append((str(chain).strip(), int(resnr)))
                except Exception:
                    pass

    return hits

HB_TEXT_PATTERNS = [
    # pick fragments such as "C:ARG2325" or "Chain C Res 2325 ARG" from text
    re.compile(r"\b([A-Za-z0-9])\s*[:]\s*([A-Z]{3})\s*(\d+)\b"),
    re.compile(r"\bchain\s*([A-Za-z0-9])\b.*?\b(\d+)\b.*?\b([A-Z]{3})\b", re.IGNORECASE),
    re.compile(r"\bres(?:idue)?\s*(\d+)\b.*?\bchain\s*([A-Za-z0-9])\b", re.IGNORECASE),
]

def parse_plip_txt(path: Path):
    txt = path.read_text(encoding="utf-8", errors="ignore")
    hits = []

    # restrict to the hydrogen-bond section where possible
    # scan the whole file if no section header is found
    chunk = txt
    m = re.search(r"hydrogen\s*bonds?.*?\n(-+|\=+)", txt, flags=re.IGNORECASE)
    if m:
        chunk = txt[m.start():]

    for pat in HB_TEXT_PATTERNS:
        for mm in pat.finditer(chunk):
            groups = mm.groups()
            # group order differs between patterns
            if len(groups) == 3 and groups[0] and groups[2]:
                # pattern1: chain, aa3, resnr
                chain = groups[0]
                resnr = groups[2]
            elif len(groups) == 3:
                # pattern2: chain, resnr, aa3
                chain = groups[0]
                resnr = groups[1]
            elif len(groups) == 2:
                # pattern3: resnr, chain
                resnr = groups[0]
                chain = groups[1]
            else:
                continue
            try:
                hits.append((str(chain).strip(), int(resnr)))
            except Exception:
                pass

    return hits

def parse_plip_file(path: Path):
    suf = path.suffix.lower()
    if suf in [".json"]:
        return parse_plip_json(path)
    if suf in [".xml"]:
        return parse_plip_xml(path)
    # fallback: text
    return parse_plip_txt(path)

# ------------------------ matching / output ------------------------

def find_plip_for_csv(csv_path: Path, plip_dir: Path, exts=(".xml", ".json", ".txt")):
    """
    Find the PLIP file in plip_dir matching the CSV base name.
    Rule: stem prefix match (starting with csv_stem), preferring xml > json > txt.
    """
    base = csv_path.stem
    cand = []
    for ext in exts:
        cand.extend(sorted(plip_dir.glob(f"{base}*{ext}")))
    # prefer XML/JSON
    cand = sorted(cand, key=lambda p: {".xml":0, ".json":1, ".txt":2}.get(p.suffix.lower(), 9))
    return cand[0] if cand else None

def write_outputs(df, out_tsv: Path, out_fasta: Path, star=False):
    df = df.copy()
    df["aa1"] = df["residue_name"].apply(aa3_to_aa1)
    df = df.sort_values(["chain", "residue_label"])

    cols = ["chain", "residue_label", "residue_name", "aa1", "pocket", "hb"]
    for c in ["score", "zscore", "probability", "hb_count", "hb_sources"]:
        if c in df.columns:
            cols.append(c)
    df[cols].to_csv(out_tsv, sep="\t", index=False)

    if out_fasta is not None:
        lines = []
        for chain, g in df.groupby("chain", sort=False):
            # pocket sequence: hydrogen-bond sites marked with * or lower case
            seq = []
            for _, r in g.iterrows():
                aa = r["aa1"]
                if bool(r["hb"]):
                    aa = (aa.lower() if not star else aa + "*")
                seq.append(aa)
            seq = "".join(seq)
            resi = ",".join(str(x) for x in g["residue_label"].tolist())
            header = f">{out_tsv.stem}|chain={chain}|n={len(g)}|residues={resi}"
            lines.append(header)
            for i in range(0, len(seq), 60):
                lines.append(seq[i:i+60])
        out_fasta.write_text("\n".join(lines) + "\n", encoding="utf-8")

def main():
    ap = argparse.ArgumentParser(
        description="Extract residues of selected pockets from pocket-prediction CSVs and mark PLIP hydrogen-bond residues."
    )
    ap.add_argument("-i", "--inputs", nargs="+", required=True,
                    help="Pocket-prediction CSVs (several or wildcards), e.g. pockets/*.csv")
    ap.add_argument("-p", "--pockets", type=parse_pockets, required=True,
                    help="Pockets to keep, e.g. 1,2")
    ap.add_argument("--plip", nargs="*", default=None,
                    help="PLIP output files (xml/json/txt). One file applies to all CSVs; otherwise give one per CSV.")
    ap.add_argument("--plip_dir", default=None,
                    help="PLIP output directory; files are matched by name prefix (recommended for batches).")
    ap.add_argument("-o", "--outdir", default="pocket_hb_out",
                    help="Output directory")
    ap.add_argument("--format", choices=["tsv", "both"], default="both",
                    help="Write tsv or tsv+fasta (default: both)")
    ap.add_argument("--star", action="store_true",
                    help="Mark hydrogen-bond residues as AA* in FASTA; otherwise as lower-case aa")
    ap.add_argument("--dedup", action="store_true",
                    help="De-duplicate (chain, residue_label)")
    ap.add_argument("--min_prob", type=float, default=None,
                    help="Optional: keep probability >= min_prob")
    ap.add_argument("--min_score", type=float, default=None,
                    help="Optional: keep score >= min_score")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    csv_files = iter_files(args.inputs)
    if not csv_files:
        raise SystemExit("No CSV input found.")

    plip_files = iter_files(args.plip) if args.plip else []
    plip_dir = Path(args.plip_dir) if args.plip_dir else None

    # pairing:
    # 1) plip_dir: automatic matching
    # 2) a single --plip: shared by all CSVs
    # 3) as many --plip as CSVs: paired in order
    # otherwise: error
    pair_mode = None
    if plip_dir is not None:
        pair_mode = "dir"
    elif len(plip_files) == 1:
        pair_mode = "one"
    elif len(plip_files) == len(csv_files) and len(plip_files) > 1:
        pair_mode = "pairwise"
    elif len(plip_files) == 0:
        pair_mode = "none"
    else:
        raise SystemExit("PLIP file count mismatch: use --plip_dir, a single --plip, or one --plip per CSV.")

    for idx, csvf in enumerate(csv_files):
        csv_path = Path(csvf)
        base = csv_path.stem
        pocket_tag = "_".join(map(str, args.pockets))

        df = normalize_pocket_csv(pd.read_csv(csv_path))

        # pocket filter
        df = df[df["pocket"].isin(args.pockets)].copy()

        # threshold filter
        if args.min_prob is not None and "probability" in df.columns:
            df["probability"] = pd.to_numeric(df["probability"], errors="coerce")
            df = df[df["probability"] >= args.min_prob].copy()
        if args.min_score is not None and "score" in df.columns:
            df["score"] = pd.to_numeric(df["score"], errors="coerce")
            df = df[df["score"] >= args.min_score].copy()

        if args.dedup:
            df = df.drop_duplicates(subset=["chain", "residue_label"])

        if df.empty:
            sys.stderr.write(f"[WARN] {csv_path.name}: pocket={args.pockets} leaves no residues\n")
            continue

        # matching PLIP file
        plip_path = None
        if pair_mode == "dir":
            plip_path = find_plip_for_csv(csv_path, plip_dir)
        elif pair_mode == "one":
            plip_path = Path(plip_files[0])
        elif pair_mode == "pairwise":
            plip_path = Path(plip_files[idx])

        # parse the set of PLIP hydrogen-bond residues
        hb_set = set()
        hb_count = {}  # (chain,resnr) -> count
        hb_sources = {}  # (chain,resnr) -> set(paths)

        if plip_path is not None and plip_path.exists():
            hits = parse_plip_file(plip_path)
            for chain, resnr in hits:
                key = (str(chain).strip(), int(resnr))
                hb_set.add(key)
                hb_count[key] = hb_count.get(key, 0) + 1
                hb_sources.setdefault(key, set()).add(plip_path.name)

        # mark residues
        df["hb"] = df.apply(lambda r: (r["chain"], r["residue_label"]) in hb_set, axis=1)

        # extra summary columns (when a PLIP file exists)
        if hb_set:
            df["hb_count"] = df.apply(lambda r: hb_count.get((r["chain"], r["residue_label"]), 0), axis=1)
            df["hb_sources"] = df.apply(
                lambda r: ",".join(sorted(hb_sources.get((r["chain"], r["residue_label"]), set()))),
                axis=1
            )

        out_tsv = outdir / f"{base}.pocket{pocket_tag}.hb.tsv"
        out_fa = None
        if args.format == "both":
            out_fa = outdir / f"{base}.pocket{pocket_tag}.hb.fasta"

        write_outputs(df, out_tsv, out_fa, star=args.star)

        print(f"[OK] {csv_path.name} -> {out_tsv.name}" + (f" + {out_fa.name}" if out_fa else ""))

if __name__ == "__main__":
    import sys
    main()
