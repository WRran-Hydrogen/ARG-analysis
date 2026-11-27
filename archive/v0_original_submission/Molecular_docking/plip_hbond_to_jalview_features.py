#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, glob
from pathlib import Path
import pandas as pd
import xml.etree.ElementTree as ET

# Interaction element tags commonly found in PLIP XML reports
# hydrophobic_interaction / hydrogen_bond / salt_bridge entries carry complete fields.
PLIP_TAG_TO_FEATURETYPE = {
    "hydrogen_bond": "HBOND",
    "hydrophobic_interaction": "HYDROPHOBIC",
    "salt_bridge": "SALT_BRIDGE",
    "water_bridge": "WATER_BRIDGE",
    "pi_stack": "PI_STACK",
    "pi_cation_interaction": "PI_CATION",
    "halogen_bond": "HALOGEN_BOND",
    "metal_complex": "METAL_COMPLEX",
}

# Colour definitions for the Jalview features file (FeatureType -> color)
# Jalview lets each FeatureType define its display style at the top of the file as "type<TAB>color".
DEFAULT_COLORS = {
    "HBOND": "128,128,128",           # grey
    "HYDROPHOBIC": "255,200,0",       # yellow
    "SALT_BRIDGE": "255,0,0",         # red
    "WATER_BRIDGE": "0,180,255",      # blue
    "PI_STACK": "160,0,255",          # purple
    "PI_CATION": "255,0,160",         # pink
    "HALOGEN_BOND": "0,200,0",        # green
    "METAL_COMPLEX": "120,80,40",     # brown
}

def iter_files(items):
    files=[]
    for it in items:
        it=os.path.expanduser(it)
        m=glob.glob(it)
        if m: files.extend(m)
        elif os.path.exists(it): files.append(it)
    seen=set(); out=[]
    for f in files:
        if f not in seen:
            out.append(f); seen.add(f)
    return out

def load_seqmap(path: Path):
    """
    seqmap.tsv: seq_pos, chain, resnr, icode, resname3 ...
    Returns: (chain,resnr) -> seq_pos (first occurrence kept on duplicates)
    """
    df=pd.read_csv(path, sep="\t")
    mp={}
    for _,r in df.iterrows():
        key=(str(r["chain"]).strip(), int(r["resnr"]))
        mp.setdefault(key, int(r["seq_pos"]))
    return mp

def parse_plip_xml_all(xml_path: Path, enabled_featuretypes: set):
    """
    Parse a PLIP XML report and extract all interacting residues:
      Returns: list of dicts [{feature_type, chain, resnr, restype, extra_desc}, ...]
    Each PLIP interaction record usually has <resnr>, <reschain>, <restype> fields (hydrogen bonds, hydrophobic contacts and salt bridges alike).
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    hits = []
    for elem in root.iter():
        tag = elem.tag.lower().strip()
        if tag not in PLIP_TAG_TO_FEATURETYPE:
            continue
        ftype = PLIP_TAG_TO_FEATURETYPE[tag]
        if enabled_featuretypes and ftype not in enabled_featuretypes:
            continue

        # read child-node fields
        sub = {c.tag.lower().strip(): (c.text.strip() if c.text else "") for c in list(elem)}
        chain = sub.get("reschain")
        resnr = sub.get("resnr")
        restype = sub.get("restype")  # can be used in the description

        if not chain or not resnr:
            continue
        try:
            resnr_i = int(resnr)
        except:
            continue

        # extra information goes into the description (optional)
        # e.g. hbond has dist_h-a / don_angle; hydrophobic contacts and salt bridges have dist.
        extra = []
        for k in ("dist", "dist_h-a", "dist_d-a", "don_angle", "protisdon"):
            if k in sub and sub[k]:
                extra.append(f"{k}={sub[k]}")
        extra_desc = ";".join(extra) if extra else ""

        hits.append({
            "feature_type": ftype,
            "chain": chain.strip(),
            "resnr": resnr_i,
            "restype": (restype.strip() if restype else ""),
            "extra_desc": extra_desc,
        })

    return hits

def choose_enabled_types(types_arg: str):
    """
    --types accepts:
      all
      HBOND,HYDROPHOBIC,SALT_BRIDGE,...
      or tag names: hbond,hydrophobic,salt_bridge...
    """
    if not types_arg or types_arg.lower() == "all":
        return set(DEFAULT_COLORS.keys())  # all
    parts = [x.strip() for x in types_arg.split(",") if x.strip()]
    enabled=set()
    # accept lower-case short names
    alias = {
        "hbond":"HBOND",
        "hydrogen":"HBOND",
        "hydrophobic":"HYDROPHOBIC",
        "salt":"SALT_BRIDGE",
        "salt_bridge":"SALT_BRIDGE",
        "water":"WATER_BRIDGE",
        "water_bridge":"WATER_BRIDGE",
        "pi_stack":"PI_STACK",
        "pistack":"PI_STACK",
        "pi_cation":"PI_CATION",
        "pi_cation_interaction":"PI_CATION",
        "halogen":"HALOGEN_BOND",
        "halogen_bond":"HALOGEN_BOND",
        "metal":"METAL_COMPLEX",
        "metal_complex":"METAL_COMPLEX",
    }
    for p in parts:
        key = alias.get(p.lower(), p.upper())
        enabled.add(key)
    return enabled

def main():
    ap = argparse.ArgumentParser(
        description="Extract all (or selected) non-covalent interaction residues from PLIP XML, map them to sequence positions (seq_pos) and write a Jalview features file."
    )
    ap.add_argument("--seqmaps", nargs="+", required=True, help="seqmap.tsv files (from pdb_to_fasta_seqmap.py)")
    ap.add_argument("--plip", nargs="+", required=True, help="PLIP XML reports (wildcards allowed), e.g. *_report.xml")
    ap.add_argument("-o","--out", default="plip_noncovalent.features", help="Output Jalview features file")
    ap.add_argument("--match_by_prefix", action="store_true",
                    help="Pair files by name: a seqmap stem and a PLIP stem match when one contains the other")
    ap.add_argument("--types", default="all",
                    help="Interaction types to write: all, or comma-separated, e.g. HBOND,HYDROPHOBIC,SALT_BRIDGE")
    ap.add_argument("--merge_same_residue", action="store_true",
                    help="Merge several interactions of one residue into a single NONCOV feature (types listed in the description)")
    ap.add_argument("--merged_type_name", default="NONCOV", help="FeatureType name in merged mode")
    ap.add_argument("--merged_color", default="0,0,0", help="Colour in merged mode (default black)")
    args = ap.parse_args()

    enabled = choose_enabled_types(args.types)

    seqmap_files=[Path(x) for x in iter_files(args.seqmaps)]
    plip_files=[Path(x) for x in iter_files(args.plip)]

    seqmap_by={f.stem: f for f in seqmap_files}
    plip_by={f.stem: f for f in plip_files}

    pairs=[]
    if args.match_by_prefix:
        for sm_stem, sm in seqmap_by.items():
            pp = next((f for pstem,f in plip_by.items() if (sm_stem in pstem) or (pstem in sm_stem)), None)
            if pp:
                pairs.append((sm, pp))
    else:
        if len(seqmap_files) != len(plip_files):
            raise SystemExit("Without --match_by_prefix, --seqmaps and --plip must have the same number of files (paired in order).")
        pairs=list(zip(seqmap_files, plip_files))

    if not pairs:
        raise SystemExit("No seqmap/PLIP pairs matched; check the file names or use --match_by_prefix.")

    out_path=Path(args.out)
    lines=[]

    # 1) colour definition block (a Jalview features file may define each FeatureType style first).
    if args.merge_same_residue:
        lines.append(f"{args.merged_type_name}\t{args.merged_color}")
    else:
        for ft in sorted(enabled):
            color = DEFAULT_COLORS.get(ft, "128,128,128")
            lines.append(f"{ft}\t{color}")

    # 2) feature lines
    for sm_path, plip_path in pairs:
        # seq_id must match the first token of the alignment FASTA header, otherwise Jalview does not show it.
        seq_id = sm_path.stem.replace(".seqmap","")
        mp = load_seqmap(sm_path)

        if plip_path.suffix.lower() != ".xml":
            raise SystemExit(f"Only XML reports are supported: {plip_path} (TXT parsing of all types is unreliable)")

        hits = parse_plip_xml_all(plip_path, enabled_featuretypes=(set() if args.merge_same_residue else enabled))

        if args.merge_same_residue:
            # (chain,resnr)-> set(types)
            bucket={}
            extra_info={}
            for h in hits:
                key=(h["chain"], h["resnr"])
                bucket.setdefault(key, set()).add(h["feature_type"])
                if h["extra_desc"]:
                    extra_info.setdefault(key, set()).add(h["extra_desc"])
            for (chain,resnr), typeset in sorted(bucket.items()):
                if (chain,resnr) not in mp:
                    continue
                pos = mp[(chain,resnr)]
                desc = f"PLIP_NONCOV;chain={chain};resnr={resnr};types={','.join(sorted(typeset))}"
                if extra_info.get((chain,resnr)):
                    desc += f";info={'|'.join(sorted(extra_info[(chain,resnr)]))}"
                # Jalview feature line:description sequenceId sequenceIndex start end featureType [score]
                lines.append(f"{desc}\t{seq_id}\t-1\t{pos}\t{pos}\t{args.merged_type_name}")
        else:
            # one line per interaction type
            for h in hits:
                chain,resnr,ftype = h["chain"], h["resnr"], h["feature_type"]
                if ftype not in enabled:
                    continue
                if (chain,resnr) not in mp:
                    continue
                pos = mp[(chain,resnr)]
                restype = h["restype"]
                extra = h["extra_desc"]
                desc = f"PLIP_{ftype};chain={chain};resnr={resnr}"
                if restype:
                    desc += f";aa={restype}"
                if extra:
                    desc += f";{extra}"
                lines.append(f"{desc}\t{seq_id}\t-1\t{pos}\t{pos}\t{ftype}")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[OK] wrote Jalview features -> {out_path}")

if __name__=="__main__":
    main()
