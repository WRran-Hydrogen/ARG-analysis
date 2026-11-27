#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, glob, os
from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser

AA3_TO_AA1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "SEC":"U","PYL":"O","ASX":"B","GLX":"Z","XLE":"J","UNK":"X"
}

def aa3_to_aa1(resname: str) -> str:
    return AA3_TO_AA1.get(resname.strip().upper(), "X")

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

def main():
    ap = argparse.ArgumentParser(description="Extract protein chain sequences from PDB files to FASTA and write a sequence-position to PDB residue-number map (seqmap.tsv).")
    ap.add_argument("-i","--inputs",nargs="+",required=True,help="Input PDB files (wildcards allowed)")
    ap.add_argument("-c","--chain",default=None,help="Chain ID (e.g. A/B). If omitted, every chain is written separately")
    ap.add_argument("-o","--outdir",default="seqmaps",help="Output directory")
    ap.add_argument("--keep_x",action="store_true",help="Keep unknown residues as X (default: drop X)")
    args = ap.parse_args()

    outdir=Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    parser=PDBParser(QUIET=True)

    for pdb in iter_files(args.inputs):
        pdbp=Path(pdb)
        struct=parser.get_structure(pdbp.stem, str(pdbp))
        model=struct[0]

        chains = [model[args.chain]] if args.chain else list(model.get_chains())

        for ch in chains:
            chain_id=ch.id
            seq=[]
            rows=[]
            seq_pos=0

            for res in ch.get_residues():
                hetflag, resseq, icode = res.id
                if hetflag.strip():  # skip hetero residues / water
                    continue
                aa1=aa3_to_aa1(res.get_resname())
                if aa1=="X" and not args.keep_x:
                    continue
                seq_pos += 1
                seq.append(aa1)
                rows.append({
                    "seq_pos": seq_pos,
                    "chain": chain_id,
                    "resnr": int(resseq),
                    "icode": (icode.strip() if isinstance(icode,str) else ""),
                    "resname3": res.get_resname().strip().upper()
                })

            if not seq:
                continue

            # keep IDs free of spaces; Jalview uses the first header word as sequenceId
            seq_id = f"{pdbp.stem}_chain{chain_id}"
            fasta_path = outdir / f"{seq_id}.fasta"
            map_path   = outdir / f"{seq_id}.seqmap.tsv"

            fasta_path.write_text(f">{seq_id}\n{''.join(seq)}\n", encoding="utf-8")
            pd.DataFrame(rows).to_csv(map_path, sep="\t", index=False)
            print(f"[OK] {pdbp.name} chain {chain_id} -> {fasta_path.name}, {map_path.name}")

if __name__=="__main__":
    main()
