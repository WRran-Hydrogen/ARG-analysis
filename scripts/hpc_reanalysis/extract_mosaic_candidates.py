#!/usr/bin/env python3
"""Extract unique CL0013/MOSAiC SBL candidate proteins from mosaic_virus.faa."""

from pathlib import Path

import hpc_config as hpc
import sys

MOSAIC_FAA = hpc.RAW_INPUT_DIR / "mosaic_virus.faa"
ORIG_TSV = hpc.RAW_INPUT_DIR / "sbl_mosaic.tsv"
ORIG_FAA = hpc.RAW_INPUT_DIR / "mosaic_virus_pf00144.faa"
WORK = hpc.BASE / "01_curated_datasets" / "mosaic_audit"
DOMTBL = WORK / "mosaic_cl0013.domtbl"
OUT_FAA = WORK / "mosaic_cl0013_candidates.faa"


def read_fa_dict(path):
    recs = {}
    cur_id, cur_hdr, buf = None, None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if cur_id is not None:
                    recs[cur_id] = (cur_hdr, "".join(buf))
                cur_hdr = s[1:].strip()
                cur_id = cur_hdr.split()[0]
                buf = []
            else:
                buf.append(s)
        if cur_id is not None:
            recs[cur_id] = (cur_hdr, "".join(buf))
    return recs


def main():
    cl_ids = set()
    with open(DOMTBL, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                cl_ids.add(line.split()[0])

    orig_tsv_ids = set()
    with open(ORIG_TSV, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                orig_tsv_ids.add(line.split()[0])

    orig_fa = read_fa_dict(ORIG_FAA)
    target_ids = cl_ids | orig_tsv_ids | set(orig_fa.keys())
    print(f"CL0013 GA hits: {len(cl_ids)}, orig_tsv: {len(orig_tsv_ids)}, orig_faa: {len(orig_fa)}, union: {len(target_ids)}")
    if len(target_ids) > 5000:
        print(f"ABORT_GT_5000_HITS: {len(target_ids)} > 5000", file=sys.stderr)
        sys.exit(2)

    found = {}
    cur_id, cur_hdr, buf = None, None, []
    keep = False
    with open(MOSAIC_FAA, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if keep and cur_id is not None:
                    found[cur_id] = (cur_hdr, "".join(buf))
                cur_hdr = s[1:].strip()
                cur_id = cur_hdr.split()[0]
                keep = cur_id in target_ids
                buf = []
            elif keep:
                buf.append(s)
        if keep and cur_id is not None:
            found[cur_id] = (cur_hdr, "".join(buf))

    for sid, val in orig_fa.items():
        if sid not in found:
            found[sid] = val

    with open(OUT_FAA, "w", encoding="utf-8") as f:
        for sid in sorted(found):
            hdr, seq = found[sid]
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")

    missing = target_ids - set(found)
    print(f"Extracted {len(found)} unique candidate sequences -> {OUT_FAA} (missing={len(missing)})")


if __name__ == "__main__":
    main()
