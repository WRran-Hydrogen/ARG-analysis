#!/usr/bin/env python3
"""Selection analysis, step 5 (sensitivity): MEME and FEL on GARD partitions.

GARD breakpoints are found on nucleotides, so each breakpoint is moved to the nearest codon boundary; a
segment that shrinks to nothing is merged into its left neighbour. Each partition gets the GARD segment
topology (branch lengths and node labels removed; HyPhy re-fits branch lengths).

  build <set_dir>   write <set_dir>/gard_partitioned.nex from codon.fna and gard.json
  tasks             write gard_partition_tasks.tsv (subclade sets whose GARD run is 'complete')
  summarize         compare single-tree and partitioned MEME/FEL site calls and write
                    GARD_Partitioned_MEME_FEL_Sensitivity_2026.tsv (per set) and
                    GARD_Partitioned_MEME_FEL_Sites_2026.tsv (per site)
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

OUT_DIR = hpc.BASE / "05_selection_pressure"
TASKS_TSV = OUT_DIR / "hyphy_tasks.tsv"
PART_TASKS_TSV = OUT_DIR / "gard_partition_tasks.tsv"
SET_TSV = OUT_DIR / "GARD_Partitioned_MEME_FEL_Sensitivity_2026.tsv"
SITE_TSV = OUT_DIR / "GARD_Partitioned_MEME_FEL_Sites_2026.tsv"
ALPHA = 0.05


def read_fasta(path: Path) -> dict[str, str]:
    seqs, key = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            key = line[1:].split()[0]
            seqs[key] = []
        elif key:
            seqs[key].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def topology_only(newick: str) -> str:
    t = re.sub(r"\)[A-Za-z_][A-Za-z0-9_.]*", ")", newick)  # internal node labels
    t = re.sub(r":-?[0-9.]+(?:[eE][-+]?[0-9]+)?", "", t)  # branch lengths
    return t.rstrip(";") + ";"


def codon_partitions(gard: dict, n_codons: int) -> list[tuple[int, int, str]]:
    """GARD segments as (first_codon, last_codon, tree), 1-based and contiguous."""
    segs = sorted(gard["breakpointData"].values(), key=lambda v: v["bps"][0][0])
    out: list[tuple[int, int, str]] = []
    start = 1
    for i, seg in enumerate(segs):
        end = n_codons if i == len(segs) - 1 else min(n_codons, round(seg["bps"][0][1] / 3))
        if end < start:  # segment vanished after snapping: extend the previous partition instead
            continue
        out.append((start, end, topology_only(seg["tree"])))
        start = end + 1
    if out and out[-1][1] < n_codons:
        s, _, t = out[-1]
        out[-1] = (s, n_codons, t)
    return out


def build(set_dir: Path) -> None:
    seqs = read_fasta(set_dir / "codon.fna")
    nchar = len(next(iter(seqs.values())))
    parts = codon_partitions(json.loads((set_dir / "gard.json").read_text(encoding="utf-8")), nchar // 3)
    lines = ["#NEXUS", "", "BEGIN TAXA;", f"\tDIMENSIONS NTAX = {len(seqs)};",
             "\tTAXLABELS " + " ".join(f"'{k}'" for k in seqs) + " ;", "END;", "",
             "BEGIN CHARACTERS;", f"\tDIMENSIONS NCHAR = {nchar};",
             "\tFORMAT DATATYPE = DNA GAP=- MISSING=? ;", "MATRIX"]
    lines += [f"\t'{k}' {v}" for k, v in seqs.items()]
    lines += [";", "END;", "", "BEGIN ASSUMPTIONS;"]
    lines += [f"\tCHARSET span_{i} = {3 * (s - 1) + 1}-{3 * e};" for i, (s, e, _) in enumerate(parts, 1)]
    lines += ["END;", "", "BEGIN TREES;"]
    lines += [f"\tTREE tree_{i} = {t}" for i, (_, _, t) in enumerate(parts, 1)]
    lines += ["END;", ""]
    (set_dir / "gard_partitioned.nex").write_text("\n".join(lines), encoding="utf-8")
    print(f"{set_dir.name}: {len(parts)} partitions, codons " + ", ".join(f"{s}-{e}" for s, e, _ in parts))


def write_tasks() -> None:
    with open(TASKS_TSV, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    keep = []
    for r in rows:
        st = Path(r["set_dir"]) / "gard.status"
        if r["set_group"] == "subclade" and st.is_file() and st.read_text().strip() == "complete":
            keep.append(r)
    with open(PART_TASKS_TSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(keep)
    print(f"{len(keep)} subclade sets with a complete GARD run -> {PART_TASKS_TSV}")


def site_table(path: Path) -> list[dict]:
    """MEME/FEL MLE rows in alignment order (partitions are contiguous, so they are concatenated)."""
    data = json.loads(path.read_text(encoding="utf-8"))
    headers = [h[0] for h in data["MLE"]["headers"]]
    content = data["MLE"]["content"]
    rows = []
    for key in sorted(content, key=int):
        rows += [dict(zip(headers, map(float, r))) for r in content[key]]
    return rows


def calls(rows: list[dict], kind: str) -> dict[str, set[int]]:
    if kind == "meme":
        return {"pos": {i for i, r in enumerate(rows, 1) if r["p-value"] < ALPHA}}
    return {"pos": {i for i, r in enumerate(rows, 1) if r["p-value"] < ALPHA and r["beta"] > r["alpha"]},
            "pur": {i for i, r in enumerate(rows, 1) if r["p-value"] < ALPHA and r["alpha"] > r["beta"]}}


def summarize() -> None:
    with open(PART_TASKS_TSV, encoding="utf-8") as fh:
        tasks = list(csv.DictReader(fh, delimiter="\t"))
    set_rows, site_rows = [], []
    for t in tasks:
        d = Path(t["set_dir"])
        pairs = {k: (d / f"{k}.json", d / f"{k}_gard.json") for k in ("meme", "fel")}
        ok = {k: all(p.is_file() and p.stat().st_size for p in ab) for k, ab in pairs.items()}
        row = {"set_id": t["set_id"], "family": t["family"], "n_codons": t["n_codons"],
               "meme_status": "done" if ok["meme"] else "failed", "fel_status": "done" if ok["fel"] else "failed"}
        if not any(ok.values()):
            set_rows.append(row)
            continue
        row["n_partitions"] = len(json.loads(pairs["fel" if ok["fel"] else "meme"][1].read_text(encoding="utf-8"))["MLE"]["content"])
        tabs = {}
        for kind, (single_p, part_p) in pairs.items():
            if not ok[kind]:
                continue
            single, part = site_table(single_p), site_table(part_p)
            if len(single) != len(part):
                raise SystemExit(f"{t['set_id']} {kind}: {len(single)} vs {len(part)} sites")
            tabs[kind] = (single, part)
            cs, cp = calls(single, kind), calls(part, kind)
            for lab in cs:
                a, b = cs[lab], cp[lab]
                row[f"{kind}_{lab}_single"] = len(a)
                row[f"{kind}_{lab}_partitioned"] = len(b)
                row[f"{kind}_{lab}_shared"] = len(a & b)
                row[f"{kind}_{lab}_retained_pct"] = round(100.0 * len(a & b) / len(a), 1) if a else ""
        set_rows.append(row)
        ms, mp = tabs.get("meme", (None, None))
        fs, fp = tabs.get("fel", (None, None))
        n_sites = len(ms or fs)
        for i in range(n_sites):
            site_rows.append({"set_id": t["set_id"], "family": t["family"], "hyphy_site_1based": i + 1,
                              "meme_p_single": ms[i]["p-value"] if ms else "",
                              "meme_p_partitioned": mp[i]["p-value"] if mp else "",
                              "fel_alpha_single": fs[i]["alpha"] if fs else "",
                              "fel_beta_single": fs[i]["beta"] if fs else "",
                              "fel_p_single": fs[i]["p-value"] if fs else "",
                              "fel_alpha_partitioned": fp[i]["alpha"] if fp else "",
                              "fel_beta_partitioned": fp[i]["beta"] if fp else "",
                              "fel_p_partitioned": fp[i]["p-value"] if fp else ""})
    fields = list(dict.fromkeys(k for r in set_rows for k in r))
    for path, rows, fn in ((SET_TSV, set_rows, fields), (SITE_TSV, site_rows, list(site_rows[0]) if site_rows else [])):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fn, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(rows)
    for r in set_rows:
        print("\t".join(str(r.get(k, "")) for k in fields))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "build" and len(sys.argv) == 3:
        build(Path(sys.argv[2]))
    elif cmd == "tasks":
        write_tasks()
    elif cmd == "summarize":
        summarize()
    else:
        raise SystemExit(__doc__)
