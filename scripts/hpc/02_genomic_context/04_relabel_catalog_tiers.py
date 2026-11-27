#!/usr/bin/env python3
"""Update `viral_origin_tier` in the curated viral catalogs from the current Phase 2 table, keeping
every other column, the viral nr85 representatives and all FASTA exports unchanged.

Used when the tier rule changes after the trees have been built (2026-11-25): the nr85 representative
of each cluster was chosen with tier priority (01_curation/09 --stage finalize), so re-running finalize
could change the representatives and therefore the trees. This script relabels only, and reports how
many clusters would have chosen a representative of a better tier under the new labels.

Catalogs updated in place (a copy of each original is written next to it with --backup-suffix):
  $SBL_HPC_BASE/01_final_curated/Final_Viral_SBL_Catalog.tsv
  $SBL_HPC_BASE/01_final_curated/Viral_Enzymatic_Function_Verification_2026.tsv
Report: $SBL_HPC_BASE/02_genomic_context/Catalog_Tier_Relabel_Report.txt
"""
import argparse
import csv
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

csv.field_size_limit(2**31 - 1)
P2 = hpc.BASE / "02_genomic_context" / "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"
CATALOGS = (hpc.BASE / "01_final_curated" / "Final_Viral_SBL_Catalog.tsv",
            hpc.BASE / "01_final_curated" / "Viral_Enzymatic_Function_Verification_2026.tsv")
PRIORITY = {"Tier1_HighConfidence_Viral": 0, "MOSAiC_Proteome_Verified_vMAG": 1,
            "Tier2_Supported_Viral": 2, "Tier3_Fragment_or_Borderline": 3}
SHORT = {"Tier1_HighConfidence_Viral": "T1", "Tier2_Supported_Viral": "T2",
         "Tier3_Fragment_or_Borderline": "T3", "MOSAiC_Proteome_Verified_vMAG": "MOSAiC"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--phase2", default=str(P2))
    ap.add_argument("--backup-suffix", default=".before_tier_relabel")
    ap.add_argument("--report", default=str(hpc.BASE / "02_genomic_context" / "Catalog_Tier_Relabel_Report.txt"))
    a = ap.parse_args()

    with open(a.phase2, encoding="utf-8") as f:
        p2 = list(csv.DictReader(f, delimiter="\t"))
    if not p2 or "viral_origin_tier_local_w10" not in p2[0]:
        sys.exit(f"{a.phase2} is not a Phase 2 table with gene-origin evidence tiers; run 01_build_neighborhoods.py first")
    new = {r["uid"]: r["viral_origin_tier"] for r in p2}

    p2_path = Path(a.phase2)
    try:
        p2_label = "$SBL_HPC_BASE/" + p2_path.resolve().relative_to(hpc.BASE.resolve()).as_posix()
    except ValueError:
        p2_label = p2_path.name
    lines = [f"Phase 2 table: {p2_label}"]
    for path in CATALOGS:
        with open(path, encoding="utf-8") as f:
            rd = csv.DictReader(f, delimiter="\t")
            fields, rows = rd.fieldnames, list(rd)
        backup = path.with_name(path.name + a.backup_suffix)
        if not backup.exists():
            shutil.copy2(path, backup)
        shift, missing = Counter(), 0
        for r in rows:
            old = r["viral_origin_tier"]
            nt = new.get(r["uid"])
            if nt is None:
                missing += 1
                continue
            r["viral_origin_tier"] = nt
            shift[(old, nt)] += 1
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")   # csv default line ending, as 01_curation/09 writes
            w.writeheader()
            w.writerows(rows)
        lines += ["", f"{path.name}: {len(rows):,} rows; relabelled {sum(shift.values()):,}; not in Phase 2 (kept) {missing:,}; "
                  f"backup {backup.name}", "  old -> new: " + ", ".join(
                      f"{SHORT.get(o, o)}->{SHORT.get(n, n)} {c:,}" for (o, n), c in sorted(shift.items()) if c)]

        if "viral_nr85_rep_uid" in fields:
            members = defaultdict(list)
            for r in rows:
                if r.get("viral_nr85_rep_uid"):
                    members[r["viral_nr85_rep_uid"]].append(r)
            tier = {r["uid"]: r["viral_origin_tier"] for r in rows}
            worse = [rep for rep, ms in members.items()
                     if PRIORITY.get(tier.get(rep, ""), 9) > min(PRIORITY.get(m["viral_origin_tier"], 9) for m in ms)]
            reps = [r for r in rows if r.get("in_viral_nr85_rep", "").lower() == "true"]
            lines += [f"  nr85 clusters: {len(members):,}; representatives kept fixed; clusters where another member now has a "
                      f"better tier than the representative: {len(worse):,}",
                      "  representative tiers (single-family exports): " + ", ".join(
                          f"{SHORT.get(k, k)} {v:,}" for k, v in sorted(Counter(r['viral_origin_tier'] for r in reps).items()))]
    rep = "\n".join(lines) + "\n"
    Path(a.report).write_text(rep, encoding="utf-8")
    print(rep)


if __name__ == "__main__":
    main()
