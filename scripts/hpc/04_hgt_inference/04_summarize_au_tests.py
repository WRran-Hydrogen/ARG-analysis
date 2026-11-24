#!/usr/bin/env python3
"""Collect AU-test results (03_run_au_test.sh) and add them to the primary HGT catalog.

For each task in au_test_tasks.tsv, reads au_tests/{event_id}.au.iqtree (tree 1 = unconstrained ML
tree, tree 2 = constrained "not nested" tree) and reports logL, deltaL and p-AU of the constrained
tree. The alternative is rejected when p-AU < 0.05.

Outputs ($SBL_HPC_BASE/04_hgt_inference/):
  AU_Test_Results.tsv
  Primary_HGT_Events_Catalog_2026.tsv (columns au_p_alternative, au_alternative_rejected, au_status)
  Stage3_HGT_Summary_Report.txt (AU section appended)
"""

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

csv.field_size_limit(10_000_000)
OUT_DIR = hpc.path_setting("SBL_HPC_HGT_DIR", hpc.BASE / "04_hgt_inference")
AU_DIR = OUT_DIR / "au_tests"
ALPHA = 0.05


def parse_au(path: Path) -> dict:
    """USER TREES table of an IQ-TREE -z ... -au run: per tree logL, deltaL, p-AU."""
    if not path.is_file():
        return {}
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    try:
        start = next(i for i, l in enumerate(lines) if "p-AU" in l and l.strip().startswith("Tree"))
    except StopIteration:
        return {}
    trees = {}
    for l in lines[start + 2:]:
        if not re.match(r"^\s*\d+\s", l):
            break
        nums = [t for t in l.split() if t not in ("+", "-")]
        trees[int(nums[0])] = {"logL": float(nums[1]), "deltaL": float(nums[2]), "p_AU": float(nums[-1])}
    return trees


def main() -> int:
    tasks = list(csv.DictReader(open(OUT_DIR / "au_test_tasks.tsv", encoding="utf-8"), delimiter="\t"))
    results = {}
    rows = []
    for t in tasks:
        res = parse_au(AU_DIR / f"{t['event_id']}.au.iqtree")
        if 1 in res and 2 in res:
            p = res[2]["p_AU"]
            row = {"event_id": t["event_id"], "family": t["family"], "status": "done",
                   "logL_ml": res[1]["logL"], "logL_constrained": res[2]["logL"], "deltaL_constrained": res[2]["deltaL"],
                   "p_AU_constrained": p, "alternative_rejected": p < ALPHA}
        else:
            row = {"event_id": t["event_id"], "family": t["family"], "status": "missing",
                   "logL_ml": "", "logL_constrained": "", "deltaL_constrained": "", "p_AU_constrained": "",
                   "alternative_rejected": ""}
        rows.append(row)
        results[t["event_id"]] = row

    with open(OUT_DIR / "AU_Test_Results.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["event_id"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    cat_path = OUT_DIR / "Primary_HGT_Events_Catalog_2026.tsv"
    cat = list(csv.DictReader(open(cat_path, encoding="utf-8"), delimiter="\t"))
    if cat:
        for r in cat:
            res = results.get(r["event_id"])
            r["au_status"] = res["status"] if res else ("not_tested" if r.get("is_tier1_hgt_prereg") == "True" else "")
            r["au_p_alternative"] = res["p_AU_constrained"] if res else ""
            r["au_alternative_rejected"] = res["alternative_rejected"] if res else ""
        with open(cat_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(cat[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(cat)

    done = [r for r in rows if r["status"] == "done"]
    rej = [r for r in done if r["alternative_rejected"]]
    lines = ["", "--- AU tests of pre-registered Tier-1 events (constrained: focal tips not nested in the donor background) ---",
             f"Tasks: {len(rows)}, finished: {len(done)}, alternative rejected (p-AU < {ALPHA}): {len(rej)}"]
    for r in rows:
        lines.append(f"  {r['event_id']}: {r['status']} p-AU={r['p_AU_constrained']} deltaL={r['deltaL_constrained']}")
    rep = OUT_DIR / "Stage3_HGT_Summary_Report.txt"
    txt = rep.read_text(encoding="utf-8") if rep.is_file() else ""
    txt = txt.split("\n--- AU tests of pre-registered Tier-1 events")[0].rstrip("\n") + "\n"
    rep.write_text(txt + "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
