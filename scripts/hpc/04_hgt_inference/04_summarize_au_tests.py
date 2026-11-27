#!/usr/bin/env python3
"""Collect AU-test results (03_run_au_test.sh) and add them to the primary HGT catalog.

For each task in au_test_tasks.tsv, reads au_tests/{event_id}.au.iqtree (tree 1 = unconstrained ML
tree, tree 2 = constrained "not nested" tree) and reports logL, deltaL and p-AU of the constrained
tree. The alternative is rejected when p-AU < 0.05 and the ML tree is more likely than the constrained
tree; a constrained tree with the higher likelihood (also as the score of its own search) is reported as
'ml_tree_suboptimal' (the unconstrained search missed a better tree; the nesting is not supported).

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


def best_score(log: Path) -> float | None:
    """Final 'BEST SCORE FOUND' of an IQ-TREE search log."""
    if not log.is_file():
        return None
    m = re.findall(r"BEST SCORE FOUND\s*:\s*(-?[0-9.]+)", log.read_text(encoding="utf-8", errors="ignore"))
    return float(m[-1]) if m else None


AU_FIELDS = ["status", "logL_ml", "logL_constrained", "logL_constrained_search", "deltaL_constrained",
             "p_AU_constrained", "alternative_rejected"]


def au_result(prefix: Path) -> dict:
    """Result of one AU test ({prefix}.au.iqtree: tree 1 = unconstrained ML tree, tree 2 = constrained tree).

    deltaL_constrained = logL(ML tree) - logL(constrained tree), signed. When the constrained tree has the
    higher likelihood (in the AU table, or as the score of its own constrained search), the unconstrained
    search did not find the ML tree: status 'ml_tree_suboptimal' and the alternative is not rejected.
    """
    res = parse_au(Path(f"{prefix}.au.iqtree"))
    if 1 not in res or 2 not in res:
        return {k: "" for k in AU_FIELDS} | {"status": "missing"}
    ml, co = res[1]["logL"], res[2]["logL"]
    cs = best_score(Path(f"{prefix}.constr.log"))
    subopt = co > ml or (cs is not None and cs > ml)
    p = res[2]["p_AU"]
    return {"status": "ml_tree_suboptimal" if subopt else "done", "logL_ml": ml, "logL_constrained": co,
            "logL_constrained_search": cs if cs is not None else "", "deltaL_constrained": round(ml - co, 4),
            "p_AU_constrained": p, "alternative_rejected": (not subopt) and p < ALPHA}


def main() -> int:
    tasks = list(csv.DictReader(open(OUT_DIR / "au_test_tasks.tsv", encoding="utf-8"), delimiter="\t"))
    results = {}
    rows = []
    for t in tasks:
        row = {"event_id": t["event_id"], "family": t["family"]} | au_result(AU_DIR / t["event_id"])
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

    done = [r for r in rows if r["status"] in ("done", "ml_tree_suboptimal")]
    rej = [r for r in done if r["alternative_rejected"]]
    sub = [r for r in done if r["status"] == "ml_tree_suboptimal"]
    lines = ["", "--- AU tests of pre-registered Tier-1 events (constrained: focal tips not nested in the donor background) ---",
             f"Tasks: {len(rows)}, finished: {len(done)}, alternative rejected (p-AU < {ALPHA}): {len(rej)}, "
             f"constrained tree more likely than the ML tree (ml_tree_suboptimal, nesting not supported): {len(sub)}",
             "deltaL = logL(ML tree) - logL(constrained tree); model parameters fitted on the ML tree"]
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
