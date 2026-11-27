#!/usr/bin/env python3
"""Protein-architecture classes of the curated viral SBLs (Fig. 1C input).

Each curated viral SBL (HMM coverage >= 65 %, FIMO motif pass) is assigned one class:
  Verified hydrolytic SBL      enzymatic_validation_class == Verified_Hydrolytic_SBL
  Signalling / sensor          BlaR1 sensor class, or a co-domain of that module
  Virion structure             co-domain of a tail spike / tail fibre / tape measure protein
  Lysis / peptidoglycan hydrolysis
  Cell-surface binding
  Cell-wall synthesis (PBP)    complete multi-domain penicillin-binding protein
  Other fusion                 any other co-domain (DUFs, unrelated enzymes)
  SBL domain only              no co-domain in the protein
Co-domains are the Pfam domains other than the SBL domain itself listed in domain_architecture; their module comes
from data/metadata/SBL_CoDomain_Module_Map.tsv (Linker domains are ignored). With co-domains of several modules,
PRIORITY decides (virion structure > lysis > cell-wall synthesis > signalling > cell-surface binding > other).
Truncated proteins and pseudogenes keep their architecture class.

Outputs (results/reanalysis_2026/06_structure_and_ecology/)
  Viral_SBL_Architecture_Classes_2026.tsv   one row per curated viral SBL
  Viral_SBL_Architecture_Summary_2026.tsv   counts per family x class (main_text = False for PF00768)
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

if sys.platform.startswith("win"):
    try:
        import pandas  # noqa: F401
    except ImportError:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import sbl_config as _cfg

        sys.exit(subprocess.run(_cfg.linux_cmd(_cfg.WSL_PY, _cfg.to_wsl_path(Path(__file__).resolve()),
                                               *sys.argv[1:])).returncode)

import pandas as pd

BASE = Path(__file__).resolve().parents[2]
RES = BASE / "results" / "reanalysis_2026"
VIRAL_CAT = RES / "01_curation_and_context" / "Final_Viral_SBL_Catalog.tsv"
MODULE_MAP = BASE / "data" / "metadata" / "SBL_CoDomain_Module_Map.tsv"
OUT_DIR = RES / "06_structure_and_ecology"
FAMILIES = ["PF00144", "PF00905", "PF13354", "PF00768"]
MAIN_FAMILIES = ["PF00144", "PF00905", "PF13354"]
CLASS_ORDER = ["Verified hydrolytic SBL", "Signalling / sensor", "Virion structure", "Lysis / peptidoglycan hydrolysis",
               "Cell-surface binding", "Cell-wall synthesis (PBP)", "Other fusion", "SBL domain only"]
# module precedence when a protein carries co-domains of several modules (a kinase next to a tape measure protein
# stays a virion protein; a surface Ig domain on a complete PBP stays a PBP)
PRIORITY = ["Virion structure", "Lysis / peptidoglycan hydrolysis", "Cell-wall synthesis (PBP)", "Signalling / sensor",
            "Cell-surface binding", "Other fusion"]


def classify(fam: str, arch: str, enz_class: str, modules: dict[str, str]) -> tuple[str, str]:
    accs = sorted(set(re.findall(r"(PF\d{5})", arch or "")) - {fam})
    mods = {modules.get(a, "Other fusion") for a in accs} - {"Linker"}
    if enz_class == "Verified_Hydrolytic_SBL":
        return CLASS_ORDER[0], ";".join(accs)
    if enz_class == "BlaR1_Sensor_Receptor":
        return "Signalling / sensor", ";".join(accs)
    for c in PRIORITY:
        if c in mods:
            return c, ";".join(accs)
    return ("SBL domain only" if not mods else "Other fusion"), ";".join(accs)


def main() -> None:
    modules = dict(pd.read_csv(MODULE_MAP, sep="\t")[["acc", "module"]].itertuples(index=False))
    cat = pd.read_csv(VIRAL_CAT, sep="\t", low_memory=False)
    cur = cat[(cat["pass_cov65"] == True) & (cat["pass_fimo_1e5"] == True)  # noqa: E712
              & cat["winning_acc"].isin(FAMILIES)].drop_duplicates("uid").copy()
    res = [classify(f, a, e, modules) for f, a, e in
           zip(cur["winning_acc"], cur["domain_architecture"].fillna(""), cur["enzymatic_validation_class"].fillna(""))]
    cur["architecture_class"] = [r[0] for r in res]
    cur["co_domains"] = [r[1] for r in res]
    out = cur[["uid", "source", "winning_acc", "viral_origin_tier", "enzymatic_validation_class", "domain_architecture",
               "co_domains", "architecture_class"]].rename(columns={"winning_acc": "family"})
    out["main_text"] = out["family"] != "PF00768"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_DIR / "Viral_SBL_Architecture_Classes_2026.tsv", sep="\t", index=False)
    summ = (out.groupby(["family", "architecture_class"]).size().rename("n").reset_index()
              .merge(out.groupby("family").size().rename("total").reset_index(), on="family"))
    summ["share"] = (summ["n"] / summ["total"]).round(4)
    summ["main_text"] = summ["family"] != "PF00768"
    summ.to_csv(OUT_DIR / "Viral_SBL_Architecture_Summary_2026.tsv", sep="\t", index=False)
    tab = pd.crosstab(out["architecture_class"], out["family"]).reindex(CLASS_ORDER).fillna(0).astype(int)
    print(tab[FAMILIES].to_string())
    unmapped = sorted({a for s in out["co_domains"] for a in s.split(";") if a and a not in modules})
    print(f"\nco-domains without a module (counted as Other fusion): {len(unmapped)}")


if __name__ == "__main__":
    main()
