#!/usr/bin/env python3
"""Habitat enrichment of SBL-carrying viruses by viral phylum, and habitat -> viral phylum -> host flows
(Fig. 1C and 1D inputs).

Habitats: the IMG/VR v4.1 ecosystem path (first three levels) is collapsed into the 14 categories used in the
original submission (habitat_category()). Viral phylum and host phylum are the IMG/VR v4.1 geNomad taxonomy and
host prediction; GTDB r207 host phyla are given their current names (Firmicutes* -> Bacillota,
Proteobacteria -> Pseudomonadota, Cyanobacteria -> Cyanobacteriota).

Enrichment (one 2x2 table per viral phylum x habitat, UViG level, high-confidence UViGs only, because the
background table holds IMG/VR high-confidence UViGs):
    a = SBL-carrying UViGs of the phylum in the habitat      b = SBL-carrying UViGs of the phylum elsewhere
    c = other UViGs of the phylum in the habitat             d = other UViGs of the phylum elsewhere
    log2 enrichment = log2( (a / (a + b)) / ((a + c) / (a + b + c + d)) ), left empty when a = 0 (no SBL-carrying
    UViG observed; the expected count is reported instead, because a pseudocount would turn small expected
    counts into apparent enrichment); two-sided Fisher exact test; Benjamini-Hochberg across all tested cells of the pooled analysis, and
    separately within each Pfam family. Group 'All_main' pools UViGs carrying any of the three main-text families
    (PF00144, PF00905, PF13354); 'Hydrolytic_main' pools UViGs carrying a verified hydrolytic SBL of those families;
    PF00768 rows (main_text = False) are kept for the supplement.
A UViG is SBL-carrying (for a family) when it holds at least one curated viral SBL row (HMM coverage >= 65 %,
FIMO motif pass) of that family.
Architecture states (Fig. 1B): the six protein-architecture states of the ancestral state reconstruction, from
Viral_SBL_Architecture_Classes_2026.tsv; each state is tested against the IMG/VR background (UViG level, and with at
most one UViG per metagenome) and against the other main-family SBL-carrying UViGs of the same phylum.
"Respiratory" (ecosystem level 3 = Respiratory system) is split from "Other host-associated" because the tail-spike
fusions concentrate there.

Inputs
  results/reanalysis_2026/01_curation_and_context/Final_Viral_SBL_Catalog.tsv
  data/metadata/information.xlsx (sheet ALL_information: IMG/VR v4.1 records of the SBL-carrying UViGs)
  data/metadata/IMGVR_v4.1_HighConfidence_Phylum_Ecosystem_Counts.tsv
      UViG counts of IMGVR_all_Sequence_information-high_confidence.tsv by viral phylum and ecosystem levels 1-3
  data/metadata/GOLD_TaxonOID_Coordinates.tsv.gz (from build_sample_geography_2026.py; habitat of MOSAiC vMAGs)
  results/reanalysis_2026/01_curation_and_context/Phase2_Viral_Genomic_Neighborhoods_Full.tsv (MOSAiC phylum)
  results/reanalysis_2026/06_structure_and_ecology/Viral_SBL_Architecture_Classes_2026.tsv
      (from classify_sbl_architecture_2026.py)
Outputs (results/reanalysis_2026/06_structure_and_ecology/)
  Habitat_Enrichment_By_Viral_Phylum_2026.tsv
  Habitat_Enrichment_By_Architecture_State_2026.tsv   per state: comparison (vs_background / vs_other_SBL) x unit
  Habitat_ViralPhylum_Host_Flows_2026.tsv
  Family_Function_Host_Composition_2026.tsv   functional class of the curated rows (Other = candidate, BlaR1
      sensor, truncated) and host phylum of the SBL-carrying UViGs, per family
  Hydrolytic_SBL_Profile_2026.tsv   one row per verified hydrolytic SBL (main families): CARD enzyme, viral origin
      tier, IMG/VR confidence, viral phylum, host phylum and class, habitat
"""

from __future__ import annotations

import math
import re
import subprocess
import sys
from pathlib import Path

if sys.platform.startswith("win"):
    try:
        import openpyxl  # noqa: F401
        import pandas  # noqa: F401
        import scipy.stats  # noqa: F401
    except ImportError:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import sbl_config as _cfg

        sys.exit(subprocess.run(_cfg.linux_cmd(_cfg.WSL_PY, _cfg.to_wsl_path(Path(__file__).resolve()),
                                               *sys.argv[1:])).returncode)

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

BASE = Path(__file__).resolve().parents[2]
META = BASE / "data" / "metadata"
RES = BASE / "results" / "reanalysis_2026"
VIRAL_CAT = RES / "01_curation_and_context" / "Final_Viral_SBL_Catalog.tsv"
INFO_XLSX = META / "information.xlsx"
BG_COUNTS = META / "IMGVR_v4.1_HighConfidence_Phylum_Ecosystem_Counts.tsv"
GOLD_CACHE = META / "GOLD_TaxonOID_Coordinates.tsv.gz"  # written by build_sample_geography_2026.py
NEIGHBORHOODS = RES / "01_curation_and_context" / "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"
OUT_DIR = RES / "06_structure_and_ecology"
ARCH_CLASSES = OUT_DIR / "Viral_SBL_Architecture_Classes_2026.tsv"  # written by classify_sbl_architecture_2026.py
# eight architecture classes -> the six states of the ancestral state reconstruction (asr_architecture_tips_2026.py)
ARCH_STATE = {"Verified hydrolytic SBL": "Hydrolytic", "Virion structure": "Virion_fusion",
              "Lysis / peptidoglycan hydrolysis": "Surface_lysis_fusion", "Cell-surface binding": "Surface_lysis_fusion",
              "Cell-wall synthesis (PBP)": "Complete_PBP", "SBL domain only": "SBL_only",
              "Signalling / sensor": "Other", "Other fusion": "Other"}
STATES = ["Hydrolytic", "Virion_fusion", "Surface_lysis_fusion", "Complete_PBP", "SBL_only", "Other"]
FAMILIES = ["PF00144", "PF00905", "PF13354", "PF00768"]
# Main-text families; PF00768 (homologous PBPs split from PF13354 by the clan bitscore assignment) is reported in
# the supplement only. "All_main" pools the three main families.
MAIN_FAMILIES = ["PF00144", "PF00905", "PF13354"]
HYDROLYTIC = "Verified_Hydrolytic_SBL"
FUNC_GROUP = {"Verified_Hydrolytic_SBL": "Verified hydrolytic SBL",
              "PBP_Esterase_or_PhageLysin_Homolog": "PBP / esterase / lysin homolog"}  # the rest: "Other"
MAJOR_HOSTS = ["Bacillota", "Actinobacteriota", "Bacteroidota", "Pseudomonadota", "Cyanobacteriota"]
HABITATS = ["Digestive", "Respiratory", "Other host-associated", "Plants", "Microbial", "Freshwater", "Marine", "Thermal", "Soil",
            "Other environmental", "Bioreactor", "Wastewater", "Solid waste", "Other engineered", "Unclassified"]
# Catch-all categories are reported but not shown as rows of the enrichment heat map in the original figure
SHOWN_HABITATS = [h for h in HABITATS if not h.startswith("Other") and h != "Unclassified"]


def habitat_category(path: object) -> str:
    """IMG/VR ecosystem path (level1;level2;level3[;...]) -> one of HABITATS."""
    p = (str(path).split(";") + ["", "", ""])[:3] if isinstance(path, str) else ["", "", ""]
    l1, l2, l3 = (s.strip() for s in p)
    if not l1:
        return "Unclassified"
    if l1 == "Host-associated":
        if "Digestive system" in l3:
            return "Digestive"
        if "Respiratory system" in l3:
            return "Respiratory"
        if l2 == "Plants":
            return "Plants"
        if l2 == "Microbial":
            return "Microbial"
        return "Other host-associated"
    if l1 == "Environmental":
        if l2 == "Aquatic":
            return {"Freshwater": "Freshwater", "Marine": "Marine", "Thermal springs": "Thermal"}.get(l3, "Other environmental")
        if l2 == "Terrestrial" and l3 == "Soil":
            return "Soil"
        return "Other environmental"
    if l1 == "Engineered":
        if l2 == "Bioreactor":
            return "Bioreactor"
        if l2 in ("Wastewater", "WWTP"):
            return "Wastewater"
        if l2 == "Solid waste":
            return "Solid waste"
        if l2 == "Artificial ecosystem" and re.search(r"soil|plant", l3, re.I):
            return "Soil"
        return "Other engineered"
    return "Unclassified"


def rank(lineage: object, r: str) -> str:
    m = re.search(rf"(?:^|;){r}__([^;]+)", str(lineage))
    return m.group(1).strip() if m else ""


def genomad_phylum(lineage: object) -> str:
    """Phylum (4th rank) of a geNomad 'Viruses;realm;kingdom;phylum;class;order;family' string; vMAG labels such as
    '..._caudoviricetes' map Caudoviricetes to Uroviricota."""
    s = str(lineage)
    p = s.split(";")
    if s.startswith("Viruses") and len(p) > 3 and p[3]:
        return p[3]
    if s.lower().endswith("caudoviricetes"):
        return "Uroviricota"
    return "Unassigned"


def card_enzyme(card_family: str) -> str:
    """Short enzyme label from a CARD gene family ('OXA-23-like beta-lactamase' -> 'OXA')."""
    name = re.sub(r"\s*beta-lactamase.*$", "", str(card_family), flags=re.I).strip()
    return "OXA" if name.upper().startswith("OXA") else (name or "No CARD hit")


def host_phylum(lineage: object) -> str:
    p = rank(lineage, "p")
    if not p:
        return "Unknown"
    if p.startswith("Firmicutes"):
        return "Bacillota"
    return {"Proteobacteria": "Pseudomonadota", "Cyanobacteria": "Cyanobacteriota"}.get(p, p)


def bh(p) -> np.ndarray:
    p = np.asarray(p, float)
    if len(p) == 0:
        return p
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def metagenome_id(uvig: str) -> str:
    """IMG taxon OID of the metagenome a UViG was assembled from (IMGVR_UViG_3300012345_000001,
    IMGVR_UViG_GVMAG-S-3300012345-93_000001); the UViG itself when no OID is found."""
    m = re.search(r"(\d{10})", str(uvig))
    return m.group(1) if m else str(uvig)


def fisher_block(fgt: pd.DataFrame, reft: pd.DataFrame, base: dict, ref_is_background: bool) -> list[dict]:
    """2x2 Fisher tests for every viral phylum x habitat of the foreground counts fgt against the reference counts
    reft (phylum x habitat). With ref_is_background the reference contains the foreground (IMG/VR background, so
    c = reference in habitat - a); otherwise the reference is a disjoint comparison set. BH within the block."""
    block = []
    for ph in fgt.index:
        n_fg = int(fgt.loc[ph].sum())
        n_ref = int(reft.loc[ph].sum()) if ph in reft.index else 0
        n_all = n_ref if ref_is_background else n_ref + n_fg
        for h in HABITATS:
            a = int(fgt.at[ph, h])
            ref_h = int(reft.at[ph, h]) if ph in reft.index else 0
            tot_h = ref_h if ref_is_background else ref_h + a
            b, c = n_fg - a, max(tot_h - a, 0)
            d = max(n_all - n_fg - c, 0)
            exp_a = n_fg * tot_h / n_all if n_all else math.nan
            log2e = math.log2((a / n_fg) / (tot_h / n_all)) if a and tot_h else math.nan
            odds, p = fisher_exact([[a, b], [c, d]]) if n_all > n_fg else (math.nan, math.nan)
            block.append({**base, "viral_phylum": ph, "habitat": h, "shown_in_heatmap": h in SHOWN_HABITATS,
                          "sbl_uvig_in_habitat": a, "sbl_uvig_phylum_total": n_fg,
                          "background_uvig_in_habitat": tot_h, "background_uvig_phylum_total": n_all,
                          "sbl_share": round(a / n_fg, 5) if n_fg else "",
                          "background_share": round(tot_h / n_all, 5) if n_all else "",
                          "expected_sbl_uvig_in_habitat": round(exp_a, 3),
                          "log2_enrichment": round(log2e, 4), "odds_ratio": odds, "fisher_p": p})
    ok = [i for i, r in enumerate(block) if not math.isnan(r["fisher_p"])]
    q = bh([block[i]["fisher_p"] for i in ok])
    for r in block:
        r["BH_q"], r["signif"] = math.nan, "ns"
    for i, qq in zip(ok, q):
        block[i]["BH_q"] = qq
        block[i]["signif"] = "***" if qq < 0.001 else "**" if qq < 0.01 else "*" if qq < 0.05 else "ns"
    return block


def counts(df: pd.DataFrame) -> pd.DataFrame:
    return df.groupby(["viral_phylum", "habitat"]).size().unstack(fill_value=0).reindex(columns=HABITATS, fill_value=0)


def main() -> None:
    cat = pd.read_csv(VIRAL_CAT, sep="\t", low_memory=False)
    cur = cat[(cat["pass_cov65"] == True) & (cat["pass_fimo_1e5"] == True)].drop_duplicates("uid").copy()  # noqa: E712
    cur["uvig"] = cur["uid"].astype(str).str.split("|").str[0]
    fam_by_uvig = cur.groupby("uvig")["winning_acc"].agg(lambda s: set(s))
    src = cur.groupby("uvig")["source"].first()

    info = pd.read_excel(INFO_XLSX, sheet_name="ALL_information", dtype=str)
    info = info.drop_duplicates("UVIG").set_index("UVIG")
    u = pd.DataFrame({"uvig": fam_by_uvig.index, "families": fam_by_uvig.values, "source": src.reindex(fam_by_uvig.index).values})
    u["in_imgvr_record"] = u["uvig"].isin(info.index)
    rec = info.reindex(u["uvig"])
    u["confidence"] = rec["Confidence"].fillna("").values
    u["habitat"] = [habitat_category(x) if ok else "Unclassified" for x, ok in zip(rec["Ecosystem classification"].values, u["in_imgvr_record"])]
    u["viral_phylum"] = [rank(x, "p") or "Unassigned" for x in rec["Taxonomic classification"].values]
    u["host_phylum"] = [host_phylum(x) for x in rec["Host taxonomy prediction"].values]
    u["host_method"] = rec["Host prediction method"].fillna("").values
    # MOSAiC vMAGs have no IMG/VR record: habitat from the GOLD biosample ecosystem of their metagenome,
    # viral phylum from the geNomad taxonomy of the genomic-context table
    mos = u["source"] == "MOSAiC"
    if mos.any():
        oid = u.loc[mos, "uvig"].str.split("_").str[1]
        gold = pd.read_csv(GOLD_CACHE, sep="\t", dtype=str).drop_duplicates("taxon_oid").set_index("taxon_oid")
        u.loc[mos, "habitat"] = oid.map(gold["gold_ecosystem_path"]).map(habitat_category).values
        nb = pd.read_csv(NEIGHBORHOODS, sep="\t", usecols=["uid", "genomad_taxonomy"], dtype=str)
        nb["uvig"] = nb["uid"].str.split("|").str[0]
        tax = nb.drop_duplicates("uvig").set_index("uvig")["genomad_taxonomy"]
        u.loc[mos, "viral_phylum"] = [genomad_phylum(tax.get(x, "")) for x in u.loc[mos, "uvig"]]
    for f in FAMILIES:
        u[f] = u["families"].map(lambda s, f=f: f in s)
    u["carries_main_family"] = u[MAIN_FAMILIES].any(axis=1)
    main_rows = cur[cur["winning_acc"].isin(MAIN_FAMILIES)]
    u["hydrolytic"] = u["uvig"].isin(main_rows.loc[main_rows["enzymatic_validation_class"] == HYDROLYTIC, "uvig"])

    # ---- Fig 1D flows: every curated SBL-carrying UViG (all confidence levels); the figure keeps carries_main_family
    flows = (u.assign(n_families=u["families"].map(len))
               .groupby(["habitat", "viral_phylum", "host_phylum", "confidence", "carries_main_family"], dropna=False)
               .agg(n_uvig=("uvig", "size"), **{f"n_{f}": (f, "sum") for f in FAMILIES})
               .reset_index().sort_values("n_uvig", ascending=False))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    flows.to_csv(OUT_DIR / "Habitat_ViralPhylum_Host_Flows_2026.tsv", sep="\t", index=False)

    # ---- Fig 1C enrichment: high-confidence UViGs against the IMG/VR high-confidence background
    bg = pd.read_csv(BG_COUNTS, sep="\t", dtype={"ecosystem_l1_l3": str})
    bg["habitat"] = bg["ecosystem_l1_l3"].fillna("").map(habitat_category)
    bgt = bg.groupby(["viral_phylum", "habitat"])["n_uvig"].sum().unstack(fill_value=0).reindex(columns=HABITATS, fill_value=0)
    hc = u[u["confidence"] == "High-confidence"].copy()
    rows = []
    for group in ["All_main", "Hydrolytic_main"] + FAMILIES:
        fg = (hc[hc["carries_main_family"]] if group == "All_main" else hc[hc["hydrolytic"]] if group == "Hydrolytic_main"
              else hc[hc[group]])
        rows += fisher_block(counts(fg), bgt, {"group": group, "main_text": group != "PF00768", "unit": "UViG"}, True)
    enr = pd.DataFrame(rows)
    enr.to_csv(OUT_DIR / "Habitat_Enrichment_By_Viral_Phylum_2026.tsv", sep="\t", index=False)

    # ---- Fig 1B: habitat enrichment per protein-architecture state (main families). A UViG belongs to a state when
    # it carries at least one curated main-family SBL of that state. Three analyses:
    #   vs_background   state UViGs against the IMG/VR high-confidence background of the same phylum (UViG level)
    #   vs_background, unit = metagenome   the same with at most one UViG per metagenome and state (one study or
    #                   sample cannot drive a result through many co-assembled UViGs)
    #   vs_other_SBL    state UViGs against the other main-family SBL-carrying UViGs of the phylum (does the
    #                   architecture prefer the habitat beyond what SBL carriage already shows?)
    arch = pd.read_csv(ARCH_CLASSES, sep="\t")
    arch = arch[arch["family"].isin(MAIN_FAMILIES)].assign(
        state=lambda d: d["architecture_class"].map(ARCH_STATE),
        uvig=lambda d: d["uid"].astype(str).str.split("|").str[0])
    states_by_uvig = arch.groupby("uvig")["state"].agg(set)
    hc["states"] = hc["uvig"].map(states_by_uvig).map(lambda s: s if isinstance(s, set) else set())
    hc["metagenome"] = hc["uvig"].map(metagenome_id)
    main_hc = hc[hc["carries_main_family"]]
    srows = []
    for st in ["All_SBL"] + STATES:  # All_SBL: every main-family SBL-carrying UViG (reference row of the heat map)
        has = main_hc["states"].map(lambda s, st=st: st == "All_SBL" or st in s)
        fg, other = main_hc[has], main_hc[~has]
        base = {"group": st, "main_text": True}
        srows += fisher_block(counts(fg), bgt, {**base, "comparison": "vs_background", "unit": "UViG"}, True)
        srows += fisher_block(counts(fg.drop_duplicates("metagenome")), bgt,
                              {**base, "comparison": "vs_background", "unit": "metagenome"}, True)
        if len(other):
            srows += fisher_block(counts(fg), counts(other), {**base, "comparison": "vs_other_SBL", "unit": "UViG"},
                                  False)
    senr = pd.DataFrame(srows)
    senr.to_csv(OUT_DIR / "Habitat_Enrichment_By_Architecture_State_2026.tsv", sep="\t", index=False)

    # ---- Fig 1C composition: functional class of the curated rows and host phylum of the UViGs, per family
    comp = []
    for fam in MAIN_FAMILIES + ["PF00768"]:
        fr = cur[cur["winning_acc"] == fam]
        fc = fr["enzymatic_validation_class"].map(lambda c: FUNC_GROUP.get(c, "Other"))
        for cat_, n in fc.value_counts().items():
            comp.append({"group": fam, "panel": "function", "category": cat_, "n": int(n), "total": len(fr)})
    for grp in MAIN_FAMILIES + ["PF00768", "Hydrolytic_main"]:
        fu = u[u["hydrolytic"]] if grp == "Hydrolytic_main" else u[u[grp]]
        hp = fu["host_phylum"].where(fu["host_phylum"].isin(MAJOR_HOSTS + ["Unknown"]), "Others")
        for cat_, n in hp.value_counts().items():
            comp.append({"group": grp, "panel": "host", "category": cat_, "n": int(n), "total": len(fu)})
    pd.DataFrame(comp).assign(main_text=lambda d: d["group"] != "PF00768").to_csv(
        OUT_DIR / "Family_Function_Host_Composition_2026.tsv", sep="\t", index=False)

    # ---- Fig 1C profile: one row per verified hydrolytic SBL of the main families
    hyd = main_rows[main_rows["enzymatic_validation_class"] == HYDROLYTIC].copy()
    ur = u.set_index("uvig")
    hrec = info.reindex(hyd["uvig"])
    hyd_prof = pd.DataFrame({
        "uid": hyd["uid"].values, "family": hyd["winning_acc"].values,
        "card_gene_family": hyd["card_gene_family"].fillna("").values,
        "card_enzyme": [card_enzyme(x) for x in hyd["card_gene_family"].fillna("")],
        "card_pident": hyd["card_pident"].values,
        "viral_origin_tier": hyd["viral_origin_tier"].values,
        "imgvr_confidence": ur.loc[hyd["uvig"], "confidence"].values,
        "viral_phylum": ur.loc[hyd["uvig"], "viral_phylum"].values,
        "host_phylum": ur.loc[hyd["uvig"], "host_phylum"].values,
        "host_class": [rank(x, "c") or "Unknown" for x in hrec["Host taxonomy prediction"].values],
        "host_method": ur.loc[hyd["uvig"], "host_method"].values,
        "habitat": ur.loc[hyd["uvig"], "habitat"].values})
    hyd_prof.to_csv(OUT_DIR / "Hydrolytic_SBL_Profile_2026.tsv", sep="\t", index=False)

    print(f"curated SBL-carrying UViGs {len(u):,} (IMG/VR record {int(u['in_imgvr_record'].sum()):,}; "
          f"high-confidence {len(hc):,}); background high-confidence UViGs {int(bgt.values.sum()):,}")
    print("habitat x viral phylum (all curated UViGs):")
    print(pd.crosstab(u["habitat"], u["viral_phylum"], margins=True).to_string())
    print("\nhost phylum (all curated UViGs):", u["host_phylum"].value_counts().to_dict())
    show = enr[(enr["group"] == "All_main") & enr["shown_in_heatmap"] & (enr["sbl_uvig_phylum_total"] >= 5)]
    print("\nAll_main enrichment, high-confidence (log2 [a/n]; signif):")
    piv = show.assign(cell=show["log2_enrichment"].round(1).astype(str) + " [" + show["sbl_uvig_in_habitat"].astype(str)
                      + "/" + show["expected_sbl_uvig_in_habitat"].round(1).astype(str) + "] " + show["signif"]
                      ).pivot(index="habitat", columns="viral_phylum", values="cell")
    print(piv.reindex(SHOWN_HABITATS).to_string())
    print("\nArchitecture states, Uroviricota (log2 [a/expected] signif):")
    for (cmp_, unit), d in senr[(senr["viral_phylum"] == "Uroviricota") & senr["shown_in_heatmap"]].groupby(
            ["comparison", "unit"], sort=False):
        print(f"-- {cmp_} / {unit}")
        print(d.assign(cell=d["log2_enrichment"].round(1).astype(str) + " [" + d["sbl_uvig_in_habitat"].astype(str)
                       + "/" + d["expected_sbl_uvig_in_habitat"].round(1).astype(str) + "] " + d["signif"])
              .pivot(index="habitat", columns="group", values="cell").reindex(SHOWN_HABITATS)[STATES].to_string())


if __name__ == "__main__":
    main()
