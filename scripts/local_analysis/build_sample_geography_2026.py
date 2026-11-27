#!/usr/bin/env python3
"""Sample geography of IMG/VR v4.1 and climate-zone prevalence of viral SBL carriers (Fig. 1 inputs).

Every IMG/VR sample is an IMG taxon OID. GOLD links it to coordinates:
  Analysis Project (AP IMG TAXON ID) -> Sequencing Project (AP PROJECT GOLD IDS)
  -> Biosample (metagenomes) or Organism (isolates) latitude / longitude.
Biosample coordinates are used when present, otherwise organism coordinates. Coordinates outside the valid
range and the (0, 0) placeholder are discarded.

Climate zones by absolute latitude: Tropical <= 23.44, Temperate 23.44-66.56, Polar > 66.56.
A sample is positive for a Pfam family when it carries at least one curated viral SBL row of that family
(Final_Viral_SBL_Catalog.tsv, HMM coverage >= 65 % and FIMO motif pass). Prevalence is reported per sample
(taxon OID) and per site (unique coordinates rounded to 0.01 degree), with Wilson 95 % intervals and pairwise
Fisher exact tests between zones (Benjamini-Hochberg within each family, unit and scope). The pooled group
'any_main' counts samples carrying any of the three main-text families (PF00144, PF00905, PF13354); PF00768 rows
(main_text = False) are kept for the supplement.

Two sample universes (column 'scope' of the prevalence table):
  IMGVR_only          IMG/VR v4.1 samples (102,253 taxon OIDs)
  IMGVR_plus_MOSAiC   plus the 270 MOSAiC metagenomes whose viral proteomes were screened (central Arctic,
                      Polar zone); 8 of them are also IMG/VR samples and are counted once. MOSAiC viruses were
                      identified with a different pipeline (geNomad + vMAG binning) from IMG/VR, so this scope is
                      reported as a sensitivity analysis of the Polar estimate.

Inputs
  data/metadata/goldData.xlsx                                GOLD public export (not tracked; see .gitignore)
  data/metadata/GOLD_TaxonOID_Coordinates.tsv.gz             cache of load_gold(), rebuilt when goldData.xlsx is newer
  data/metadata/IMGVR_v4.1_HighConfidence_Taxon_Ecosystem.tsv.gz
      unique (Taxon_oid, Ecosystem classification) of IMGVR_all_Sequence_information-high_confidence.tsv
  data/metadata/MOSAiC_Viral_Metagenome_TaxonOIDs.tsv        MOSAiC metagenomes in mosaic_virus.faa (viral proteins each)
  results/reanalysis_2026/01_curation_and_context/Final_Viral_SBL_Catalog.tsv
Habitat: the IMG/VR ecosystem path, or the GOLD biosample ecosystem path for samples outside IMG/VR (MOSAiC),
collapsed with habitat_enrichment_2026.habitat_category().

Outputs (results/reanalysis_2026/06_structure_and_ecology/)
  Sample_Geography_2026.tsv
  Climate_Zone_Prevalence_2026.tsv
"""

from __future__ import annotations

import math
import subprocess
import sys
from itertools import combinations
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from habitat_enrichment_2026 import habitat_category  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
META = BASE / "data" / "metadata"
RES = BASE / "results" / "reanalysis_2026"
GOLD_XLSX = META / "goldData.xlsx"
GOLD_CACHE = META / "GOLD_TaxonOID_Coordinates.tsv.gz"
IMGVR_TAXA = META / "IMGVR_v4.1_HighConfidence_Taxon_Ecosystem.tsv.gz"
MOSAIC_TAXA = META / "MOSAiC_Viral_Metagenome_TaxonOIDs.tsv"
VIRAL_CAT = RES / "01_curation_and_context" / "Final_Viral_SBL_Catalog.tsv"
OUT_DIR = RES / "06_structure_and_ecology"
FAMILIES = ["PF00144", "PF00905", "PF13354", "PF00768"]
# Main-text families. PF00768 (homologous PBPs separated from PF13354 by the clan bitscore assignment) is kept in
# the tables for the supplement; the pooled group "any_main" pools the three main families only.
MAIN_FAMILIES = ["PF00144", "PF00905", "PF13354"]
ZONES = ["Tropical", "Temperate", "Polar"]
TROPIC, POLAR = 23.44, 66.56


def climate_zone(lat: float) -> str:
    a = abs(lat)
    return "Tropical" if a <= TROPIC else ("Temperate" if a <= POLAR else "Polar")


def to_float(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.strip().replace({"": np.nan, "None": np.nan}), errors="coerce")


def valid_coords(lat: pd.Series, lon: pd.Series) -> pd.Series:
    return lat.between(-90, 90) & lon.between(-180, 180) & ~((lat == 0) & (lon == 0))


def load_gold() -> pd.DataFrame:
    """One row per IMG taxon OID with GOLD coordinates, from the cache when it is up to date."""
    fresh = GOLD_CACHE.is_file() and (not GOLD_XLSX.is_file() or GOLD_CACHE.stat().st_mtime > GOLD_XLSX.stat().st_mtime)
    if fresh and "gold_ecosystem_path" in pd.read_csv(GOLD_CACHE, sep="\t", nrows=0).columns:
        return pd.read_csv(GOLD_CACHE, sep="\t", dtype=str).assign(
            latitude=lambda d: pd.to_numeric(d["latitude"]), longitude=lambda d: pd.to_numeric(d["longitude"]),
            coord_source=lambda d: d["coord_source"].fillna(""))
    df = parse_gold_xlsx()
    df.to_csv(GOLD_CACHE, sep="\t", index=False)
    return df


def parse_gold_xlsx() -> pd.DataFrame:
    """One row per IMG taxon OID with GOLD coordinates (biosample first, organism second)."""
    cols = {
        "Analysis Project": ["AP GOLD ID", "AP TYPE", "AP IMG TAXON ID", "AP PROJECT GOLD IDS"],
        "Sequencing Project": ["PROJECT GOLD ID", "BIOSAMPLE GOLD ID", "ORGANISM GOLD ID"],
        "Biosample": ["BIOSAMPLE GOLD ID", "BIOSAMPLE LATITUDE", "BIOSAMPLE LONGITUDE", "BIOSAMPLE GEOGRAPHIC LOCATION",
                      "BIOSAMPLE ECOSYSTEM", "BIOSAMPLE ECOSYSTEM CATEGORY", "BIOSAMPLE ECOSYSTEM TYPE"],
        "Organism": ["ORGANISM GOLD ID", "ORGANISM LATITUDE", "ORGANISM LONGITUDE", "ORGANISM GEOGRAPHIC LOCATION"],
    }
    sh = {k: pd.read_excel(GOLD_XLSX, sheet_name=k, usecols=v, dtype=str) for k, v in cols.items()}
    ap = sh["Analysis Project"].dropna(subset=["AP IMG TAXON ID"]).copy()
    ap["taxon_oid"] = ap["AP IMG TAXON ID"].str.strip().str.replace(r"\.0$", "", regex=True)
    ap["project"] = ap["AP PROJECT GOLD IDS"].fillna("").str.split(r"[|,;\s]+")
    ap = ap.explode("project")
    ap = ap[ap["project"].str.len() > 0]
    sp = sh["Sequencing Project"].rename(columns={"PROJECT GOLD ID": "project"})
    bs, og = sh["Biosample"], sh["Organism"]
    bs["b_lat"], bs["b_lon"] = to_float(bs["BIOSAMPLE LATITUDE"]), to_float(bs["BIOSAMPLE LONGITUDE"])
    bs["gold_ecosystem_path"] = (bs[["BIOSAMPLE ECOSYSTEM", "BIOSAMPLE ECOSYSTEM CATEGORY", "BIOSAMPLE ECOSYSTEM TYPE"]]
                                 .fillna("").agg(";".join, axis=1))
    og["o_lat"], og["o_lon"] = to_float(og["ORGANISM LATITUDE"]), to_float(og["ORGANISM LONGITUDE"])
    df = (ap.merge(sp, on="project", how="left")
            .merge(bs[["BIOSAMPLE GOLD ID", "b_lat", "b_lon", "BIOSAMPLE GEOGRAPHIC LOCATION", "gold_ecosystem_path"]], on="BIOSAMPLE GOLD ID", how="left")
            .merge(og[["ORGANISM GOLD ID", "o_lat", "o_lon", "ORGANISM GEOGRAPHIC LOCATION"]], on="ORGANISM GOLD ID", how="left"))
    b_ok, o_ok = valid_coords(df["b_lat"], df["b_lon"]), valid_coords(df["o_lat"], df["o_lon"])
    df["latitude"] = np.where(b_ok, df["b_lat"], np.where(o_ok, df["o_lat"], np.nan))
    df["longitude"] = np.where(b_ok, df["b_lon"], np.where(o_ok, df["o_lon"], np.nan))
    df["coord_source"] = np.where(b_ok, "GOLD_Biosample", np.where(o_ok, "GOLD_Organism", ""))
    df["geographic_location"] = df["BIOSAMPLE GEOGRAPHIC LOCATION"].where(b_ok, df["ORGANISM GEOGRAPHIC LOCATION"])
    df["has_coord"] = df["coord_source"] != ""
    df = df.sort_values(["taxon_oid", "has_coord"], ascending=[True, False]).drop_duplicates("taxon_oid")
    return df[["taxon_oid", "AP GOLD ID", "AP TYPE", "project", "BIOSAMPLE GOLD ID", "ORGANISM GOLD ID",
               "latitude", "longitude", "coord_source", "geographic_location", "gold_ecosystem_path"]].rename(
        columns={"AP GOLD ID": "ap_gold_id", "AP TYPE": "ap_type", "project": "project_gold_id",
                 "BIOSAMPLE GOLD ID": "biosample_gold_id", "ORGANISM GOLD ID": "organism_gold_id"})


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def bh(p: list[float]) -> list[float]:
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out.tolist()


def prevalence_rows(gz: pd.DataFrame, scope: str) -> list[dict]:
    rows = []
    for unit in ("sample", "site"):
        for fam in FAMILIES + ["any_main"]:
            col = f"n_{fam}"
            if unit == "sample":
                tot = gz.groupby("climate_zone").size()
                pos = gz[gz[col] > 0].groupby("climate_zone").size()
            else:
                tot = gz.groupby("climate_zone")["site_id"].nunique()
                pos = gz[gz[col] > 0].groupby("climate_zone")["site_id"].nunique()
            kn = {z: (int(pos.get(z, 0)), int(tot.get(z, 0))) for z in ZONES}
            pairs = list(combinations(ZONES, 2))
            ps = [fisher_exact([[kn[a][0], kn[a][1] - kn[a][0]], [kn[b][0], kn[b][1] - kn[b][0]]])[1] for a, b in pairs]
            qs = dict(zip(pairs, bh(ps)))
            for z in ZONES:
                k, n = kn[z]
                lo, hi = wilson(k, n)
                rows.append({"scope": scope, "unit": unit, "family": fam, "main_text": fam != "PF00768", "climate_zone": z, "positive": k, "total": n,
                             "prevalence": round(k / n, 5) if n else "", "ci95_low": round(lo, 5), "ci95_high": round(hi, 5),
                             **{f"BH_q_vs_{b}": round(qs[(z, b)] if (z, b) in qs else qs[(b, z)], 6) for b in ZONES if b != z}})
    return rows


def main() -> None:
    taxa = pd.read_csv(IMGVR_TAXA, sep="\t", header=None, names=["taxon_oid", "imgvr_ecosystem"], dtype=str)
    taxa = taxa.drop_duplicates("taxon_oid").assign(in_imgvr=True)
    mos = pd.read_csv(MOSAIC_TAXA, sep="\t", dtype={"taxon_oid": str})[["taxon_oid"]].assign(in_mosaic=True)
    samples = taxa.merge(mos, on="taxon_oid", how="outer")
    for c in ("in_imgvr", "in_mosaic"):
        samples[c] = samples[c].eq(True)
    geo = samples.merge(load_gold(), on="taxon_oid", how="left")
    geo["coord_source"] = geo["coord_source"].fillna("")
    eco = geo["imgvr_ecosystem"].where(geo["imgvr_ecosystem"].fillna("").str.strip(";").str.len() > 0,
                                       geo["gold_ecosystem_path"])
    geo["habitat"] = eco.map(habitat_category)
    has = geo["coord_source"] != ""
    geo["climate_zone"] = ""
    geo.loc[has, "climate_zone"] = geo.loc[has, "latitude"].map(climate_zone)
    geo["site_id"] = ""
    geo.loc[has, "site_id"] = (geo.loc[has, "latitude"].round(2).map("{:.2f}".format) + "_"
                               + geo.loc[has, "longitude"].round(2).map("{:.2f}".format))

    cat = pd.read_csv(VIRAL_CAT, sep="\t", low_memory=False, dtype={"uid": str})
    cur = cat[(cat["pass_cov65"] == True) & (cat["pass_fimo_1e5"] == True)].drop_duplicates("uid").copy()  # noqa: E712
    is_img = cur["source"] == "IMG_VR"
    cur["taxon_oid"] = np.where(is_img, cur["uid"].str.split("|").str[1], cur["uid"].str.split("_").str[1])
    for fam in FAMILIES:
        n = cur[cur["winning_acc"] == fam].groupby("taxon_oid").size()
        geo[f"n_{fam}"] = geo["taxon_oid"].map(n).fillna(0).astype(int)
    geo["n_any_main"] = geo[[f"n_{f}" for f in MAIN_FAMILIES]].sum(axis=1)
    hyd = cur[cur["winning_acc"].isin(MAIN_FAMILIES) & (cur["enzymatic_validation_class"] == "Verified_Hydrolytic_SBL")]
    geo["n_hydrolytic_main"] = geo["taxon_oid"].map(hyd.groupby("taxon_oid").size()).fillna(0).astype(int)
    for fam in MAIN_FAMILIES:
        geo[f"n_hyd_{fam}"] = geo["taxon_oid"].map(hyd[hyd["winning_acc"] == fam].groupby("taxon_oid").size()).fillna(0).astype(int)
    missing = set(cur["taxon_oid"]) - set(geo["taxon_oid"])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    geo.to_csv(OUT_DIR / "Sample_Geography_2026.tsv", sep="\t", index=False)

    rows = prevalence_rows(geo[has & geo["in_imgvr"]], "IMGVR_only") + prevalence_rows(geo[has], "IMGVR_plus_MOSAiC")
    pd.DataFrame(rows).to_csv(OUT_DIR / "Climate_Zone_Prevalence_2026.tsv", sep="\t", index=False)

    for lab, m in (("IMG/VR", geo["in_imgvr"]), ("MOSAiC", geo["in_mosaic"])):
        print(f"{lab} samples {int(m.sum()):,}; with coordinates {int((m & has).sum()):,}; "
              f"sites {geo.loc[m & has, 'site_id'].nunique():,}; SBL-positive (main families) {int((m & (geo['n_any_main'] > 0)).sum()):,} "
              f"(with coordinates {int((m & has & (geo['n_any_main'] > 0)).sum()):,})")
    print(f"IMG/VR and MOSAiC overlap: {int((geo['in_imgvr'] & geo['in_mosaic']).sum())}; "
          f"curated SBL taxon OIDs not in either sample list: {len(missing)}")
    out = pd.DataFrame(rows)
    print(out[out["climate_zone"] == "Polar"][["scope", "unit", "family", "positive", "total", "prevalence",
                                                 "BH_q_vs_Temperate", "BH_q_vs_Tropical"]].to_string(index=False))


if __name__ == "__main__":
    main()
