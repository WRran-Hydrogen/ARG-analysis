#!/usr/bin/env python3
"""Audit of the CheckV host-gene check used by the viral-origin confidence tiers.

Read-only with respect to 01_build_neighborhoods.py outputs. For every primary IMG/VR locus it
  1. finds the CheckV gene (tmp/gene_features.tsv) that overlaps the focal SBL/PBP ORF and
     records its CheckV category (hmm_cat: 0 = no hit, 1 = viral, -1 = host/microbial);
  2. re-derives the current tier from the Phase 2 table (consistency check against the stored
     `viral_origin_tier`);
  3. evaluates two alternative host-gene rules:
       focal_corrected   host check `viral_genes >= host_genes - f`, where f = 1 only when the
                         focal ORF itself is a CheckV host gene inside the evaluated region; applied
                         to every locus, with no flank-based exemption;
       focal_corrected_hq  as above, plus whole-contig viruses that CheckV rates Complete or
                         High-quality with 0 % contamination (no host region detected) pass.

Outputs ($SBL_HPC_BASE/02_genomic_context/):
  Phase2_Host_Gene_Check_Audit.tsv
  Phase2_Host_Gene_Check_Audit_Summary.txt
"""

import csv
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

csv.field_size_limit(10_000_000)

OUT_DIR = hpc.BASE / "02_genomic_context"
PHASE2_TSV = OUT_DIR / "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"
CHECKV_DIR = hpc.BASE / "00_checkv" / "out"
GENE_FEATURES = CHECKV_DIR / "tmp" / "gene_features.tsv"
CONTAMINATION = CHECKV_DIR / "contamination.tsv"
BOUNDARY_MIN_BP = 1000
TIERS = ("Tier1_HighConfidence_Viral", "Tier2_Supported_Viral", "Tier3_Fragment_or_Borderline")


def overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def as_int(x, default=0):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return default


def as_float(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def truthy(x):
    return str(x).strip().lower() in ("true", "1", "yes")


def checkv_regions(row):
    """Viral regions of a CheckV provirus as (start, end, viral_genes, host_genes)."""
    types = (row.get("region_types") or "").split(",")
    coords = (row.get("region_coords_bp") or "").split(",")
    vgs = (row.get("region_viral_genes") or "").split(",")
    hgs = (row.get("region_host_genes") or "").split(",")
    out = []
    for i, rt in enumerate(types):
        if rt.strip().lower() != "viral" or i >= len(coords):
            continue
        m = re.match(r"^(\d+)-(\d+)$", coords[i].strip())
        if m:
            vg = int(vgs[i]) if i < len(vgs) and vgs[i].strip().isdigit() else 0
            hg = int(hgs[i]) if i < len(hgs) and hgs[i].strip().isdigit() else 0
            out.append((int(m.group(1)), int(m.group(2)), vg, hg))
    return out


def tier_from_flags(boundary_ok, n_hm, n_rep, gn_present, gn_score, gn_nhm, vg, t1_host_ok, t2_host_ok):
    if boundary_ok and n_hm >= 1 and t1_host_ok:
        return TIERS[0]
    if gn_present and gn_score >= 0.70 and t2_host_ok and (n_rep >= 1 or vg >= 1 or n_hm >= 1 or gn_nhm >= 1):
        return TIERS[1]
    return TIERS[2]


def main():
    rows = []
    with open(PHASE2_TSV, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["source"] != "MOSAiC" and truthy(r["in_primary_cov65_fimo1e5"]) and r["host_gene_check_basis"] != "NA":
                rows.append(r)
    uvigs = {r["uvig"] for r in rows}
    print(f"Primary IMG/VR loci with genomic context: {len(rows):,}")

    genes = {}
    with open(GENE_FEATURES, encoding="utf-8") as fh:
        for g in csv.DictReader(fh, delimiter="\t"):
            if g["contig_id"] in uvigs:
                genes.setdefault(g["contig_id"], []).append(
                    (int(g["start"]), int(g["end"]), as_int(g["hmm_cat"]), g["hmm_db"], g["hmm_name"]))

    regions = {}
    with open(CONTAMINATION, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["contig_id"] in uvigs and r.get("provirus", "No") == "Yes":
                regions[r["contig_id"]] = checkv_regions(r)

    out, mismatches = [], 0
    for r in rows:
        uvig, st, en = r["uvig"], as_int(r["sbl_start_bp"]), as_int(r["sbl_end_bp"])
        basis = r["host_gene_check_basis"]

        # focal ORF -> CheckV gene with the largest overlap (>= 50 % of the shorter gene)
        focal_cat, focal_db, focal_hmm = "unmapped", "NA", "NA"
        if st > 0:
            best, best_ov = None, 0
            for g in genes.get(uvig, []):
                ov = overlap(st, en, g[0], g[1])
                if ov > best_ov:
                    best, best_ov = g, ov
            if best is not None and best_ov >= 0.5 * min(en - st + 1, best[1] - best[0] + 1):
                focal_cat = {0: "no_hit", 1: "viral", -1: "host"}.get(best[2], str(best[2]))
                focal_db, focal_hmm = best[3], best[4]

        vg_c, hg_c = as_int(r["checkv_viral_genes"]), as_int(r["checkv_host_genes"])
        if basis == "checkv_viral_region":
            reg = [x for x in regions.get(uvig, []) if overlap(st, en, x[0], x[1]) > 0]
            _, _, vg, hg = max(reg, key=lambda x: overlap(st, en, x[0], x[1])) if reg else (0, 0, vg_c, hg_c)
        else:
            vg, hg = vg_c, hg_c

        topo = r["genomad_topology"]
        gn_present = topo != "Not_in_virus_summary"
        boundary_ok = gn_present and (topo in ("DTR", "ITR") or as_float(r["min_dist_to_viral_boundary_bp"]) >= BOUNDARY_MIN_BP)
        n_hm, n_rep = as_int(r["n_viral_hallmark_in_region"]), as_int(r["n_viral_replic_lysis_in_region"])
        gn_score, gn_nhm = as_float(r["genomad_virus_score"]), as_int(r["genomad_n_hallmarks"])
        flank_hm2 = truthy(r["flanked_strict_structural_hallmark_both"])
        flank_sr2 = truthy(r["flanked_viral_structural_or_replic_both"])
        raw_ok = vg >= hg
        skip = basis == "genomad_provirus_skip"

        # current rule (01_build_neighborhoods.py)
        if skip:
            t1, t2 = flank_sr2 and vg >= 1, True
        else:
            t1 = vg >= 1 and (raw_ok or (flank_hm2 and vg >= hg - 1))
            t2 = raw_ok or flank_hm2
        tier_cur = tier_from_flags(boundary_ok, n_hm, n_rep, gn_present, gn_score, gn_nhm, vg, t1, t2)
        if tier_cur != r["viral_origin_tier"]:
            mismatches += 1

        # alternative A: focal-gene correction for every locus, no flank exemption
        f = 1 if (focal_cat == "host" and not skip) else 0
        ok_a = vg >= hg - f
        t1a, t2a = (flank_sr2 and vg >= 1, True) if skip else (vg >= 1 and ok_a, ok_a)
        tier_a = tier_from_flags(boundary_ok, n_hm, n_rep, gn_present, gn_score, gn_nhm, vg, t1a, t2a)

        # alternative B: A + CheckV Complete/High-quality whole contigs with no host region
        hq_clean = (basis == "whole_contig" and r["checkv_quality"] in ("Complete", "High-quality")
                    and as_float(r["checkv_contamination_pct"]) == 0.0)
        ok_b = ok_a or hq_clean
        t1b, t2b = (flank_sr2 and vg >= 1, True) if skip else (vg >= 1 and ok_b, ok_b)
        tier_b = tier_from_flags(boundary_ok, n_hm, n_rep, gn_present, gn_score, gn_nhm, vg, t1b, t2b)

        out.append({
            "uid": r["uid"], "uvig": uvig, "winning_acc": r["winning_acc"], "host_gene_check_basis": basis,
            "eval_viral_genes": vg, "eval_host_genes": hg, "focal_checkv_category": focal_cat,
            "focal_checkv_hmm_db": focal_db, "focal_checkv_hmm": focal_hmm,
            "flanked_strict_structural_hallmark_both": flank_hm2, "checkv_quality": r["checkv_quality"],
            "checkv_contamination_pct": r["checkv_contamination_pct"],
            "tier_current": r["viral_origin_tier"], "tier_rederived": tier_cur,
            "tier_focal_corrected": tier_a, "tier_focal_corrected_hq": tier_b,
        })

    out_tsv = OUT_DIR / "Phase2_Host_Gene_Check_Audit.tsv"
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    def short(t):
        return t.split("_")[0]

    lines = ["=== CheckV host-gene check audit (primary IMG/VR loci) ===",
             f"Loci: {len(out):,} | re-derived tier differs from stored tier: {mismatches}", "",
             "Focal ORF CheckV category by host-check basis:"]
    for (b, c), n in sorted(Counter((x["host_gene_check_basis"], x["focal_checkv_category"]) for x in out).items()):
        lines.append(f"  {b:24s} {c:9s}: {n:,}")
    lines += ["", "Focal ORF CheckV category by family:"]
    for (fam, c), n in sorted(Counter((x["winning_acc"], x["focal_checkv_category"]) for x in out).items()):
        lines.append(f"  {fam:8s} {c:9s}: {n:,}")
    rescued = [x for x in out if x["eval_viral_genes"] < x["eval_host_genes"] and short(x["tier_current"]) in ("Tier1", "Tier2")
               and x["host_gene_check_basis"] != "genomad_provirus_skip"]
    lines += ["", f"Loci in Tier 1/2 that fail the raw check (flank exemption): {len(rescued)}",
              "  focal category: " + ", ".join(f"{k}={v}" for k, v in sorted(Counter(x['focal_checkv_category'] for x in rescued).items())),
              "  host - viral genes: " + ", ".join(f"{k}={v}" for k, v in sorted(Counter(min(x['eval_host_genes'] - x['eval_viral_genes'], 10) for x in rescued).items()))
              + " (10 = >=10)"]
    for col in ("tier_focal_corrected", "tier_focal_corrected_hq"):
        lines += ["", f"Current tier (rows) -> {col} (columns):", f"  {'':30s}" + "".join(f"{short(t):>9s}" for t in TIERS)]
        for a in TIERS:
            lines.append(f"  {a:30s}" + "".join(f"{sum(1 for x in out if x['tier_current'] == a and x[col] == b):>9,}" for b in TIERS))
        lines.append("  by family (T1/T2/T3): " + "; ".join(
            f"{fam} " + "/".join(str(sum(1 for x in out if x["winning_acc"] == fam and x[col] == t)) for t in TIERS)
            for fam in sorted({x["winning_acc"] for x in out})))
    rep = "\n".join(lines) + "\n"
    (OUT_DIR / "Phase2_Host_Gene_Check_Audit_Summary.txt").write_text(rep, encoding="utf-8")
    print(rep)


if __name__ == "__main__":
    sys.exit(main())
