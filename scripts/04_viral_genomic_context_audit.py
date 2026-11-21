"""Module 02: Locus-Level Viral Origin Verification, Flanking Gene Neighborhoods, and Prophage Boundary Audit.

Addresses Reviewer Comment 9 (following Martin et al., Nature Microbiology 2025) and Comment 1 (part 1):
1. Documents for every viral SBL locus (n=5,701 across 5,618 contigs) its coordinates, CheckV/geNomad quality metrics,
   geNomad viral hallmark count, and exact distance/relative position within contig/prophage boundaries.
2. Reconstructs full upstream and downstream ORF genomic neighborhoods across all annotated high-quality viral contigs
   (merging HGT and non-HGT ORF/Pfam datasets), recording neighboring gene identities, Pfam domains, functional categories
   (Viral Hallmark/Structural, Viral Replication/Lysis, Recombination/Lysogeny, Auxiliary), relative ORF positions,
   and physical bp distances from each SBL gene.
3. Stratifies all 5,701 loci into explicit Viral Origin Confidence Tiers and generates a multi-panel verification figure.
"""

from pathlib import Path
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sbl_config as cfg

TABLES_DIR = str(cfg.TABLES_DIR)
SUPP_DIR = str(cfg.SUPP_DIR)
FIG_DIR = str(cfg.FIG_DIR)

M01_RAW = str(cfg.SEQ_DIR)
META_DIR = cfg.META_DIR
INFO_XLSX = str(META_DIR / "information.xlsx")
PFAM_CLANS = str(META_DIR / "Pfam-A.clans.tsv")
SUM_HIGH_FAA = str(META_DIR / "sum_high.faa")
SUM_HIGH_PFAM = str(META_DIR / "sum_high.pfam")
NO_HGT_FAA = str(META_DIR / "no_HGT.faa")
NO_HGT_PFAM = str(META_DIR / "no_HGT.pfam")
INPHAGE_XLSX = str(META_DIR / "Inphage_gggenome.xlsx")
PROPHAGE_BLAST = str(META_DIR / "Prophage_identify.blastnout")

FAMILIES = ["PF00144", "PF00905", "PF13354"]

RE_HALLMARK = re.compile(
    r"(capsid|portal|terminase|tail|baseplate|virion|sheath|tube|neck|head|spike|tape_?measure|"
    r"scaffolding|prohead|coat|spike|wedge|collar|connector|fiber|fibritin|bacteriophage_head|"
    r"phage_tail|phage_capsid|gp[0-9]|major_capsid|minor_capsid|packaging)",
    re.IGNORECASE,
)
RE_REPLIC_LYSIS = re.compile(
    r"(holin|endolysin|lysozyme|amidase|spanin|muramidase|peptidoglycan|phage_lysin|"
    r"dna_pol|helicase|primase|replicat|ssb|antirepressor|repressor|cro|ci_|cII|"
    r"transcription_factor|antiterminat|excisionase|xis|nuclease|exonuclease|endonuclease|"
    r"methylase|methyltransferase|thymidylate|ribonucleotide_reductase|dutpase)",
    re.IGNORECASE,
)
RE_RECOMB = re.compile(
    r"(integrase|resolvase|recombinase|transposase|insertion_element|rve\b|phage_integrase|"
    r"phage_int|pin\b|hin\b|gin\b|xer[cd]|rec[ao]|rad5[12]|orfB_IS)",
    re.IGNORECASE,
)


def classify_pfam(pfam_id, short_name, desc):
    base_pf = pfam_id.split(".")[0]
    if base_pf in FAMILIES:
        return "SBL_Target"
    if base_pf in ("PF00589", "PF00239", "PF02899", "PF09003", "PF13102", "PF13495", "PF00665", "PF01527", "PF07508"):
        return "Recombination_Lysogeny"
    text = f"{short_name} {desc}"
    if RE_HALLMARK.search(text):
        return "Viral_Hallmark_Structural"
    if RE_RECOMB.search(text):
        return "Recombination_Lysogeny"
    if RE_REPLIC_LYSIS.search(text):
        return "Viral_Replication_Lysis"
    return "Other_Auxiliary"


def load_pfam_clans():
    clans = {}
    with open(PFAM_CLANS, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                clans[parts[0]] = (parts[3], parts[4])
    return clans


def parse_sbl_loci():
    """Parse all 5,701 final viral SBL loci with Prodigal ORF coordinates."""
    loci = []
    hdr_re = re.compile(r"^(\S+)\s*#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?\d+)\s*#\s*ID=([^;\s]+)")
    for fam in FAMILIES:
        path = os.path.join(M01_RAW, f"{fam}_Viral_Final.faa")
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                hdr = line[1:].strip()
                m = hdr_re.match(hdr)
                if m:
                    sid, start, end, strand, prod_id = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)), m.group(5)
                else:
                    sid = hdr.split()[0]
                    start, end, strand, prod_id = np.nan, np.nan, np.nan, ""
                uvig = sid.split("|")[0]
                # Extract trailing ORF number from scaffold name if available (e.g., ..._35)
                orf_num = np.nan
                tail_part = sid.split("|")[-1]
                m_orf = re.search(r"_(\d+)$", tail_part)
                if m_orf:
                    orf_num = int(m_orf.group(1))
                loci.append({
                    "sbl_protein_id": sid,
                    "uvig": uvig,
                    "sbl_family": fam,
                    "sbl_start": start,
                    "sbl_end": end,
                    "sbl_strand": strand,
                    "prodigal_id": prod_id,
                    "sbl_orf_num": orf_num,
                })
    return pd.DataFrame(loci)


def parse_orf_faa(path, orfs_by_uvig):
    """Parse Prodigal .faa headers to get all ORFs and their coordinates on each contig."""
    hdr_re = re.compile(r"^(\S+)\s*#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?\d+)\s*#")
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            m = hdr_re.match(line[1:].strip())
            if not m:
                continue
            feat_id, start, end, strand = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
            uvig = feat_id.split("|")[0]
            if "|" in feat_id and feat_id.split("|")[-1].isdigit():
                orf_idx = int(feat_id.split("|")[-1])
            else:
                orf_idx = None
            orfs_by_uvig.setdefault(uvig, {})[feat_id] = {
                "feat_id": feat_id,
                "uvig": uvig,
                "start": start,
                "end": end,
                "strand": strand,
                "orf_idx": orf_idx,
                "pfams": [],
            }


def parse_pfam_hits(path, orfs_by_uvig, clans):
    """Attach Pfam-A domain annotations to ORFs."""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            feat_id = parts[0]
            uvig = feat_id.split("|")[0]
            # Handle different pfam_scan column layouts
            pfam_acc = ""
            pfam_name = ""
            for i, tok in enumerate(parts):
                if tok.startswith("PF") and re.match(r"^PF\d{5}(\.\d+)?$", tok):
                    pfam_acc = tok.split(".")[0]
                    if i + 1 < len(parts):
                        pfam_name = parts[i + 1]
                    break
            if not pfam_acc:
                continue
            short_n, desc = clans.get(pfam_acc, (pfam_name, pfam_name))
            cat = classify_pfam(pfam_acc, short_n or pfam_name, desc)
            if uvig in orfs_by_uvig and feat_id in orfs_by_uvig[uvig]:
                orfs_by_uvig[uvig][feat_id]["pfams"].append({
                    "pfam": pfam_acc,
                    "name": short_n or pfam_name,
                    "desc": desc,
                    "category": cat,
                })


def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(SUPP_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    print("Loading Pfam clans and 5,701 SBL loci...")
    clans = load_pfam_clans()
    sbl_df = parse_sbl_loci()

    # Load ALL_information (n=5,618)
    info = pd.read_excel(INFO_XLSX, sheet_name="ALL_information")
    info["UVIG"] = info["UVIG"].astype(str)
    # Parse Gene content (total genes;cds;tRNA;geNomad marker)
    gc_split = info["Gene content (total genes;cds;tRNA;geNomad marker)"].astype(str).str.split(";", expand=True)
    info["total_genes"] = pd.to_numeric(gc_split[0], errors="coerce")
    info["total_cds"] = pd.to_numeric(gc_split[1], errors="coerce")
    info["total_trna"] = pd.to_numeric(gc_split[2], errors="coerce")
    info["genomad_marker_genes"] = pd.to_numeric(gc_split[3], errors="coerce")
    info = info.rename(columns={
        "Coordinates ('whole' if the UViG is the entire contig)": "uvig_coords",
        "Length": "contig_length",
        "geNomad score": "genomad_score",
        "Estimated completeness": "checkv_completeness",
        "Estimated contamination": "checkv_contamination",
        "MIUViG quality": "miuvig_quality",
        "Taxonomic classification": "viral_taxonomy",
        "Host taxonomy prediction": "host_taxonomy",
        "Host prediction method": "host_method",
        "Ecosystem classification": "ecosystem",
    })

    # Load Prophage RefSeq matches
    prophage_hits = {}
    if os.path.exists(PROPHAGE_BLAST):
        with open(PROPHAGE_BLAST, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 6:
                    prophage_hits[p[0]] = f"{p[1]}:{p[4]}-{p[5]}"

    inph_bac_map = {}
    if os.path.exists(INPHAGE_XLSX):
        inph_bac = pd.read_excel(INPHAGE_XLSX, sheet_name="Bac")
        for _, r in inph_bac.iterrows():
            uv = str(r.get("Unnamed: 6", ""))
            if uv.startswith("IMGVR"):
                inph_bac_map[uv] = f"{r['seq_id']}:{r['start']}-{r['end']}(len={r['length']})"

    # Load ORF neighborhoods from sum_high + no_HGT
    print("Parsing ORF coordinates and Pfam annotations across high-quality neighborhoods...")
    orfs_by_uvig = {}
    parse_orf_faa(SUM_HIGH_FAA, orfs_by_uvig)
    parse_orf_faa(NO_HGT_FAA, orfs_by_uvig)
    parse_pfam_hits(SUM_HIGH_PFAM, orfs_by_uvig, clans)
    parse_pfam_hits(NO_HGT_PFAM, orfs_by_uvig, clans)
    print(f"Loaded multi-ORF Pfam neighborhoods for {len(orfs_by_uvig)} unique viral contigs.")

    # Merge SBL loci with contig metadata
    merged = sbl_df.merge(info, left_on="uvig", right_on="UVIG", how="left")

    # Compute prophage boundary metrics and flanking neighbor summaries per SBL locus
    records = []
    for _, row in merged.iterrows():
        uvig = row["uvig"]
        s_start = row["sbl_start"]
        s_end = row["sbl_end"]
        c_len = row["contig_length"]
        coords = str(row["uvig_coords"])

        # Prophage / contig boundary calculation
        is_prophage_sub = (coords != "whole") and ("-" in coords)
        if pd.notna(s_start) and pd.notna(s_end) and pd.notna(c_len) and c_len > 0:
            s_mid = 0.5 * (s_start + s_end)
            dist_5p = max(0.0, s_start - 1)
            dist_3p = max(0.0, c_len - s_end)
            min_dist_edge = min(dist_5p, dist_3p)
            rel_pos = round(s_mid / c_len, 4)
        else:
            dist_5p, dist_3p, min_dist_edge, rel_pos = np.nan, np.nan, np.nan, np.nan

        # Neighborhood analysis if contig is in orfs_by_uvig
        orf_dict = orfs_by_uvig.get(uvig, {})
        has_orf_map = len(orf_dict) > 0

        up_ann_str = ""
        dn_ann_str = ""
        up_viral_str = ""
        dn_viral_str = ""
        up5_str = ""
        dn5_str = ""
        n_hallmark = 0
        n_replic = 0
        n_recomb = 0
        has_up_viral_or_recomb = False
        has_dn_viral_or_recomb = False
        nearest_hallmark_dist_bp = np.nan
        nearest_recomb_dist_bp = np.nan

        if has_orf_map:
            orf_list = sorted(orf_dict.values(), key=lambda x: (x["start"], x["end"]))
            # Match SBL ORF by overlapping coordinates
            sbl_idx = None
            if pd.notna(s_start) and pd.notna(s_end):
                best_ov = -1
                for idx, o in enumerate(orf_list):
                    ov = min(s_end, o["end"]) - max(s_start, o["start"])
                    if ov > best_ov:
                        best_ov = ov
                        sbl_idx = idx
            if sbl_idx is None:
                # fallback: look for SBL Pfam in orf_list
                for idx, o in enumerate(orf_list):
                    if any(p["category"] == "SBL_Target" for p in o["pfams"]):
                        sbl_idx = idx
                        break

            if sbl_idx is not None:
                sbl_orf = orf_list[sbl_idx]
                s_mid_orf = 0.5 * (sbl_orf["start"] + sbl_orf["end"])

                def summarize_orf(idx):
                    o = orf_list[idx]
                    d_orf = idx - sbl_idx
                    # gap distance in bp (0 if adjacent/overlapping)
                    if idx < sbl_idx:
                        d_bp = max(0, sbl_orf["start"] - o["end"])
                    elif idx > sbl_idx:
                        d_bp = max(0, o["start"] - sbl_orf["end"])
                    else:
                        d_bp = 0
                    if o["pfams"]:
                        p0 = o["pfams"][0]
                        pf_label = "+".join(p["pfam"] + ":" + p["name"] for p in o["pfams"][:2])
                        cats = [p["category"] for p in o["pfams"]]
                        # prioritize viral hallmark > recomb > replic > aux
                        if "Viral_Hallmark_Structural" in cats:
                            main_cat = "Viral_Hallmark_Structural"
                        elif "Recombination_Lysogeny" in cats:
                            main_cat = "Recombination_Lysogeny"
                        elif "Viral_Replication_Lysis" in cats:
                            main_cat = "Viral_Replication_Lysis"
                        elif "SBL_Target" in cats:
                            main_cat = "SBL_Target"
                        else:
                            main_cat = "Other_Auxiliary"
                    else:
                        pf_label = "No_Pfam_Hit"
                        main_cat = "Unannotated_ORF"
                    return d_orf, d_bp, pf_label, main_cat, o["feat_id"]

                up_indices = list(range(sbl_idx - 1, -1, -1))
                dn_indices = list(range(sbl_idx + 1, len(orf_list)))

                up_summaries = [summarize_orf(i) for i in up_indices]
                dn_summaries = [summarize_orf(i) for i in dn_indices]

                # Count categories on contig
                for d_orf, d_bp, pfl, cat, fid in up_summaries + dn_summaries:
                    if cat == "Viral_Hallmark_Structural":
                        n_hallmark += 1
                        if pd.isna(nearest_hallmark_dist_bp) or d_bp < nearest_hallmark_dist_bp:
                            nearest_hallmark_dist_bp = d_bp
                    elif cat == "Viral_Replication_Lysis":
                        n_replic += 1
                    elif cat == "Recombination_Lysogeny":
                        n_recomb += 1
                        if pd.isna(nearest_recomb_dist_bp) or d_bp < nearest_recomb_dist_bp:
                            nearest_recomb_dist_bp = d_bp

                # Nearest annotated upstream & downstream
                up_ann = [x for x in up_summaries if x[3] != "Unannotated_ORF"]
                dn_ann = [x for x in dn_summaries if x[3] != "Unannotated_ORF"]
                if up_ann:
                    u = up_ann[0]
                    up_ann_str = f"{u[4]} (dORF={u[0]}, {u[1]}bp, {u[2]} [{u[3]}])"
                if dn_ann:
                    d = dn_ann[0]
                    dn_ann_str = f"{d[4]} (dORF={d[0]}, {d[1]}bp, {d[2]} [{d[3]}])"

                viral_cats = {"Viral_Hallmark_Structural", "Viral_Replication_Lysis", "Recombination_Lysogeny"}
                up_vir = [x for x in up_summaries if x[3] in viral_cats]
                dn_vir = [x for x in dn_summaries if x[3] in viral_cats]
                if up_vir:
                    u = up_vir[0]
                    up_viral_str = f"{u[4]} (dORF={u[0]}, {u[1]}bp, {u[2]} [{u[3]}])"
                    has_up_viral_or_recomb = True
                if dn_vir:
                    d = dn_vir[0]
                    dn_viral_str = f"{d[4]} (dORF={d[0]}, {d[1]}bp, {d[2]} [{d[3]}])"
                    has_dn_viral_or_recomb = True

                up5_str = "; ".join(f"{x[0]}:{x[1]}bp:{x[2]}[{x[3]}]" for x in up_summaries[:5])
                dn5_str = "; ".join(f"{x[0]}:{x[1]}bp:{x[2]}[{x[3]}]" for x in dn_summaries[:5])

        # Confidence tier classification
        g_score = row["genomad_score"] if pd.notna(row["genomad_score"]) else 0.0
        g_mark = row["genomad_marker_genes"] if pd.notna(row["genomad_marker_genes"]) else 0
        ckv_comp = row["checkv_completeness"] if pd.notna(row["checkv_completeness"]) else 0.0
        ckv_cont = row["checkv_contamination"] if pd.notna(row["checkv_contamination"]) else 100.0

        flanked_both = has_up_viral_or_recomb and has_dn_viral_or_recomb
        flanked_any = has_up_viral_or_recomb or has_dn_viral_or_recomb

        if (ckv_comp >= 90 or flanked_both or (g_mark >= 3 and flanked_any)) and ckv_cont <= 10.0:
            tier = "Tier1_HighConfidence_Viral"
        elif g_score >= 0.80 and g_mark >= 1 and ckv_cont <= 20.0:
            tier = "Tier2_Supported_Viral"
        else:
            tier = "Tier3_Fragment_or_Borderline"

        records.append({
            "sbl_protein_id": row["sbl_protein_id"],
            "uvig": uvig,
            "sbl_family": row["sbl_family"],
            "viral_origin_tier": tier,
            "miuvig_quality": row["miuvig_quality"],
            "contig_length_bp": c_len,
            "total_genes": row["total_genes"],
            "genomad_score": g_score,
            "genomad_marker_genes": g_mark,
            "checkv_completeness": row["checkv_completeness"],
            "checkv_contamination": row["checkv_contamination"],
            "uvig_coordinates": coords,
            "is_prophage_subregion": is_prophage_sub,
            "sbl_start_bp": s_start,
            "sbl_end_bp": s_end,
            "sbl_strand": row["sbl_strand"],
            "dist_to_5prime_boundary_bp": dist_5p,
            "dist_to_3prime_boundary_bp": dist_3p,
            "min_dist_to_boundary_bp": min_dist_edge,
            "relative_pos_in_uvig": rel_pos,
            "has_full_orf_pfam_map": has_orf_map,
            "n_viral_hallmark_on_contig": n_hallmark,
            "n_viral_replic_lysis_on_contig": n_replic,
            "n_recomb_lysogeny_on_contig": n_recomb,
            "nearest_hallmark_dist_bp": nearest_hallmark_dist_bp,
            "nearest_recomb_dist_bp": nearest_recomb_dist_bp,
            "flanked_both_sides_viral_or_recomb": flanked_both,
            "flanked_at_least_one_side_viral_or_recomb": flanked_any,
            "nearest_upstream_annotated_gene": up_ann_str,
            "nearest_downstream_annotated_gene": dn_ann_str,
            "nearest_upstream_viral_or_recomb_gene": up_viral_str,
            "nearest_downstream_viral_or_recomb_gene": dn_viral_str,
            "upstream_5_orfs": up5_str,
            "downstream_5_orfs": dn5_str,
            "refseq_prophage_hit": prophage_hits.get(uvig, inph_bac_map.get(uvig, "")),
            "viral_taxonomy": row["viral_taxonomy"],
            "host_taxonomy": row["host_taxonomy"],
            "host_method": row["host_method"],
            "ecosystem": row["ecosystem"],
        })

    out_df = pd.DataFrame(records)
    out_tsv = os.path.join(SUPP_DIR, "Supplementary_Data_Viral_Genomic_Neighborhoods.tsv")
    out_df.to_csv(out_tsv, sep="\t", index=False)

    # Summary statistics table
    tier_rows = []
    for tier_name, sub in [
        ("All_Viral_SBL_Loci", out_df),
        ("Tier1_HighConfidence_Viral", out_df[out_df["viral_origin_tier"] == "Tier1_HighConfidence_Viral"]),
        ("Tier2_Supported_Viral", out_df[out_df["viral_origin_tier"] == "Tier2_Supported_Viral"]),
        ("Tier3_Fragment_or_Borderline", out_df[out_df["viral_origin_tier"] == "Tier3_Fragment_or_Borderline"]),
        ("ORF_Neighborhood_Mapped_Subset", out_df[out_df["has_full_orf_pfam_map"]]),
        ("Prophage_Subregion_Subset", out_df[out_df["is_prophage_subregion"]]),
    ]:
        n_loci = len(sub)
        n_contigs = sub["uvig"].nunique()
        has_map = sub[sub["has_full_orf_pfam_map"]]
        tier_rows.append({
            "stratum": tier_name,
            "n_sbl_loci": n_loci,
            "pct_of_total_loci": round(100 * n_loci / len(out_df), 2),
            "n_unique_contigs": n_contigs,
            "median_genomad_score": round(float(sub["genomad_score"].median()), 4),
            "median_genomad_markers": round(float(sub["genomad_marker_genes"].median()), 1),
            "pct_with_ge1_genomad_marker": round(100 * float((sub["genomad_marker_genes"] >= 1).mean()), 2),
            "median_checkv_contamination": round(float(sub["checkv_contamination"].median()), 2),
            "median_min_dist_to_boundary_bp": round(float(sub["min_dist_to_boundary_bp"].median()), 1),
            "n_with_orf_neighborhood_map": len(has_map),
            "pct_mapped_with_viral_hallmark_on_contig": round(100 * float((has_map["n_viral_hallmark_on_contig"] >= 1).mean()), 2) if len(has_map) else np.nan,
            "pct_mapped_flanked_any_viral_recomb": round(100 * float(has_map["flanked_at_least_one_side_viral_or_recomb"].mean()), 2) if len(has_map) else np.nan,
            "pct_mapped_flanked_both_sides": round(100 * float(has_map["flanked_both_sides_viral_or_recomb"].mean()), 2) if len(has_map) else np.nan,
        })

    tier_df = pd.DataFrame(tier_rows)
    tier_df.to_csv(os.path.join(TABLES_DIR, "Table_S04_Viral_Origin_Confidence_Tiers.csv"), index=False)

    print("\n=== Viral Origin & Genomic Neighborhood Confidence Summary ===")
    print(tier_df.to_string(index=False))

    # Generate Figure: Viral Genomic Context Verification
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), dpi=300)

    # Panel A: geNomad viral hallmark marker genes per contig across tiers
    ax = axes[0]
    for t_name, col, lbl in [
        ("Tier1_HighConfidence_Viral", "#1B9E77", "Tier 1 (High-Conf)"),
        ("Tier2_Supported_Viral", "#386CB0", "Tier 2 (Supported)"),
        ("Tier3_Fragment_or_Borderline", "#E6AB02", "Tier 3 (Fragment)"),
    ]:
        vals = np.clip(out_df[out_df["viral_origin_tier"] == t_name]["genomad_marker_genes"].dropna().values, 0, 100)
        if len(vals) > 0:
            ax.hist(vals, bins=np.arange(0, 105, 5), alpha=0.65, color=col, label=f"{lbl} (n={len(vals)})", edgecolor="white")
    ax.set_xlabel("geNomad Viral Marker Genes per Contig (capped at 100)", fontweight="bold")
    ax.set_ylabel("Number of SBL Loci", fontweight="bold")
    ax.set_title("(A) Viral Hallmark Marker Density", fontweight="bold", loc="left")
    ax.legend(fontsize=8)

    # Panel B: Relative position of SBL inside prophage vs whole viral contigs
    ax = axes[1]
    pro_pos = out_df[out_df["is_prophage_subregion"]]["relative_pos_in_uvig"].dropna()
    whole_pos = out_df[~out_df["is_prophage_subregion"]]["relative_pos_in_uvig"].dropna()
    ax.hist(whole_pos, bins=20, density=True, alpha=0.55, color="#386CB0", label=f"Whole Viral Contigs (n={len(whole_pos)})", edgecolor="white")
    ax.hist(pro_pos, bins=20, density=True, alpha=0.65, color="#D95F02", label=f"Prophage Subregions (n={len(pro_pos)})", edgecolor="white")
    ax.axvline(0.5, color="#333333", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Normalized Position Inside Viral Contig / Prophage (0 = 5', 1 = 3')", fontweight="bold")
    ax.set_ylabel("Density", fontweight="bold")
    ax.set_title("(B) SBL Position Relative to Boundaries", fontweight="bold", loc="left")
    ax.legend(fontsize=8)

    # Panel C: Physical distance (kb) from SBL to nearest flanking viral hallmark vs recombination genes
    ax = axes[2]
    mapped = out_df[out_df["has_full_orf_pfam_map"]]
    d_hall = mapped["nearest_hallmark_dist_bp"].dropna() / 1000.0
    d_rec = mapped["nearest_recomb_dist_bp"].dropna() / 1000.0
    bp = ax.boxplot([d_hall, d_rec], tick_labels=[f"Viral Hallmark/\nStructural\n(n={len(d_hall)})", f"Integrase/\nRecombinase\n(n={len(d_rec)})"],
                    patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], ["#1B9E77", "#7570B3"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
    ax.set_ylabel("Genomic Distance from SBL Gene (kb)", fontweight="bold")
    ax.set_title("(C) Proximity to Flanking Viral Modules", fontweight="bold", loc="left")

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "Fig_S01_Viral_Genomic_Context_Verification.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(FIG_DIR, "Fig_S01_Viral_Genomic_Context_Verification.png"), bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
