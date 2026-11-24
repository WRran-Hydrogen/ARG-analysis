#!/usr/bin/env python3
"""Phase 2: 100%-Coverage Viral Genomic Neighborhood & Strict Structural Confidence Tier Integration.

Integrates for all curated IMG/VR viral SBL/PRE loci:
- Prodigal meta ORFs (all_orfs.faa)
- Full-genome Pfam-A (--cut_ga) annotations (orfs_pfam.domtbl)
- Full-genome eggNOG-mapper annotations (sbl_neighbourhood.emapper.annotations.gz)
- geNomad v1.7 summary, hallmarks, and provirus boundaries
- CheckV v1.0.3 (DB v1.5) quality and viral/host region coordinates

Tier rules:
- Strictly separates viral structural hallmark genes (capsid, portal, terminase, tail, baseplate,
  virion coat/sheath/fiber/spike, prohead/scaffolding) from viral replication/lysis genes and
  from recombination/integrase/transposase genes.
- Tier1_HighConfidence_Viral: SBL >= 1 kb inside the viral region (or DTR/ITR
  complete genome) + >=1 ORF-verified strict viral structural hallmark inside the viral region
  (`n_viral_hallmark_in_region >= 1`; raw geNomad `gn_nhm` is restricted to Tier 2 and never
  counting integrase/recombinase toward Tier 1) + CheckV viral_genes >= 1 and host-gene check:
  `viral_genes >= host_genes` on the CheckV viral region containing the SBL (or whole contig for
  non-provirus contigs), or `viral_genes >= host_genes - 1` when flanked on BOTH sides by strict
  viral structural hallmarks (accounting for CheckV self-counting the focal SBL/PBP as 1 host gene);
  for `genomad_provirus_skip` intervals, Tier 1 requires double-sided viral structural or
  replication/lysis flanking (`flanked_viral_structural_or_replic_both`) plus `viral_genes >= 1`.
  `host_gene_check_pass` preserves the raw `vg_eval >= hg_eval` status for audit transparency.
- geNomad hallmark counts are taken from the geNomad entry (provirus region or whole contig) that
  contains the SBL ORF.
- Also reports, per locus, the tier under a lenient comparison rule (>= 2 hallmarks including
  holin/endolysin, 2 kb boundary) and a lenient -> strict tier shift table.

Run after 01_curation/09_build_curated_datasets.py --stage curate and before --stage finalize.
"""

import argparse
import csv
import gzip
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

BASE = hpc.BASE
P1_OUT = BASE / "01_final_curated"
P1_DEC = BASE / "01_curated_datasets" / "Phase1_Viral_Decontamination_Audit.tsv"
CHECKV_DIR = BASE / "00_checkv" / "out"
GENOMAD_DIR = hpc.GENOMAD_DIR
PRODIGAL_ORFS = hpc.PRODIGAL_DIR / "all_orfs.faa"
ORFS_PFAM = hpc.ORFS_PFAM_DIR / "orfs_pfam.domtbl"
EGGNOG_GZ = hpc.EGGNOG_DIR / "sbl_neighbourhood.emapper.annotations.gz"

OUT_DIR = BASE / "02_genomic_context"
BOUNDARY_MIN_BP = 1000   # minimum distance from the SBL ORF to the viral region boundary for Tier 1

RE_STRICT_STRUCTURAL = re.compile(
    r"(terminase|portal|capsid|head[- _]protein|tail|baseplate|tape[- _]?measure|virion|"
    r"phage[- _]?(?:coat|sheath|neck|fiber|spike|tube|wedge|hub|structural)|"
    r"prohead|scaffolding|DTR|MCP|major[- _]?capsid|minor[- _]?capsid|connector|collar|"
    r"gp\d+.*(?:tail|head|capsid|portal))",
    re.I,
)
RE_REPLIC_LYSIS = re.compile(
    r"(holin|endolysin|lysozyme|muramidase|amidase|peptidoglycan|phage[- _]?lysin|spanin|lysis|"
    r"dna[- _]?polymerase|rna[- _]?polymerase|helicase|primase|exonuclease|endonuclease|"
    r"dutpase|ribonucleotide[- _]?reductase|thymidylate|single[- _]?strand.*binding|"
    r"sigma[- _]?70|anti[- _]?repressor|replication|ssb|recT|recE|erf)",
    re.I,
)
# Lenient-rule hallmark regex (includes holin/endolysin), used only for the tier comparison
RE_HALLMARK_lenient = re.compile(
    r"(terminase|portal|capsid|head[- _]protein|tail|baseplate|tape[- _]?measure|virion|"
    r"phage[- _]?(?:coat|sheath|neck|fiber|spike|tube|wedge|hub|lysozyme|structural)|"
    r"prohead|scaffolding|holin|endolysin|muramidase|DTR|MCP|major[- _]?capsid|"
    r"gp\d+.*(?:tail|head|capsid|portal))",
    re.I,
)
RE_RECOMB = re.compile(
    r"(integrase|recombinase|resolvase|excisionase|transposase|Phage_integrase|"
    r"\bHTH_Tnp|\brve\b|site[- _]?specific.*recomb|invertase|insertion[- _]?sequence|"
    r"prophage[- _]?repressor|\bcro\b|lytic[- _]?repressor)",
    re.I,
)


def open_text(path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8", errors="ignore") if p.endswith(".gz") else open(p, "rt", encoding="utf-8", errors="ignore")


def interval_overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def parse_checkv_viral_intervals(region_types_str, region_coords_str):
    out = []
    if not region_types_str or region_types_str == "NA" or not region_coords_str or region_coords_str == "NA":
        return out
    for rt, rc in zip(region_types_str.split(","), region_coords_str.split(",")):
        if rt.strip().lower() == "viral":
            m = re.match(r"^(\d+)-(\d+)$", rc.strip())
            if m:
                out.append((int(m.group(1)), int(m.group(2))))
    return out


def parse_checkv_viral_regions_with_genes(row):
    """CheckV contamination.tsv viral regions as (start, end, viral_genes, host_genes)."""
    types = (row.get("region_types") or "").split(",")
    coords = (row.get("region_coords_bp") or "").split(",")
    vgs = (row.get("region_viral_genes") or "").split(",")
    hgs = (row.get("region_host_genes") or "").split(",")
    out = []
    for i, rt in enumerate(types):
        if rt.strip().lower() != "viral" or i >= len(coords):
            continue
        m = re.match(r"^(\d+)-(\d+)$", coords[i].strip())
        if not m:
            continue
        vg = int(vgs[i]) if i < len(vgs) and vgs[i].strip().isdigit() else 0
        hg = int(hgs[i]) if i < len(hgs) and hgs[i].strip().isdigit() else 0
        out.append((int(m.group(1)), int(m.group(2)), vg, hg))
    return out


def classify_viral_origin_tier(
    is_provirus_boundary_ok,
    n_strict_hm_reg,
    n_rep_reg,
    gn_present,
    gn_score,
    gn_nhm,
    cv_vg,
    cv_hg,
    host_ok,
    host_basis,
    flanked_strict_both,
    flanked_sr_both,
):
    """Strict structural-hallmark-gated Viral Origin Confidence Tier assignment.

    1. Tier 1 requires strictly `n_strict_hm_reg >= 1` (ORF-verified non-recombinase viral
       structural hallmark inside the viral region; raw geNomad `gn_nhm` is NOT sufficient for Tier 1).
    2. For `genomad_provirus_skip` intervals (where CheckV did not split the provirus and contig
       host_genes are skipped), Tier 1 requires double-sided viral structural/replication flanking
       (`flanked_sr_both`) and `cv_vg >= 1` to rule out boundary over-extension into host flanks.
    3. For loci flanked on BOTH sides by strict viral structural hallmarks (`flanked_strict_both`),
       CheckV's self-counting of the focal SBL/PBP as 1 host gene is accounted for:
       `cv_vg >= 1 and cv_vg >= cv_hg - 1` qualifies for Tier 1 (when boundary_ok), and
       any double-structurally-flanked locus with `gn_score >= 0.70` qualifies for at least Tier 2.
    """
    if host_basis == "genomad_provirus_skip":
        t1_host_ok = flanked_sr_both and (cv_vg >= 1)
        t2_host_ok = True
    else:
        t1_host_ok = (cv_vg >= 1) and (host_ok or (flanked_strict_both and cv_vg >= cv_hg - 1))
        t2_host_ok = host_ok or flanked_strict_both

    if is_provirus_boundary_ok and (n_strict_hm_reg >= 1) and t1_host_ok:
        return "Tier1_HighConfidence_Viral"
    has_any_hm = (n_strict_hm_reg >= 1) or (gn_nhm >= 1)
    if gn_present and gn_score >= 0.70 and t2_host_ok and (n_rep_reg >= 1 or cv_vg >= 1 or has_any_hm):
        return "Tier2_Supported_Viral"
    return "Tier3_Fragment_or_Borderline"


def classify_tier_lenient_rule(gn_present, gn_score, gn_topo, cv_qual, n_hm_lenient, gn_nhm_best, min_dist, cv_vg, n_rep_reg):
    """Lenient tier rule (>=2 hallmarks incl. holin/endolysin, 2 kb boundary). Comparison only."""
    is_complete_topo = gn_topo in ("DTR", "ITR") or cv_qual in ("Complete", "High-quality", "Medium-quality")
    has_strong_hm = (n_hm_lenient >= 2) or (gn_nhm_best >= 2)
    if (gn_score >= 0.90 or is_complete_topo) and has_strong_hm and (min_dist >= 2000 or gn_topo in ("DTR", "ITR")):
        return "Tier1_HighConfidence_Viral"
    if gn_present and (n_hm_lenient >= 1 or gn_nhm_best >= 1 or cv_vg >= 1 or n_rep_reg >= 1):
        return "Tier2_Supported_Viral"
    return "Tier3_Fragment_or_Borderline"


def genomad_entry_span(row):
    """Return (start, end) of a geNomad virus_summary entry, or None for whole-contig entries."""
    m = re.search(r"\|provirus_(\d+)_(\d+)$", row.get("seq_name", ""))
    if m:
        return int(m.group(1)), int(m.group(2))
    coords = (row.get("coordinates") or "NA").strip()
    m = re.match(r"^(\d+)-(\d+)$", coords)
    return (int(m.group(1)), int(m.group(2))) if m else None


def pick_genomad_entry(entries, st, en):
    """Choose the geNomad entry containing the SBL ORF: overlapping provirus > whole contig > best score."""
    if not entries:
        return {}
    if st > 0:
        overl = [(interval_overlap(st, en, *sp), e) for e in entries if (sp := genomad_entry_span(e))]
        overl = [x for x in overl if x[0] > 0]
        if overl:
            return max(overl, key=lambda x: (x[0], float(x[1]["virus_score"])))[1]
    whole = [e for e in entries if genomad_entry_span(e) is None]
    if whole:
        return max(whole, key=lambda e: float(e["virus_score"]))
    return {}


def main():
    parser = argparse.ArgumentParser(description="Phase 2 strict viral genomic neighborhood & Tier assignment")
    parser.add_argument("--catalog", default=str(P1_OUT / "Final_Viral_SBL_Catalog.tsv"))
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/7] Loading viral catalog & decontamination audit coordinates...")
    decontam = {}
    with open(P1_DEC, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            decontam[r["sbl_protein_id"]] = r

    viral_rows = []
    with open(args.catalog, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cat = r["pfam_category"]
            if not (cat.startswith("focal_") or cat.startswith("aux_")):
                continue
            viral_rows.append(r)

    target_uvigs = {decontam[r["orig_uid"]]["uvig"] for r in viral_rows if r["source"] == "IMG_VR" and r["orig_uid"] in decontam}
    print(f"  Candidate non-drop viral loci: {len(viral_rows)} across {len(target_uvigs)} IMG/VR contigs (+ MOSAiC={sum(1 for r in viral_rows if r['source']=='MOSAiC')}).")

    print("[2/7] Indexing Prodigal ORFs on target contigs...")
    contig_orfs = defaultdict(list)
    target_orf_ids = set()
    with open(PRODIGAL_ORFS, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            hdr = line[1:].strip()
            p = hdr.split(" # ")
            oid = p[0].strip()
            uvig = oid.rsplit("_", 1)[0]
            if uvig not in target_uvigs or len(p) < 4:
                continue
            st, en, strand = int(p[1]), int(p[2]), int(p[3])
            m_gc = re.search(r"gc_cont=([0-9.]+)", p[4]) if len(p) > 4 else None
            gc = m_gc.group(1) if m_gc else "NA"
            contig_orfs[uvig].append({
                "orf_id": oid,
                "start": st,
                "end": en,
                "strand": strand,
                "gc": gc,
                "pfams": [],
                "eggnog": None,
                "gn_hallmark": False,
                "gn_plasmid": False,
                "gn_desc": "",
            })
            target_orf_ids.add(oid)
    for uvig in contig_orfs:
        contig_orfs[uvig].sort(key=lambda o: (o["start"], o["end"]))
    print(f"  Indexed {len(target_orf_ids):,} ORFs across {len(contig_orfs):,} target contigs.")

    orf_by_id = {o["orf_id"]: o for lst in contig_orfs.values() for o in lst}

    print("[3/7] Loading full-genome Pfam-A (--cut_ga) annotations...")
    n_pfam = 0
    with open(ORFS_PFAM, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 22:
                continue
            oid = p[0]
            if oid not in target_orf_ids:
                continue
            acc = p[4].split(".")[0]
            name = p[3]
            sc = float(p[13])
            st, en = int(p[19]), int(p[20])
            orf_by_id[oid]["pfams"].append((st, en, sc, f"{acc}:{name}"))
            n_pfam += 1
    print(f"  Loaded {n_pfam:,} Pfam domain hits on target contig ORFs.")

    print("[4/7] Loading full-genome eggNOG-mapper annotations...")
    n_egg = 0
    with open_text(EGGNOG_GZ) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            oid = p[0]
            if oid not in target_orf_ids:
                continue
            cog = p[6] if p[6] != "-" else ""
            desc = p[7] if p[7] != "-" else ""
            pref = p[8] if len(p) > 8 and p[8] != "-" else ""
            epf = p[20] if len(p) > 20 and p[20] != "-" else ""
            orf_by_id[oid]["eggnog"] = (cog, pref, desc, epf)
            n_egg += 1
    print(f"  Loaded {n_egg:,} eggNOG annotations on target contig ORFs.")

    print("[5/7] Loading geNomad v1.7 gene hallmarks & summaries...")
    gn_vs = {}                       # best-scoring entry per contig (lenient-rule comparison only)
    gn_entries = defaultdict(list)   # all entries per contig (provirus regions + whole contig)
    with open_text(GENOMAD_DIR / "sbl_contigs_virus_summary.tsv.gz") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            uv = r["seq_name"].split("|")[0]
            if uv in target_uvigs:
                gn_entries[uv].append(r)
                if uv not in gn_vs or float(r["virus_score"]) > float(gn_vs[uv]["virus_score"]):
                    gn_vs[uv] = r

    gn_pv = defaultdict(list)
    with open_text(GENOMAD_DIR / "sbl_contigs_provirus.tsv.gz") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["source_seq"] in target_uvigs:
                gn_pv[r["source_seq"]].append((int(r["start"]), int(r["end"])))

    gn_genes_file = GENOMAD_DIR / "genomad_out" / "sbl_contigs_annotate" / "sbl_contigs_genes.tsv"
    n_gn_map = 0
    with open(gn_genes_file, "r", encoding="utf-8", errors="ignore") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gid = r["gene"]
            uv = gid.rsplit("_", 1)[0]
            if uv not in target_uvigs:
                continue
            vh = r.get("virus_hallmark", "0") == "1"
            ph = r.get("plasmid_hallmark", "0") == "1"
            desc = r.get("annotation_description", "NA")
            acc = r.get("annotation_accessions", "NA")
            if desc == "NA" and acc == "NA" and not vh and not ph:
                continue
            txt = f"{acc}:{desc}" if desc != "NA" else acc
            o = orf_by_id.get(gid)
            if o is None:
                gs, ge = int(r["start"]), int(r["end"])
                best_ov, best_o = 0, None
                for cand in contig_orfs.get(uv, []):
                    ov = interval_overlap(gs, ge, cand["start"], cand["end"])
                    if ov > best_ov:
                        best_ov = ov
                        best_o = cand
                if best_o is not None and best_ov >= 0.5 * min(ge - gs + 1, best_o["end"] - best_o["start"] + 1):
                    o = best_o
            if o is not None:
                o["gn_hallmark"] = o["gn_hallmark"] or vh
                o["gn_plasmid"] = o["gn_plasmid"] or ph
                if txt != "NA":
                    o["gn_desc"] = txt if not o["gn_desc"] else f"{o['gn_desc']};{txt}"
                n_gn_map += 1
    print(f"  Mapped {n_gn_map:,} annotated geNomad genes onto target contig ORFs.")

    print("[6/7] Loading CheckV v1.0.3 quality & viral boundaries...")
    cv_qs = {}
    with open(CHECKV_DIR / "quality_summary.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["contig_id"] in target_uvigs:
                cv_qs[r["contig_id"]] = r

    cv_v_regions = {}
    cv_v_regions_genes = {}
    with open(CHECKV_DIR / "contamination.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = r["contig_id"]
            if cid in target_uvigs and r.get("provirus", "No") == "Yes":
                ivs = parse_checkv_viral_intervals(r.get("region_types", ""), r.get("region_coords_bp", ""))
                if ivs:
                    cv_v_regions[cid] = ivs
                    cv_v_regions_genes[cid] = parse_checkv_viral_regions_with_genes(r)

    # Classify each ORF strictly across the 3 orthogonal categories
    for uvig, orfs in contig_orfs.items():
        for o in orfs:
            pf_sorted = sorted(o["pfams"], key=lambda x: (-x[2], x[0]))
            kept_pf = []
            for st, en, sc, label in pf_sorted:
                if not any(interval_overlap(st, en, kst, ken) >= 20 for kst, ken, _ in kept_pf):
                    kept_pf.append((st, en, label))
            pf_str = "+".join(x[2] for x in sorted(kept_pf))
            egg = o["eggnog"]
            egg_str = "|".join(x for x in (egg[1], egg[2], egg[3]) if x and x != "NA") if egg else ""
            combined_txt = " ".join(x for x in (pf_str, egg_str, o["gn_desc"]) if x)
            o["label"] = pf_str or (egg[1] if egg and egg[1] else "") or (o["gn_desc"][:40] if o["gn_desc"] else "hypothetical")
            o["is_recomb"] = bool(RE_RECOMB.search(combined_txt))
            o["is_replic"] = bool(RE_REPLIC_LYSIS.search(combined_txt))
            # Strict viral structural hallmark: explicit structural annotation OR non-recomb geNomad virus hallmark
            is_strict_txt = bool(RE_STRICT_STRUCTURAL.search(combined_txt))
            o["is_strict_structural_hallmark"] = is_strict_txt or (o["gn_hallmark"] and not o["is_recomb"])
            o["is_hallmark_lenient"] = o["gn_hallmark"] or bool(RE_HALLMARK_lenient.search(combined_txt))
            o["is_structural_or_replic"] = o["is_strict_structural_hallmark"] or o["is_replic"]
            o["has_any_annot"] = bool(combined_txt.strip())

    print("[7/7] Computing locus-level neighborhood features & strict confidence tiers...")
    out_rows = []
    tier_counts = Counter()
    tier_shift = Counter()   # (lenient-rule tier, strict tier) on the primary IMG/VR cohort

    for r in viral_rows:
        uid = r["uid"]
        orig_uid = r["orig_uid"]
        src = r["source"]
        fam = r["winning_acc"]
        cat = r["pfam_category"]
        in_primary = (r.get("pass_fimo_1e5") == "True") and (r.get("pass_cov65") == "True")

        if src != "IMG_VR" or orig_uid not in decontam:
            tier = "MOSAiC_Proteome_Verified_vMAG"
            if in_primary:
                tier_counts[(cat, tier)] += 1
            out_rows.append({
                "uid": uid, "orig_uid": orig_uid, "source": src, "uvig": "NA",
                "pfam_category": cat, "winning_acc": fam,
                "dual_PF13354_PF00768_GA": r.get("dual_PF13354_PF00768_GA", "NA"),
                "in_primary_cov65_fimo1e5": in_primary,
                "contig_context_available": False, "viral_origin_tier": tier,
                "viral_origin_tier_lenient_rule": tier, "provirus_interval_used": "NA",
                "host_gene_check_basis": "NA", "host_gene_check_pass": "NA",
                "contig_length_bp": "NA", "total_orfs_on_contig": "NA",
                "sbl_start_bp": "NA", "sbl_end_bp": "NA", "sbl_strand": "NA",
                "dist_to_5p_viral_bp": "NA", "dist_to_3p_viral_bp": "NA", "min_dist_to_viral_boundary_bp": "NA",
                "genomad_virus_score": "NA", "genomad_topology": "NA", "genomad_n_hallmarks": "NA",
                "genomad_taxonomy": r.get("vir_sub_info", "NA"), "checkv_quality": "NA", "checkv_miuvig_quality": "NA",
                "checkv_viral_genes": "NA", "checkv_host_genes": "NA", "checkv_contamination_pct": "NA",
                "n_orfs_in_viral_region": "NA", "pct_orfs_annotated_pfam_or_eggnog": "NA",
                "n_viral_hallmark_in_region": "NA", "n_viral_hallmark_10kb": "NA",
                "n_viral_replic_lysis_in_region": "NA", "n_recomb_lysogeny_in_region": "NA",
                "flanked_strict_structural_hallmark_any": "NA",
                "flanked_strict_structural_hallmark_both": "NA",
                "flanked_viral_structural_or_replic_any": "NA",
                "flanked_viral_structural_or_replic_both": "NA",
                "flanked_recomb_lysogeny_any": "NA",
                "nearest_hallmark_dist_bp": "NA", "nearest_recomb_dist_bp": "NA",
                "upstream_5_orfs": "NA", "downstream_5_orfs": "NA",
            })
            continue

        dec = decontam[orig_uid]
        uvig = dec["uvig"]
        st = int(dec["orf_start"]) if dec["orf_start"] != "NA" else 0
        en = int(dec["orf_end"]) if dec["orf_end"] != "NA" else 0

        cq = cv_qs.get(uvig, {})
        gn = pick_genomad_entry(gn_entries.get(uvig, []), st, en)   # entry that contains this SBL
        gn_best = gn_vs.get(uvig, {})                                # lenient-rule comparison only
        orfs = contig_orfs.get(uvig, [])

        clen = int(cq.get("contig_length", 0) or 0)
        if clen == 0 and orfs:
            clen = max(o["end"] for o in orfs)

        v_st, v_en = 1, clen
        has_pv_interval = False
        if uvig in gn_pv and st > 0:
            cands = [(ps, pe) for ps, pe in gn_pv[uvig] if interval_overlap(st, en, ps, pe) > 0]
            if cands:
                best_gn = max(cands, key=lambda x: interval_overlap(st, en, x[0], x[1]))
                v_st = max(v_st, best_gn[0])
                v_en = min(v_en, best_gn[1])
                has_pv_interval = True
        if uvig in cv_v_regions and st > 0:
            cands = [(ps, pe) for ps, pe in cv_v_regions[uvig] if interval_overlap(st, en, ps, pe) > 0]
            if cands:
                best_cv = max(cands, key=lambda x: interval_overlap(st, en, x[0], x[1]))
                v_st = max(v_st, best_cv[0])
                v_en = min(v_en, best_cv[1])
                has_pv_interval = True
        if v_st > v_en:
            v_st, v_en = 1, clen

        dist_5p = max(0, st - v_st) if st > 0 else 0
        dist_3p = max(0, v_en - en) if en > 0 else 0
        min_dist = min(dist_5p, dist_3p) if st > 0 else 0

        reg_orfs = [o for o in orfs if interval_overlap(o["start"], o["end"], v_st, v_en) > 0]
        sbl_idx = None
        sbl_strand = 0
        if st > 0 and reg_orfs:
            best_i, best_ov = None, 0
            for idx, o in enumerate(reg_orfs):
                ov = interval_overlap(st, en, o["start"], o["end"])
                if ov > best_ov:
                    best_ov = ov
                    best_i = idx
            if best_i is not None:
                sbl_idx = best_i
                sbl_strand = reg_orfs[best_i]["strand"]

        n_annot = sum(1 for o in reg_orfs if o["has_any_annot"])
        pct_annot = round(100.0 * n_annot / max(len(reg_orfs), 1), 2)

        neigh = [o for idx, o in enumerate(reg_orfs) if idx != sbl_idx]
        n_hm_reg = sum(1 for o in neigh if o["is_strict_structural_hallmark"])
        n_rep_reg = sum(1 for o in neigh if o["is_replic"])
        n_rec_reg = sum(1 for o in neigh if o["is_recomb"])
        n_hm_reg_lenient = sum(1 for o in neigh if o["is_hallmark_lenient"])

        def orf_dist(o):
            if st == 0:
                return 9999999
            if o["end"] < st:
                return st - o["end"]
            if o["start"] > en:
                return o["start"] - en
            return 0

        hm_dists = [orf_dist(o) for o in neigh if o["is_strict_structural_hallmark"]]
        rec_dists = [orf_dist(o) for o in neigh if o["is_recomb"]]
        n_hm_10k = sum(1 for d in hm_dists if d <= 10000)
        nearest_hm = min(hm_dists) if hm_dists else "NA"
        nearest_rec = min(rec_dists) if rec_dists else "NA"

        if sbl_idx is not None:
            left_orfs = reg_orfs[:sbl_idx]
            right_orfs = reg_orfs[sbl_idx + 1 :]
            left_hm = any(o["is_strict_structural_hallmark"] for o in left_orfs)
            right_hm = any(o["is_strict_structural_hallmark"] for o in right_orfs)
            left_sr = any(o["is_structural_or_replic"] for o in left_orfs)
            right_sr = any(o["is_structural_or_replic"] for o in right_orfs)
            any_rec = any(o["is_recomb"] for o in neigh)
            up5 = "; ".join(f"{o['orf_id'].rsplit('_',1)[-1]}({o['label']})" for o in left_orfs[-5:]) or "NONE"
            dn5 = "; ".join(f"{o['orf_id'].rsplit('_',1)[-1]}({o['label']})" for o in right_orfs[:5]) or "NONE"
        else:
            left_hm = right_hm = left_sr = right_sr = any_rec = False
            up5, dn5 = "UNMAPPED", "UNMAPPED"

        gn_score = float(gn.get("virus_score") or dec.get("genomad_virus_score") or 0.0)
        gn_nhm = int(gn.get("n_hallmarks") or 0)
        gn_topo = gn.get("topology", "Not_in_virus_summary")
        gn_tax = gn.get("taxonomy", "NA")
        cv_qual = cq.get("checkv_quality", "NA")
        cv_miuvig = cq.get("miuvig_quality", "NA")
        cv_vg = int(float(cq.get("viral_genes", 0) or 0))
        cv_hg = int(float(cq.get("host_genes", 0) or 0))
        cv_ctm = float(cq.get("contamination", 0.0) or 0.0)

        # The SBL must sit >= 1 kb inside the viral region (provirus interval or, for whole-contig
        # viruses, the contig ends); contig length alone is no longer accepted as boundary evidence.
        is_boundary_ok = bool(gn) and (gn_topo in ("DTR", "ITR") or min_dist >= BOUNDARY_MIN_BP)

        # Host-vs-viral gene check on the region that actually contains the SBL:
        #  - inside a CheckV viral region -> that region's own viral/host gene counts;
        #  - inside a geNomad-only provirus -> contig-level counts include the host flank, so the
        #    comparison is not informative and is skipped (viral genes >= 1 still required);
        #  - whole-contig virus -> contig-level CheckV counts.
        cv_reg = [x for x in cv_v_regions_genes.get(uvig, []) if st > 0 and interval_overlap(st, en, x[0], x[1]) > 0]
        if cv_reg:
            _, _, reg_vg, reg_hg = max(cv_reg, key=lambda x: interval_overlap(st, en, x[0], x[1]))
            vg_eval, hg_eval, host_ok, host_basis = reg_vg, reg_hg, reg_vg >= reg_hg, "checkv_viral_region"
        elif has_pv_interval:
            vg_eval, hg_eval, host_ok, host_basis = cv_vg, cv_hg, True, "genomad_provirus_skip"
        else:
            vg_eval, hg_eval, host_ok, host_basis = cv_vg, cv_hg, cv_vg >= cv_hg, "whole_contig"

        tier = classify_viral_origin_tier(
            is_provirus_boundary_ok=is_boundary_ok,
            n_strict_hm_reg=n_hm_reg,
            n_rep_reg=n_rep_reg,
            gn_present=bool(gn),
            gn_score=gn_score,
            gn_nhm=gn_nhm,
            cv_vg=vg_eval,
            cv_hg=hg_eval,
            host_ok=host_ok,
            host_basis=host_basis,
            flanked_strict_both=(left_hm and right_hm),
            flanked_sr_both=(left_sr and right_sr),
        )

        tier_lenient = classify_tier_lenient_rule(
            gn_present=bool(gn_best),
            gn_score=float(gn_best.get("virus_score") or dec.get("genomad_virus_score") or 0.0),
            gn_topo=gn_best.get("topology", "Not_in_virus_summary"),
            cv_qual=cv_qual,
            n_hm_lenient=n_hm_reg_lenient,
            gn_nhm_best=int(gn_best.get("n_hallmarks") or 0),
            min_dist=min_dist,
            cv_vg=cv_vg,
            n_rep_reg=n_rep_reg,
        )

        if in_primary:
            tier_counts[(cat, tier)] += 1
            tier_shift[(tier_lenient, tier)] += 1

        out_rows.append({
            "uid": uid,
            "orig_uid": orig_uid,
            "source": src,
            "uvig": uvig,
            "pfam_category": cat,
            "winning_acc": fam,
            "dual_PF13354_PF00768_GA": r.get("dual_PF13354_PF00768_GA", "NA"),
            "in_primary_cov65_fimo1e5": in_primary,
            "contig_context_available": True,
            "viral_origin_tier": tier,
            "viral_origin_tier_lenient_rule": tier_lenient,
            "provirus_interval_used": has_pv_interval,
            "host_gene_check_basis": host_basis,
            "host_gene_check_pass": host_ok,
            "contig_length_bp": clen,
            "total_orfs_on_contig": len(orfs),
            "sbl_start_bp": st,
            "sbl_end_bp": en,
            "sbl_strand": sbl_strand,
            "dist_to_5p_viral_bp": dist_5p,
            "dist_to_3p_viral_bp": dist_3p,
            "min_dist_to_viral_boundary_bp": min_dist,
            "genomad_virus_score": gn_score,
            "genomad_topology": gn_topo,
            "genomad_n_hallmarks": gn_nhm,
            "genomad_taxonomy": gn_tax,
            "checkv_quality": cv_qual,
            "checkv_miuvig_quality": cv_miuvig,
            "checkv_viral_genes": cv_vg,
            "checkv_host_genes": cv_hg,
            "checkv_contamination_pct": cv_ctm,
            "n_orfs_in_viral_region": len(reg_orfs),
            "pct_orfs_annotated_pfam_or_eggnog": pct_annot,
            "n_viral_hallmark_in_region": n_hm_reg,
            "n_viral_hallmark_10kb": n_hm_10k,
            "n_viral_replic_lysis_in_region": n_rep_reg,
            "n_recomb_lysogeny_in_region": n_rec_reg,
            "flanked_strict_structural_hallmark_any": left_hm or right_hm,
            "flanked_strict_structural_hallmark_both": left_hm and right_hm,
            "flanked_viral_structural_or_replic_any": left_sr or right_sr,
            "flanked_viral_structural_or_replic_both": left_sr and right_sr,
            "flanked_recomb_lysogeny_any": any_rec,
            "nearest_hallmark_dist_bp": nearest_hm,
            "nearest_recomb_dist_bp": nearest_rec,
            "upstream_5_orfs": up5,
            "downstream_5_orfs": dn5,
        })

    out_tsv = OUT_DIR / "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"
    with open(out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    img_prim = [x for x in out_rows if x["contig_context_available"] and x["in_primary_cov65_fimo1e5"]]
    focal_img = [x for x in img_prim if x["pfam_category"].startswith("focal_")]
    lines = [
        "=== Phase 2: 100%-Coverage Viral Genomic Neighborhood Summary (Strict Structural Hallmarks) ===",
        f"Total primary curated viral loci (cov>=0.65 & FIMO<=1e-5): {sum(1 for x in out_rows if x['in_primary_cov65_fimo1e5']):,} (IMG/VR: {len(img_prim):,}, MOSAiC: {sum(1 for x in out_rows if x['in_primary_cov65_fimo1e5'] and not x['contig_context_available'])})",
        f"  Focal IMG/VR (PF00144+PF00905+PF13354): {len(focal_img):,} | Sister PRE IMG/VR (PF00768): {len(img_prim)-len(focal_img):,}",
        f"Median % ORFs annotated per viral contig (Pfam-A + eggNOG): {statistics.median(x['pct_orfs_annotated_pfam_or_eggnog'] for x in img_prim):.2f}%",
        f"Strict viral structural hallmark >=1 in viral region (capsid/portal/terminase/tail/geNomad hallmark): {sum(1 for x in img_prim if x['n_viral_hallmark_in_region'] >= 1):,} / {len(img_prim):,} ({100*sum(1 for x in img_prim if x['n_viral_hallmark_in_region'] >= 1)/max(len(img_prim),1):.2f}%)",
        f"Flanked by STRICT viral structural hallmark on >=1 side: {sum(1 for x in img_prim if x['flanked_strict_structural_hallmark_any']):,} / {len(img_prim):,} ({100*sum(1 for x in img_prim if x['flanked_strict_structural_hallmark_any'])/max(len(img_prim),1):.2f}%)",
        f"Flanked by STRICT viral structural hallmark on BOTH sides: {sum(1 for x in img_prim if x['flanked_strict_structural_hallmark_both']):,} / {len(img_prim):,} ({100*sum(1 for x in img_prim if x['flanked_strict_structural_hallmark_both'])/max(len(img_prim),1):.2f}%)",
        f"Flanked by viral structural OR replication/lysis (NO integrase/recomb) on >=1 side: {sum(1 for x in img_prim if x['flanked_viral_structural_or_replic_any']):,} / {len(img_prim):,} ({100*sum(1 for x in img_prim if x['flanked_viral_structural_or_replic_any'])/max(len(img_prim),1):.2f}%)",
        f"Flanked by viral structural OR replication/lysis (NO integrase/recomb) on BOTH sides: {sum(1 for x in img_prim if x['flanked_viral_structural_or_replic_both']):,} / {len(img_prim):,} ({100*sum(1 for x in img_prim if x['flanked_viral_structural_or_replic_both'])/max(len(img_prim),1):.2f}%)",
        f"Recombination / integrase / lysogeny gene in region (reported separately as MGE background): {sum(1 for x in img_prim if x['flanked_recomb_lysogeny_any']):,} / {len(img_prim):,} ({100*sum(1 for x in img_prim if x['flanked_recomb_lysogeny_any'])/max(len(img_prim),1):.2f}%)",
        f"Median min distance to viral region boundary (IMG/VR): {statistics.median(x['min_dist_to_viral_boundary_bp'] for x in img_prim):.1f} bp",
        "\nBreakdown by (Pfam Category, Strict Viral Origin Confidence Tier) on Primary Curated Cohort:",
    ]
    for (cat, tr), cnt in sorted(tier_counts.items()):
        lines.append(f"  {cat:16s} | {tr:30s}: {cnt:,}")
    basis = Counter((x["host_gene_check_basis"], x["host_gene_check_pass"]) for x in img_prim)
    lines.append("\nHost-vs-viral gene check (basis, pass) on primary IMG/VR loci:")
    for (bs, ok), cnt in sorted(basis.items(), key=lambda kv: str(kv[0])):
        lines.append(f"  {bs:24s} pass={ok!s:5s}: {cnt:,}")
    lines.append("\nTier shift, lenient rule (rows) -> strict rule (columns), primary IMG/VR loci:")
    tiers = ["Tier1_HighConfidence_Viral", "Tier2_Supported_Viral", "Tier3_Fragment_or_Borderline"]
    lines.append("  " + " " * 30 + "".join(f"{t.split('_')[0]:>10s}" for t in tiers))
    for old in tiers:
        lines.append(f"  {old:30s}" + "".join(f"{tier_shift.get((old, new), 0):>10,}" for new in tiers))
    rep = "\n".join(lines) + "\n"
    (OUT_DIR / "Phase2_Genomic_Context_Summary.txt").write_text(rep, encoding="utf-8")
    print(rep)


if __name__ == "__main__":
    main()
