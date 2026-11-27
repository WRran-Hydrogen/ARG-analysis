#!/usr/bin/env python3
"""Phase 2: 100%-Coverage Viral Genomic Neighborhood & Strict Structural Confidence Tier Integration.

Integrates for all curated IMG/VR viral SBL/PRE loci:
- Prodigal meta ORFs (all_orfs.faa)
- Full-genome Pfam-A (--cut_ga) annotations (orfs_pfam.domtbl)
- Full-genome eggNOG-mapper annotations (sbl_neighbourhood.emapper.annotations.gz)
- geNomad v1.7 summary, hallmarks, and provirus boundaries
- CheckV v1.0.3 (DB v1.5) quality and viral/host region coordinates

Gene origin evidence (per ORF, no free-text matching):
- viral evidence: CheckV HMM category viral (tmp/gene_features.tsv hmm_cat = 1), a geNomad marker
  whose specificity class is virus-first (VV/VC/VP, genomad_marker_metadata.tsv), an eggNOG seed
  ortholog in Viruses (eggNOG_OGs level 10239), or a geNomad virus hallmark;
- host evidence: CheckV microbial (hmm_cat = -1) or a chromosome-first geNomad marker (CC/CP), and
  no viral evidence;
- viral hallmark: geNomad virus hallmark, excluding recombination/transposition markers (marker
  description: integrase, recombinase, resolvase, invertase, excisionase, transposase).

Tier rules (primary, `viral_origin_tier`):
- Tier1_HighConfidence_Viral: SBL >= 1 kb inside the viral region (or DTR/ITR complete genome)
  + >= 1 viral hallmark in the viral region + CheckV viral_genes >= 1 + host-gene check
  `viral_genes >= host_genes - f` on the CheckV viral region containing the SBL (whole contig for
  non-provirus contigs), f = 1 when CheckV calls the focal ORF itself a host gene, else 0; for
  `genomad_provirus_skip` intervals (contig-level counts include the host flank) viral-evidence
  ORFs on both sides of the SBL within the interval plus viral_genes >= 1 replace the count check.
- Tier2_Supported_Viral: geNomad score >= 0.70 + host-gene check as above (not required for
  `genomad_provirus_skip`) + any viral evidence (viral-evidence ORF, CheckV viral gene, hallmark).
- Tier3_Fragment_or_Borderline: all other loci.
Sensitivity (`viral_origin_tier_local_w10`): Tier 1 = boundary + hallmark + viral-evidence ORF on
both sides within 10 ORFs + viral >= host evidence in that window; Tier 2 = geNomad score >= 0.70
+ >= 1 viral-evidence ORF within 10 ORFs and viral >= host evidence there.
Comparison columns: `viral_origin_tier_text_rule` (rule used until 2026-11-25, annotation-text
categories and a flank-based host-gene allowance) and `viral_origin_tier_lenient_rule` (>= 2
hallmarks including holin/endolysin, 2 kb boundary). Text categories (`*_text_*` columns) are kept
for labels only. geNomad hallmark counts per entry come from the geNomad entry containing the SBL.

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

CHECKV_GENES = CHECKV_DIR / "tmp" / "gene_features.tsv"
GENOMAD_MARKER_META = hpc.DB_DIR / "genomad_db" / "genomad_marker_metadata.tsv"

OUT_DIR = BASE / "02_genomic_context"
BOUNDARY_MIN_BP = 1000   # minimum distance from the SBL ORF to the viral region boundary for Tier 1
LOCAL_WINDOW_ORFS = 10   # ORFs on each side of the SBL for the local-evidence sensitivity tier
VIRAL_MARKER_CLASSES = ("VV", "VC", "VP")
HOST_MARKER_CLASSES = ("CC", "CP")
RE_RECOMB_MARKER = re.compile(r"\b(integrase|recombinase|resolvase|invertase|excisionase|transposase)\b", re.I)

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


TIER1, TIER2, TIER3 = "Tier1_HighConfidence_Viral", "Tier2_Supported_Viral", "Tier3_Fragment_or_Borderline"


def classify_viral_origin_tier(boundary_ok, n_hallmark, n_viral_ev, gn_present, gn_score, gn_nhm,
                               cv_vg, cv_hg, focal_is_host_gene, host_basis, flanked_viral_ev_both):
    """Primary viral-origin tier (gene-origin evidence, focal-corrected host-gene check).

    `n_hallmark`: geNomad virus hallmarks (no recombination markers) in the viral region, focal excluded.
    `n_viral_ev`: viral-evidence ORFs in the viral region, focal excluded.
    Host-gene check: `cv_vg >= cv_hg - f` with f = 1 when CheckV calls the focal ORF a host gene;
    skipped for `genomad_provirus_skip` (contig-level counts include the host flank), where Tier 1
    instead needs viral-evidence ORFs on both sides of the SBL and `cv_vg >= 1`.
    """
    if host_basis == "genomad_provirus_skip":
        t1_host_ok = flanked_viral_ev_both and cv_vg >= 1
        t2_host_ok = True
    else:
        ok = cv_vg >= cv_hg - (1 if focal_is_host_gene else 0)
        t1_host_ok = cv_vg >= 1 and ok
        t2_host_ok = ok
    if boundary_ok and n_hallmark >= 1 and t1_host_ok:
        return TIER1
    if gn_present and gn_score >= 0.70 and t2_host_ok and (n_viral_ev >= 1 or cv_vg >= 1 or n_hallmark >= 1 or gn_nhm >= 1):
        return TIER2
    return TIER3


def classify_tier_local(boundary_ok, n_hallmark, gn_present, gn_score, left, right, w=LOCAL_WINDOW_ORFS):
    """Sensitivity tier from the w ORFs on each side of the SBL (left/right ordered outward from the SBL)."""
    lw, rw = left[:w], right[:w]
    n_v = sum(o["viral_ev"] for o in lw + rw)
    n_h = sum(o["host_ev"] for o in lw + rw)
    both = any(o["viral_ev"] for o in lw) and any(o["viral_ev"] for o in rw)
    if boundary_ok and n_hallmark >= 1 and both and n_v >= n_h:
        return TIER1, n_v, n_h, both
    if gn_present and gn_score >= 0.70 and n_v >= 1 and n_v >= n_h:
        return TIER2, n_v, n_h, both
    return TIER3, n_v, n_h, both


def classify_tier_text_rule(
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
    """Tier rule used until 2026-11-25 (comparison only): hallmark / replication-lysis categories from
    annotation text, and a `cv_vg >= cv_hg - 1` allowance for loci flanked on both sides by text
    structural hallmarks."""
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
    parser.add_argument("--out-dir", default=str(OUT_DIR), help="output folder (default: the Phase 2 folder)")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

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
                "gn_hallmark_recomb": False,
                "gn_plasmid": False,
                "gn_desc": "",
                "gn_marker_class": "none",
                "egg_domain": "no_hit",
                "cv_cat": None,
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
            # domain of the seed ortholog's lineage (eggNOG_OGs lists one OG per level of that lineage)
            ogs = p[4]
            orf_by_id[oid]["egg_domain"] = ("Viruses" if "@10239|Viruses" in ogs else "Bacteria" if "@2|Bacteria" in ogs
                                            else "Archaea" if "@2157|Archaea" in ogs else "Eukaryota" if "@2759|Eukaryota" in ogs
                                            else "root_only")
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

    def map_to_orf(uv, gs, ge, gid=None):
        """Prodigal ORF for a gene called by another tool: same id and coordinates, else the ORF
        with the largest overlap covering >= 50 % of the shorter of the two."""
        o = orf_by_id.get(gid) if gid else None
        if o is not None and o["start"] == gs and o["end"] == ge:
            return o
        best_ov, best_o = 0, None
        for cand in contig_orfs.get(uv, []):
            ov = interval_overlap(gs, ge, cand["start"], cand["end"])
            if ov > best_ov:
                best_ov, best_o = ov, cand
        if best_o is not None and best_ov >= 0.5 * min(ge - gs + 1, best_o["end"] - best_o["start"] + 1):
            return best_o
        return o if gid else None

    marker_class, recomb_markers = {}, set()
    with open_text(GENOMAD_MARKER_META) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            marker_class[r["MARKER"]] = r["SPECIFICITY_CLASS"]
            if r.get("VIRUS_HALLMARK") == "1" and RE_RECOMB_MARKER.search(r.get("ANNOTATION_DESCRIPTION") or ""):
                recomb_markers.add(r["MARKER"])
    print(f"  geNomad markers: {len(marker_class):,} ({len(recomb_markers)} recombination/transposition hallmarks excluded)")

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
            mk = r.get("marker", "NA")
            desc = r.get("annotation_description", "NA")
            acc = r.get("annotation_accessions", "NA")
            if desc == "NA" and acc == "NA" and not vh and not ph and mk in ("NA", ""):
                continue
            txt = f"{acc}:{desc}" if desc != "NA" else acc
            o = map_to_orf(uv, int(r["start"]), int(r["end"]), gid)
            if o is not None:
                o["gn_hallmark"] = o["gn_hallmark"] or vh
                o["gn_hallmark_recomb"] = o["gn_hallmark_recomb"] or (vh and mk in recomb_markers)
                o["gn_plasmid"] = o["gn_plasmid"] or ph
                if mk not in ("NA", "") and o["gn_marker_class"] == "none":
                    o["gn_marker_class"] = marker_class.get(mk, "none")
                if txt != "NA":
                    o["gn_desc"] = txt if not o["gn_desc"] else f"{o['gn_desc']};{txt}"
                n_gn_map += 1
    print(f"  Mapped {n_gn_map:,} geNomad genes (annotation, hallmark or marker) onto target contig ORFs.")

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

    n_cv_map = 0
    with open(CHECKV_GENES, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = r["contig_id"]
            if cid in target_uvigs and r["hmm_cat"] not in ("NA", ""):
                o = map_to_orf(cid, int(r["start"]), int(r["end"]))
                if o is not None:
                    o["cv_cat"] = int(r["hmm_cat"])   # 1 viral, -1 microbial, 0 no HMM hit
                    n_cv_map += 1
    print(f"  Mapped {n_cv_map:,} CheckV gene calls onto target contig ORFs.")

    # Gene-origin evidence (tiers) and annotation-text categories (labels, comparison rule only)
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
            o["hallmark"] = o["gn_hallmark"] and not o["gn_hallmark_recomb"]
            o["viral_ev"] = (o["cv_cat"] == 1 or o["gn_marker_class"] in VIRAL_MARKER_CLASSES
                             or o["egg_domain"] == "Viruses" or o["hallmark"])
            o["host_ev"] = (o["cv_cat"] == -1 or o["gn_marker_class"] in HOST_MARKER_CLASSES) and not o["viral_ev"]

    print("[7/7] Computing locus-level neighborhood features & strict confidence tiers...")
    out_rows = []
    tier_counts = Counter()
    tier_shift = Counter()        # (lenient-rule tier, primary tier) on the primary IMG/VR cohort
    tier_shift_text = Counter()   # (text-rule tier, primary tier)
    tier_shift_local = Counter()  # (primary tier, local w10 tier)

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
            out_rows.append({   # other columns are written as NA
                "uid": uid, "orig_uid": orig_uid, "source": src, "uvig": "NA",
                "pfam_category": cat, "winning_acc": fam,
                "dual_PF13354_PF00768_GA": r.get("dual_PF13354_PF00768_GA", "NA"),
                "in_primary_cov65_fimo1e5": in_primary,
                "contig_context_available": False, "viral_origin_tier": tier,
                "viral_origin_tier_local_w10": tier, "viral_origin_tier_text_rule": tier,
                "viral_origin_tier_lenient_rule": tier, "genomad_taxonomy": r.get("vir_sub_info", "NA"),
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
        n_gh_reg = sum(1 for o in neigh if o["hallmark"])
        n_vev_reg = sum(1 for o in neigh if o["viral_ev"])
        n_hev_reg = sum(1 for o in neigh if o["host_ev"])

        def orf_dist(o):
            if st == 0:
                return 9999999
            if o["end"] < st:
                return st - o["end"]
            if o["start"] > en:
                return o["start"] - en
            return 0

        gh_dists = [orf_dist(o) for o in neigh if o["hallmark"]]
        rec_dists = [orf_dist(o) for o in neigh if o["is_recomb"]]
        n_gh_10k = sum(1 for d in gh_dists if d <= 10000)
        nearest_gh = min(gh_dists) if gh_dists else "NA"
        nearest_rec = min(rec_dists) if rec_dists else "NA"

        if sbl_idx is not None:
            focal = reg_orfs[sbl_idx]
            left_orfs = reg_orfs[:sbl_idx]
            right_orfs = reg_orfs[sbl_idx + 1 :]
            left_out = left_orfs[::-1]   # ordered outward from the SBL
            left_hm = any(o["is_strict_structural_hallmark"] for o in left_orfs)
            right_hm = any(o["is_strict_structural_hallmark"] for o in right_orfs)
            left_sr = any(o["is_structural_or_replic"] for o in left_orfs)
            right_sr = any(o["is_structural_or_replic"] for o in right_orfs)
            left_gh, right_gh = any(o["hallmark"] for o in left_orfs), any(o["hallmark"] for o in right_orfs)
            left_vev, right_vev = any(o["viral_ev"] for o in left_orfs), any(o["viral_ev"] for o in right_orfs)
            any_rec = any(o["is_recomb"] for o in neigh)
            up5 = "; ".join(f"{o['orf_id'].rsplit('_',1)[-1]}({o['label']})" for o in left_orfs[-5:]) or "NONE"
            dn5 = "; ".join(f"{o['orf_id'].rsplit('_',1)[-1]}({o['label']})" for o in right_orfs[:5]) or "NONE"
        else:
            focal = None
            left_out, right_orfs = [], []
            left_hm = right_hm = left_sr = right_sr = any_rec = False
            left_gh = right_gh = left_vev = right_vev = False
            up5, dn5 = "UNMAPPED", "UNMAPPED"
        focal_cv = {1: "viral", -1: "host", 0: "no_hit"}.get(focal["cv_cat"], "unmapped") if focal else "unmapped"

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

        focal_host = focal_cv == "host" and host_basis != "genomad_provirus_skip"
        host_ok_corrected = (vg_eval >= hg_eval - (1 if focal_host else 0)) if host_basis != "genomad_provirus_skip" else True
        tier = classify_viral_origin_tier(
            boundary_ok=is_boundary_ok,
            n_hallmark=n_gh_reg,
            n_viral_ev=n_vev_reg,
            gn_present=bool(gn),
            gn_score=gn_score,
            gn_nhm=gn_nhm,
            cv_vg=vg_eval,
            cv_hg=hg_eval,
            focal_is_host_gene=focal_host,
            host_basis=host_basis,
            flanked_viral_ev_both=(left_vev and right_vev),
        )
        tier_local, n_vev_w, n_hev_w, vev_both_w = classify_tier_local(
            is_boundary_ok, n_gh_reg, bool(gn), gn_score, left_out, right_orfs)

        tier_text = classify_tier_text_rule(
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
            tier_shift_text[(tier_text, tier)] += 1
            tier_shift_local[(tier, tier_local)] += 1

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
            "viral_origin_tier_local_w10": tier_local,
            "viral_origin_tier_text_rule": tier_text,
            "viral_origin_tier_lenient_rule": tier_lenient,
            "provirus_interval_used": has_pv_interval,
            "host_gene_check_basis": host_basis,
            "host_gene_check_viral_genes": vg_eval,
            "host_gene_check_host_genes": hg_eval,
            "focal_checkv_category": focal_cv,
            "host_gene_check_pass": host_ok_corrected,
            "host_gene_check_pass_raw": host_ok,
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
            "n_genomad_hallmark_in_region": n_gh_reg,
            "n_genomad_hallmark_10kb": n_gh_10k,
            "n_viral_evidence_in_region": n_vev_reg,
            "n_host_evidence_in_region": n_hev_reg,
            "flanked_genomad_hallmark_any": left_gh or right_gh,
            "flanked_genomad_hallmark_both": left_gh and right_gh,
            "flanked_viral_evidence_any": left_vev or right_vev,
            "flanked_viral_evidence_both": left_vev and right_vev,
            "n_viral_evidence_w10": n_vev_w,
            "n_host_evidence_w10": n_hev_w,
            "flanked_viral_evidence_both_w10": vev_both_w,
            "nearest_hallmark_dist_bp": nearest_gh,
            "n_text_hallmark_in_region": n_hm_reg,
            "n_text_replic_lysis_in_region": n_rep_reg,
            "n_recomb_lysogeny_in_region": n_rec_reg,
            "flanked_text_structural_both": left_hm and right_hm,
            "flanked_text_structural_or_replic_both": left_sr and right_sr,
            "flanked_recomb_lysogeny_any": any_rec,
            "nearest_recomb_dist_bp": nearest_rec,
            "upstream_5_orfs": up5,
            "downstream_5_orfs": dn5,
        })

    out_tsv = out_dir / "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"
    fields = list(next(x for x in out_rows if x["contig_context_available"]).keys())
    with open(out_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", restval="NA")
        w.writeheader()
        w.writerows(out_rows)

    img_prim = [x for x in out_rows if x["contig_context_available"] and x["in_primary_cov65_fimo1e5"]]
    focal_img = [x for x in img_prim if x["pfam_category"].startswith("focal_")]
    n_img = max(len(img_prim), 1)

    def share(key, test=bool):
        k = sum(1 for x in img_prim if test(x[key]))
        return f"{k:,} / {len(img_prim):,} ({100 * k / n_img:.2f}%)"

    orf_ids = {o["orf_id"] for x in img_prim for o in contig_orfs.get(x["uvig"], [])}
    all_o = [orf_by_id[i] for i in orf_ids]
    n_o = max(len(all_o), 1)
    lines = [
        "=== Phase 2: viral genomic neighbourhoods and viral-origin tiers ===",
        f"Total primary curated viral loci (cov>=0.65 & FIMO<=1e-5): {sum(1 for x in out_rows if x['in_primary_cov65_fimo1e5']):,} (IMG/VR: {len(img_prim):,}, MOSAiC: {sum(1 for x in out_rows if x['in_primary_cov65_fimo1e5'] and not x['contig_context_available'])})",
        f"  Focal IMG/VR (PF00144+PF00905+PF13354): {len(focal_img):,} | Sister PRE IMG/VR (PF00768): {len(img_prim)-len(focal_img):,}",
        f"Median % ORFs annotated per viral contig (Pfam-A + eggNOG): {statistics.median(x['pct_orfs_annotated_pfam_or_eggnog'] for x in img_prim):.2f}%",
        f"Gene-origin evidence on ORFs of these contigs ({len(all_o):,} ORFs): viral {sum(o['viral_ev'] for o in all_o) / n_o:.1%}, "
        f"host only {sum(o['host_ev'] for o in all_o) / n_o:.1%}; CheckV viral {sum(o['cv_cat'] == 1 for o in all_o) / n_o:.1%}, "
        f"CheckV microbial {sum(o['cv_cat'] == -1 for o in all_o) / n_o:.1%}, geNomad virus-first marker "
        f"{sum(o['gn_marker_class'] in VIRAL_MARKER_CLASSES for o in all_o) / n_o:.1%}, eggNOG Viruses seed "
        f"{sum(o['egg_domain'] == 'Viruses' for o in all_o) / n_o:.1%}",
        f">= 1 geNomad virus hallmark (no recombination markers) in viral region: {share('n_genomad_hallmark_in_region', lambda v: v >= 1)}",
        f"geNomad hallmark on >= 1 side: {share('flanked_genomad_hallmark_any')}; on BOTH sides: {share('flanked_genomad_hallmark_both')}",
        f"Viral-evidence ORF on >= 1 side: {share('flanked_viral_evidence_any')}; on BOTH sides: {share('flanked_viral_evidence_both')}; "
        f"on both sides within {LOCAL_WINDOW_ORFS} ORFs: {share('flanked_viral_evidence_both_w10')}",
        f"Recombination / integrase / lysogeny gene in region (annotation text, MGE background only): {share('flanked_recomb_lysogeny_any')}",
        f"Median min distance to viral region boundary (IMG/VR): {statistics.median(x['min_dist_to_viral_boundary_bp'] for x in img_prim):.1f} bp",
        "\nBreakdown by (Pfam category, primary viral-origin tier) on the primary curated cohort:",
    ]
    for (cat, tr), cnt in sorted(tier_counts.items()):
        lines.append(f"  {cat:16s} | {tr:30s}: {cnt:,}")
    lines.append("\nPrimary-tier counts by family (T1 / T2 / T3): primary | local w10 | text rule")
    for fam in sorted({x["winning_acc"] for x in img_prim}):
        rs = [x for x in img_prim if x["winning_acc"] == fam]
        lines.append(f"  {fam:8s} n={len(rs):>5,}  " + "  |  ".join(
            "/".join(str(sum(1 for x in rs if x[col] == t)) for t in (TIER1, TIER2, TIER3))
            for col in ("viral_origin_tier", "viral_origin_tier_local_w10", "viral_origin_tier_text_rule")))
    basis = Counter((x["host_gene_check_basis"], x["host_gene_check_pass"], x["focal_checkv_category"]) for x in img_prim)
    lines.append("\nHost-gene check (basis, pass with focal correction, focal ORF CheckV category):")
    for (bs, ok, fc), cnt in sorted(basis.items(), key=lambda kv: str(kv[0])):
        lines.append(f"  {bs:24s} pass={ok!s:5s} focal={fc:8s}: {cnt:,}")
    tiers = [TIER1, TIER2, TIER3]
    for title, shift in (("text rule (rows) -> primary (columns)", tier_shift_text),
                         ("primary (rows) -> local w10 sensitivity (columns)", tier_shift_local),
                         ("lenient rule (rows) -> primary (columns)", tier_shift)):
        lines.append(f"\nTier shift, {title}, primary IMG/VR loci:")
        lines.append("  " + " " * 30 + "".join(f"{t.split('_')[0]:>10s}" for t in tiers))
        for old in tiers:
            lines.append(f"  {old:30s}" + "".join(f"{shift.get((old, new), 0):>10,}" for new in tiers))
    rep = "\n".join(lines) + "\n"
    (out_dir / "Phase2_Genomic_Context_Summary.txt").write_text(rep, encoding="utf-8")
    print(rep)


if __name__ == "__main__":
    main()
