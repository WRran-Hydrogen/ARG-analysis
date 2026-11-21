#!/usr/bin/env python3
"""Phase 2: 100%-Coverage Viral Genomic Neighborhood & Confidence Tier Integration.

Integrates for all curated IMG/VR viral SBL loci:
- Prodigal meta ORFs (51cc034a/all_orfs.faa)
- Full-genome Pfam-A (--cut_ga) annotations (5113374/orfs_pfam.domtbl, 201,317 hits)
- Full-genome eggNOG-mapper annotations (5113454/sbl_neighbourhood.emapper.annotations.gz, 103,110 ORFs)
- geNomad v1.7 summary, hallmarks, and provirus boundaries (5113378/6ca11f3e)
- CheckV v1.0.3 (DB v1.5) quality and viral/host region coordinates (5113663)
"""

import csv
import gzip
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import hpc_config as hpc

BASE = hpc.BASE
P1_OUT = BASE / "01_final_curated"
P1_DEC = BASE / "01_curated_datasets" / "Phase1_Viral_Decontamination_Audit.tsv"
CHECKV_DIR = BASE / "00_checkv" / "out"
GENOMAD_DIR = hpc.GENOMAD_DIR
PRODIGAL_ORFS = hpc.PRODIGAL_DIR / "all_orfs.faa"
ORFS_PFAM = hpc.ORFS_PFAM_DIR / "orfs_pfam.domtbl"
EGGNOG_GZ = hpc.EGGNOG_DIR / "sbl_neighbourhood.emapper.annotations.gz"

OUT_DIR = BASE / "02_genomic_context"

RE_HALLMARK = re.compile(
    r"(terminase|portal|capsid|head[- _]protein|tail|baseplate|tape[- _]?measure|virion|"
    r"phage[- _]?(?:coat|sheath|neck|fiber|spike|tube|wedge|hub|lysozyme|structural)|"
    r"prohead|scaffolding|holin|endolysin|muramidase|DTR|MCP|major[- _]?capsid|"
    r"gp\d+.*(?:tail|head|capsid|portal))",
    re.I,
)
RE_REPLIC = re.compile(
    r"(dna[- _]?polymerase|rna[- _]?polymerase|helicase|primase|exonuclease|endonuclease|"
    r"dutpase|ribonucleotide[- _]?reductase|thymidylate|single[- _]?strand.*binding|"
    r"sigma[- _]?70|anti[- _]?repressor|lysis|spanin|replication|ssb|recT|recE|erf)",
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


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load Final Viral SBL Catalog & Decontamination coordinates
    print("[1/7] Loading final curated viral SBL catalog...")
    decontam = {}
    with open(P1_DEC, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            decontam[r["sbl_protein_id"]] = r

    viral_rows = []
    with open(P1_OUT / "Final_Viral_SBL_Catalog.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cat = r["pfam_category"]
            if cat.startswith("focal_") and r["pass_fimo_1e5"] != "True":
                continue
            if not (cat.startswith("focal_") or cat.startswith("aux_")):
                continue
            viral_rows.append(r)

    target_uvigs = {decontam[r["orig_uid"]]["uvig"] for r in viral_rows if r["source"] == "IMG_VR" and r["orig_uid"] in decontam}
    print(f"  Final viral loci: {len(viral_rows)} across {len(target_uvigs)} IMG/VR contigs (+ MOSAiC={sum(1 for r in viral_rows if r['source']=='MOSAiC')}).")

    # 2. Load Prodigal ORFs for target contigs
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

    # Quick lookup maps
    orf_by_id = {}
    for uvig, lst in contig_orfs.items():
        for o in lst:
            orf_by_id[o["orf_id"]] = o

    # 3. Parse full-genome Pfam-A annotations (5113374)
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

    # 4. Parse full-genome eggNOG-mapper annotations (5113454)
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

    # 5. Parse geNomad genes, virus summary, and provirus boundaries (5113378)
    print("[5/7] Loading geNomad v1.7 gene hallmarks & summaries...")
    gn_vs = {}
    with open_text(GENOMAD_DIR / "sbl_contigs_virus_summary.tsv.gz") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            uv = r["seq_name"].split("|")[0]
            if uv in target_uvigs:
                # Keep best virus_score entry if contig has multiple proviruses
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
                # Fallback by coordinate overlap >= 50%
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

    # 6. Load CheckV v1.0.3 quality_summary & contamination boundaries
    print("[6/7] Loading CheckV v1.0.3 quality & viral boundaries...")
    cv_qs = {}
    with open(CHECKV_DIR / "quality_summary.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["contig_id"] in target_uvigs:
                cv_qs[r["contig_id"]] = r

    cv_v_regions = {}
    with open(CHECKV_DIR / "contamination.tsv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = r["contig_id"]
            if cid in target_uvigs and r.get("provirus", "No") == "Yes":
                ivs = parse_checkv_viral_intervals(r.get("region_types", ""), r.get("region_coords_bp", ""))
                if ivs:
                    cv_v_regions[cid] = ivs

    # Pre-classify every ORF's functional roles on all target contigs
    for uvig, orfs in contig_orfs.items():
        for o in orfs:
            # Resolve non-overlapping top Pfams for readable label
            pf_sorted = sorted(o["pfams"], key=lambda x: (-x[2], x[0]))
            kept_pf = []
            for st, en, sc, label in pf_sorted:
                if not any(interval_overlap(st, en, kst, ken) >= 20 for kst, ken, _ in kept_pf):
                    kept_pf.append((st, en, label))
            pf_str = "+".join(x[2] for x in sorted(kept_pf))
            egg = o["eggnog"]
            egg_str = ""
            if egg:
                parts = [x for x in (egg[1], egg[2], egg[3]) if x and x != "NA"]
                egg_str = "|".join(parts[:2])
            combined_txt = " ".join(x for x in (pf_str, egg_str, o["gn_desc"]) if x)
            o["label"] = pf_str or (egg[1] if egg and egg[1] else "") or (o["gn_desc"][:40] if o["gn_desc"] else "hypothetical")
            o["is_hallmark"] = o["gn_hallmark"] or bool(RE_HALLMARK.search(combined_txt))
            o["is_replic"] = bool(RE_REPLIC.search(combined_txt))
            o["is_recomb"] = bool(RE_RECOMB.search(combined_txt))
            o["is_viral_or_recomb"] = o["is_hallmark"] or o["is_replic"] or o["is_recomb"]
            o["has_any_annot"] = bool(combined_txt.strip())

    # 7. Assemble per-locus genomic neighborhood report & Confidence Tiers
    print("[7/7] Computing locus-level neighborhood features & confidence tiers...")
    out_rows = []
    tier_counts = Counter()

    for r in viral_rows:
        uid = r["uid"]
        orig_uid = r["orig_uid"]
        src = r["source"]
        fam = r["winning_acc"]
        cat = r["pfam_category"]

        if src != "IMG_VR" or orig_uid not in decontam:
            tier = "MOSAiC_Proteome_Verified_vMAG"
            tier_counts[(cat, tier)] += 1
            out_rows.append({
                "uid": uid, "orig_uid": orig_uid, "source": src, "uvig": "NA",
                "pfam_category": cat, "winning_acc": fam,
                "dual_PF13354_PF00768_GA": r["dual_PF13354_PF00768_GA"],
                "contig_context_available": False, "viral_origin_tier": tier,
                "contig_length_bp": "NA", "total_orfs_on_contig": "NA",
                "sbl_start_bp": "NA", "sbl_end_bp": "NA", "sbl_strand": "NA",
                "dist_to_5p_viral_bp": "NA", "dist_to_3p_viral_bp": "NA", "min_dist_to_viral_boundary_bp": "NA",
                "genomad_virus_score": "NA", "genomad_topology": "NA", "genomad_n_hallmarks": "NA",
                "genomad_taxonomy": r["vir_sub_info"], "checkv_quality": "NA", "checkv_miuvig_quality": "NA",
                "checkv_viral_genes": "NA", "checkv_host_genes": "NA", "checkv_contamination_pct": "NA",
                "n_orfs_in_viral_region": "NA", "pct_orfs_annotated_pfam_or_eggnog": "NA",
                "n_viral_hallmark_in_region": "NA", "n_viral_hallmark_10kb": "NA",
                "n_viral_replic_lysis_in_region": "NA", "n_recomb_lysogeny_in_region": "NA",
                "flanked_both_sides_viral_or_recomb": "NA", "flanked_any_viral_or_recomb": "NA",
                "nearest_hallmark_dist_bp": "NA", "nearest_recomb_dist_bp": "NA",
                "upstream_5_orfs": "NA", "downstream_5_orfs": "NA",
            })
            continue

        dec = decontam[orig_uid]
        uvig = dec["uvig"]
        st = int(dec["orf_start"]) if dec["orf_start"] != "NA" else 0
        en = int(dec["orf_end"]) if dec["orf_end"] != "NA" else 0

        cq = cv_qs.get(uvig, {})
        gn = gn_vs.get(uvig, {})
        orfs = contig_orfs.get(uvig, [])

        clen = int(cq.get("contig_length", 0) or 0)
        if clen == 0 and orfs:
            clen = max(o["end"] for o in orfs)

        # Effective viral interval on contig containing SBL
        v_st, v_en = 1, clen
        if uvig in gn_pv and st > 0:
            cands = [(ps, pe) for ps, pe in gn_pv[uvig] if interval_overlap(st, en, ps, pe) > 0]
            if cands:
                best_gn = max(cands, key=lambda x: interval_overlap(st, en, x[0], x[1]))
                v_st = max(v_st, best_gn[0])
                v_en = min(v_en, best_gn[1])
        if uvig in cv_v_regions and st > 0:
            cands = [(ps, pe) for ps, pe in cv_v_regions[uvig] if interval_overlap(st, en, ps, pe) > 0]
            if cands:
                best_cv = max(cands, key=lambda x: interval_overlap(st, en, x[0], x[1]))
                v_st = max(v_st, best_cv[0])
                v_en = min(v_en, best_cv[1])
        if v_st > v_en:
            v_st, v_en = 1, clen

        dist_5p = max(0, st - v_st) if st > 0 else 0
        dist_3p = max(0, v_en - en) if en > 0 else 0
        min_dist = min(dist_5p, dist_3p) if st > 0 else 0

        # ORFs within the viral region
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

        # Exclude SBL itself when counting neighboring hallmarks/replic/recomb
        neigh = [o for idx, o in enumerate(reg_orfs) if idx != sbl_idx]
        n_hm_reg = sum(1 for o in neigh if o["is_hallmark"])
        n_rep_reg = sum(1 for o in neigh if o["is_replic"])
        n_rec_reg = sum(1 for o in neigh if o["is_recomb"])

        def orf_dist(o):
            if st == 0:
                return 9999999
            if o["end"] < st:
                return st - o["end"]
            if o["start"] > en:
                return o["start"] - en
            return 0

        hm_dists = [orf_dist(o) for o in neigh if o["is_hallmark"]]
        rec_dists = [orf_dist(o) for o in neigh if o["is_recomb"]]
        n_hm_10k = sum(1 for d in hm_dists if d <= 10000)
        nearest_hm = min(hm_dists) if hm_dists else "NA"
        nearest_rec = min(rec_dists) if rec_dists else "NA"

        if sbl_idx is not None:
            left_orfs = reg_orfs[:sbl_idx]
            right_orfs = reg_orfs[sbl_idx + 1 :]
            left_vr = any(o["is_viral_or_recomb"] for o in left_orfs)
            right_vr = any(o["is_viral_or_recomb"] for o in right_orfs)
            up5 = "; ".join(f"{o['orf_id'].rsplit('_',1)[-1]}({o['label']})" for o in left_orfs[-5:]) or "NONE"
            dn5 = "; ".join(f"{o['orf_id'].rsplit('_',1)[-1]}({o['label']})" for o in right_orfs[:5]) or "NONE"
        else:
            left_vr, right_vr = False, False
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

        # Upgraded Tier 1 / Tier 2 / Tier 3 classification
        is_complete_topo = gn_topo in ("DTR", "ITR") or cv_qual in ("Complete", "High-quality", "Medium-quality")
        has_strong_hm = (n_hm_reg >= 2) or (gn_nhm >= 2)
        if (gn_score >= 0.90 or is_complete_topo) and has_strong_hm and (min_dist >= 2000 or gn_topo in ("DTR", "ITR")):
            tier = "Tier1_HighConfidence_Viral"
        elif gn and (n_hm_reg >= 1 or gn_nhm >= 1 or cv_vg >= 1 or n_rep_reg >= 1):
            tier = "Tier2_Supported_Viral"
        else:
            tier = "Tier3_Fragment_or_Borderline"

        tier_counts[(cat, tier)] += 1
        out_rows.append({
            "uid": uid,
            "orig_uid": orig_uid,
            "source": src,
            "uvig": uvig,
            "pfam_category": cat,
            "winning_acc": fam,
            "dual_PF13354_PF00768_GA": r["dual_PF13354_PF00768_GA"],
            "contig_context_available": True,
            "viral_origin_tier": tier,
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
            "flanked_both_sides_viral_or_recomb": left_vr and right_vr,
            "flanked_any_viral_or_recomb": left_vr or right_vr,
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

    # Summary statistics
    img_rows = [x for x in out_rows if x["contig_context_available"]]
    focal_img = [x for x in img_rows if x["pfam_category"].startswith("focal_")]
    lines = [
        "=== Phase 2: 100%-Coverage Viral Genomic Neighborhood Summary ===",
        f"Total final curated viral SBL loci: {len(out_rows):,} (IMG/VR w/ full contig ORFs: {len(img_rows):,}, MOSAiC: {len(out_rows)-len(img_rows)})",
        f"Focal IMG/VR loci (PF00144+PF00905+PF13354): {len(focal_img):,}",
        f"Median % ORFs annotated per viral contig (Pfam-A + eggNOG): {statistics.median(x['pct_orfs_annotated_pfam_or_eggnog'] for x in img_rows):.2f}%",
        f"Contigs with >=1 viral hallmark gene inside viral region: {sum(1 for x in img_rows if x['n_viral_hallmark_in_region'] >= 1):,} / {len(img_rows):,} ({100*sum(1 for x in img_rows if x['n_viral_hallmark_in_region'] >= 1)/max(len(img_rows),1):.2f}%)",
        f"Flanked on at least one side by viral/recomb gene: {sum(1 for x in img_rows if x['flanked_any_viral_or_recomb']):,} / {len(img_rows):,} ({100*sum(1 for x in img_rows if x['flanked_any_viral_or_recomb'])/max(len(img_rows),1):.2f}%)",
        f"Flanked on BOTH sides by viral/recomb genes: {sum(1 for x in img_rows if x['flanked_both_sides_viral_or_recomb']):,} / {len(img_rows):,} ({100*sum(1 for x in img_rows if x['flanked_both_sides_viral_or_recomb'])/max(len(img_rows),1):.2f}%)",
        f"Median min distance to viral region boundary (IMG/VR): {statistics.median(x['min_dist_to_viral_boundary_bp'] for x in img_rows):.1f} bp",
        "\nBreakdown by (Pfam Category, Viral Origin Confidence Tier):",
    ]
    for (cat, tr), cnt in sorted(tier_counts.items()):
        lines.append(f"  {cat:16s} | {tr:30s}: {cnt:,}")
    rep = "\n".join(lines) + "\n"
    (OUT_DIR / "Phase2_Genomic_Context_Summary.txt").write_text(rep, encoding="utf-8")
    print(rep)


if __name__ == "__main__":
    main()
