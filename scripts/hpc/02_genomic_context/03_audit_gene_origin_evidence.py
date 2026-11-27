#!/usr/bin/env python3
"""Read-only audit of gene-origin evidence for the ORFs around viral SBL/PRE loci.

Evidence per ORF: eggNOG seed-ortholog domain (from eggNOG_OGs), CheckV HMM category (viral / microbial),
geNomad marker specificity class and virus hallmark flag. Compared with the annotation-text categories of
the tier rule used until 2026-11-25 (regular expressions read from 01_build_neighborhoods.py), which this
audit led to replace. Writes only to --out (default $SBL_HPC_BASE/02_genomic_context/origin_evidence_audit).
"""
import argparse
import csv
import gzip
import re
import statistics
from collections import Counter, defaultdict
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

csv.field_size_limit(2**31 - 1)
P2 = hpc.BASE / "02_genomic_context"
ap = argparse.ArgumentParser()
ap.add_argument("--nb-script", default=str(Path(__file__).resolve().parent / "01_build_neighborhoods.py"))
ap.add_argument("--phase2", default=str(P2 / "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"))
ap.add_argument("--orfs", default=str(hpc.PRODIGAL_DIR / "all_orfs.faa"))
ap.add_argument("--pfam", default=str(hpc.ORFS_PFAM_DIR / "orfs_pfam.domtbl"))
ap.add_argument("--eggnog", default=str(hpc.EGGNOG_DIR / "sbl_neighbourhood.emapper.annotations.gz"))
ap.add_argument("--genomad-genes", default=str(hpc.GENOMAD_DIR / "genomad_out" / "sbl_contigs_annotate" / "sbl_contigs_genes.tsv"))
ap.add_argument("--marker-meta", default=str(hpc.DB_DIR / "genomad_db" / "genomad_marker_metadata.tsv"))
ap.add_argument("--checkv-genes", default=str(hpc.BASE / "00_checkv" / "out" / "tmp" / "gene_features.tsv"))
ap.add_argument("--out", default=str(P2 / "origin_evidence_audit"))
a = ap.parse_args()
out = Path(a.out)
out.mkdir(parents=True, exist_ok=True)

src = Path(a.nb_script).read_text(encoding="utf-8")
ns = {}
exec(src[src.index("RE_STRICT_STRUCTURAL ="):src.index("def open_text")], {"re": re}, ns)
RE_S, RE_R, RE_REC = ns["RE_STRICT_STRUCTURAL"], ns["RE_REPLIC_LYSIS"], ns["RE_RECOMB"]


def opn(p):
    return gzip.open(p, "rt", encoding="utf-8", errors="ignore") if str(p).endswith(".gz") else open(p, encoding="utf-8", errors="ignore")


def ov(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def pct(n, d):
    return f"{n:,} ({100.0 * n / d:.1f}%)" if d else f"{n:,} (NA)"


# ---- loci
loci = []
with opn(a.phase2) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["source"] != "IMG_VR" or r["in_primary_cov65_fimo1e5"] != "True" or r["contig_context_available"] != "True":
            continue
        try:
            st, en = int(float(r["sbl_start_bp"])), int(float(r["sbl_end_bp"]))
            d5, d3 = int(float(r["dist_to_5p_viral_bp"])), int(float(r["dist_to_3p_viral_bp"]))
        except ValueError:
            continue
        if st <= 0:
            continue
        loci.append(dict(uid=r["uid"], uvig=r["uvig"], fam=r["winning_acc"], tier=r.get("viral_origin_tier_text_rule", r["viral_origin_tier"]),
                         basis=r["host_gene_check_basis"], st=st, en=en, vst=st - d5, ven=en + d3,
                         old_sr_both=r.get("flanked_text_structural_or_replic_both", r.get("flanked_viral_structural_or_replic_both")) == "True",
                         old_hm_both=r.get("flanked_text_structural_both", r.get("flanked_strict_structural_hallmark_both")) == "True"))
uvigs = {x["uvig"] for x in loci}
print(len(loci), "loci on", len(uvigs), "contigs")

# ---- ORFs
orfs = defaultdict(list)
by_id = {}
with open(a.orfs, encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line[0] != ">":
            continue
        p = line[1:].split(" # ")
        oid = p[0].strip()
        u = oid.rsplit("_", 1)[0]
        if u not in uvigs or len(p) < 4:
            continue
        o = dict(id=oid, st=int(p[1]), en=int(p[2]), pf=[], egg=None, cv=None, gn=None)
        orfs[u].append(o)
        by_id[oid] = o
for u in orfs:
    orfs[u].sort(key=lambda o: (o["st"], o["en"]))
print(len(by_id), "ORFs")

with open(a.pfam, encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line[0] == "#":
            continue
        p = line.split()
        if len(p) >= 22 and p[0] in by_id:
            by_id[p[0]]["pf"].append(f"{p[4].split('.')[0]}:{p[3]}")

DOM = (("10239|Viruses", "Viruses"), ("2|Bacteria", "Bacteria"), ("2157|Archaea", "Archaea"), ("2759|Eukaryota", "Eukaryota"))
with opn(a.eggnog) as f:
    for line in f:
        if line[0] == "#":
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[0] not in by_id:
            continue
        dom = next((d for tag, d in DOM if "@" + tag in p[4]), "root_only")
        by_id[p[0]]["egg"] = dict(dom=dom, cog=p[6], desc=p[7], pref=p[8], pfams=p[20] if len(p) > 20 else "-",
                                  seed_taxid=p[1].split(".")[0])


def map_genes(contig, gs, ge, gid=None):
    if gid and gid in by_id and by_id[gid]["st"] == gs and by_id[gid]["en"] == ge:
        return by_id[gid]
    best, bo = None, 0
    for o in orfs.get(contig, []):
        x = ov(gs, ge, o["st"], o["en"])
        if x > bo:
            best, bo = o, x
    if best is not None and bo >= 0.5 * min(ge - gs + 1, best["en"] - best["st"] + 1):
        return best
    return None


meta = {}
with opn(a.marker_meta) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        meta[r["MARKER"]] = (r["SPECIFICITY_CLASS"], r["SPM_V"])
n_gn = 0
with opn(a.genomad_genes) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        gid = r["gene"]
        contig = gid.rsplit("_", 1)[0].split("|")[0]
        if contig not in uvigs:
            continue
        o = map_genes(contig, int(r["start"]), int(r["end"]), gid)
        if o is None:
            continue
        mk = r["marker"]
        cls = meta.get(mk, ("none", ""))[0] if mk not in ("NA", "") else "none"
        o["gn"] = dict(cls=cls, vh=r["virus_hallmark"] == "1", desc=r.get("annotation_description", "NA"))
        n_gn += 1
n_cv = 0
with opn(a.checkv_genes) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["contig_id"] not in uvigs:
            continue
        o = map_genes(r["contig_id"], int(r["start"]), int(r["end"]))
        if o is not None:
            o["cv"] = int(r["hmm_cat"]) if r["hmm_cat"] not in ("NA", "") else 0
            n_cv += 1
print("mapped geNomad", n_gn, "CheckV", n_cv)


def classify(o):
    e = o["egg"]
    egg_str = "|".join(x for x in ((e["pref"], e["desc"], e["pfams"]) if e else ()) if x and x != "-")
    gdesc = o["gn"]["desc"] if o["gn"] and o["gn"]["desc"] != "NA" else ""
    txt = " ".join(x for x in ("+".join(o["pf"]), egg_str, gdesc) if x)
    rec = RE_REC.search(txt)
    s = RE_S.search(txt)
    r = RE_R.search(txt)
    gvh = bool(o["gn"] and o["gn"]["vh"])
    cls = o["gn"]["cls"] if o["gn"] else "none"
    o.update(
        txt=txt,
        rx_struct_txt=s.group(0).lower() if s else "",
        rx_struct=bool(s) or (gvh and not rec),
        rx_replic=r.group(0).lower() if r else "",
        rx_recomb=bool(rec),
        egg_dom=e["dom"] if e else "no_hit",
        cv_cat={1: "viral", -1: "microbial", 0: "no_hmm_hit"}.get(o["cv"], "not_called") if o["cv"] is not None else "not_called",
        gn_cls=cls,
        gn_grp=("virus" if cls[0] == "V" else "chromosome" if cls[0] == "C" else "plasmid" if cls[0] == "P" else "none") if cls != "none" else "none",
        gn_vh=gvh,
    )
    o["viral_ev"] = o["cv_cat"] == "viral" or o["gn_grp"] == "virus" or o["egg_dom"] == "Viruses" or gvh
    o["host_ev"] = (o["cv_cat"] == "microbial" or o["gn_grp"] == "chromosome") and not o["viral_ev"]


for o in by_id.values():
    classify(o)

# ---- neighbourhood ORFs per locus
nb_ids = set()
loc_rows = []
for L in loci:
    reg = [o for o in orfs.get(L["uvig"], []) if ov(o["st"], o["en"], L["vst"], L["ven"]) > 0]
    if not reg:
        continue
    fi = max(range(len(reg)), key=lambda i: ov(L["st"], L["en"], reg[i]["st"], reg[i]["en"]))
    left, right = reg[:fi][::-1], reg[fi + 1:]
    for o in left + right:
        nb_ids.add(o["id"])
    row = dict(uid=L["uid"], family=L["fam"], tier=L["tier"], host_basis=L["basis"], n_region_orfs=len(reg) - 1,
               old_flank_sr_both=L["old_sr_both"], old_flank_struct_both=L["old_hm_both"],
               focal_checkv=reg[fi]["cv_cat"], focal_egg_domain=reg[fi]["egg_dom"], focal_genomad=reg[fi]["gn_cls"])
    nb = left + right
    row["n_viral_ev"] = sum(o["viral_ev"] for o in nb)
    row["n_host_ev"] = sum(o["host_ev"] for o in nb)
    for k in (5, 10):
        lv, rv = sum(o["viral_ev"] for o in left[:k]), sum(o["viral_ev"] for o in right[:k])
        lh, rh = sum(o["host_ev"] for o in left[:k]), sum(o["host_ev"] for o in right[:k])
        row[f"left{k}_viral"], row[f"right{k}_viral"], row[f"left{k}_host"], row[f"right{k}_host"] = lv, rv, lh, rh
        row[f"new_flank_viral_both_{k}"] = lv >= 1 and rv >= 1
        row[f"new_flank_hallmark_both_{k}"] = any(o["gn_vh"] for o in left[:k]) and any(o["gn_vh"] for o in right[:k])
    row["new_flank_viral_both_region"] = any(o["viral_ev"] for o in left) and any(o["viral_ev"] for o in right)
    row["new_flank_hallmark_both_region"] = any(o["gn_vh"] for o in left) and any(o["gn_vh"] for o in right)
    loc_rows.append(row)

with open(out / "Locus_Origin_Evidence.tsv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(loc_rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(loc_rows)
NB = [by_id[i] for i in sorted(nb_ids)]
with open(out / "Neighbourhood_ORF_Origin_Evidence.tsv", "w", newline="", encoding="utf-8") as f:
    cols = ["id", "egg_dom", "cv_cat", "gn_cls", "gn_vh", "viral_ev", "host_ev", "rx_struct", "rx_struct_txt", "rx_replic", "rx_recomb", "txt"]
    w = csv.writer(f, delimiter="\t")
    w.writerow(cols)
    for o in NB:
        w.writerow([o[c] for c in cols])

# ---- summary
S = []
n = len(NB)
S.append(f"Loci: {len(loc_rows):,} primary IMG/VR loci; neighbourhood ORFs (viral region, focal excluded, unique): {n:,}")
for key, title in (("egg_dom", "eggNOG seed-ortholog domain"), ("cv_cat", "CheckV HMM category"), ("gn_cls", "geNomad marker specificity class")):
    c = Counter(o[key] for o in NB)
    S.append(f"\n--- {title} ---")
    for k, v in c.most_common():
        S.append(f"  {k:14s} {pct(v, n)}")
S.append(f"\ngeNomad virus hallmark: {pct(sum(o['gn_vh'] for o in NB), n)}")
S.append(f"Any viral evidence (CheckV viral | geNomad V* marker | eggNOG Viruses | geNomad hallmark): {pct(sum(o['viral_ev'] for o in NB), n)}")
S.append(f"Host evidence only (CheckV microbial | geNomad C* marker, no viral evidence): {pct(sum(o['host_ev'] for o in NB), n)}")
S.append(f"Neither: {pct(sum(not o['viral_ev'] and not o['host_ev'] for o in NB), n)}")


def xtab(rows, k1, k2, title):
    S.append(f"\n--- {title} (rows {k1}, cols {k2}) ---")
    c = Counter((o[k1], o[k2]) for o in rows)
    cs = sorted({o[k2] for o in rows})
    S.append("  " + "".join(f"{x:>14s}" for x in [""] + cs))
    for r_ in sorted({o[k1] for o in rows}):
        S.append("  " + f"{r_:>14s}" + "".join(f"{c[(r_, x)]:>14,}" for x in cs))


xtab(NB, "cv_cat", "gn_grp", "CheckV vs geNomad marker group")
xtab(NB, "cv_cat", "egg_dom", "CheckV vs eggNOG domain")

for flag, title in (("rx_replic", "regex replication/lysis"), ("rx_struct_txt", "regex structural (text only)")):
    hit = [o for o in NB if o[flag]]
    S.append(f"\n--- ORFs flagged by {title}: {len(hit):,} ---")
    S.append(f"  with viral evidence: {pct(sum(o['viral_ev'] for o in hit), len(hit))}; host evidence only: {pct(sum(o['host_ev'] for o in hit), len(hit))}")
    tok = defaultdict(list)
    for o in hit:
        tok[o[flag]].append(o)
    S.append(f"  {'matched token':24s}{'n':>8s}{'viral_ev':>10s}{'host_only':>10s}{'egg_Bact':>10s}")
    for t, lst in sorted(tok.items(), key=lambda x: -len(x[1]))[:30]:
        S.append(f"  {t[:24]:24s}{len(lst):>8,}{sum(o['viral_ev'] for o in lst):>10,}{sum(o['host_ev'] for o in lst):>10,}"
                 f"{sum(o['egg_dom'] == 'Bacteria' for o in lst):>10,}")

S.append("\n--- Locus level: old regex flank vs new evidence-based flank ---")
for basis in sorted({r["host_basis"] for r in loc_rows}) + ["ALL"]:
    rs = [r for r in loc_rows if basis == "ALL" or r["host_basis"] == basis]
    S.append(f"  [{basis}] n={len(rs):,}")
    S.append(f"    old flanked_structural_or_replic_both: {pct(sum(r['old_flank_sr_both'] for r in rs), len(rs))}")
    for k in ("5", "10", "region"):
        key = f"new_flank_viral_both_{k}"
        S.append(f"    new viral-evidence both sides (window {k}): {pct(sum(r[key] for r in rs), len(rs))}; "
                 f"old pass & new fail: {sum(r['old_flank_sr_both'] and not r[key] for r in rs):,}")
    S.append(f"    old flanked_strict_structural_both: {pct(sum(r['old_flank_struct_both'] for r in rs), len(rs))}; "
             f"geNomad hallmark both sides (region): {pct(sum(r['new_flank_hallmark_both_region'] for r in rs), len(rs))}")

S.append("\n--- Locus level by text-rule tier: median viral / host evidence ORFs in region; focal ORF CheckV call ---")
for t in sorted({r["tier"] for r in loc_rows}):
    rs = [r for r in loc_rows if r["tier"] == t]
    fc = Counter(r["focal_checkv"] for r in rs)
    S.append(f"  {t:30s} n={len(rs):>5,}  viral_ev median {statistics.median(r['n_viral_ev'] for r in rs):>5}  "
             f"host_ev median {statistics.median(r['n_host_ev'] for r in rs):>5}  "
             f"both-side viral (10 ORFs) {pct(sum(r['new_flank_viral_both_10'] for r in rs), len(rs))}  focal CheckV {dict(fc)}")
S.append("\n--- Focal ORF eggNOG domain by family ---")
for fam in sorted({r["family"] for r in loc_rows}):
    rs = [r for r in loc_rows if r["family"] == fam]
    S.append(f"  {fam}: {dict(Counter(r['focal_egg_domain'] for r in rs))}")

txt = "\n".join(S)
(out / "Protein_Origin_Audit_Summary.txt").write_text(txt + "\n", encoding="utf-8")
print(txt)
