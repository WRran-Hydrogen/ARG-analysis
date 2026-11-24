#!/usr/bin/env python3
"""Re-stratify viral & bacterial CL0013 enzymatic functions in the local curated catalogs.

Applies the biochemical rules of scripts/hpc/01_curation/09_build_curated_datasets.py to the local
catalogs (Final_Viral_SBL_Catalog.tsv, Final_Bacterial_TripleSource_Catalog.tsv) and writes
Viral_Enzymatic_Function_Verification_2026.tsv. The HPC catalogs (01_final_curated/) already carry
these tiers, so the HGT run uses them directly.

  python reclassify_enzymatic_function_2026.py              # catalogs + verification table
  python reclassify_enzymatic_function_2026.py --sync-hgt   # also rewrite the functional columns
                                                            # of HGT tables made before the fix
--sync-hgt regenerates HGT_Summary_By_Family_And_Function_2026.tsv and Stage3_HGT_Summary_Report.txt
in its own (older) format, so do not use it on output of the current 04_hgt_inference/01 script.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Auto-dispatch under WSL Python on Windows so /usr/bin/hmmalign and pandas/numpy run natively
if sys.platform.startswith("win"):
    _in_wsl = False
else:
    _in_wsl = True

if not _in_wsl:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import sbl_config as _cfg

    cmd = _cfg.linux_cmd(_cfg.WSL_PY, _cfg.to_wsl_path(Path(__file__).resolve()), *sys.argv[1:])
    proc = subprocess.run(cmd)
    sys.exit(proc.returncode)

import pandas as pd

csv.field_size_limit(10_000_000)

ROOT = Path(__file__).resolve().parents[2]
CUR_DIR = ROOT / "data" / "curated_2026"
HMM_DIR = ROOT / "data" / "hmm_models"
RES_DIR = ROOT / "results" / "reanalysis_2026"
C1_DIR = RES_DIR / "01_curation_and_context"
HGT_DIR = RES_DIR / "04_hgt_inference"

HMMALIGN_BIN = shutil.which("hmmalign") or "/usr/bin/hmmalign"

CL0013_HMM_LENS = {
    "PF00144": 329,
    "PF00905": 300,
    "PF13354": 256,
    "PF00768": 239,
    "PF02113": 445,
}

MOTIF_WINDOWS = {
    "PF00144": {"m1": (49, 62), "m2": (137, 154), "omega_e": None, "omega_w": None, "m3": (269, 300)},
    "PF13354": {"m1": (37, 49), "m2": (98, 110), "omega_e": (126, 149), "omega_w": None, "m3": (198, 214)},
    "PF00905": {"m1": (37, 52), "m2": (119, 140), "omega_e": None, "omega_w": (119, 175), "m3": (229, 252)},
}

CARD_PBP_RE = re.compile(
    r"(penicillin[-_ ]?binding|\bPBP\d*[a-z]?\b|\bpbp\d*[a-z]?\b|"
    r"\bmec[ACRI]\d*\b|\bbla[RI]\d*\b|\bamp[RDG]\b|\bpenA\b|ftsI|"
    r"gene[-_ ]modulating[-_ ]beta-lactam[-_ ]resistance)",
    re.I,
)
OMEGA_E166_RE = re.compile(r"E[A-Z]{2}[LIVMFCA]N|TE[A-Z]{2}[MLIVFA][SNHRKAE]")
CLASSD_TRPCAP_RE = re.compile(r"[FYW]WL|W[A-Z][EDQNG][A-Z]{2}[LIVMF]")
NON_SBL_SENSOR_RE = re.compile(r"(PF05569|\bblaR\d*\b|\bmecR\d*\b|modulating beta-lactam resistance)", re.I)
PHAGE_FUSION_DOMS = (
    "TMP_3", "Endolysin", "SLT", "VanY", "Peptidase_M23", "Amidase",
    "SH3", "Phage_spike", "Phage_tail", "SLH", "Gly_rich_2", "Choline_bind",
)
PBP_SYNTH_DOMS = (
    "Transgly", "PASTA", "PCB_OB", "PBP_middle",
    "PBP_1A_C", "PBP_dimer", "BiPBP_C", "PBP5_C",
)


def as_bool(val: object) -> bool:
    return str(val).strip().lower() in ("true", "t", "1", "yes")


def parse_fasta_str(text: str) -> dict[str, str]:
    """Parse FASTA records from a string into {seq_id: sequence}."""
    seqs: dict[str, str] = {}
    hdr: str | None = None
    buf: list[str] = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            if hdr is not None:
                seqs[hdr] = "".join(buf)
            hdr = s[1:].split()[0]
            buf = []
        else:
            buf.append(s)
    if hdr is not None:
        seqs[hdr] = "".join(buf)
    return seqs


def split_a2m_match_states(a2m_row: str) -> tuple[list[str], list[list[str]]]:
    """Split an A2M row into HMM match states and following lowercase insertions."""
    states: list[str] = []
    inserts: list[list[str]] = []
    for ch in a2m_row:
        if ch == ".":
            continue
        if ch.isupper() or ch == "-":
            states.append(ch)
            inserts.append([])
        elif ch.islower() and inserts:
            inserts[-1].append(ch.upper())
    return states, inserts


def load_domain_sequences() -> dict[str, dict[str, str]]:
    """Load all local domain FASTAs per Pfam family (preferring unclustered *_dom.faa)."""
    by_fam: dict[str, dict[str, str]] = defaultdict(dict)
    for fam in ("PF00144", "PF00905", "PF13354", "PF00768"):
        candidates = [
            CUR_DIR / f"{fam}_viral_dom.faa",
            CUR_DIR / f"{fam}_bacterial_nr85_dom.faa",
            CUR_DIR / f"{fam}_combined_tree_cov60_dom.faa",
            CUR_DIR / f"{fam}_combined_tree_dom.faa",
            CUR_DIR / f"{fam}_combined_tree_phylo_dom.faa",
        ]
        for p in candidates:
            if p.is_file():
                parsed = parse_fasta_str(p.read_text(encoding="utf-8", errors="ignore"))
                for uid, seq in parsed.items():
                    by_fam[fam].setdefault(uid, seq)
    return by_fam


def compute_a2m_motif_flags(fam: str, seqs: dict[str, str]) -> dict[str, dict]:
    """Run hmmalign --outformat A2M against official Pfam HMM and extract motif flags."""
    if not seqs or fam not in MOTIF_WINDOWS:
        return {}
    hmm_path = HMM_DIR / f"{fam}_official.hmm"
    if not hmm_path.is_file():
        raise FileNotFoundError(f"Missing official HMM: {hmm_path}")

    win_cfg = MOTIF_WINDOWS[fam]
    hmm_len = CL0013_HMM_LENS[fam]

    # Map to compact safe IDs so hmmalign never truncates long IMG/VR headers
    sid_to_uid: dict[str, str] = {}
    fa_chunks: list[str] = []
    for idx, (uid, seq) in enumerate(seqs.items()):
        sid = f"S_{idx:06d}"
        sid_to_uid[sid] = uid
        fa_chunks.append(f">{sid}\n{seq}\n")

    proc = subprocess.run(
        [HMMALIGN_BIN, "--outformat", "A2M", str(hmm_path), "-"],
        input="".join(fa_chunks),
        text=True,
        capture_output=True,
        check=True,
    )
    a2m_map = parse_fasta_str(proc.stdout)
    out: dict[str, dict] = {}

    for sid, row in a2m_map.items():
        uid = sid_to_uid[sid]
        states, inserts = split_a2m_match_states(row)
        if len(states) != hmm_len:
            raise RuntimeError(f"{fam}: {uid} yielded {len(states)} match states != LENG={hmm_len}")
        mseq = "".join(states)

        s1, e1 = win_cfg["m1"]
        seg1 = mseq[s1 - 1 : e1].replace("-", "")
        has_m1 = bool(re.search(r"S[A-Z]{2}K", seg1))
        has_m1_clsD = bool(re.search(r"S[TFASV][FYWLIVM]K", seg1)) if fam == "PF00905" else False

        raw_om_w = ""
        has_trp_cap = False
        if win_cfg.get("omega_w") is not None:
            ww_lo, ww_hi = win_cfg["omega_w"]
            raw_om_w = "".join(
                (states[k - 1] if states[k - 1] != "-" else "") + "".join(inserts[k - 1])
                for k in range(max(1, ww_lo), min(hmm_len, ww_hi) + 1)
            )
            has_trp_cap = bool(CLASSD_TRPCAP_RE.search(raw_om_w))

        s2, e2 = win_cfg["m2"]
        seg2 = mseq[s2 - 1 : e2].replace("-", "")
        if fam == "PF00905":
            if re.search(r"[YF]GN", raw_om_w):
                m2_cls = "YGN_ClassD"
            elif re.search(r"[SY][A-Z]N", seg2):
                m2_cls = "SXN_PBP"
            else:
                m2_cls = "NONE"
        elif fam == "PF13354":
            if re.search(r"[SY][DN]N", seg2):
                m2_cls = "SDN_ClassA"
            elif re.search(r"[SY][A-Z]N", seg2):
                m2_cls = "SXN_PBP"
            else:
                m2_cls = "NONE"
        else:
            if re.search(r"Y[AS]N", seg2):
                m2_cls = "YAN_ClassC"
            elif re.search(r"[SY][A-Z]N", seg2):
                m2_cls = "SXN_PBP"
            else:
                m2_cls = "NONE"

        s3, e3 = win_cfg["m3"]
        seg3 = mseq[s3 - 1 : e3].replace("-", "")
        if fam == "PF13354":
            has_m3 = bool(re.search(r"[KR][TS][GA]", seg3))
        else:
            has_m3 = bool(re.search(r"K[TS]G", seg3))

        has_e166 = False
        if win_cfg.get("omega_e") is not None:
            w_lo, w_hi = win_cfg["omega_e"]
            raw_om = "".join(
                (states[k - 1] if states[k - 1] != "-" else "") + "".join(inserts[k - 1])
                for k in range(max(1, w_lo), min(hmm_len, w_hi) + 1)
            )
            has_e166 = bool(OMEGA_E166_RE.search(raw_om))

        out[uid] = {
            "motif1_SXXK_hmm": has_m1,
            "has_m1_clsD": has_m1_clsD,
            "motif2_class": m2_cls,
            "motif3_KTG_hmm": has_m3,
            "has_omega_E166": has_e166,
            "has_trp_cap": has_trp_cap,
            "_triad_aligned": True,
        }
    return out


def classify_row_enzymatic(row: dict, flags: dict | None) -> dict:
    """Apply corrected four-tier enzymatic rules to a viral or bacterial catalog row."""
    fam = row.get("winning_acc", "")
    if fam not in ("PF00144", "PF00905", "PF13354", "PF00768", "PF02113"):
        row["enzymatic_validation_class"] = "CL0013_NonPRE_Outgroup"
        return row

    if flags is not None:
        m1 = bool(flags["motif1_SXXK_hmm"])
        m1_d = bool(flags["has_m1_clsD"])
        m2 = str(flags["motif2_class"])
        m3 = bool(flags["motif3_KTG_hmm"])
        om_e = bool(flags["has_omega_E166"])
        trp_cap = bool(flags["has_trp_cap"])
        triad_aligned = bool(flags.get("_triad_aligned", True))
        row["motif1_SXXK_hmm"] = str(m1)
        row["motif2_class"] = m2
        row["motif3_KTG_hmm"] = str(m3)
        row["has_omega_E166"] = str(om_e)
    else:
        # Unaligned fragments / PF00768 / PF02113: disable HTG false-positives & mark unaligned
        m1 = as_bool(row.get("motif1_SXXK_hmm", False))
        m1_d = False
        m2 = str(row.get("motif2_class", "NONE"))
        m3 = False if fam in MOTIF_WINDOWS else as_bool(row.get("motif3_KTG_hmm", False))
        om_e = as_bool(row.get("has_omega_E166", False))
        trp_cap = False
        triad_aligned = False
        if fam in MOTIF_WINDOWS:
            row["motif3_KTG_hmm"] = "False"

    arch = str(row.get("domain_architecture", ""))
    card_hit = str(row.get("best_card_hit", "NONE"))
    card_fam = str(row.get("card_gene_family", "NONE"))
    card_text = f"{card_hit} {card_fam}"
    pident = float(row.get("card_pident") or 0.0)
    cov = float(row.get("card_cov") or 0.0)

    is_sensor = bool(NON_SBL_SENSOR_RE.search(f"{arch} {card_text}"))
    has_phage_fusion = any(d in arch for d in PHAGE_FUSION_DOMS)
    has_pbp_fusion = any(d in arch for d in PBP_SYNTH_DOMS)
    card_ok = (pident >= 40.0) and (cov >= 70.0) and (not is_sensor) and (not bool(CARD_PBP_RE.search(card_text)))

    if not triad_aligned:
        triad_intact = False
    elif fam == "PF13354":
        triad_intact = m1 and (m2 in ("SDN_ClassA", "SXN_PBP")) and m3
    elif fam == "PF00144":
        triad_intact = m1 and (m2 == "YAN_ClassC") and m3
    elif fam == "PF00905":
        triad_intact = m1 and m3
    else:
        triad_intact = False

    is_cls_a = (fam == "PF13354") and triad_aligned and (not is_sensor) and (not has_phage_fusion) and m1 and (m2 == "SDN_ClassA") and om_e and m3
    is_cls_c = (fam == "PF00144") and triad_aligned and (not is_sensor) and (not has_phage_fusion) and m1 and (m2 == "YAN_ClassC") and m3
    is_cls_d = (fam == "PF00905") and triad_aligned and (not is_sensor) and (not has_phage_fusion) and (not has_pbp_fusion) and m1_d and (m2 == "YGN_ClassD") and trp_cap and m3

    if (card_ok and triad_intact) or is_cls_a or is_cls_c or is_cls_d:
        row["functional_tier"] = "Putative_Hydrolytic_BetaLactamase"
        row["enzymatic_validation_class"] = "Verified_Hydrolytic_SBL" if card_ok else "Candidate_Divergent_SBL"
        ev = []
        if card_ok:
            ev.append("CARD_Match")
        if is_cls_a:
            ev.append("ClassA_SXXK+SDN+E166+KTG")
        if is_cls_c:
            ev.append("ClassC_SXXK+YAN+KTG_Lys315")
        if is_cls_d:
            ev.append("ClassD_OXA_SXXK+YFGN+TrpCap+KTG")
        row["functional_evidence"] = ";".join(ev)
    elif card_ok and not triad_intact:
        row["functional_tier"] = "PBP_or_PRE_Homolog"
        row["enzymatic_validation_class"] = "Truncated_SBL_Pseudogene"
        row["functional_evidence"] = "CARD_Homolog_Missing_Catalytic_Motif"
    else:
        row["functional_tier"] = "PBP_or_PRE_Homolog"
        row["enzymatic_validation_class"] = "BlaR1_Sensor_Receptor" if is_sensor else "PBP_Esterase_or_PhageLysin_Homolog"
        row["functional_evidence"] = "BlaR1_MecR1_Sensor" if is_sensor else f"PBP_PRE_{fam}"

    return row


def update_catalog(tsv_path: Path, flags_by_fam: dict[str, dict[str, dict]], is_bacterial: bool = False) -> list[dict]:
    """Read catalog TSV, re-evaluate enzymatic classification, and write back preserving columns."""
    with open(tsv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        orig_fields = list(reader.fieldnames or [])
        rows = list(reader)

    # Inherit nr85 representative motif flags for redundant non-rep bacterial members
    cluster_flags: dict[str, dict] = {}
    if is_bacterial:
        for r in rows:
            uid = r["uid"]
            fam = r.get("winning_acc", "")
            cid = r.get("mmseqs85_cluster_id", "")
            if cid and cid != "NA" and uid in flags_by_fam.get(fam, {}):
                cluster_flags[cid] = flags_by_fam[fam][uid]

    updated: list[dict] = []
    for r in rows:
        fam = r.get("winning_acc", "")
        uid = r["uid"]
        fl = flags_by_fam.get(fam, {}).get(uid)
        if fl is None and is_bacterial:
            cid = r.get("mmseqs85_cluster_id", "")
            if cid and cid != "NA" and cid in cluster_flags and as_bool(r.get("pass_cov65", False)) and as_bool(r.get("pass_fimo_1e5", False)):
                fl = cluster_flags[cid]
        updated.append(classify_row_enzymatic(r, fl))

    if "enzymatic_validation_class" not in orig_fields:
        idx = orig_fields.index("functional_evidence") if "functional_evidence" in orig_fields else len(orig_fields)
        out_fields = orig_fields[:idx] + ["enzymatic_validation_class"] + orig_fields[idx:]
    else:
        out_fields = orig_fields

    with open(tsv_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(updated)
    return updated


def update_phase1_summary(v_rows: list[dict], b_rows: list[dict]) -> None:
    """Update SBL counts inside Final_Phase1_Dataset_Summary.txt."""
    sum_path = C1_DIR / "Final_Phase1_Dataset_Summary.txt"
    if not sum_path.is_file():
        return

    def count_sbl(accs: tuple[str, ...]) -> tuple[int, int, int]:
        v_sbl = sum(
            1 for r in v_rows
            if r["winning_acc"] in accs
            and as_bool(r["pass_cov65"]) and as_bool(r["pass_fimo_1e5"])
            and r["functional_tier"] == "Putative_Hydrolytic_BetaLactamase"
        )
        v85_sbl = sum(
            1 for r in v_rows
            if r["winning_acc"] in accs
            and as_bool(r["pass_cov65"]) and as_bool(r["pass_fimo_1e5"])
            and as_bool(r["in_viral_nr85_rep"])
            and r["functional_tier"] == "Putative_Hydrolytic_BetaLactamase"
        )
        b85_sbl = sum(
            1 for r in b_rows
            if r["winning_acc"] in accs
            and as_bool(r["pass_cov65"]) and as_bool(r["pass_fimo_1e5"])
            and not as_bool(r["is_viral_contig_in_genome"])
            and as_bool(r["in_nr85_rep"])
            and r["functional_tier"] == "Putative_Hydrolytic_BetaLactamase"
        )
        return v_sbl, v85_sbl, b85_sbl

    tag_map = {
        "PF00144": ("PF00144",),
        "PF00905": ("PF00905",),
        "PF13354": ("PF13354",),
        "PF00768": ("PF00768",),
        "PF13354_PF00768": ("PF13354", "PF00768"),
    }
    lines = sum_path.read_text(encoding="utf-8").splitlines()
    new_lines: list[str] = []
    for line in lines:
        m = re.match(r"^\s*(PF00144|PF00905|PF13354_PF00768|PF13354|PF00768)\s*:", line)
        if m:
            tag = m.group(1)
            vs, v85s, b85s = count_sbl(tag_map[tag])
            line = re.sub(r"(\(SBL=)[0-9,]+(\) -> Viral nr85=[0-9,]+ \(SBL=)[0-9,]+(, MOSAiC=)", rf"\g<1>{vs:,}\g<2>{v85s:,}\g<3>", line)
            line = re.sub(r"(Bac nr85=[0-9,]+ \(SBL=)[0-9,]+(\))", rf"\g<1>{b85s:,}\g<2>", line)
        new_lines.append(line)
    sum_path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")


def export_viral_verification_table(v_rows: list[dict]) -> None:
    """Write dedicated 4,203-row strict viral enzymatic verification table."""
    cols = [
        "uid", "source", "winning_acc", "viral_origin_tier", "seq_len", "hmm_cov",
        "in_viral_nr85_rep", "viral_cluster_size",
        "motif1_SXXK_hmm", "motif2_class", "motif3_KTG_hmm", "has_omega_E166",
        "best_card_hit", "card_gene_family", "card_pident", "card_cov",
        "functional_tier", "enzymatic_validation_class", "functional_evidence", "domain_architecture",
    ]
    strict_rows = [
        {c: r.get(c, "") for c in cols}
        for r in v_rows
        if as_bool(r.get("pass_cov65")) and as_bool(r.get("pass_fimo_1e5"))
    ]
    out_p = C1_DIR / "Viral_Enzymatic_Function_Verification_2026.tsv"
    with open(out_p, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(strict_rows)
    print(f"  -> Exported {len(strict_rows):,} rows to {out_p.name}")


def recompute_stage3_hgt_tables(v_rows: list[dict], b_rows: list[dict]) -> None:
    """Recompute functional columns across All_Trees, Primary, OAT_Matrix, and HGT Summary."""
    meta_by_uid: dict[str, dict] = {}
    for r in v_rows + b_rows:
        meta_by_uid[r["uid"]] = {
            "functional_tier": r["functional_tier"],
            "card_pident": float(r.get("card_pident") or 0.0),
            "card_cov": float(r.get("card_cov") or 0.0),
            "card_gene_family": r.get("card_gene_family", "NONE"),
        }

    def refresh_event_file(path: Path) -> pd.DataFrame:
        df = pd.read_csv(path, sep="\t", low_memory=False)
        n_hyd_list, n_pbp_list, n_c40_list, cfams_list, fcat_list = [], [], [], [], []
        for _, row in df.iterrows():
            tips = [t.strip() for t in str(row.get("all_recipient_tips", "")).split(";") if t.strip() and t.strip() != "nan"]
            n_hyd = 0
            n_pbp = 0
            n_c40 = 0
            cfams: set[str] = set()
            for uid in tips:
                info = meta_by_uid.get(uid, {})
                ft = info.get("functional_tier", "PBP_or_PRE_Homolog")
                if ft == "Putative_Hydrolytic_BetaLactamase":
                    n_hyd += 1
                    if info.get("card_pident", 0.0) >= 40.0 and info.get("card_cov", 0.0) >= 70.0:
                        n_c40 += 1
                        cf = str(info.get("card_gene_family", "NONE"))
                        if cf and cf != "NONE":
                            cfams.add(cf.split(";")[0].strip())
                else:
                    n_pbp += 1
            fcat = "Putative_Hydrolytic_BetaLactamase" if (n_hyd >= n_pbp and n_hyd > 0) else "PBP_or_PRE_Homolog"
            n_hyd_list.append(n_hyd)
            n_pbp_list.append(n_pbp)
            n_c40_list.append(n_c40)
            cfams_list.append(";".join(sorted(cfams)[:5]))
            fcat_list.append(fcat)

        df["n_hydrolytic_sbl_tips"] = n_hyd_list
        df["n_pbp_pre_tips"] = n_pbp_list
        df["n_card_40_70_tips"] = n_c40_list
        df["card_families"] = cfams_list
        df["functional_category"] = fcat_list
        df.to_csv(path, sep="\t", index=False)
        return df

    df_all = refresh_event_file(HGT_DIR / "All_Trees_HGT_Events_Long.tsv")
    df_prim = refresh_event_file(HGT_DIR / "Primary_HGT_Events_Catalog_2026.tsv")

    # Update tier1_hydrolytic_sbl_events in OAT_Sensitivity_Summary_Matrix_2026.tsv
    mat_path = HGT_DIR / "OAT_Sensitivity_Summary_Matrix_2026.tsv"
    df_mat = pd.read_csv(mat_path, sep="\t", low_memory=False)
    t1_sbl_counts = (
        df_all[df_all["is_tier1_hgt"].astype(str).str.lower().isin(["true", "1"]) & (df_all["functional_category"] == "Putative_Hydrolytic_BetaLactamase")]
        .groupby(["family", "tree_label", "rooting_method"])
        .size()
        .to_dict()
    )
    df_mat["tier1_hydrolytic_sbl_events"] = [
        int(t1_sbl_counts.get((r["family"], r["tree_label"], r["rooting_method"]), 0))
        for _, r in df_mat.iterrows()
    ]
    df_mat.to_csv(mat_path, sep="\t", index=False)

    # Regenerate HGT_Summary_By_Family_And_Function_2026.tsv
    families = ("PF00144", "PF00905", "PF13354", "PF00768", "PF13354_PF00768")
    sum_rows = []
    for fam in families:
        sub_f = df_prim[df_prim["family"] == fam]
        for d in ("Bacterial_to_Viral", "Viral_to_Bacterial"):
            for fc in ("Putative_Hydrolytic_BetaLactamase", "PBP_or_PRE_Homolog"):
                s = sub_f[(sub_f["direction"] == d) & (sub_f["functional_category"] == fc)]
                t1 = s[s["is_tier1_hgt"].astype(str).str.lower().isin(["true", "1"])]
                t2 = s[s["is_tier2_hgt"].astype(str).str.lower().isin(["true", "1"])]
                pb = s[s["is_prophage_borne"].astype(str).str.lower().isin(["true", "1"])]
                rc = s[s["robust_core_hgt"].astype(str).str.lower().isin(["true", "1"])]
                ri = t1[t1["root_Invariant"].astype(str).str.lower().isin(["true", "1"])]
                sum_rows.append({
                    "family": fam,
                    "direction": d,
                    "functional_category": fc,
                    "all_nested_concordant": len(s),
                    "tier1_events": len(t1),
                    "tier1_root_invariant": len(ri),
                    "tier1_robust_core": len(t1[t1["robust_core_hgt"].astype(str).str.lower().isin(["true", "1"])]),
                    "tier1_recip_nr85_tips": int(t1["n_focal_recipient_tips"].sum()),
                    "tier1_expanded_viral_seqs": int(t1["expanded_viral_cluster_seqs"].sum()),
                    "tier2_events": len(t2),
                    "robust_core_total_t1_t2": len(rc),
                    "prophage_borne_events": len(pb),
                })
    df_sum = pd.DataFrame(sum_rows)
    df_sum.to_csv(HGT_DIR / "HGT_Summary_By_Family_And_Function_2026.tsv", sep="\t", index=False)

    # Regenerate Stage3_HGT_Summary_Report.txt
    cl_path = HGT_DIR / "CL0013_Superfamily_Rooting_Validation.tsv"
    df_cl = pd.read_csv(cl_path, sep="\t") if cl_path.is_file() else pd.DataFrame()
    lines = [
        "=== Stage 3 (Phase 4): Pre-Locked Primary HGT Inference & OAT Robustness Summary ===",
        f"Trees analyzed across all families/perturbations: {len(df_mat)}",
        "",
        "--- Primary Nested Concordant HGT Summary (by Family x Direction x Function) ---",
        df_sum[df_sum["all_nested_concordant"] > 0].to_string(index=False) if len(df_sum) else "(none)",
        "",
        "--- Primary Trees Overview ---",
        df_mat[df_mat["oat_key"] == "primary_out"][[
            "family", "total_tips", "viral_tips", "bacterial_tips",
            "nested_fitch_ml_events", "tier1_hgt_events", "tier1_B_to_V_events",
            "tier1_V_to_B_events", "tier1_hydrolytic_sbl_events",
            "tier1_B_to_V_viral_nr85_tips", "tier1_B_to_V_expanded_seqs",
            "tier2_hgt_events", "prophage_borne_events",
        ]].to_string(index=False),
    ]
    if len(df_cl) > 0:
        lines.extend(["", "--- CL0013 Superfamily & LOO Rooting Validation ---", df_cl.to_string(index=False)])
    (HGT_DIR / "Stage3_HGT_Summary_Report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  -> Refreshed Stage 3 HGT tables across {len(df_mat)} trees.")


def main() -> int:
    print("[1/4] Loading local domain FASTAs & aligning against official Pfam HMMs (A2M)...")
    seqs_by_fam = load_domain_sequences()
    flags_by_fam: dict[str, dict[str, dict]] = {}
    for fam in ("PF00144", "PF00905", "PF13354"):
        flags_by_fam[fam] = compute_a2m_motif_flags(fam, seqs_by_fam[fam])
        print(f"  {fam}: A2M-aligned {len(flags_by_fam[fam]):,} sequences")

    print("[2/4] Updating Final_Viral_SBL_Catalog.tsv & Final_Bacterial_TripleSource_Catalog.tsv...")
    v_rows = update_catalog(C1_DIR / "Final_Viral_SBL_Catalog.tsv", flags_by_fam, is_bacterial=False)
    b_rows = update_catalog(C1_DIR / "Final_Bacterial_TripleSource_Catalog.tsv", flags_by_fam, is_bacterial=True)
    update_phase1_summary(v_rows, b_rows)
    export_viral_verification_table(v_rows)

    strict_v = [r for r in v_rows if as_bool(r["pass_cov65"]) and as_bool(r["pass_fimo_1e5"])]
    print("  Strict curated viral breakdown (n=%d):" % len(strict_v))
    for (fam, ft, vc), cnt in sorted(Counter((r["winning_acc"], r["functional_tier"], r["enzymatic_validation_class"]) for r in strict_v).items()):
        print(f"    {fam:8s} | {ft:34s} | {vc:34s} -> {cnt}")

    if "--sync-hgt" in sys.argv[1:]:
        print("[3/4] Recomputing HGT functional annotations (--sync-hgt)...")
        recompute_stage3_hgt_tables(v_rows, b_rows)
    else:
        print("[3/4] HGT tables left unchanged (HPC catalogs already carry these tiers; see --sync-hgt)")
    print("[4/4] Done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
