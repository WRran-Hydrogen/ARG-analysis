#!/usr/bin/env python3
"""Selection analysis, step 1: verified codon alignments, dS saturation & Xia's test,
domain-anchored purifying selection tests, and HyPhy dataset construction.

Builds:
  1. 1:1 verified codon alignments (hmmalign match-state anchored + MAFFT/pal2nal)
     for PF00144, PF00905, PF13354, PF00768.
  2. Formal dS saturation diagnostics (pairwise Nei-Gojobori pS / dS, Xia et al. 2003
     substitution saturation index I_SS vs simulated symmetric & asymmetric critical
     thresholds I_SS.c, and Bruen's Pairwise Homoplasy Index PHI).
  3. Per-site HMM match-column codon dN/dS (omega_NG), amino-acid Shannon entropy H_aa,
     nucleotide 1st+2nd vs 3rd position entropy, and catalytic motif purifying selection
     tests (Mann-Whitney, Fisher, binomial, and 10,000-iteration permutation tests).
  4. HyPhy datasets (codon.fna + pruned rooted tree.nwk + id_map + col_map):
       - Group A: 4-family matched N=90 benchmark (45 viral + 45 bacterial, hmm track)
       - Group B: PF00144 2x2 factorial benchmark ([hmm vs mafft] x [N=90 vs N=240])
       - Group C: Up to 4 low-divergence HGT & monophyletic subclades per family (8 <= N <= 55)
"""

from __future__ import annotations

import csv
import math
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from ete4 import Tree
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/hpc/hpc_config.py
import hpc_config as hpc  # noqa: E402

csv.field_size_limit(10_000_000)

CURATED_DIR = hpc.BASE / "01_final_curated"
INPUTS_DIR = hpc.BASE / "inputs"
P3_ALN_DIR = hpc.PHASE3_DIR / "alignments"
P3_TREES_DIR = hpc.PHASE3_DIR / "trees"
HGT_DIR = hpc.BASE / "04_hgt_inference"
OUT_DIR = hpc.BASE / "05_selection_pressure"
CODON_DIR = OUT_DIR / "codon_alignments"
HYPHY_DIR = OUT_DIR / "hyphy_sets"

PAL2NAL_BIN = hpc.PAL2NAL_BIN
ESL_REFORMAT_BIN = hpc.setting("SBL_HPC_ESL_REFORMAT", str(Path(hpc.HMMALIGN_BIN).with_name("esl-reformat")))
FAMILIES = ("PF00144", "PF00905", "PF13354", "PF00768")
THREADS = 16

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
SENSE_CODONS = sorted(c for c, aa in CODON_TABLE.items() if aa != "*")
CODON_TO_IDX = {c: i for i, c in enumerate(SENSE_CODONS)}
BASES = ("A", "C", "G", "T")


def _build_ng_tables():
    """Precompute Nei-Gojobori (1986) site counts and pairwise N/S differences for 61 sense codons."""
    n_c = len(SENSE_CODONS)
    site_n = np.zeros(n_c, dtype=np.float64)
    site_s = np.zeros(n_c, dtype=np.float64)
    for i, c in enumerate(SENSE_CODONS):
        sc = 0.0
        nc = 0.0
        for pos in range(3):
            for b in BASES:
                if b == c[pos]:
                    continue
                c2 = c[:pos] + b + c[pos + 1:]
                aa2 = CODON_TABLE[c2]
                if aa2 == "*":
                    continue
                if aa2 == CODON_TABLE[c]:
                    sc += 1.0 / 3.0
                else:
                    nc += 1.0 / 3.0
        site_n[i] = nc
        site_s[i] = sc

    diff_n = np.zeros((n_c, n_c), dtype=np.float64)
    diff_s = np.zeros((n_c, n_c), dtype=np.float64)
    for i, c1 in enumerate(SENSE_CODONS):
        for j in range(i + 1, n_c):
            c2 = SENSE_CODONS[j]
            diffs = [p for p in range(3) if c1[p] != c2[p]]
            nd = len(diffs)
            if nd == 1:
                if CODON_TABLE[c1] == CODON_TABLE[c2]:
                    sd, ndiff = 1.0, 0.0
                else:
                    sd, ndiff = 0.0, 1.0
            elif nd == 2:
                paths_s = 0.0
                paths_n = 0.0
                valid_p = 0
                for first, second in ((diffs[0], diffs[1]), (diffs[1], diffs[0])):
                    mid = c1[:first] + c2[first] + c1[first + 1:]
                    if CODON_TABLE[mid] == "*":
                        continue
                    valid_p += 1
                    paths_s += 1.0 if CODON_TABLE[c1] == CODON_TABLE[mid] else 0.0
                    paths_n += 0.0 if CODON_TABLE[c1] == CODON_TABLE[mid] else 1.0
                    paths_s += 1.0 if CODON_TABLE[mid] == CODON_TABLE[c2] else 0.0
                    paths_n += 0.0 if CODON_TABLE[mid] == CODON_TABLE[c2] else 1.0
                if valid_p > 0:
                    sd, ndiff = paths_s / valid_p, paths_n / valid_p
                else:
                    sd, ndiff = 0.0, 2.0
            else:
                # 3-nt difference: average over non-stop 6 pathways
                paths_s = 0.0
                paths_n = 0.0
                valid_p = 0
                for p0 in range(3):
                    for p1 in range(3):
                        if p1 == p0:
                            continue
                        p2 = 3 - p0 - p1
                        m1 = c1[:p0] + c2[p0] + c1[p0 + 1:]
                        m2 = m1[:p1] + c2[p1] + m1[p1 + 1:]
                        if CODON_TABLE[m1] == "*" or CODON_TABLE[m2] == "*":
                            continue
                        valid_p += 1
                        for a, b in ((c1, m1), (m1, m2), (m2, c2)):
                            if CODON_TABLE[a] == CODON_TABLE[b]:
                                paths_s += 1.0
                            else:
                                paths_n += 1.0
                if valid_p > 0:
                    sd, ndiff = paths_s / valid_p, paths_n / valid_p
                else:
                    sd, ndiff = 0.25, 2.75
            diff_n[i, j] = diff_n[j, i] = ndiff
            diff_s[i, j] = diff_s[j, i] = sd
    return site_n, site_s, diff_n, diff_s


SITE_N, SITE_S, DIFF_N, DIFF_S = _build_ng_tables()


# ---------------------------------------------------------------------------
# FASTA / Subprocess Utilities
# ---------------------------------------------------------------------------
def run(cmd: str, stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    else:
        with open(stdout_path, "w", encoding="utf-8") as fh:
            p = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {cmd}\n{p.stderr}")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, str]] = []
    name: str | None = None
    buf: list[str] = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path: Path, recs: list[tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def hmm_leng(hmm_path: Path) -> int:
    with open(hmm_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("LENG "):
                return int(line.split()[1])
    raise ValueError(f"No LENG in {hmm_path}")


def clean_verified_pairs(
    aa_recs: list[tuple[str, str]], cds_map: dict[str, str]
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """Keep sequences with exact 1:1 CDS translation match and zero internal stop codons."""
    ok_aa: list[tuple[str, str]] = []
    ok_cds: list[tuple[str, str]] = []
    for uid, aa in aa_recs:
        nt = cds_map.get(uid, "").upper().replace("U", "T")
        if not nt:
            continue
        # Handle optional trailing stop codon if CDS is 3*(len(aa)+1)
        if len(nt) == 3 * (len(aa) + 1) and nt[-3:] in ("TAA", "TAG", "TGA"):
            nt = nt[:-3]
        if len(nt) != 3 * len(aa):
            continue
        cods = [nt[i:i + 3] for i in range(0, len(nt), 3)]
        # Mask trailing stop if present, reject internal stops
        if cods and cods[-1] in ("TAA", "TAG", "TGA"):
            continue
        has_internal_stop = any(c in ("TAA", "TAG", "TGA") for c in cods)
        if has_internal_stop:
            continue
        # Verify residue agreement for unambiguous standard codons
        mismatch = False
        clean_cods = []
        clean_aa = []
        for a, c in zip(aa.upper(), cods):
            if any(b not in BASES for b in c):
                clean_cods.append("---")
                clean_aa.append(a if a != "*" else "X")
                continue
            trans = CODON_TABLE.get(c, "X")
            if a not in ("X", "B", "Z", "J", "U", "O") and trans != a:
                mismatch = True
                break
            clean_cods.append(c)
            clean_aa.append(trans if a in ("U", "O", "*") else a)
        if mismatch:
            continue
        ok_aa.append((uid, "".join(clean_aa)))
        ok_cds.append((uid, "".join(clean_cods)))
    return ok_aa, ok_cds


def backtranslate_msa_direct(msa_aa: list[tuple[str, str]], cds_map: dict[str, str]) -> list[tuple[str, str]]:
    """Exact 1:1 back-translation of an AA MSA using verified codon strings (handles masked --- codons)."""
    out: list[tuple[str, str]] = []
    for uid, arow in msa_aa:
        nt = cds_map[uid]
        cods = [nt[i:i + 3] for i in range(0, len(nt), 3)]
        idx = 0
        buf = []
        for ch in arow:
            if ch in ("-", "."):
                buf.append("---")
            else:
                buf.append(cods[idx] if idx < len(cods) else "---")
                idx += 1
        out.append((uid, "".join(buf)))
    return out


def build_codon_alignments_for_family(
    fam: str, aa_recs: list[tuple[str, str]], cds_recs: list[tuple[str, str]]
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    fdir = CODON_DIR / fam
    fdir.mkdir(parents=True, exist_ok=True)
    cds_map = dict(cds_recs)

    unaligned_aa = fdir / f"{fam}.verified.faa"
    unaligned_cds = fdir / f"{fam}.verified.fna"
    write_fasta(unaligned_aa, aa_recs)
    write_fasta(unaligned_cds, cds_recs)

    # Track 1: hmmalign (A2M -> extract exact LENG HMM match-state codons)
    hmm = INPUTS_DIR / f"{fam}_official.hmm"
    leng = hmm_leng(hmm)
    a2m_path = fdir / f"{fam}.hmm.a2m"
    afa_path = fdir / f"{fam}.hmm.afa"
    run(f"{hpc.HMMALIGN_BIN} --outformat A2M {hmm} {unaligned_aa}", stdout_path=a2m_path)
    run(f"{ESL_REFORMAT_BIN} afa {a2m_path}", stdout_path=afa_path)

    # Identify match columns in the rectangular AFA (columns without lowercase or '.')
    afa_recs = read_fasta(afa_path)
    aln_len = len(afa_recs[0][1])
    is_match_col = [True] * aln_len
    for _, seq in afa_recs:
        for j, ch in enumerate(seq):
            if ch == "." or ch.islower():
                is_match_col[j] = False
    mcols = [j for j, ok in enumerate(is_match_col) if ok]
    if len(mcols) != leng:
        raise RuntimeError(f"{fam}: expected {leng} match cols in AFA, got {len(mcols)}")

    # Back-translate full AFA (also run pal2nal cross-check where possible, fallback to direct)
    upper_afa = [(u, s.replace(".", "-").upper()) for u, s in afa_recs]
    upper_afa_path = fdir / f"{fam}.hmm.upper.afa"
    write_fasta(upper_afa_path, upper_afa)

    hmm_full_codon_raw = fdir / f"{fam}.hmm.pal2nal_raw.fna"
    pal_ok = False
    try:
        run(f"{PAL2NAL_BIN} {upper_afa_path} {unaligned_cds} -output fasta", stdout_path=hmm_full_codon_raw)
        pal_recs = read_fasta(hmm_full_codon_raw)
        if len(pal_recs) == len(upper_afa) and len(pal_recs[0][1]) == 3 * aln_len:
            pal_ok = True
    except Exception:
        pal_ok = False

    if pal_ok:
        full_codon_recs = read_fasta(hmm_full_codon_raw)
    else:
        full_codon_recs = backtranslate_msa_direct(upper_afa, cds_map)

    # Slice exact LENG match-state codons (3 * LENG nt)
    hmm_match_codons = []
    for uid, nt_row in full_codon_recs:
        sliced = "".join(nt_row[3 * j:3 * j + 3] for j in mcols)
        hmm_match_codons.append((uid, sliced))
    hmm_out_path = fdir / f"{fam}.full.hmm_codon.fna"
    write_fasta(hmm_out_path, hmm_match_codons)

    # Track 2: MAFFT --auto -> pal2nal.pl -> codon-triplet <50% gap trimming
    mafft_aa_path = fdir / f"{fam}.mafft.faa"
    mafft_raw_codon = fdir / f"{fam}.mafft.pal2nal_raw.fna"
    run(f"{hpc.MAFFT_BIN} --auto --anysymbol --thread {THREADS} {unaligned_aa}", stdout_path=mafft_aa_path)
    mafft_aa = read_fasta(mafft_aa_path)
    m_len = len(mafft_aa[0][1])

    pal_m_ok = False
    try:
        run(f"{PAL2NAL_BIN} {mafft_aa_path} {unaligned_cds} -output fasta", stdout_path=mafft_raw_codon)
        pal_m = read_fasta(mafft_raw_codon)
        if len(pal_m) == len(mafft_aa) and len(pal_m[0][1]) == 3 * m_len:
            pal_m_ok = True
    except Exception:
        pal_m_ok = False

    mafft_codons = read_fasta(mafft_raw_codon) if pal_m_ok else backtranslate_msa_direct(mafft_aa, cds_map)
    n_seq = len(mafft_codons)
    keep_mcols = [
        j for j in range(m_len)
        if (sum(1 for _, s in mafft_codons if "-" in s[3 * j:3 * j + 3]) / max(n_seq, 1)) < 0.50
    ]
    mafft_trimmed = [
        (uid, "".join(s[3 * j:3 * j + 3] for j in keep_mcols))
        for uid, s in mafft_codons
    ]
    mafft_out_path = fdir / f"{fam}.full.mafft_codon.fna"
    write_fasta(mafft_out_path, mafft_trimmed)

    print(
        f"[{fam}] Codon alignments built: n={n_seq} (hmm={leng} codons, "
        f"mafft={len(keep_mcols)} codons; pal2nal={pal_ok}/{pal_m_ok})",
        flush=True,
    )
    return hmm_match_codons, mafft_trimmed


# ---------------------------------------------------------------------------
# Per-Site Codon Metrics, Purifying Selection Tests, dS Saturation, Xia's Test & PHI
# ---------------------------------------------------------------------------
def shannon_entropy(chars: list[str], ignore: tuple[str, ...] = ("-", "X", "*")) -> float:
    valid = [c for c in chars if c not in ignore]
    if not valid:
        return 0.0
    n = len(valid)
    cnt = Counter(valid)
    return -sum((v / n) * math.log2(v / n) for v in cnt.values())


def nei_gojobori_column_omega(col_codons: list[str], min_seqs: int = 10) -> float:
    valid = [c for c in col_codons if c in CODON_TO_IDX]
    if len(valid) < min_seqs:
        return float("nan")
    cnt = Counter(valid)
    items = [(CODON_TO_IDX[c], w) for c, w in cnt.items()]
    tot_nd = tot_sd = tot_en = tot_es = 0.0
    for i in range(len(items)):
        i1, w1 = items[i]
        n1, s1 = SITE_N[i1], SITE_S[i1]
        tot_en += w1 * (w1 - 1) * n1
        tot_es += w1 * (w1 - 1) * s1
        for k in range(i + 1, len(items)):
            i2, w2 = items[k]
            w = 2.0 * w1 * w2
            tot_en += w * 0.5 * (n1 + SITE_N[i2])
            tot_es += w * 0.5 * (s1 + SITE_S[i2])
            tot_nd += w * DIFF_N[i1, i2]
            tot_sd += w * DIFF_S[i1, i2]
    dn = tot_nd / max(tot_en, 1e-6)
    ds = tot_sd / max(tot_es, 1e-6)
    if ds < 1e-4:
        return float("nan")
    return float(min(dn / ds, 10.0))


def compute_persite_metrics(fam: str, codon_recs: list[tuple[str, str]]) -> tuple[pd.DataFrame, dict]:
    n_seq = len(codon_recs)
    n_cod = len(codon_recs[0][1]) // 3
    cod_mat = [[s[3 * j:3 * j + 3].upper() for j in range(n_cod)] for _, s in codon_recs]
    aa_seqs = ["".join(CODON_TABLE.get(c, "-") for c in row) for row in cod_mat]

    # Identify exact modal HMM columns for Motif1 (SXXK), Motif2, Motif3 (KTG), and OmegaLoop (PF13354)
    motif_defs = [
        ("Motif1_SXXK", r"S[A-Z]{2}K", 4),
        ("Motif2_SDN_YAN_SXN", r"[SY][A-Z]N", 3),
        ("Motif3_KTG", r"[KH][TS]G", 3),
    ]
    if fam == "PF13354":
        motif_defs.append(("OmegaLoop_E166", r"E[A-Z]{2}[LIVM]N", 5))

    motif_map: dict[int, str] = {}
    ext_set: set[int] = set()
    for mname, pat, mlen in motif_defs:
        cnt = Counter()
        for aas in aa_seqs:
            for m in re.finditer(pat, aas):
                cnt[m.start()] += 1
        if cnt:
            top_col = cnt.most_common(1)[0][0]
            for j in range(top_col, min(n_cod, top_col + mlen)):
                motif_map[j] = mname
            for j in range(max(0, top_col - 1), min(n_cod, top_col + mlen + 1)):
                ext_set.add(j)

    rows = []
    for j in range(n_cod):
        col_c = [row[j] for row in cod_mat]
        col_a = [CODON_TABLE.get(c, "-") for c in col_c]
        occ = sum(1 for c in col_c if c in CODON_TO_IDX) / max(n_seq, 1)
        h_aa = shannon_entropy(col_a)
        h_p1 = shannon_entropy([c[0] for c in col_c if c in CODON_TO_IDX], ignore=())
        h_p2 = shannon_entropy([c[1] for c in col_c if c in CODON_TO_IDX], ignore=())
        h_p3 = shannon_entropy([c[2] for c in col_c if c in CODON_TO_IDX], ignore=())
        om_ng = nei_gojobori_column_omega(col_c)
        mrole = motif_map.get(j, "Background")
        rows.append({
            "family": fam,
            "hmm_col": j + 1,
            "motif_role": mrole,
            "in_catalytic_motif_core": mrole != "Background",
            "in_catalytic_pocket_ext": j in ext_set,
            "occupancy": round(occ, 4),
            "H_aa": round(h_aa, 4),
            "H_nt_pos12": round(0.5 * (h_p1 + h_p2), 4),
            "H_nt_pos3": round(h_p3, 4),
            "codon_omega_ng": round(om_ng, 5) if not math.isnan(om_ng) else float("nan"),
        })
    df = pd.DataFrame(rows)

    # Purifying selection test comparing motif core vs background
    mot = df[df["in_catalytic_motif_core"]]
    bkg = df[~df["in_catalytic_motif_core"]]
    mot_om = mot["codon_omega_ng"].dropna().values
    bkg_om = bkg["codon_omega_ng"].dropna().values
    mot_h = mot["H_aa"].dropna().values
    bkg_h = bkg["H_aa"].dropna().values

    _, p_mw_om = stats.mannwhitneyu(mot_om, bkg_om, alternative="less") if (len(mot_om) and len(bkg_om)) else (0, np.nan)
    _, p_mw_h = stats.mannwhitneyu(mot_h, bkg_h, alternative="less") if (len(mot_h) and len(bkg_h)) else (0, np.nan)

    # Strongly constrained sites: omega_NG < 0.25 (or H_aa < 1.5 bits when synonymous ds~0)
    mot_con = int(((mot["codon_omega_ng"] < 0.25) | ((mot["H_aa"] < 1.5) & mot["codon_omega_ng"].isna())).sum())
    bkg_con = int(((bkg["codon_omega_ng"] < 0.25) | ((bkg["H_aa"] < 1.5) & bkg["codon_omega_ng"].isna())).sum())
    or_con, p_fish_con = stats.fisher_exact(
        [[mot_con, len(mot) - mot_con], [bkg_con, len(bkg) - bkg_con]],
        alternative="greater",
    ) if (len(mot) and len(bkg)) else (np.nan, np.nan)
    p_binom = float(stats.binomtest(mot_con, len(mot), p=max(bkg_con / max(len(bkg), 1), 1e-6), alternative="greater").pvalue) if len(mot) else np.nan

    # 10,000-iteration permutation test on mean H_aa and mean omega_NG
    rng = np.random.default_rng(20261121)
    all_h = df["H_aa"].dropna().values
    k_m = len(mot_h)
    obs_dh = float(np.mean(bkg_h) - np.mean(mot_h)) if k_m else np.nan
    cnt_perm_h = 0
    if k_m:
        for _ in range(10000):
            perm = rng.permutation(all_h)
            if (np.mean(perm[k_m:]) - np.mean(perm[:k_m])) >= obs_dh:
                cnt_perm_h += 1
        p_perm_h = (cnt_perm_h + 1) / 10001.0
    else:
        p_perm_h = np.nan

    proof = {
        "family": fam,
        "n_hmm_cols": len(df),
        "n_motif_core_cols": len(mot),
        "n_bg_cols": len(bkg),
        "motif_median_omega_ng": round(float(np.median(mot_om)), 4) if len(mot_om) else np.nan,
        "bg_median_omega_ng": round(float(np.median(bkg_om)), 4) if len(bkg_om) else np.nan,
        "MW_p_omega_motif_lt_bg": float(p_mw_om),
        "motif_median_H_aa": round(float(np.median(mot_h)), 4) if len(mot_h) else np.nan,
        "bg_median_H_aa": round(float(np.median(bkg_h)), 4) if len(bkg_h) else np.nan,
        "MW_p_Haa_motif_lt_bg": float(p_mw_h),
        "Perm10k_p_Haa_motif_lt_bg": round(float(p_perm_h), 6) if pd.notna(p_perm_h) else np.nan,
        "motif_constrained_pct": round(100.0 * mot_con / max(len(mot), 1), 2),
        "bg_constrained_pct": round(100.0 * bkg_con / max(len(bkg), 1), 2),
        "Fisher_OR_constrained": round(float(or_con), 3) if pd.notna(or_con) else np.nan,
        "Fisher_p_constrained": float(p_fish_con),
        "Binomial_p_constrained": p_binom,
    }
    return df, proof


def _sim_hky_site_entropies(
    k: int, n_sites: int, pi: np.ndarray, kappa: float, tree_type: str, rng: np.random.Generator
) -> np.ndarray:
    """Simulate critical saturation entropy I_SS.c on symmetric vs asymmetric k-OTU trees (Xia et al. 2003)."""
    # Build HKY85 transition matrix for branch length t_crit representing saturation regime
    # Total tree height ~2.0 subst/site across k taxa; internal vs terminal scaling
    q = np.outer(np.ones(4), pi).astype(np.float64)
    # Transitions: A<->G (0<->2), C<->T (1<->3)
    for i, j in ((0, 2), (2, 0), (1, 3), (3, 1)):
        q[i, j] *= kappa
    np.fill_diagonal(q, 0.0)
    np.fill_diagonal(q, -q.sum(axis=1))
    rate_scale = -float(np.sum(pi * np.diag(q)))
    q /= max(rate_scale, 1e-12)

    def expm_hky(t: float) -> np.ndarray:
        vals, vecs = np.linalg.eig(q * t)
        P = (vecs * np.exp(vals)) @ np.linalg.inv(vecs)
        P = np.real(P)
        P = np.clip(P, 0.0, 1.0)
        P /= P.sum(axis=1, keepdims=True)
        return P

    root_states = rng.choice(4, size=n_sites, p=pi)
    if tree_type == "sym":
        # Balanced binary tree of depth >= log2(k), total root-to-tip ~1.0 (saturation limit)
        depth = max(2, int(math.ceil(math.log2(max(k, 2)))))
        b_len = 1.0 / depth
        P_step = expm_hky(b_len)
        layer = [root_states]
        for _ in range(depth):
            nxt = []
            for st in layer:
                for _c in range(2):
                    u = rng.random(n_sites)
                    cdf = np.cumsum(P_step[st], axis=1)
                    child = (u[:, None] > cdf).sum(axis=1)
                    nxt.append(np.clip(child, 0, 3))
            layer = nxt
        tips = np.stack(layer[:k], axis=0)
    else:
        # Asymmetric pectinate tree of k OTUs (Xia et al. 2003 worst-case topology)
        b_int = 1.0 / max(k - 1, 1)
        b_tip = 0.75
        P_int = expm_hky(b_int)
        P_tip = expm_hky(b_tip)
        tips_list = []
        cur = root_states
        for step in range(k - 1):
            u = rng.random(n_sites)
            cdf = np.cumsum(P_tip[cur], axis=1)
            leaf_st = np.clip((u[:, None] > cdf).sum(axis=1), 0, 3)
            tips_list.append(leaf_st)
            u2 = rng.random(n_sites)
            cdf2 = np.cumsum(P_int[cur], axis=1)
            cur = np.clip((u2[:, None] > cdf2).sum(axis=1), 0, 3)
        tips_list.append(cur)
        tips = np.stack(tips_list, axis=0)

    # Compute mean site Shannon entropy (bits) across simulated k taxa
    ents = np.zeros(n_sites, dtype=np.float64)
    for b in range(4):
        p = (tips == b).mean(axis=0)
        mask = p > 0
        ents[mask] -= p[mask] * np.log2(p[mask])
    return ents


def compute_saturation_and_phi(
    label: str, family: str, set_group: str, codon_recs: list[tuple[str, str]]
) -> dict:
    """Compute pairwise Nei-Gojobori dS/dN, Xia's I_SS vs I_SS.c, and Bruen's PHI."""
    rng = np.random.default_rng(20261121)
    n_seq = len(codon_recs)
    n_cod = len(codon_recs[0][1]) // 3

    # Subsample up to 200 sequences evenly for pairwise dS distribution
    if n_seq > 200:
        idxs = np.linspace(0, n_seq - 1, 200, dtype=int)
        sub_recs = [codon_recs[i] for i in idxs]
    else:
        sub_recs = codon_recs

    m = len(sub_recs)
    cmat = np.full((m, n_cod), -1, dtype=np.int32)
    for i, (_, s) in enumerate(sub_recs):
        su = s.upper()
        for j in range(n_cod):
            c = su[3 * j:3 * j + 3]
            if c in CODON_TO_IDX:
                cmat[i, j] = CODON_TO_IDX[c]

    # Pairwise NG pS, dS, pN, dN (chunked vectorized lookup)
    ps_list = []
    ds_list = []
    pn_list = []
    dn_list = []
    sat_cnt = 0
    tot_pairs = 0

    valid_mask = cmat >= 0
    for i in range(m):
        ci = cmat[i]
        vi = valid_mask[i]
        for j in range(i + 1, m):
            both = vi & valid_mask[j]
            if int(both.sum()) < 30:
                continue
            a = ci[both]
            b = cmat[j][both]
            tot_pairs += 1
            en = 0.5 * (SITE_N[a].sum() + SITE_N[b].sum())
            es = 0.5 * (SITE_S[a].sum() + SITE_S[b].sum())
            nd = DIFF_N[a, b].sum()
            sd = DIFF_S[a, b].sum()
            ps = sd / max(es, 1e-6)
            pn = nd / max(en, 1e-6)
            ps_list.append(ps)
            pn_list.append(pn)
            if ps >= 0.74:
                ds_val = 3.0
                sat_cnt += 1
            else:
                ds_val = -0.75 * math.log(1.0 - (4.0 / 3.0) * ps)
                if ds_val >= 2.5:
                    sat_cnt += 1
            ds_list.append(ds_val)
            if pn >= 0.74:
                dn_list.append(3.0)
            else:
                dn_list.append(-0.75 * math.log(1.0 - (4.0 / 3.0) * pn))

    # Xia's test on 3rd codon positions & 1st+2nd codon positions (sites with >=50% ACGT)
    nt_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    full_nt = np.full((n_seq, 3 * n_cod), -1, dtype=np.int8)
    for i, (_, s) in enumerate(codon_recs):
        for j, ch in enumerate(s.upper()):
            if ch in nt_map:
                full_nt[i, j] = nt_map[ch]

    def xia_for_positions(pos_cols: np.ndarray) -> tuple[float, float, float, float, float, int]:
        sub = full_nt[:, pos_cols]
        occ = (sub >= 0).mean(axis=0)
        sub = sub[:, occ >= 0.50]
        ns = sub.shape[1]
        if ns < 10 or n_seq < 4:
            return (np.nan, np.nan, np.nan, np.nan, np.nan, ns)

        # Global base frequencies pi and rough transition/transversion ratio kappa
        flat = sub[sub >= 0]
        counts = np.bincount(flat, minlength=4).astype(np.float64)
        pi = np.clip(counts / max(counts.sum(), 1.0), 0.01, 0.97)
        pi /= pi.sum()

        # Estimate empirical kappa from 100 random sequence pairs
        ti = tv = 0
        for _ in range(min(200, n_seq * (n_seq - 1) // 2)):
            a, b = rng.choice(n_seq, size=2, replace=False)
            ok = (sub[a] >= 0) & (sub[b] >= 0)
            va, vb = sub[a, ok], sub[b, ok]
            diff = va != vb
            xa, xb = va[diff], vb[diff]
            is_ti = ((xa ^ xb) == 2)  # A(0)<->G(2), C(1)<->T(3)
            ti += int(is_ti.sum())
            tv += int((~is_ti).sum())
        kappa = max(0.5, min(10.0, (2.0 * ti) / max(tv, 1)))

        k_otu = min(32, n_seq)
        # Finite-sample full-saturation expected entropy H_FSS(k_otu, pi)
        fss_draws = rng.multinomial(k_otu, pi, size=1000) / k_otu
        h_fss_vals = np.zeros(1000, dtype=np.float64)
        for b in range(4):
            p = fss_draws[:, b]
            mk = p > 0
            h_fss_vals[mk] -= p[mk] * np.log2(p[mk])
        h_fss = float(np.mean(h_fss_vals))

        # 100 replicates of k_otu taxa for I_SS
        iss_reps = []
        for _ in range(100):
            r_idx = rng.choice(n_seq, size=k_otu, replace=False)
            block = sub[r_idx]
            h_sites = np.zeros(ns, dtype=np.float64)
            valid_c = (block >= 0).sum(axis=0)
            for b in range(4):
                p = (block == b).sum(axis=0) / np.maximum(valid_c, 1)
                mk = (p > 0) & (valid_c >= 3)
                h_sites[mk] -= p[mk] * np.log2(p[mk])
            iss_reps.append(float(np.mean(h_sites[valid_c >= 3])) / max(h_fss, 1e-6))

        iss_arr = np.array(iss_reps)
        # Critical thresholds from HKY85 simulation under symmetric & asymmetric trees
        h_sym = _sim_hky_site_entropies(k_otu, min(ns, 400), pi, kappa, "sym", rng)
        h_asym = _sim_hky_site_entropies(k_otu, min(ns, 400), pi, kappa, "asym", rng)
        iss_c_sym = float(np.mean(h_sym)) / max(h_fss, 1e-6)
        iss_c_asym = float(np.mean(h_asym)) / max(h_fss, 1e-6)

        _, p_sym = stats.ttest_1samp(iss_arr, iss_c_sym)
        _, p_asym = stats.ttest_1samp(iss_arr, iss_c_asym)
        return (
            round(float(np.mean(iss_arr)), 4),
            round(iss_c_sym, 4),
            round(iss_c_asym, 4),
            float(p_sym),
            float(p_asym),
            ns,
        )

    pos3_cols = np.arange(2, 3 * n_cod, 3)
    pos12_cols = np.array([j for j in range(3 * n_cod) if j % 3 != 2])
    iss3, iss3_cs, iss3_ca, p3_sym, p3_asym, n_pos3 = xia_for_positions(pos3_cols)
    iss12, iss12_cs, iss12_ca, _, _, n_pos12 = xia_for_positions(pos12_cols)

    # Bruen PHI test (same algorithm as scripts/07_selection_meme_and_catalytic_motifs.py)
    nt_seqs_sub = [s for _, s in sub_recs]
    phi_obs, phi_null, p_phi, phi_inf = compute_phi_homoplasy(nt_seqs_sub, n_perm=300, seed=42)

    return {
        "dataset": label,
        "family": family,
        "set_group": set_group,
        "n_seqs": n_seq,
        "n_codons": n_cod,
        "n_pairs_evaluated": tot_pairs,
        "median_pS": round(float(np.median(ps_list)), 4) if ps_list else np.nan,
        "median_dS": round(float(np.median(ds_list)), 4) if ds_list else np.nan,
        "iqr_dS": f"{np.percentile(ds_list, 25):.3f}-{np.percentile(ds_list, 75):.3f}" if ds_list else "NA",
        "pct_saturated_pairs_pS75_or_dS25": round(100.0 * sat_cnt / max(tot_pairs, 1), 2),
        "median_dN": round(float(np.median(dn_list)), 4) if dn_list else np.nan,
        "pos3_n_sites": n_pos3,
        "pos3_Iss": iss3,
        "pos3_Iss_c_sym": iss3_cs,
        "pos3_Iss_c_asym": iss3_ca,
        "pos3_t_p_sym": p3_sym,
        "pos3_t_p_asym": p3_asym,
        "pos3_unsaturated_sym": bool(pd.notna(iss3) and iss3 < iss3_cs and p3_sym < 0.05),
        "pos12_Iss": iss12,
        "pos12_Iss_c_sym": iss12_cs,
        "pos12_unsaturated_sym": bool(pd.notna(iss12) and iss12 < iss12_cs),
        "PHI_obs": round(phi_obs, 4) if pd.notna(phi_obs) else np.nan,
        "PHI_p": round(p_phi, 4) if pd.notna(p_phi) else np.nan,
        "PHI_n_informative_sites": phi_inf,
    }


def compute_phi_homoplasy(nt_seqs: list[str], n_perm: int = 300, seed: int = 42):
    rng = np.random.default_rng(seed)
    idx = np.linspace(0, len(nt_seqs) - 1, min(80, len(nt_seqs)), dtype=int)
    seqs = [nt_seqs[i].upper() for i in idx]
    L = len(seqs[0])
    inf_cols = []
    for col_i in range(L):
        col = [s[col_i] for s in seqs]
        counts = {}
        for ch in col:
            if ch in BASES:
                counts[ch] = counts.get(ch, 0) + 1
        if sum(1 for v in counts.values() if v >= 2) >= 2:
            inf_cols.append(col)
    m = len(inf_cols)
    if m < 10:
        return np.nan, np.nan, np.nan, m
    if m > 200:
        sub_idx = np.linspace(0, m - 1, 200, dtype=int)
        inf_cols = [inf_cols[k] for k in sub_idx]
        m = len(inf_cols)

    top2_masks = []
    for col in inf_cols:
        counts = {}
        for i, ch in enumerate(col):
            if ch in BASES:
                counts.setdefault(ch, set()).add(i)
        ordered = sorted(counts.values(), key=lambda s: -len(s))
        top2_masks.append((ordered[0], ordered[1]))

    full_mat = np.zeros((m, m), dtype=float)
    for j1 in range(m):
        a0, a1 = top2_masks[j1]
        for j2 in range(j1 + 1, m):
            b0, b1 = top2_masks[j2]
            gam = int(bool(a0 & b0)) + int(bool(a0 & b1)) + int(bool(a1 & b0)) + int(bool(a1 & b1))
            if gam == 4:
                full_mat[j1, j2] = 1.0
                full_mat[j2, j1] = 1.0

    w = max(5, min(50, m // 3))

    def phi_of_order(order):
        sub = full_mat[np.ix_(order, order)]
        s = 0.0
        c = 0
        for d in range(1, w + 1):
            diag = np.diagonal(sub, offset=d)
            s += diag.sum()
            c += len(diag)
        return s / max(c, 1)

    obs_phi = phi_of_order(np.arange(m))
    base_ord = np.arange(m)
    null_phis = np.array([phi_of_order(rng.permutation(base_ord)) for _ in range(n_perm)])
    p_phi = float((np.sum(null_phis <= obs_phi) + 1) / (n_perm + 1))
    return float(obs_phi), float(np.mean(null_phis)), p_phi, m


# ---------------------------------------------------------------------------
# Tree Rooting & HyPhy Set Extraction
# ---------------------------------------------------------------------------
def load_rooted_ingroup_tree(fam: str) -> Tree:
    tpath = P3_TREES_DIR / fam / f"{fam}.primary.hmm_out.treefile"
    og_path = P3_ALN_DIR / fam / "outgroup_taxa.txt"
    og_ids = {l.strip() for l in og_path.read_text(encoding="utf-8").splitlines() if l.strip()}
    t = Tree(tpath.read_text(encoding="utf-8"), parser=1)
    all_leaves = {l.name for l in t.leaves()}
    og_present = og_ids & all_leaves
    in_leaves = all_leaves - og_present

    n_og = len(og_present)
    n_tot = len(all_leaves)
    for nd in t.traverse("postorder"):
        if nd.is_leaf:
            nd.add_prop("_og", 1 if nd.name in og_present else 0)
            nd.add_prop("_sz", 1)
        else:
            nd.add_prop("_og", sum(c.props["_og"] for c in nd.children))
            nd.add_prop("_sz", sum(c.props["_sz"] for c in nd.children))

    best_nd = None
    best_f1 = -1.0
    for nd in t.traverse():
        if nd.up is None:
            continue
        k_og, k_sz = nd.props["_og"], nd.props["_sz"]
        if k_og > 0:
            f1 = 2.0 * k_og / (k_sz + n_og)
            if f1 > best_f1:
                best_f1, best_nd = f1, nd
        if n_og - k_og > 0:
            f1 = 2.0 * (n_og - k_og) / ((n_tot - k_sz) + n_og)
            if f1 > best_f1:
                best_f1, best_nd = f1, nd
    if best_nd is not None:
        t.set_outgroup(best_nd)
    t.prune(sorted(in_leaves), preserve_branch_length=True)
    return t


def pick_balanced_taxa(
    tree: Tree,
    valid_ids: set[str],
    v_meta: dict[str, dict],
    hgt_priority_ids: set[str],
    n_vir: int,
    n_bac: int,
) -> list[str]:
    """Select balanced viral + bacterial taxa along tree order, prioritizing SBL/Tier1/HGT viral tips."""
    ordered = [l.name for l in tree.leaves() if l.name in valid_ids]
    v_all = [u for u in ordered if u in v_meta]
    b_all = [u for u in ordered if u not in v_meta]

    # Viral selection: prioritize hydrolytic SBL + Tier 1 + HGT recipients, then evenly spaced
    v_prio = [
        u for u in v_all
        if v_meta[u]["functional_tier"] == "Putative_Hydrolytic_BetaLactamase"
        or v_meta[u]["viral_origin_tier"] == "Tier1_HighConfidence_Viral"
        or u in hgt_priority_ids
    ]
    if len(v_prio) >= n_vir:
        idxs = np.linspace(0, len(v_prio) - 1, n_vir, dtype=int)
        v_sel = [v_prio[i] for i in idxs]
    else:
        v_set = set(v_prio)
        rem = [u for u in v_all if u not in v_set]
        need = min(n_vir - len(v_prio), len(rem))
        idxs = np.linspace(0, len(rem) - 1, need, dtype=int) if need > 0 else []
        v_sel = v_prio + [rem[i] for i in idxs]

    # Bacterial selection: prioritize HGT sisters/recipients + evenly spaced along phylogeny
    b_prio = [u for u in b_all if u in hgt_priority_ids]
    if len(b_prio) >= min(10, n_bac):
        b_prio = b_prio[:min(10, n_bac)]
    b_set = set(b_prio)
    b_rem = [u for u in b_all if u not in b_set]
    need_b = min(n_bac - len(b_prio), len(b_rem))
    idxs_b = np.linspace(0, len(b_rem) - 1, need_b, dtype=int) if need_b > 0 else []
    b_sel = b_prio + [b_rem[i] for i in idxs_b]
    return v_sel + b_sel


def select_low_divergence_subclades(
    fam: str,
    tree: Tree,
    valid_ids: set[str],
    v_meta: dict[str, dict],
    b_meta: dict[str, dict],
    hgt_df: pd.DataFrame,
    codon_map: dict[str, str],
) -> list[tuple[str, list[str], str]]:
    """Select up to 4 low-divergence HGT and monophyletic focal subclades (8 <= N <= 55)."""
    # Annotate clade domain & median pairwise nucleotide identity proxy (branch length / depth)
    for nd in tree.traverse("postorder"):
        if nd.is_leaf:
            ok = 1 if nd.name in valid_ids else 0
            is_v = 1 if (ok and nd.name in v_meta) else 0
            nd.add_prop("_ok", ok)
            nd.add_prop("_v", is_v)
            nd.add_prop("_b", ok - is_v)
            nd.add_prop("_max_h", 0.0)
        else:
            nd.add_prop("_ok", sum(c.props["_ok"] for c in nd.children))
            nd.add_prop("_v", sum(c.props["_v"] for c in nd.children))
            nd.add_prop("_b", sum(c.props["_b"] for c in nd.children))
            nd.add_prop("_max_h", max(c.props["_max_h"] + float(c.dist or 0.0) for c in nd.children))

    # Candidate 1: HGT recipient + donor sister subclades from Stage 3 catalog
    hgt_sub = hgt_df[
        (hgt_df["family"] == fam)
        & (hgt_df["is_tier1_hgt"] | hgt_df["is_tier2_hgt"] | hgt_df["robust_core_hgt"])
    ].copy()
    hgt_anchors = set()
    for s in hgt_sub["all_recipient_tips"].dropna():
        for u in str(s).split(";"):
            if u in valid_ids:
                hgt_anchors.add(u)

    def median_nt_id_and_ps(uids: list[str]) -> tuple[float, float]:
        sample = uids[:min(25, len(uids))]
        ids = []
        ps_vals = []
        for i in range(len(sample)):
            s1 = codon_map[sample[i]]
            c1 = [s1[k:k + 3].upper() for k in range(0, len(s1), 3)]
            for j in range(i + 1, len(sample)):
                s2 = codon_map[sample[j]]
                same = tot = 0
                for a, b in zip(s1, s2):
                    if a != "-" and b != "-":
                        tot += 1
                        if a == b:
                            same += 1
                if tot >= 60:
                    ids.append(same / tot)
                # Fast Nei-Gojobori pS estimate on shared sense codons
                en_es = sd_sum = es_sum = 0.0
                c2 = [s2[k:k + 3].upper() for k in range(0, len(s2), 3)]
                for a_c, b_c in zip(c1, c2):
                    ia, ib = CODON_TO_IDX.get(a_c, -1), CODON_TO_IDX.get(b_c, -1)
                    if ia >= 0 and ib >= 0:
                        es_sum += 0.5 * (SITE_S[ia] + SITE_S[ib])
                        sd_sum += DIFF_S[ia, ib]
                        en_es += 1.0
                if en_es >= 30 and es_sum > 0:
                    ps_vals.append(sd_sum / es_sum)
        return (
            float(np.median(ids)) if ids else 0.0,
            float(np.median(ps_vals)) if ps_vals else 1.0,
        )

    cands = []
    for nd in tree.traverse("postorder"):
        if nd.is_leaf or nd.up is None:
            continue
        ok = nd.props["_ok"]
        if not (8 <= ok <= 55):
            continue
        nv, nb = nd.props["_v"], nd.props["_b"]
        leaves = [l.name for l in nd.leaves() if l.name in valid_ids]
        nt_id, med_ps = median_nt_id_and_ps(leaves)
        unsat_tier = 0 if (nt_id >= 0.70 and med_ps < 0.72) else (1 if med_ps < 0.72 else 2)
        n_hgt = sum(1 for u in leaves if u in hgt_anchors)
        if nv >= 3 and nb >= 3 and n_hgt >= 1:
            kind = "HGT_clade"
            prio = (0, unsat_tier, med_ps, -nt_id, -n_hgt, -ok)
        elif nv >= 0.85 * ok and nv >= 6:
            kind = "Viral_clade"
            prio = (1, unsat_tier, med_ps, -nt_id, -nv, -ok)
        elif nb >= 0.85 * ok:
            # Prefer bacterial clades containing hydrolytic SBLs or CARD refs
            n_sbl = sum(
                1 for u in leaves
                if (v_meta.get(u) or b_meta.get(u) or {}).get("functional_tier") == "Putative_Hydrolytic_BetaLactamase"
            )
            kind = "Bacterial_clade"
            prio = (2, unsat_tier, med_ps, -nt_id, -n_sbl, -ok)
        else:
            continue
        cands.append((prio, kind, nt_id, med_ps, set(leaves), leaves))

    cands.sort(key=lambda x: x[0])
    chosen: list[tuple[str, list[str], str]] = []
    chosen_sets: list[set[str]] = []
    kind_counts: Counter = Counter()

    for _, kind, nt_id, med_ps, lset, leaves in cands:
        # Deduplicate overlapping clades (Jaccard < 0.50) and cap per kind to ensure diversity
        if any((len(lset & cs) / max(len(lset | cs), 1)) >= 0.50 for cs in chosen_sets):
            continue
        if kind_counts[kind] >= 2:
            continue
        idx = len(chosen) + 1
        sid = f"{fam}_sub{idx}_{kind}"
        chosen.append((sid, leaves, f"nt_id={nt_id:.3f},pS={med_ps:.3f}"))
        chosen_sets.append(lset)
        kind_counts[kind] += 1
        if len(chosen) >= 4:
            break
    return chosen


def write_hyphy_dataset(
    set_id: str,
    family: str,
    set_group: str,
    track: str,
    selected_ids: list[str],
    codon_map: dict[str, str],
    full_tree: Tree,
    v_meta: dict[str, dict],
) -> dict:
    sdir = HYPHY_DIR / set_id
    sdir.mkdir(parents=True, exist_ok=True)

    # Prune tree to selected_ids (unique, preserving order)
    sel_unique = list(dict.fromkeys(u for u in selected_ids if u in codon_map))
    t = full_tree.copy()
    t.prune(sel_unique, preserve_branch_length=True)

    # Ensure strictly matching tip set and strip internal support labels so HyPhy assigns clean names
    final_tips = [l.name for l in t.leaves()]
    n_cod = len(codon_map[final_tips[0]]) // 3
    # Drop all-gap or <20%-occupied codon columns while tracking col_map
    seqs_raw = [codon_map[u] for u in final_tips]
    keep_cols = []
    for j in range(n_cod):
        if any("-" not in s[3 * j:3 * j + 3] for s in seqs_raw):
            keep_cols.append(j)

    id_map = {u: f"SEQ{i + 1:04d}" for i, u in enumerate(final_tips)}
    for nd in t.traverse():
        if nd.is_leaf:
            nd.name = id_map[nd.name]
        else:
            nd.name = ""

    trimmed_recs = []
    for u in final_tips:
        s = codon_map[u]
        trimmed_recs.append((id_map[u], "".join(s[3 * j:3 * j + 3] for j in keep_cols)))

    write_fasta(sdir / "codon.fna", trimmed_recs)
    # Write unquoted Newick with branch lengths only (parser=1/format=5 style)
    nw = t.write(parser=5)
    (sdir / "tree.nwk").write_text(nw + ("\n" if not nw.endswith("\n") else ""), encoding="utf-8")

    with open(sdir / "id_map.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["hyphy_id", "orig_uid", "domain"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        for u in final_tips:
            w.writerow({"hyphy_id": id_map[u], "orig_uid": u, "domain": "Viral" if u in v_meta else "Bacterial"})

    with open(sdir / "col_map.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["hyphy_site_1based", "orig_col_1based"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        for idx, j in enumerate(keep_cols, start=1):
            w.writerow({"hyphy_site_1based": idx, "orig_col_1based": j + 1})

    sat_row = compute_saturation_and_phi(set_id, family, set_group, trimmed_recs)
    return {
        "task_row": {
            "set_id": set_id,
            "family": family,
            "set_group": set_group,
            "track": track,
            "n_seqs": len(final_tips),
            "n_codons": len(keep_cols),
            "median_dS": sat_row["median_dS"],
            "iss": sat_row["pos3_Iss"],
            "iss_c_sym": sat_row["pos3_Iss_c_sym"],
            "iss_c_asym": sat_row["pos3_Iss_c_asym"],
            "set_dir": str(sdir),
        },
        "sat_row": sat_row,
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CODON_DIR.mkdir(parents=True, exist_ok=True)
    HYPHY_DIR.mkdir(parents=True, exist_ok=True)

    # Load Viral & Bacterial catalog metadata for priority sampling
    v_meta: dict[str, dict] = {}
    with open(CURATED_DIR / "Final_Viral_SBL_Catalog.tsv", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            v_meta[row["uid"]] = row
    b_meta: dict[str, dict] = {}
    with open(CURATED_DIR / "Final_Bacterial_TripleSource_Catalog.tsv", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            b_meta[row["uid"]] = row

    hgt_path = HGT_DIR / "Primary_HGT_Events_Catalog_2026.tsv"
    hgt_df = pd.read_csv(hgt_path, sep="\t") if hgt_path.is_file() else pd.DataFrame(
        columns=["family", "is_tier1_hgt", "is_tier2_hgt", "robust_core_hgt", "all_recipient_tips"]
    )

    persite_dfs = []
    proof_rows = []
    sat_rows = []
    hyphy_tasks = []

    for fam in FAMILIES:
        v_aa = read_fasta(CURATED_DIR / f"{fam}_viral_nr85_dom.faa")
        b_aa = read_fasta(CURATED_DIR / f"{fam}_bacterial_tree_dom.faa")
        v_cds = dict(read_fasta(CURATED_DIR / f"{fam}_viral_nr85_dom_cds.fna"))
        b_cds = dict(read_fasta(CURATED_DIR / f"{fam}_bacterial_tree_dom_cds.fna"))

        ok_aa_v, ok_cds_v = clean_verified_pairs(v_aa, v_cds)
        ok_aa_b, ok_cds_b = clean_verified_pairs(b_aa, b_cds)
        all_aa = ok_aa_v + ok_aa_b
        all_cds = ok_cds_v + ok_cds_b

        hmm_codons, mafft_codons = build_codon_alignments_for_family(fam, all_aa, all_cds)
        hmm_map = dict(hmm_codons)
        mafft_map = dict(mafft_codons)
        valid_ids = set(hmm_map.keys()) & set(mafft_map.keys())

        # Per-site domain metrics & catalytic motif purifying selection test
        df_site, proof = compute_persite_metrics(fam, hmm_codons)
        persite_dfs.append(df_site)
        proof_rows.append(proof)

        # Full-family dS saturation & Xia's test
        sat_rows.append(compute_saturation_and_phi(f"{fam}_full_hmm", fam, "full_family_hmm", hmm_codons))
        if fam == "PF00144":
            sat_rows.append(compute_saturation_and_phi(f"{fam}_full_mafft", fam, "full_family_mafft", mafft_codons))

        # Load primary Sister-PRE rooted ingroup tree
        tree = load_rooted_ingroup_tree(fam)
        tree_ids = {l.name for l in tree.leaves()}
        valid_tree_ids = valid_ids & tree_ids

        # HGT priority taxa for this family
        hgt_prio = set()
        if len(hgt_df) > 0:
            sub_h = hgt_df[hgt_df["family"] == fam]
            for s in sub_h[sub_h["is_tier1_hgt"] | sub_h["robust_core_hgt"]]["all_recipient_tips"].dropna():
                hgt_prio.update(str(s).split(";"))

        # Group A: N=90 (45 viral + 45 bacterial) matched benchmark on hmm track
        sel_90 = pick_balanced_taxa(tree, valid_tree_ids, v_meta, hgt_prio, n_vir=45, n_bac=45)
        res_90 = write_hyphy_dataset(f"{fam}_N90_hmm", fam, "benchmark_N90", "hmm", sel_90, hmm_map, tree, v_meta)
        hyphy_tasks.append(res_90["task_row"])
        sat_rows.append(res_90["sat_row"])

        # Group B: PF00144 2x2 factorial ([hmm vs mafft] x [N=90 vs N=240])
        if fam == "PF00144":
            res_90_m = write_hyphy_dataset(f"{fam}_N90_mafft", fam, "benchmark_2x2", "mafft", sel_90, mafft_map, tree, v_meta)
            hyphy_tasks.append(res_90_m["task_row"])
            sat_rows.append(res_90_m["sat_row"])

            sel_240 = pick_balanced_taxa(tree, valid_tree_ids, v_meta, hgt_prio, n_vir=120, n_bac=120)
            res_240_h = write_hyphy_dataset(f"{fam}_N240_hmm", fam, "benchmark_2x2", "hmm", sel_240, hmm_map, tree, v_meta)
            res_240_m = write_hyphy_dataset(f"{fam}_N240_mafft", fam, "benchmark_2x2", "mafft", sel_240, mafft_map, tree, v_meta)
            hyphy_tasks.extend([res_240_h["task_row"], res_240_m["task_row"]])
            sat_rows.extend([res_240_h["sat_row"], res_240_m["sat_row"]])

        # Group C: Low-divergence HGT & monophyletic focal subclades (8 <= N <= 55)
        subclades = select_low_divergence_subclades(fam, tree, valid_tree_ids, v_meta, b_meta, hgt_df, hmm_map)
        for sid, s_leaves, note in subclades:
            res_sub = write_hyphy_dataset(sid, fam, "subclade", "hmm", s_leaves, hmm_map, tree, v_meta)
            hyphy_tasks.append(res_sub["task_row"])
            sat_rows.append(res_sub["sat_row"])
            print(f"  [Subclade] {sid}: n={len(s_leaves)} ({note}), median_dS={res_sub['sat_row']['median_dS']}", flush=True)

    # Pooled row in the purifying selection table
    df_persite = pd.concat(persite_dfs, ignore_index=True)
    mot_all = df_persite[df_persite["in_catalytic_motif_core"]]
    bkg_all = df_persite[~df_persite["in_catalytic_motif_core"]]
    mot_om_a = mot_all["codon_omega_ng"].dropna().values
    bkg_om_a = bkg_all["codon_omega_ng"].dropna().values
    mot_h_a = mot_all["H_aa"].dropna().values
    bkg_h_a = bkg_all["H_aa"].dropna().values
    _, p_mw_om_a = stats.mannwhitneyu(mot_om_a, bkg_om_a, alternative="less")
    _, p_mw_h_a = stats.mannwhitneyu(mot_h_a, bkg_h_a, alternative="less")
    mot_con_a = int(((mot_all["codon_omega_ng"] < 0.25) | ((mot_all["H_aa"] < 1.5) & mot_all["codon_omega_ng"].isna())).sum())
    bkg_con_a = int(((bkg_all["codon_omega_ng"] < 0.25) | ((bkg_all["H_aa"] < 1.5) & bkg_all["codon_omega_ng"].isna())).sum())
    or_a, p_fish_a = stats.fisher_exact([[mot_con_a, len(mot_all) - mot_con_a], [bkg_con_a, len(bkg_all) - bkg_con_a]], alternative="greater")
    p_bin_a = float(stats.binomtest(mot_con_a, len(mot_all), p=bkg_con_a / max(len(bkg_all), 1), alternative="greater").pvalue)
    proof_rows.append({
        "family": "All_Pooled",
        "n_hmm_cols": len(df_persite),
        "n_motif_core_cols": len(mot_all),
        "n_bg_cols": len(bkg_all),
        "motif_median_omega_ng": round(float(np.median(mot_om_a)), 4),
        "bg_median_omega_ng": round(float(np.median(bkg_om_a)), 4),
        "MW_p_omega_motif_lt_bg": float(p_mw_om_a),
        "motif_median_H_aa": round(float(np.median(mot_h_a)), 4),
        "bg_median_H_aa": round(float(np.median(bkg_h_a)), 4),
        "MW_p_Haa_motif_lt_bg": float(p_mw_h_a),
        "Perm10k_p_Haa_motif_lt_bg": np.nan,
        "motif_constrained_pct": round(100.0 * mot_con_a / max(len(mot_all), 1), 2),
        "bg_constrained_pct": round(100.0 * bkg_con_a / max(len(bkg_all), 1), 2),
        "Fisher_OR_constrained": round(float(or_a), 3),
        "Fisher_p_constrained": float(p_fish_a),
        "Binomial_p_constrained": p_bin_a,
    })

    df_persite.to_csv(OUT_DIR / "PerSite_Domain_Selection_Metrics_2026.tsv", sep="\t", index=False)
    df_proof = pd.DataFrame(proof_rows)
    df_proof.to_csv(OUT_DIR / "Catalytic_Motif_Purifying_Proof_2026.tsv", sep="\t", index=False)
    df_sat = pd.DataFrame(sat_rows)
    df_sat.to_csv(OUT_DIR / "dS_Saturation_And_Xia_Test_2026.tsv", sep="\t", index=False)
    df_tasks = pd.DataFrame(hyphy_tasks)
    df_tasks.to_csv(OUT_DIR / "hyphy_tasks.tsv", sep="\t", index=False)

    print("\n=== Catalytic Motif Purifying Selection Tests ===")
    print(df_proof[["family", "motif_median_omega_ng", "bg_median_omega_ng", "MW_p_omega_motif_lt_bg", "motif_median_H_aa", "bg_median_H_aa", "MW_p_Haa_motif_lt_bg", "Fisher_OR_constrained"]].to_string(index=False))
    print("\n=== dS Saturation & Xia's Test ===")
    print(df_sat[["dataset", "n_seqs", "median_dS", "pct_saturated_pairs_pS75_or_dS25", "pos3_Iss", "pos3_Iss_c_sym", "pos3_Iss_c_asym", "pos3_unsaturated_sym", "PHI_p"]].to_string(index=False))
    print(f"\nWrote hyphy_tasks.tsv ({len(df_tasks)} datasets)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
