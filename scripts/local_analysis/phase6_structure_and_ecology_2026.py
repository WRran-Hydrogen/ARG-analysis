#!/usr/bin/env python3
"""Stage 4B (2026 Curated SBL Re-analysis): 3D Biophysics, Confounder-Controlled
Structural/Phylogenetic Diversity, and Stratified Ecology (CMH + BH-FDR).

Outputs in results/reanalysis_2026/06_structure_and_ecology/:
  1. PerSite_3D_Biophysics_2026.tsv
  2. Table_3D_Structural_Localization_2026.tsv
  3. PDB_3D_Structural_Catalog_2026.tsv
  4. Table_Structural_Diversity_Controlled_2026.tsv
  5. Table_Ecology_CMH_Controlled_2026.tsv
  6. Stage4_Structure_And_Ecology_Summary_Report.txt
"""

from __future__ import annotations

import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Auto-dispatch under WSL Python when Bio / scipy / ete4 are absent on Windows
if sys.platform.startswith("win"):
    try:
        import Bio.PDB  # noqa: F401
        import ete4  # noqa: F401
        import scipy.stats  # noqa: F401
    except ImportError:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import sbl_config as _cfg

        cmd = _cfg.linux_cmd(_cfg.WSL_PY, _cfg.to_wsl_path(Path(__file__).resolve()), *sys.argv[1:])
        proc = subprocess.run(cmd)
        sys.exit(proc.returncode)

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
from ete4 import Tree
from scipy.stats import fisher_exact, levene, mannwhitneyu, norm

# ---------------------------------------------------------------------------
# Constants & Paths
# ---------------------------------------------------------------------------
BASE = Path(__file__).resolve().parents[2]
DATA_DIR = BASE / "data"
RES_DIR = BASE / "results" / "reanalysis_2026"
SEL_DIR = RES_DIR / "05_selection_pressure"
OUT_DIR = RES_DIR / "06_structure_and_ecology"

PDB_DIR = DATA_DIR / "structures" / "ProteinStructure"
REF_PDB_DIR = DATA_DIR / "structures" / "ref_pdbs"
MAXCC_NWK = DATA_DIR / "phylogeny" / "maxcc.nwk"
MOTIF_NODE = DATA_DIR / "motifs" / "motif_node.csv"
MOTIF_EDGE = DATA_DIR / "motifs" / "motif_edge.csv"
INFO_XLSX = DATA_DIR / "metadata" / "information.xlsx"
VIRAL_CAT = RES_DIR / "01_curation_and_context" / "Final_Viral_SBL_Catalog.tsv"
BAC_CAT = RES_DIR / "01_curation_and_context" / "Final_Bacterial_TripleSource_Catalog.tsv"
METRICS_TSV = SEL_DIR / "PerSite_Domain_Selection_Metrics_2026.tsv"

FAMILIES = ["PF00144", "PF00905", "PF13354", "PF00768"]

REF_PDBS = {
    "PF00144": ("1bls_P99.pdb", "1bls_P99"),
    "PF00905": ("1k25_pneumoniae_PF00905.pdb", "1k25_pneumoniae"),
    "PF13354": ("1bue_ENTEROBACTER_PF13354.pdb", "1bue_Enterobacter"),
}

CRYSTAL_BAC_INFO = {
    "1axb_ESCHERICHIA_PF13354": ("PF13354", "Putative_Hydrolytic_BetaLactamase"),
    "1blp_STAPHYLOCOCCUS_PF13354": ("PF13354", "Putative_Hydrolytic_BetaLactamase"),
    "1bls_P99": ("PF00144", "Putative_Hydrolytic_BetaLactamase"),
    "1bue_ENTEROBACTER_PF13354": ("PF13354", "Putative_Hydrolytic_BetaLactamase"),
    "1cef_STREPTOMYCES_PF00144": ("PF00144", "PBP_or_PRE_Homolog"),
    "1ei5_OCHROBACTRUM_PF00144": ("PF00144", "PBP_or_PRE_Homolog"),
    "1k25_pneumoniae_PF00905": ("PF00905", "PBP_or_PRE_Homolog"),
    "1mwr_methicillin_resistant_Staphylococcus_PF00905": ("PF00905", "PBP_or_PRE_Homolog"),
}

MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K",
}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

RELIABLE_HOST_METHODS = {"Isolate taxonomy", "CRISPR spacer match"}
HIGH_CONF_TIERS = {
    "Tier1_HighConfidence_Viral",
    "Tier2_Supported_Viral",
    "MOSAiC_Proteome_Verified_vMAG",
}


# ---------------------------------------------------------------------------
# Generic Statistical & Taxonomy Helpers
# ---------------------------------------------------------------------------
def bh_fdr(pvals: np.ndarray | list[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        val = min(prev, ranked[i] * n / (i + 1))
        adj[i] = val
        prev = val
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0.0, 1.0)
    return out


def cmh_test_stratified(in_feat: np.ndarray, in_hab: np.ndarray, strata: np.ndarray) -> tuple[float, float]:
    """Cochran-Mantel-Haenszel pooled odds ratio and two-sided normal-approx p."""
    num_or = 0.0
    den_or = 0.0
    sum_diff = 0.0
    sum_var = 0.0
    for s in np.unique(strata):
        m = strata == s
        f = in_feat[m]
        h = in_hab[m]
        n = len(f)
        if n < 2:
            continue
        a = float(np.sum(f & h))
        b = float(np.sum(f & ~h))
        c = float(np.sum(~f & h))
        d = float(np.sum(~f & ~h))
        r1, r0 = a + b, c + d
        c1, c0 = a + c, b + d
        if r1 == 0 or r0 == 0 or c1 == 0 or c0 == 0:
            continue
        num_or += (a * d) / n
        den_or += (b * c) / n
        e_a = (r1 * c1) / n
        var_a = (r1 * r0 * c1 * c0) / (n * n * max(n - 1, 1))
        sum_diff += a - e_a
        sum_var += var_a
    if sum_var <= 0 or den_or <= 0:
        return np.nan, np.nan
    or_cmh = num_or / den_or
    z = abs(sum_diff) / math.sqrt(sum_var)
    p = float(2.0 * (1.0 - norm.cdf(z)))
    return float(or_cmh), p


def safe_mw(x: np.ndarray, y: np.ndarray, alt: str = "two-sided") -> float:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return np.nan
    if np.all(x == x[0]) and np.all(y == y[0]) and x[0] == y[0]:
        return 1.0
    return float(mannwhitneyu(x, y, alternative=alt).pvalue)


def extract_taxon_rank(tax_str: object, prefix: str) -> str:
    if not isinstance(tax_str, str):
        return "Unassigned"
    m = re.search(rf"{prefix}__([^;]+)", tax_str)
    if m:
        val = m.group(1).strip()
        if val and val.lower() not in {"unclassified", "unknown", ""}:
            return val
    return "Unassigned"


def classify_ecosystem(eco_str: object) -> str:
    """Standardize IMG/VR Ecosystem classification into interpretable niches."""
    if not isinstance(eco_str, str):
        return "Unclassified"
    s = eco_str.strip()
    if not s or s.replace(";", "").strip() == "":
        return "Unclassified"
    sl = s.lower()
    if "human" in sl and "digestive" in sl:
        return "Human_Digestive"
    if "human" in sl and "respiratory" in sl:
        return "Human_Respiratory"
    if "human" in sl:
        return "Human_Other"
    if "digestive" in sl:
        return "NonHuman_Digestive"
    if "freshwater" in sl:
        return "Freshwater"
    if "marine" in sl:
        return "Marine"
    if "soil" in sl or "rhizosphere" in sl or "peat" in sl:
        return "Soil_Terrestrial"
    if "engineered" in sl or "wastewater" in sl or "bioreactor" in sl or "sludge" in sl:
        return "Engineered_Wastewater"
    if "host-associated" in sl:
        return "Other_Host_Associated"
    if "environmental" in sl:
        return "Other_Environmental"
    return "Unclassified"


# ---------------------------------------------------------------------------
# Alignment Translation & Pairwise HMM Mapping (Part 1)
# ---------------------------------------------------------------------------
def translate_codon_aln(fna_path: Path) -> tuple[dict[str, str], list[tuple[str, str]], str]:
    """Translate {fam}.full.hmm_codon.fna -> exact 1:1 hmm_col amino acid strings."""
    seqs_all: list[tuple[str, str]] = []
    by_id: dict[str, str] = {}
    cur_h = None
    chunks: list[str] = []
    with open(fna_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_h is not None:
                    aa = _codon_to_aa("".join(chunks))
                    seqs_all.append((cur_h, aa))
                    by_id[cur_h] = aa
                    by_id.setdefault(cur_h.split("|")[0], aa)
                cur_h = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if cur_h is not None:
            aa = _codon_to_aa("".join(chunks))
            seqs_all.append((cur_h, aa))
            by_id[cur_h] = aa
            by_id.setdefault(cur_h.split("|")[0], aa)

    if not seqs_all:
        return {}, [], ""
    ncols = len(seqs_all[0][1])
    cons_chars: list[str] = []
    for j in range(ncols):
        col = [aa[j] for _, aa in seqs_all if j < len(aa) and aa[j] not in {"-", "X"}]
        cons_chars.append(Counter(col).most_common(1)[0][0] if col else "A")
    return by_id, seqs_all, "".join(cons_chars)


def _codon_to_aa(nt: str) -> str:
    out: list[str] = []
    n = len(nt) - (len(nt) % 3)
    for i in range(0, n, 3):
        cod = nt[i:i + 3]
        if "-" in cod or "N" in cod:
            out.append("-")
        else:
            out.append(CODON_TABLE.get(cod, "X"))
    return "".join(out)


def _make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.end_deletion_score = 0.0
    aligner.end_insertion_score = 0.0
    return aligner


ALIGNER = _make_aligner()


def _align_pdb_to_hmm_msa(pdb_seq: str, msa_aa: str) -> tuple[dict[int, int], float]:
    """Align ungapped pdb_seq against ungapped positions of msa_aa, returning {pdb_idx: 1-based hmm_col}."""
    target_chars: list[str] = []
    hmm_cols: list[int] = []
    for col0, ch in enumerate(msa_aa):
        if ch not in {"-", "*", "?"}:
            ch_use = ch if ch in MAX_ASA else "X"
            target_chars.append(ch_use)
            hmm_cols.append(col0 + 1)
    if not target_chars or len(pdb_seq) < 30:
        return {}, 0.0
    target_str = "".join(target_chars)
    query_str = "".join(c if c in MAX_ASA else "X" for c in pdb_seq)
    aln = ALIGNER.align(target_str, query_str)[0]
    t_blocks, q_blocks = aln.aligned
    mapping: dict[int, int] = {}
    ident = 0
    aligned_len = 0
    for (t0, t1), (q0, q1) in zip(t_blocks, q_blocks):
        for off in range(t1 - t0):
            ti = t0 + off
            qi = q0 + off
            mapping[qi] = hmm_cols[ti]
            aligned_len += 1
            if target_str[ti] == query_str[qi]:
                ident += 1
    pident = ident / max(aligned_len, 1)
    return mapping, pident


def map_pdb_seq_to_hmm_cols(
    pdb_seq: str,
    seq_keys: list[str],
    by_id: dict[str, str],
    seqs_all: list[tuple[str, str]],
    cons_seq: str,
) -> tuple[dict[int, int], float]:
    """Map PDB sequence positions to 1-based HMM columns cleanly."""
    for sk in seq_keys:
        if sk and sk in by_id:
            m, pid = _align_pdb_to_hmm_msa(pdb_seq, by_id[sk])
            if pid >= 0.70 or len(m) >= 50:
                return m, pid

    # Crystal / missing IDs: screen top 5-mer matches across alignment + consensus
    kmers = {pdb_seq[i:i + 5] for i in range(0, max(0, len(pdb_seq) - 4), 3)}
    candidates = [("", cons_seq)]
    if seqs_all and kmers:
        scored: list[tuple[int, str, str]] = []
        for hdr, msa_aa in seqs_all[::2]:
            ung = msa_aa.replace("-", "")
            hits = sum(1 for k in kmers if k in ung)
            if hits > 0:
                scored.append((hits, hdr, msa_aa))
        scored.sort(key=lambda x: -x[0])
        for _, hdr, msa_aa in scored[:4]:
            candidates.append((hdr, msa_aa))

    best_map: dict[int, int] = {}
    best_score = -1.0
    for _, msa_aa in candidates:
        m, pid = _align_pdb_to_hmm_msa(pdb_seq, msa_aa)
        score = len(m) * pid
        if score > best_score:
            best_score = score
            best_map = m
    best_pid = best_score / max(len(best_map), 1) if best_map else 0.0
    return best_map, float(best_pid)


# ---------------------------------------------------------------------------
# 3D Biophysics Helpers (Part 1)
# ---------------------------------------------------------------------------
def classify_ss_phi_psi(phi: float | None, psi: float | None) -> str:
    if phi is None or psi is None:
        return "Loop_Coil"
    phi_d = math.degrees(phi)
    psi_d = math.degrees(psi)
    if -115.0 <= phi_d <= -30.0 and -75.0 <= psi_d <= -10.0:
        return "Helix"
    if -180.0 <= phi_d <= -45.0 and (psi_d >= 75.0 or psi_d <= -140.0):
        return "Sheet"
    return "Loop_Coil"


def extract_pdb_biophysics(pdb_path: Path) -> tuple[str, list[dict]] | None:
    """Compute per-residue RSA, phi/psi SS class, and CA coordinates on longest chain."""
    parser = PDBParser(QUIET=True)
    try:
        struct = parser.get_structure(pdb_path.stem, str(pdb_path))
        model = struct[0]
    except Exception:
        return None

    best_chain = None
    best_n = -1
    for ch in model:
        n_aa = sum(1 for r in ch if is_aa(r, standard=True) and "CA" in r)
        if n_aa > best_n:
            best_n = n_aa
            best_chain = ch
    if best_chain is None or best_n < 30:
        return None

    try:
        sr = ShrakeRupley()
        sr.compute(best_chain, level="R")
    except Exception:
        return None

    ppb = PPBuilder()
    phi_psi_map: dict[tuple, tuple] = {}
    for pp in ppb.build_peptides(best_chain):
        pps = pp.get_phi_psi_list()
        for res, (phi, psi) in zip(pp, pps):
            phi_psi_map[res.id] = (phi, psi)

    records: list[dict] = []
    aa_chars: list[str] = []
    for res in best_chain:
        if not is_aa(res, standard=True) or "CA" not in res:
            continue
        aa1 = AA3_TO_1.get(res.get_resname().strip().upper())
        if not aa1 or aa1 not in MAX_ASA:
            continue
        rsa = min(1.0, float(getattr(res, "sasa", np.nan)) / MAX_ASA[aa1])
        phi, psi = phi_psi_map.get(res.id, (None, None))
        ss = classify_ss_phi_psi(phi, psi)
        records.append({
            "pdb_idx": len(records),
            "aa1": aa1,
            "resname_resseq": f"{aa1}{res.id[1]}",
            "rsa": rsa,
            "ss": ss,
            "is_loop": int(ss == "Loop_Coil"),
            "ca_coord": res["CA"].get_coord().astype(float),
        })
        aa_chars.append(aa1)

    if len(records) < 30:
        return None
    return "".join(aa_chars), records


def find_catalytic_ser_idx(
    pdb_seq: str,
    pdb_to_hmm: dict[int, int],
    sxxk_hmm_col: int | None,
) -> int | None:
    """Identify catalytic SXXK Serine index using HMM-anchored [ST]..[KRH] motif."""
    regex_hits = [m.start() for m in re.finditer(r"[ST]..[KRH]", pdb_seq)]
    if sxxk_hmm_col is not None and pdb_to_hmm:
        best_i = None
        best_dist = 9999
        for i in regex_hits:
            hc = pdb_to_hmm.get(i)
            if hc is not None and abs(hc - sxxk_hmm_col) < best_dist:
                best_dist = abs(hc - sxxk_hmm_col)
                best_i = i
        if best_i is not None and best_dist <= 12:
            return best_i
        # Fallback: direct HMM column match (exact or +-3)
        near_hm = [
            (abs(hc - sxxk_hmm_col), i)
            for i, hc in pdb_to_hmm.items()
            if abs(hc - sxxk_hmm_col) <= 3 and pdb_seq[i] in {"S", "T", "C"}
        ]
        if near_hm:
            near_hm.sort()
            return near_hm[0][1]
    if regex_hits:
        return regex_hits[0]
    return None


def _match_viral_pdbs_to_catalog(v26: pd.DataFrame) -> dict[str, Path]:
    """Map exact curated 2026 protein uid -> PDB path (uid_key prefix match, uvig fallback)."""
    by_uvig: dict[str, list[Path]] = defaultdict(list)
    for p in PDB_DIR.glob("IMGVR*.pdb"):
        by_uvig[p.stem.split("__")[0]].append(p)
    matched: dict[str, Path] = {}
    for _, r in v26.iterrows():
        uv = r["uvig"]
        cands = by_uvig.get(uv, [])
        if not cands:
            continue
        uk = r["uid_key"]
        exact = [p for p in cands if p.stem.startswith(uk + "_") or p.stem == uk]
        if exact:
            matched[r["uid"]] = exact[0]
        elif len(cands) == 1:
            matched[r["uid"]] = cands[0]
    return matched


def select_representative_viral_pdbs(v26: pd.DataFrame, max_per_fam: int = 40) -> dict[str, list[dict]]:
    """Select diverse representative curated 2026 viral PDBs per family for full-atom 3D mapping."""
    pdb_by_uid = _match_viral_pdbs_to_catalog(v26)

    selected: dict[str, list[dict]] = {}
    for fam in FAMILIES:
        sub = v26[(v26["winning_acc"] == fam) & (v26["uid"].isin(pdb_by_uid))].copy()
        if sub.empty:
            selected[fam] = []
            continue
        sub["is_hyd"] = (sub["functional_tier"] == "Putative_Hydrolytic_BetaLactamase").astype(int)
        sub["is_highconf"] = sub["viral_origin_tier"].isin(HIGH_CONF_TIERS).astype(int)
        sub = sub.sort_values(
            ["is_hyd", "is_highconf", "hmm_cov", "dom_score"],
            ascending=[False, False, False, False],
        )
        # Prioritize up to 15 Putative_Hydrolytic_BetaLactamase + cluster-diverse PBP_or_PRE
        chosen_rows = []
        seen_clusters = set()
        for _, r in sub[sub["is_hyd"] == 1].iterrows():
            if len(chosen_rows) >= 15:
                break
            chosen_rows.append(r)
            seen_clusters.add(r["mmseqs85_cluster_id"])
        for _, r in sub.iterrows():
            if len(chosen_rows) >= max_per_fam:
                break
            cid = r["mmseqs85_cluster_id"]
            if cid not in seen_clusters:
                chosen_rows.append(r)
                seen_clusters.add(cid)
        for _, r in sub.iterrows():
            if len(chosen_rows) >= max_per_fam:
                break
            if r["uid"] not in {x["uid"] for x in chosen_rows}:
                chosen_rows.append(r)
        selected[fam] = [
            {"pdb_path": pdb_by_uid[r["uid"]], "uid": r["uid"], "uvig": r["uvig"], "tier": r["functional_tier"]}
            for r in chosen_rows
        ]
    return selected


def run_part1_3d_biophysics(v26: pd.DataFrame, metrics_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Map Ref crystal & 2026 viral ESMFold 3D biophysics onto HMM selection columns."""
    print("[Part 1] Mapping 3D crystal + viral ESMFold biophysics onto 2026 HMM selection columns...")
    has_meme = "meme_min_p" in metrics_df.columns and metrics_df["meme_min_p"].notna().any()

    # Load translated codon alignments and Motif1_SXXK Serine anchor columns
    fam_alns = {}
    sxxk_col_map: dict[str, int | None] = {}
    for fam in FAMILIES:
        fna = SEL_DIR / "codon_alignments" / f"{fam}.full.hmm_codon.fna"
        if fna.exists():
            fam_alns[fam] = translate_codon_aln(fna)
        else:
            fam_alns[fam] = ({}, [], "")
        sub_s = metrics_df[(metrics_df["family"] == fam) & (metrics_df["motif_role"] == "Motif1_SXXK")]
        sxxk_col_map[fam] = int(sub_s["hmm_col"].min()) if len(sub_s) else None

    # 1A. Ref crystal PDB mapping
    ref_site_records: dict[tuple[str, int], dict] = {}
    for fam, (fname, label) in REF_PDBS.items():
        p = REF_PDB_DIR / fname
        if not p.exists():
            continue
        parsed = extract_pdb_biophysics(p)
        if parsed is None:
            continue
        pdb_seq, res_recs = parsed
        by_id, seqs_all, cons_seq = fam_alns[fam]
        mapping, pid = map_pdb_seq_to_hmm_cols(pdb_seq, [], by_id, seqs_all, cons_seq)
        ser_idx = find_catalytic_ser_idx(pdb_seq, mapping, sxxk_col_map[fam])
        if ser_idx is None:
            print(f"  WARN: could not locate SXXK Ser in crystal {fname}")
            continue
        ser_coord = res_recs[ser_idx]["ca_coord"]
        mapped_cols = 0
        for rec in res_recs:
            hcol = mapping.get(rec["pdb_idx"])
            if hcol is None:
                continue
            dist = float(np.linalg.norm(rec["ca_coord"] - ser_coord))
            ref_site_records[(fam, hcol)] = {
                "ref_pdb_id": label,
                "ref_resname_resseq": rec["resname_resseq"],
                "ref_rsa": rec["rsa"],
                "ref_dist_ser_A": dist,
                "ref_ss": rec["ss"],
                "ref_is_loop": rec["is_loop"],
            }
            mapped_cols += 1
        print(f"  [Crystal {fam}] {label}: mapped {mapped_cols} HMM cols (id={pid:.2f}, ser={res_recs[ser_idx]['resname_resseq']})")

    # 1B. Viral ESMFold ensemble mapping (up to 40 curated PDBs per family)
    rep_viral = select_representative_viral_pdbs(v26, max_per_fam=40)
    viral_col_rsa: dict[tuple[str, int], list[float]] = defaultdict(list)
    viral_col_dist: dict[tuple[str, int], list[float]] = defaultdict(list)
    viral_col_loop: dict[tuple[str, int], list[int]] = defaultdict(list)

    for fam in FAMILIES:
        by_id, seqs_all, cons_seq = fam_alns[fam]
        n_ok = 0
        for item in rep_viral.get(fam, []):
            parsed = extract_pdb_biophysics(item["pdb_path"])
            if parsed is None:
                continue
            pdb_seq, res_recs = parsed
            mapping, pid = map_pdb_seq_to_hmm_cols(
                pdb_seq, [item["uid"], item["uvig"]], by_id, seqs_all, cons_seq
            )
            if len(mapping) < 50:
                continue
            ser_idx = find_catalytic_ser_idx(pdb_seq, mapping, sxxk_col_map[fam])
            if ser_idx is None:
                continue
            ser_coord = res_recs[ser_idx]["ca_coord"]
            for rec in res_recs:
                hcol = mapping.get(rec["pdb_idx"])
                if hcol is None:
                    continue
                key = (fam, hcol)
                viral_col_rsa[key].append(rec["rsa"])
                viral_col_dist[key].append(float(np.linalg.norm(rec["ca_coord"] - ser_coord)))
                viral_col_loop[key].append(rec["is_loop"])
            n_ok += 1
        print(f"  [Viral ESMFold {fam}] mapped {n_ok}/{len(rep_viral.get(fam, []))} representative curated 2026 PDBs")

    # Merge onto PerSite_Domain_Selection_Metrics_2026.tsv
    merged = metrics_df.copy()
    for col in ["ref_pdb_id", "ref_resname_resseq", "ref_ss"]:
        merged[col] = [ref_site_records.get((r["family"], int(r["hmm_col"])), {}).get(col, "") for _, r in merged.iterrows()]
    for col in ["ref_rsa", "ref_dist_ser_A", "ref_is_loop"]:
        merged[col] = [ref_site_records.get((r["family"], int(r["hmm_col"])), {}).get(col, np.nan) for _, r in merged.iterrows()]

    v_n, v_rsa_mean, v_rsa_med, v_dist_mean, v_dist_med, v_loop_pct = [], [], [], [], [], []
    for _, r in merged.iterrows():
        key = (r["family"], int(r["hmm_col"]))
        rsas = viral_col_rsa.get(key, [])
        dists = viral_col_dist.get(key, [])
        loops = viral_col_loop.get(key, [])
        v_n.append(len(rsas))
        v_rsa_mean.append(float(np.mean(rsas)) if rsas else np.nan)
        v_rsa_med.append(float(np.median(rsas)) if rsas else np.nan)
        v_dist_mean.append(float(np.mean(dists)) if dists else np.nan)
        v_dist_med.append(float(np.median(dists)) if dists else np.nan)
        v_loop_pct.append(float(np.mean(loops) * 100.0) if loops else np.nan)
    merged["n_viral_pdbs_mapped"] = v_n
    merged["viral_rsa_mean"] = v_rsa_mean
    merged["viral_rsa_median"] = v_rsa_med
    merged["viral_dist_ser_mean"] = v_dist_mean
    merged["viral_dist_ser_median"] = v_dist_med
    merged["viral_loop_pct"] = v_loop_pct

    # Classify each HMM column per family into 1_Catalytic_Motif_Core / 2_High_Diversification_Sites / 3_Background_Core
    merged["site_class"] = "4_Low_Occupancy"
    crit_used = "meme_min_p<0.05 & omega_ng>=0.50" if has_meme else "family_topQ25_omega_ng & >=median_H_aa"
    for fam in FAMILIES:
        core_mask = (merged["family"] == fam) & (merged["occupancy"] >= 0.30)
        bg_vals = merged.loc[core_mask & (~merged["in_catalytic_motif_core"])]
        if bg_vals.empty:
            continue
        q75_omega = float(np.nanquantile(bg_vals["codon_omega_ng"], 0.75))
        q50_haa = float(np.nanquantile(bg_vals["H_aa"], 0.50))

        m_cat = (merged["family"] == fam) & (merged["in_catalytic_motif_core"])
        if has_meme:
            m_div = core_mask & (~merged["in_catalytic_motif_core"]) & (
                ((merged["meme_min_p"] < 0.05) & (merged["codon_omega_ng"] >= 0.50))
                | ((merged["codon_omega_ng"] >= q75_omega) & (merged["H_aa"] >= q50_haa))
            )
        else:
            m_div = core_mask & (~merged["in_catalytic_motif_core"]) & (
                (merged["codon_omega_ng"] >= q75_omega) & (merged["H_aa"] >= q50_haa)
            )
        m_bg = core_mask & (~merged["in_catalytic_motif_core"]) & (~m_div)
        merged.loc[m_bg, "site_class"] = "3_Background_Core"
        merged.loc[m_div, "site_class"] = "2_High_Diversification_Sites"
        merged.loc[m_cat, "site_class"] = "1_Catalytic_Motif_Core"

    merged["diversification_criterion"] = crit_used
    merged.to_csv(OUT_DIR / "PerSite_3D_Biophysics_2026.tsv", sep="\t", index=False)

    # Summarize structural localization tests (Ref Crystal + Viral ESMFold per family & pooled)
    test_rows = []
    sources = [("Ref_Crystal", "ref_dist_ser_A", "ref_rsa", "ref_is_loop", False)]
    sources.append(("Viral_ESMFold_2026", "viral_dist_ser_median", "viral_rsa_median", "viral_loop_pct", True))

    for src_name, dist_c, rsa_c, loop_c, loop_is_pct in sources:
        fams_to_run = [f for f in FAMILIES if (src_name != "Ref_Crystal" or f in REF_PDBS)] + [f"All_{src_name}_Pooled"]
        for fam_lbl in fams_to_run:
            if fam_lbl.startswith("All_"):
                sub = merged[merged[dist_c].notna() & merged["site_class"].isin(["1_Catalytic_Motif_Core", "2_High_Diversification_Sites", "3_Background_Core"])].copy()
            else:
                sub = merged[(merged["family"] == fam_lbl) & merged[dist_c].notna() & merged["site_class"].isin(["1_Catalytic_Motif_Core", "2_High_Diversification_Sites", "3_Background_Core"])].copy()
            if sub.empty:
                continue

            cat = sub[sub["site_class"] == "1_Catalytic_Motif_Core"]
            div = sub[sub["site_class"] == "2_High_Diversification_Sites"]
            bg = sub[sub["site_class"] == "3_Background_Core"]

            div_loop = (div[loop_c] >= 50.0) if loop_is_pct else (div[loop_c] == 1)
            bg_loop = (bg[loop_c] >= 50.0) if loop_is_pct else (bg[loop_c] == 1)
            cat_loop = (cat[loop_c] >= 50.0) if loop_is_pct else (cat[loop_c] == 1)

            # Fisher OR: high diversification sites distal (>=10 A & RSA >= 0.20) vs background
            div_distal = (div[dist_c] >= 10.0) & (div[rsa_c] >= 0.20)
            bg_distal = (bg[dist_c] >= 10.0) & (bg[rsa_c] >= 0.20)
            or_distal, p_distal = fisher_exact([[int(div_distal.sum()), int((~div_distal).sum())], [int(bg_distal.sum()), int((~bg_distal).sum())]]) if len(div) and len(bg) else (np.nan, np.nan)
            or_loop, p_loop = fisher_exact([[int(div_loop.sum()), int((~div_loop).sum())], [int(bg_loop.sum()), int((~bg_loop).sum())]]) if len(div) and len(bg) else (np.nan, np.nan)

            # Spatial partition test: Active center (<10 A) vs Distal solvent-exposed (>=10 A & RSA >= 0.20)
            act_ctr = sub[sub[dist_c] < 10.0]
            exp_surf = sub[(sub[dist_c] >= 10.0) & (sub[rsa_c] >= 0.20)]

            test_rows.append({
                "structure_source": src_name,
                "family": fam_lbl,
                "n_mapped_sites": len(sub),
                "n_cat": len(cat),
                "n_div": len(div),
                "n_bg": len(bg),
                "cat_median_dist_A": round(float(cat[dist_c].median()), 3) if len(cat) else np.nan,
                "div_median_dist_A": round(float(div[dist_c].median()), 3) if len(div) else np.nan,
                "bg_median_dist_A": round(float(bg[dist_c].median()), 3) if len(bg) else np.nan,
                "cat_median_RSA": round(float(cat[rsa_c].median()), 4) if len(cat) else np.nan,
                "div_median_RSA": round(float(div[rsa_c].median()), 4) if len(div) else np.nan,
                "bg_median_RSA": round(float(bg[rsa_c].median()), 4) if len(bg) else np.nan,
                "cat_loop_pct": round(float(cat_loop.mean() * 100.0), 1) if len(cat) else np.nan,
                "div_loop_pct": round(float(div_loop.mean() * 100.0), 1) if len(div) else np.nan,
                "bg_loop_pct": round(float(bg_loop.mean() * 100.0), 1) if len(bg) else np.nan,
                "MW_p_div_gt_cat_dist": safe_mw(div[dist_c].values, cat[dist_c].values, "greater"),
                "MW_p_div_gt_cat_RSA": safe_mw(div[rsa_c].values, cat[rsa_c].values, "greater"),
                "MW_p_div_gt_bg_dist": safe_mw(div[dist_c].values, bg[dist_c].values, "greater"),
                "MW_p_div_gt_bg_RSA": safe_mw(div[rsa_c].values, bg[rsa_c].values, "greater"),
                "Fisher_OR_div_vs_bg_distal_exposed": round(float(or_distal), 3) if pd.notna(or_distal) else np.nan,
                "Fisher_p_div_vs_bg_distal_exposed": float(p_distal) if pd.notna(p_distal) else np.nan,
                "Fisher_OR_div_vs_bg_loop": round(float(or_loop), 3) if pd.notna(or_loop) else np.nan,
                "Fisher_p_div_vs_bg_loop": float(p_loop) if pd.notna(p_loop) else np.nan,
                "n_active_center_lt10A": len(act_ctr),
                "n_distal_exposed_ge10A_RSA20": len(exp_surf),
                "active_center_median_omega_ng": round(float(act_ctr["codon_omega_ng"].median()), 4) if len(act_ctr) else np.nan,
                "distal_exposed_median_omega_ng": round(float(exp_surf["codon_omega_ng"].median()), 4) if len(exp_surf) else np.nan,
                "MW_p_distal_gt_active_omega_ng": safe_mw(exp_surf["codon_omega_ng"].values, act_ctr["codon_omega_ng"].values, "greater"),
                "active_center_median_H_aa": round(float(act_ctr["H_aa"].median()), 4) if len(act_ctr) else np.nan,
                "distal_exposed_median_H_aa": round(float(exp_surf["H_aa"].median()), 4) if len(exp_surf) else np.nan,
                "MW_p_distal_gt_active_H_aa": safe_mw(exp_surf["H_aa"].values, act_ctr["H_aa"].values, "greater"),
            })

    df_loc = pd.DataFrame(test_rows)
    df_loc.to_csv(OUT_DIR / "Table_3D_Structural_Localization_2026.tsv", sep="\t", index=False)
    return merged, df_loc


# ---------------------------------------------------------------------------
# Fast 3D Geometry + Phylogenetic Diversity (Part 2)
# ---------------------------------------------------------------------------
def fast_pdb_ca_features(pdb_path: Path) -> tuple[int, float, float] | None:
    """Fast CA-only extraction of residue count, length-normalized Rg, and asphericity."""
    coords = []
    seen = set()
    try:
        with open(pdb_path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
                    key = (line[21], line[22:27])
                    if key not in seen:
                        seen.add(key)
                        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    except Exception:
        return None
    if len(coords) < 20:
        return None
    pts = np.asarray(coords, dtype=float)
    center = pts.mean(axis=0)
    d = pts - center
    rg = float(np.sqrt(np.mean(np.sum(d * d, axis=1))))
    norm_rg = rg / (len(pts) ** (1.0 / 3.0))
    cov = np.cov(d.T)
    eigvals = np.maximum(np.linalg.eigvalsh(cov), 0.0)
    lam_sum = float(np.sum(eigvals))
    if lam_sum > 0:
        asph = 1.5 * float(np.sum((eigvals - lam_sum / 3.0) ** 2)) / (lam_sum * lam_sum)
    else:
        asph = 0.0
    return len(pts), float(norm_rg), float(asph)


def rarefied_tree_diversity(
    tree: Tree,
    viral_tips: list[str],
    bac_tips: list[str],
    n_iter: int = 200,
    seed: int = 42,
) -> dict[str, float]:
    """Subsample k = min(|V|, |B|, 50) leaves and compute rarefied Faith's PD and MPD."""
    rng = np.random.default_rng(seed)
    nv, nb = len(viral_tips), len(bac_tips)
    k = min(nv, nb, 50)
    if k < 4:
        return {
            "rarefied_k": 0, "viral_PD_mean": np.nan, "viral_PD_sd": np.nan,
            "bac_PD_mean": np.nan, "bac_PD_sd": np.nan, "PD_MW_p": np.nan,
            "viral_MPD_mean": np.nan, "bac_MPD_mean": np.nan, "MPD_MW_p": np.nan,
        }

    # Precompute leaf-to-ancestors path dictionary for fast exact Faith's PD
    leaf_nodes = {leaf.name: leaf for leaf in tree.leaves()}
    anc_dist: dict[str, tuple[list[tuple[int, float]], dict[int, float]]] = {}
    for nm in set(viral_tips) | set(bac_tips):
        nd = leaf_nodes.get(nm)
        if nd is None:
            continue
        path = []
        dmap = {}
        d_acc = 0.0
        cur = nd
        while cur is not None:
            dmap[id(cur)] = d_acc
            if cur.up is not None:
                path.append((id(cur), float(cur.dist or 0.0)))
                d_acc += float(cur.dist or 0.0)
            cur = cur.up
        anc_dist[nm] = (path, dmap)

    def _pd_mpd(sample_names: np.ndarray) -> tuple[float, float]:
        # Faith's PD = sum of unique branches above the MRCA of the sample
        # Approximate exact union up to root (difference cancels identically when comparing clades or prune to MRCA)
        lca_node = tree.common_ancestor([leaf_nodes[nm] for nm in sample_names])
        stop_id = id(lca_node)
        seen_edges = set()
        pd_sum = 0.0
        for nm in sample_names:
            for eid, bl in anc_dist[nm][0]:
                if eid == stop_id:
                    break
                if eid not in seen_edges:
                    seen_edges.add(eid)
                    pd_sum += bl
        # Exact MPD from root-to-node depths & pairwise LCA depths via fast ancestor sweep
        pw_dists = []
        n_s = len(sample_names)
        for i in range(n_s):
            _, dmap_i = anc_dist[sample_names[i]]
            for j in range(i + 1, n_s):
                _, dmap_j = anc_dist[sample_names[j]]
                for anc_id, di in dmap_i.items():
                    if anc_id in dmap_j:
                        pw_dists.append(di + dmap_j[anc_id])
                        break
        return pd_sum, float(np.mean(pw_dists)) if pw_dists else np.nan

    v_pds, b_pds, v_mpds, b_mpds = [], [], [], []
    v_arr = np.asarray([t for t in viral_tips if t in anc_dist])
    b_arr = np.asarray([t for t in bac_tips if t in anc_dist])
    k = min(len(v_arr), len(b_arr), 50)
    if k < 4:
        return {
            "rarefied_k": 0, "viral_PD_mean": np.nan, "viral_PD_sd": np.nan,
            "bac_PD_mean": np.nan, "bac_PD_sd": np.nan, "PD_MW_p": np.nan,
            "viral_MPD_mean": np.nan, "bac_MPD_mean": np.nan, "MPD_MW_p": np.nan,
        }
    for _ in range(n_iter):
        sv = rng.choice(v_arr, size=k, replace=False)
        sb = rng.choice(b_arr, size=k, replace=False)
        vp, vm = _pd_mpd(sv)
        bp, bm = _pd_mpd(sb)
        v_pds.append(vp)
        b_pds.append(bp)
        v_mpds.append(vm)
        b_mpds.append(bm)

    return {
        "rarefied_k": int(k),
        "viral_PD_mean": round(float(np.mean(v_pds)), 4),
        "viral_PD_sd": round(float(np.std(v_pds, ddof=1)), 4),
        "bac_PD_mean": round(float(np.mean(b_pds)), 4),
        "bac_PD_sd": round(float(np.std(b_pds, ddof=1)), 4),
        "PD_MW_p": safe_mw(np.asarray(v_pds), np.asarray(b_pds)),
        "viral_MPD_mean": round(float(np.mean(v_mpds)), 4),
        "bac_MPD_mean": round(float(np.mean(b_mpds)), 4),
        "MPD_MW_p": safe_mw(np.asarray(v_mpds), np.asarray(b_mpds)),
    }


def _classify_bac_pdb_from_seq(pdb_path: Path, fam: str) -> str:
    """Classify a bacterial PDB sequence into Putative_Hydrolytic_BetaLactamase vs PBP_or_PRE_Homolog."""
    aa_chars: list[str] = []
    seen = set()
    try:
        with open(pdb_path, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM  ") and line[12:16].strip() == "CA":
                    key = (line[21], line[22:27])
                    if key not in seen:
                        seen.add(key)
                        aa_chars.append(AA3_TO_1.get(line[17:20].strip(), "X"))
    except Exception:
        return "PBP_or_PRE_Homolog"
    seq = "".join(aa_chars)
    if fam == "PF00144":
        ok = bool(re.search(r"S[A-Z]{2}K", seq) and re.search(r"Y[AS]N", seq) and re.search(r"K[TS]G", seq))
    elif fam == "PF00905":
        m_span = re.search(r"_(\d+)-(\d+)$", pdb_path.stem)
        end_pos = int(m_span.group(2)) if m_span else len(seq)
        ok = bool(
            len(seq) < 400
            and end_pos < 400
            and re.search(r"S[TFASV][FYWLIVM]K", seq)
            and re.search(r"[YF]GN", seq)
            and re.search(r"[FYW]WL|W[A-Z][EDQNG]", seq)
            and re.search(r"K[TS]G", seq)
        )
    elif fam == "PF13354":
        ok = bool(
            re.search(r"S[A-Z]{2}K", seq)
            and re.search(r"[SY][DN]N", seq)
            and re.search(r"E[A-Z]{2}[LIVMFCA]N|TE[A-Z]{2}[MLIVFA]", seq)
            and re.search(r"[KR][TS][GA]", seq)
        )
    else:
        ok = False
    return "Putative_Hydrolytic_BetaLactamase" if ok else "PBP_or_PRE_Homolog"


def _load_bac_catalog_lookup() -> dict[str, tuple[str, str]]:
    """Index Final_Bacterial_TripleSource_Catalog.tsv by normalized orig_uid for PDB lookup."""
    if not BAC_CAT.is_file():
        return {}
    bcat = pd.read_csv(BAC_CAT, sep="\t", usecols=["orig_uid", "winning_acc", "functional_tier"], low_memory=False)
    lookup: dict[str, tuple[str, str]] = {}
    for _, r in bcat.iterrows():
        u = str(r["orig_uid"])
        val = (str(r["winning_acc"]), str(r["functional_tier"]))
        lookup[u] = val
        u_base = re.sub(r"\.p\d+$", "", u)
        lookup[re.sub(r"[^A-Za-z0-9]+", "_", u_base)] = val
        u_nover = re.sub(r"^([A-Z]+_[A-Z0-9]+)\.\d+\.", r"\1_", u_base)
        lookup[re.sub(r"[^A-Za-z0-9]+", "_", u_nover)] = val
    return lookup


def run_part2_structural_diversity(
    v26_meta: pd.DataFrame,
    mn: pd.DataFrame,
    me: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare Viral vs Bacterial 3D fold descriptors and rarefied phylogenetic diversity."""
    print("[Part 2] Computing fast CA geometry across all 1,507 PDBs & maxcc.nwk rarefied diversity...")
    tree = Tree(str(MAXCC_NWK), parser=0)
    tips_strip_map = {re.sub(r"_A$", "", lf.name): lf.name for lf in tree.leaves()}

    v_pdb_map = {str(p): uid for uid, p in _match_viral_pdbs_to_catalog(v26_meta).items()}
    v_uid_lookup = v26_meta.set_index("uid")
    mn_bac = mn[mn["taxonomy"] == "BACTERIAL"].set_index("Id")
    bcat_lookup = _load_bac_catalog_lookup()

    rows = []
    for p in sorted(PDB_DIR.glob("*.pdb")):
        stem = p.stem
        feats = fast_pdb_ca_features(p)
        if feats is None:
            continue
        n_res, norm_rg, asph = feats
        tip_name = tips_strip_map.get(stem, "")
        if stem.startswith("IMGVR"):
            matched_uid = v_pdb_map.get(str(p))
            in_2026 = matched_uid is not None
            if in_2026:
                r = v_uid_lookup.loc[matched_uid]
                fam = str(r["winning_acc"])
                ftier = str(r["functional_tier"])
                otier = str(r["viral_origin_tier"])
                vphy = str(r["viral_phylum"])
                hphy = str(r["host_phylum"])
                hmet = str(r["host_method"])
            else:
                fam, ftier, otier, vphy, hphy, hmet = "FilteredOut_2026", "FilteredOut_2026", "FilteredOut_2026", "Unassigned", "Unassigned", "NA"
            rows.append({
                "pdb_stem": stem,
                "domain_group": "Viral",
                "in_curated_2026": in_2026,
                "sbl_family": fam,
                "functional_tier": ftier,
                "viral_origin_tier": otier,
                "viral_phylum": vphy,
                "host_phylum": hphy,
                "host_method": hmet,
                "reliable_host": hmet in RELIABLE_HOST_METHODS,
                "high_conf_viral": otier in HIGH_CONF_TIERS,
                "tree_tip": tip_name,
                "n_ca_res": n_res,
                "norm_rg": norm_rg,
                "asphericity": asph,
            })
        else:
            if stem in CRYSTAL_BAC_INFO:
                fam, ftier = CRYSTAL_BAC_INFO[stem]
            else:
                stem_base = re.sub(r"_\d+-\d+$", "", stem)
                stem_no1 = re.sub(r"_1$", "", stem_base)
                b_hit = bcat_lookup.get(stem) or bcat_lookup.get(stem_base) or bcat_lookup.get(stem_no1)
                if b_hit is not None:
                    fam, ftier = b_hit
                elif stem in mn_bac.index:
                    fam = str(mn_bac.loc[stem, "pfam"])
                    ftier = _classify_bac_pdb_from_seq(p, fam)
                else:
                    fam, ftier = "Unknown", "PBP_or_PRE_Homolog"
            rows.append({
                "pdb_stem": stem,
                "domain_group": "Bacterial",
                "in_curated_2026": True,
                "sbl_family": fam,
                "functional_tier": ftier,
                "viral_origin_tier": "Bacterial_Reference",
                "viral_phylum": "Bacteria",
                "host_phylum": "Bacteria",
                "host_method": "Bacterial_Reference",
                "reliable_host": True,
                "high_conf_viral": False,
                "tree_tip": tip_name,
                "n_ca_res": n_res,
                "norm_rg": norm_rg,
                "asphericity": asph,
            })

    df_pdb = pd.DataFrame(rows)
    df_pdb.to_csv(OUT_DIR / "PDB_3D_Structural_Catalog_2026.tsv", sep="\t", index=False)

    # Evaluate 10 controlled strata
    cur_v = (df_pdb["domain_group"] == "Viral") & (df_pdb["in_curated_2026"])
    all_v = df_pdb["domain_group"] == "Viral"
    all_b = df_pdb["domain_group"] == "Bacterial"

    strata_defs = [
        ("1_Global_Baseline_Unfiltered", all_v, all_b),
        ("2_Global_Curated_2026", cur_v, all_b),
        ("3_Within_PF00144", cur_v & (df_pdb["sbl_family"] == "PF00144"), all_b & (df_pdb["sbl_family"] == "PF00144")),
        ("3_Within_PF00905", cur_v & (df_pdb["sbl_family"] == "PF00905"), all_b & (df_pdb["sbl_family"] == "PF00905")),
        ("3_Within_PF13354", cur_v & (df_pdb["sbl_family"] == "PF13354"), all_b & (df_pdb["sbl_family"] == "PF13354")),
        ("3_Within_PF00768_ViralOnly", cur_v & (df_pdb["sbl_family"] == "PF00768"), pd.Series(False, index=df_pdb.index)),
        ("4_Tier_PBP_or_PRE_Homolog", cur_v & (df_pdb["functional_tier"] == "PBP_or_PRE_Homolog"), all_b & (df_pdb["functional_tier"] == "PBP_or_PRE_Homolog")),
        ("4_Tier_Putative_Hydrolytic_BetaLactamase", cur_v & (df_pdb["functional_tier"] == "Putative_Hydrolytic_BetaLactamase"), all_b & (df_pdb["functional_tier"] == "Putative_Hydrolytic_BetaLactamase")),
        ("5_ViralOrigin_Tier1_Tier2_HighConf", cur_v & df_pdb["high_conf_viral"], all_b),
        ("6_ReliableHost_Isolate_CRISPR", cur_v & df_pdb["reliable_host"], all_b),
    ]

    div_rows = []
    for sname, m_v, m_b in strata_defs:
        vg = df_pdb[m_v]
        bg = df_pdb[m_b]
        vtips = [t for t in vg["tree_tip"] if t]
        btips = [t for t in bg["tree_tip"] if t]
        phylo = rarefied_tree_diversity(tree, vtips, btips) if (len(vtips) >= 4 and len(btips) >= 4) else {
            "rarefied_k": 0, "viral_PD_mean": np.nan, "viral_PD_sd": np.nan,
            "bac_PD_mean": np.nan, "bac_PD_sd": np.nan, "PD_MW_p": np.nan,
            "viral_MPD_mean": np.nan, "bac_MPD_mean": np.nan, "MPD_MW_p": np.nan,
        }
        rg_p = safe_mw(vg["norm_rg"].values, bg["norm_rg"].values) if len(vg) and len(bg) else np.nan
        asp_p = safe_mw(vg["asphericity"].values, bg["asphericity"].values) if len(vg) and len(bg) else np.nan
        lev_rg_p = float(levene(vg["norm_rg"].values, bg["norm_rg"].values).pvalue) if len(vg) >= 3 and len(bg) >= 3 else np.nan

        div_rows.append({
            "stratum": sname,
            "n_viral_pdbs": len(vg),
            "n_bac_pdbs": len(bg),
            "n_viral_tree_tips": len(vtips),
            "n_bac_tree_tips": len(btips),
            **phylo,
            "viral_norm_rg_median": round(float(vg["norm_rg"].median()), 4) if len(vg) else np.nan,
            "viral_norm_rg_iqr": f"{vg['norm_rg'].quantile(0.25):.4f}-{vg['norm_rg'].quantile(0.75):.4f}" if len(vg) else "",
            "bac_norm_rg_median": round(float(bg["norm_rg"].median()), 4) if len(bg) else np.nan,
            "bac_norm_rg_iqr": f"{bg['norm_rg'].quantile(0.25):.4f}-{bg['norm_rg'].quantile(0.75):.4f}" if len(bg) else "",
            "norm_rg_MW_p": rg_p,
            "norm_rg_Levene_var_p": lev_rg_p,
            "viral_asphericity_median": round(float(vg["asphericity"].median()), 4) if len(vg) else np.nan,
            "bac_asphericity_median": round(float(bg["asphericity"].median()), 4) if len(bg) else np.nan,
            "asphericity_MW_p": asp_p,
        })

    df_div = pd.DataFrame(div_rows)
    df_div.to_csv(OUT_DIR / "Table_Structural_Diversity_Controlled_2026.tsv", sep="\t", index=False)
    print(f"  [Part 2 Done] 2026 curated viral PDBs={int(cur_v.sum())}, baseline viral={int(all_v.sum())}, bacterial={int(all_b.sum())}")
    return df_pdb, df_div


# ---------------------------------------------------------------------------
# Confounder-Controlled Ecological Stratification (Part 3)
# ---------------------------------------------------------------------------
def _evaluate_cmh_battery(
    df_sub: pd.DataFrame,
    features: list[tuple[str, str, np.ndarray]],
    habitat_col: str,
    cohort_name: str,
    min_hab_n: int = 20,
    min_feat_n: int = 20,
) -> list[dict]:
    """Run global Fisher + BH-FDR and multi-strata CMH controls over one ecological cohort."""
    valid_habs = [
        h for h, cnt in df_sub[habitat_col].value_counts().items()
        if h not in {"Unclassified", "Unclassfied", "NA", ""} and cnt >= min_hab_n
    ]
    raw_hits = []
    for feat_name, feat_cat, feat_mask in features:
        if int(feat_mask.sum()) < min_feat_n:
            continue
        for hab in valid_habs:
            h_mask = (df_sub[habitat_col] == hab).values
            a = int(np.sum(feat_mask & h_mask))
            b = int(np.sum(feat_mask & ~h_mask))
            c = int(np.sum(~feat_mask & h_mask))
            d = int(np.sum(~feat_mask & ~h_mask))
            if a < 3:
                continue
            or_val, p_val = fisher_exact([[a, b], [c, d]], alternative="greater")
            raw_hits.append({
                "cohort": cohort_name,
                "feature": feat_name,
                "feature_category": feat_cat,
                "habitat": hab,
                "n_total": len(df_sub),
                "a_feat_in_hab": a,
                "b_feat_other": b,
                "n_feat": a + b,
                "n_hab": int(h_mask.sum()),
                "pct_in_hab": round(100.0 * a / max(a + b, 1), 2),
                "pct_hab_among_nonfeat": round(100.0 * c / max(c + d, 1), 2),
                "global_OR": round(float(or_val), 3) if np.isfinite(or_val) else np.nan,
                "global_p": float(p_val),
                "_fmask": feat_mask,
                "_hmask": h_mask,
            })

    if not raw_hits:
        return []
    padjs = bh_fdr([r["global_p"] for r in raw_hits])
    for r, adj in zip(raw_hits, padjs):
        r["global_BH_padj"] = float(adj)

    # Keep all globally significant enrichments (BH < 0.05) + any nominal hits (p < 0.01) for audit transparency
    selected = [r for r in raw_hits if r["global_BH_padj"] < 0.05 or r["global_p"] < 0.01]
    selected.sort(key=lambda x: x["global_p"])

    fams = df_sub["winning_acc"].astype(str).values
    vphys = df_sub["viral_phylum"].astype(str).values
    otiers = df_sub["viral_origin_tier"].astype(str).values
    hphys = df_sub["host_phylum"].astype(str).values
    hmets = df_sub["host_method"].astype(str).values

    fam_x_phy = np.array([f"{x}|{y}" for x, y in zip(fams, vphys)])
    fam_x_tier = np.array([f"{x}|{y}" for x, y in zip(fams, otiers)])
    m_hc = np.isin(otiers, list(HIGH_CONF_TIERS))
    m_uro = vphys == "Uroviricota"
    m_rel = np.isin(hmets, list(RELIABLE_HOST_METHODS)) & (hphys != "Unassigned")

    out_rows = []
    for r in selected:
        fm = r.pop("_fmask")
        hm = r.pop("_hmask")
        cmh_fam_or, cmh_fam_p = cmh_test_stratified(fm, hm, fams)
        cmh_phy_or, cmh_phy_p = cmh_test_stratified(fm, hm, vphys)
        cmh_fp_or, cmh_fp_p = cmh_test_stratified(fm, hm, fam_x_phy)
        cmh_ft_or, cmh_ft_p = cmh_test_stratified(fm, hm, fam_x_tier)
        cmh_host_or, cmh_host_p = cmh_test_stratified(fm[m_rel], hm[m_rel], hphys[m_rel]) if int(m_rel.sum()) >= 20 else (np.nan, np.nan)

        # Tier1/Tier2 subset exact test
        if int(m_hc.sum()) >= 20 and int(np.sum(fm & m_hc)) >= 2 and int(np.sum(hm & m_hc)) >= 5:
            hc_or, hc_p = fisher_exact([
                [int(np.sum(fm & hm & m_hc)), int(np.sum(fm & (~hm) & m_hc))],
                [int(np.sum((~fm) & hm & m_hc)), int(np.sum((~fm) & (~hm) & m_hc))],
            ], alternative="greater")
        else:
            hc_or, hc_p = np.nan, np.nan

        # Uroviricota-only exact test
        if int(m_uro.sum()) >= 20 and int(np.sum(fm & m_uro)) >= 2 and int(np.sum(hm & m_uro)) >= 5:
            uro_or, uro_p = fisher_exact([
                [int(np.sum(fm & hm & m_uro)), int(np.sum(fm & (~hm) & m_uro))],
                [int(np.sum((~fm) & hm & m_uro)), int(np.sum((~fm) & (~hm) & m_uro))],
            ], alternative="greater")
        else:
            uro_or, uro_p = np.nan, np.nan

        # Reliable host subset exact test
        if int(m_rel.sum()) >= 20 and int(np.sum(fm & m_rel)) >= 2 and int(np.sum(hm & m_rel)) >= 5:
            rh_or, rh_p = fisher_exact([
                [int(np.sum(fm & hm & m_rel)), int(np.sum(fm & (~hm) & m_rel))],
                [int(np.sum((~fm) & hm & m_rel)), int(np.sum((~fm) & (~hm) & m_rel))],
            ], alternative="greater")
        else:
            rh_or, rh_p = np.nan, np.nan

        out_rows.append({
            **r,
            "CMH_OR_ctrl_Family": round(cmh_fam_or, 3) if pd.notna(cmh_fam_or) else np.nan,
            "CMH_p_ctrl_Family": cmh_fam_p,
            "CMH_OR_ctrl_ViralPhylum": round(cmh_phy_or, 3) if pd.notna(cmh_phy_or) else np.nan,
            "CMH_p_ctrl_ViralPhylum": cmh_phy_p,
            "CMH_OR_ctrl_FamilyXPhylum": round(cmh_fp_or, 3) if pd.notna(cmh_fp_or) else np.nan,
            "CMH_p_ctrl_FamilyXPhylum": cmh_fp_p,
            "CMH_OR_ctrl_FamilyXOriginTier": round(cmh_ft_or, 3) if pd.notna(cmh_ft_or) else np.nan,
            "CMH_p_ctrl_FamilyXOriginTier": cmh_ft_p,
            "Tier1_Tier2_n_in_hab": f"{int(np.sum(fm & hm & m_hc))}/{int(np.sum(hm & m_hc))}",
            "Tier1_Tier2_OR": round(float(hc_or), 3) if pd.notna(hc_or) and np.isfinite(hc_or) else np.nan,
            "Tier1_Tier2_p": float(hc_p) if pd.notna(hc_p) else np.nan,
            "Uroviricota_OR": round(float(uro_or), 3) if pd.notna(uro_or) and np.isfinite(uro_or) else np.nan,
            "Uroviricota_p": float(uro_p) if pd.notna(uro_p) else np.nan,
            "ReliableHost_OR": round(float(rh_or), 3) if pd.notna(rh_or) and np.isfinite(rh_or) else np.nan,
            "ReliableHost_p": float(rh_p) if pd.notna(rh_p) else np.nan,
            "CMH_OR_ctrl_HostPhylum": round(cmh_host_or, 3) if pd.notna(cmh_host_or) else np.nan,
            "CMH_p_ctrl_HostPhylum": cmh_host_p,
        })

    if out_rows:
        df_c = pd.DataFrame(out_rows)
        df_c["CMH_BH_padj_ctrl_FamilyXPhylum"] = bh_fdr(df_c["CMH_p_ctrl_FamilyXPhylum"].fillna(1.0).values)
        df_c["persists_after_FamilyXPhylum_FDR05"] = (
            (df_c["CMH_BH_padj_ctrl_FamilyXPhylum"] < 0.05)
            & (df_c["CMH_OR_ctrl_FamilyXPhylum"].fillna(0.0) > 1.0)
        )
        return df_c.to_dict("records")
    return []


def run_part3_ecology_cmh(
    v26_meta: pd.DataFrame,
    mn: pd.DataFrame,
    me: pd.DataFrame,
) -> pd.DataFrame:
    """Evaluate ecological niche enrichments before and after CMH lineage/family control."""
    print("[Part 3] Running confounder-controlled CMH ecological analyses (Full 2026 IMG/VR + Curated Motif Network)...")
    # Cohort A: Full 2026 IMG/VR curated catalog with classified Ecosystems
    df_a = v26_meta[v26_meta["ecosystem_cat"] != "Unclassified"].copy().reset_index(drop=True)
    feats_a = [
        ("PF00144_ClassC", "SBL_Family", (df_a["winning_acc"] == "PF00144").values),
        ("PF00905_Transpeptidase", "SBL_Family", (df_a["winning_acc"] == "PF00905").values),
        ("PF13354_ClassA_D", "SBL_Family", (df_a["winning_acc"] == "PF13354").values),
        ("PF00768_DD_Peptidase", "SBL_Family", (df_a["winning_acc"] == "PF00768").values),
        ("Putative_Hydrolytic_BetaLactamase", "Functional_Tier", (df_a["functional_tier"] == "Putative_Hydrolytic_BetaLactamase").values),
        ("motif1_SXXK_hmm", "Catalytic_Motif", df_a["motif1_SXXK_hmm"].astype(bool).values),
        ("motif2_SXN_PBP", "Catalytic_Motif", (df_a["motif2_class"] == "SXN_PBP").values),
        ("motif2_YAN_ClassC", "Catalytic_Motif", (df_a["motif2_class"] == "YAN_ClassC").values),
        ("motif2_SDN_ClassA", "Catalytic_Motif", (df_a["motif2_class"] == "SDN_ClassA").values),
        ("motif3_KTG_hmm", "Catalytic_Motif", df_a["motif3_KTG_hmm"].astype(bool).values),
        ("Omega_Loop_E166", "Catalytic_Motif", df_a["has_omega_E166"].astype(bool).values),
    ]
    rows_a = _evaluate_cmh_battery(
        df_a, feats_a, "ecosystem_cat", "A_Full_Curated_2026_IMGVR", min_hab_n=30, min_feat_n=25
    )

    # Cohort B: Curated 2026 Motif Network subset (motif_node.csv filtered to 2026 curated UVIGs)
    mn_v = mn[~mn["taxonomy"].isin(["M1", "M2", "M3", "BACTERIAL"])].copy()
    mn_v["uvig"] = mn_v["Id"].astype(str).str.split("__").str[0]
    df_b = mn_v.merge(v26_meta.drop_duplicates("uvig"), on="uvig", how="inner", suffixes=("_mn", ""))
    df_b = df_b[~df_b["habitat"].isin(["Unclassfied", "Unclassified", "NA"]) & df_b["habitat"].notna()].copy().reset_index(drop=True)

    # Build motif sets from motif_edge.csv
    seq_motifs: dict[str, set[str]] = defaultdict(set)
    for _, r in me.iterrows():
        seq_motifs[r["Target"]].add(r["Source"])
    motif_cats = mn[mn["taxonomy"].isin(["M1", "M2", "M3"])].set_index("Id")["taxonomy"].to_dict()
    motif_counts = Counter()
    for sid in df_b["Id"]:
        for m in seq_motifs.get(sid, ()):
            if m in motif_cats:
                motif_counts[m] += 1

    feats_b = [
        (m, motif_cats[m], np.array([m in seq_motifs.get(sid, ()) for sid in df_b["Id"]]))
        for m, cnt in motif_counts.items() if cnt >= 20
    ]
    rows_b = _evaluate_cmh_battery(
        df_b, feats_b, "habitat", "B_Curated_2026_MotifNetwork", min_hab_n=20, min_feat_n=20
    )

    df_eco = pd.DataFrame(rows_a + rows_b)
    df_eco.to_csv(OUT_DIR / "Table_Ecology_CMH_Controlled_2026.tsv", sep="\t", index=False)
    print(f"  [Part 3 Done] Cohort A (n={len(df_a)}) hits={len(rows_a)}, Cohort B (n={len(df_b)}) hits={len(rows_b)}")
    return df_eco


# ---------------------------------------------------------------------------
# Summary Report & Main Entrypoint
# ---------------------------------------------------------------------------
def write_summary_report(
    v26_meta: pd.DataFrame,
    df_loc: pd.DataFrame,
    df_div: pd.DataFrame,
    df_eco: pd.DataFrame,
) -> None:
    lines = [
        "========================================================================",
        "Stage 4B Report: 3D Biophysics, Controlled Diversity & Stratified Ecology",
        "2026 Strictly Curated CL0013 Viral Cohort (pass_cov65 & pass_fimo_1e5)",
        "========================================================================",
        f"Curated viral cohort size: n={len(v26_meta)} (matched info.xlsx: {int(v26_meta['in_info'].sum())})",
        f"Reliable host subset (Isolate/CRISPR): n={int(v26_meta['host_method'].isin(RELIABLE_HOST_METHODS).sum())}",
        f"High-confidence viral origin (Tier1+Tier2+MOSAiC): n={int(v26_meta['viral_origin_tier'].isin(HIGH_CONF_TIERS).sum())}",
        "",
        "--- 1. 3D Structural Localization of Constraint vs Diversification ---",
    ]
    for _, r in df_loc.iterrows():
        lines.append(
            f"  {r['structure_source']:18s} | {r['family']:24s} | n={int(r['n_mapped_sites']):3d} "
            f"| dist(cat/div/bg)={r['cat_median_dist_A']}/{r['div_median_dist_A']}/{r['bg_median_dist_A']} A "
            f"(p_div>cat={r['MW_p_div_gt_cat_dist']:.2e}) "
            f"| RSA(cat/div/bg)={r['cat_median_RSA']}/{r['div_median_RSA']}/{r['bg_median_RSA']} "
            f"(p_div>cat={r['MW_p_div_gt_cat_RSA']:.2e}) "
            f"| omega_ng(<10A vs >=10A&RSA>=0.20)={r['active_center_median_omega_ng']} vs {r['distal_exposed_median_omega_ng']} "
            f"(p={r['MW_p_distal_gt_active_omega_ng']:.2e})"
        )

    lines.extend(["", "--- 2. Confounder-Controlled Structural & Phylogenetic Diversity ---"])
    for _, r in df_div.iterrows():
        lines.append(
            f"  {r['stratum']:40s} | V/B PDB={int(r['n_viral_pdbs']):4d}/{int(r['n_bac_pdbs']):3d} "
            f"| rarefied_k={int(r['rarefied_k']):2d} PD(V/B)={r['viral_PD_mean']}/{r['bac_PD_mean']} (p={r['PD_MW_p']:.2e}) "
            f"| norm_Rg(V/B)={r['viral_norm_rg_median']}/{r['bac_norm_rg_median']} (MW_p={r['norm_rg_MW_p']:.2e})"
        )

    lines.extend(["", "--- 3. Confounder-Controlled Ecological Enrichment (CMH Summary) ---"])
    if not df_eco.empty:
        for coh, sub in df_eco.groupby("cohort"):
            n_g = int((sub["global_BH_padj"] < 0.05).sum())
            n_p = int(sub["persists_after_FamilyXPhylum_FDR05"].sum())
            lines.append(f"  [{coh}] evaluated={len(sub)}, global_BH<0.05={n_g}, persists_after_FamilyXPhylum_FDR05={n_p}")
            for _, r in sub.head(12).iterrows():
                lines.append(
                    f"    {r['feature']:32s} @ {r['habitat']:20s} | a={int(r['a_feat_in_hab']):4d} "
                    f"| gOR={r['global_OR']} (gFDR={r['global_BH_padj']:.2e}) "
                    f"| CMH_FamXPhy_OR={r['CMH_OR_ctrl_FamilyXPhylum']} (cFDR={r['CMH_BH_padj_ctrl_FamilyXPhylum']:.2e}) "
                    f"| persists={r['persists_after_FamilyXPhylum_FDR05']}"
                )
    (OUT_DIR / "Stage4_Structure_And_Ecology_Summary_Report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    vcat = pd.read_csv(VIRAL_CAT, sep="\t", low_memory=False)
    v26 = vcat[vcat["pass_cov65"] & vcat["pass_fimo_1e5"]].copy()
    v26["uvig"] = v26["uid"].astype(str).str.split("|").str[0]
    v26["uid_key"] = v26["uid"].astype(str).str.replace("|", "__", regex=False)
    v26 = v26.drop_duplicates("uid").reset_index(drop=True)

    info = pd.read_excel(INFO_XLSX, sheet_name="ALL_information")
    info["uvig"] = info["UVIG"].astype(str)
    info["viral_phylum"] = info["Taxonomic classification"].apply(lambda x: extract_taxon_rank(x, "p"))
    info["host_phylum"] = info["Host taxonomy prediction"].apply(lambda x: extract_taxon_rank(x, "p"))
    info["host_method"] = info["Host prediction method"].fillna("NA").astype(str)
    info["ecosystem_cat"] = info["Ecosystem classification"].apply(classify_ecosystem)
    info_sub = info[["uvig", "viral_phylum", "host_phylum", "host_method", "ecosystem_cat"]].drop_duplicates("uvig")

    v26_meta = v26.merge(info_sub, on="uvig", how="left")
    v26_meta["viral_phylum"] = v26_meta["viral_phylum"].fillna("Unassigned")
    v26_meta["host_phylum"] = v26_meta["host_phylum"].fillna("Unassigned")
    v26_meta["host_method"] = v26_meta["host_method"].fillna("NA")
    v26_meta["ecosystem_cat"] = v26_meta["ecosystem_cat"].fillna("Unclassified")
    v26_meta["in_info"] = v26_meta["uvig"].isin(set(info_sub["uvig"]))

    metrics_df = pd.read_csv(METRICS_TSV, sep="\t")
    mn = pd.read_csv(MOTIF_NODE)
    me = pd.read_csv(MOTIF_EDGE)

    _, df_loc = run_part1_3d_biophysics(v26_meta, metrics_df)
    _, df_div = run_part2_structural_diversity(v26_meta, mn, me)
    df_eco = run_part3_ecology_cmh(v26_meta, mn, me)
    write_summary_report(v26_meta, df_loc, df_div, df_eco)
    print(f"[DONE] Stage 4B outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
