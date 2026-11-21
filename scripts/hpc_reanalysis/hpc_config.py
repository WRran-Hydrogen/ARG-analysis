"""Path and tool configuration for the 2026 HPC re-analysis scripts (Slurm cluster).

Values are resolved in this order:
1. environment variables,
2. ``hpc.env`` next to this file (handy when the folder is copied to the cluster alone),
3. ``config/hpc.env`` in the repository (git-ignored; copy ``config/hpc.env.example``),
4. the defaults below.

Only ``SBL_HPC_BASE`` normally needs to be set; everything else defaults to a
location under it. Upstream job folders (geNomad, Prodigal, Pfam scans, eggNOG)
can be pointed anywhere with their own variables.
"""

from __future__ import annotations

import os
from pathlib import Path

_HERE = Path(__file__).resolve().parent


def _load_env_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("export "):
            line = line[len("export "):]
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        values[key.strip()] = os.path.expandvars(val.strip().strip('"').strip("'"))
    return values


_FILE_VALUES: dict[str, str] = {}
for _candidate in (_HERE.parents[1] / "config" / "hpc.env", _HERE / "hpc.env"):
    _FILE_VALUES.update(_load_env_file(_candidate))


def setting(key: str, default: str | Path = "") -> str:
    return os.environ.get(key) or _FILE_VALUES.get(key) or str(default)


def path_setting(key: str, default: str | Path) -> Path:
    return Path(setting(key, default)).expanduser()


def _tool(bin_dir: str, name: str) -> str:
    return f"{bin_dir.rstrip('/')}/{name}" if bin_dir else name


# ---------------------------------------------------------------------------
# Roots
# ---------------------------------------------------------------------------
BASE = path_setting("SBL_HPC_BASE", Path.cwd() / "sbl_reanalysis")   # pipeline working root
DB_DIR = path_setting("SBL_HPC_DB_DIR", BASE / "sbl_db")               # Pfam-A.hmm, Pfam-A.clans.tsv
RAW_INPUT_DIR = path_setting("SBL_HPC_RAW_INPUT_DIR", BASE / "raw_inputs")      # MOSAiC FASTA/TSV, vir.info, euk_info.tsv
MARINE_MAGS_DIR = path_setting("SBL_HPC_MARINE_MAGS_DIR", BASE / "raw_inputs" / "marine_mags")
JOBS_DIR = path_setting("SBL_HPC_JOBS_DIR", BASE / "upstream_jobs")

# Upstream job outputs (one folder per completed upstream job)
GENOMAD_DIR = path_setting("SBL_HPC_GENOMAD_DIR", JOBS_DIR / "genomad")
PRODIGAL_DIR = path_setting("SBL_HPC_PRODIGAL_DIR", JOBS_DIR / "prodigal")
VIRAL_PFAM_DIR = path_setting("SBL_HPC_VIRAL_PFAM_DIR", JOBS_DIR / "viral_sbl_pfam")   # sbl_proteins.faa + domtbl
STEP1_DIR = path_setting("SBL_HPC_STEP1_DIR", JOBS_DIR / "step1")
ORFS_PFAM_DIR = path_setting("SBL_HPC_ORFS_PFAM_DIR", JOBS_DIR / "orfs_pfam")
EGGNOG_DIR = path_setting("SBL_HPC_EGGNOG_DIR", JOBS_DIR / "eggnog")

PFAM_ALL_HMM = DB_DIR / "Pfam-A.hmm"
PFAM_CLANS_TSV = DB_DIR / "Pfam-A.clans.tsv"

# ---------------------------------------------------------------------------
# Tools (empty bin dir -> found on PATH)
# ---------------------------------------------------------------------------
MEME_BIN_DIR = setting("SBL_HPC_MEME_BIN_DIR")
MEME_BIN = _tool(MEME_BIN_DIR, "meme")
FIMO_BIN = _tool(MEME_BIN_DIR, "fimo")
HMMSEARCH_BIN = setting("SBL_HPC_HMMSEARCH", "hmmsearch")
HMMALIGN_BIN = setting("SBL_HPC_HMMALIGN", str(Path(HMMSEARCH_BIN).with_name("hmmalign")) if "/" in HMMSEARCH_BIN else "hmmalign")
MMSEQS_BIN = setting("SBL_HPC_MMSEQS", "mmseqs")

# Phase 3 (alignment & trees)
MAFFT_BIN = setting("SBL_HPC_MAFFT", "mafft")
CLIPKIT_BIN = setting("SBL_HPC_CLIPKIT", "clipkit")
IQTREE_BIN = setting("SBL_HPC_IQTREE", "iqtree2")
PHASE3_DIR = path_setting("SBL_HPC_PHASE3_DIR", BASE / "03_alignment_and_trees")
