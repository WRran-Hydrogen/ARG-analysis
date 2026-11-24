# Global Diversity and Functional Evolution of Viral Serine-β-Lactamase Genes (`SBL_research`)

Reproducible bioinformatics codebase and curated datasets supporting the manuscript:

> **Global analyses reveal the diversity and functional evolution of viral serine-β-lactamase genes linked to antibiotic resistance**

This repository reconstructs and extends the complete serine-β-lactamase (SBL) analytical pipeline centered on the three focal families—**`PF00144`** (Beta-lactamase, Class A/C-like), **`PF00905`** (Transpeptidase / Penicillin-binding protein domain), and **`PF13354`** (Beta-lactamase2)—while integrating the broader **Pfam Clan `CL0013` evolutionary context** (`PF00768` [`Peptidase_S11`], `PF02113` [`Peptidase_S13`], and non-penicillin-recognizing outgroup families `PF04960`, `PF06089`, `PF27310`) as auxiliary validation. The workflow implements a five-layer credibility framework covering official Pfam `--cut_ga` profiling, unified MEME/FIMO motif filtering, bacterial reference prophage decontamination, locus-level viral genomic neighborhood verification, deterministic phylogenetic horizontal gene transfer (HGT) inference, lineage-resolved and 3D-localized selection analyses, and confounder-controlled ecological modeling (full methods: [`docs/METHODS_AND_PIPELINE.md`](docs/METHODS_AND_PIPELINE.md)).

---

## Five-Layer Credibility Architecture

1. **Layer 1 — Official Pfam `CL0013` Profiling & Unified Motif Filtering (`01`, `02`)**:
   - Evaluates the 3 focal families (`PF00144.31` GA=27.0, `PF00905.29` GA=22.1, `PF13354.14` GA=22.0) alongside all other Clan `CL0013` members using official Pfam-A `--cut_ga` thresholds and rebuilt seed HMMs (`MUSCLE v5.3` + `HMMER v3.4`).
   - Relies on strictly uniform two-round `MEME/FIMO` motif filtering ($p \le 10^{-5}$ and $p \le 10^{-10}$) across viral and bacterial sets to eliminate false positives lacking conserved catalytic motifs, without enforcing a hard-coded `SXXK` regex pre-filter (`SXXK` conservation is tracked descriptively).
   - Auxiliary `CL0013` profiling demonstrates that 100.00% (`1,370/1,370`) of final viral `PF13354` SBLs belong to the penicillin-recognizing enzyme (`PRE`) subset (`88.54%` pass `PF13354` `--cut_ga`; `64.31%` simultaneously pass both `PF13354` and `PF00768:Peptidase_S11` `--cut_ga`), with `0` hits to non-penicillin-recognizing `CL0013` enzymes (`PF04960`, `PF06089`, `PF27310`).
2. **Layer 2 — Bacterial Reference Decontamination & Viral Genomic Context (`03`, `04`)**:
   - Audits all `2,220` bacterial reference sequences for mislabeled viral entries (`61` `SARG` `Viruses`), artificial/synthetic vectors (`153`), phage-annotated headers (`223`), and high-identity (`>=90%` aa identity) `NCBI nr` cryptic prophages (`393`), exporting `CleanStrict` (`n=1,961`) and `CleanCore` (`n=1,775`) subsets.
   - Stratifies all `5,701` viral SBL loci (`5,618` unique contigs) into explicit viral-origin confidence tiers (`Tier 1`: `1,254`; `Tier 2`: `3,238`; `Tier 3`: `1,209`) with prophage boundary distances and full upstream/downstream ORF Pfam neighborhoods (`n=879` mapped loci; `86.35%` carrying viral structural hallmarks on the same contig, `76.45%` flanked on both sides), following Martin et al. (*Nat. Microbiol.* 2025).
3. **Layer 3 — Deterministic Phylogenetic HGT Inference & Sensitivity Pruning (`05`)**:
   - Implements reproducible phylogenetic HGT detection modeled on Irwin et al. (*Nat. Microbiol.* 2022) and Dmitrijeva et al. (*Nat. Ecol. Evol.* 2024): empirical Youden-*J* inter-domain backbone rooting, zero-length polytomy collapse, two-pass Fitch parsimony ancestral state reconstruction, double-outgroup topological nesting (`>=80%` recipient purity, `<25%` tree tips, `>=75%` donor outgroup fraction across sister `S1` and cousin `S2`), `UFBoot` support stratification (`>=70/80/90%`), and sensitivity pruning across `CleanStrict`, `CleanCore`, and `Viral Tier 1+2`.
4. **Layer 4 — Lineage-Resolved Selection, MEME Diagnostics & 3D Structural Biophysics (`06`, `07`, `08`)**:
   - Normalizes HyPhy `aBSREL` branch selection rates against tested background branches (`Internal` vs. `Terminal`; `Viral` vs. `Bacterial` vs. `Mixed_Nested`) and tests paired root-to-tip $\omega$ transitions along HGT and non-HGT lineages.
   - Diagnoses `PF00144` MEME positive selection (`41.87%` at $p < 0.10$ vs. `29.76%` at $p < 0.05$ and `15.22%` at BH-FDR $q < 0.05$) across threshold sensitivity, branch power (`1,382` branches), sequence divergence/alignment shifts (`hmmalign` locks `SXXK` at `89.94%` vs. `3.15%` in unanchored `codon_trim.fna`), and intragenic recombination (Bruen's PHI $p > 0.61$).
   - Tests catalytic motif purifying selection via enrichment-depletion statistics (`Fisher`, `Binomial`, `10,000-permutation`), domain-anchored (`hmmalign` back-translated) Nei–Gojobori $d_N/d_S$, and residue-level 3D mapping (`Shrake–Rupley` RSA, $\phi/\psi$ secondary structure, and $C_\alpha$ distance to catalytic Serine across `110` crystal + ESMFold PDB structures).
5. **Layer 5 — Confounder-Controlled Structural Diversity & Habitat–Motif Ecology (`09`)**:
   - Controls structural diversity comparisons for SBL family (`PF00144`, `PF00905`, `PF13354`), viral phylum (`Uroviricota`), equal-sample-size rarefaction (`k=50` / `k=26`), and reliable host assignments (`Isolate + CRISPR`).
   - Applies Cochran–Mantel–Haenszel (`CMH`) stratified tests controlling for `SBL Family x Viral Phylum` and `Host Phylum` across `43` globally enriched habitat–motif pairs, distinguishing `13` confounder-robust habitat-specific motif signals from `30` taxonomically confounded associations.

---

## Repository Structure

```text
SBL_research/
├── README.md
├── run_pipeline.py              # Verification & execution entrypoint (--verify / --run-fast / --run-all)
├── environment.yml              # Conda environment (Python + HMMER/MEME/HyPhy/BLAST/IQ-TREE)
├── config/
│   ├── local.env.example        # Template: WSL/Linux tool locations for scripts 01-09 -> copy to local.env
│   └── hpc.env.example          # Template: cluster paths for scripts/hpc -> copy to hpc.env
├── docs/
│   └── METHODS_AND_PIPELINE.md  # Full methods, software versions & HPC Phase 0-6 pipeline description
├── data/                        # INPUTS (read-only for the pipeline)
│   ├── sequences/               # Original cohort: seed, initial viral (11,368), FIMO-r1, final viral (5,701), bacterial refs
│   ├── curated_2026/            # HPC re-analysis curated protein, domain-core, CDS & CL0013 FASTAs
│   ├── hmm_models/              # Official Pfam-A / Clan CL0013 HMMs (*_official.hmm)
│   ├── metadata/                # information.xlsx, Pfam-A.clans.tsv, prophage BLAST, contig ORF/Pfam annotations
│   ├── motifs/                  # Motif graph node/edge tables
│   ├── phylogeny/               # Family ML trees (PF*.nwk) & structure-guided tree (maxcc.nwk)
│   ├── selection/               # Codon MSAs, Newick trees, HyPhy aBSREL/MEME/MG94 JSONs, branch metadata
│   ├── structures/
│   │   ├── ref_pdbs/            # Reference crystal structures (1BLS, 1K25, 1BUE)
│   │   └── ProteinStructure/    # 1,508 predicted SBL structures - NOT in git, download from Figshare (see below)
│   └── processed/               # Regenerated by the pipeline
│       ├── hmm_rebuilt/         # Rebuilt seed HMMs + MUSCLE MSAs (script 01)
│       └── bacterial_clean/     # CleanStrict / CleanCore bacterial reference FASTAs (script 03)
├── scripts/
│   ├── sbl_config.py            # Single source of truth for repository paths & tool locations
│   ├── 01_build_and_verify_hmm.py ... 09_confounder_controlled_ecology.py
│   ├── 10_plot_phylogenetic_trees.R   # Circular figures of the primary HPC trees
│   ├── run_hmm_motif_pipeline.sh
│   └── hpc/                     # Slurm pipeline, numbered by stage (see "HPC pipeline" below)
│       ├── hpc_config.py        # Cluster path configuration (reads config/hpc.env or SBL_HPC_* variables)
│       ├── run_curation_and_context.sh
│       ├── 01_curation/         # 01-09: viral decontamination, MOSAiC, GTDB, marine MAGs, curated datasets
│       ├── 02_genomic_context/  # 01 neighborhoods & viral-origin tiers, 02 CheckV host-gene audit
│       ├── 03_alignment_and_trees/  # 01-03 primary trees, 04-06 CL0013 & PMSF, 07-08 structure-guided CL0013
│       ├── 04_hgt_inference/    # 01-02 rooting, ancestral states, HGT events, OAT; 03-05 AU tests
│       └── 05_selection/        # 01-04 codon alignments, saturation, HyPhy, summaries
├── results/
│   ├── tables/                  # Table_S01-Table_S13
│   ├── supplementary_data/      # Per-sequence, per-contig, HGT catalog & per-site TSVs (7 datasets)
│   └── reanalysis_2026/         # HPC re-analysis catalogs & summaries; primary_trees/ (trees, rooted trees, reports); hgt_inference/
├── figures/                     # Fig_S01-Fig_S05 (.pdf + .png); phylogenetic_trees/ (script 10)
└── legacy/                      # Original analysis code (R/Python/visualisation) from the first release
```

| Script | Reads | Writes |
|---|---|---|
| `01` | `data/sequences/`, `data/hmm_models/` | `data/processed/hmm_rebuilt/`, `Table_S01`, `PerSequence_HMM_Scores.tsv` |
| `02` | `data/sequences/` | `Table_S02`, `PerSequence_FIMO_Sensitivity.tsv` |
| `03` | `data/sequences/`, `data/metadata/` | `data/processed/bacterial_clean/`, `Table_S03`, `PerSequence_Bacterial_Audit.tsv` |
| `04` | `data/sequences/`, `data/metadata/` | `Table_S04`, `Supplementary_Data_Viral_Genomic_Neighborhoods.tsv`, `Fig_S01` |
| `05` | `data/phylogeny/`, outputs of `03`/`04` | `Table_S05`, `Table_S06`, `Supplementary_Data_HGT_Events_Catalog.tsv` |
| `06` | `data/selection/` | `Table_S07`, `Table_S08`, `Fig_S02` |
| `07` | `data/selection/`, `data/processed/hmm_rebuilt/` | `Table_S09`, `Table_S10`, `PerSite_Selection_Metrics.tsv`, `Fig_S03` |
| `08` | `data/structures/`, output of `07` | `Table_S11`, `PerSite_Selection_and_3D_Biophysics.tsv`, `Fig_S04` |
| `09` | `data/structures/`, `data/motifs/`, `data/phylogeny/maxcc.nwk`, `data/metadata/` | `Table_S12`, `Table_S13`, `Fig_S05` |
| `10` (R) | `results/reanalysis_2026/primary_trees/`, `results/reanalysis_2026/*Catalog.tsv`, Phase 2 table | `figures/phylogenetic_trees/`, `Tree_Tip_Annotations_Merged.tsv` |

---

## Core Cohort & Official Pfam `CL0013` Summary

- **Focal Viral SBL Cohort (`n = 5,701` loci across `5,618` contigs)**:
  - `PF00144` (`2,524` final from `2,540` initial): **`98.34%` (`2,482`)** pass official `PF00144.31` `--cut_ga`, `98.34%` pass unified Round-2 FIMO $p \le 10^{-5}$.
  - `PF00905` (`1,807` final from `2,356` initial): **`94.30%` (`1,704`)** pass official `PF00905.29` `--cut_ga` (enriched from `73.13%` pre-FIMO, confirming motif-filtering efficacy), `98.34%` pass Round-2 FIMO $p \le 10^{-5}$.
  - `PF13354` (`1,370` final from `6,472` initial): **`88.54%` (`1,213`)** pass official `PF13354.14` `--cut_ga`, **`100.00%` (`1,370`)** pass official Clan `CL0013` penicillin-recognizing enzyme (`PRE`) `--cut_ga` (`881` dual-matching both `PF13354` and `PF00768:Peptidase_S11`, `0` matching non-PRE `CL0013` families), and `99.20%` (`1,359`) pass Round-2 FIMO $p \le 10^{-5}$.
- **Bacterial Reference Cohort (`n = 2,220` sequences)**:
  - `PF00144` (`485`), `PF00905` (`505`), and `PF13354` (`1,230`) all achieve `100.00%` official self `--cut_ga` validation (`Table_S01`). Auditing yields **`1,961` `CleanStrict` (Tier A)** and **`1,775` `CleanCore` (Tier B)** sensitivity reference sets (`Table_S03`).

---

## Reproducibility & Quick Start

### 1. Environment Setup

```bash
conda env create -f environment.yml
conda activate sbl_virome_evolution
```

### 2. Configure Tool Locations

Scripts `01`-`09` call HMMER, MUSCLE, MEME/FIMO and a Python interpreter with `ete4` inside WSL (on Windows) or directly (on Linux/macOS). If these tools are not on `PATH`, copy the template and edit it:

```bash
cp config/local.env.example config/local.env   # SBL_WSL_PYTHON, SBL_MEME_BIN_DIR, SBL_MUSCLE, SBL_USE_WSL
```

For the cluster scripts in `scripts/hpc/`, copy `config/hpc.env.example` to `config/hpc.env` (or to `scripts/hpc/hpc.env` when copying that folder to the cluster on its own) and set `SBL_HPC_BASE` plus the upstream job folders. Both `*.env` files are git-ignored; any value can also be exported as an environment variable.

Scripts `08` and `09` additionally need the predicted structures: download **Protein Structure Data Parts 1-4** from the Figshare collection below and extract all structure files into `data/structures/ProteinStructure/`.

### 3. Verify Repository Integrity

Confirm all required inputs, official Pfam HMMs, scripts (`01`–`09`), documentation, supplementary tables (`Table_S01–S13`), supplementary TSVs (`7` datasets), and figures (`Fig_S01–S05`) are present, non-empty, and within GitHub's `<100 MB` file limit:

```bash
python run_pipeline.py --verify
```

### 4. Re-Run Analyses

```bash
# Fast re-execution of downstream statistical & structural modules (scripts 03-09)
python run_pipeline.py --run-fast

# Complete pipeline execution from official Pfam CL0013 / seed HMMs & FIMO sensitivity (scripts 01-09)
python run_pipeline.py --run-all
```

Or run individual modules sequentially:

```bash
python scripts/01_build_and_verify_hmm.py
python scripts/02_fimo_threshold_sensitivity.py
python scripts/03_bacterial_prophage_audit.py
python scripts/04_viral_genomic_context_audit.py
python scripts/05_phylogenetic_hgt_inference.py
python scripts/06_selection_absrel_and_lineages.py
python scripts/07_selection_meme_and_catalytic_motifs.py
python scripts/08_structural_localization_3d.py
python scripts/09_confounder_controlled_ecology.py
Rscript scripts/10_plot_phylogenetic_trees.R   # needs ape, phangorn, ggtree, ggtreeExtra, ggnewscale, cowplot, ragg
```

> **Note on large HPC intermediates**: Multi-gigabyte raw HPC intermediates (`bac120_metadata.tsv.gz`, full `Pfam-A.domtbl`, raw bacterial genomes, geNomad/CheckV/eggNOG job outputs) remain on the HPC cluster; the curated sequence sets (`data/curated_2026/`), audit catalogs (`results/reanalysis_2026/`) and the scripts that produced them (`scripts/hpc/`) are included here.

### 5. HPC Pipeline (`scripts/hpc/`)

Each stage folder matches an output folder under `$SBL_HPC_BASE`; scripts are numbered in run order. Python
scripts find `hpc_config.py` one level up; shell scripts read `config/hpc.env` or `scripts/hpc/hpc.env`.
Submit each driver from its own folder (Slurm uses the submit directory to locate the scripts).

| Stage | Order | Entry point |
|---|---|---|
| `01_curation/` | `01` viral decontamination; `02`-`04` MOSAiC; `05`-`07` GTDB genomes, SBL scan, contig shards for geNomad; `08` marine MAGs | `sbatch 01_viral_decontamination.sh`, then the Python scripts in order |
| curated datasets + context | `01_curation/09 --stage curate` -> `02_genomic_context/01` -> `01_curation/09 --stage finalize` | `sbatch run_curation_and_context.sh` (from `scripts/hpc/`) |
| `02_genomic_context/` | `02` CheckV host-gene audit (read-only, writes `Phase2_Host_Gene_Check_Audit*`) | `python 02_audit_host_gene_check.py` |
| `03_alignment_and_trees/` | `01`-`03` alignments + wave-1/wave-2 IQ-TREE; `04`-`06` CL0013 superfamily + PMSF C20; `07`-`08` structure-guided CL0013 + clean GTDB outgroups | `bash 03_submit_primary_trees.sh`, `bash 06_submit_cl0013_and_pmsf.sh`, `08_run_cl0013_c20.sh` |
| `04_hgt_inference/` | `01` rooting, ancestral states, HGT events, OAT robustness, AU-test constraints; `03`-`04` AU tests of Tier-1 events | `sbatch 02_run_hgt_inference.sh`, then `bash 05_submit_au_tests.sh` |
| `05_selection/` | `01` codon alignments & saturation; `02` HyPhy; `03` summaries | `bash 04_submit_selection.sh` |

---

## Data Availability

Large supporting datasets are archived in the **SBL-virus Supporting Data** Figshare collection - DOI [10.6084/m9.figshare.c.8403348](https://doi.org/10.6084/m9.figshare.c.8403348):

| Dataset | DOI |
|---|---|
| Global Distribution Data | [10.6084/m9.figshare.31942896](https://doi.org/10.6084/m9.figshare.31942896) |
| Motif Data | [10.6084/m9.figshare.31942899](https://doi.org/10.6084/m9.figshare.31942899) |
| Phylogeny Data | [10.6084/m9.figshare.31942905](https://doi.org/10.6084/m9.figshare.31942905) |
| Prophage Data | [10.6084/m9.figshare.31942911](https://doi.org/10.6084/m9.figshare.31942911) |
| DomainDistance Data | [10.6084/m9.figshare.31942914](https://doi.org/10.6084/m9.figshare.31942914) |
| OriginalFile Data | [10.6084/m9.figshare.31942917](https://doi.org/10.6084/m9.figshare.31942917) |
| SelectionPressure Data | [10.6084/m9.figshare.31942920](https://doi.org/10.6084/m9.figshare.31942920) |
| MolecularDocking Data | [10.6084/m9.figshare.31946697](https://doi.org/10.6084/m9.figshare.31946697) |
| Gephi Network Data | [10.6084/m9.figshare.31942926](https://doi.org/10.6084/m9.figshare.31942926) |
| Protein Structure Data (Part 1) | [10.6084/m9.figshare.31942962](https://doi.org/10.6084/m9.figshare.31942962) |
| Protein Structure Data (Part 2) | [10.6084/m9.figshare.31943334](https://doi.org/10.6084/m9.figshare.31943334) |
| Protein Structure Data (Part 3) | [10.6084/m9.figshare.31943400](https://doi.org/10.6084/m9.figshare.31943400) |
| Protein Structure Data (Part 4) | [10.6084/m9.figshare.31943583](https://doi.org/10.6084/m9.figshare.31943583) |

## Legacy Code

`legacy/` preserves the analysis and visualisation code of the first release (R plotting, motif enrichment, selection-pressure plots, molecular-docking helpers). It is kept for provenance and is not called by the pipeline above. The earlier Observable/JavaScript visualisations (bubble world map, Sankey graph) were removed; all figures are now made with Python and R; they remain available in the git history (branch `master`, commit `925e189`).

## License

MIT License
