# Global Diversity and Functional Evolution of Viral Serine-β-Lactamase Genes (`SBL_research`)

Reproducible bioinformatics codebase and curated datasets supporting the manuscript:

> **Global analyses reveal the diversity and functional evolution of viral serine-β-lactamase genes linked to antibiotic resistance**

This repository reconstructs and extends the complete serine-β-lactamase (SBL) analytical pipeline centered on the three focal families—**`PF00144`** (Beta-lactamase, Class A/C-like), **`PF00905`** (Transpeptidase / Penicillin-binding protein domain), and **`PF13354`** (Beta-lactamase2)—while integrating the broader **Pfam Clan `CL0013` evolutionary context** (`PF00768` [`Peptidase_S11`], `PF02113` [`Peptidase_S13`], and non-penicillin-recognizing outgroup families `PF04960`, `PF06089`, `PF27310`) as auxiliary validation. The workflow implements a five-layer credibility framework covering official Pfam `--cut_ga` profiling, unified MEME/FIMO motif filtering, bacterial reference prophage decontamination, locus-level viral genomic neighborhood verification, deterministic phylogenetic horizontal gene transfer (HGT) inference, lineage-resolved and 3D-localized selection analyses, and confounder-controlled ecological modeling (full methods: [`docs/METHODS_AND_PIPELINE.md`](docs/METHODS_AND_PIPELINE.md)).

---

## Five-Layer Credibility Architecture

1. **Layer 1 — Official Pfam `CL0013` Profiling, Unified Motif Filtering & Strict Biochemical Stratification (`01_curation_and_context`)**:
   - Evaluates the three focal families (`PF00144.31`, `PF00905.29`, `PF13354.14`) alongside sister PRE families (`PF00768`, `PF02113`) and non-PRE `CL0013` outgroups (`PF04960`, `PF06089`, `PF27310`) using official Pfam-A `--cut_ga` competition, two-round `MEME/FIMO 5.5.9` ($p \le 10^{-5}$), `hmm_cov >= 0.65`, and independent domain-slice `MMseqs2` `nr85` clustering across `n = 4,203` strictly curated viral `CL0013` loci (`PF00144=1,876`, `PF00905=1,242`, `PF13354=248`, `PF00768=837`; viral `nr85`: `416`, `215`, `63`, `195`) and a triple-source bacterial reference (`GTDB r220 + CARD/SARG + Marine MAGs`).
   - Applies strict CARD + `hmmalign` A2M catalytic-triad verification (`Viral_Enzymatic_Function_Verification_2026.tsv`) excluding PBP-resistance entries, `BlaR1`/`MecR1` (`PF05569`) sensor receptors, `PF00144` `H[TS]G` R61 DD-peptidases/esterases, and truncated fragments without an intact triad, plus Class A Ω-loop (`E166-x-x-[LIVM]-N`) and Class D `[YF]GN` + `TrpCap` verification: **`155 / 4,203` (`3.69%`)** viral sequences are intact **`Putative_Hydrolytic_BetaLactamase`** (`PF13354=145` [`144` Verified + `1` Candidate], `PF00905=9` Class D `OXA`, `PF00144=1` `YRC-1`, `PF00768=0`), whereas **`4,048 / 4,203` (`96.31%`)** are **`PBP_or_PRE_Homolog`** (`4,042` PBP/esterase/phage-lysin homologs, `3` `BlaR1_Sensor_Receptor`, and `3` `Truncated_SBL_Pseudogene`).
2. **Layer 2 — Bacterial Reference Decontamination & Viral Genomic Context (`01_curation_and_context`)**:
   - Filters whole viral contigs embedded in bacterial assemblies/MAGs (`is_viral_contig_in_genome`), tags chromosomal prophage loci for OAT sensitivity exclusion, and stratifies `4,193` IMG/VR + `10` MOSAiC viral loci across `298,926` Prodigal ORFs (`Phase2_Viral_Genomic_Neighborhoods_Full.tsv`), strictly separating recombination/integrase markers (`44.96%`) from structural virion hallmarks (`Tier 1` high-confidence viral origin: `1,649` IMG/VR + `10` MOSAiC).
3. **Layer 3 — Phylogenetic HGT Inference with Rooting Checks & OAT Robustness (`03_primary_trees`, `04_hgt_inference`)**:
   - Infers HGT on `5` primary trees (four single families + joint `PF13354+PF00768`, the latter a sensitivity analysis whose events repeating a single-family event are not added to totals) with concordant Fitch parsimony + F81 marginal ML (`PastML`-style) ancestral states, double-outgroup nesting and parent/grandparent support thresholds, and tracks every event across `85` OAT trees (`MAD`, `MinVar`, `Midpoint`, clean `12`-taxon chromosomal GTDB outgroups, legacy 24-taxon outgroups, `LG+C20` PMSF, MAFFT alignments, `rep2..rep5`/`phylo_rep` subsampling, `without_prophage`, `cov60`/`cov70`).
   - **Rooting**: crown projection from the structure-guided `CL0013.struct_c20.no_PF27310` (`FoldMason + hmmalign`, `LG+C20+F+G`) superfamily tree is used only for a family that is a supported clade with a supported root position (UFBoot `>= 95`) in every independent CL0013 taxon sample; no family met this (root-split UFBoot `34–97`, never `>= 95` in both samples), so the primary root is **MAD** for every tree and all comparators are reported (`Rooting_Decisions_And_Agreement.tsv`). Every tree is exported as `{fam}.primary.rooted.treefile`, so R figures and Python statistics share one root.
   - **Results** (viral-evidence gate: `>= 1` Tier 1/Tier 2/MOSAiC viral locus on the viral side of each event): single-family trees **Tier 1 = `5`**, all single-tip and all `Root_Invariant` — `PF00144` `3` B→V + `1` V→B (`PBP_or_PRE_Homolog`), `PF13354` `1` B→V hydrolytic SBL (BlaZ-like); `PF00905` and `PF00768` `0`; Tier 2 = `46`; `Prophage_Borne = 0`. Joint tree: Tier 1 = `3` hydrolytic SBL events (BlaZ-like and AER-like repeat single-family events, CfxA-like V→B is joint-only). Pre-registered Tier-1 events are tested with AU tests, and single-tip recipients are checked for long terminal branches. These criteria detect recent transfers only; deep transfers are listed in `All_Trees_HGT_Events_Long.tsv` but cannot reach Tier 1 (see `docs/METHODS_AND_PIPELINE.md`).
   - **Ancient HGT (exploratory, not pre-registered)**: deep viral clades on the primary trees and on reduced LG+C60 PMSF backbone trees (UFBoot/SH-aLRT and TBE support), with the GTDB donor lineage of each clade's bacterial sister compared with its IMG/VR host lineage (permutation tests), a root-free host-proximity test and AU tests of deep-clade placement (`docs/METHODS_AND_PIPELINE.md`, Section 4.3).
4. **Layer 4 — Codon Saturation Diagnostics, Catalytic Purifying Selection & 3D Structural Biophysics (`05_selection_pressure`, `06_structure_and_ecology`)**:
   - Documents family-wide synonymous saturation (`dS_Saturation_And_Xia_Test_2026.tsv`; `median dS = 3.0`, `67.6–78.8%` saturated pairs, `pos3_unsaturated_sym = False` on full-family MSAs) while evaluating `21` completed core HyPhy `GARD`/`aBSREL`/`MEME`/`FEL` datasets (`5` `N=90` benchmarks + `16` subclades, including unsaturated viral/HGT subclades such as `PF00144_sub4_Viral_clade` `median dS = 0.7239` and `PF00768_sub1_HGT_clade` `median dS = 0.7608`, `pos3_unsaturated_sym = True`) plus `PF00144_N240` (`hmm`/`mafft`) `MEME`/`FEL` site benchmarks.
   - Proves catalytic motif purifying constraint across `1,124` HMM columns (`Catalytic_Motif_Purifying_Proof_2026.tsv`: motif core median $\omega_{\text{NG}} = 0.3997$ vs. background $0.8282$, `MW p = 3.20e-14`; `42.22%` vs. `3.61%` constrained columns, `Fisher OR = 19.49`, `p = 1.45e-14`) and 3D proximity in pooled crystal (`n = 809`, `omega_ng <10 Å` vs. surface `0.6088 vs. 0.9162`, `p = 1.76e-15`) and viral ESMFold structures (`n = 1,064`, `0.5929 vs. 0.8930`, `p = 4.74e-20`).
5. **Layer 5 — Confounder-Controlled Structural Diversity & Habitat–Motif Ecology (`06_structure_and_ecology`)**:
   - Evaluates `825` curated 2026 viral and `170` bacterial PDB structures (`PDB_3D_Structural_Catalog_2026.tsv` and `Table_Structural_Diversity_Controlled_2026.tsv`). Crucially, once stratified by verified enzymatic tier, true hydrolytic viral and bacterial SBLs exhibit virtually identical normalized radius of gyration (`4_Tier_Putative_Hydrolytic_BetaLactamase`, V/B = `19/137`, median `norm_Rg` `2.6885` vs. `2.6862`, `MW p = 0.192`), and PBP/PRE homologs also match bacterial PBP controls (`4_Tier_PBP_or_PRE_Homolog`, V/B = `806/33`, median `norm_Rg` `2.7638` vs. `2.8009`, `MW p = 0.512`).
   - Applies Cochran–Mantel–Haenszel (`CMH`) tests controlling for `SBL Family × Viral Phylum`, finding `13 / 38` Cohort A (`n = 3,956`) and `11 / 16` Cohort B (`n = 813`) globally enriched habitat–motif associations remain significant after taxonomic confounding adjustment.

---

## Repository Structure

```text
SBL_research/
├── README.md
├── environment.yml                         # Conda environment (local tools, R figure packages)
├── envs/hpc.yml                            # HPC environment specification
├── config/
│   ├── local.env.example                   # Template -> copy to config/local.env (git-ignored)
│   └── hpc.env.example                     # Template -> copy to config/hpc.env (git-ignored)
├── docs/
│   └── METHODS_AND_PIPELINE.md             # Full methods, software versions & pipeline specification
├── data/                                   # INPUTS (read-only for the pipeline)
│   ├── curated_2026/                       # Curated protein, domain-slice, CDS & CL0013 FASTAs
│   ├── sequences/                          # Official Pfam seed alignments (*_seed.faa)
│   ├── hmm_models/                         # Official Pfam-A / Clan CL0013 HMMs (*_official.hmm)
│   ├── metadata/                           # information.xlsx, Pfam-A.clans.tsv, Color_Mapping_Guide.txt
│   ├── motifs/                             # Motif graph node/edge tables
│   ├── phylogeny/                          # Structure-guided tree (maxcc.nwk)
│   ├── structures/
│   │   ├── ref_pdbs/                       # Reference crystal structures + cl0013_crystal_anchors/ (14 FoldMason seeds)
│   │   └── ProteinStructure/               # 1,508 predicted SBL structures (Figshare, git-ignored)
│   └── legacy_2023/                        # Archived 2023 baseline inputs (selection, processed, sequences, phylogeny, metadata)
├── scripts/
│   ├── sbl_config.py                       # Shared local/WSL path & tool configuration (with legacy_2023 fallback)
│   ├── hpc/                                # Slurm pipeline, one folder per stage (see "HPC pipeline" below)
│   │   ├── hpc_config.py                   # Cluster paths (config/hpc.env or SBL_HPC_* variables)
│   │   ├── run_curation_and_context.sh
│   │   ├── 01_curation/                    # 01-09: viral decontamination, MOSAiC, GTDB, marine MAGs, curated datasets
│   │   ├── 02_genomic_context/             # 01 neighborhoods & viral-origin tiers, 02 CheckV host-gene audit, 03 gene-origin evidence audit
│   │   ├── 03_alignment_and_trees/         # 01-03 primary trees, 04-06 CL0013 & PMSF, 07-09 structure-guided CL0013 & joint tree
│   │   ├── 04_hgt_inference/               # 01-02 rooting, ancestral states, HGT events, OAT; 03-05 AU tests; 06-11 ancient HGT
│   │   └── 05_selection/                   # 01-04 codon alignments, saturation, HyPhy, summaries
│   ├── local_analysis/                     # Enzymatic re-stratification & structure/ecology analyses (run locally)
│   └── visualization/                      # R figure scripts (trees, CL0013 superfamily, HGT robustness)
├── results/
│   └── reanalysis_2026/                    # Stage-numbered results
│       ├── 01_curation_and_context/        # Curated catalogs, Viral_Enzymatic_Function_Verification_2026.tsv, Phase 2 tables & upstream_audits/
│       ├── 03_primary_trees/               # Primary trees, rooted trees, IQ-TREE reports, CL0013 superfamily trees & tip annotations
│       ├── 04_hgt_inference/               # HGT event catalogs, rooting decisions, OAT matrix & summary report
│       ├── 05_selection_pressure/          # Codon saturation/Xia tests, motif constraint & HyPhy selection
│       └── 06_structure_and_ecology/       # 3D crystal/ESMFold biophysics, diversity & CMH ecology tables
├── figures/
│   └── trees_2026_R/                       # Phylogenetic, domain-context, CL0013 and HGT robustness figures
└── archive/
    ├── v0_original_submission/             # First-release R/Python/docking code (provenance only)
    └── v1_baseline_5701_audit/             # Baseline 5,701-cohort audit (run_pipeline.py, scripts 01-09, Table_S01-S13, Fig_S01-S05)
```

---

## Reproducibility & Quick Start

### 1. Environment & Tool Setup

```bash
conda env create -f environment.yml
conda activate sbl_virome_evolution
cp config/local.env.example config/local.env   # Set SBL_WSL_PYTHON on Windows/WSL
cp config/hpc.env.example config/hpc.env       # Set SBL_HPC_BASE for Slurm runs
```

### 2. HPC Pipeline (`scripts/hpc/`)

Each stage folder matches an output folder under `$SBL_HPC_BASE`; scripts are numbered in run order. Python
scripts find `hpc_config.py` one level up; shell scripts read `config/hpc.env` or `scripts/hpc/hpc.env`.
Submit each driver from its own folder (Slurm uses the submit directory to locate the scripts).

| Stage | Order | Entry point |
|---|---|---|
| `01_curation/` | `01` viral decontamination; `02`-`04` MOSAiC; `05`-`07` GTDB genomes, SBL scan, contig shards for geNomad; `08` marine MAGs | `sbatch 01_viral_decontamination.sh`, then the Python scripts in order |
| curated datasets + context | `01_curation/09 --stage curate` -> `02_genomic_context/01` -> `01_curation/09 --stage finalize` | `sbatch run_curation_and_context.sh` (from `scripts/hpc/`) |
| `02_genomic_context/` | `02` CheckV host-gene audit (`Phase2_Host_Gene_Check_Audit*`); `03` gene-origin evidence audit (`origin_evidence_audit/`); both read-only | `python 02_audit_host_gene_check.py`, `python 03_audit_gene_origin_evidence.py` |
| `03_alignment_and_trees/` | `01`-`03` alignments + wave-1/wave-2 IQ-TREE; `04`-`06` CL0013 superfamily + PMSF C20; `07`-`08` structure-guided CL0013 + clean GTDB outgroups; `09` structure-guided joint tree | `bash 03_submit_primary_trees.sh`, `bash 06_submit_cl0013_and_pmsf.sh`, `08_run_cl0013_c20.sh`, `bash 09_submit_joint_struct_tree.sh` |
| `04_hgt_inference/` | `01` rooting, ancestral states, HGT events, OAT robustness, AU-test constraints; `03`-`04` AU tests of Tier-1 events | `sbatch 02_run_hgt_inference.sh`, then `bash 05_submit_au_tests.sh` |
| `04_hgt_inference/` (ancient HGT) | `06` IMG/VR host and GTDB lineages; `07` backbone alignments; `08` C60 backbone trees (UFBoot, TBE); `10` deep clades, donor lineage, host concordance, AU constraints; `11` AU tests | `export SBL_HPC_IMGVR_HOST_TABLE=<IMG/VR host table>`, then `bash 09_submit_ancient_hgt.sh` |
| `05_selection/` | `01` codon alignments & saturation; `02` HyPhy; `03` summaries | `bash 04_submit_selection.sh` |

### 3. Local Analyses & Figures

- **Enzymatic function verification** (local catalogs + `Viral_Enzymatic_Function_Verification_2026.tsv`; the HPC catalogs already carry these tiers):
  ```bash
  python scripts/local_analysis/reclassify_enzymatic_function_2026.py
  ```
- **Structure, diversity & CMH ecology**: download **Protein Structure Data Parts 1–4** from the Figshare collection below into `data/structures/ProteinStructure/`, then run:
  ```bash
  python scripts/local_analysis/phase6_structure_and_ecology_2026.py
  ```
- **Tree & HGT figures (R)**:
  ```bash
  Rscript scripts/visualization/plot_phylogenetic_trees_2026.R
  Rscript scripts/visualization/plot_phylogeny_cl0013_and_hgt_2026.R
  ```

> **Note on large HPC files**: Multi-gigabyte raw HPC intermediates (`bac120_metadata.tsv.gz`, full `Pfam-A.domtbl`, raw bacterial genomes, geNomad/CheckV/eggNOG job outputs) remain on the cluster; earlier first-release and 5,701-baseline code/tables are preserved in `archive/`, and superseded 2023 input datasets in `data/legacy_2023/`.

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

## Legacy & Archived Files

- **`archive/v0_original_submission/`**: Preserves the analysis and visualization code released with the first submission (R plotting, motif enrichment, selection-pressure plots, molecular-docking helpers). The first submission's Observable/JavaScript visualizations remain available in git history (branch `master`, commit `925e189`).
- **`archive/v1_baseline_5701_audit/`**: Preserves the baseline 5,701-cohort audit (`run_pipeline.py`, scripts `01–09`, `Table_S01–S13`, `Fig_S01–S05`).
- **`data/legacy_2023/`**: Preserves superseded 2023 baseline input datasets (`selection/`, `processed/`, legacy FASTAs, trees, and BLAST/Pfam metadata), isolated from active 2026 pipeline inputs while `scripts/sbl_config.py` retains a fallback resolver for provenance reproducibility.

## License

MIT License
