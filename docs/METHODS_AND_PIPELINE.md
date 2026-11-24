# Viral Serine-β-Lactamase / PRE (`CL0013`) Re-Analysis: Methods, Checkpoints & Session Reference

> **Purpose**: Authoritative, reproducible Materials & Methods specification **and** complete HPC/local execution handoff (`Memorize`) for the re-analysis of:
> *"Global analyses reveal the diversity and functional evolution of viral serine-β-lactamase / penicillin-recognizing enzymes (`CL0013`)"*
> Last updated: `2026-11-24` (**Stage 1 Completed**: Strict CL0013 bitscore taxonomy, own-seed FIMO, `hmm_cov >= 0.65`, independent domain-slice `nr85`, CARD + catalytic functional stratification, tagged bacterial prophages, strict structural hallmark Phase 2 `Tier 1`, and locked Primary HGT protocol).
>
> ✅ **Stage 1 re-run with the fixes below on 2026-11-21** (Slurm jobs `5120440`, `5120495`; `5120556` adds merged multi-piece HMM coverage; **final run `5120575`** restricts merge candidates, counts unchanged). Sections 2.3–3 show the counts of job `5120575`, which are the Phase 3 inputs. Fixes applied: (i) CARD reference set now excludes PBP-type entries (mecA/mecC, PBP1a/2b/2x/3); (ii) Class A Ω-loop requires the `E166-x-x-L-N` signature (previously any glutamate in the window) with A2M match-state verification; (iii) whole viral contigs inside bacterial assemblies/MAGs are tagged `viral_contig` and excluded from bacterial sets (no longer treated as prophages); (iv) PF13354/PF00768 trees are rooted only with PF00144 + PF00905; (v) `cov60` sensitivity applied to bacteria too; (vi) strict Phase 2 `Tier 1` (see Section 3); (vii) HMM coverage merges co-linear same-family hits (see Section 2.3). Reproduce with `scripts/hpc/run_curation_and_context.sh` (`--stage curate` → Phase 2 → `--stage finalize`).

---

## 1. Environment, Compute Infrastructure & Repository Layout

### 1.1 Compute Infrastructure & Verified Software Versions
- **HPC Cluster**: institutional Slurm cluster (`compute` partition; outbound HTTPS downloads from login node via `aria2c 1.37.0` + `curl`).
- **Dual Environment Specification**:
  - **Local baseline audit environment** ([`environment.yml`](../environment.yml)): `python=3.11`, `meme=5.0.5`, `hmmer=3.4`, `muscle=5.3`, `hyphy=2.5.93`, `iqtree>=2.2`, `clipkit>=2.10`.
  - **HPC re-analysis environment** ([`envs/hpc.yml`](../envs/hpc.yml)):
    - **HMMER `3.4`**, **MMseqs2 `18.8cc5c`**, **MEME/FIMO `5.5.9`**, **DIAMOND `2.1.24`**
    - **CheckV `v1.0.3`** (`checkv-db-v1.5`, 13 GB)
    - **geNomad `1.7`** (`genomad_db` v1.7 — IMG/VR viral screening) **and `geNomad 1.8.1`** (bacterial 8-shard screening & Marine MAGs)
    - **eggNOG-mapper `2.1.12`**, **Prodigal `2.6.3` (`-p meta`)**
    - **MAFFT `7.526`**, **MUSCLE `5.3`**, **ClipKIT `2.14.0`**, **trimAl `1.5.0`**, **PAL2NAL `14.1`**
    - **IQ-TREE `2.4.0`** (with fallback `/Phylo/bin/iqtree2` `2.3.6`), **PastML `1.9.42`** / `ete4`, **HyPhy `2.5.93(MP)`**, **Foldseek `10.941cd33`**, **FoldMason** (build `4bd8b8d`), **Python `3.13`**

### 1.2 Self-Contained Repository & HPC Directory Map
- **GitHub Repository Root**: configures local paths via [`config/local.env.example`](../config/local.env.example) / [`scripts/sbl_config.py`](../scripts/sbl_config.py) and HPC paths via [`config/hpc.env.example`](../config/hpc.env.example) / [`scripts/hpc/hpc_config.py`](../scripts/hpc/hpc_config.py).
  - `data/curated_2026/` — Stage 1 curated 2026 protein, domain-slice (`env ± 10 aa`), 1:1 verified CDS (`*.fna`), viral & bacterial `nr85` clusters, `rep1..rep5` + `phylo_rep` combined tree inputs, `cov60`/`cov70` sensitivity inputs, Sister-PRE outgroups, and `CL0013` non-PRE biological outgroups (`PF04960`, `PF06089`, `PF27310`).
  - `data/hmm_models/` & `data/sequences/` — Official Pfam HMMs (`CL0013_official.hmm` + all 8 single-family HMMs) and official seed alignments (including `PF00768_seed.faa` & `PF02113_seed.faa` fetched from InterPro REST API).
  - `results/reanalysis_2026/` — Authoritative Stage 1/2 summary reports and per-sequence / per-locus TSV catalogs (`Final_Viral_SBL_Catalog.tsv`, `Final_Bacterial_TripleSource_Catalog.tsv`, `Phase2_Viral_Genomic_Neighborhoods_Full.tsv`).
  - `results/tables/`, `results/supplementary_data/`, `figures/` — Currently hold baseline `5,701` outputs (`Table_S01–S13`, `Fig_S01–S05`), to be migrated to the **2026 curated `nr85` primary cohort** in Stage 4 upon completion of Stages 2–4.
  - `scripts/hpc/` — Complete HPC workflow, one numbered folder per stage, matching the HPC output folders: `01_curation/` (`01` viral decontamination … `09` curated datasets), `02_genomic_context/` (`01` neighborhoods & tiers, `02` CheckV host-gene audit), `03_alignment_and_trees/` (`01`–`08`), `04_hgt_inference/` (`01`–`02`), `05_selection/` (`01`–`04`); `run_curation_and_context.sh` runs curate → Phase 2 → finalize; `hpc_config.py` holds the paths.
  - `scripts/10_plot_phylogenetic_trees.R` — circular figures of the five primary trees (`figures/phylogenetic_trees/`), rooted on the Sister-PRE outgroup.

---

## 2. Stage 1 (Completed): Phase 1/1B Decontamination, Strict Clan Taxonomy, Functional Stratification & Domain-Slice `nr85`

### 2.1 Viral Sources, Decontamination & 1:1 CDS Verification
1. **IMG/VR v4 Baseline (`n = 5,701` across `5,618` contigs)**:
   - 1:1 translation-verified CDS mapping (`99.98%`, `5,700 / 5,701`) against 16-part Prodigal `-p meta` ORFs + `ARG_cds.fna`.
   - **Scheme B Joint Decontamination (`geNomad v1.7` + `CheckV v1.0.3` / DB v1.5)**: Excluded `8` geNomad plasmids, `620` host-flank provirus contaminants (`126` geNomad + `494` CheckV), `163` chromosome/plasmid-dominant non-viral contigs, and `215` unverified non-summary fragments (rescuing only `4` fragments with `checkv_viral_genes >= 1 & host_genes == 0`), leaving **`4,695`** surviving IMG/VR proteins (`drop_no_pfam_ga = 82`, `drop_non_CL0013_only = 46`).
2. **Arctic MOSAiC Viral Addition (`clean_mosaic_euk.py`)**:
   - Re-scanned full `mosaic_virus.faa` (650 MB) with `CL0013_official.hmm --cut_ga` + full `Pfam-A.hmm --cut_ga` competition (`83` -> `31` PRE hits, excluding `49` `PF04960` + `2` `PF06089`), then purged `6` contigs overlapping eukaryotic MAGs in `euk_info.tsv` -> **`25` clean MOSAiC viral PREs** (`22 PF00144`, `1 PF00905`, `2 PF00768`; protein phylogenies only, `contig_context_available = False`).

### 2.2 Verification of All 8 Clan `CL0013` Members (Including `PF27310`)
InterPro REST Pfam Stockholm headers directly verify all 8 official Clan `CL0013` (`#=GF CL CL0013`) entries:
- **Penicillin-Recognizing Enzymes (PRE)**: `PF00144.31` (`Beta-lactamase`, LENG=329, GA=27.0), `PF00905.29` (`Transpeptidase`, LENG=300, GA=22.1), `PF13354.14` (`Beta-lactamase2`, LENG=256, GA=22.0), `PF00768.26` (`Peptidase_S11`, LENG=239, GA=27.6), `PF02113.22` (`Peptidase_S13`, LENG=445, GA=22.8).
- **Non-PRE Outgroup Enzymes**: `PF04960.22` (`Glutaminase`, LENG=286, GA=27.0), `PF06089.19` (`Asparaginase_II`, LENG=322, GA=25.0), `PF27310.2` (`HrcA_C_epsilon`, `#=GF CL CL0013`, LENG=188, GA=27.0).

### 2.3 Strict Clan Bitscore Assignment, Own-Seed FIMO, Coverage Filter & Tagged Prophages
1. **Strict Clan Highest-Bitscore Assignment (`resolve_strict_cl0013`)**:
   - Within overlapping `Pfam-A --cut_ga` domain clusters (`ov >= 30 aa`), sequences are assigned **strictly to the `CL0013` family with highest domain bitscore**, resolving the artificial `PF13354`/`PF00768` boundary conflict (`dual_PF13354_PF00768_GA` remains logged).
   - Previously, `764` viral `PF00768`-highest-scoring proteins (incl. `294 / 295` FIMO-rejected candidates carrying `SXXK` with median length `237 aa`) were misassigned to `PF13354`; under strict Clan bitscore, viral GA survivors cleanly split into **`PF00144 = 2,105`**, **`PF00905 = 1,331`**, **`PF13354 = 261`**, and sister PBP **`PF00768 = 895`**.
2. **Symmetric Own-Family Official-Seed `FIMO 5.5.9 (p <= 1e-5)`**:
   - Every sequence is filtered exclusively against its **own assigned family's** official Pfam seed MEME model across domain slices (`env_from..env_to ± 10 aa`), completely eliminating the prior `71% vs 88%` asymmetry:
     - `PF00144`: Viral `2,086 / 2,105 (99.10%)` | Bacterial `13,225 / 13,312 (99.35%)`
     - `PF00905`: Viral `1,323 / 1,331 (99.40%)` | Bacterial `16,080 / 16,203 (99.24%)`
     - `PF13354`: Viral **`260 / 261 (99.62%)`** | Bacterial `3,386 / 3,456 (97.97%)`
     - `PF00768`: Viral **`895 / 895 (100.00%)`** | Bacterial `3,889 / 3,896 (99.82%)`
     - `PF02113`: Bacterial `1,990 / 2,014 (98.81%)` (Viral `0`)
3. **HMM Coverage Filter (`hmm_cov >= 0.65` Primary; `>=0.60` / `>=0.70` Tracked)**:
   - `hmm_cov` = number of HMM match states covered by the winning family's domain / HMM `LENG`. hmmsearch can report one domain as several hits of the same family (a divergent, poorly scoring middle region, or a genuine internal insertion), so **co-linear same-family hits are merged** (no ≥30 aa overlap in the protein, same order in protein and HMM, each piece adding ≥50% new match states; candidates only from the winning hit's own envelope cluster or from clusters the same family wins, never from regions won by another family) and coverage is the union of their HMM intervals; the domain slice spans all pieces (insertions become lowercase insert states in hmmalign). A real internal deletion still lowers coverage. `hmm_cov_best_piece` and `n_domain_pieces` are kept per sequence.
   - Merging rescued `131` viral (`PF00768 116`, `PF00144 9`, `PF00905 3`, `PF13354 3`; median merged coverage `0.84`) and `461` bacterial sequences (`PF00905 279`, `PF13354 141`, `PF00144 27`, `PF00768 14`); multi-piece domains still `< 0.65`: `6` viral, `6` bacterial.
   - Viral sequences passing FIMO & `hmm_cov >= 0.65`: `PF00144 1,876`; `PF00905 1,242`; `PF13354 248`; `PF00768 837`.
4. **Triple-Source Bacterial Pool with Tagged Prophages (`in_bacterial_prophage`)**:
   - Combines `3,000` phylogenetically balanced GTDB/RefSeq representatives (`46` phyla, `837` families; `99.80%` 1:1 native CDS), `259` dereplicated Marine MAGs (`100%` 1:1 CDS), and `1,219` `CARD/SARG` clinical gold standards.
   - Loci inside geNomad provirus intervals are **retained with `in_bacterial_prophage = True` (`replicon = prophage`)** for paired `with_prophage` (Primary) vs. `without_prophage` (Sensitivity) testing.
   - Whole contigs called viral by geNomad inside a GTDB assembly or Marine MAG (most likely free phage genomes binned with the host, not integrated prophages) are tagged **`replicon = viral_contig` / `is_viral_contig_in_genome = True`**, kept in the catalog, and **excluded from every bacterial tree, HyPhy and outgroup set**.

### 2.4 Independent Domain-Slice `nr85` Clustering & SBL vs. PBP/PRE Functional Stratification
- **Independent `MMseqs2` (`--min-seq-id 0.85 -c 0.80 --cov-mode 0`)**:
  - Viral and bacterial domain slices are clustered **separately** (never mixed) to preserve recent `>=85%` cross-domain HGT pairs while removing clonal viral inflation; viral `nr85` rep selection prioritizes Phase 2 Strict `Tier 1 > MOSAiC > Tier 2 > Tier 3` and logs `viral_cluster_size` strictly for reporting (never leaf-weighted).
- **Functional Stratification (`CARD + hmmalign` Catalytic Triad)**:
  - Searches domain slices via `DIAMOND` against **CARD β-lactamases after removing PBP-type entries** (mecA/mecC PBP2a, PBP1a/2b/2x/3 resistance variants; removed headers listed in `work/card_excluded_PBP_entries.txt`) (`pident >= 40% & cov >= 70%`), alongside calibrated `hmmalign` match-state positions for **`SXXK`**, **`SDN` (Class A) / `YAN` (Class C) / `YGN` (Class D) / `SXN` (PBP)**, **`[KH][TS]G`**, and the **Class A Ω-loop signature `E166-x-x-[LIVM]-N`** within the calibrated window ± 3 match states.
  - Alignments are read in **A2M** format; each sequence must yield exactly `LENG` match states of its family HMM, otherwise the run stops (guards against mis-placed motif windows).
  - Separates **`Putative_Hydrolytic_BetaLactamase`** (`100%` recall on CARD/SARG Class A `PF13354` `632/632` and Class C `PF00144` `199/199`, while flagging `126` `SARG PBP2/MecA` non-SBLs as `PBP_or_PRE_Homolog`!) from **`PBP_or_PRE_Homolog`**:

| Family Group (`cov>=0.65 & FIMO<=1e-5`) | Viral Total (`SBL`) | Viral `nr85` (`SBL`, `MOSAiC`) | Bacterial Total | Bacterial `nr85` (`SBL`) | Tree Balanced `rep1` (`prophage`) |
|---|---:|---:|---:|---:|---:|
| **`PF00144`** | `1,876` (`44`) | **`416`** (`27`, `4`) | `12,851` | `12,149` (`1,008`) | `1,500` (`3`) |
| **`PF00905`** | `1,242` (`12`) | **`215`** (`6`, `1`) | `15,797` | `14,339` (`377`) | `1,500` (`1`) |
| **`PF13354`** | `248` (`146`) | **`63`** (`22`, `0`) | `2,674` | `2,001` (`476`) | `1,500` (`0`) |
| **`PF00768` (Sister PBP)** | `837` (`0`) | **`195`** (`0`, `1`) | `3,872` | `3,681` (`0`) | `1,500` (`4`) |
| **Joint `PF13354+PF00768`** | `1,085` (`146`) | **`258`** (`22`, `1`) | `6,546` | `5,681` (`476`) | `1,500` (`4`) |

- **Calibration of the stricter Ω-loop rule**: `E166-x-x-[LIVM]-N` is detected in `577 / 631` (`91.4%`) CARD/SARG Class A (`PF13354`) references; all `631 / 631` remain `Putative_Hydrolytic_BetaLactamase` (CARD match). Among non-CARD bacterial `PF13354`, `771` carry `SXXK + SDN + KTG` without the Ω-loop signature or a CARD match and are now `PBP_or_PRE_Homolog` (bacterial `nr85` SBL `930 → 476`).
- **CARD PBP filter**: `1` PBP-type entry removed from the CARD reference set (`4,895` kept).
- **Whole viral contigs in bacterial genomes** (`is_viral_contig_in_genome`, excluded from bacterial sets): `3` (`GTDB PF00905`, `GTDB PF13354`, `Marine MAG PF00905`).

- **Additional Sensitivity & Outgroup Exports**: `rep2..rep5` (4 taxonomy-stratified random seeds), `phylo_rep` (`MMseqs2 65%` deep-homology balance `<=1,500`), `rep1_no_prophage`, `cov60` & `cov70` sets, Sister-PRE outgroups, `PF02113` (`100` 60%-reps), and **all 3 non-PRE `CL0013` outgroups (`PF04960: 24`, `PF06089: 24`, `PF27310: 24`)**.

---

## 3. Stage 1 Phase 2 (Completed): Strict Structural Hallmark Genomic Context & Strict `Tier 1`

Integrates `298,926` Prodigal ORFs, `171,901` full-genome `Pfam-A --cut_ga` hits, `89,824` `eggNOG` annotations, `62,864` `geNomad v1.7` genes, and `CheckV v1.0.3` across `100%` (`4,193` IMG/VR + `10` MOSAiC) of primary curated viral loci (`Phase2_Viral_Genomic_Neighborhoods_Full.tsv`):
- **Recombinase/Integrase Strictly Separated from Viral Hallmarks**:
  - **Strict viral structural hallmarks** (capsid, portal, terminase, tail, baseplate, virion coat/sheath/fiber, prohead/scaffolding, `geNomad` non-recomb hallmarks): **`54.54%` (`2,287 / 4,193`)** in viral region; flanked on `>=1` side by strict structural hallmark **`54.47%` (`2,284 / 4,193`)**, on **BOTH sides `21.04%` (`882 / 4,193`)**.
  - **Viral structural OR replication/lysis (NO integrase/recomb)**: flanked on `>=1` side **`77.61%` (`3,254 / 4,193`)**, on **BOTH sides `44.69%` (`1,874 / 4,193`)**.
  - **Recombination / integrase / lysogeny (reported separately as MGE background)**: **`44.96%` (`1,885 / 4,193`)** — NEVER counted toward viral origin confidence.
  - Median min distance to viral/provirus boundary: **`5,160.0 bp`**; median contig ORF annotation coverage: **`47.00%`**.
- **`Tier1_HighConfidence_Viral`** (all required): SBL ORF `>= 1 kb` inside the viral region (geNomad/CheckV provirus interval, or contig ends for whole-contig viruses; DTR/ITR complete genomes exempt — contig length alone no longer counts) + `>= 1` ORF-verified non-recombinase strict viral structural hallmark inside the viral region (`n_viral_hallmark_in_region >= 1`; raw geNomad `n_hallmarks` is restricted to `Tier 2`) + CheckV `viral_genes >= 1` and host-gene check: (i) `viral_genes >= host_genes` on the CheckV viral region containing the SBL (or whole contig for non-provirus contigs), or (ii) `viral_genes >= host_genes - 1` when flanked on **BOTH** sides by strict viral structural hallmarks (accounting for CheckV self-counting the focal SBL/PBP as 1 host gene), or (iii) for `genomad_provirus_skip` intervals lacking a CheckV-split region, double-sided viral structural/replication flanking (`flanked_viral_structural_or_replic_both`) plus `viral_genes >= 1`. `Tier 2` likewise enforces the host-gene check (with double-sided strict structural flanking qualifying for `Tier 2` when `genomad_virus_score >= 0.70`). geNomad hallmark counts come from the geNomad entry that **contains the SBL ORF** (not the contig's best-scoring provirus). Each locus also reports its tier under the lenient comparison rule (`>= 2` hallmarks including holin/endolysin, `2 kb` boundary; column `viral_origin_tier_lenient_rule`, named `viral_origin_tier_2023rule` in outputs produced before 2026-11-24), and the summary prints a lenient → strict tier shift table.
- **Host-vs-viral gene check basis**: the CheckV viral region containing the SBL when one exists (`102` loci; `96` pass raw `vg >= hg`), skipped for geNomad-only provirus intervals whose contig-level counts include the host flank (`239`, requiring double-sided viral flanking for `Tier 1`), and contig-level for whole-contig viruses (`3,852`; **`1,634` fail** raw `vg >= hg`, i.e. CheckV finds more host than viral genes on the contig, unless rescued by double-sided strict structural flanking).
- **CheckV host-gene audit** (`scripts/hpc/02_genomic_context/02_audit_host_gene_check.py`, outputs `Phase2_Host_Gene_Check_Audit*.tsv/txt`): CheckV's own gene calls (`tmp/gene_features.tsv`, `hmm_cat` −1 host / 1 viral / 0 no hit) show that the focal SBL/PBP ORF is itself counted as a host gene in `2,852 / 4,193` loci (`68%`; PF00905 `1,065/1,241`, PF00768 `810/836`, PF13354 `247/248`, PF00144 `730/1,868`). The `host_genes - 1` allowance in rule (ii) above was set after the Stage 1 counts had been seen and is applied only to loci with double-sided strict structural flanking; `38` loci pass Tier 1/2 only through this flank exemption (`5` Tier 1, `33` Tier 2; `37/38` focal ORFs are CheckV host genes, and `19` have `>= 10` more host than viral genes). Alternative rules on the same loci (re-derived tiers reproduce the stored tiers exactly):
  - *Focal-corrected* (subtract the focal ORF only when CheckV calls it a host gene, for every locus, no flank exemption): Tier 1 `1,701`, Tier 2 `1,006`, Tier 3 `1,486` (current rule: `1,649` / `786` / `1,758`); PF00905 `136 / 299 / 806`.
  - *Focal-corrected + CheckV Complete/High-quality with 0 % contamination passes*: Tier 1 `1,764`, Tier 2 `1,022`, Tier 3 `1,407`; PF00905 `189 / 314 / 738`.
  - The primary tiers used by Stages 2–4 are still the current rule; switching rules changes the Tier-1 gate used by the CL0013 sampling (`03/07`), HGT singleton gating (`04`) and HyPhy set selection (`05/01`).
- Re-run counts (primary cohort, `cov >= 0.65 & FIMO <= 1e-5`):
  - `PF00144`: `Tier1 = 1,164`, `Tier2 = 325`, `Tier3 = 379`, `MOSAiC = 8`
  - `PF00905`: `Tier1 = 128`, `Tier2 = 215`, `Tier3 = 898`, `MOSAiC = 1` ← **`72.4%` of IMG/VR PF00905 loci fall in `Tier 3`** (mostly whole contigs where CheckV host genes > viral genes)
  - `PF13354`: `Tier1 = 125`, `Tier2 = 50`, `Tier3 = 73`
  - `PF00768`: `Tier1 = 232`, `Tier2 = 196`, `Tier3 = 408`, `MOSAiC = 1`
- **Lenient rule → strict rule** (primary IMG/VR loci): of `1,539` old `Tier 1`, `1,363` stay `Tier 1`, `52` → `Tier 2`, `124` → `Tier 3`; of `1,993` old `Tier 2`, `286` → `Tier 1`, `734` stay, `973` → `Tier 3`; all `661` old `Tier 3` stay `Tier 3`.

---

## 4. Pre-Locked Analysis Specification & Execution Roadmap (Stages 2–4)

### 4.1 Pre-Locked Primary HGT Analysis & OAT Robustness Protocol (Anti-Cherry-Picking Standard)
1. **Primary Analysis Configuration (Pre-Registered)**:
   - **MSA**: `hmmalign` (official family `.hmm`, lowercase insertions stripped) on independent viral `nr85` + bacterial `nr85` `rep1` (`hmm_cov >= 0.65`, tagged bacterial prophages retained).
   - **Primary Rooting**: **Sister-PRE biological outgroup** (verified against the **`CL0013` PRE superfamily tree** rooted on non-PRE outgroups `PF04960+PF06089+PF27310` with leave-one-outgroup-family-out checks & structure-informed alignment).
   - **Tree Support & HGT Inference (`IQ-TREE 2.4.0`, `-m MFP -mset LG,WAG,Q.pfam,VT,JTT -B 1000 --bnni -alrt 1000`)**:
     - Double-outgroup topological nesting (Irwin et al. 2022) + **concordant two-pass Fitch parsimony AND `PastML` marginal ML ancestral state reconstruction**;
     - Strict **`Tier-1 HGT`**: `UFBoot >= 95%` AND `SH-aLRT >= 80%` on **`min(P, G)`**;
     - Any **singleton-tip HGT** (`focal_is_singleton_tip`) **MUST pass Phase 2 Strict `Tier1_HighConfidence_Viral`** (`>=1` strict viral structural hallmark + CheckV viral gene + boundary verified).
2. **Unbiased Three-Family Rooting Cross-Validation Rule**:
   - *(1) Sister-PRE Outgroup* (Primary), *(2) MAD* on ingroup-only tree, *(3) Clock-like `MinVar`/`Midpoint` (counted as ONE clock-dependent family)*. **`Youden's J` is EXCLUDED** from primary cross-validation due to circularity (archived in Appendix only as a biased lower-bound comparator). Events corroborated under Sister-PRE + MAD + MinVar are classified `Root_Invariant`.
3. **One-Factor-at-a-Time (OAT) Robustness & Compute Budget**:
   - Sensitivity perturbations test: alternative MSA (`MAFFT 7.526 + ClipKIT 2.14.0`), alternative rooting (`MAD`, `MinVar`), alternative bacterial sampling (`rep2..rep5`, `phylo_rep` fixing `-m <best_from_rep1>` + Fast/ultrafast branch supports), `without_prophage`, `cov >= 0.60 / 0.70`, joint `PF13354_PF00768` tree, and **Functional Stratification (`Putative_Hydrolytic_BetaLactamase` vs. `PBP_or_PRE_Homolog`)**.
   - Each Primary HGT clade is tracked across OAT trees via tip-set Jaccard overlap (`>= 0.70`); events achieving **OAT Recurrence Rate `>= 80%`** are designated **`Robust_Core_HGT`**.
4. **Prophage-Borne Event Rule (Pre-Registered)**:
   - A `Viral → Bacterial` event whose bacterial recipient tips are **all** `in_bacterial_prophage = True` is reported as **`Prophage_Borne`** (the gene is carried by a viral genome integrated in a bacterium), **not** as transfer into the bacterial chromosome, and is excluded from the primary V→B count. Events with at least one chromosomal/plasmid recipient are counted normally, with the prophage fraction of recipients reported.
   - `is_viral_contig_in_genome` sequences never enter the trees (Section 2.3).
5. **Rooting Constraint for the PF13354/PF00768 Continuum**: PF13354, PF00768 and the joint `PF13354_PF00768` trees use **PF00144 + PF00905** Sister-PRE outgroups only; PF13354 and PF00768 are never used to root each other.

### 4.2 Remaining Stages
- **Stage 2 (`03_alignment_and_trees/` on HPC)** (inputs = Stage 1 job `5120575`; scripts `scripts/hpc/03_alignment_and_trees/01`–`08`):
  - **Alignments (job `5120670`)**: primary track = `hmmalign` against the family's official Pfam HMM with insert (lowercase) states removed (every row = `LENG` match columns); sensitivity track = `MAFFT --auto` + `ClipKIT -m gappy -g 0.5`. The joint `PF13354_PF00768` tree has no single family HMM, so its **primary track is MAFFT + ClipKIT**. Every tree alignment includes the Sister-PRE outgroup (`24` seqs); ingroup-only alignments drop those rows and all-gap columns without re-alignment (for MAD/MinVar).
  - `ClipKIT kpic-smart-gap` was tested first and **rejected**: on these insertion-rich domain slices it retained 733–2,235 columns at 16–35 % occupancy, and `kpic` modes delete singleton sites (shortening terminal branches used by MAD/MinVar rooting).
  - Alignment QC (`alignment_qc.tsv`, primary sets, ingroup rows): columns `239–329` (hmm) / `243–319` (MAFFT); mean occupancy `0.88–0.93` / `0.86–0.93`; SXXK lock-in `0.96–0.98` on both tracks; KTG lock-in `0.95–0.97` (PF00905, PF00768), `0.73` (PF13354), `0.27` (PF00144, where KTG is not conserved across the family).
  - **Trees**: wave 1 (job array `5121329`, 9 trees) = primary `_out` trees for both tracks with ModelFinder (`-m MFP -mset LG,WAG,Q.pfam,VT,JTT`, `-B 1000 --bnni -alrt 1000`, seed `20261121`); wave 2 (job array `5121330`, 41 trees) = ingroup-only, `rep2..rep5`, `phylo`, `nopro`, `cov60`, `cov70` (primary track) with the model fixed from the matching wave-1 tree. Trees are left unrooted (rooting in Stage 3). The first launch of both waves (`5120747`/`5120748`) failed at start-up (CRLF line endings in the task files) and was resubmitted unchanged. All 50 trees are finished; best-fit models (BIC): `Q.pfam+I+R10` (PF00144, PF00768), `Q.pfam+R10` (PF00905, PF13354, joint PF13354_PF00768).
  - **Outgroup monophyly check (primary trees)**: the 24-sequence Sister-PRE outgroup is monophyletic in the unrooted tree only for PF00905 (largest monophyletic outgroup subset `22/24`); PF00144 `8/24`, PF13354 `6/24`, PF00768 `5/24`, joint PF13354_PF00768 `8/24`. Sister-PRE rooting of those four trees is therefore taken from the largest monophyletic outgroup subset in the figures, and the rooting cross-validation (MAD, MinVar, CL0013 superfamily, clean GTDB outgroups) carries the weight for them.
  - Rooting/LBA checks: `PMSF (LG+C20+F+G)` on the five primary trees and the `CL0013` superfamily + leave-one-outgroup-family-out trees are finished (scripts `04`–`06`); the structure-guided (FoldMason) `CL0013` C20 trees and the clean GTDB sister-outgroup trees (scripts `07`–`08`) were running on 2026-11-24.
- **Stage 3 (`04_hgt_inference/` on HPC)**: Execute pre-locked Primary Fitch + `PastML` ML HGT pipeline, OAT Robustness Recurrence scoring, singleton strict `Tier 1` gating, and separate reporting of Hydrolytic SBL vs. PBP/PRE HGT counts and expanded `viral_cluster_size`.
- **Stage 4 (`05_selection_pressure/`, `06_structure_and_ecology/`, & full local repo migration)**: `PAL2NAL` codon alignments on verified CDS -> formal **`dS` saturation diagnostics** (`dS` distribution & Xia's test) -> `GARD` recombination partitioning & `HyPhy` (`aBSREL`/`MEME`/`FEL`) within low-divergence HGT/monophyletic subclades (`>=70%` nt id) + `PF00144` 2×2 downsampling benchmark -> 3D PDB biophysics & stratified ecology -> migrate all local `scripts/01–09`, `Table_S01–S13`, `Fig_S01–S05`, and `README.md` cleanly to the 2026 cohort.
