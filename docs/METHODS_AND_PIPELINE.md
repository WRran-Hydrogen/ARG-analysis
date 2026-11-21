# Viral Serine-β-Lactamase (SBL) Full Re-Analysis: Methods, Checkpoints & Session Memorize Reference

> **Purpose**: Authoritative, reproducible Materials & Methods specification **and** complete HPC/local execution handoff (`Memorize`) for the full upgrade of:
> *"Global analyses reveal the diversity and functional evolution of viral serine-β-lactamase genes linked to antibiotic resistance"*
> Last updated: `2026-11-21` (Phases 0–2 completed and verified on the HPC cluster; Phases 3–6 staged and ready for execution).

---

## 1. Environment, Compute Infrastructure & File Map

### 1.1 Compute Infrastructure & Verified Software Versions
- **HPC Cluster**: institutional Slurm cluster (partitions: `compute` [77 nodes, 7-day limit]; default QOS). Outbound HTTPS downloads (NCBI RefSeq / GTDB) executed from login node via `aria2c 1.37.0` + `curl`; all HMM/geNomad/CheckV/eggNOG/MMseqs2/tree/HyPhy workloads executed on `compute` nodes.
- **Software Versions & Install Locations on the HPC Cluster**:
  - **HMMER `3.4`**: `$HPC_HOME/.pixi/bin/hmmsearch`, `hmmalign`
  - **MMseqs2 `18.8cc5c`**: `$HPC_HOME/.pixi/bin/mmseqs`
  - **MEME Suite / FIMO `5.5.9`**: `$HPC_HOME/micromamba/envs/sbl/bin/meme`, `fimo`
  - **CheckV `v1.0.3`** (DB `checkv-db-v1.5`, 13 GB): `$HPC_HOME/micromamba/envs/checkV/bin/checkv`
  - **geNomad `1.7`** (`genomad_db` v1.7, 2.8 GB, IMG/VR viral screening) & **`1.8.1`** (bacterial 8-shard screening): `$HPC_HOME/micromamba/envs/genomad/bin/genomad`
  - **eggNOG-mapper** (`emapper.py -m diamond`, DB `$HPC_HOME/eggnog_db`): `$HPC_HOME/micromamba/envs/eggnog/bin/emapper.py`
  - **Prodigal `-p meta`** (`51cc034a`): `$HPC_HOME/.pixi/bin/prodigal`
  - **Phylogenetics & Selection Tools (Phases 3–5)**: `MAFFT 7.526` (`module load MAFFT/7.526`), `MUSCLE 5.3`, `ClipKIT 2.14.0`, `trimAl v1.5.rev0`, `pal2nal.pl`, `seqkit v2.14.0`, `IQ-TREE 2.4.0`, `HyPhy 2.5.93(MP)`, `Foldseek 10.941cd33`, `Python 3.13` (`envs/sbl`).

### 1.2 Master Workspace Directory Tree (HPC + Local)
- **HPC Root**: `$SBL_HPC_BASE/`
  - `inputs/` — Uploaded `PerSequence_Bacterial_Audit.tsv`, `PF*_SARG.faa`, `PF*_seed.faa`, `CL0013_official.hmm`, `PF*_official.hmm`, viral fallback `ARG_cds.fna`
  - `00_checkv/out/` — Fresh CheckV v1.0.3 (DB v1.5) output across all `5,618` viral contigs (`Job 5113663`)
  - `01_curated_datasets/` — Initial Scheme B viral decontamination (`Job 5117517`), `r1_meme/meme.txt`, and `mosaic_audit/` (`Job 5113697` + `clean_mosaic_euk.py`)
  - `01b_bacterial_ref_db/` — GTDB `3,000` representative bacterial genomes (`Job 5113993`), 8-shard `geNomad` (`Job 5117565`), and `marine_mags/` (`259` MAGs, `Job 5113811`)
  - `01_final_curated/` — **Authoritative final Phase 1 protein, domain-core, and CDS outputs** (`Job 5118649`)
  - `02_genomic_context/` — **100%-coverage Phase 2 viral genomic neighborhood & confidence tier catalog** (`Job 5118649`)
  - `03_alignment_and_trees/` .. `06_structure_and_ecology/` — Staged for Phases 3–6
- **Upstream HPC Annotation Jobs**:
  - `51cc034a` (Job `5113325`): Prodigal `-p meta` across all `5,618` `sbl_contigs.fna` (`340,422` ORFs)
  - `6ca11f3e` (Job `5113378`): `geNomad v1.7` end-to-end across all `5,618` `sbl_contigs.fna`
  - `de0940b0` (Job `5113374`, `5h44m42s`): Full-genome `Pfam-A.hmm` (53,275 HMMs) `--cut_ga` on `340,422` viral ORFs (`201,317` hits)
  - `8ee5379a` (Job `5113454`, `9h48m24s`): Full-genome `eggNOG-mapper` on `340,422` viral ORFs (`103,146` queries scanned, `103,110` annotated)
- **Local Self-Contained Master Repository (this repository)**:
  - `docs/METHODS_AND_PIPELINE.md` — Master methods & pipeline document
  - `docs/private/` — Peer-review correspondence (kept locally, not version-controlled)
  - `scripts/` (`01`–`09`, `run_pipeline.py --verify`) & `scripts/hpc_reanalysis/` (all 11 HPC Phase 1–2 pipeline scripts)
  - `data/curated_2026/` — All 39 tree-ready (`*_dom.faa`, `*_combined_tree_dom.faa`), 1:1 CDS-verified (`*_cds.fna`, `*_dom_cds.fna`), and `CL0013` auxiliary/outgroup FASTAs from `01_final_curated/`
  - `data/metadata/` (`information.xlsx`, `Pfam-A.clans.tsv`, `sum_high.*`, `no_HGT.*`, `Prophage_*`) & `data/structures/` (`ref_pdbs/` + `ProteinStructure/` 1,507 PDBs + `1ewz` CIF)
  - `results/reanalysis_2026/` — All 17 Phase 1/1B/2 audit catalogs, summary tables, and `Phase2_Viral_Genomic_Neighborhoods_Full.tsv`
  - `results/tables/`, `results/supplementary_data/`, `figures/` — Baseline `Table_S01–S13` and `Fig_S01–S05`.

---

## 2. Phase 1: Viral Decontamination & Domain-Level Curation (Completed)

### 2.1 Viral Sources & 1:1 Nucleotide CDS Mapping
1. **IMG/VR v4 Curated Baseline**: `5,701` viral SBL proteins encoded on `5,618` unique contigs (`sbl_contigs.fna`, covering `100.0%` of baseline final loci). Each protein was mapped back to its nucleotide CDS and genomic coordinates via:
   - Exact/substring amino-acid match against 16-part Prodigal `-p meta` ORFs (`51cc034a/orf/*.fna`, Table 11 + alternative start codons `GTG/TTG/ATT/ATC/CTG`): `5,358` exact + `314` substring;
   - Strict 1:1 translation-verified fallback against `ARG_cds.fna` (capturing `28` non-standard/GVMAG ORFs);
   - **Total verified viral CDS mapping**: **`5,700 / 5,701 (99.98%)`**.
2. **Arctic MOSAiC Viral Proteome Addition & Eukaryotic Purge**:
   - Original MOSAiC search (`$SBL_HPC_RAW_INPUT_DIR/sbl.hmm`, `-E 1e-5` vs `mosaic_virus.faa`, `650 MB`) yielded `sbl_mosaic.tsv` (`28 PF00144 + 1 PF00905`) and `mosaic_virus_pf00144.faa` (`29` seqs; union = `30`, zero overlap with IMG/VR `5,701`).
   - Re-scanning the full `mosaic_virus.faa` with official `CL0013_official.hmm --cut_ga` (Job `5113697`, union = `83` candidates) followed by full `Pfam-A.hmm` (53,275 HMMs) `--cut_ga` overlapping-envelope competition eliminated **`51` non-penicillin-recognizing `CL0013` enzymes** (`49` `PF04960:Glutaminase` — 43 from *Phycodnaviridae* — and `2` `PF06089:Asparaginase_II`) plus `1` below-GA hit, leaving `31` focal/aux hits (`28 PF00144`, `1 PF00905`, `2 PF00768`).
   - **Cross-domain co-binning purge (`clean_mosaic_euk.py`)**: Querying all `31` inferred MOSAiC contigs against `$SBL_HPC_RAW_INPUT_DIR/euk_info.tsv` detected and purged **`6` `PF00144` `unknown`-taxonomic entries** overlapping eukaryotic MAGs (`bacillariophyta`, `chrysophyceae`, `amoebozoa`, `euglenozoa`, `mbat_multi.43`), yielding **`25` clean MOSAiC viral SBLs** (`22 PF00144` [incl. `15` `caudoviricetes` vMAGs], `1 PF00905`, `2 PF00768`). Because MOSAiC contig `.fna` bins are not staged in `$SBL_HPC_RAW_INPUT_DIR`, MOSAiC sequences are included in Phase 3/4 protein phylogenies & HGT analyses (`contig_context_available = False`) and excluded from nucleic-acid-only HyPhy runs.

### 2.2 Scheme B Joint Viral Decontamination (`geNomad v1.7` + `CheckV v1.0.3`)
To strictly prevent bacterial chromosomal flanking contamination or plasmid false positives (Comments #9, #10), every IMG/VR viral SBL locus (`n = 5,701`) was filtered at **gene-coordinate resolution** using both `geNomad v1.7` (`6ca11f3e`) and fresh `CheckV v1.0.3` (`checkv-db-v1.5`, `5113663`):
1. **`DROP_genomad_plasmid` (`n = 8`)**: Contig classified in `sbl_contigs_plasmid_summary.tsv` (`plasmid_score = 0.704–0.960`).
2. **`DROP_genomad_host_flank` (`n = 126`)**: Contig contains a `geNomad` provirus region (`401` contigs / `406` SBLs), but the SBL ORF coordinates fall strictly **outside** the provirus boundaries (in the bacterial host flank).
3. **`DROP_checkv_host_flank` (`n = 494`)**: Contig flagged `provirus == Yes` in CheckV `contamination.tsv`, and the SBL ORF has zero or minority overlap with CheckV `viral` regions (residing in a CheckV `host` region). Together, geNomad + CheckV excised **`620` bacterial host-flank SBL loci**.
4. **`DROP_nonviral_chrom_or_plasmid` (`n = 163`)**: Contig absent from `geNomad virus_summary` and highest aggregated geNomad score is `chromosome_score` or `plasmid_score`.
5. **Scheme B Fragment Rescue Rule** for contigs absent from `virus_summary` whose highest score is `virus_score`:
   - **`KEEP_Tier3_CheckV_Rescued` (`n = 4`)**: Retained ONLY IF `checkv_viral_genes >= 1` AND `checkv_host_genes == 0`;
   - **`DROP_nonviral_unverified_fragment` (`n = 215`)**: Otherwise excluded.
- **Total IMG/VR viral SBLs surviving Scheme B**: **`4,695`** (`4,691 KEEP_VIRAL` + `4 Tier3_Rescued`, out of `5,701`; `1,006` excluded).

### 2.3 Official Pfam-A `--cut_ga` Overlapping-Envelope & Focal-Family Priority Rule
- All surviving viral (`4,695` IMG/VR + `25` MOSAiC) and bacterial proteins were searched against full `Pfam-A.hmm` (53,275 HMMs) under official gathering thresholds (`--cut_ga`).
- Domain envelopes overlapping by `>= 30 aa` were clustered into connected components (`resolve_focal_priority`), **preserving genuine non-overlapping multi-domain fusions** (e.g., `PF00912:Transgly` + `PF00905:Transpeptidase`, phage tail/CBM fusions + `PF00144`) while isolating non-`CL0013` false hits (`drop_no_pfam_ga: 82`, `drop_non_CL0013_only: 46` on IMG/VR).
- **Focal-Family Priority within Clan `CL0013` clusters**: Any `CL0013`-winning cluster containing at least one `--cut_ga`-passing hit to the three focal SBL families (`PF00144.31` GA=27.0, `PF00905.29` GA=22.1, `PF13354.14` GA=22.0) was retained in its focal family (`focal_{acc}`), recording dual-threshold hits to `PF00768.26:Peptidase_S11` (GA=27.6) in **`dual_PF13354_PF00768_GA`**.
  - Result on `4,695` IMG/VR viral survivors: `focal_PF00144 = 2,083`, `focal_PF00905 = 1,330`, `focal_PF13354 = 1,017` (of which **`764 / 1,017 = 75.12%`** pass BOTH `PF13354` and `PF00768` `--cut_ga`), `aux_PF00768 = 137` (passing `PF00768` GA only, routed to `CL0013` auxiliary group).
  - Result on `25` MOSAiC survivors: `focal_PF00144 = 22`, `focal_PF00905 = 1`, `aux_PF00768 = 2`.
- Catalytic domain cores (`*_dom.faa`) and codon-frame nucleotide slices (`*_dom_cds.fna`) were extracted uniformly using **`[max(1, env_from - 10) .. min(seq_len, env_to + 10)]`**.

---

## 3. Phase 1B: Triple-Source Phylogenetically Balanced Bacterial Reference Library (Completed)

To simultaneously resolve Comment #10 (exact prophage removal), Comment #7 (phylum-balanced bacterial host spectrum), Comment #8 (100% symmetric HMM+FIMO curation), and verified bacterial CDS availability for Phase 5 HyPhy, legacy `NCBI nr` bacterial fragments were replaced by a **Triple-Source Bacterial Reference Pool**:

### 3.1 Source 1 — GTDB × RefSeq `3,000` Species-Representative Bacterial Genomes
1. **Selection (`select_gtdb_genomes.py`)**: From GTDB latest `bac120_metadata.tsv.gz` joined to NCBI RefSeq `assembly_summary.txt` (`39,085` `RS_GCF_` bacterial species representatives -> `32,861` passing `completeness >= 95%` and `contamination <= 2%`), two-stage balanced sampling selected **`3,000` distinct bacterial genera** spanning **`46` phyla, `115` classes, `333` orders, and `837` families** (`1,830 Complete Genome + 71 Chromosome + 741 Scaffold + 358 Contig`).
   - Top phyla: `Pseudomonadota` (637), `Bacillota` (633 + `Bacillota_I` 77), `Bacteroidota` (474), `Actinomycetota` (472), `Cyanobacteriota` (177), `Desulfobacterota` (108), `Planctomycetota` (53), `Verrucomicrobiota` (52), `Campylobacterota` (46), `Spirochaetota` (34), `Chloroflexota` (23), `Acidobacteriota` (22), `Fusobacteriota` (21).
2. **Proteome & CDS Verification (`scan_gtdb_sbls.py`, Job `5113993`)**:
   - Scanned **`10,280,376` proteins across `3,000 / 3,000` proteomes** (`CL0013_official.hmm --cut_ga` -> `37,241` hits -> full `Pfam-A --cut_ga` across `2,959` positive genomes);
   - **1:1 translation-verified CDS match**: **`37,163 / 37,239 (99.80%)`** directly matched native `_cds_from_genomic.fna.gz`.
3. **8-Shard Parallel `geNomad 1.8.1` Prophage Excision (`Job 5117565_0..7`)**:
   - Downloaded `2,959 / 2,959` `_genomic.fna.gz` assemblies; extracted all **`10,082 / 10,082` SBL-bearing contigs (`9.61 Gb`)** into 8 length-balanced shards (`1,201.58 Mb` each).
   - `geNomad 1.8.1` identified `1,430` provirus-bearing contigs, `2` full-viral contigs, and `147` plasmid contigs; coordinate intersection **physically excised `13` bacterial loci inside viral/prophage regions**.
   - Focal-priority surviving counts: `focal_PF00144 = 11,689`, `focal_PF00905 = 14,829`, `focal_PF13354 = 6,359`, `aux_PF00768 = 638`, `aux_PF02113 = 1,639`, biological outgroups (`PF04960 = 1,191`, `PF06089 = 843`, `PF27310 = 38`).

### 3.2 Source 2 — `259` Dereplicated High-Quality Marine Bacterial MAGs (`Job 5113811`)
- From `$SBL_HPC_MARINE_MAGS_DIR/` (`CheckM2` completeness `>=90%`, contamination `<5%`, 99% ANI dereplicated), `2,692` raw candidates were audited under full `Pfam-A --cut_ga`, yielding **`2,659 / 2,659 (100.00%)`** native contig 1:1 CDS verification.
- Excising `2` full-viral `geNomad` contigs and `162` edge/incomplete partial ORFs retained **`2,495` clean marine bacterial SBLs** (`PF00144 = 1,422`, `PF00905 = 995`, `PF13354 = 78`, with `chromosome`/`plasmid` tags).

### 3.3 Source 3 — `CARD/SARG` Clinical Gold Standards (`n = 1,219`)
- Retained `1,219` non-`NCBI_nr` `CleanCore` clinical references (`PF00144 = 199`, `PF00905 = 388`, `PF13354 = 632`) with explicit `plasmid_or_MGE` / clinical β-lactamase tags (so genuine plasmid-borne clinical β-lactamases such as `TEM/SHV/CTX-M/OXA/CMY/FOX/DHA` remain represented).

---

## 4. Motif Filtering Design & Final Deduplicated Counts (`build_final_phase1_datasets.py`, Job `5118649`)

### 4.1 Why One Round vs. Two Rounds of `MEME/FIMO`? (Empirical Methodological Note)
- **Original rationale & empirical diagnosis**: Previously, Round 1 ran MEME (`-protein -nmotifs 5 -minw 6 -maxw 20 -mod zoops -maxsize 1000000`) on official Pfam seed sequences (`PF00144_seed.faa`: 127; `PF00905_seed.faa`: 40; `PF13354_seed.faa`: 76), and Round 2 re-trained MEME on selected viral + bacterial Round-1 hits (`p <= 1e-5`) to expand motifs onto divergent environmental homologs.
- **Empirical test on our decontaminated dataset revealed**:
  1. For **`PF00144`** and **`PF00905`**, **a single round (`Round 1` official-seed `FIMO p <= 1e-5`) is already completely sufficient and avoids Round-2 overfitting** (`PF00144`: `99.10%` viral / `99.35%` bacterial pass R1, whereas naive top-bitscore Round-2 selection overfit a narrow sub-clade down to `82` viral hits).
  2. For **`PF13354`**, `722 / 1,017 (70.99%)` viral and `6,201 / 7,069 (87.72%)` bacterial hits pass Round 1 at `p <= 1e-5`; the remaining **`295` viral `PF13354` `--cut_ga` hits** cluster just above the cutoff (`median p = 1.49e-5`, IQR `1.16e-5–3.44e-5`), because `294/295` belong to the **`dual_PF13354_PF00768_GA`** PBP-bridging clade (not covered by Pfam's 76 classical bacterial seeds) despite `293/295 (99.3%)` carrying canonical `SXXK`.
- **Uniform Policy Applied in `01_final_curated/`**: Both viral and Triple-Source bacterial focal pools were filtered identically by **single-round official-Pfam-seed `FIMO (p <= 1e-5)`** (no hard-coded `SXXK` regex pre-filter), while logging `p <= 1e-10` and `PF13354` `dual_PF00768` borderline candidates for optional sensitivity checks.

### 4.2 Final Exported Phase 1 Counts & Bacterial `MMseqs2` Deduplication
Bacterial FIMO passers were clustered with `MMseqs2 18.8cc5c` (`mmseqs easy-cluster --min-seq-id 0.85 -c 0.80 --cov-mode 0`) selecting 1 representative per 85% cluster (`CARD/SARG` > `GTDB` complete/CDS-verified > `MarineMAG` CDS-verified, also exporting dedicated 100%-CDS-verified `*_hyphy_dom_cds.fna` cluster reps + `tree_balanced` sets capped at `1,500` preserving all `CARD/SARG` reps and all bacterial families/phyla):

| Family | Viral Clean (IMG/VR + MOSAiC) | Verified Viral CDS | Bacterial Raw FIMO `<=1e-5` | Bacterial `nr85` Clusters | Tree-Balanced Bacterial (`phyla`) | HyPhy CDS-Verified in Balanced (`%`) |
|---|---:|---:|---:|---:|---:|---:|
| **`PF00144`** | **`2,086`** (`2,065` + `21`) | `2,065 / 2,065` | `13,223` (`99.35%`) | `12,516` | **`1,500`** (`39` phyla) | `1,483` (`98.9%`) |
| **`PF00905`** | **`1,323`** (`1,322` + `1`) | `1,322 / 1,322` | `16,089` (`99.24%`) | `14,594` | **`1,500`** (`57` `p__` taxa) | `1,454` (`96.9%`) |
| **`PF13354`** | **`722`** (`722` + `0`; `470` dual `PF00768`) | `721 / 722` | `6,201` (`87.72%`) | `5,337` | **`1,500`** (`36` phyla) | `1,457` (`97.1%`) |
| **Focal Sum** | **`4,131`** (`4,109` + `22`) | **`4,108 / 4,109`** | **`35,513`** | **`32,447`** | **`4,500`** | **`4,394 / 4,500`** |
| **`CL0013` Aux (`PF00768`/`PF02113`)** | **`139`** (`137` + `2` `PF00768`) | — | `2,277` (`638+1,639`) | `300` (60%-rep) | `300` | — |
| **Outgroup (`PF04960:Glutaminase`)** | `0` | — | `1,191` | `24` (50%-rep) | `24` | — |

---

## 5. Phase 2: 100%-Coverage Viral Genomic Neighborhoods (Completed, `Job 5118649`)

Script `$SBL_HPC_BASE/02_genomic_context/build_phase2_neighborhoods.py` integrates **all five genome-wide evidence layers** across `100%` of final curated IMG/VR viral loci (`n = 4,246` = `4,109` focal + `137` `aux_PF00768`, plus `24` MOSAiC labeled `contig_context_available = False`; outputs `Phase2_Viral_Genomic_Neighborhoods_Full.tsv` & `Phase2_Genomic_Context_Summary.txt`):
- Mapped `163,012` full-genome `Pfam-A --cut_ga` domain hits, `85,869` `eggNOG` functional annotations, and `57,703` `geNomad` gene hallmarks onto target contig Prodigal ORFs;
- **Median contig ORF annotation coverage**: **`46.79%`** (`Pfam-A` + `eggNOG`);
- **Viral hallmark & flanking context**: **`57.49%` (`2,441 / 4,246`)** encode `>=1` viral structural hallmark (capsid, portal, terminase, tail, holin/endolysin) inside the viral region; **`81.28%` (`3,451 / 4,246`)** are flanked on at least one side and **`47.36%` (`2,011 / 4,246`)** on BOTH sides by viral/recombination genes; median minimum distance to the viral/provirus boundary is **`5,042.0 bp`**;
- **Upgraded Viral Origin Confidence Tiers**:
  - `PF00144` (`2,065` IMG/VR + `21` MOSAiC): `Tier1 = 1,055`, `Tier2 = 833`, `Tier3 = 177`
  - `PF00905` (`1,322` IMG/VR + `1` MOSAiC): `Tier1 = 209`, `Tier2 = 669`, `Tier3 = 444`
  - `PF13354` (`722` IMG/VR): `Tier1 = 245`, `Tier2 = 395`, `Tier3 = 82`
  - `aux_PF00768` (`137` IMG/VR + `2` MOSAiC): `Tier1 = 17`, `Tier2 = 78`, `Tier3 = 42`

---

## 6. Execution Blueprint & Design Rules for Phases 3–6 (Next Steps)

1. **Phase 3 — Alignment Benchmarking & Phylogenetics (`03_alignment_and_trees/`)**:
   - **Sequence Input**: Catalytic domain cores (`*_dom.faa`, `env_from..env_to ± 10 aa`) for each focal family (`Viral Clean + Bacterial Tree-Balanced`). Exclude extreme incomplete fragments (`<50%` aligned core coverage).
   - **MSA & Trimming Benchmark**: Compare `hmmalign` (using official family `.hmm`, **stripping lowercase insertion columns**) vs. `MAFFT 7.526 --auto` / `--localpair`, crossed with `ClipKIT 2.14.0 (-m smart-gap)` vs. `trimAl v1.5.rev0 (-gappyout / -gt 0.5)`; verify catalytic `SXXK` / `SXN` / `KTG` column conservation (`>=85%`) before selecting primary MSA.
   - **Tree Inference (`IQ-TREE 2.4.0`)**: Run BOTH:
     - *(a) Ingroup-only trees* (for unbiased **MAD** & **Youden's *J*** statistical rooting without LBA from distant outgroups);
     - *(b) Biological outgroup-rooted family trees* (with `PF00768`/`PF02113` or `PF04960` outgroup) + **`CL0013` superfamily auxiliary tree** (`PF00144 + PF00905 + PF13354 + PF00768 + PF02113`, rooted on `PF04960:Glutaminase`);
     - Standard model parameters: `-m MFP -mset LG,WAG,Q.pfam,VT,JTT -B 1000 --bnni -alrt 1000 -T 32 --seed 20261120 --safe`, PLUS **PMSF LBA sensitivity check** (`-m Q.pfam+C20+F+G -ft guide.treefile`).
2. **Phase 4 — Upgraded HGT Inference (`04_hgt_inference/`)**:
   - Apply the **4 approved upgrades** over baseline `05_phylogenetic_hgt_inference.py`:
     1. **Dual Support + Dual Bounding Node Strict Tier (`Tier-1 HGT`)**: Require **`UFBoot >= 95%` AND `SH-aLRT >= 80%`** simultaneously on **`min(P, G)`** (`P` and `G` wrapping nodes), alongside moderate (`max >= 70%`) sensitivity;
     2. **Singleton-Tip × Phase 2 Genomic Flank Linkage**: Automatically cross-validate every `focal_is_singleton_tip` event against `Phase2_Viral_Genomic_Neighborhoods_Full.tsv` (`Tier1/2` hallmarks & boundary distance) and flag unsupported singletons;
     3. **Multi-Rooting Cross-Validation**: Cross-check Biological Outgroup vs. **MAD (Minimum Ancestor Deviation)** vs. **Youden's *J***;
     4. **Methods Citation Precision**: Document explicitly as **Irwin et al. (2022) double-outgroup nesting + two-pass Fitch maximum parsimony** aligned with **Dmitrijeva et al. (2024)** independent-event counting logic (without misclaiming RANGER-DTL).
3. **Phase 5 — Co-Evolution & Selection Pressure (`05_selection_pressure/`)**:
   - Use **100%-translation-verified CDSs** (`*_viral_dom_cds.fna` + `*_bacterial_nr85_hyphy_dom_cds.fna`), deduplicated via `MMseqs2` (`80%–85%` aa identity, retaining HGT clade reps), mapped 1:1 to `hmmalign` match-state codon alignments (strip lowercase insertions cleanly per codon triplet, benchmark `trimAl`), matching `IQ-TREE` subtrees.
   - Run `HyPhy 2.5.93(MP)`: `GARD` / `PHI`, `aBSREL` (background-normalized & paired root-to-tip trajectories), `MEME`, formal `FEL`, and the **2×2 orthogonal `PF00144` MEME benchmark** (downsampling `N=155 × 3` replicates × `hmmalign` realignment) to partition statistical branch-power vs. alignment effects without presupposing a single driver.
4. **Phase 6 — 3D Structural Biophysics & Confounder-Controlled Ecology (`06_structure_and_ecology/`)**:
   - Re-map newly identified `hmmalign`-anchored `MEME` & `FEL` sites onto `110` PDB structures (RSA, loop/coil fraction, distance to catalytic Ser), recompute stratified Faith's `PD` / 3D dispersion / `CMH` habitat–motif tests, and sync final tables/figures back into this repository.
