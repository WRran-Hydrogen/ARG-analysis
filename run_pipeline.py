#!/usr/bin/env python3
"""Master verification and execution entrypoint for SBL_research.

Usage:
  python run_pipeline.py --verify     # Verify repository layout, data, scripts, tables, and figures
  python run_pipeline.py --run-fast   # Re-run downstream scripts 03-09 (~2-3 min)
  python run_pipeline.py --run-all    # Re-run complete pipeline scripts 01-09 (includes MUSCLE/FIMO)
"""

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

REQUIRED_DATA = [
    "data/processed/hmm_rebuilt/PF00144.hmm",
    "data/processed/hmm_rebuilt/PF00905.hmm",
    "data/processed/hmm_rebuilt/PF13354.hmm",
    "data/hmm_models/CL0013_official.hmm",
    "data/hmm_models/PF00144_official.hmm",
    "data/hmm_models/PF00905_official.hmm",
    "data/hmm_models/PF13354_official.hmm",
    "data/hmm_models/PF00768_official.hmm",
    "data/hmm_models/PF02113_official.hmm",
    "data/hmm_models/PF04960_official.hmm",
    "data/hmm_models/PF06089_official.hmm",
    "data/hmm_models/PF27310_official.hmm",
    "data/processed/hmm_rebuilt/PF00144_seed.aln.faa",
    "data/processed/hmm_rebuilt/PF00905_seed.aln.faa",
    "data/processed/hmm_rebuilt/PF13354_seed.aln.faa",
    "data/phylogeny/PF00144.nwk",
    "data/phylogeny/PF00905.nwk",
    "data/phylogeny/PF13354.nwk",
    "data/phylogeny/maxcc.nwk",
    "data/motifs/motif_node.csv",
    "data/motifs/motif_edge.csv",
    "data/selection/pf00144.codon_trim.fna",
    "data/selection/pf00144.nwk",
    "data/selection/pf00144_info.tsv",
    "data/selection/pf00905.codon_trim.fna",
    "data/selection/pf00905.nwk",
    "data/selection/pf00905_info.tsv",
    "data/selection/pf13354.codon_trim.fna",
    "data/selection/pf13354.nwk",
    "data/selection/pf13354_info.tsv",
]

for fam in ("PF00144", "PF00905", "PF13354"):
    REQUIRED_DATA.extend([
        f"data/sequences/{fam}_seed.faa",
        f"data/sequences/{fam}_Viral.faa",
        f"data/sequences/{fam}_Viral_FIMO_r1.faa",
        f"data/sequences/{fam}_Viral_Final.faa",
        f"data/sequences/{fam}_Bacterial.faa",
        f"data/processed/bacterial_clean/{fam}_Bacterial_CleanStrict.faa",
        f"data/processed/bacterial_clean/{fam}_Bacterial_CleanCore.faa",
        f"data/curated_2026/{fam}_combined_tree_cov60_dom.faa",
        f"data/curated_2026/{fam}_combined_tree_cov70_dom.faa",
    ])

REQUIRED_DATA.extend([
    "data/sequences/PF00768_seed.faa",
    "data/sequences/PF02113_seed.faa",
])

for tag in ("PF00144", "PF00905", "PF13354", "PF00768", "PF13354_PF00768"):
    REQUIRED_DATA.extend([
        f"data/curated_2026/{tag}_viral_clean.faa",
        f"data/curated_2026/{tag}_viral_dom.faa",
        f"data/curated_2026/{tag}_viral_cds.fna",
        f"data/curated_2026/{tag}_viral_dom_cds.fna",
        f"data/curated_2026/{tag}_viral_nr85_clean.faa",
        f"data/curated_2026/{tag}_viral_nr85_dom.faa",
        f"data/curated_2026/{tag}_viral_nr85_cds.fna",
        f"data/curated_2026/{tag}_viral_nr85_dom_cds.fna",
        f"data/curated_2026/{tag}_bacterial_nr85_dom.faa",
        f"data/curated_2026/{tag}_bacterial_nr85_dom_cds.fna",
        f"data/curated_2026/{tag}_bacterial_nr85_hyphy_dom.faa",
        f"data/curated_2026/{tag}_bacterial_nr85_hyphy_dom_cds.fna",
        f"data/curated_2026/{tag}_bacterial_tree_clean.faa",
        f"data/curated_2026/{tag}_bacterial_tree_dom.faa",
        f"data/curated_2026/{tag}_bacterial_tree_cds.fna",
        f"data/curated_2026/{tag}_bacterial_tree_dom_cds.fna",
        f"data/curated_2026/{tag}_combined_tree_dom.faa",
        f"data/curated_2026/{tag}_combined_tree_nopro_dom.faa",
        f"data/curated_2026/{tag}_combined_tree_phylo_dom.faa",
        f"data/curated_2026/{tag}_combined_tree_rep2_dom.faa",
        f"data/curated_2026/{tag}_combined_tree_rep3_dom.faa",
        f"data/curated_2026/{tag}_combined_tree_rep4_dom.faa",
        f"data/curated_2026/{tag}_combined_tree_rep5_dom.faa",
        f"data/curated_2026/{tag}_sister_pre_outgroup_dom.faa",
    ])

REQUIRED_DATA.extend([
    "data/curated_2026/CL0013_aux_PF02113_bacterial_rep_dom.faa",
    "data/curated_2026/CL0013_outgroup_PF04960_clean.faa",
    "data/curated_2026/CL0013_outgroup_PF04960_dom.faa",
    "data/curated_2026/CL0013_outgroup_PF06089_clean.faa",
    "data/curated_2026/CL0013_outgroup_PF06089_dom.faa",
    "data/curated_2026/CL0013_outgroup_PF27310_clean.faa",
    "data/curated_2026/CL0013_outgroup_PF27310_dom.faa",
    "data/metadata/information.xlsx",
    "data/metadata/Color_Mapping_Guide.txt",
    "data/metadata/Prophage_identify.blastnout",
    "data/metadata/Inphage_gggenome.xlsx",
    "data/metadata/Pfam-A.clans.tsv",
    "data/metadata/sum_high.faa",
    "data/metadata/sum_high.pfam",
    "data/metadata/no_HGT.faa",
    "data/metadata/no_HGT.pfam",
    "data/structures/ref_pdbs/1bls_P99.pdb",
    "data/structures/ref_pdbs/1k25_pneumoniae_PF00905.pdb",
    "data/structures/ref_pdbs/1bue_ENTEROBACTER_PF13354.pdb",
])

# Large inputs distributed via Figshare instead of git (needed only by scripts 08 and 09).
OPTIONAL_EXTERNAL = [
    "data/structures/ProteinStructure",
]

SCRIPTS_ALL = [
    "scripts/01_build_and_verify_hmm.py",
    "scripts/02_fimo_threshold_sensitivity.py",
    "scripts/03_bacterial_prophage_audit.py",
    "scripts/04_viral_genomic_context_audit.py",
    "scripts/05_phylogenetic_hgt_inference.py",
    "scripts/06_selection_absrel_and_lineages.py",
    "scripts/07_selection_meme_and_catalytic_motifs.py",
    "scripts/08_structural_localization_3d.py",
    "scripts/09_confounder_controlled_ecology.py",
]

SCRIPTS_HPC = [
    "scripts/hpc/01_curation/01_viral_decontamination.py",
    "scripts/hpc/01_curation/01_viral_decontamination.sh",
    "scripts/hpc/01_curation/02_mosaic_extract_candidates.py",
    "scripts/hpc/01_curation/03_mosaic_parse_audit.py",
    "scripts/hpc/01_curation/04_mosaic_clean_euk.py",
    "scripts/hpc/01_curation/05_gtdb_select_genomes.py",
    "scripts/hpc/01_curation/06_gtdb_scan_sbls.py",
    "scripts/hpc/01_curation/07_gtdb_extract_contigs_and_shard.py",
    "scripts/hpc/01_curation/08_marine_mags_process.py",
    "scripts/hpc/01_curation/09_build_curated_datasets.py",
    "scripts/hpc/02_genomic_context/01_build_neighborhoods.py",
    "scripts/hpc/02_genomic_context/02_audit_host_gene_check.py",
    "scripts/hpc/03_alignment_and_trees/01_build_alignments.py",
    "scripts/hpc/03_alignment_and_trees/02_run_iqtree.sh",
    "scripts/hpc/03_alignment_and_trees/03_submit_primary_trees.sh",
    "scripts/hpc/03_alignment_and_trees/04_build_cl0013_and_pmsf_tasks.py",
    "scripts/hpc/03_alignment_and_trees/05_run_pmsf.sh",
    "scripts/hpc/03_alignment_and_trees/06_submit_cl0013_and_pmsf.sh",
    "scripts/hpc/03_alignment_and_trees/07_build_structural_cl0013_and_clean_outgroups.py",
    "scripts/hpc/03_alignment_and_trees/08_run_cl0013_c20.sh",
    "scripts/hpc/04_hgt_inference/01_hgt_inference.py",
    "scripts/hpc/04_hgt_inference/02_run_hgt_inference.sh",
    "scripts/hpc/04_hgt_inference/03_run_au_test.sh",
    "scripts/hpc/04_hgt_inference/04_summarize_au_tests.py",
    "scripts/hpc/04_hgt_inference/05_submit_au_tests.sh",
    "scripts/hpc/05_selection/01_build_codon_and_saturation.py",
    "scripts/hpc/05_selection/02_run_hyphy.sh",
    "scripts/hpc/05_selection/03_summarize_hyphy.py",
    "scripts/hpc/05_selection/04_submit_selection.sh",
    "scripts/hpc/hpc_config.py",
    "scripts/hpc/run_curation_and_context.sh",
]

REQUIRED_REANALYSIS_2026 = [
    "results/reanalysis_2026/Final_Phase1_Dataset_Summary.txt",
    "results/reanalysis_2026/Phase2_Genomic_Context_Summary.txt",
    "results/reanalysis_2026/Final_Viral_SBL_Catalog.tsv",
    "results/reanalysis_2026/Final_Bacterial_TripleSource_Catalog.tsv",
    "results/reanalysis_2026/Phase2_Viral_Genomic_Neighborhoods_Full.tsv",
    "results/reanalysis_2026/Phase1_Viral_Decontamination_Audit.tsv",
    "results/reanalysis_2026/Phase1_Summary_Report.txt",
    "results/reanalysis_2026/Phase1_Motif_Filter_Audit.tsv",
    "results/reanalysis_2026/Phase1_MultiDomain_Architecture_Catalog.tsv",
    "results/reanalysis_2026/Phase1_Bacterial_Metadata.tsv",
    "results/reanalysis_2026/MOSAiC_SBL_Clean_Audit.tsv",
    "results/reanalysis_2026/MOSAiC_Clean_Summary.txt",
    "results/reanalysis_2026/Marine_MAGs_SBL_Audit.tsv",
    "results/reanalysis_2026/Marine_MAGs_Summary.txt",
    "results/reanalysis_2026/bac_gtdb_scan_summary.txt",
    "results/reanalysis_2026/selection_summary.txt",
    "results/reanalysis_2026/selected_3000_genomes.tsv",
    "results/reanalysis_2026/primary_trees/PF00144.primary.hmm_out.treefile",
    "results/reanalysis_2026/primary_trees/PF00144.primary.hmm_out.iqtree",
    "results/reanalysis_2026/primary_trees/PF00905.primary.hmm_out.treefile",
    "results/reanalysis_2026/primary_trees/PF00905.primary.hmm_out.iqtree",
    "results/reanalysis_2026/primary_trees/PF13354.primary.hmm_out.treefile",
    "results/reanalysis_2026/primary_trees/PF13354.primary.hmm_out.iqtree",
    "results/reanalysis_2026/primary_trees/PF00768.primary.hmm_out.treefile",
    "results/reanalysis_2026/primary_trees/PF00768.primary.hmm_out.iqtree",
    "results/reanalysis_2026/primary_trees/PF13354_PF00768.primary.mafft_out.treefile",
    "results/reanalysis_2026/primary_trees/PF13354_PF00768.primary.mafft_out.iqtree",
    "results/reanalysis_2026/primary_trees/PF00144_outgroup.txt",
    "results/reanalysis_2026/primary_trees/PF00905_outgroup.txt",
    "results/reanalysis_2026/primary_trees/PF13354_outgroup.txt",
    "results/reanalysis_2026/primary_trees/PF00768_outgroup.txt",
    "results/reanalysis_2026/primary_trees/PF13354_PF00768_outgroup.txt",
    "results/reanalysis_2026/primary_trees/Tree_Tip_Annotations_Merged.tsv",
    "results/reanalysis_2026/hgt_inference/Primary_HGT_Events_Catalog_2026.tsv",
    "results/reanalysis_2026/hgt_inference/OAT_Sensitivity_Summary_Matrix_2026.tsv",
    "results/reanalysis_2026/hgt_inference/HGT_Summary_By_Family_And_Function_2026.tsv",
    "results/reanalysis_2026/hgt_inference/Rooting_Decisions_And_Agreement.tsv",
    "results/reanalysis_2026/hgt_inference/CL0013_Superfamily_Rooting_Validation.tsv",
    "results/reanalysis_2026/hgt_inference/Stage3_HGT_Summary_Report.txt",
    "results/reanalysis_2026/hgt_inference/au_test_tasks.tsv",
    "results/reanalysis_2026/primary_trees/PF00144.primary.rooted.treefile",
    "results/reanalysis_2026/primary_trees/PF00905.primary.rooted.treefile",
    "results/reanalysis_2026/primary_trees/PF13354.primary.rooted.treefile",
    "results/reanalysis_2026/primary_trees/PF00768.primary.rooted.treefile",
    "results/reanalysis_2026/primary_trees/PF13354_PF00768.primary.rooted.treefile",
    "figures/phylogenetic_trees/PF00144_primary_tree.pdf",
    "figures/phylogenetic_trees/PF00905_primary_tree.pdf",
    "figures/phylogenetic_trees/PF13354_primary_tree.pdf",
    "figures/phylogenetic_trees/PF00768_primary_tree.pdf",
    "figures/phylogenetic_trees/PF13354_PF00768_joint_tree.pdf",
    "figures/phylogenetic_trees/Overview_4Panel_trees.pdf",
]

REQUIRED_RESULTS = [
    "results/tables/Table_S01_HMM_Profile_Summary.csv",
    "results/tables/Table_S02_Filtering_Threshold_Sensitivity.csv",
    "results/tables/Table_S03_Bacterial_Reference_Prophage_Audit.csv",
    "results/tables/Table_S04_Viral_Origin_Confidence_Tiers.csv",
    "results/tables/Table_S05_HGT_Summary_By_Criteria.csv",
    "results/tables/Table_S06_HGT_Prophage_Pruning_Sensitivity.csv",
    "results/tables/Table_S07_aBSREL_Background_Normalized_Frequencies.csv",
    "results/tables/Table_S08_Lineage_Resolved_Selection_Trajectories.csv",
    "results/tables/Table_S09_MEME_Divergence_Alignment_Recombination_Diagnostics.csv",
    "results/tables/Table_S10_Catalytic_Motif_Purifying_Selection_Proof.csv",
    "results/tables/Table_S11_Structural_Assessment_Selected_Sites.csv",
    "results/tables/Table_S12_Stratified_Structural_Diversity.csv",
    "results/tables/Table_S13_Stratified_Habitat_Motif_Enrichment.csv",
    "results/supplementary_data/PerSequence_HMM_Scores.tsv",
    "results/supplementary_data/PerSequence_FIMO_Sensitivity.tsv",
    "results/supplementary_data/PerSequence_Bacterial_Audit.tsv",
    "results/supplementary_data/Supplementary_Data_Viral_Genomic_Neighborhoods.tsv",
    "results/supplementary_data/Supplementary_Data_HGT_Events_Catalog.tsv",
    "results/supplementary_data/PerSite_Selection_Metrics.tsv",
    "results/supplementary_data/PerSite_Selection_and_3D_Biophysics.tsv",
    "figures/Fig_S01_Viral_Genomic_Context_Verification.pdf",
    "figures/Fig_S01_Viral_Genomic_Context_Verification.png",
    "figures/Fig_S02_Lineage_Selection_Trajectories.pdf",
    "figures/Fig_S02_Lineage_Selection_Trajectories.png",
    "figures/Fig_S03_MEME_Diagnostics_and_Motif_Purifying.pdf",
    "figures/Fig_S03_MEME_Diagnostics_and_Motif_Purifying.png",
    "figures/Fig_S04_Structural_RSA_Distance_Selection.pdf",
    "figures/Fig_S04_Structural_RSA_Distance_Selection.png",
    "figures/Fig_S05_Confounder_Controlled_Structure_and_Motifs.pdf",
    "figures/Fig_S05_Confounder_Controlled_Structure_and_Motifs.png",
]


def verify_repo() -> int:
    """Verify repository structure, file non-emptiness, and GitHub <100 MB compatibility (git-tracked files)."""
    required = REQUIRED_DATA + SCRIPTS_ALL + SCRIPTS_HPC + [
        "scripts/10_plot_phylogenetic_trees.R",
        "scripts/run_hmm_motif_pipeline.sh",
        "scripts/sbl_config.py",
        "README.md",
        "environment.yml",
        "envs/hpc.yml",
        ".gitignore",
        "config/local.env.example",
        "config/hpc.env.example",
        "docs/METHODS_AND_PIPELINE.md",
        "data/selection/pf00144_absrel.json",
        "data/selection/pf00905_absrel.json",
        "data/selection/pf13354_absrel.json",
        "data/selection/pf00144_meme.json",
        "data/selection/pf00905_meme.json",
        "data/selection/pf13354_meme.json",
        "data/selection/PF00144_rename.faa",
        "data/selection/PF00905_rename.faa",
        "data/selection/PF13354_rename.faa",
    ] + REQUIRED_RESULTS + REQUIRED_REANALYSIS_2026

    missing = []
    empty = []
    for rel in required:
        p = REPO_ROOT / rel
        if not p.exists():
            missing.append(rel)
        elif p.stat().st_size == 0:
            empty.append(rel)

    external_missing = [rel for rel in OPTIONAL_EXTERNAL if not (REPO_ROOT / rel).exists()]

    skip_parts = {".git", ".wsl_stage", "private", "ProteinStructure"}
    oversized = []
    total_bytes = 0
    file_count = 0
    for p in REPO_ROOT.rglob("*"):
        if p.is_file() and not skip_parts.intersection(p.parts):
            sz = p.stat().st_size
            total_bytes += sz
            file_count += 1
            if sz >= 100 * 1024 * 1024:
                oversized.append((str(p.relative_to(REPO_ROOT)), round(sz / (1024 * 1024), 2)))

    print(f"Repository root : {REPO_ROOT}")
    print(f"Checked items   : {len(required)}")
    print(f"Total repo files: {file_count} ({total_bytes / (1024 * 1024):.2f} MB)")
    print(f"Missing items   : {len(missing)}")
    for m in missing:
        print(f"  [MISSING] {m}")
    print(f"Empty items     : {len(empty)}")
    for e in empty:
        print(f"  [EMPTY]   {e}")
    print(f"External (Figshare) inputs missing: {len(external_missing)}")
    for rel in external_missing:
        print(f"  [OPTIONAL] {rel}  (needed by scripts 08/09; see README > Data availability)")
    print(f"Files >= 100 MB : {len(oversized)}")
    for rel, mb in oversized:
        print(f"  [>100MB]  {rel} ({mb} MB)")

    if missing or empty or oversized:
        print("VERIFICATION FAILED.")
        return 1
    print("VERIFICATION PASSED: All data, scripts, docs, tables, supplementary TSVs, and figures verified.")
    return 0


def run_scripts(scripts: list[str]) -> int:
    """Execute pipeline scripts sequentially after checking the Linux tool environment."""
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    import sbl_config

    problems = sbl_config.preflight()
    if problems:
        print("PREFLIGHT FAILED:")
        for msg in problems:
            print(f"  - {msg}")
        return 1
    for rel in scripts:
        p = REPO_ROOT / rel
        print(f"\n>>> Running {rel} ...")
        res = subprocess.run([sys.executable, str(p)], cwd=str(REPO_ROOT))
        if res.returncode != 0:
            print(f"ERROR: {rel} exited with code {res.returncode}")
            return res.returncode
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="SBL_research reproducible bioinformatics pipeline")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--verify", action="store_true", help="Verify all files, outputs, and GitHub limits")
    group.add_argument("--run-fast", action="store_true", help="Run downstream analysis scripts (03-09)")
    group.add_argument("--run-all", action="store_true", help="Run complete analysis pipeline (01-09)")
    args = parser.parse_args()

    if args.verify:
        return verify_repo()
    if args.run_fast:
        rc = run_scripts(SCRIPTS_ALL[2:])
        return rc if rc != 0 else verify_repo()
    if args.run_all:
        rc = run_scripts(SCRIPTS_ALL)
        return rc if rc != 0 else verify_repo()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
