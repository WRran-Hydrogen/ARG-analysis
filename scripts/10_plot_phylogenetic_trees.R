#!/usr/bin/env Rscript
# ==============================================================================
# Circular phylogenetic tree figures for the primary IQ-TREE trees (PDF + LZW-TIFF + PNG)
#
# Usage:
#   1) Command line: Rscript scripts/10_plot_phylogenetic_trees.R
#   2) RStudio: open this script and Source it (the repository root is detected automatically)
#
# Rooting: ROOTING = "hgt" (default) draws the rooted trees written by the HGT analysis
#   ({fam}.primary.rooted.treefile: crown projection or MAD, see Rooting_Decisions_And_Agreement.tsv),
#   so figures and HGT tables share one root. "outgroup" roots on the 24-sequence Sister-PRE outgroup
#   (largest monophyletic subset) and "midpoint" gives the midpoint view.
# Branch lengths beyond the CAP_DEPTH_QUANTILE tip depth are shortened for display only; the
#   subtitle states when this is applied.
#
# Inputs:
#   - results/reanalysis_2026/primary_trees/*.treefile, *.primary.rooted.treefile, *.iqtree & *_outgroup.txt
#   - results/reanalysis_2026/hgt_inference/Rooting_Decisions_And_Agreement.tsv
#   - results/reanalysis_2026/Final_Viral_SBL_Catalog.tsv
#   - results/reanalysis_2026/Final_Bacterial_TripleSource_Catalog.tsv
#   - results/reanalysis_2026/Phase2_Viral_Genomic_Neighborhoods_Full.tsv
#
# Outputs:
#   - results/reanalysis_2026/primary_trees/Tree_Tip_Annotations_Merged.tsv
#   - figures/phylogenetic_trees/{PF00144,PF00905,PF13354,PF00768,PF13354_PF00768,Overview_4Panel}.{pdf,tiff,png}
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(cowplot)
  library(ragg)
})

# ==============================================================================
# SECTION 1: Adjustable parameters and colour palettes (tune the figures here only)
# ==============================================================================

# 1.1 Locate the repository root
locate_repo_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
    cand <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
    if (dir.exists(file.path(cand, "results", "reanalysis_2026"))) return(cand)
  }
  for (cand in c(getwd(), normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE))) {
    if (dir.exists(file.path(cand, "results", "reanalysis_2026"))) return(normalizePath(cand, winslash = "/"))
  }
  stop("Cannot locate the SBL_research root; run from the repository root or from scripts/.")
}

REPO_ROOT <- locate_repo_root()
TREE_DIR  <- file.path(REPO_ROOT, "results", "reanalysis_2026", "primary_trees")
META_DIR  <- file.path(REPO_ROOT, "results", "reanalysis_2026")
OUT_DIR   <- file.path(REPO_ROOT, "figures", "phylogenetic_trees")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 1.2 Layout and resolution
ROOTING            <- "hgt"       # "hgt" (same root as the HGT analysis), "outgroup" or "midpoint"
OPEN_ANGLE         <- 16      # fan opening angle (degrees); leaves room for the R1/R2/R3 ring labels
CAP_DEPTH_QUANTILE <- 0.985   # shorten the longest 1.5% of tips for display (1.0 = no capping)
SHOW_BOOT_NODES    <- TRUE    # mark well-supported nodes (UFBoot >= 95 & SH-aLRT >= 80)
SHOW_SBL_LABELS    <- TRUE    # label representative viral CARD-homolog SBLs on the outer ring
LABEL_ENV_SBL      <- FALSE   # also print "Motif-SBL" labels (default: purple diamond marker only)
MIN_LABEL_TIP_GAP  <- 18      # minimum number of tips between outer-ring labels (avoids overlap)
FIG_WIDTH          <- 13.0    # single-tree figure width (inches)
FIG_HEIGHT         <- 10.2    # single-tree figure height (inches)
TIFF_DPI           <- 300     # 600 for final production files; 300 dpi balances sharpness and file size
EXPORT_PNG_PREVIEW <- TRUE    # also write a PNG preview

# 1.3 Branch colours
PAL_BRANCH <- c(
  "Viral clade (>=80% V)"         = "#D62728",
  "Hydrolytic SBL clade (>=80%)"  = "#0D9488",
  "Bacterial / PBP backbone"      = "#94A3B8"
)
LW_BRANCH <- c(
  "Viral clade (>=80% V)"         = 0.42,
  "Hydrolytic SBL clade (>=80%)"  = 0.36,
  "Bacterial / PBP backbone"      = 0.18
)

# 1.4 Ring 1: origin and viral-origin tier
PAL_RING1 <- c(
  "Viral Tier 1 (High-conf.)"      = "#DC2626",
  "Viral Tier 2 (Supported)"       = "#F97316",
  "Viral Tier 3 (Fragment)"        = "#FACC15",
  "Viral MOSAiC (Arctic vMAG)"     = "#9333EA",
  "Bacterial Prophage"             = "#0284C7",
  "Bacterial Reference"            = "#E2E8F0"
)

# 1.5 Ring 2: catalytic function (CARD-homolog SBL vs motif-defined SBL vs PBP/PRE)
# CARD homolog = DIAMOND vs CARD (PBP entries removed) >= 40 % identity & >= 70 % coverage;
# motif-defined = class A (SXXK + SDN + E166xxLN + KTG) or class C (SXXK + YAN + KTG) without a
# CARD hit (class D always requires a CARD hit, see 09_build_curated_datasets.py).
LBL_V_CARD  <- "Viral CARD-homolog SBL"
LBL_V_MOTIF <- "Viral motif-defined SBL (class A/C)"
PAL_RING2 <- c(
  "Viral CARD-homolog SBL"                  = "#DB2777",
  "Viral motif-defined SBL (class A/C)"     = "#7E22CE",
  "Bacterial CARD-homolog SBL"              = "#059669",
  "Bacterial motif-defined SBL (class A/C)" = "#34D399",
  "PBP / PRE Homolog"                       = "#F1F5F9"
)

# 1.6 Ring 3: domain architecture (single-family trees) or Pfam family (joint tree)
PAL_RING3_ARCH <- c(
  "Viral Phage_spike/tail_CBM"       = "#2563EB",
  "Viral SLH+PF00144"                = "#7C3AED",
  "Transgly+PF00905+PASTA/C"         = "#2563EB",
  "Transgly+PF00905"                 = "#38BDF8",
  "PBP_dimer+PF00905"                = "#CBD5E1",
  "Viral TMP_3/Endolysin+PF13354"    = "#2563EB",
  "Viral/Bac Double-S11/PBP5_C"      = "#2563EB",
  "Other multi-domain"               = "#94A3B8",
  "Single-domain"                    = "#F8FAFC"
)

PAL_RING3_JOINT <- c(
  "PF13354 (Beta-lactamase2)" = "#4F46E5",
  "PF00768 (Peptidase_S11)"   = "#06B6D4"
)

# ==============================================================================
# SECTION 2: Read and merge the tip annotation tables
# ==============================================================================

clean_card_family <- function(card_family, card_pident, card_cov, func_tier) {
  ok <- !is.na(card_pident) & !is.na(card_cov) & card_pident >= 40 & card_cov >= 70
  out <- rep("NONE", length(card_family))
  if (any(ok)) {
    lbl <- card_family[ok]
    lbl <- sub(";.*$", "", lbl)
    lbl <- gsub("\\s*beta-lactamases?.*$", "", lbl, ignore.case = TRUE)
    lbl <- gsub("\\s*pending classification.*$", "", lbl, ignore.case = TRUE)
    lbl <- gsub("class A Bacillus cereus Bc", "Bc-like", lbl, ignore.case = TRUE)
    lbl <- gsub("-like$", "", lbl)
    lbl <- trimws(lbl)
    lbl[lbl == "" | lbl == "NONE"] <- "CARD-SBL"
    out[ok] <- lbl
  }
  env_idx <- !ok & func_tier == "Putative_Hydrolytic_BetaLactamase"
  out[env_idx] <- "Motif-SBL"
  out
}

classify_arch_group <- function(arch_str, origin, fam) {
  arch_clean <- gsub("\\(\\d+-\\d+\\)", "", ifelse(is.na(arch_str), "", arch_str))
  is_multi <- grepl("\\+", arch_clean)
  if (fam == "PF00144") {
    if (origin == "Viral" && grepl("Phage_spike_CBM|Phage_tail_CBM", arch_clean)) return("Viral Phage_spike/tail_CBM")
    if (origin == "Viral" && grepl("\\bSLH\\b", arch_clean)) return("Viral SLH+PF00144")
    return(ifelse(is_multi, "Other multi-domain", "Single-domain"))
  }
  if (fam == "PF00905") {
    if (grepl("Transgly", arch_clean) && grepl("PASTA|PBP_1A_C|PBP_middle|PCB_OB", arch_clean)) return("Transgly+PF00905+PASTA/C")
    if (grepl("Transgly", arch_clean)) return("Transgly+PF00905")
    if (grepl("PBP_dimer", arch_clean)) return("PBP_dimer+PF00905")
    return(ifelse(is_multi, "Other multi-domain", "Single-domain"))
  }
  if (fam == "PF13354") {
    if (origin == "Viral" && grepl("TMP_3|Endolysin|SLT|Peptidase_M23|VanY", arch_clean)) return("Viral TMP_3/Endolysin+PF13354")
    return(ifelse(is_multi, "Other multi-domain", "Single-domain"))
  }
  if (fam == "PF00768") {
    if (grepl("PBP5_C|Peptidase_S11\\+Peptidase_S11", arch_clean)) return("Viral/Bac Double-S11/PBP5_C")
    return(ifelse(is_multi, "Other multi-domain", "Single-domain"))
  }
  ifelse(is_multi, "Other multi-domain", "Single-domain")
}

as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

# Ring 3 of the joint tree only: PF13354 vs PF00768; other families get NA
joint_family_label <- function(acc) {
  case_when(acc == "PF13354" ~ "PF13354 (Beta-lactamase2)",
            acc == "PF00768" ~ "PF00768 (Peptidase_S11)",
            TRUE ~ NA_character_)
}

build_merged_metadata <- function() {
  message("[Metadata] Loading curated viral, bacterial, and Phase 2 neighborhood catalogs...")
  v_raw <- read_tsv(file.path(META_DIR, "Final_Viral_SBL_Catalog.tsv"), show_col_types = FALSE)
  b_raw <- read_tsv(file.path(META_DIR, "Final_Bacterial_TripleSource_Catalog.tsv"), show_col_types = FALSE)
  p2_raw <- read_tsv(
    file.path(META_DIR, "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"),
    col_select = c(
      "uid", "contig_length_bp", "genomad_virus_score", "checkv_viral_genes",
      "checkv_host_genes", "n_viral_hallmark_in_region",
      "flanked_strict_structural_hallmark_both", "upstream_5_orfs", "downstream_5_orfs"
    ),
    show_col_types = FALSE
  )

  v_mge <- p2_raw |>
    mutate(
      nb = paste(coalesce(upstream_5_orfs, ""), coalesce(downstream_5_orfs, ""), sep = " ; "),
      has_plasmid_rep = grepl("IncF|repA|Rep_1|Rep_2|Rep_3|Rep_trans|StbA|ParA|ParB", nb, ignore.case = TRUE),
      has_mob_relax   = grepl("Relaxase|Mob_Pre|MobC|MobA|MbeA|TrwC|oriT|Conjugal", nb, ignore.case = TRUE),
      viral_mge_context = case_when(
        has_plasmid_rep ~ "Phage_Plasmid_Hybrid",
        has_mob_relax   ~ "Phage_Mobilizable_Cassette",
        TRUE            ~ "Canonical_Phage"
      )
    ) |>
    select(uid, contig_length_bp, checkv_viral_genes, checkv_host_genes,
           n_viral_hallmark_in_region, flanked_strict_structural_hallmark_both, viral_mge_context)

  v_df <- v_raw |>
    filter(as_bool(pass_cov65), as_bool(pass_fimo_1e5)) |>
    left_join(v_mge, by = "uid") |>
    mutate(
      origin = "Viral",
      card_verified = (as.numeric(card_pident) >= 40 & as.numeric(card_cov) >= 70),
      ring1_origin = case_when(
        source == "MOSAiC" | grepl("MOSAiC", viral_origin_tier) ~ "Viral MOSAiC (Arctic vMAG)",
        viral_origin_tier == "Tier1_HighConfidence_Viral"       ~ "Viral Tier 1 (High-conf.)",
        viral_origin_tier == "Tier2_Supported_Viral"            ~ "Viral Tier 2 (Supported)",
        TRUE                                                    ~ "Viral Tier 3 (Fragment)"
      ),
      ring2_function = case_when(
        functional_tier == "Putative_Hydrolytic_BetaLactamase" & card_verified ~ LBL_V_CARD,
        functional_tier == "Putative_Hydrolytic_BetaLactamase"                 ~ LBL_V_MOTIF,
        TRUE                                                                   ~ "PBP / PRE Homolog"
      ),
      ring3_joint_fam = joint_family_label(winning_acc),
      sbl_display_label = clean_card_family(card_gene_family, as.numeric(card_pident), as.numeric(card_cov), functional_tier),
      taxon_group = coalesce(viral_mge_context, "MOSAiC_vMAG")
    )
  v_df$ring3_arch <- mapply(classify_arch_group, v_df$domain_architecture, "Viral", v_df$winning_acc, USE.NAMES = FALSE)

  b_df <- b_raw |>
    filter(as_bool(pass_cov65), as_bool(pass_fimo_1e5), !as_bool(is_viral_contig_in_genome)) |>
    mutate(
      origin = "Bacterial",
      card_verified = (as.numeric(card_pident) >= 40 & as.numeric(card_cov) >= 70),
      ring1_origin = ifelse(as_bool(in_bacterial_prophage), "Bacterial Prophage", "Bacterial Reference"),
      ring2_function = case_when(
        functional_tier == "Putative_Hydrolytic_BetaLactamase" & card_verified ~ "Bacterial CARD-homolog SBL",
        functional_tier == "Putative_Hydrolytic_BetaLactamase"                 ~ "Bacterial motif-defined SBL (class A/C)",
        TRUE                                                                   ~ "PBP / PRE Homolog"
      ),
      ring3_joint_fam = joint_family_label(winning_acc),
      sbl_display_label = clean_card_family(card_gene_family, as.numeric(card_pident), as.numeric(card_cov), functional_tier),
      taxon_group = phylum,
      viral_mge_context = NA_character_,
      viral_origin_tier = NA_character_,
      viral_cluster_size = NA_integer_,
      contig_length_bp = NA_real_,
      checkv_viral_genes = NA_real_,
      checkv_host_genes = NA_real_,
      n_viral_hallmark_in_region = NA_real_,
      flanked_strict_structural_hallmark_both = NA
    )
  b_df$ring3_arch <- mapply(classify_arch_group, b_df$domain_architecture, "Bacterial", b_df$winning_acc, USE.NAMES = FALSE)

  keep_cols <- c(
    "uid", "origin", "winning_acc", "functional_tier", "functional_evidence",
    "card_verified", "card_gene_family", "card_pident", "card_cov", "sbl_display_label",
    "ring1_origin", "ring2_function", "ring3_arch", "ring3_joint_fam",
    "taxon_group", "viral_origin_tier", "viral_cluster_size", "viral_mge_context",
    "domain_architecture", "hmm_cov", "motif1_SXXK_hmm", "motif2_class", "motif3_KTG_hmm", "has_omega_E166"
  )
  merged <- bind_rows(select(v_df, all_of(keep_cols)), select(b_df, all_of(keep_cols)))
  # guard against duplicated uids
  merged <- distinct(merged, uid, .keep_all = TRUE)

  out_tsv <- file.path(TREE_DIR, "Tree_Tip_Annotations_Merged.tsv")
  write_tsv(merged, out_tsv)
  message("  -> Saved unified tip metadata table: ", out_tsv, " (n=", nrow(merged), ")")
  merged
}

# ==============================================================================
# SECTION 3: Rooting, display-only branch capping and clade annotation
# ==============================================================================

prepare_tree_and_nodes <- function(tree_path, outgroup_path, meta_df, rooted_path = NULL, rooted_note = "") {
  use_hgt_root <- ROOTING == "hgt" && !is.null(rooted_path) && file.exists(rooted_path)
  tr <- ape::read.tree(if (use_hgt_root) rooted_path else tree_path)
  og <- read_lines(outgroup_path)
  og <- trimws(og[nzchar(trimws(og))])

  drop_ids <- intersect(tr$tip.label, og)
  root_note <- "midpoint-rooted"
  if (use_hgt_root) {
    if (!ape::is.rooted(tr)) stop("Expected a rooted tree: ", rooted_path)
    root_note <- rooted_note
  } else if (ROOTING %in% c("outgroup", "hgt") && length(drop_ids) > 0) {
    rooted <- root_on_outgroup(tr, drop_ids)
    tr <- rooted$tree
    root_note <- paste0(rooted$note, "; outgroup tips removed")
  }
  if (length(drop_ids) > 0) {
    tr <- ape::drop.tip(tr, drop_ids)
  }
  if (!use_hgt_root && (ROOTING == "midpoint" || length(drop_ids) == 0)) {
    tr <- phangorn::midpoint(tr)
  }
  tr <- ape::ladderize(tr, right = FALSE)

  # Shorten the longest 1.5% of tips (display only; parent depth <= child depth keeps branch lengths non-negative)
  n_capped <- 0L
  if (CAP_DEPTH_QUANTILE < 1.0 && !is.null(tr$edge.length)) {
    depths <- ape::node.depth.edgelength(tr)
    tip_depths <- depths[seq_len(length(tr$tip.label))]
    q_cap <- as.numeric(quantile(tip_depths, probs = CAP_DEPTH_QUANTILE, na.rm = TRUE))
    if (is.finite(q_cap) && q_cap > 0) {
      n_capped <- sum(tip_depths > q_cap)
      d_cap <- pmin(depths, q_cap)
      tr$edge.length <- pmax(0, d_cap[tr$edge[, 2]] - d_cap[tr$edge[, 1]])
    }
  }

  # Parse IQ-TREE SH-aLRT / UFBoot support values
  n_tips <- length(tr$tip.label)
  n_nodes <- tr$Nnode
  tot_nodes <- n_tips + n_nodes

  sh_val <- rep(NA_real_, tot_nodes)
  uf_val <- rep(NA_real_, tot_nodes)
  if (!is.null(tr$node.label) && length(tr$node.label) == n_nodes) {
    lbls <- tr$node.label
    has_slash <- grepl("/", lbls)
    if (any(has_slash)) {
      parts <- str_split_fixed(lbls[has_slash], "/", 2)
      idx_n <- n_tips + which(has_slash)
      sh_val[idx_n] <- suppressWarnings(as.numeric(parts[, 1]))
      uf_val[idx_n] <- suppressWarnings(as.numeric(parts[, 2]))
    }
  }

  # Bottom-up count of viral and hydrolytic SBL descendants per node
  tip_meta <- tibble(uid = tr$tip.label) |> left_join(meta_df, by = "uid")
  n_desc <- integer(tot_nodes)
  n_v    <- integer(tot_nodes)
  n_sbl  <- integer(tot_nodes)

  n_desc[seq_len(n_tips)] <- 1L
  n_v[seq_len(n_tips)]    <- ifelse(!is.na(tip_meta$origin) & tip_meta$origin == "Viral", 1L, 0L)
  n_sbl[seq_len(n_tips)]  <- ifelse(!is.na(tip_meta$functional_tier) & tip_meta$functional_tier == "Putative_Hydrolytic_BetaLactamase", 1L, 0L)

  tr_po <- ape::reorder.phylo(tr, "postorder")
  for (i in seq_len(nrow(tr_po$edge))) {
    p <- tr_po$edge[i, 1]
    c <- tr_po$edge[i, 2]
    n_desc[p] <- n_desc[p] + n_desc[c]
    n_v[p]    <- n_v[p] + n_v[c]
    n_sbl[p]  <- n_sbl[p] + n_sbl[c]
  }

  branch_group <- ifelse(
    (n_v / pmax(n_desc, 1)) >= 0.80,
    "Viral clade (>=80% V)",
    ifelse((n_sbl / pmax(n_desc, 1)) >= 0.80 & n_desc >= 5, "Hydrolytic SBL clade (>=80%)", "Bacterial / PBP backbone")
  )

  node_df <- tibble(
    node = seq_len(tot_nodes),
    n_desc = n_desc,
    n_viral_desc = n_v,
    n_sbl_desc = n_sbl,
    sh_alrt = sh_val,
    ufboot = uf_val,
    high_support = (seq_len(tot_nodes) > n_tips) & !is.na(uf_val) & !is.na(sh_val) & uf_val >= 95 & sh_val >= 80 & n_desc >= 5,
    branch_group = factor(branch_group, levels = names(PAL_BRANCH))
  )

  # tip_df keeps only the columns used by the outer rings, avoiding clashes with p$data / node_df
  tip_full <- tip_meta |>
    mutate(
      id = uid,
      ring1_origin = factor(ring1_origin, levels = names(PAL_RING1)),
      ring2_function = factor(ring2_function, levels = names(PAL_RING2)),
      ring3_arch = factor(ring3_arch, levels = names(PAL_RING3_ARCH)),
      ring3_joint_fam = factor(ring3_joint_fam, levels = names(PAL_RING3_JOINT))
    ) |>
    select(
      id, origin, functional_tier, card_verified, sbl_display_label,
      ring1_origin, ring2_function, ring3_arch, ring3_joint_fam
    )

  list(tree = tr, node_df = node_df, tip_df = tip_full, root_note = root_note, n_capped = n_capped)
}

# Root on the Sister-PRE outgroup. If the outgroup is not monophyletic in the unrooted tree,
# root on its largest monophyletic subset and report that in the subtitle.
root_on_outgroup <- function(tr, og) {
  res <- tryCatch(ape::root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE), error = function(e) NULL)
  if (!is.null(res) && ape::is.monophyletic(res, og)) {
    return(list(tree = res, note = sprintf("outgroup-rooted (%d/%d outgroup tips monophyletic)", length(og), length(og))))
  }
  ing <- setdiff(tr$tip.label, og)
  tmp <- ape::root(tr, outgroup = ing[1], resolve.root = TRUE, edgelabel = TRUE)
  n_t <- ape::Ntip(tmp)
  clades <- ape::prop.part(tmp)   # tip sets of every internal node, in node order
  best <- integer(0)
  for (cl in clades) {
    tips <- attr(clades, "labels")[cl]
    if (all(tips %in% og) && length(tips) > length(best)) best <- tips
  }
  if (length(best) == 0) best <- og[1]
  list(
    tree = ape::root(tr, outgroup = best, resolve.root = TRUE, edgelabel = TRUE),
    note = sprintf("outgroup-rooted on largest monophyletic subset (%d/%d outgroup tips)", length(best), length(og))
  )
}

# Root used by the HGT analysis for one family (Rooting_Decisions_And_Agreement.tsv)
hgt_root_note <- function(tag) {
  f <- file.path(META_DIR, "hgt_inference", "Rooting_Decisions_And_Agreement.tsv")
  if (!file.exists(f)) return("rooted as in the HGT analysis")
  d <- read_tsv(f, show_col_types = FALSE)
  r <- d[d$family == tag, ]
  if (nrow(r) == 0) return("rooted as in the HGT analysis")
  m <- r$primary_rooting[1]
  agree <- if ("root_agreement_primary_in_MAD" %in% names(r)) r$root_agreement_primary_in_MAD[1] else NA
  sprintf("rooted as in the HGT analysis (%s%s)", gsub("_", " ", m),
          if (!is.na(agree) && m != "MAD") sprintf("; root agreement with MAD = %.2f", agree) else "")
}

# Best-fit model from the IQ-TREE report next to the tree file (falls back to spec$model)
read_iqtree_model <- function(tree_path, fallback) {
  rep <- sub("\\.treefile$", ".iqtree", tree_path)
  if (!file.exists(rep)) return(paste(fallback, "(model not verified)"))
  hit <- grep("Best-fit model according to BIC:", readLines(rep, warn = FALSE), value = TRUE)
  if (length(hit) == 0) return(paste(fallback, "(model not verified)"))
  trimws(sub(".*BIC:", "", hit[1]))
}

# Pick outer-ring viral SBL labels at least MIN_LABEL_TIP_GAP tips apart (avoids overlap)
select_spaced_sbl_labels <- function(tip_coords_df, min_gap = MIN_LABEL_TIP_GAP) {
  cand <- tip_coords_df |>
    filter(origin == "Viral", functional_tier == "Putative_Hydrolytic_BetaLactamase")
  if (!LABEL_ENV_SBL) {
    cand <- cand |> filter(card_verified == TRUE)
  }
  cand <- cand |> arrange(y)
  if (nrow(cand) == 0) return(cand)

  keep_idx <- integer(0)
  last_y <- -Inf
  for (i in seq_len(nrow(cand))) {
    if ((cand$y[i] - last_y) >= min_gap) {
      keep_idx <- c(keep_idx, i)
      last_y <- cand$y[i]
    } else if (length(keep_idx) > 0) {
      prev_i <- keep_idx[length(keep_idx)]
      if (!cand$card_verified[prev_i] && cand$card_verified[i]) {
        keep_idx[length(keep_idx)] <- i
        last_y <- cand$y[i]
      }
    }
  }
  cand[keep_idx, , drop = FALSE]
}

# ==============================================================================
# SECTION 4: ggtree + ggtreeExtra plotting
# ==============================================================================

build_ggtree_plot <- function(spec, meta_df, compact_mode = FALSE) {
  tree_path <- file.path(TREE_DIR, spec$treefile)
  og_path   <- file.path(TREE_DIR, spec$outgroup)
  rooted_path <- file.path(TREE_DIR, paste0(spec$tag, ".primary.rooted.treefile"))
  prep <- prepare_tree_and_nodes(tree_path, og_path, meta_df, rooted_path, hgt_root_note(spec$tag))
  tr      <- prep$tree
  node_df <- prep$node_df
  tip_df  <- prep$tip_df

  n_tot  <- nrow(tip_df)
  n_v    <- sum(tip_df$origin == "Viral", na.rm = TRUE)
  n_b    <- n_tot - n_v
  n_vsbl_clin <- sum(tip_df$origin == "Viral" & tip_df$functional_tier == "Putative_Hydrolytic_BetaLactamase" & tip_df$card_verified, na.rm = TRUE)
  n_vsbl_env  <- sum(tip_df$origin == "Viral" & tip_df$functional_tier == "Putative_Hydrolytic_BetaLactamase" & !tip_df$card_verified, na.rm = TRUE)
  n_bsbl      <- sum(tip_df$origin == "Bacterial" & tip_df$functional_tier == "Putative_Hydrolytic_BetaLactamase", na.rm = TRUE)

  p <- ggtree(tr, layout = "fan", open.angle = OPEN_ANGLE, aes(color = branch_group, linewidth = branch_group)) %<+% node_df +
    scale_color_manual(name = "Branches", values = PAL_BRANCH, drop = FALSE) +
    scale_linewidth_manual(name = "Branches", values = if (compact_mode) LW_BRANCH * 0.72 else LW_BRANCH, drop = FALSE)

  if (SHOW_BOOT_NODES && !compact_mode) {
    p <- p + geom_point2(aes(subset = high_support), size = 0.75, color = "#334155", alpha = 0.55, show.legend = FALSE)
  }

  max_x <- max(p$data$x, na.rm = TRUE)
  ring_w <- 0.068 * max_x
  gap_w  <- 0.015 * max_x

  # Ring 1: Origin & Viral Tier
  p <- p +
    geom_fruit(
      data = tip_df,
      geom = geom_tile,
      mapping = aes(y = id, fill = ring1_origin),
      width = ring_w,
      offset = 0.04,
      color = NA
    ) +
    scale_fill_manual(
      name = "Ring 1: Origin & Viral Tier",
      values = PAL_RING1,
      drop = TRUE,
      na.translate = FALSE,
      guide = guide_legend(order = 2)
    )

  # Ring 2: 3-Tier Catalytic Function (Mechanistic + CARD)
  p <- p +
    new_scale_fill() +
    geom_fruit(
      data = tip_df,
      geom = geom_tile,
      mapping = aes(y = id, fill = ring2_function),
      width = ring_w,
      offset = 0.082,
      color = NA
    ) +
    scale_fill_manual(
      name = "Ring 2: Catalytic Function",
      values = PAL_RING2,
      drop = TRUE,
      na.translate = FALSE,
      guide = guide_legend(order = 3)
    )

  # Ring 3: domain architecture (single-family trees) or Pfam family (joint tree)
  is_joint <- identical(spec$tag, "PF13354_PF00768")
  if (is_joint) {
    p <- p +
      new_scale_fill() +
      geom_fruit(
        data = tip_df,
        geom = geom_tile,
        mapping = aes(y = id, fill = ring3_joint_fam),
        width = ring_w,
        offset = 0.082,
        color = NA
      ) +
      scale_fill_manual(
        name = "Ring 3: Pfam Family",
        values = PAL_RING3_JOINT,
        drop = TRUE,
        na.translate = FALSE,
        guide = guide_legend(order = 4)
      )
  } else {
    p <- p +
      new_scale_fill() +
      geom_fruit(
        data = tip_df,
        geom = geom_tile,
        mapping = aes(y = id, fill = ring3_arch),
        width = ring_w,
        offset = 0.082,
        color = NA
      ) +
      scale_fill_manual(
        name = "Ring 3: Domain Architecture",
        values = PAL_RING3_ARCH,
        drop = TRUE,
        na.translate = FALSE,
        guide = guide_legend(order = 4)
      )
  }

  # Outer markers: viral hydrolytic SBLs (CARD homolog vs motif-defined)
  tip_coords <- p$data |>
    filter(isTip == TRUE) |>
    select(id = label, y, angle) |>
    left_join(tip_df, by = "id")

  vsbl_pts <- tip_coords |>
    filter(origin == "Viral", functional_tier == "Putative_Hydrolytic_BetaLactamase") |>
    mutate(sbl_marker = ifelse(card_verified, LBL_V_CARD, LBL_V_MOTIF))

  r_marker <- max_x * 1.29
  if (nrow(vsbl_pts) > 0) {
    p <- p +
      new_scale_color() +
      geom_point(
        data = vsbl_pts,
        aes(x = r_marker, y = y, color = sbl_marker, shape = sbl_marker),
        size = if (compact_mode) 1.6 else 2.4,
        stroke = 0.45,
        inherit.aes = FALSE
      ) +
      scale_color_manual(
        name = "Outer Markers (Viral SBL)",
        values = setNames(c("#DB2777", "#7E22CE"), c(LBL_V_CARD, LBL_V_MOTIF)),
        drop = TRUE,
        guide = guide_legend(order = 5)
      ) +
      scale_shape_manual(
        name = "Outer Markers (Viral SBL)",
        values = setNames(c(8, 18), c(LBL_V_CARD, LBL_V_MOTIF)),   # asterisk / diamond
        drop = TRUE,
        guide = guide_legend(order = 5)
      )
  }

  # Sparse outer-ring gene labels (non-overlapping)
  if (SHOW_SBL_LABELS && !compact_mode) {
    lbl_df <- select_spaced_sbl_labels(tip_coords, min_gap = max(14, round(n_tot / 95)))
    if (nrow(lbl_df) > 0) {
      # In the fan layout angle is in [0, 360]; labels on the left half (90-270) are flipped by 180 degrees to stay readable
      lbl_df <- lbl_df |>
        mutate(
          ang = (angle %% 360),
          flip = ang > 90 & ang < 270,
          text_angle = ifelse(flip, ang - 180, ang),
          text_hjust = ifelse(flip, 1, 0),
          lbl_col = ifelse(card_verified, "#9D174D", "#581C87")
        )
      p <- p +
        geom_text(
          data = lbl_df,
          aes(x = max_x * 1.325, y = y, label = sbl_display_label, angle = text_angle, hjust = text_hjust),
          color = lbl_df$lbl_col,
          size = 2.55,
          fontface = "bold",
          inherit.aes = FALSE
        )
    }
  }

  # R1 / R2 / R3 ring labels inside the fan opening
  ring_centers <- max_x * c(1.04 + 0.034, 1.122 + 0.034, 1.204 + 0.034)
  r_lbl_df <- tibble(
    x = ring_centers,
    y = max(p$data$y, na.rm = TRUE) + max(8, n_tot * 0.012),
    lab = c("R1", "R2", "R3")
  )
  p <- p +
    geom_text(
      data = r_lbl_df,
      aes(x = x, y = y, label = lab),
      size = if (compact_mode) 2.2 else 2.8,
      fontface = "bold",
      color = "#334155",
      inherit.aes = FALSE
    )

  # Widen the radial limit so outer labels are not clipped
  x_lim_hi <- if (SHOW_SBL_LABELS && !compact_mode) max_x * 1.50 else max_x * 1.35
  sub_txt <- sprintf(
    "Ingroup n=%s (Viral nr85=%s, Bacterial=%s) | Viral SBL: CARD-homolog=%d, motif-defined=%d | Bac SBL=%d | %s",
    format(n_tot, big.mark = ","), format(n_v, big.mark = ","), format(n_b, big.mark = ","),
    n_vsbl_clin, n_vsbl_env, n_bsbl,
    paste(read_iqtree_model(tree_path, spec$model), spec$track, sep = " | ")
  )
  sub_txt <- paste0(
    sub_txt, "
", prep$root_note,
    if (prep$n_capped > 0) sprintf("; %d longest tips shortened to the %.1fth percentile depth (display only)",
                                   prep$n_capped, 100 * CAP_DEPTH_QUANTILE) else ""
  )

  p <- p +
    xlim(0, x_lim_hi) +
    labs(
      title = if (compact_mode) spec$short else spec$title,
      subtitle = if (compact_mode) sprintf("n=%s (V=%s, B=%s) | Viral SBL: %d CARD + %d motif", format(n_tot, big.mark=","), n_v, n_b, n_vsbl_clin, n_vsbl_env) else sub_txt
    ) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title       = element_text(size = if (compact_mode) 11 else 14, face = "bold", color = "#0F172A", hjust = 0),
      plot.subtitle    = element_text(size = if (compact_mode) 8.5 else 9.5, color = "#475569", hjust = 0),
      legend.position  = if (compact_mode) "none" else "right",
      legend.title     = element_text(size = 9.0, face = "bold", color = "#0F172A"),
      legend.text      = element_text(size = 8.0, color = "#334155"),
      legend.key.size  = unit(0.42, "cm"),
      legend.spacing.y = unit(0.18, "cm"),
      plot.margin      = margin(8, 8, 8, 8)
    )

  p
}

# ==============================================================================
# SECTION 5: Export PDF (vector) + LZW-TIFF (raster) + PNG preview
# ==============================================================================

save_multi_format <- function(plot_obj, stem, width, height, dpi = TIFF_DPI) {
  pdf_path  <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  tiff_path <- file.path(OUT_DIR, paste0(stem, ".tiff"))
  png_path  <- file.path(OUT_DIR, paste0(stem, ".png"))

  # PDF (cairo_pdf when available, otherwise the standard pdf device)
  tryCatch(
    ggsave(pdf_path, plot = plot_obj, device = grDevices::cairo_pdf, width = width, height = height, units = "in"),
    error = function(e) ggsave(pdf_path, plot = plot_obj, device = "pdf", width = width, height = height, units = "in")
  )

  # LZW-compressed TIFF
  ragg::agg_tiff(tiff_path, width = width, height = height, units = "in", res = dpi, compression = "lzw", background = "white")
  print(plot_obj)
  dev.off()

  if (EXPORT_PNG_PREVIEW) {
    ragg::agg_png(png_path, width = width, height = height, units = "in", res = 300, background = "white")
    print(plot_obj)
    dev.off()
  }
  message("  [Exported] ", stem, " -> .pdf / .tiff", ifelse(EXPORT_PNG_PREVIEW, " / .png", ""))
}

main <- function() {
  meta_df <- build_merged_metadata()

  specs <- list(
    list(tag = "PF00144",         stem = "PF00144_primary_tree",         title = "PF00144 (Beta-lactamase / Class C-like)",                  short = "A. PF00144 (Class C-like)",          treefile = "PF00144.primary.hmm_out.treefile",         outgroup = "PF00144_outgroup.txt",         model = "Q.pfam+I+R10", track = "hmmalign"),
    list(tag = "PF00905",         stem = "PF00905_primary_tree",         title = "PF00905 (Transpeptidase / HMW-PBP & Class D-like)",        short = "B. PF00905 (Transpeptidase)",        treefile = "PF00905.primary.hmm_out.treefile",         outgroup = "PF00905_outgroup.txt",         model = "Q.pfam+R10", track = "hmmalign"),
    list(tag = "PF13354",         stem = "PF13354_primary_tree",         title = "PF13354 (Beta-lactamase2 / Class A-like)",                 short = "C. PF13354 (Class A-like)",          treefile = "PF13354.primary.hmm_out.treefile",         outgroup = "PF13354_outgroup.txt",         model = "Q.pfam+R10", track = "hmmalign"),
    list(tag = "PF00768",         stem = "PF00768_primary_tree",         title = "PF00768 (Peptidase_S11 / Sister LMW-PBP)",                 short = "PF00768 (Peptidase_S11)",            treefile = "PF00768.primary.hmm_out.treefile",         outgroup = "PF00768_outgroup.txt",         model = "Q.pfam+I+R10", track = "hmmalign"),
    list(tag = "PF13354_PF00768", stem = "PF13354_PF00768_joint_tree",   title = "Joint PF13354 + PF00768 (Beta-lactamase2 vs. Peptidase_S11)", short = "D. Joint PF13354 + PF00768",     treefile = "PF13354_PF00768.primary.mafft_out.treefile", outgroup = "PF13354_PF00768_outgroup.txt", model = "Q.pfam+R10", track = "MAFFT+ClipKIT")
  )

  for (sp in specs) {
    message("[Plotting] ", sp$stem, " ...")
    p <- build_ggtree_plot(sp, meta_df, compact_mode = FALSE)
    save_multi_format(p, sp$stem, width = FIG_WIDTH, height = FIG_HEIGHT)
  }

  # 2x2 overview panel (PF00144, PF00905, PF13354, PF13354_PF00768)
  message("[Plotting] Overview_4Panel_trees ...")
  ov_specs <- specs[c(1, 2, 3, 5)]
  ov_plots <- lapply(ov_specs, function(sp) build_ggtree_plot(sp, meta_df, compact_mode = TRUE))
  p_for_leg <- build_ggtree_plot(specs[[1]], meta_df, compact_mode = FALSE) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal")
  shared_leg <- cowplot::get_legend(p_for_leg)

  # Extra legend for Ring 3 of panel D (PF13354 vs PF00768)
  d_leg_df <- tibble(
    x = c(1, 2),
    y = c(1, 1),
    fam = factor(names(PAL_RING3_JOINT), levels = names(PAL_RING3_JOINT))
  )
  p_joint_leg <- ggplot(d_leg_df, aes(x, y, fill = fam)) +
    geom_tile() +
    scale_fill_manual(name = "Ring 3 (Panel D): Joint Family", values = PAL_RING3_JOINT) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 9.0, face = "bold", color = "#0F172A"),
      legend.text = element_text(size = 8.0, color = "#334155"),
      legend.key.size = unit(0.42, "cm")
    )
  joint_leg <- cowplot::get_legend(p_joint_leg)

  combined_leg <- cowplot::plot_grid(shared_leg, joint_leg, ncol = 1, rel_heights = c(1, 0.22))
  grid_2x2 <- cowplot::plot_grid(plotlist = ov_plots, ncol = 2, align = "hv")
  ov_full  <- cowplot::plot_grid(grid_2x2, combined_leg, ncol = 1, rel_heights = c(1, 0.20))

  save_multi_format(ov_full, "Overview_4Panel_trees", width = 14.5, height = 15.8)
  message("[Done] All R tree figures saved to: ", OUT_DIR)
}

main()
