#!/usr/bin/env Rscript
# ==============================================================================
# Circular phylogenetic trees (domain-based functional annotation) and rectangular
# domain-context supplementary figures. Trees are drawn with the root used by the HGT
# analysis ({fam}.primary.rooted.treefile; see 04_hgt_inference/Rooting_Decisions_And_Agreement.tsv).
#
# Usage:
#   1) Command line: Rscript scripts/visualization/plot_phylogenetic_trees_2026.R
#   2) RStudio: open this script and Source it (Ctrl+Alt+R); the repository root is detected automatically
#
# Inputs:
#   - results/reanalysis_2026/03_primary_trees/*.treefile & *_outgroup.txt
#   - results/reanalysis_2026/01_curation_and_context/Final_Viral_SBL_Catalog.tsv
#   - results/reanalysis_2026/01_curation_and_context/Final_Bacterial_TripleSource_Catalog.tsv
#   - results/reanalysis_2026/01_curation_and_context/Phase2_Viral_Genomic_Neighborhoods_Full.tsv
#
# Outputs:
#   - results/reanalysis_2026/03_primary_trees/Tree_Tip_Annotations_Merged.tsv
#   - figures/trees_2026_R/{PF00144_primary_tree,PF00905_primary_tree,PF13354_primary_tree,PF00768_primary_tree,PF13354_PF00768_joint_tree,Fig_Main_PF00144_PF00905,FigS_PF13354_PF00768_Joint}.{pdf,tiff,png}
#   - figures/trees_2026_R/FigS_DomainContext_{PF00144,PF00905,PF13354,PF00768}.{pdf,tiff,png}
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggtree)
  library(ggnewscale)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(patchwork)
  library(cowplot)
  library(ragg)
})

# ==============================================================================
# SECTION 1: Adjustable parameters and colour palettes (edit here when fine-tuning figures)
# ==============================================================================

# 1.1 Locate the repository root
locate_repo_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
    for (rel in c("..", file.path("..", ".."))) {
      cand <- normalizePath(file.path(dirname(script_path), rel), winslash = "/", mustWork = FALSE)
      if (dir.exists(file.path(cand, "results", "reanalysis_2026"))) return(cand)
    }
  }
  for (cand in c(
    getwd(),
    normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = FALSE)
  )) {
    if (dir.exists(file.path(cand, "results", "reanalysis_2026"))) return(normalizePath(cand, winslash = "/"))
  }
  stop("Cannot locate the SBL_research root; run from the repository root or scripts/visualization/.")
}

REPO_ROOT <- locate_repo_root()
TREE_DIR  <- file.path(REPO_ROOT, "results", "reanalysis_2026", "03_primary_trees")
META_DIR  <- file.path(REPO_ROOT, "results", "reanalysis_2026", "01_curation_and_context")
OUT_DIR   <- file.path(REPO_ROOT, "figures", "trees_2026_R")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
HGT_DIR   <- file.path(REPO_ROOT, "results", "reanalysis_2026", "04_hgt_inference")

# Subtitle "model | alignment | root" from the HGT rooting decisions and the IQ-TREE report
tree_method_label <- function(tag) {
  root_txt <- "rooted as in the HGT analysis"
  track <- if (tag == "PF13354_PF00768") "mafft" else "hmm"
  f <- file.path(HGT_DIR, "Rooting_Decisions_And_Agreement.tsv")
  if (file.exists(f)) {
    r <- suppressMessages(readr::read_tsv(f, show_col_types = FALSE))
    r <- r[r$family == tag, , drop = FALSE]
    if (nrow(r) == 1) {
      root_txt <- paste0(gsub("_", " ", r$primary_rooting[1]), " root")
      if ("primary_alignment_track" %in% names(r) && !is.na(r$primary_alignment_track[1])) {
        track <- r$primary_alignment_track[1]
      }
    }
  }
  aln_txt <- c(hmm = "hmmalign", mafft = "MAFFT+ClipKIT", struct = "FoldMason-guided hmmalign")[[track]]
  model <- "model n/a"
  iq <- file.path(TREE_DIR, sprintf("%s.primary.%s_out.iqtree", tag, track))
  if (file.exists(iq)) {
    hit <- grep("^(Best-fit model according to BIC|Model of substitution):", readLines(iq, warn = FALSE), value = TRUE)
    if (length(hit) > 0) model <- trimws(sub("^[^:]+:", "", hit[1]))
  }
  paste(model, aln_txt, root_txt, sep = " | ")
}

# 1.2 Layout, thresholds and resolution
OPEN_ANGLE         <- 10      # Opening angle of the circular layout (degrees)
CAP_DEPTH_QUANTILE <- 0.985   # Cap the longest 1.5% of branches to expand the ingroup topology (1.0 = no capping)
SHOW_BOOT_NODES    <- TRUE    # Mark well-supported nodes (UFBoot >= 95 & SH-aLRT >= 80)
SHOW_SBL_LABELS    <- TRUE    # Label representative viral CARD resistance genes on the outer ring
MIN_LABEL_TIP_GAP  <- 18      # Minimum tip spacing between outer-ring gene labels (avoids overlaps)
FIG_WIDTH          <- 13.0    # Single-tree figure width (inches)
FIG_HEIGHT         <- 10.2    # Single-tree figure height (inches)
SUPP_ROW_IN        <- 0.011   # Row height per tip in the rectangular domain-context figure (inches)
SUPP_WIDTH         <- 14.0    # Width of the rectangular domain-context figure (inches)
TIFF_DPI           <- 300     # Use 600 for final submission (300 dpi balances sharpness and file size)
EXPORT_PNG_PREVIEW <- TRUE    # Also write PNG previews

# 1.3 Annotation thresholds
VIRAL_CLADE_FRAC   <- 0.80                      # Minimum viral descendant fraction for a viral branch
MOSAIC_RING1_LEVEL <- "Viral Tier 2 (supported)" # MOSAiC environmental vMAGs are grouped with supported viral loci
CARD_MIN_PIDENT    <- 40                        # Minimum identity to a CARD reference resistance gene (%)
CARD_MIN_COV       <- 70                        # Minimum coverage of a CARD reference resistance gene (%)
ARC_MIN_TIPS       <- 5                         # Minimum tips for a functional clade arc
ARC_MIN_TIPS_PHAGE <- 3                         # Minimum tips for a phage tail-fibre / lysin fusion clade arc
ARC_MIN_FRAC       <- 0.80                      # Minimum fraction of one functional class within a clade arc

# 1.4 Branch colours (viral vs bacterial lineages only)
PAL_BRANCH <- c(
  "Viral lineage"     = "#D62728",
  "Bacterial lineage" = "#94A3B8"
)
LW_BRANCH <- c(
  "Viral lineage"     = 0.40,
  "Bacterial lineage" = 0.18
)

# 1.5 Ring 1: viral confidence tier and bacterial reference (4 levels)
PAL_RING1 <- c(
  "Viral Tier 1 (high-confidence)" = "#B91C1C",
  "Viral Tier 2 (supported)"       = "#F97316",
  "Viral Tier 3 (fragment)"        = "#FCD34D",
  "Bacterial genome"               = "#E2E8F0"
)

# 1.6 Outer CARD ring (known or close CARD reference genes only)
PAL_CARD <- c(
  "Bacterial CARD \u03b2-lactamase" = "#334155",
  "Viral CARD \u03b2-lactamase"     = "#E11D48"
)

# 1.7 Joint-tree ring: Pfam family
PAL_RING3_JOINT <- c(
  "PF13354 (Beta-lactamase2)" = "#4F46E5",
  "PF00768 (Peptidase_S11)"   = "#06B6D4"
)

# 1.8 Pfam domain vocabulary and mapping
CATALYTIC_DOMAINS <- c("Beta-lactamase", "Beta-lactamase2", "Transpeptidase", "Peptidase_S11")
DOM_PHAGE_STRUCT  <- c("Phage_tail_CBM", "Phage_spike_CBM", "TMP_3", "Gly_rich_2")
DOM_LYSIN         <- c(
  "Endolysin", "Phage_lysozyme", "VanY", "Peptidase_M23",
  "Amidase_2", "Amidase_3", "Amidase_6", "Glyco_hydro_25",
  "Glucosaminidase", "SLT", "CHAP", "Peptidase_M15", "Polysacc_deac_1"
)
DOM_CW_BIND <- c(
  "LysM", "PG_binding_1", "PG_binding_3", "PG_binding_4", "SLH", "SH3_3",
  "Choline_bind_1", "Choline_bind_2", "Choline_bind_3", "GW", "ChW", "AMIN", "Dockerin_1"
)
DOM_PG_RECYC <- c("Glyco_hydro_3", "Glyco_hydro_3_C", "NamZ_N", "NamZ_C")
DOM_TG       <- c("Transgly")
DOM_CLASSB   <- c("PBP_dimer", "PBP_dimer_2", "MrdA_TM_N", "MecA_N", "Pbp2x_N")
DOM_CTERM    <- c("PASTA", "PBP_1A_C", "PCB_OB", "PBP_middle", "BiPBP_C")
DOM_LMW      <- c("PBP5_C", "PBP4_C", "SPOR")
DOM_SENSOR   <- c("Peptidase_M56")

# 1.9 Main-figure outer arcs & supplementary Function strip colours (function inferred from domain context)
PAL_FUNC <- c(
  "Phage tail/lysin fusion (putative PG-lytic, non-\u03b2-lactamase)" = "#9333EA",
  "Cell-wall binding/hydrolase fusion"  = "#84CC16",
  "PG-recycling fusion (GH3/NamZ)"      = "#0EA5E9",
  "Class A HMW-PBP (TG + TP)"           = "#2563EB",
  "Class B HMW-PBP (dimer + TP)"        = "#10B981",
  "LMW-PBP (PBP5_C/SPOR)"               = "#F59E0B"
)

# 1.10 Per-row domain group colours in the rectangular supplementary figure
PAL_DOMGRP <- c(
  "SBL-fold catalytic domain"            = "#475569",
  "Transglycosylase (TG)"                = "#2563EB",
  "PBP dimerisation / N-term (class B)"  = "#10B981",
  "PBP C-terminal ext. (PASTA etc.)"     = "#FCD34D",
  "LMW-PBP C-term (PBP5_C/SPOR)"         = "#F59E0B",
  "Phage tail/spike CBM, TMP"            = "#9333EA",
  "PG hydrolase / lysin"                 = "#E879F9",
  "Cell-wall binding"                    = "#84CC16",
  "PG recycling (GH3/NamZ)"              = "#0EA5E9",
  "Other / DUF"                          = "#CBD5E1"
)

domain_group <- function(name) {
  dplyr::case_when(
    name %in% CATALYTIC_DOMAINS ~ "SBL-fold catalytic domain",
    name %in% DOM_TG            ~ "Transglycosylase (TG)",
    name %in% DOM_CLASSB        ~ "PBP dimerisation / N-term (class B)",
    name %in% DOM_CTERM         ~ "PBP C-terminal ext. (PASTA etc.)",
    name %in% DOM_LMW           ~ "LMW-PBP C-term (PBP5_C/SPOR)",
    name %in% DOM_PHAGE_STRUCT  ~ "Phage tail/spike CBM, TMP",
    name %in% DOM_LYSIN         ~ "PG hydrolase / lysin",
    name %in% DOM_CW_BIND       ~ "Cell-wall binding",
    name %in% DOM_PG_RECYC      ~ "PG recycling (GH3/NamZ)",
    name %in% DOM_SENSOR        ~ "Other / DUF",
    TRUE                        ~ "Other / DUF"
  )
}

# ==============================================================================
# SECTION 2: Read and build one normalised annotation metadata table (reusable in R)
# ==============================================================================

clean_card_family <- function(card_family, card_pident, card_cov, func_tier) {
  ok <- !is.na(card_pident) & !is.na(card_cov) &
        card_pident >= CARD_MIN_PIDENT & card_cov >= CARD_MIN_COV &
        func_tier == "Putative_Hydrolytic_BetaLactamase"
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
  out[env_idx] <- "Env-SBL (Triad)"
  out
}

parse_domains <- function(arch_str) {
  if (is.na(arch_str) || !nzchar(trimws(arch_str))) {
    return(tibble(name = character(0), start = numeric(0), end = numeric(0)))
  }
  m <- stringr::str_match_all(arch_str, "(PF\\d+):([^()+]+)\\((\\d+)-(\\d+)\\)")[[1]]
  if (nrow(m) == 0) {
    return(tibble(name = character(0), start = numeric(0), end = numeric(0)))
  }
  tibble(
    name  = m[, 3],
    start = as.numeric(m[, 4]),
    end   = as.numeric(m[, 5])
  )
}

classify_func_class <- function(arch_str, origin, func_tier) {
  d <- parse_domains(arch_str)$name
  if (length(d) == 0) return(NA_character_)
  non_hydrolytic <- is.na(func_tier) || !identical(func_tier, "Putative_Hydrolytic_BetaLactamase")
  if (identical(origin, "Viral") && non_hydrolytic && any(d %in% c(DOM_PHAGE_STRUCT, DOM_LYSIN, DOM_CW_BIND))) {
    return("Phage tail/lysin fusion (putative PG-lytic, non-\u03b2-lactamase)")
  }
  if (any(d %in% c(DOM_LYSIN, DOM_CW_BIND, DOM_PHAGE_STRUCT))) {
    return("Cell-wall binding/hydrolase fusion")
  }
  if (any(d %in% DOM_PG_RECYC)) return("PG-recycling fusion (GH3/NamZ)")
  if (any(d %in% DOM_TG))       return("Class A HMW-PBP (TG + TP)")
  if (any(d %in% DOM_CLASSB))   return("Class B HMW-PBP (dimer + TP)")
  if (any(d %in% DOM_LMW))      return("LMW-PBP (PBP5_C/SPOR)")
  NA_character_
}

as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

build_merged_metadata <- function() {
  message("[Metadata] Loading 2026 curated viral, bacterial, and Phase 2 neighborhood catalogs...")
  v_raw <- read_tsv(file.path(META_DIR, "Final_Viral_SBL_Catalog.tsv"), show_col_types = FALSE)
  b_raw <- read_tsv(file.path(META_DIR, "Final_Bacterial_TripleSource_Catalog.tsv"), show_col_types = FALSE)
  # Phase 2 tables written before 2026-11-25 name the hallmark columns n_viral_hallmark_in_region and
  # flanked_strict_structural_hallmark_both (annotation-text hallmarks); both layouts are accepted.
  p2_raw <- read_tsv(
    file.path(META_DIR, "Phase2_Viral_Genomic_Neighborhoods_Full.tsv"),
    col_select = any_of(c(
      "uid", "contig_length_bp", "genomad_virus_score", "checkv_viral_genes",
      "checkv_host_genes", "n_genomad_hallmark_in_region", "n_viral_hallmark_in_region",
      "flanked_genomad_hallmark_both", "flanked_strict_structural_hallmark_both",
      "upstream_5_orfs", "downstream_5_orfs"
    )),
    show_col_types = FALSE
  ) |>
    rename(any_of(c(n_genomad_hallmark_in_region = "n_viral_hallmark_in_region",
                    flanked_genomad_hallmark_both = "flanked_strict_structural_hallmark_both")))

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
           n_genomad_hallmark_in_region, flanked_genomad_hallmark_both, viral_mge_context)

  v_df <- v_raw |>
    filter(as_bool(pass_cov65), as_bool(pass_fimo_1e5)) |>
    left_join(v_mge, by = "uid") |>
    mutate(
      origin = "Viral",
      seq_len = as.numeric(seq_len),
      card_verified = !is.na(card_pident) & !is.na(card_cov) &
                      as.numeric(card_pident) >= CARD_MIN_PIDENT &
                      as.numeric(card_cov) >= CARD_MIN_COV &
                      functional_tier == "Putative_Hydrolytic_BetaLactamase",
      ring1_origin = case_when(
        source == "MOSAiC" | grepl("MOSAiC", viral_origin_tier) ~ MOSAIC_RING1_LEVEL,
        viral_origin_tier == "Tier1_HighConfidence_Viral"       ~ "Viral Tier 1 (high-confidence)",
        viral_origin_tier == "Tier2_Supported_Viral"            ~ "Viral Tier 2 (supported)",
        TRUE                                                    ~ "Viral Tier 3 (fragment)"
      ),
      ring3_joint_fam = ifelse(winning_acc == "PF13354", "PF13354 (Beta-lactamase2)", "PF00768 (Peptidase_S11)"),
      sbl_display_label = clean_card_family(card_gene_family, as.numeric(card_pident), as.numeric(card_cov), functional_tier),
      taxon_group = coalesce(viral_mge_context, "MOSAiC_vMAG")
    )

  b_df <- b_raw |>
    filter(as_bool(pass_cov65), as_bool(pass_fimo_1e5), !as_bool(is_viral_contig_in_genome)) |>
    mutate(
      origin = "Bacterial",
      seq_len = as.numeric(seq_len),
      card_verified = !is.na(card_pident) & !is.na(card_cov) &
                      as.numeric(card_pident) >= CARD_MIN_PIDENT &
                      as.numeric(card_cov) >= CARD_MIN_COV &
                      functional_tier == "Putative_Hydrolytic_BetaLactamase",
      ring1_origin = "Bacterial genome",
      ring3_joint_fam = ifelse(winning_acc == "PF13354", "PF13354 (Beta-lactamase2)", "PF00768 (Peptidase_S11)"),
      sbl_display_label = clean_card_family(card_gene_family, as.numeric(card_pident), as.numeric(card_cov), functional_tier),
      taxon_group = phylum,
      viral_mge_context = NA_character_,
      viral_origin_tier = NA_character_,
      viral_cluster_size = NA_integer_,
      contig_length_bp = NA_real_,
      checkv_viral_genes = NA_real_,
      checkv_host_genes = NA_real_,
      n_genomad_hallmark_in_region = NA_real_,
      flanked_genomad_hallmark_both = NA
    )

  keep_cols <- c(
    "uid", "origin", "winning_acc", "functional_tier", "functional_evidence",
    "card_verified", "card_gene_family", "card_pident", "card_cov", "sbl_display_label",
    "ring1_origin", "ring3_joint_fam",
    "taxon_group", "viral_origin_tier", "viral_cluster_size", "viral_mge_context",
    "domain_architecture", "seq_len", "hmm_cov", "motif1_SXXK_hmm", "motif2_class", "motif3_KTG_hmm", "has_omega_E166"
  )
  merged <- bind_rows(select(v_df, all_of(keep_cols)), select(b_df, all_of(keep_cols))) |>
    distinct(uid, .keep_all = TRUE) |>
    mutate(
      func_class = mapply(classify_func_class, domain_architecture, origin, functional_tier, USE.NAMES = FALSE)
    )

  out_tsv <- file.path(TREE_DIR, "Tree_Tip_Annotations_Merged.tsv")
  write_tsv(merged, out_tsv)
  message("  -> Saved unified tip metadata table: ", out_tsv, " (n=", nrow(merged), ")")
  merged
}

# ==============================================================================
# SECTION 3: Rooting, long-branch capping, branch states and domain-based functional clade arcs
# ==============================================================================

find_func_clades <- function(tr, tip_class, min_frac = ARC_MIN_FRAC) {
  n_tips <- length(tr$tip.label)
  tot <- n_tips + tr$Nnode
  parent <- integer(tot)
  parent[tr$edge[, 2]] <- tr$edge[, 1]
  po <- ape::reorder.phylo(tr, "postorder")$edge
  n_desc <- integer(tot)
  n_desc[seq_len(n_tips)] <- 1L
  for (i in seq_len(nrow(po))) {
    n_desc[po[i, 1]] <- n_desc[po[i, 1]] + n_desc[po[i, 2]]
  }
  out <- list()
  for (cl in unique(na.omit(tip_class))) {
    mt <- if (startsWith(cl, "Phage tail/lysin")) ARC_MIN_TIPS_PHAGE else ARC_MIN_TIPS
    cnt <- integer(tot)
    cnt[seq_len(n_tips)] <- as.integer(!is.na(tip_class) & tip_class == cl)
    for (i in seq_len(nrow(po))) {
      cnt[po[i, 1]] <- cnt[po[i, 1]] + cnt[po[i, 2]]
    }
    q <- seq_len(tot) > n_tips & n_desc >= mt & (cnt / pmax(n_desc, 1)) >= min_frac
    qn <- which(q)
    keep <- vapply(qn, function(nd) {
      a <- parent[nd]
      while (a > 0) {
        if (q[a]) return(FALSE)
        a <- parent[a]
      }
      TRUE
    }, logical(1))
    if (any(keep)) {
      out[[cl]] <- tibble(
        node       = qn[keep],
        func_class = cl,
        n_tips     = n_desc[qn[keep]],
        n_class    = cnt[qn[keep]]
      )
    }
  }
  bind_rows(out)
}

prepare_tree_and_nodes <- function(tree_path, outgroup_path, meta_df, tag) {
  tr <- ape::read.tree(tree_path)
  if (!is.null(outgroup_path) && file.exists(outgroup_path)) {
    og <- read_lines(outgroup_path)
    og <- trimws(og[nzchar(trimws(og))])
    drop_ids <- intersect(tr$tip.label, og)
    if (length(drop_ids) > 0) {
      tr <- ape::drop.tip(tr, drop_ids)
    }
  }
  if (!ape::is.rooted(tr)) {
    warning("[Root fallback] ", basename(tree_path), " is unrooted; falling back to midpoint().")
    tr <- phangorn::midpoint(tr)
  }
  tr <- ape::ladderize(tr, right = FALSE)

  # Cap the longest 1.5% of branches (parent depth <= child depth, so branch lengths stay non-negative)
  if (CAP_DEPTH_QUANTILE < 1.0 && !is.null(tr$edge.length)) {
    depths <- ape::node.depth.edgelength(tr)
    tip_depths <- depths[seq_len(length(tr$tip.label))]
    q_cap <- as.numeric(quantile(tip_depths, probs = CAP_DEPTH_QUANTILE, na.rm = TRUE))
    if (is.finite(q_cap) && q_cap > 0) {
      d_cap <- pmin(depths, q_cap)
      tr$edge.length <- pmax(0, d_cap[tr$edge[, 2]] - d_cap[tr$edge[, 1]])
    }
  }

  # Parse IQ-TREE SH-aLRT / UFBoot support
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

  # Count viral descendants of every branch bottom-up
  tip_meta <- tibble(uid = tr$tip.label) |> left_join(meta_df, by = "uid")
  n_desc <- integer(tot_nodes)
  n_v    <- integer(tot_nodes)

  n_desc[seq_len(n_tips)] <- 1L
  n_v[seq_len(n_tips)]    <- ifelse(!is.na(tip_meta$origin) & tip_meta$origin == "Viral", 1L, 0L)

  tr_po <- ape::reorder.phylo(tr, "postorder")
  for (i in seq_len(nrow(tr_po$edge))) {
    p <- tr_po$edge[i, 1]
    c <- tr_po$edge[i, 2]
    n_desc[p] <- n_desc[p] + n_desc[c]
    n_v[p]    <- n_v[p] + n_v[c]
  }

  branch_group <- ifelse(
    (n_v / pmax(n_desc, 1)) >= VIRAL_CLADE_FRAC,
    "Viral lineage",
    "Bacterial lineage"
  )

  node_df <- tibble(
    node         = seq_len(tot_nodes),
    n_desc       = n_desc,
    n_viral_desc = n_v,
    sh_alrt      = sh_val,
    ufboot       = uf_val,
    high_support = (seq_len(tot_nodes) > n_tips) & !is.na(uf_val) & !is.na(sh_val) & uf_val >= 95 & sh_val >= 80 & n_desc >= 5,
    branch_group = factor(branch_group, levels = names(PAL_BRANCH))
  )

  tip_class <- tip_meta$func_class
  arc_df <- find_func_clades(tr, tip_class)

  if (nrow(arc_df) > 0) {
    arc_summary <- arc_df |>
      group_by(func_class) |>
      summarise(n_arcs = n(), tips_covered = sum(n_class), .groups = "drop")
    msg_parts <- sprintf("%s: %d arc(s), %d tips", arc_summary$func_class, arc_summary$n_arcs, arc_summary$tips_covered)
    message("  [Func clades | ", tag, "] ", paste(msg_parts, collapse = " ; "))
  } else {
    message("  [Func clades | ", tag, "] none (>= ", ARC_MIN_TIPS, " tips / ", ARC_MIN_TIPS_PHAGE, " phage, >= ", round(ARC_MIN_FRAC * 100), "%)")
  }

  tip_full <- tip_meta |>
    mutate(
      id              = uid,
      func_class      = factor(func_class, levels = names(PAL_FUNC)),
      ring1_origin    = factor(ring1_origin, levels = names(PAL_RING1)),
      card_ring       = factor(
        ifelse(!is.na(card_verified) & card_verified, paste0(origin, " CARD \u03b2-lactamase"), NA_character_),
        levels = names(PAL_CARD)
      ),
      ring3_joint_fam = factor(ring3_joint_fam, levels = names(PAL_RING3_JOINT))
    ) |>
    select(
      id, origin, functional_tier, card_verified, sbl_display_label,
      winning_acc, domain_architecture, seq_len, func_class,
      ring1_origin, card_ring, ring3_joint_fam
    )

  list(tree = tr, node_df = node_df, tip_df = tip_full, arc_df = arc_df)
}

# Choose outer-ring viral CARD gene labels (at least MIN_LABEL_TIP_GAP tips apart in tree order, to avoid overlaps)
select_spaced_sbl_labels <- function(tip_coords_df, min_gap = MIN_LABEL_TIP_GAP) {
  cand <- tip_coords_df |>
    filter(origin == "Viral", !is.na(card_verified) & card_verified) |>
    arrange(y)
  if (nrow(cand) == 0) return(cand)

  keep_idx <- integer(0)
  last_y <- -Inf
  for (i in seq_len(nrow(cand))) {
    if ((cand$y[i] - last_y) >= min_gap) {
      keep_idx <- c(keep_idx, i)
      last_y <- cand$y[i]
    }
  }
  cand[keep_idx, , drop = FALSE]
}

# ==============================================================================
# SECTION 4: Main circular trees (with functional clade arcs inferred from domain architecture)
# ==============================================================================

build_tree_panel <- function(spec, meta_df, mode = c("single", "composite"), panel_letter = NULL) {
  mode <- match.arg(mode)
  is_comp  <- (mode == "composite")
  is_joint <- identical(spec$tag, "PF13354_PF00768")

  tree_path <- file.path(TREE_DIR, spec$treefile)
  og_path   <- file.path(TREE_DIR, spec$outgroup)
  prep <- prepare_tree_and_nodes(tree_path, og_path, meta_df, spec$tag)
  tr      <- prep$tree
  node_df <- prep$node_df
  tip_df  <- prep$tip_df
  arc_df  <- prep$arc_df

  n_tot   <- nrow(tip_df)
  n_v     <- sum(tip_df$origin == "Viral", na.rm = TRUE)
  n_b     <- n_tot - n_v
  n_vcard <- sum(tip_df$origin == "Viral" & !is.na(tip_df$card_verified) & tip_df$card_verified, na.rm = TRUE)
  n_bcard <- sum(tip_df$origin == "Bacterial" & !is.na(tip_df$card_verified) & tip_df$card_verified, na.rm = TRUE)

  p <- ggtree(tr, layout = "fan", open.angle = OPEN_ANGLE, aes(color = branch_group, linewidth = branch_group)) %<+% node_df +
    scale_color_manual(name = "Branches", values = PAL_BRANCH, drop = FALSE, guide = guide_legend(order = 1)) +
    scale_linewidth_manual(
      name   = "Branches",
      values = if (is_comp) LW_BRANCH * 0.78 else LW_BRANCH,
      drop   = FALSE,
      guide  = guide_legend(order = 1)
    )

  if (SHOW_BOOT_NODES) {
    p <- p + geom_point2(aes(subset = high_support), size = 0.5, color = "#334155", alpha = 0.5, show.legend = FALSE)
  }

  M <- max(p$data$x, na.rm = TRUE)
  tip_coords <- p$data |>
    filter(isTip == TRUE) |>
    select(id = label, y, angle) |>
    left_join(tip_df, by = "id")

  # Ring 1: Viral origin confidence (4 levels)
  p <- p +
    new_scale_fill() +
    geom_tile(
      data = tip_coords,
      aes(x = M * 1.05, y = y, fill = ring1_origin),
      width = M * 0.06,
      height = 1,
      color = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name         = "Viral origin confidence",
      values       = PAL_RING1,
      drop         = FALSE,
      na.translate = FALSE,
      guide        = guide_legend(order = 2)
    )

  # Joint tree extra ring: Pfam family (placed between Ring 1 and CARD ring)
  if (is_joint) {
    p <- p +
      new_scale_fill() +
      geom_tile(
        data = tip_coords,
        aes(x = M * 1.115, y = y, fill = ring3_joint_fam),
        width = M * 0.04,
        height = 1,
        color = NA,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(
        name         = "Pfam family",
        values       = PAL_RING3_JOINT,
        drop         = FALSE,
        na.translate = FALSE,
        guide        = guide_legend(order = 4)
      )
  }

  # CARD reference ring (drawn only at CARD-verified tips)
  card_x <- if (is_joint) M * 1.165 else M * 1.115
  card_tips <- tip_coords |> filter(!is.na(card_ring))
  p <- p +
    new_scale_fill() +
    geom_tile(
      data = card_tips,
      aes(x = card_x, y = y, fill = card_ring),
      width = M * 0.04,
      height = 1,
      color = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name         = "CARD reference (\u226540% id, \u226570% cov)",
      values       = PAL_CARD,
      drop         = FALSE,
      na.translate = FALSE,
      guide        = guide_legend(order = 3)
    )

  # Domain-inferred function arcs
  arc_x0 <- card_x + M * 0.04
  arc_x1 <- arc_x0 + M * 0.022

  if (nrow(arc_df) > 0) {
    arc_rects <- lapply(seq_len(nrow(arc_df)), function(i) {
      nd <- arc_df$node[i]
      cl <- arc_df$func_class[i]
      desc_nodes <- phangorn::Descendants(tr, nd, "tips")[[1]]
      ys <- p$data$y[match(desc_nodes, p$data$node)]
      tibble(
        xmin       = arc_x0,
        xmax       = arc_x1,
        ymin       = min(ys, na.rm = TRUE) - 0.45,
        ymax       = max(ys, na.rm = TRUE) + 0.45,
        y_mid      = (min(ys, na.rm = TRUE) + max(ys, na.rm = TRUE)) / 2,
        n_tips     = arc_df$n_tips[i],
        func_class = factor(cl, levels = names(PAL_FUNC))
      )
    }) |> bind_rows()
  } else {
    arc_rects <- tibble(
      xmin       = numeric(0),
      xmax       = numeric(0),
      ymin       = numeric(0),
      ymax       = numeric(0),
      y_mid      = numeric(0),
      n_tips     = integer(0),
      func_class = factor(character(0), levels = names(PAL_FUNC))
    )
  }

  p <- p +
    new_scale_fill() +
    geom_rect(
      data = arc_rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = func_class),
      colour = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name         = "Domain-inferred function",
      values       = PAL_FUNC,
      drop         = FALSE,
      na.translate = FALSE,
      guide        = guide_legend(order = 5, ncol = 2)
    )

  # Text callout ONLY for phage tail/lysin arcs (merge adjacent arcs within 3% of tips so label is never doubled)
  phage_arcs <- arc_rects |>
    filter(startsWith(as.character(func_class), "Phage tail/lysin")) |>
    arrange(y_mid)
  if (nrow(phage_arcs) > 0) {
    gap_thr <- max(40, round(n_tot * 0.03))
    grp_id <- integer(nrow(phage_arcs))
    grp_id[1] <- 1L
    if (nrow(phage_arcs) > 1) {
      for (i in 2:nrow(phage_arcs)) {
        grp_id[i] <- grp_id[i - 1] + as.integer((phage_arcs$ymin[i] - phage_arcs$ymax[i - 1]) > gap_thr)
      }
    }
    phage_lbls <- phage_arcs |>
      mutate(grp = grp_id, span = ymax - ymin) |>
      group_by(grp) |>
      summarise(y_mid = weighted.mean(y_mid, w = pmax(span, 1)), .groups = "drop") |>
      mutate(
        nearest_idx = vapply(y_mid, function(ym) which.min(abs(tip_coords$y - ym)), integer(1)),
        ang         = tip_coords$angle[nearest_idx] %% 360,
        flip        = ang > 90 & ang < 270,
        text_angle  = ifelse(flip, ang - 180, ang),
        text_hjust  = ifelse(flip, 1, 0),
        callout     = "Phage tail/lysin"
      )
    p <- p +
      geom_text(
        data = phage_lbls,
        aes(x = arc_x1 + M * 0.03, y = y_mid, label = callout, angle = text_angle, hjust = text_hjust),
        color = "#7E22CE",
        size = if (is_comp) 2.6 else 3.0,
        fontface = "bold",
        inherit.aes = FALSE
      )
  }

  # Viral CARD stars outside the domain arcs
  star_x <- arc_x1 + M * 0.06
  vsbl_pts <- tip_coords |> filter(origin == "Viral", !is.na(card_verified) & card_verified)
  if (nrow(vsbl_pts) > 0) {
    p <- p +
      geom_point(
        data = vsbl_pts,
        aes(x = star_x, y = y),
        shape = 8,
        color = "#E11D48",
        size = if (is_comp) 1.8 else 2.3,
        stroke = 0.45,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

  # Sparse viral CARD gene-family labels
  if (SHOW_SBL_LABELS) {
    lbl_df <- select_spaced_sbl_labels(tip_coords, min_gap = max(14, round(n_tot / 95)))
    if (nrow(lbl_df) > 0) {
      lbl_df <- lbl_df |>
        mutate(
          ang        = (angle %% 360),
          flip       = ang > 90 & ang < 270,
          text_angle = ifelse(flip, ang - 180, ang),
          text_hjust = ifelse(flip, 1, 0)
        )
      p <- p +
        geom_text(
          data = lbl_df,
          aes(x = arc_x1 + M * 0.085, y = y, label = sbl_display_label, angle = text_angle, hjust = text_hjust),
          color = "#9D174D",
          size = if (is_comp) 2.3 else 2.6,
          fontface = "bold",
          inherit.aes = FALSE
        )
    }
  }

  base_sub <- sprintf(
    "n = %s tips (viral %s | bacterial %s) \u00b7 CARD hits: viral %d | bacterial %d",
    format(n_tot, big.mark = ","),
    format(n_v, big.mark = ","),
    format(n_b, big.mark = ","),
    n_vcard, n_bcard
  )
  title_txt <- if (is_comp && !is.null(panel_letter)) paste0(panel_letter, "  ", spec$short_title) else spec$title
  sub_txt   <- if (is_comp) base_sub else paste0(base_sub, " \u00b7 ", spec$model)

  p <- p +
    xlim(0, M * 1.50) +
    labs(title = title_txt, subtitle = sub_txt) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title       = element_text(size = if (is_comp) 11.5 else 14, face = "bold", color = "#0F172A", hjust = 0),
      plot.subtitle    = element_text(size = if (is_comp) 8.5 else 9.5, color = "#475569", hjust = 0),
      legend.position  = if (is_comp) "none" else "right",
      legend.title     = element_text(size = 8.8, face = "bold", color = "#0F172A"),
      legend.text      = element_text(size = 7.8, color = "#334155"),
      legend.key.size  = unit(0.40, "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      plot.margin      = margin(8, 8, 6, 8)
    )

  p
}

# ==============================================================================
# SECTION 5: Supplementary figures (rectangular tree + per-row domain-context map)
# ==============================================================================

build_domain_context_figure <- function(spec, meta_df) {
  stem <- paste0("FigS_DomainContext_", spec$tag)
  message("[Supplementary Domain Context] Building ", stem, " ...")

  tree_path <- file.path(TREE_DIR, spec$treefile)
  og_path   <- file.path(TREE_DIR, spec$outgroup)
  prep <- prepare_tree_and_nodes(tree_path, og_path, meta_df, spec$tag)
  tr      <- prep$tree
  node_df <- prep$node_df
  tip_df  <- prep$tip_df
  n       <- nrow(tip_df)

  top_pad <- max(48, round(n * 0.032))
  y_lims  <- c(0.5, n + top_pad)

  base_tree <- ggtree(tr, aes(color = branch_group), linewidth = 0.15) %<+% node_df +
    scale_color_manual(name = "Branches", values = PAL_BRANCH, drop = FALSE, guide = guide_legend(order = 1))

  M <- max(base_tree$data$x, na.rm = TRUE)
  tip_coords <- base_tree$data |>
    filter(isTip == TRUE) |>
    select(id = label, y) |>
    left_join(tip_df, by = "id")

  strip_hdr <- tibble(
    x   = c(M * 1.04, M * 1.08, M * 1.12),
    y   = n + 5,
    lab = c("Origin", "CARD", "Function")
  )

  tree_p <- base_tree +
    new_scale_fill() +
    geom_tile(
      data = tip_coords,
      aes(x = M * 1.04, y = y, fill = ring1_origin),
      width = M * 0.035,
      height = 1,
      colour = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name         = "Viral origin confidence",
      values       = PAL_RING1,
      drop         = FALSE,
      na.translate = FALSE,
      guide        = guide_legend(order = 2)
    ) +
    new_scale_fill() +
    geom_tile(
      data = tip_coords |> filter(!is.na(card_ring)),
      aes(x = M * 1.08, y = y, fill = card_ring),
      width = M * 0.035,
      height = 1,
      colour = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name         = "CARD reference (\u226540% id, \u226570% cov)",
      values       = PAL_CARD,
      drop         = FALSE,
      na.translate = FALSE,
      guide        = guide_legend(order = 3)
    ) +
    new_scale_fill() +
    geom_tile(
      data = tip_coords |> filter(!is.na(func_class)),
      aes(x = M * 1.12, y = y, fill = func_class),
      width = M * 0.035,
      height = 1,
      colour = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name         = "Domain-inferred function",
      values       = PAL_FUNC,
      drop         = FALSE,
      na.translate = FALSE,
      guide        = guide_legend(order = 4, ncol = 2)
    ) +
    geom_text(
      data = strip_hdr,
      aes(x = x, y = y, label = lab),
      size = 2.4,
      angle = 90,
      hjust = 0,
      fontface = "bold",
      color = "#0F172A",
      inherit.aes = FALSE
    ) +
    scale_y_continuous(limits = y_lims, expand = c(0, 0)) +
    xlim(0, M * 1.16) +
    theme_tree() +
    theme(plot.margin = margin(2, 2, 4, 6))

  # Domain panel: per-tip backbone + Pfam domain boxes
  XMAX <- ceiling(as.numeric(quantile(tip_coords$seq_len, 0.99, na.rm = TRUE)) / 100) * 100
  bb_df <- tip_coords |>
    transmute(y = y, L = pmin(seq_len, XMAX), truncated = seq_len > XMAX)

  dom_list <- lapply(seq_len(nrow(tip_coords)), function(i) {
    dt <- parse_domains(tip_coords$domain_architecture[i])
    if (nrow(dt) == 0) return(NULL)
    dt$y <- tip_coords$y[i]
    dt
  })
  dom_df <- bind_rows(dom_list) |>
    filter(start <= XMAX) |>
    mutate(
      end = pmin(end, XMAX),
      grp = factor(domain_group(name), levels = names(PAL_DOMGRP))
    )

  trunc_df <- bb_df |> filter(truncated)

  dom_p <- ggplot() +
    geom_segment(
      data = bb_df,
      aes(x = 1, xend = L, y = y, yend = y),
      colour = "#94A3B8",
      linewidth = 0.12
    ) +
    geom_rect(
      data = dom_df,
      aes(xmin = start, xmax = end, ymin = y - 0.42, ymax = y + 0.42, fill = grp),
      colour = NA
    ) +
    geom_point(
      data = trunc_df,
      aes(x = XMAX, y = y),
      shape = "|",
      size = 0.6,
      colour = "#0F172A"
    ) +
    scale_fill_manual(
      name   = "Domain group",
      values = PAL_DOMGRP,
      drop   = TRUE,
      guide  = guide_legend(order = 1, ncol = 4)
    ) +
    scale_x_continuous(
      name     = "Residue position (aa, N\u2192C)",
      limits   = c(0, XMAX),
      expand   = c(0.01, 0),
      position = "top"
    ) +
    scale_y_continuous(limits = y_lims, expand = c(0, 0)) +
    theme_minimal(base_size = 8) +
    theme(
      axis.title.y       = element_blank(),
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = "#EEF2F6", linewidth = 0.35),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.margin        = margin(2, 8, 4, 2)
    )

  # Align panels and title with patchwork; keep the full ggnewscale bottom legend via cowplot
  leg_theme <- theme(
    legend.position      = "bottom",
    legend.box           = "vertical",
    legend.direction     = "horizontal",
    legend.justification = "center",
    legend.box.just      = "center",
    legend.title         = element_text(size = 8.0, face = "bold", color = "#0F172A"),
    legend.text          = element_text(size = 7.0, color = "#334155"),
    legend.key.size      = unit(0.32, "cm"),
    legend.spacing.y     = unit(0.04, "cm"),
    legend.box.spacing   = unit(2, "pt"),
    legend.margin        = margin(1, 2, 1, 2)
  )

  leg_tree <- suppressWarnings(cowplot::get_legend(tree_p + leg_theme))
  leg_dom  <- suppressWarnings(cowplot::get_legend(dom_p + leg_theme))
  leg_all  <- cowplot::plot_grid(leg_tree, leg_dom, ncol = 1, rel_heights = c(1.5, 1))

  row_pw <- (tree_p + theme(legend.position = "none") | dom_p + theme(legend.position = "none")) +
    patchwork::plot_layout(widths = c(1.15, 1)) +
    patchwork::plot_annotation(
      title    = sprintf("%s: domain context of every tip (rectangular ML tree)", spec$title),
      subtitle = sprintf(
        "n = %s tips; each row = one protein drawn from residue 1; proteins > %d aa truncated (|).",
        format(n, big.mark = ","), XMAX
      ),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold", color = "#0F172A", hjust = 0),
        plot.subtitle = element_text(size = 9.0, color = "#475569", hjust = 0),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )

  fig_h <- max(8.0, n * SUPP_ROW_IN + 2.5)
  fig <- cowplot::plot_grid(
    patchwork::patchworkGrob(row_pw),
    leg_all,
    ncol = 1,
    rel_heights = c(fig_h - 2.3, 2.3)
  )

  save_multi_format(fig, stem, width = SUPP_WIDTH, height = fig_h)
}

# ==============================================================================
# SECTION 6: Generate and export main, composite and supplementary domain-context figures (PDF + LZW-TIFF + PNG)
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

  # LZW-compressed TIFF (journal submission standard)
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

extract_bottom_legend <- function(spec, meta_df) {
  leg_panel <- build_tree_panel(spec, meta_df, mode = "single") +
    theme(
      legend.position      = "bottom",
      legend.box           = "vertical",
      legend.direction     = "horizontal",
      legend.justification = "center",
      legend.box.just      = "center",
      legend.title         = element_text(size = 8.2, face = "bold", color = "#0F172A"),
      legend.text          = element_text(size = 7.4, color = "#334155"),
      legend.key.size      = unit(0.32, "cm"),
      legend.spacing.y     = unit(0.03, "cm"),
      legend.box.spacing   = unit(2, "pt"),
      legend.margin        = margin(1, 2, 1, 2)
    )
  leg <- suppressWarnings(cowplot::get_legend(leg_panel))
  if (is.null(leg) || length(leg$grobs) == 0) {
    comps <- suppressWarnings(cowplot::get_plot_component(leg_panel, "guide-box-bottom", return_all = TRUE))
    if (length(comps) > 0) leg <- comps[[1]]
  }
  leg
}

make_composite <- function(specs_sub, letters_sub, legend_from_idx, meta_df, stem, width, height, leg_rel_h = 0.17) {
  message("[Composite] Building ", stem, " ...")
  panels <- lapply(seq_along(specs_sub), function(i) {
    build_tree_panel(specs_sub[[i]], meta_df, mode = "composite", panel_letter = letters_sub[i])
  })
  row_grid   <- cowplot::plot_grid(plotlist = panels, nrow = 1, align = "h", axis = "tb")
  shared_leg <- extract_bottom_legend(specs_sub[[legend_from_idx]], meta_df)
  full_fig   <- cowplot::plot_grid(row_grid, shared_leg, ncol = 1, rel_heights = c(1, leg_rel_h))
  save_multi_format(full_fig, stem, width = width, height = height)
}

main <- function() {
  meta_df <- build_merged_metadata()

  specs <- list(
    list(
      tag = "PF00144", stem = "PF00144_primary_tree",
      title = "PF00144 (Beta-lactamase / Class C-like)", short_title = "PF00144 (class C-like)",
      treefile = "PF00144.primary.rooted.treefile", outgroup = "PF00144_outgroup.txt",
      model = tree_method_label("PF00144")
    ),
    list(
      tag = "PF00905", stem = "PF00905_primary_tree",
      title = "PF00905 (Transpeptidase / HMW-PBP & Class D-like)", short_title = "PF00905 (transpeptidase)",
      treefile = "PF00905.primary.rooted.treefile", outgroup = "PF00905_outgroup.txt",
      model = tree_method_label("PF00905")
    ),
    list(
      tag = "PF13354", stem = "PF13354_primary_tree",
      title = "PF13354 (Beta-lactamase2 / Class A-like)", short_title = "PF13354 (class A-like)",
      treefile = "PF13354.primary.rooted.treefile", outgroup = "PF13354_outgroup.txt",
      model = tree_method_label("PF13354")
    ),
    list(
      tag = "PF00768", stem = "PF00768_primary_tree",
      title = "PF00768 (Peptidase_S11 / Sister LMW-PBP)", short_title = "PF00768 (Peptidase_S11)",
      treefile = "PF00768.primary.rooted.treefile", outgroup = "PF00768_outgroup.txt",
      model = tree_method_label("PF00768")
    ),
    list(
      tag = "PF13354_PF00768", stem = "PF13354_PF00768_joint_tree",
      title = "Joint PF13354 + PF00768 (Beta-lactamase2 vs. Peptidase_S11)", short_title = "Joint PF13354 + PF00768",
      treefile = "PF13354_PF00768.primary.rooted.treefile", outgroup = "PF13354_PF00768_outgroup.txt",
      model = tree_method_label("PF13354_PF00768")
    )
  )

  # 1) Single-family circular trees
  for (sp in specs) {
    message("[Plotting] ", sp$stem, " ...")
    p <- build_tree_panel(sp, meta_df, mode = "single")
    save_multi_format(p, sp$stem, width = FIG_WIDTH, height = FIG_HEIGHT)
  }

  # 2) Main figure: PF00144 + PF00905 side by side
  make_composite(
    specs_sub = specs[c(1, 2)],
    letters_sub = c("A", "B"),
    legend_from_idx = 1,
    meta_df = meta_df,
    stem = "Fig_Main_PF00144_PF00905",
    width = 15.0,
    height = 9.2,
    leg_rel_h = 0.17
  )

  # 3) Supplementary circular composite: PF13354 + PF00768 + joint PF13354+PF00768
  make_composite(
    specs_sub = specs[c(3, 4, 5)],
    letters_sub = c("A", "B", "C"),
    legend_from_idx = 3,
    meta_df = meta_df,
    stem = "FigS_PF13354_PF00768_Joint",
    width = 21.0,
    height = 9.3,
    leg_rel_h = 0.21
  )

  # 4) Supplementary domain-context figures: rectangular tree + per-sequence domain map
  for (sp in specs[1:4]) {
    build_domain_context_figure(sp, meta_df)
  }

  message("[Done] All R tree figures saved to: ", OUT_DIR)
}

main()
