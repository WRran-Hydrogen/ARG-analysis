#!/usr/bin/env Rscript
# ==============================================================================
# CL0013 superfamily + leave-one-outgroup-out validation and HGT/OAT figures
#
# Usage:
#   Rscript scripts/visualization/plot_phylogeny_cl0013_and_hgt_2026.R
#
# Outputs in figures/trees_2026_R/ (.pdf, .tiff @ 300 dpi, .png):
#   1. Fig_CL0013_Superfamily_and_LOO
#   2. Fig_Phylogenetic_HGT_and_OAT_Robustness
# ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggtree)
  library(ggnewscale)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(cowplot)
  library(ragg)
})

# ==============================================================================
# 1. Path & Parameter Configuration
# ==============================================================================
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
  stop("Could not locate SBL_research repository root.")
}

REPO_ROOT  <- locate_repo_root()
CL_DIR     <- file.path(REPO_ROOT, "results", "reanalysis_2026", "03_primary_trees", "cl0013_superfamily")
HGT_DIR    <- file.path(REPO_ROOT, "results", "reanalysis_2026", "04_hgt_inference")
OUT_DIR    <- file.path(REPO_ROOT, "figures", "trees_2026_R")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TIFF_DPI   <- 300

# All 5 PMSF C20 jobs evaluated when completed and verified on HPC; keep character(0)
PENDING_PMSF_FAMS <- character(0)

# High-contrast, colorblind-safe CL0013 8-family palette
FAM_ORDER <- c("PF00144", "PF00905", "PF13354", "PF00768", "PF02113", "PF04960", "PF06089", "PF27310")
FAM_LABELS <- c(
  "PF00144" = "PF00144 (Beta-lactamase, Class C)",
  "PF00905" = "PF00905 (Transpeptidase, HMW-PBP)",
  "PF13354" = "PF13354 (Beta-lactamase2, Class A)",
  "PF00768" = "PF00768 (Peptidase_S11, LMW-PBP)",
  "PF02113" = "PF02113 (Peptidase_S13, PRE)",
  "PF04960" = "PF04960 (Glutaminase, Outgroup)",
  "PF06089" = "PF06089 (Asparaginase_II, Outgroup)",
  "PF27310" = "PF27310 (HrcA_C_epsilon, Outgroup)"
)
FAM_SHORT <- c(
  "PF00144" = "PF00144", "PF00905" = "PF00905", "PF13354" = "PF13354",
  "PF00768" = "PF00768", "PF02113" = "PF02113", "PF04960" = "PF04960",
  "PF06089" = "PF06089", "PF27310" = "PF27310"
)
PAL_CL0013_BRANCH <- c(
  "PF00144" = "#DC2626",
  "PF00905" = "#2563EB",
  "PF13354" = "#7C3AED",
  "PF00768" = "#0891B2",
  "PF02113" = "#D97706",
  "PF04960" = "#334155",
  "PF06089" = "#64748B",
  "PF27310" = "#94A3B8",
  "Backbone" = "#475569"
)
PAL_CL0013_RING <- PAL_CL0013_BRANCH[FAM_ORDER]

PAL_ORIGIN_FUNC <- c(
  "Viral PRE / SBL"                  = "#E11D48",
  "Bacterial Hydrolytic SBL (CARD)"  = "#1E293B",
  "Bacterial PBP / PRE"              = "#CBD5E1",
  "Non-PRE Outgroup"                 = "#F1F5F9"
)

PAL_DIR <- c(
  "Bacterial -> Viral" = "#E11D48",
  "Viral -> Bacterial" = "#0284C7"
)

# ==============================================================================
# 2. Shared Export & Tree Utilities
# ==============================================================================
as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

save_multi_format <- function(plot_obj, stem, width, height) {
  pdf_path  <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  tiff_path <- file.path(OUT_DIR, paste0(stem, ".tiff"))
  png_path  <- file.path(OUT_DIR, paste0(stem, ".png"))

  cairo_pdf(pdf_path, width = width, height = height, family = "sans")
  print(plot_obj)
  dev.off()

  ragg::agg_tiff(tiff_path, width = width, height = height, units = "in", res = TIFF_DPI,
                 background = "white", compression = "lzw")
  print(plot_obj)
  dev.off()

  ragg::agg_png(png_path, width = width, height = height, units = "in", res = TIFF_DPI, background = "white")
  print(plot_obj)
  dev.off()

  message("[Exported] ", stem, " -> .pdf / .tiff / .png (", width, "x", height, " in)")
}

parse_support <- function(lbl) {
  if (is.null(lbl) || is.na(lbl) || !nzchar(lbl)) return(c(sh = NA_real_, uf = NA_real_))
  parts <- strsplit(lbl, "/")[[1]]
  if (length(parts) >= 2) {
    return(c(sh = suppressWarnings(as.numeric(parts[1])), uf = suppressWarnings(as.numeric(parts[2]))))
  }
  c(sh = NA_real_, uf = suppressWarnings(as.numeric(parts[1])))
}

root_cl0013_safely <- function(tr, og_tips, meta_df) {
  og <- intersect(og_tips, tr$tip.label)
  if (length(og) > 0 && is.monophyletic(tr, og)) {
    return(ladderize(root(tr, outgroup = og, resolve.root = TRUE)))
  }
  # Fallback when a LOO outgroup split is non-monophyletic (e.g., loo_no_PF06089):
  # root on a single monophyletic outgroup family (prefer PF04960, then PF27310)
  for (f_try in c("PF04960", "PF27310", "PF06089")) {
    sub_og <- intersect(meta_df$uid[meta_df$cl0013_family == f_try], tr$tip.label)
    if (length(sub_og) >= 5 && is.monophyletic(tr, sub_og)) {
      return(ladderize(root(tr, outgroup = sub_og, resolve.root = TRUE)))
    }
  }
  ladderize(phangorn::midpoint(tr))
}

prepare_cl0013_tree <- function(tree_file, outgroup_file, meta_df) {
  tr <- read.tree(tree_file)
  og <- readLines(outgroup_file, warn = FALSE)
  og <- og[nzchar(trimws(og))]
  tr <- root_cl0013_safely(tr, og, meta_df)

  n_tips <- length(tr$tip.label)
  n_nodes <- n_tips + tr$Nnode
  desc_list <- phangorn::Descendants(tr, seq_len(n_nodes), type = "tips")

  tip_meta <- tibble(uid = tr$tip.label) |>
    left_join(meta_df, by = "uid") |>
    mutate(
      cl0013_family = factor(cl0013_family, levels = FAM_ORDER),
      origin_func = case_when(
        functional_tier == "CL0013_NonPRE_Outgroup"            ~ "Non-PRE Outgroup",
        origin == "Viral"                                      ~ "Viral PRE / SBL",
        functional_tier == "Putative_Hydrolytic_BetaLactamase" ~ "Bacterial Hydrolytic SBL (CARD)",
        TRUE                                                   ~ "Bacterial PBP / PRE"
      ),
      origin_func = factor(origin_func, levels = names(PAL_ORIGIN_FUNC))
    )

  fam_vec <- as.character(tip_meta$cl0013_family)
  branch_group <- character(n_nodes)
  for (i in seq_len(n_nodes)) {
    if (i <= n_tips) {
      branch_group[i] <- fam_vec[i]
    } else {
      f_sub <- fam_vec[desc_list[[i]]]
      tab <- sort(table(f_sub), decreasing = TRUE)
      if (length(tab) > 0 && (tab[1] / length(f_sub)) >= 0.70) {
        branch_group[i] <- names(tab)[1]
      } else {
        branch_group[i] <- "Backbone"
      }
    }
  }

  sup_mat <- t(sapply(c(rep("", n_tips), tr$node.label), parse_support))
  node_df <- tibble(
    node = seq_len(n_nodes),
    branch_group = factor(branch_group, levels = names(PAL_CL0013_BRANCH)),
    sh_alrt = sup_mat[, "sh"],
    ufboot  = sup_mat[, "uf"],
    high_support = (!is.na(sh_alrt) & !is.na(ufboot) & sh_alrt >= 80 & ufboot >= 95 & node > n_tips)
  )

  list(tree = tr, node_df = node_df, tip_df = tip_meta)
}

# Compute angular centroid of the contiguous core block of each family to avoid wraparound arcs
compute_core_family_labels <- function(tip_coords, fam_list) {
  res <- lapply(fam_list, function(f) {
    sub <- tip_coords |> filter(as.character(cl0013_family) == f) |> arrange(y)
    if (nrow(sub) == 0) return(NULL)
    # Use the tightest sliding 80% window along y to avoid 1-2 outlier tips skewing the midpoint
    k <- max(1L, floor(nrow(sub) * 0.80))
    if (nrow(sub) > k) {
      spans <- sapply(seq_len(nrow(sub) - k + 1L), function(s) sub$y[s + k - 1L] - sub$y[s])
      s0 <- which.min(spans)
      sub <- sub[s0:(s0 + k - 1L), ]
    }
    tibble(
      family = f,
      label = ifelse(f %in% c("PF04960", "PF06089", "PF27310"), paste0(f, "\n(Outgroup)"), f),
      y_mid = median(sub$y),
      angle_mid = median(sub$angle),
      n = nrow(tip_coords |> filter(as.character(cl0013_family) == f))
    )
  })
  bind_rows(res)
}

# ==============================================================================
# 3. Figure 1: CL0013 Superfamily Tree & Leave-One-Outgroup-Out (LOO) Validation
# ==============================================================================
build_cl0013_main_tree <- function(meta_df) {
  struct_tree <- file.path(CL_DIR, "CL0013.struct_c20.no_PF27310.treefile")
  if (file.exists(struct_tree) && file.exists(file.path(CL_DIR, "outgroup_struct_no_PF27310.txt"))) {
    use_tree <- struct_tree
    use_og <- file.path(CL_DIR, "outgroup_struct_no_PF27310.txt")
    main_title <- "A   CL0013 Superfamily Phylogeny (FoldMason + hmmalign, LG+C20+F+G, n = 448)"
    main_sub <- "5 PRE families (400 seqs) rooted on full-length non-PRE outgroups PF04960 + PF06089 (n = 48; short PF27310 removed)"
    active_fams <- setdiff(FAM_ORDER, "PF27310")
  } else {
    use_tree <- file.path(CL_DIR, "CL0013.all3_out.mafft_out.treefile")
    use_og <- file.path(CL_DIR, "outgroup_taxa_all3.txt")
    main_title <- "A   CL0013 Superfamily Phylogeny (All 3 Non-PRE Outgroups, n = 472)"
    main_sub <- "5 PRE families (400 seqs) cleanly resolved from 3 non-PRE CL0013 outgroups (PF04960, PF06089, PF27310; n = 72) | MAFFT + ClipKIT"
    active_fams <- FAM_ORDER
  }
  prep <- prepare_cl0013_tree(use_tree, use_og, meta_df)
  tr      <- prep$tree
  node_df <- prep$node_df
  tip_df  <- prep$tip_df |> rename(id = uid)

  p <- ggtree(tr, layout = "fan", open.angle = 14, aes(color = branch_group), linewidth = 0.42) %<+% node_df +
    scale_color_manual(
      name   = "CL0013 Family",
      values = PAL_CL0013_BRANCH,
      breaks = FAM_ORDER,
      labels = FAM_LABELS[FAM_ORDER],
      drop   = FALSE,
      guide  = guide_legend(order = 1, ncol = 2, override.aes = list(linewidth = 1.6))
    ) +
    geom_point2(aes(subset = high_support), size = 1.15, color = "#0F172A", alpha = 0.72, show.legend = FALSE)

  M <- max(p$data$x, na.rm = TRUE)
  tip_coords <- p$data |>
    filter(isTip == TRUE) |>
    select(id = label, y, angle) |>
    left_join(tip_df, by = "id")

  # Ring 1: CL0013 Family
  p <- p +
    new_scale_fill() +
    geom_tile(
      data = tip_coords,
      aes(x = M * 1.045, y = y, fill = cl0013_family),
      width = M * 0.055,
      height = 1,
      color = NA,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = PAL_CL0013_RING, drop = FALSE)

  # Ring 2: Viral PRE vs Bacterial SBL / PBP vs Non-PRE Outgroup
  p <- p +
    new_scale_fill() +
    geom_tile(
      data = tip_coords,
      aes(x = M * 1.105, y = y, fill = origin_func),
      width = M * 0.045,
      height = 1,
      color = NA,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      name   = "Origin / Function Ring",
      values = PAL_ORIGIN_FUNC,
      drop   = FALSE,
      guide  = guide_legend(order = 2, ncol = 2)
    )

  # Crisp outer family labels placed well outside Ring 2
  lbl_df <- compute_core_family_labels(tip_coords, active_fams) |>
    mutate(
      ang_norm = ((angle_mid %% 360) + 360) %% 360,
      flip = ang_norm > 90 & ang_norm < 270,
      text_ang = ifelse(flip, ang_norm + 180, ang_norm),
      hjust_val = ifelse(flip, 1, 0),
      col = unname(PAL_CL0013_RING[family])
    )

  p <- p +
    geom_text(
      data = lbl_df,
      aes(x = M * 1.175, y = y_mid, label = label, angle = text_ang, hjust = hjust_val, color = family),
      size = 3.1,
      fontface = "bold",
      lineheight = 0.88,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    xlim(0, M * 1.42) +
    labs(
      title = main_title,
      subtitle = main_sub
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 11.5, face = "bold", color = "#0F172A", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 8.4, color = "#475569", margin = margin(b = 4)),
      plot.margin = margin(6, 6, 4, 6),
      legend.position = "none"
    )

  # Build clean standalone legend panel separately
  leg_p <- ggtree(tr, layout = "fan", open.angle = 14, aes(color = branch_group), linewidth = 0.5) %<+% node_df +
    scale_color_manual(
      name = "CL0013 Family Branches & Inner Ring",
      values = PAL_CL0013_BRANCH,
      breaks = FAM_ORDER,
      labels = FAM_LABELS[FAM_ORDER],
      drop = FALSE,
      guide = guide_legend(order = 1, ncol = 4, override.aes = list(linewidth = 2.0))
    ) +
    new_scale_fill() +
    geom_tile(data = tip_coords, aes(x = M * 1.05, y = y, fill = origin_func), inherit.aes = FALSE) +
    scale_fill_manual(
      name = "Outer Origin / Function Ring",
      values = PAL_ORIGIN_FUNC,
      drop = FALSE,
      guide = guide_legend(order = 2, ncol = 4)
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 8.5, face = "bold", color = "#0F172A"),
      legend.text  = element_text(size = 7.8, color = "#334155"),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.y = unit(0.05, "cm"),
      legend.margin = margin(2, 4, 2, 4)
    )
  shared_leg <- cowplot::get_legend(leg_p)
  list(plot = p, legend = shared_leg)
}

build_cl0013_loo_heatmap <- function() {
  df <- read_tsv(file.path(HGT_DIR, "CL0013_Superfamily_Rooting_Validation.tsv"), show_col_types = FALSE)
  tree_lbls <- c(
    "CL0013.struct_c20.no_PF27310"     = "FoldMason+C20\nPF04960+PF06089",
    "CL0013.struct_c20.only_PF04960"   = "FoldMason+C20\nOnly PF04960",
    "CL0013.struct_c20.only_PF06089"   = "FoldMason+C20\nOnly PF06089",
    "CL0013.all3_out.mafft_out"        = "MAFFT Baseline\nAll 3 (+PF27310)",
    "CL0013.loo_no_PF04960.mafft_out"  = "MAFFT Baseline\n-PF04960",
    "CL0013.loo_no_PF06089.mafft_out"  = "MAFFT Baseline\n-PF06089",
    "CL0013.loo_no_PF27310.mafft_out"  = "MAFFT Baseline\n-PF27310"
  )

  # Reshape Split F1 + 5 PRE family F1/supports into clean matrix
  rows <- list()
  present_lbls <- character()
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    t_lab <- ifelse(r$cl0013_tree %in% names(tree_lbls), tree_lbls[r$cl0013_tree], r$cl0013_tree)
    present_lbls <- c(present_lbls, t_lab)
    rows[[length(rows) + 1]] <- tibble(
      tree = t_lab, metric = "Non-PRE\nOutgroup Split",
      f1 = r$nonpre_outgroup_split_F1,
      cell_txt = sprintf("F1 = %.3f\n(%d out)", r$nonpre_outgroup_split_F1, r$n_outgroup_nonpre)
    )
    for (fam in c("PF00144", "PF00905", "PF13354", "PF00768", "PF02113")) {
      f1_v <- as.numeric(r[[paste0(fam, "_clade_F1")]])
      sh_v <- as.numeric(r[[paste0(fam, "_SH_aLRT")]])
      uf_v <- as.numeric(r[[paste0(fam, "_UFBoot")]])
      rows[[length(rows) + 1]] <- tibble(
        tree = t_lab, metric = fam, f1 = f1_v,
        cell_txt = sprintf("F1 = %.3f\nSH %.0f / UF %.0f", f1_v, sh_v, uf_v)
      )
    }
  }
  f1_floor <- min(0.65, floor(min(vapply(rows, function(x) x$f1, numeric(1)), na.rm = TRUE) * 20) / 20)
  mat <- bind_rows(rows) |>
    mutate(
      tree   = factor(tree, levels = rev(unique(present_lbls))),
      metric = factor(metric, levels = c("Non-PRE\nOutgroup Split", "PF00144", "PF00905", "PF13354", "PF00768", "PF02113")),
      txt_col = ifelse(f1 >= 0.90, "#FFFFFF", "#0F172A")
    )

  ggplot(mat, aes(x = metric, y = tree, fill = f1)) +
    geom_tile(color = "#FFFFFF", linewidth = 0.9) +
    geom_text(aes(label = cell_txt, color = txt_col), size = 2.55, fontface = "bold", lineheight = 0.90) +
    scale_color_identity() +
    scale_fill_gradientn(
      colours = c("#FDE68A", "#38BDF8", "#0284C7", "#0F766E"),
      values = rescale(c(f1_floor, 0.85, 0.95, 1.00)),
      limits = c(f1_floor, 1.00),
      name = "Clade F1"
    ) +
    labs(
      title = "B   CL0013 FoldMason+C20 vs. Baseline LOO Monophyly",
      subtitle = "Removing short PF27310 (~200 aa) & adding FoldMason+hmmalign C20 across all PRE families",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 11.5, face = "bold", color = "#0F172A", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 8.0, color = "#475569", margin = margin(b = 6)),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 8.4, face = "bold", color = "#1E293B"),
      axis.text.y = element_text(size = 8.4, face = "bold", color = "#1E293B"),
      legend.position = "right",
      legend.key.height = unit(0.85, "cm"),
      legend.key.width  = unit(0.32, "cm"),
      legend.title = element_text(size = 8.0, face = "bold"),
      plot.margin = margin(6, 8, 6, 4)
    )
}

build_cl0013_mini_loo_trees <- function(meta_df) {
  val_df <- read_tsv(file.path(HGT_DIR, "CL0013_Superfamily_Rooting_Validation.tsv"), show_col_types = FALSE)
  get_split_sub <- function(tname) {
    v <- val_df$nonpre_outgroup_split_F1[val_df$cl0013_tree == tname]
    if (length(v) == 0 || is.na(v[1])) return("Non-PRE Split F1 = NA")
    sprintf("Non-PRE Split F1 = %.3f", as.numeric(v[1]))
  }

  has_struct_loo <- file.exists(file.path(CL_DIR, "CL0013.struct_c20.only_PF04960.treefile")) &&
                    file.exists(file.path(CL_DIR, "CL0013.struct_c20.only_PF06089.treefile")) &&
                    file.exists(file.path(CL_DIR, "outgroup_struct_only_PF04960.txt")) &&
                    file.exists(file.path(CL_DIR, "outgroup_struct_only_PF06089.txt"))

  if (has_struct_loo) {
    specs <- list(
      list(
        tfile = file.path(CL_DIR, "CL0013.struct_c20.only_PF04960.treefile"),
        ogfile = file.path(CL_DIR, "outgroup_struct_only_PF04960.txt"),
        stem = "struct_only_PF04960",
        title = "C1  FoldMason+C20: Only PF04960",
        sub = get_split_sub("CL0013.struct_c20.only_PF04960")
      ),
      list(
        tfile = file.path(CL_DIR, "CL0013.struct_c20.only_PF06089.treefile"),
        ogfile = file.path(CL_DIR, "outgroup_struct_only_PF06089.txt"),
        stem = "struct_only_PF06089",
        title = "C2  FoldMason+C20: Only PF06089",
        sub = get_split_sub("CL0013.struct_c20.only_PF06089")
      ),
      list(
        tfile = file.path(CL_DIR, "CL0013.all3_out.mafft_out.treefile"),
        ogfile = file.path(CL_DIR, "outgroup_taxa_all3.txt"),
        stem = "baseline_all3",
        title = "C3  MAFFT Baseline (+PF27310)",
        sub = get_split_sub("CL0013.all3_out.mafft_out")
      )
    )
    hdr_txt <- "C   Single-Outgroup FoldMason+C20 LOO vs. Baseline (+PF27310) Topologies"
  } else {
    specs <- list(
      list(
        tfile = file.path(CL_DIR, "CL0013.loo_no_PF04960.mafft_out.treefile"),
        ogfile = file.path(CL_DIR, "outgroup_taxa_loo_no_PF04960.txt"),
        stem = "loo_no_PF04960",
        title = "C1  LOO: -PF04960",
        sub = get_split_sub("CL0013.loo_no_PF04960.mafft_out")
      ),
      list(
        tfile = file.path(CL_DIR, "CL0013.loo_no_PF06089.mafft_out.treefile"),
        ogfile = file.path(CL_DIR, "outgroup_taxa_loo_no_PF06089.txt"),
        stem = "loo_no_PF06089",
        title = "C2  LOO: -PF06089",
        sub = get_split_sub("CL0013.loo_no_PF06089.mafft_out")
      ),
      list(
        tfile = file.path(CL_DIR, "CL0013.loo_no_PF27310.mafft_out.treefile"),
        ogfile = file.path(CL_DIR, "outgroup_taxa_loo_no_PF27310.txt"),
        stem = "loo_no_PF27310",
        title = "C3  LOO: -PF27310",
        sub = get_split_sub("CL0013.loo_no_PF27310.mafft_out")
      )
    )
    hdr_txt <- "C   Leave-One-Outgroup-Out Tree Topologies (Identical 8-Family Coloring)"
  }

  plots <- lapply(specs, function(sp) {
    prep <- prepare_cl0013_tree(sp$tfile, sp$ogfile, meta_df)
    p <- ggtree(prep$tree, layout = "fan", open.angle = 12, aes(color = branch_group), linewidth = 0.28) %<+% prep$node_df +
      scale_color_manual(values = PAL_CL0013_BRANCH, drop = FALSE, guide = "none")
    M <- max(p$data$x, na.rm = TRUE)
    tip_c <- p$data |>
      filter(isTip == TRUE) |>
      select(id = label, y) |>
      left_join(prep$tip_df |> rename(id = uid), by = "id")
    is_warn <- grepl("NA$", sp$sub) || (!grepl("F1 = 1\\.000", sp$sub) && sp$stem == "loo_no_PF06089")
    p +
      geom_tile(
        data = tip_c,
        aes(x = M * 1.05, y = y, fill = cl0013_family),
        width = M * 0.08, height = 1, color = NA, inherit.aes = FALSE, show.legend = FALSE
      ) +
      scale_fill_manual(values = PAL_CL0013_RING, drop = FALSE) +
      xlim(0, M * 1.15) +
      labs(title = sp$title, subtitle = sp$sub) +
      theme_void() +
      theme(
        plot.title    = element_text(size = 9.2, face = "bold", color = "#0F172A", hjust = 0.5, margin = margin(b = 1)),
        plot.subtitle = element_text(size = 7.6, color = ifelse(is_warn, "#B45309", "#0F766E"),
                                     face = "bold", hjust = 0.5, margin = margin(b = 2)),
        plot.margin   = margin(2, 2, 2, 2)
      )
  })
  row_p <- cowplot::plot_grid(plotlist = plots, nrow = 1)
  hdr <- ggdraw() +
    draw_label(
      hdr_txt,
      x = 0.01, hjust = 0, fontface = "bold", size = 11, color = "#0F172A"
    )
  cowplot::plot_grid(hdr, row_p, ncol = 1, rel_heights = c(0.11, 0.89))
}

render_cl0013_figure <- function() {
  message("[Figure 1] Rendering Fig_CL0013_Superfamily_and_LOO ...")
  meta_df <- read_tsv(file.path(CL_DIR, "cl0013_taxa_metadata.tsv"), show_col_types = FALSE)
  main_obj <- build_cl0013_main_tree(meta_df)
  loo_heat <- build_cl0013_loo_heatmap()
  loo_mini <- build_cl0013_mini_loo_trees(meta_df)

  right_col <- cowplot::plot_grid(loo_heat, loo_mini, ncol = 1, rel_heights = c(0.49, 0.51))
  body_grid <- cowplot::plot_grid(main_obj$plot, right_col, nrow = 1, rel_widths = c(1.08, 1.0))
  full_fig  <- cowplot::plot_grid(body_grid, main_obj$legend, ncol = 1, rel_heights = c(1.0, 0.13))
  save_multi_format(full_fig, "Fig_CL0013_Superfamily_and_LOO", width = 16.5, height = 9.6)
}

# ==============================================================================
# 4. Figure 2: Stage 3 Stratified HGT & OAT Robustness Sensitivity Matrix
# ==============================================================================
FAM_FIVE <- c("PF00144", "PF00905", "PF13354", "PF00768", "PF13354_PF00768")
FAM_FIVE_LBL <- c(
  "PF00144"         = "PF00144\n(Class C-like)",
  "PF00905"         = "PF00905\n(Transpeptidase)",
  "PF13354"         = "PF13354\n(Class A-like)",
  "PF00768"         = "PF00768\n(Peptidase_S11)",
  "PF13354_PF00768" = "Joint\nPF13354+PF00768"
)

build_hgt_stratified_bars <- function() {
  df <- read_tsv(file.path(HGT_DIR, "HGT_Summary_By_Family_And_Function_2026.tsv"), show_col_types = FALSE) |>
    mutate(
      family_lbl = factor(FAM_FIVE_LBL[family], levels = rev(FAM_FIVE_LBL[FAM_FIVE])),
      dir_lbl = ifelse(direction == "Bacterial_to_Viral", "Bacterial -> Viral", "Viral -> Bacterial"),
      func_short = ifelse(functional_category == "Putative_Hydrolytic_BetaLactamase", "Hydrolytic SBL", "PBP / PRE Homolog"),
      stratum = factor(
        paste(dir_lbl, func_short, sep = " | "),
        levels = c(
          "Bacterial -> Viral | Hydrolytic SBL",
          "Bacterial -> Viral | PBP / PRE Homolog",
          "Viral -> Bacterial | Hydrolytic SBL",
          "Viral -> Bacterial | PBP / PRE Homolog"
        )
      )
    )

  long <- df |>
    select(family_lbl, stratum, dir_lbl, func_short, tier1_events, tier2_events, tier1_recip_nr85_tips) |>
    pivot_longer(cols = c("tier1_events", "tier2_events"), names_to = "tier", values_to = "count") |>
    mutate(
      tier_lbl = factor(ifelse(tier == "tier1_events", "Tier 1 (UFBoot>=95 & SH>=80)", "Tier 2 (Supported)"),
                        levels = c("Tier 2 (Supported)", "Tier 1 (UFBoot>=95 & SH>=80)"))
    )

  # Direct end-of-bar summary annotations
  bar_tot <- df |>
    mutate(tot = tier1_events + tier2_events) |>
    filter(tot > 0) |>
    mutate(
      lbl = ifelse(
        tier1_recip_nr85_tips > 0,
        sprintf("T1:%d | T1+2:%d (%d recipient tips)", tier1_events, tot, tier1_recip_nr85_tips),
        sprintf("T1:%d | T1+2:%d", tier1_events, tot)
      )
    )

  pal_tier <- c(
    "Tier 1 (UFBoot>=95 & SH>=80)" = "#BE123C",
    "Tier 2 (Supported)"           = "#FDBA74"
  )

  ggplot(long, aes(x = count, y = stratum, fill = tier_lbl)) +
    geom_col(width = 0.68, color = "#1E293B", linewidth = 0.25) +
    geom_text(
      data = bar_tot,
      aes(x = tot + 0.35, y = stratum, label = lbl),
      hjust = 0, size = 2.6, fontface = "bold", color = "#334155", inherit.aes = FALSE
    ) +
    facet_grid(family_lbl ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = pal_tier, name = "Confidence Tier") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.36))) +
    labs(
      title = "A   Primary HGT Counts by Family, Direction & Functional Tier",
      subtitle = "Tier 1 vs. Tier 2 nested Fitch+PastML events (Prophage-borne = 0 across all primary trees)",
      x = "Number of Independent HGT Clades",
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 11, face = "bold", color = "#0F172A", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 8.0, color = "#475569", margin = margin(b = 4)),
      strip.text.y = element_text(size = 8.0, face = "bold", color = "#0F172A", angle = 0),
      strip.background = element_rect(fill = "#F1F5F9", color = "#CBD5E1"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 7.6, color = "#1E293B"),
      legend.position = "bottom",
      legend.title = element_text(size = 8.2, face = "bold"),
      legend.text = element_text(size = 7.8),
      plot.margin = margin(6, 10, 4, 6)
    )
}

build_event_robustness_scatter <- function() {
  ev <- read_tsv(file.path(HGT_DIR, "Primary_HGT_Events_Catalog_2026.tsv"), show_col_types = FALSE) |>
    filter(as_bool(is_tier1_hgt) | as_bool(is_tier2_hgt)) |>
    mutate(
      oat_pct = oat_recurrence_rate * 100.0,
      dir_lbl = ifelse(direction == "Bacterial_to_Viral", "Bacterial -> Viral", "Viral -> Bacterial"),
      func_lbl = ifelse(functional_category == "Putative_Hydrolytic_BetaLactamase", "Hydrolytic SBL", "PBP / PRE Homolog"),
      is_t1 = as_bool(is_tier1_hgt),
      is_ri = as_bool(root_Invariant),
      size_seqs = ifelse(
        direction == "Bacterial_to_Viral",
        pmax(coalesce(expanded_viral_cluster_seqs, 0), n_focal_recipient_tips, 1),
        pmax(n_focal_recipient_tips, 1)
      ),
      short_id = paste0(
        ifelse(family == "PF13354_PF00768", "Joint", family),
        " ", sub(".*__(HGT\\d+)$", "\\1", event_id),
        ifelse(size_seqs > 1, sprintf(" (n=%d)", size_seqs), "")
      )
    )

  t1_ev <- ev |> filter(is_t1)

  # Subtle deterministic jitter so overlapping integer (OAT%, UFBoot%) points remain countable
  set.seed(2026)
  ev <- ev |>
    mutate(
      x_jit = pmin(100.6, pmax(-0.6, oat_pct + runif(n(), -1.4, 1.4))),
      y_jit = pmin(100.4, min_PG_ufboot + runif(n(), -0.45, 0.45))
    )
  t1_lbl_df <- ev |> filter(is_t1)

  y_min <- min(40, floor(min(ev$min_PG_ufboot, na.rm = TRUE) / 5) * 5)

  ggplot(ev, aes(x = x_jit, y = y_jit)) +
    annotate("rect", xmin = 80, xmax = 102.5, ymin = 95, ymax = 108.5, fill = "#FEF2F2", alpha = 0.85) +
    annotate("rect", xmin = 80, xmax = 102.5, ymin = y_min - 1, ymax = 95, fill = "#F8FAFC", alpha = 0.85) +
    geom_vline(xintercept = 80, linetype = "dashed", color = "#DC2626", linewidth = 0.55) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "#334155", linewidth = 0.55) +
    annotate("text", x = 80.8, y = 108.0, label = "Tier 1 + Robust Core", hjust = 0, vjust = 1,
             size = 2.75, fontface = "bold", color = "#B91C1C") +
    annotate("text", x = 80.8, y = y_min + 1.0, label = "OAT >= 80% Core", hjust = 0, vjust = 0,
             size = 2.6, fontface = "italic", color = "#475569") +
    # Root-invariant outer halo
    geom_point(
      data = filter(ev, is_ri),
      aes(size = log1p(size_seqs)),
      shape = 21, fill = NA, color = "#0F172A", stroke = 0.9, alpha = 0.55
    ) +
    geom_point(
      aes(color = dir_lbl, shape = func_lbl, size = log1p(size_seqs)),
      alpha = 0.86
    ) +
    ggrepel::geom_text_repel(
      data = t1_lbl_df,
      aes(label = short_id, color = dir_lbl),
      size = 2.35, fontface = "bold", force = 2.8, max.overlaps = Inf,
      box.padding = 0.35, point.padding = 0.2,
      xlim = c(-2, 102), ylim = c(85, 106.2),
      min.segment.length = 0.1, segment.size = 0.25, show.legend = FALSE
    ) +
    scale_color_manual(values = PAL_DIR, name = "HGT Direction") +
    scale_shape_manual(values = c("Hydrolytic SBL" = 17, "PBP / PRE Homolog" = 16), name = "Function") +
    scale_size_continuous(range = c(2.0, 6.2), guide = "none") +
    scale_x_continuous(limits = c(-3, 103), breaks = seq(0, 100, 20)) +
    scale_y_continuous(limits = c(y_min - 1, 108.8), breaks = seq(y_min, 100, 10)) +
    labs(
      title = sprintf("B   Primary HGT Robustness Landscape (n = %d Tier 1 + Tier 2 Events)", nrow(ev)),
      subtitle = sprintf("%d Tier 1 | %d Robust Core (OAT >= 80%%) | %d Root-Invariant (dark outline); point size ~ viral seqs",
                         nrow(t1_ev), sum(as_bool(ev$robust_core_hgt)), sum(ev$is_ri)),
      x = "OAT Recurrence Rate Across Perturbations (%)",
      y = "min(Parent, Grandparent) UFBoot (%)"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 11, face = "bold", color = "#0F172A", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.9, color = "#475569", margin = margin(b = 4)),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 8.2, face = "bold"),
      legend.text = element_text(size = 7.8),
      plot.margin = margin(6, 10, 4, 6)
    )
}

build_oat_heatmap <- function() {
  mat <- read_tsv(file.path(HGT_DIR, "OAT_Sensitivity_Summary_Matrix_2026.tsv"), show_col_types = FALSE)

  oat_order <- c(
    "primary_out", "primary_in_MAD", "primary_in_MinVar", "primary_in_Midpoint"
  )
  if (any(mat$oat_key == "clean_prof_out")) {
    oat_order <- c(oat_order, "clean_prof_out")
  }
  if (any(mat$oat_key == "legacy_sarg24_out")) {
    oat_order <- c(oat_order, "legacy_sarg24_out")
  }
  oat_order <- c(
    oat_order,
    "primary_in_YoudenJ", "mafft_out", "pmsf_c20",
    "rep2", "rep3", "rep4", "rep5", "phylo", "nopro",
    "cov60", "cov70"
  )

  prim_method <- unique(na.omit(mat$rooting_method[mat$oat_key == "primary_out"]))
  prim_label <- paste0("Primary\n(", if (length(prim_method)) paste(gsub("_", " ", prim_method), collapse = "/") else "primary root", ")")
  oat_labels <- c(
    "primary_out"         = prim_label,
    "primary_in_MAD"      = "Ingroup\nMAD",
    "primary_in_MinVar"   = "Ingroup\nMinVar",
    "primary_in_Midpoint" = "Ingroup\nMidpoint",
    "clean_prof_out"      = "Clean GTDB\n(--addprofile)",
    "legacy_sarg24_out"   = "Legacy 24-OG\n(SARG-mixed)",
    "primary_in_YoudenJ"  = "Youden J\n(App.)",
    "mafft_out"           = "MAFFT +\nClipKIT",
    "pmsf_c20"            = "PMSF\n(LG+C20)",
    "rep2"                = "Sample\nrep2",
    "rep3"                = "Sample\nrep3",
    "rep4"                = "Sample\nrep4",
    "rep5"                = "Sample\nrep5",
    "phylo"               = "Phylo-65%\nrep",
    "nopro"               = "Without\nProphage",
    "cov60"               = "HMM Cov\n>= 60%",
    "cov70"               = "HMM Cov\n>= 70%"
  )

  n_root_cols <- sum(oat_order %in% c(
    "primary_out", "primary_in_MAD", "primary_in_MinVar", "primary_in_Midpoint",
    "clean_prof_out", "legacy_sarg24_out", "primary_in_YoudenJ"
  ))

  fam_lbl_heat <- c(
    "PF00144"         = "PF00144 (Class C)",
    "PF00905"         = "PF00905 (Transpep.)",
    "PF13354"         = "PF13354 (Class A)",
    "PF00768"         = "PF00768 (Pep. S11)",
    "PF13354_PF00768" = "Joint PF13354+00768"
  )

  # Build complete 5 x N grid so intentionally skipped cov60/cov70 on PF00768 & Joint
  # and any pending OAT jobs show explicitly rather than blank gaps
  full_grid <- expand_grid(
    family  = FAM_FIVE,
    oat_key = oat_order
  ) |>
    left_join(mat, by = c("family", "oat_key")) |>
    mutate(
      is_pending = ((oat_key == "pmsf_c20" & family %in% PENDING_PMSF_FAMS) |
                    (oat_key == "clean_prof_out" & is.na(tier1_hgt_events))),
      tot_t12 = ifelse(is_pending | is.na(tier1_hgt_events), NA_real_, tier1_hgt_events + tier2_hgt_events),
      fill_val = ifelse(is.na(tot_t12), NA_real_, log10(tot_t12 + 1)),
      tile_txt = case_when(
        is_pending                  ~ "Pending\n(HPC)",
        is.na(tier1_hgt_events)     ~ "N/A",
        TRUE                        ~ sprintf("T1: %d (%d/%d)\nT1+2: %d",
                                              tier1_hgt_events, tier1_B_to_V_events, tier1_V_to_B_events, tot_t12)
      ),
      txt_col = case_when(
        is_pending | is.na(tier1_hgt_events) ~ "#475569",
        tot_t12 >= 18                        ~ "#FFFFFF",
        TRUE                                 ~ "#0F172A"
      ),
      fam_fac = factor(fam_lbl_heat[family], levels = rev(fam_lbl_heat[FAM_FIVE])),
      oat_fac = factor(oat_labels[oat_key], levels = oat_labels[oat_order])
    )

  ggplot(full_grid, aes(x = oat_fac, y = fam_fac, fill = fill_val)) +
    geom_tile(color = "#FFFFFF", linewidth = 0.8) +
    geom_text(aes(label = tile_txt, color = txt_col), size = 2.40, fontface = "bold", lineheight = 0.92) +
    # Vertical group separators: Rooting | MSA+PMSF | Sampling+Prophage | Coverage
    geom_vline(xintercept = c(n_root_cols + 0.5, n_root_cols + 2.5, n_root_cols + 8.5),
               color = "#334155", linewidth = 0.65, linetype = "solid") +
    scale_color_identity() +
    scale_fill_gradientn(
      colours = c("#F8FAFC", "#BAE6FD", "#38BDF8", "#0284C7", "#1E3A8A"),
      na.value = "#E2E8F0",
      breaks = log10(c(0, 3, 8, 15, 30) + 1),
      labels = c("0", "3", "8", "15", "30"),
      name = "Total HGT\n(T1 + T2)"
    ) +
    labs(
      title = sprintf("C   One-Factor-at-a-Time (OAT) & Multi-Rooting Sensitivity Matrix (%d Completed Trees)", nrow(mat)),
      subtitle = "Tier 1 (B->V / V->B) over total Tier 1 + 2 HGT; primary root (crown projection where supported, otherwise MAD) vs. MinVar/Midpoint, clean/legacy outgroups & LG+C20 PMSF",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 11, face = "bold", color = "#0F172A", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 8.0, color = "#475569", margin = margin(b = 6)),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 7.7, face = "bold", color = "#1E293B"),
      axis.text.y = element_text(size = 8.4, face = "bold", color = "#0F172A"),
      legend.position = "right",
      legend.key.height = unit(0.95, "cm"),
      legend.key.width  = unit(0.32, "cm"),
      legend.title = element_text(size = 7.8, face = "bold"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

render_hgt_figure <- function() {
  message("[Figure 2] Rendering Fig_Phylogenetic_HGT_and_OAT_Robustness ...")
  pa <- build_hgt_stratified_bars()
  pb <- build_event_robustness_scatter()
  pc <- build_oat_heatmap()

  top_row <- cowplot::plot_grid(pa, pb, nrow = 1, rel_widths = c(1.05, 1.0))
  full_fig <- cowplot::plot_grid(top_row, pc, ncol = 1, rel_heights = c(1.15, 0.85))
  save_multi_format(full_fig, "Fig_Phylogenetic_HGT_and_OAT_Robustness", width = 16.8, height = 11.4)
}

main <- function() {
  render_cl0013_figure()
  render_hgt_figure()
  message("[Done] Both figures written to: ", OUT_DIR)
}

main()
