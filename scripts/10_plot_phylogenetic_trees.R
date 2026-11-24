#!/usr/bin/env Rscript
# ==============================================================================
# Circular phylogenetic tree figures for the primary IQ-TREE trees (PDF + LZW-TIFF + PNG)
# 主树环形系统发育树作图（PDF + LZW-TIFF + PNG 输出）
#
# Usage / 运行方式:
#   1) 命令行: Rscript scripts/10_plot_phylogenetic_trees.R
#   2) RStudio: 直接打开本脚本并 Source / Ctrl+Alt+R 运行（自动识别仓库根目录）
#
# Rooting / 定根: Sister-PRE outgroup (pre-registered primary rooting); the outgroup tips are
#   removed after rooting. Set ROOTING <- "midpoint" for the midpoint sensitivity view.
# Branch lengths beyond the CAP_DEPTH_QUANTILE tip depth are shortened for display only; the
#   subtitle states when this is applied.
#
# Inputs / 输入文件:
#   - results/reanalysis_2026/primary_trees/*.treefile, *.iqtree & *_outgroup.txt
#   - results/reanalysis_2026/Final_Viral_SBL_Catalog.tsv
#   - results/reanalysis_2026/Final_Bacterial_TripleSource_Catalog.tsv
#   - results/reanalysis_2026/Phase2_Viral_Genomic_Neighborhoods_Full.tsv
#
# Outputs / 输出文件:
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
# SECTION 1: 可自定义参数与配色区（后续微调出刊图时只需修改这里）
# ==============================================================================

# 1.1 自动识别项目根目录
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
  stop("无法定位 SBL_research 根目录，请在仓库根目录或 scripts/ 目录下运行。")
}

REPO_ROOT <- locate_repo_root()
TREE_DIR  <- file.path(REPO_ROOT, "results", "reanalysis_2026", "primary_trees")
META_DIR  <- file.path(REPO_ROOT, "results", "reanalysis_2026")
OUT_DIR   <- file.path(REPO_ROOT, "figures", "phylogenetic_trees")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 1.2 图形版式与分辨率控制
ROOTING            <- "outgroup"  # "outgroup"（Sister-PRE 外群定根，主分析）或 "midpoint"（敏感性）
OPEN_ANGLE         <- 16      # 环形缺口角度（度），用于放置 R1/R2/R3 环标签
CAP_DEPTH_QUANTILE <- 0.985   # 截断前 1.5% 极端长枝以展开内群拓扑（设为 1.0 则不截断）
SHOW_BOOT_NODES    <- TRUE    # 是否标出高支持度节点（UFBoot >= 95 & SH-aLRT >= 80）
SHOW_SBL_LABELS    <- TRUE    # 是否在最外圈标注代表性病毒临床同源 SBL 基因名
LABEL_ENV_SBL      <- FALSE   # 是否重复标出 "Motif-SBL" 文字（默认仅用紫色菱形 ◆ 表示）
MIN_LABEL_TIP_GAP  <- 18      # 防止外圈基因标签重叠的最小叶节点间隔数
FIG_WIDTH          <- 13.0    # 单树图宽（英寸）
FIG_HEIGHT         <- 10.2    # 单树图高（英寸）
TIFF_DPI           <- 300     # 出刊定稿可改为 600（当前 300 dpi 兼顾清晰度与文件体积）
EXPORT_PNG_PREVIEW <- TRUE    # 同时输出 PNG 方便快速浏览

# 1.3 分支配色（Branches）
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

# 1.4 外圈 Ring 1：来源与病毒可信度分层（Origin & Viral Tier）
PAL_RING1 <- c(
  "Viral Tier 1 (High-conf.)"      = "#DC2626",
  "Viral Tier 2 (Supported)"       = "#F97316",
  "Viral Tier 3 (Fragment)"        = "#FACC15",
  "Viral MOSAiC (Arctic vMAG)"     = "#9333EA",
  "Bacterial Prophage"             = "#0284C7",
  "Bacterial Reference"            = "#E2E8F0"
)

# 1.5 外圈 Ring 2：催化功能分层（CARD 同源 SBL vs 仅基序定义 SBL vs PBP/PRE）
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

# 1.6 外圈 Ring 3：家族特征结构域架构（单家族树）与家族归属（联立树）
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
# SECTION 2: 读取并构建统一规范化注释元数据表（R 中可直接复用）
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
  # 防止重复 uid（若某些条目同时存在于单家族目录）
  merged <- distinct(merged, uid, .keep_all = TRUE)

  out_tsv <- file.path(TREE_DIR, "Tree_Tip_Annotations_Merged.tsv")
  write_tsv(merged, out_tsv)
  message("  -> Saved unified tip metadata table: ", out_tsv, " (n=", nrow(merged), ")")
  merged
}

# ==============================================================================
# SECTION 3: 定根、极端长枝截断与枝状态注释
# ==============================================================================

prepare_tree_and_nodes <- function(tree_path, outgroup_path, meta_df) {
  tr <- ape::read.tree(tree_path)
  og <- read_lines(outgroup_path)
  og <- trimws(og[nzchar(trimws(og))])

  drop_ids <- intersect(tr$tip.label, og)
  root_note <- "midpoint-rooted"
  if (ROOTING == "outgroup" && length(drop_ids) > 0) {
    rooted <- root_on_outgroup(tr, drop_ids)
    tr <- rooted$tree
    root_note <- rooted$note
  }
  if (length(drop_ids) > 0) {
    tr <- ape::drop.tip(tr, drop_ids)
  }
  if (ROOTING != "outgroup" || length(drop_ids) == 0) {
    tr <- phangorn::midpoint(tr)
  }
  tr <- ape::ladderize(tr, right = FALSE)

  # 截断极端前 1.5% 离群长枝（仅用于显示；保持父节点深度 <= 子节点深度，保证非负枝长）
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

  # 解析 IQ-TREE SH-aLRT / UFBoot 支持度
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

  # 自底向上快速累计每个分支后代的病毒数与真水解型 SBL 数
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

  # tip_df 仅保留外圈 geom_fruit 与标注所需列，避免与 p$data / node_df 列名冲突
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

# Best-fit model from the IQ-TREE report next to the tree file (falls back to spec$model)
read_iqtree_model <- function(tree_path, fallback) {
  rep <- sub("\\.treefile$", ".iqtree", tree_path)
  if (!file.exists(rep)) return(paste(fallback, "(model not verified)"))
  hit <- grep("Best-fit model according to BIC:", readLines(rep, warn = FALSE), value = TRUE)
  if (length(hit) == 0) return(paste(fallback, "(model not verified)"))
  trimws(sub(".*BIC:", "", hit[1]))
}

# 选择外圈病毒 SBL 标签（保证按树的叶顺序至少间隔 MIN_LABEL_TIP_GAP，避免重叠）
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
# SECTION 4: ggtree + ggtreeExtra 绘图核心函数
# ==============================================================================

build_ggtree_plot <- function(spec, meta_df, compact_mode = FALSE) {
  tree_path <- file.path(TREE_DIR, spec$treefile)
  og_path   <- file.path(TREE_DIR, spec$outgroup)
  prep <- prepare_tree_and_nodes(tree_path, og_path, meta_df)
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

  # Ring 3: Domain Architecture (单家族) OR Pfam Family (联立树)
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

  # 最外圈叠加：病毒水解型 SBL 外标记点（CARD 同源 vs 仅基序定义）
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
        values = setNames(c(8, 18), c(LBL_V_CARD, LBL_V_MOTIF)),   # 星号 * / 菱形 ◆
        drop = TRUE,
        guide = guide_legend(order = 5)
      )
  }

  # 外圈稀疏基因标签（避免互相重叠）
  if (SHOW_SBL_LABELS && !compact_mode) {
    lbl_df <- select_spaced_sbl_labels(tip_coords, min_gap = max(14, round(n_tot / 95)))
    if (nrow(lbl_df) > 0) {
      # 在 ggtree fan 布局中，angle 位于 [0, 360]；左半圆 [90, 270] 自动翻转 180° 以保证正向阅读
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

  # 在环形缺口内标出 R1 / R2 / R3 环序号
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

  # 扩大径向边界避免最外圈文字被裁剪
  x_lim_hi <- if (SHOW_SBL_LABELS && !compact_mode) max_x * 1.50 else max_x * 1.35
  sub_txt <- sprintf(
    "Ingroup n=%s (Viral nr85=%s, Bacterial=%s) | Viral SBL: CARD-homolog=%d, motif-defined=%d | Bac SBL=%d | %s",
    format(n_tot, big.mark = ","), format(n_v, big.mark = ","), format(n_b, big.mark = ","),
    n_vsbl_clin, n_vsbl_env, n_bsbl,
    paste(read_iqtree_model(tree_path, spec$model), spec$track, sep = " | ")
  )
  sub_txt <- paste0(
    sub_txt, "
", prep$root_note, "; outgroup tips removed",
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
# SECTION 5: 批量生成并导出 PDF（矢量）+ LZW-TIFF（投稿位图）+ PNG 预览
# ==============================================================================

save_multi_format <- function(plot_obj, stem, width, height, dpi = TIFF_DPI) {
  pdf_path  <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  tiff_path <- file.path(OUT_DIR, paste0(stem, ".tiff"))
  png_path  <- file.path(OUT_DIR, paste0(stem, ".png"))

  # PDF（优先使用 cairo_pdf，若系统缺少 Cairo 则回退标准 pdf）
  tryCatch(
    ggsave(pdf_path, plot = plot_obj, device = grDevices::cairo_pdf, width = width, height = height, units = "in"),
    error = function(e) ggsave(pdf_path, plot = plot_obj, device = "pdf", width = width, height = height, units = "in")
  )

  # LZW 压缩 TIFF（期刊投稿标准）
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

  # 生成 2x2 四图全景对比图（PF00144, PF00905, PF13354, PF13354_PF00768）
  message("[Plotting] Overview_4Panel_trees ...")
  ov_specs <- specs[c(1, 2, 3, 5)]
  ov_plots <- lapply(ov_specs, function(sp) build_ggtree_plot(sp, meta_df, compact_mode = TRUE))
  p_for_leg <- build_ggtree_plot(specs[[1]], meta_df, compact_mode = FALSE) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal")
  shared_leg <- cowplot::get_legend(p_for_leg)

  # 额外补充 Panel D 的 Ring 3 (PF13354 vs PF00768) 家族配色图例
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
