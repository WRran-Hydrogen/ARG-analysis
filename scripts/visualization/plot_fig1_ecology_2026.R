#!/usr/bin/env Rscript
# ==============================================================================
# Figure 1 - global distribution and habitat associations of viral SBLs and homologs
# Main-text families PF00144, PF00905 and PF13354 (PF00768 is reported in the supplement).
#
#   A  Sites of SBL-carrying samples (IMG/VR v4.1 + MOSAiC), shape = Pfam family, colour = habitat,
#      size = number of positive samples at the site; a black outline marks a verified hydrolytic SBL
#   B  log2 habitat enrichment of Uroviricota UViGs per protein-architecture state (all SBL carriers on top) against
#      the IMG/VR high-confidence background; stars need significance also with one UViG per metagenome
#   C  Protein architecture of the curated viral SBLs: share of six architecture states per family (the states of the
#      ancestral state reconstruction) and representative domain architectures drawn to scale
# Supplementary figure (figures/supplementary_2026/FigS_Climate_Zone_Prevalence_and_Habitat_Flows):
#   A  Prevalence by climate zone (positive sites / total sites, Wilson 95 % CI, BH-adjusted Fisher tests)
#   B  Habitat -> viral phylum -> host phylum flows of SBL-carrying UViGs (minor categories pooled)
#
# Inputs (results/reanalysis_2026/06_structure_and_ecology/):
#   Sample_Geography_2026.tsv, Climate_Zone_Prevalence_2026.tsv,
#   Habitat_Enrichment_By_Architecture_State_2026.tsv, Habitat_ViralPhylum_Host_Flows_2026.tsv, Viral_SBL_Architecture_Summary_2026.tsv,
#   Viral_SBL_Architecture_Classes_2026.tsv; seq_len from 01_curation_and_context/Final_Viral_SBL_Catalog.tsv
# Output: figures/main_2026/Fig1_Global_Distribution_and_Habitat.{pdf,png,tiff}
#
# Usage: Rscript scripts/visualization/plot_fig1_ecology_2026.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(sf)
  library(maps); library(ggforce); library(patchwork); library(scales)
})
sf_use_s2(FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
BASE <- normalizePath(file.path(script_dir, "..", ".."))
ECO <- file.path(BASE, "results", "reanalysis_2026", "06_structure_and_ecology")
OUT <- file.path(BASE, "figures", "main_2026")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- palettes (data/metadata/Color_Mapping_Guide.txt; labels with corrected spelling) -------------------------
HAB_COL <- c(
  "Digestive" = "#F16913", "Respiratory" = "#F768A1", "Other host-associated" = "#FDBE85", "Plants" = "#31A354", "Microbial" = "#BF812D",
  "Freshwater" = "#41B6C4", "Marine" = "#08519C", "Thermal" = "#800026", "Soil" = "#8C510A",
  "Other environmental" = "#A6BDDB", "Bioreactor" = "#525252", "Wastewater" = "#756BB1",
  "Solid waste" = "#BA68C8", "Other engineered" = "#808080", "Unclassified" = "#D9D9D9")
VIR_COL <- c("Uroviricota" = "#B3CDE3", "Nucleocytoviricota" = "#FBB4AE", "Phixviricota" = "#FFF2AE",
             "Preplasmiviricota" = "#CCEBC5", "Hofneiviricota" = "#DECBE4", "Unassigned" = "#ECEFF1")
HOST_COL <- c("Bacillota" = "#5C6BC0", "Pseudomonadota" = "#00897B", "Actinobacteriota" = "#C0CA33",
              "Bacteroidota" = "#E57373", "Cyanobacteriota" = "#00E676", "Others" = "#8E24AA", "Unknown" = "#ECEFF1")
# Main-text families only; PF00768 (homologous PBPs split from PF13354 by the clan bitscore assignment) is shown
# in the supplement
FAM_LAB <- c("PF00144" = "PF00144 (β-lactamase)", "PF00905" = "PF00905 (Transpeptidase)",
             "PF13354" = "PF13354 (β-lactamase2)")
FAM_SHAPE <- c("PF00144" = 17, "PF00905" = 16, "PF13354" = 15)
ZONES <- c("Tropical", "Temperate", "Polar")

theme_fig <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey90", colour = "grey40"),
        strip.text = element_text(size = 8, face = "bold"), legend.key.size = unit(3.5, "mm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"),
        plot.tag = element_text(size = 14, face = "bold"))

# ---- A: map ------------------------------------------------------------------------------------------------
geo <- read_tsv(file.path(ECO, "Sample_Geography_2026.tsv"), col_types = cols(.default = "c"), progress = FALSE) %>%
  mutate(across(c(latitude, longitude), as.numeric), across(starts_with("n_"), as.integer))
pos <- geo %>% filter(!is.na(latitude)) %>%
  pivot_longer(c(n_PF00144, n_PF00905, n_PF13354), names_to = "family", values_to = "n",
               names_prefix = "n_") %>% filter(n > 0) %>%
  mutate(n_hyd = case_when(family == "PF00144" ~ n_hyd_PF00144, family == "PF00905" ~ n_hyd_PF00905,
                           TRUE ~ n_hyd_PF13354))
site_pts <- pos %>% group_by(site_id, family) %>%
  summarise(latitude = mean(latitude), longitude = mean(longitude), n_samples = n(),
            hydrolytic = any(n_hyd > 0),
            habitat = names(sort(table(habitat), decreasing = TRUE))[1], .groups = "drop") %>%
  mutate(habitat = factor(habitat, levels = names(HAB_COL)), family = factor(family, levels = names(FAM_LAB)),
         outline = factor(ifelse(hydrolytic, "Verified hydrolytic SBL", "Homologs only"),
                          levels = c("Verified hydrolytic SBL", "Homologs only"))) %>%
  arrange(hydrolytic)  # hydrolytic sites drawn last, on top

robin <- "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"
world <- st_as_sf(maps::map("world", fill = TRUE, plot = FALSE, wrap = c(-180, 180))) %>% st_make_valid() %>%
  st_transform(robin)
edge <- seq(-180, 180, 1)
outline <- st_sfc(st_polygon(list(cbind(
  c(rep(-180, 181), edge, rep(180, 181), rev(edge)),
  c(seq(-90, 90), rep(90, length(edge)), seq(90, -90), rep(-90, length(edge)))))), crs = 4326) %>%
  st_transform(robin)
lat_lines <- lapply(c(23.44, -23.44, 66.56, -66.56), function(l)
  st_linestring(cbind(seq(-180, 180, 1), l))) %>% st_sfc(crs = 4326) %>% st_transform(robin)
pts <- st_as_sf(site_pts, coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(robin)
FAM_SHAPE_FILL <- c("PF00144" = 24, "PF00905" = 21, "PF13354" = 22)
bb_map <- st_bbox(outline)
KEY_XY <- c(bb_map[["xmin"]] + 4e5, bb_map[["ymax"]] - 6e5)  # top-left corner, outside the globe

# fill = habitat; a black outline marks family x site combinations with a verified hydrolytic SBL
pA <- ggplot() +
  geom_sf(data = outline, fill = "white", colour = "grey30", linewidth = 0.3) +
  geom_sf(data = world, fill = "grey60", colour = NA) +
  geom_sf(data = lat_lines, linetype = "dashed", linewidth = 0.25, colour = "grey20") +
  geom_sf(data = filter(pts, !hydrolytic), aes(shape = family, fill = habitat, colour = outline, size = n_samples),
          stroke = 0, alpha = 0.9) +
  geom_sf(data = filter(pts, hydrolytic), aes(shape = family, fill = habitat, colour = outline, size = n_samples),
          stroke = 0.6, alpha = 0.95) +
  scale_shape_manual(values = FAM_SHAPE_FILL, labels = FAM_LAB, name = "Pfam family") +
  scale_fill_manual(values = HAB_COL, drop = TRUE, name = "Habitat") +
  scale_colour_manual(values = c("Verified hydrolytic SBL" = "black", "Homologs only" = "transparent"),
                      breaks = c("Verified hydrolytic SBL", "Homologs only"), name = "Outline", guide = "none") +
  annotate("point", x = KEY_XY[1], y = KEY_XY[2], shape = 22, size = 2.8, fill = "grey80", colour = "black", stroke = 0.7) +
  annotate("text", x = KEY_XY[1] + 4.5e5, y = KEY_XY[2], hjust = 0, size = 2.1, lineheight = 0.9,
           label = "Black outline: site with a
verified hydrolytic SBL") +
  scale_size(range = c(1.3, 4.2), breaks = c(1, 5, 10, 20), name = "Samples") +
  coord_sf(crs = robin, xlim = st_bbox(outline)[c(1, 3)] * 1.005, ylim = st_bbox(outline)[c(2, 4)] * 1.01,
           expand = FALSE, datum = NA) +
  guides(shape = guide_legend(override.aes = list(size = 2.6, fill = "grey45", colour = "transparent"), ncol = 1, order = 1),
         colour = "none",
         fill = guide_legend(override.aes = list(shape = 22, size = 2.8, colour = "transparent"), nrow = 3, order = 3),
         size = guide_legend(override.aes = list(shape = 21, fill = "grey45", colour = "transparent"), ncol = 1, order = 4,
                             position = "inside")) +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.title.position = "top",
        legend.position.inside = c(0.005, 0.01), legend.justification.inside = c(0, 0),
        legend.key.size = unit(2.9, "mm"), legend.text = element_text(size = 6.2),
        legend.title = element_text(size = 7, face = "bold"), legend.box.spacing = unit(1, "mm"),
        legend.spacing.x = unit(1.5, "mm"), plot.tag = element_text(size = 14, face = "bold"),
        plot.margin = margin(2, 2, 2, 2))

# ---- Supplementary A: climate-zone prevalence (sites, IMG/VR + MOSAiC) -------------------------------------
prev <- read_tsv(file.path(ECO, "Climate_Zone_Prevalence_2026.tsv"), show_col_types = FALSE) %>%
  filter(scope == "IMGVR_plus_MOSAiC", unit == "site", family %in% names(FAM_LAB)) %>%
  mutate(family = factor(family, levels = names(FAM_LAB), labels = sub(" (", "\n(", FAM_LAB, fixed = TRUE)),
         zone_lab = sprintf("%s\n(%d/%d)", climate_zone, positive, total),
         climate_zone = factor(climate_zone, levels = ZONES))
stars <- function(q) ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "ns")))
brk <- prev %>% group_by(family) %>% summarise(
  top = max(ci95_high), .groups = "drop") %>%
  left_join(prev %>% filter(climate_zone == "Tropical") %>%
              transmute(family, q_TT = BH_q_vs_Temperate, q_TP = BH_q_vs_Polar), by = "family") %>%
  left_join(prev %>% filter(climate_zone == "Temperate") %>% transmute(family, q_MP = BH_q_vs_Polar), by = "family")
brk_long <- bind_rows(
  brk %>% transmute(family, x1 = 1, x2 = 2, y = top * 1.08, lab = stars(q_TT)),
  brk %>% transmute(family, x1 = 2, x2 = 3, y = top * 1.08, lab = stars(q_MP)),
  brk %>% transmute(family, x1 = 1, x2 = 3, y = top * 1.22, lab = stars(q_TP)))
pZone <- ggplot(prev, aes(climate_zone, prevalence)) +
  geom_errorbar(aes(ymin = ci95_low, ymax = ci95_high), width = 0.15, linewidth = 0.4) +
  geom_point(size = 2) +
  geom_segment(data = brk_long, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = brk_long, aes(x = (x1 + x2) / 2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.3, size = 2.6) +
  facet_wrap(~family, nrow = 1, scales = "free_y") +
  scale_x_discrete(labels = function(z) z) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1), expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Prevalence\n(positive sites / total sites)") +
  theme_fig + theme(axis.text.x = element_text(size = 7))
zone_labels <- prev %>% select(family, climate_zone, zone_lab)
pZone <- pZone + geom_text(data = zone_labels, aes(x = climate_zone, y = 0, label = zone_lab), inherit.aes = FALSE,
                     vjust = 1.25, size = 1.85, colour = "grey20", lineheight = 0.9) +
  coord_cartesian(clip = "off") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                        plot.margin = margin(5, 5, 20, 5), panel.spacing.x = unit(5, "mm"))

# ---- C: protein architecture - six states per family, plus representative domain architectures --------------
# The eight classes of classify_sbl_architecture_2026.py are collapsed to the six states used for the ancestral
# state reconstruction (asr_architecture_tips_2026.py), so the figure and the reconstruction share one coding.
STATE_LAB <- c("Hydrolytic" = "Verified hydrolytic SBL", "Virion_fusion" = "Virion-structure fusion",
               "Surface_lysis_fusion" = "Cell-surface / lysis fusion", "Complete_PBP" = "Complete PBP",
               "Other" = "Other fusion", "SBL_only" = "No fusion detected")
STATE_COL <- c("Hydrolytic" = "#B2182B", "Virion_fusion" = "#2166AC", "Surface_lysis_fusion" = "#4DAF4A",
               "Complete_PBP" = "#762A83", "Other" = "#969696", "SBL_only" = "#DDDDDD")
COLLAPSE <- c("Verified hydrolytic SBL" = "Hydrolytic", "Virion structure" = "Virion_fusion",
              "Lysis / peptidoglycan hydrolysis" = "Surface_lysis_fusion", "Cell-surface binding" = "Surface_lysis_fusion",
              "Cell-wall synthesis (PBP)" = "Complete_PBP", "SBL domain only" = "SBL_only",
              "Signalling / sensor" = "Other", "Other fusion" = "Other")
arch_sum <- read_tsv(file.path(ECO, "Viral_SBL_Architecture_Summary_2026.tsv"), show_col_types = FALSE) %>%
  filter(family %in% names(FAM_LAB)) %>%
  mutate(state = COLLAPSE[architecture_class]) %>%
  group_by(family, state) %>% summarise(n = sum(n), total = first(total), .groups = "drop") %>%
  mutate(share = n / total, state = factor(state, levels = names(STATE_LAB)),
         family = factor(family, levels = rev(names(FAM_LAB))))
stopifnot(!anyNA(arch_sum$state))
seg_lab <- arch_sum %>% arrange(family, state) %>% group_by(family) %>%
  mutate(xmid = cumsum(share) - share / 2) %>% ungroup() %>%
  filter(share >= 0.07) %>%
  mutate(lab = comma(n), txt = ifelse(state %in% c("SBL_only", "Other"), "grey15", "white"))
fam_tot <- arch_sum %>% distinct(family, total)
# representative architectures (drawn beside the matching heat-map row of panel B): median protein length and
# median domain envelopes over all sequences that share the exact domain order (domain_architecture coordinates are
# envelope positions in aa); n = curated viral sequences with that architecture
arch_seq <- read_tsv(file.path(ECO, "Viral_SBL_Architecture_Classes_2026.tsv"), show_col_types = FALSE) %>%
  filter(main_text) %>%
  left_join(read_tsv(file.path(BASE, "results", "reanalysis_2026", "01_curation_and_context",
                               "Final_Viral_SBL_Catalog.tsv"), show_col_types = FALSE,
                     col_select = c(uid, seq_len)) %>% distinct(uid, .keep_all = TRUE), by = "uid") %>%
  mutate(sig = vapply(strsplit(domain_architecture, "+", fixed = TRUE), function(p)
    paste(sub("\\(.*", "", sub("^PF\\d+:", "", p)), collapse = "+"), ""))
REPS <- tribble(
  ~family,   ~class,                      ~sig,                                        ~title,
  "PF13354", "Verified hydrolytic SBL",   "Beta-lactamase2",                           "Free-standing SBL",
  "PF00144", "Virion structure",          "Phage_spike_CBM+Gly_rich_2+Beta-lactamase", "Tail spike",
  "PF00144", "Virion structure",          "Beta-lactamase+Phage_tail_CBM",             "Tail fibre",
  "PF13354", "Virion structure",          "TMP_3+Endolysin+Beta-lactamase2",           "Tape measure",
  "PF00144", "Cell-surface binding",      "SLH+Beta-lactamase",                        "S-layer anchored",
  "PF13354", "Lysis / peptidoglycan hydrolysis", "VanY+Beta-lactamase2",               "VanY fusion",
  "PF00905", "Cell-wall synthesis (PBP)", "Transgly+Transpeptidase",                   "Class A PBP",
  "PF00905", "Cell-wall synthesis (PBP)", "Transgly+Transpeptidase+PASTA",             "Class A PBP + PASTA",
  "PF00144", "SBL domain only",           "Beta-lactamase",                            "SBL domain alone")
DOM_LAB <- c("Beta-lactamase" = "SBL", "Beta-lactamase2" = "SBL", "Transpeptidase" = "SBL",
             "Phage_spike_CBM" = "CBM", "Gly_rich_2" = "G", "Phage_tail_CBM" = "CBM", "TMP_3" = "TMP",
             "Endolysin" = "Lys", "SLH" = "SLH", "VanY" = "VanY", "Transgly" = "TG", "PASTA" = "P")
rep_dom <- bind_rows(lapply(seq_len(nrow(REPS)), function(i) {
  r <- REPS[i, ]
  s <- filter(arch_seq, family == r$family, architecture_class == r$class, sig == r$sig)
  stopifnot(nrow(s) > 0)
  parts <- strsplit(s$domain_architecture, "+", fixed = TRUE)
  co <- lapply(seq_along(parts[[1]]), function(j) {
    pj <- vapply(parts, `[`, "", j)
    m <- regmatches(pj, regexec("^PF\\d+:([^()]+)\\((\\d+)-(\\d+)\\)", pj))
    tibble(dom = m[[1]][2], start = median(as.numeric(vapply(m, `[`, "", 3))),
           end = median(as.numeric(vapply(m, `[`, "", 4))))
  })
  bind_rows(co) %>% mutate(row = i, len = median(s$seq_len), n = nrow(s), family = r$family, class = r$class,
                           title = r$title)
})) %>%
  mutate(state = COLLAPSE[class],
         role = case_when(dom %in% c("Beta-lactamase", "Beta-lactamase2", "Transpeptidase") ~ "sbl",
                          dom == "Gly_rich_2" ~ "linker", TRUE ~ "partner"),
         fill = case_when(role == "sbl" & state == "Hydrolytic" ~ STATE_COL["Hydrolytic"],
                          role == "sbl" ~ "#525252", role == "linker" ~ "#E6E6E6", TRUE ~ STATE_COL[state]),
         txt = ifelse(role == "linker", "grey20", "white"), lab = DOM_LAB[dom])
stopifnot(!anyNA(rep_dom$lab))
# envelopes of neighbouring domains can overlap by a few residues; trim so boxes abut instead of stacking
rep_dom <- rep_dom %>% group_by(row) %>% arrange(start, .by_group = TRUE) %>%
  mutate(end = pmin(end, lead(start, default = Inf) - 1)) %>% ungroup()
# domains keep their true length; an undomained stretch longer than GAP_MAX aa is drawn GAP_MAX long with a break
# mark, so the 3,000-aa tape measure protein does not compress the other tracks
GAP_MAX <- 250
squeeze <- function(d) {
  len <- d$len[1]
  knots <- sort(unique(c(1, d$start, d$end, len)))
  in_dom <- vapply(seq_len(length(knots) - 1), function(i)
    any(d$start <= knots[i] & d$end >= knots[i + 1]), logical(1))
  w <- diff(knots); w_disp <- ifelse(!in_dom & w > GAP_MAX, GAP_MAX, w)
  disp <- c(0, cumsum(w_disp))
  brk <- (knots[-1] + knots[-length(knots)]) / 2
  list(map = function(x) approx(knots, disp, xout = x)$y, len_disp = max(disp),
       breaks = approx(knots, disp, xout = brk[!in_dom & w > GAP_MAX])$y)
}
sq <- lapply(split(rep_dom, rep_dom$row), squeeze)
rep_dom <- rep_dom %>% group_by(row) %>%
  mutate(xs = sq[[as.character(row[1])]]$map(start), xe = sq[[as.character(row[1])]]$map(end)) %>% ungroup()

# ---- B: one row per protein-architecture state ------------------------------------------------------------------
#   left    curated viral sequences per family (circle area; the top row gives the family totals)
#   middle  habitat enrichment of Uroviricota UViGs against the IMG/VR high-confidence background. Colour: log2
#           enrichment at UViG level. Stars: BH q of the UViG-level test, shown only when the test with at most one
#           UViG per metagenome is also significant (the larger q sets the stars). x: no UViG of the state observed in
#           the habitat although the background predicts some (significant depletion), drawn at the scale minimum.
#   right   representative domain architectures of the state, to scale (one shared aa scale)
# 'Other fusion' (40 sequences, 20 high-confidence UViGs) is too small to test and is left out.
SHOWN <- c("Digestive", "Respiratory", "Other host-associated", "Plants", "Microbial", "Freshwater", "Marine",
           "Thermal", "Soil", "Bioreactor", "Wastewater", "Solid waste")
B_ROWS <- c("All_SBL", "Hydrolytic", "Virion_fusion", "Surface_lysis_fusion", "Complete_PBP", "SBL_only")
B_LAB <- c("All_SBL" = "All SBL carriers", STATE_LAB)
NR <- length(B_ROWS)
row_y <- function(g) NR + 1 - match(g, B_ROWS)
FAMS <- names(FAM_LAB)
FX <- c(1.05, 2.15, 3.2)              # family columns of the circle matrix
HX0 <- 4.6                            # first heat-map column
X0 <- HX0 + length(SHOWN) + 0.15      # schematic origin
X_W <- 6.2                            # drawn width of the longest protein (column units)

senr <- read_tsv(file.path(ECO, "Habitat_Enrichment_By_Architecture_State_2026.tsv"), show_col_types = FALSE) %>%
  filter(comparison == "vs_background", viral_phylum == "Uroviricota", habitat %in% SHOWN, group %in% B_ROWS)
enr <- senr %>% filter(unit == "UViG") %>%
  left_join(senr %>% filter(unit == "metagenome") %>% select(group, habitat, q_mg = BH_q), by = c("group", "habitat")) %>%
  mutate(q = pmax(BH_q, q_mg, na.rm = FALSE),
         stars = case_when(is.na(q) ~ "", q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ ""),
         zero = sbl_uvig_in_habitat == 0 & !is.na(BH_q) & BH_q < 0.05,
         val = ifelse(zero, -6, pmax(pmin(log2_enrichment, 6), -6)),
         x = HX0 + match(habitat, SHOWN) - 1, y = row_y(group))
row_key <- enr %>% distinct(group, y, n = sbl_uvig_phylum_total) %>%
  mutate(lab = sprintf("%s\n%s UViGs", B_LAB[group], comma(n)), col = c(All_SBL = "grey35", STATE_COL)[group])

# sequences per family x state (circles); top row = family totals
bub <- bind_rows(
  arch_sum %>% mutate(state = as.character(state), family = as.character(family)) %>%
    filter(state %in% B_ROWS) %>% select(family, state, n),
  arch_sum %>% distinct(family = as.character(family), n = total) %>% mutate(state = "All_SBL")) %>%
  mutate(x = FX[match(family, FAMS)], y = row_y(state), col = c(All_SBL = "grey35", STATE_COL)[state],
         txt = ifelse(state %in% c("SBL_only"), "grey10", "white"), big = n >= 250)

# schematic tracks, stacked within the row of their state
f_aa <- X_W / max(vapply(sq, `[[`, numeric(1), "len_disp"))
TRACK_SP <- 0.27; BOX_H <- 0.095
track <- rep_dom %>% distinct(row, state, title, family, n, len) %>%
  mutate(len_disp = vapply(as.character(row), function(r) sq[[r]]$len_disp, numeric(1))) %>%
  group_by(state) %>% mutate(k = row_number(), nk = n()) %>% ungroup() %>%
  mutate(y = row_y(state) + ((nk - 1) / 2 - (k - 1)) * TRACK_SP, x_end = X0 + len_disp * f_aa,
         lab = sprintf("%s (%s, %s)", title, family, comma(n)))
stopifnot(!anyNA(track$y))
seg_dom <- rep_dom %>% left_join(track %>% select(row, y), by = "row") %>%
  mutate(x1 = X0 + xs * f_aa, x2 = X0 + xe * f_aa, inside = x2 - x1 >= 0.11 * nchar(lab) + 0.05 | role == "sbl")
seg_brk <- bind_rows(lapply(names(sq), function(r) if (length(sq[[r]]$breaks))
  tibble(row = as.integer(r), xb = sq[[r]]$breaks))) %>%
  left_join(track %>% select(row, y), by = "row") %>% mutate(x = X0 + xb * f_aa)

# hand-drawn colour bar and keys in the strip under the rows
CB_X <- HX0 - 0.5 + seq(0, 4.8, length.out = 49); CB_V <- seq(-6, 6, length.out = 49)
cbar <- tibble(x = CB_X, val = CB_V)
KEY_Y <- 0.18
pB <- ggplot() +
  # row guides
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = NR - 0.5, ymax = NR + 0.5, fill = "grey96") +
  geom_point(data = row_key, aes(x = 0.2, y = y), shape = 22, size = 2.7, stroke = 0.3, fill = row_key$col,
             colour = "grey30") +
  # family composition
  geom_point(data = bub, aes(x, y, size = n), shape = 21, fill = bub$col, colour = "grey25", stroke = 0.2) +
  geom_text(data = filter(bub, big), aes(x, y, label = comma(n)), colour = filter(bub, big)$txt, size = 1.75,
            fontface = "bold") +
  geom_text(data = filter(bub, !big), aes(x, y - 0.3, label = comma(n)), size = 1.7, colour = "grey15") +
  scale_size_area(max_size = 8.2, guide = "none") +
  # enrichment heat map
  geom_tile(data = enr, aes(x, y, fill = val), colour = "white", linewidth = 0.55, height = 0.94, width = 0.97) +
  geom_text(data = filter(enr, !zero), aes(x, y, label = stars), size = 2.2, vjust = 0.75) +
  geom_text(data = filter(enr, zero), aes(x, y, label = paste0("×", stars)), size = 2.2, colour = "white",
            vjust = 0.6) +
  # domain architectures
  geom_segment(data = track, aes(x = X0, xend = x_end, y = y, yend = y), linewidth = 0.4, colour = "grey50") +
  geom_text(data = seg_brk, aes(x = x, y = y), label = "//", size = 1.6, colour = "grey30", fontface = "bold") +
  geom_rect(data = seg_dom, aes(xmin = x1, xmax = x2, ymin = y - BOX_H, ymax = y + BOX_H), fill = seg_dom$fill,
            colour = "grey15", linewidth = 0.15) +
  geom_text(data = filter(seg_dom, inside), aes(x = (x1 + x2) / 2, y = y, label = lab),
            colour = filter(seg_dom, inside)$txt, size = 1.5, fontface = "bold") +
  geom_text(data = filter(seg_dom, !inside), aes(x = (x1 + x2) / 2, y = y + BOX_H + 0.01, label = lab),
            vjust = 0, size = 1.4, colour = "grey20") +
  geom_text(data = track, aes(x = x_end + 0.1, y = y, label = lab), hjust = 0, size = 1.7, colour = "grey15") +
  annotate("text", x = X0, y = NR + 0.12, hjust = 0, size = 2.05, fontface = "bold", colour = "grey15",
           label = "Representative domain architectures") +
  annotate("text", x = X0, y = NR - 0.2, hjust = 0, size = 1.75, colour = "grey30",
           label = "to scale; // shortened linker; (family, sequences)") +
  # keys under the rows
  geom_tile(data = cbar, aes(x = x, y = KEY_Y, fill = val), height = 0.16, width = 0.11) +
  annotate("text", x = HX0 - 0.5 + c(0, 2.4, 4.8), y = KEY_Y - 0.17, label = c("≤ -6", "0", "≥ 6"),
           size = 1.75, vjust = 1) +
  annotate("text", x = HX0 - 0.65, y = KEY_Y, label = "log2 enrichment", hjust = 1, size = 1.85, fontface = "bold") +
  annotate("tile", x = HX0 + 5.1, y = KEY_Y, width = 0.32, height = 0.2, fill = "#2166AC") +
  annotate("text", x = HX0 + 5.1, y = KEY_Y, label = "×", colour = "white", size = 2) +
  annotate("text", x = HX0 + 5.35, y = KEY_Y, hjust = 0, size = 1.75,
           label = "none observed (significant depletion);  * q < 0.05, ** < 0.01, *** < 0.001") +
  annotate("segment", x = X0 + 4.2, xend = X0 + 4.2 + 200 * f_aa, y = KEY_Y, yend = KEY_Y, linewidth = 0.4) +
  annotate("text", x = X0 + 4.2 + 200 * f_aa + 0.1, y = KEY_Y, label = "200 aa", hjust = 0, size = 1.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-6, 6),
                       na.value = "grey91", guide = "none") +
  scale_x_continuous(breaks = c(FX, HX0 + seq_along(SHOWN) - 1),
                     labels = c(FAMS, sub("Other host-associated", "Other host-assoc.", SHOWN)),
                     position = "top", expand = expansion(add = 0)) +
  scale_y_continuous(breaks = row_key$y, labels = row_key$lab, expand = expansion(add = 0)) +
  coord_cartesian(xlim = c(0, X0 + X_W + 4.15), ylim = c(-0.08, NR + 0.5), clip = "off") +
  labs(x = NULL, y = NULL) + theme_void(base_size = 9) +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 6.6, margin = margin(b = 1)),
        axis.text.y = element_text(size = 6.6, hjust = 1, lineheight = 0.9, margin = margin(r = 1)),
        plot.tag = element_text(size = 14, face = "bold"), plot.margin = margin(2, 4, 2, 2))

# ---- Supplementary B: habitat -> viral phylum -> host phylum ------------------------------------------------
# Parallel sets (ggforce): gaps between nodes, flows coloured by habitat, labels outside the axes.
# Small categories are pooled so every node can carry a readable label.
D_HAB <- c(HAB_COL[c("Digestive", "Respiratory", "Other host-associated", "Plants", "Microbial", "Freshwater", "Marine", "Soil",
                     "Other environmental")], "Engineered" = "#525252", "Unclassified" = "#D9D9D9")
D_VIR <- c(VIR_COL[c("Uroviricota", "Nucleocytoviricota")], "Other / unassigned" = "#DECBE4")
D_HOST <- c(HOST_COL[c("Bacillota", "Actinobacteriota", "Bacteroidota", "Pseudomonadota", "Cyanobacteriota", "Others")],
            "Unknown" = "#CFD8DC")
flows <- read_tsv(file.path(ECO, "Habitat_ViralPhylum_Host_Flows_2026.tsv"), show_col_types = FALSE) %>%
  filter(carries_main_family) %>%
  mutate(habitat = case_when(habitat %in% c("Bioreactor", "Wastewater", "Solid waste", "Other engineered") ~ "Engineered",
                             habitat == "Thermal" ~ "Other environmental", TRUE ~ habitat),
         viral_phylum = ifelse(viral_phylum %in% names(D_VIR), viral_phylum, "Other / unassigned"),
         host_phylum = ifelse(host_phylum %in% names(D_HOST), host_phylum, "Others")) %>%
  group_by(habitat, viral_phylum, host_phylum) %>% summarise(n = sum(n_uvig), .groups = "drop") %>%
  rename(Habitat = habitat, `Viral phylum` = viral_phylum, `Viral host` = host_phylum)
SEP <- 0.012
ps <- gather_set_data(flows, 1:3) %>%
  mutate(x = factor(x, levels = c("Habitat", "Viral phylum", "Viral host")),
         y = factor(y, levels = c(names(D_HAB), names(D_VIR), names(D_HOST))),
         Habitat = factor(Habitat, levels = names(D_HAB)))
tot <- sum(flows$n)
base_D <- ggplot(ps, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = Habitat), alpha = 0.55, axis.width = 0.1, sep = SEP, colour = NA) +
  geom_parallel_sets_axes(aes(fill = after_stat(label)), axis.width = 0.1, sep = SEP, colour = "white", linewidth = 0.2)
# label positions taken from the rendered axes so they always sit at the centre of their node
nodes <- layer_data(base_D, 2) %>%
  mutate(xn = round((xmin + xmax) / 2), y = as.character(label)) %>%
  group_by(xn, y) %>% summarise(ymid = (min(ymin) + max(ymax)) / 2, .groups = "drop") %>%
  left_join(ps %>% group_by(x, y) %>% summarise(n = sum(n), .groups = "drop") %>%
              transmute(xn = as.numeric(x), y = as.character(y), n), by = c("xn", "y")) %>%
  mutate(lab = sprintf("%s (%s)", y, comma(n)))
# host labels: keep a minimum vertical spacing (small adjacent nodes), with a short leader to each node
spread <- function(y, gap) {
  o <- order(y); v <- y[o]
  for (k in seq_along(v)[-1]) v[k] <- max(v[k], v[k - 1] + gap)
  shift <- mean(v - y[o]); v <- v - pmax(0, shift)
  for (k in rev(seq_along(v))[-1]) v[k] <- min(v[k], v[k + 1] - gap)
  out <- numeric(length(y)); out[o] <- v; out
}
host_lab <- filter(nodes, xn == 3) %>% mutate(ylab = spread(ymid, 0.035 * tot))
pFlow <- base_D +
  geom_text(data = filter(nodes, xn == 1), aes(x = xn - 0.07, y = ymid, label = lab), inherit.aes = FALSE,
            hjust = 1, size = 2.05) +
  geom_label(data = filter(nodes, xn == 2), aes(x = xn + 0.07, y = ymid, label = lab), inherit.aes = FALSE,
             hjust = 0, size = 2.05, linewidth = 0, label.padding = unit(0.6, "mm"), fill = alpha("white", 0.75)) +
  geom_segment(data = host_lab, aes(x = 3.055, xend = 3.1, y = ymid, yend = ylab), inherit.aes = FALSE,
               linewidth = 0.2, colour = "grey40") +
  geom_text(data = host_lab, aes(x = 3.11, y = ylab, label = lab), inherit.aes = FALSE, hjust = 0, size = 2.05) +
  scale_fill_manual(values = c(D_HAB, D_VIR, D_HOST), guide = "none") +
  scale_x_discrete(position = "top", expand = expansion(add = c(0.95, 0.95))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, caption = sprintf("%s SBL-carrying UViGs", comma(tot))) + theme_void(base_size = 9) +
  theme(axis.text.x = element_text(face = "bold", size = 8, margin = margin(b = 2)),
        plot.caption = element_text(size = 6.5, colour = "grey30"), plot.tag = element_text(size = 14, face = "bold"),
        plot.margin = margin(5, 5, 5, 5))

# ---- assemble ---------------------------------------------------------------------------------------------
save_fig <- function(fig, stem, h) {
  ggsave(paste0(stem, ".pdf"), fig, width = 183, height = h, units = "mm", device = cairo_pdf)
  ggsave(paste0(stem, ".png"), fig, width = 183, height = h, units = "mm", dpi = 300)
  ggsave(paste0(stem, ".tiff"), fig, width = 183, height = h, units = "mm", dpi = 300, compression = "lzw")
  cat("Written:", stem, "{.pdf,.png,.tiff}\n")
}
fig1 <- wrap_plots(wrap_elements(full = pA + labs(tag = "A")), wrap_elements(full = pB + labs(tag = "B")),
                   ncol = 1, heights = c(1.22, 1))
save_fig(fig1, file.path(OUT, "Fig1_Global_Distribution_and_Habitat"), 200)

SUPP <- file.path(BASE, "figures", "supplementary_2026")
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)
figS <- wrap_plots(wrap_elements(full = pZone + labs(tag = "A")), wrap_elements(full = pFlow + labs(tag = "B")),
                   ncol = 1, heights = c(0.8, 1.3))
save_fig(figS, file.path(SUPP, "FigS_Climate_Zone_Prevalence_and_Habitat_Flows"), 200)
