library(jsonlite)
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)

# =========================================================
# Helper functions
# =========================================================

.read_fit <- function(path) {
  j <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  j$fits$`Standard MG94`
}

.extract_global_omega <- function(fit) {
  rd <- fit$`Rate Distributions`
  ci <- fit$`Confidence Intervals`

  omega_key <- grep("non-synonymous/synonymous rate ratio", names(rd), value = TRUE)[1]
  ci_key    <- grep("non-synonymous/synonymous rate ratio", names(ci), value = TRUE)[1]

  tibble(
    group = "Global",
    omega = unname(rd[[omega_key]]),
    LB    = ci[[ci_key]]$LB,
    UB    = ci[[ci_key]]$UB
  )
}

.extract_partition_omegas <- function(fit) {
  rd <- fit$`Rate Distributions`
  ci <- fit$`Confidence Intervals`

  keys <- grep("non-synonymous/synonymous rate ratio for \\*", names(rd), value = TRUE)

  tibble(raw_key = keys) %>%
    mutate(
      group = stringr::str_match(raw_key, "for \\*(.+)\\*")[, 2],
      omega = sapply(raw_key, function(k) rd[[k]]),
      LB    = sapply(raw_key, function(k) ci[[k]]$LB),
      UB    = sapply(raw_key, function(k) ci[[k]]$UB)
    ) %>%
    select(group, omega, LB, UB)
}

.calc_lrt <- function(full_fit, reduced_fit, comparison) {
  logL_full <- full_fit$`Log Likelihood`
  k_full    <- full_fit$`estimated parameters`

  logL_red <- reduced_fit$`Log Likelihood`
  k_red    <- reduced_fit$`estimated parameters`

  lrt <- 2 * (logL_full - logL_red)
  df  <- k_full - k_red
  p   <- pchisq(lrt, df = df, lower.tail = FALSE)

  tibble(
    comparison    = comparison,
    logL_full     = logL_full,
    k_full        = k_full,
    logL_reduced  = logL_red,
    k_reduced     = k_red,
    LRT           = lrt,
    df            = df,
    p_value       = p
  )
}

.p_to_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

.format_p_label <- function(p) {
  paste0(.p_to_stars(p), "\n", "p=", format(p, scientific = TRUE, digits = 2))
}

# =========================================================
# Main function
# =========================================================
make_mg94_family_plot <- function(
  family_name,
  global_json,
  full_json,
  bg_fg_json,
  bg_nt_json,
  fg_nt_json,
  x_labels = c(
    "Background" = "Background",
    "Foreground" = "Foreground",
    "Nested"     = "Nested"
  ),
  y_limits = NULL,
  show_overall_box = TRUE
) {
  # -----------------------------
  # 1. Read fits
  # -----------------------------
  fit_global <- .read_fit(global_json)
  fit_full   <- .read_fit(full_json)
  fit_bg_fg  <- .read_fit(bg_fg_json)
  fit_bg_nt  <- .read_fit(bg_nt_json)
  fit_fg_nt  <- .read_fit(fg_nt_json)

  # -----------------------------
  # 2. Extract omega summaries
  # -----------------------------
  global_df <- .extract_global_omega(fit_global)
  full_raw  <- .extract_partition_omegas(fit_full)

  plot_df <- tibble(
    group = c("Background", "Foreground", "Nested"),
    x = c(1, 2, 3),
    label = unname(x_labels[c("Background", "Foreground", "Nested")])
  ) %>%
    left_join(full_raw, by = "group")

  if (any(is.na(plot_df$omega))) {
    stop("Missing omega estimates in full partitioned model for: ", family_name)
  }

  # -----------------------------
  # 3. Pairwise reduced-model LRT
  # -----------------------------
  pair_df <- bind_rows(
    .calc_lrt(fit_full, fit_bg_fg, "BG vs FG"),
    .calc_lrt(fit_full, fit_fg_nt, "FG vs NT"),
    .calc_lrt(fit_full, fit_bg_nt, "BG vs NT")
  ) %>%
    mutate(
      stars = .p_to_stars(p_value),
      p_label = vapply(p_value, .format_p_label, character(1))
    )

  # -----------------------------
  # 4. Overall LRT: full vs global
  # -----------------------------
  overall_lrt <- .calc_lrt(fit_full, fit_global, "Global vs Partitioned")

  overall_label <- sprintf(
    "Global vs partitioned\nLRT=%.2f, df=%d, p=%.2e",
    overall_lrt$LRT,
    overall_lrt$df,
    overall_lrt$p_value
  )

  # -----------------------------
  # 5. Significance brackets
  # -----------------------------
  ann_df <- tibble(
    comparison = c("BG vs FG", "FG vs NT", "BG vs NT"),
    x1 = c(1, 2, 1),
    x2 = c(2, 3, 3)
  ) %>%
    left_join(pair_df, by = "comparison")

  y_base <- max(plot_df$UB, na.rm = TRUE) + 0.010
  ann_df$y <- c(y_base, y_base + 0.018, y_base + 0.036)

  ymax_auto <- max(ann_df$y) + 0.045
  ymin_auto <- min(plot_df$LB, na.rm = TRUE) - 0.010

  if (is.null(y_limits)) {
    y_limits <- c(ymin_auto, ymax_auto)
  }

  # -----------------------------
  # 6. Plot
  # -----------------------------
  p <- ggplot() +
    # Global omega band
    annotate(
      "rect",
      xmin = 0.5, xmax = 3.5,
      ymin = global_df$LB, ymax = global_df$UB,
      fill = "grey70", alpha = 0.18
    ) +
    geom_hline(
      yintercept = global_df$omega,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    annotate(
      "text",
      x = 3.35, y = global_df$omega,
      label = sprintf("Global omega = %.3f", global_df$omega),
      hjust = 1, vjust = -0.35,
      size = 3.2, color = "grey30"
    ) +

    # Overall significance box
    {if (show_overall_box)
      annotate(
        "label",
        x = 0.72,
        y = y_limits[2] - 0.004,
        label = overall_label,
        hjust = 0, vjust = 1,
        size = 3.0,
        linewidth = 0.2,
        fill = "white"
      )
    } +

    # Full-model CI + points
    geom_errorbar(
      data = plot_df,
      aes(x = x, ymin = LB, ymax = UB),
      width = 0.07,
      linewidth = 0.5
    ) +
    geom_point(
      data = plot_df,
      aes(x = x, y = omega),
      size = 2.4
    ) +
    geom_text(
      data = plot_df,
      aes(x = x, y = UB + 0.0035, label = sprintf("%.3f", omega)),
      size = 3.0
    ) +

    # Pairwise significance brackets
    geom_segment(
      data = ann_df,
      aes(x = x1, xend = x1, y = y, yend = y + 0.0035),
      linewidth = 0.45
    ) +
    geom_segment(
      data = ann_df,
      aes(x = x1, xend = x2, y = y + 0.0035, yend = y + 0.0035),
      linewidth = 0.45
    ) +
    geom_segment(
      data = ann_df,
      aes(x = x2, xend = x2, y = y + 0.0035, yend = y),
      linewidth = 0.45
    ) +
    geom_text(
      data = ann_df,
      aes(x = (x1 + x2)/2, y = y + 0.0065, label = p_label),
      size = 2.8,
      lineheight = 0.9
    ) +

    scale_x_continuous(
      breaks = plot_df$x,
      labels = plot_df$label,
      limits = c(0.5, 3.5)
    ) +
    scale_y_continuous(
      name = expression(omega~"(dN/dS)"),
      limits = y_limits,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = family_name
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 18, r = 14, b = 8, l = 8)
    )

  # -----------------------------
  # 7. Return plot + metadata
  # -----------------------------
  list(
    plot = p,
    global_df = global_df,
    plot_df = plot_df,
    pair_df = pair_df,
    overall_lrt = overall_lrt
  )
}