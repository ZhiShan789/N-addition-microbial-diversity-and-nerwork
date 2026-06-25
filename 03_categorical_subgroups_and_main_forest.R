# Categorical subgroup analysis and main forest plots
# Clean GitHub code package for the N deposition soil microbiome synthesis.
# This file combines the R scripts actually used for the manuscript-level analyses.

# Sources combined: 2
# Generated on: 2026-06-25


# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/Meta????????.R
# ================================================================

# MISSING SOURCE FILE: Data_analysis/Remote_working_package/06_R_scripts/Meta????????.R

# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/plot_categorical_NC_like_forest.R
# ================================================================

rm(list = ls())
setwd("E:/BaiduSyncdisk/N_deposition1/Data_analysis")

# NC-like categorical forest plots.
# One figure per Domain and response set:
#   Bacterial/Fungi diversity: Shannon, richness, beta diversity, community structure
#   Bacterial/Fungi phylum composition: top phyla

required_packages <- c("dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

library(dplyr)
library(ggplot2)

out_dir <- "outputs"
plot_dir <- file.path(out_dir, "plots", "categorical_NC_like")
single_plot_dir <- file.path(plot_dir, "single_panels")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(single_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Knobs for single-panel outputs. Edit these if you want the exported panels
# narrower/wider or want the text columns moved farther from the CI bars.
single_panel_width <- 7.2
single_panel_right_margin <- 8
single_panel_left_label_pad <- 0.22
single_panel_right_label_pad <- 0.34
single_panel_qm_inset <- 0.02
single_panel_n_inset <- 0.02
single_panel_letter_offset <- 0.08

overall <- read.csv(file.path(out_dir, "meta_overall_summary.csv"), stringsAsFactors = FALSE)
subgroup <- read.csv(file.path(out_dir, "meta_categorical_subgroup_estimates_NC_style.csv"),
                     stringsAsFactors = FALSE)
tests <- read.csv(file.path(out_dir, "meta_categorical_moderator_tests_NC_style.csv"),
                  stringsAsFactors = FALSE)

moderator_order <- c("Ecosystem", "N_form_reclassified", "N_addition_position")
moderator_labels <- c(
  Ecosystem = "Ecosystem",
  N_form_reclassified = "N form",
  N_addition_position = "N addition position"
)

level_order <- list(
  Ecosystem = c("Farmland", "Forest", "Grassland", "Wetland", "Other"),
  N_form_reclassified = c("Ammonium_N", "Nitrate_N", "Ammonium_Nitrate", "Organic_N", "Mixed_N"),
  N_addition_position = c("Canopy", "Understory")
)

level_display_labels <- list(
  N_form_reclassified = c(
    Ammonium_N = "Ammonium",
    Nitrate_N = "Nitrate",
    Ammonium_Nitrate = "NH4NO3",
    Organic_N = "Urea",
    Mixed_N = "Mixed"
  )
)

display_level <- function(moderator, level) {
  labels <- level_display_labels[[moderator]]
  if (is.null(labels)) return(level)
  out <- labels[level]
  out[is.na(out)] <- level[is.na(out)]
  unname(out)
}

diversity_metrics <- c(
  "shannon",
  "richness",
  "RRBeta_ln_Dt_over_Dc",
  "RRStructure_ln_Db_over_meanDcDt"
)

diversity_metric_labels <- c(
  shannon = "Shannon diversity",
  richness = "species richness",
  simpson = "Simpson index",
  RRBeta_ln_Dt_over_Dc = "beta diversity",
  RRStructure_ln_Db_over_meanDcDt = "community structure"
)

alpha_metrics <- c("richness", "shannon", "simpson")
alpha_metric_labels <- c(
  richness = "species richness",
  shannon = "Shannon-Wiener",
  simpson = "Simpson index"
)

beta_metrics <- c(
  "RRBeta_ln_Dt_over_Dc",
  "RRStructure_ln_Db_over_meanDcDt"
)
beta_metric_labels <- c(
  RRBeta_ln_Dt_over_Dc = "beta diversity",
  RRStructure_ln_Db_over_meanDcDt = "community structure"
)

phylum_metrics <- list(
  Bacterial = c("Proteobacteria", "Actinobacteria", "Acidobacteria", "Firmicutes", "Chloroflexi"),
  Fungi = c("Ascomycota", "Basidiomycota", "Mortierellomycota", "Chytridiomycota", "Rozellomycota")
)

fmt_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.0001) return("P < 0.0001")
  paste0("P = ", formatC(p, format = "f", digits = ifelse(p < 0.01, 4, 3)))
}

fmt_qm <- function(qm) {
  if (is.na(qm)) return("QM = NA")
  paste0("QM = ", formatC(qm, format = "f", digits = ifelse(qm >= 100, 0, 2)))
}

make_row_layout <- function(domain, metrics, moderators) {
  rows <- data.frame(
    RowType = "Overall",
    Moderator = "Overall",
    Level = "Overall",
    RowLabel = "Overall",
    stringsAsFactors = FALSE
  )
  separators_after <- 1

  for (mod in moderators) {
    mod_levels <- unique(subgroup$Level[
      subgroup$Domain == domain &
        subgroup$Metric %in% metrics &
        subgroup$Moderator == mod
    ])
    preferred <- level_order[[mod]]
    if (!is.null(preferred)) {
      mod_levels <- c(preferred[preferred %in% mod_levels], sort(setdiff(mod_levels, preferred)))
    } else {
      mod_levels <- sort(mod_levels)
    }
    if (length(mod_levels) == 0) next

    start_i <- nrow(rows) + 1
    rows <- bind_rows(
      rows,
      data.frame(
        RowType = "Subgroup",
        Moderator = mod,
        Level = mod_levels,
        RowLabel = display_level(mod, mod_levels),
        stringsAsFactors = FALSE
      )
    )
    separators_after <- c(separators_after, nrow(rows))
  }

  rows$RowIndex <- seq_len(nrow(rows))
  rows$y <- nrow(rows) - rows$RowIndex + 1
  rows$RowKey <- paste(rows$Moderator, rows$Level, sep = "||")
  sep_df <- data.frame(y = nrow(rows) - separators_after + 0.5)
  list(rows = rows, separators = sep_df)
}

build_plot_data <- function(domain, metrics, metric_labels, letters_prefix = letters) {
  layout <- make_row_layout(domain, metrics, moderator_order)
  rows <- layout$rows

  panel_labels <- setNames(
    paste0(letters_prefix[seq_along(metrics)], "\n", domain, "\n", metric_labels[metrics]),
    metrics
  )

  plot_rows <- list()
  for (metric in metrics) {
    over <- overall[overall$Domain == domain & overall$Metric == metric, , drop = FALSE]
    if (nrow(over) > 0) {
      plot_rows[[length(plot_rows) + 1]] <- data.frame(
        Metric = metric,
        Panel = panel_labels[metric],
        Moderator = "Overall",
        Level = "Overall",
        estimate = over$estimate[1],
        ci_low = over$ci_low[1],
        ci_high = over$ci_high[1],
        p_value = over$p_value[1],
        k = over$k[1],
        studies = over$studies[1],
        letters = "",
        stringsAsFactors = FALSE
      )
    }

    sub <- subgroup[subgroup$Domain == domain & subgroup$Metric == metric &
                      subgroup$Moderator %in% moderator_order, , drop = FALSE]
    if (nrow(sub) > 0) {
      sub$Panel <- panel_labels[metric]
      plot_rows[[length(plot_rows) + 1]] <- sub[, c(
        "Metric", "Panel", "Moderator", "Level", "estimate", "ci_low", "ci_high",
        "p_value", "k", "studies", "letters"
      )]
    }
  }

  pdat <- bind_rows(plot_rows)
  pdat$RowKey <- paste(pdat$Moderator, pdat$Level, sep = "||")
  pdat <- left_join(pdat, rows[, c("RowKey", "y", "RowLabel", "RowType")], by = "RowKey")
  pdat <- pdat[!is.na(pdat$y), , drop = FALSE]
  pdat$n_label <- paste0(pdat$k, "(", pdat$studies, ")")
  pdat$letter_label <- ifelse(is.na(pdat$letters), "", pdat$letters)
  pdat$EffectClass <- ifelse(pdat$p_value < 0.05 & pdat$estimate > 0, "Positive",
                             ifelse(pdat$p_value < 0.05 & pdat$estimate < 0, "Negative", "Not significant"))
  pdat$Panel <- factor(pdat$Panel, levels = panel_labels[metrics])

  qdat <- list()
  for (metric in metrics) {
    for (mod in moderator_order) {
      td <- tests[tests$Domain == domain & tests$Metric == metric & tests$Moderator == mod, , drop = FALSE]
      if (nrow(td) == 0) next
      block_y <- rows$y[rows$Moderator == mod]
      if (length(block_y) == 0) next
      qdat[[length(qdat) + 1]] <- data.frame(
        Metric = metric,
        Panel = panel_labels[metric],
        Moderator = mod,
        y = max(block_y) - 0.35,
        label = paste0(fmt_qm(td$QM[1]), "\n", fmt_p(td$QMp[1])),
        stringsAsFactors = FALSE
      )
    }
  }
  qdat <- bind_rows(qdat)
  if (nrow(qdat) > 0) qdat$Panel <- factor(qdat$Panel, levels = panel_labels[metrics])

  list(plot = pdat, rows = rows, separators = layout$separators, q = qdat)
}

plot_nc_like <- function(domain, metrics, metric_labels, file_stub, output_dir = plot_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  built <- build_plot_data(domain, metrics, metric_labels)
  pdat <- built$plot
  rows <- built$rows
  sep_df <- built$separators
  qdat <- built$q
  if (nrow(pdat) == 0) return(invisible(NULL))

  panel_range <- pdat %>%
    group_by(Panel) %>%
    summarise(
      xmin = min(ci_low, estimate, na.rm = TRUE),
      xmax = max(ci_high, estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      span = pmax(xmax - xmin, 0.2),
      x_view_min = xmin - single_panel_left_label_pad * span,
      x_view_max = xmax + single_panel_right_label_pad * span,
      view_span = x_view_max - x_view_min,
      x_qm = x_view_min + single_panel_qm_inset * view_span,
      x_n = x_view_max - single_panel_n_inset * view_span,
      x_letter = xmax + single_panel_letter_offset * span
    )

  pdat <- left_join(
    pdat,
    panel_range[, c("Panel", "x_n", "x_letter", "x_view_min", "x_view_max")],
    by = "Panel"
  )
  pdat$ci_low_plot <- pmax(pdat$ci_low, pdat$x_view_min, na.rm = TRUE)
  pdat$ci_high_plot <- pmin(pdat$ci_high, pdat$x_view_max, na.rm = TRUE)
  if (nrow(qdat) > 0) qdat <- left_join(qdat, panel_range[, c("Panel", "x_qm")], by = "Panel")
  x_limits <- range(c(panel_range$x_view_min, panel_range$x_view_max), na.rm = TRUE)

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_hline(data = sep_df, aes(yintercept = y), linetype = "dashed",
               color = "grey55", linewidth = 0.45) +
    geom_errorbar(data = pdat, aes(y = y, xmin = ci_low_plot, xmax = ci_high_plot),
                  orientation = "y", width = 0.23, linewidth = 0.8, color = "black") +
    geom_point(data = pdat, aes(x = estimate, y = y, fill = EffectClass),
               shape = 21, size = 3.9, color = "black", stroke = 0.8) +
    geom_text(data = pdat[pdat$RowType != "Overall", , drop = FALSE],
              aes(x = x_letter, y = y, label = letter_label),
              hjust = 0.5, size = 4.3, color = "black") +
    geom_text(data = pdat, aes(x = x_n, y = y, label = n_label),
              hjust = 1, size = 3.3, color = "black") +
    geom_text(data = qdat, aes(x = x_qm, y = y, label = label),
              hjust = 0, vjust = 1, size = 3.0, color = "black") +
    facet_wrap(
      ~ Panel,
      ncol = length(metrics),
      scales = "free_x"
    ) +
    scale_y_continuous(
      breaks = rows$y,
      labels = rows$RowLabel,
      expand = expansion(add = c(0.6, 0.6))
    ) +
    scale_fill_manual(
      values = c("Positive" = "#4E9BD3", "Negative" = "#F2B84B", "Not significant" = "white")
    ) +
    scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "on") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 13, face = "bold", lineheight = 1.0),
      axis.text.y = element_text(color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 9),
      axis.ticks.y = element_line(color = "black"),
      panel.spacing.x = unit(0.7, "lines"),
      plot.margin = margin(8, single_panel_right_margin, 8, 8)
    )

  height <- max(7, min(18, 0.34 * nrow(rows) + 2.2))
  width <- if (length(metrics) == 1) single_panel_width else max(single_panel_width, 2.8 * length(metrics))
  pdf_path <- file.path(output_dir, paste0(file_stub, ".pdf"))
  png_path <- file.path(output_dir, paste0(file_stub, ".png"))
  ggsave(pdf_path, p, width = width, height = height)
  ggsave(png_path, p, width = width, height = height, dpi = 300)
  invisible(p)
}

plot_nc_like_single_metric <- function(domain, metric, metric_labels, output_dir = single_plot_dir) {
  metric_label <- metric_labels[metric]
  if (is.na(metric_label)) metric_label <- metric

  plot_nc_like(
    domain = domain,
    metrics = metric,
    metric_labels = setNames(metric_label, metric),
    file_stub = paste("categorical_NC_like", domain, metric, sep = "_"),
    output_dir = output_dir
  )
}

plot_nc_like_single_metrics <- function(domains, metrics, metric_labels, output_dir = single_plot_dir) {
  for (domain in domains) {
    for (metric in metrics) {
      plot_nc_like_single_metric(domain, metric, metric_labels, output_dir)
    }
  }
}

make_row_layout_multi_domain <- function(domains, metrics, moderators) {
  rows <- data.frame(
    RowType = "Overall",
    Moderator = "Overall",
    Level = "Overall",
    RowLabel = "Overall",
    stringsAsFactors = FALSE
  )
  separators_after <- 1

  for (mod in moderators) {
    mod_levels <- unique(subgroup$Level[
      subgroup$Domain %in% domains &
        subgroup$Metric %in% metrics &
        subgroup$Moderator == mod
    ])
    preferred <- level_order[[mod]]
    if (!is.null(preferred)) {
      mod_levels <- c(preferred[preferred %in% mod_levels], sort(setdiff(mod_levels, preferred)))
    } else {
      mod_levels <- sort(mod_levels)
    }
    if (length(mod_levels) == 0) next

    rows <- bind_rows(
      rows,
      data.frame(
        RowType = "Subgroup",
        Moderator = mod,
        Level = mod_levels,
        RowLabel = display_level(mod, mod_levels),
        stringsAsFactors = FALSE
      )
    )
    separators_after <- c(separators_after, nrow(rows))
  }

  rows$RowIndex <- seq_len(nrow(rows))
  rows$y <- nrow(rows) - rows$RowIndex + 1
  rows$RowKey <- paste(rows$Moderator, rows$Level, sep = "||")
  sep_df <- data.frame(y = nrow(rows) - separators_after + 0.5)
  list(rows = rows, separators = sep_df)
}

build_plot_data_multi_domain <- function(domains, metrics, metric_labels, letters_prefix = letters) {
  layout <- make_row_layout_multi_domain(domains, metrics, moderator_order)
  rows <- layout$rows

  panel_labels <- setNames(
    paste0(letters_prefix[seq_along(metrics)], "\n", metric_labels[metrics]),
    metrics
  )

  plot_rows <- list()
  for (domain in domains) {
    for (metric in metrics) {
      over <- overall[overall$Domain == domain & overall$Metric == metric, , drop = FALSE]
      if (nrow(over) > 0) {
        plot_rows[[length(plot_rows) + 1]] <- data.frame(
          DomainFacet = domain,
          Metric = metric,
          Panel = panel_labels[metric],
          Moderator = "Overall",
          Level = "Overall",
          estimate = over$estimate[1],
          ci_low = over$ci_low[1],
          ci_high = over$ci_high[1],
          p_value = over$p_value[1],
          k = over$k[1],
          studies = over$studies[1],
          letters = "",
          stringsAsFactors = FALSE
        )
      }

      sub <- subgroup[subgroup$Domain == domain & subgroup$Metric == metric &
                        subgroup$Moderator %in% moderator_order, , drop = FALSE]
      if (nrow(sub) > 0) {
        sub$DomainFacet <- domain
        sub$Panel <- panel_labels[metric]
        plot_rows[[length(plot_rows) + 1]] <- sub[, c(
          "DomainFacet", "Metric", "Panel", "Moderator", "Level", "estimate", "ci_low", "ci_high",
          "p_value", "k", "studies", "letters"
        )]
      }
    }
  }

  pdat <- bind_rows(plot_rows)
  pdat$RowKey <- paste(pdat$Moderator, pdat$Level, sep = "||")
  pdat <- left_join(pdat, rows[, c("RowKey", "y", "RowLabel", "RowType")], by = "RowKey")
  pdat <- pdat[!is.na(pdat$y), , drop = FALSE]
  pdat$n_label <- paste0(pdat$k, "(", pdat$studies, ")")
  pdat$letter_label <- ifelse(is.na(pdat$letters), "", pdat$letters)
  pdat$EffectClass <- ifelse(pdat$p_value < 0.05 & pdat$estimate > 0, "Positive",
                             ifelse(pdat$p_value < 0.05 & pdat$estimate < 0, "Negative", "Not significant"))
  pdat$DomainFacet <- factor(pdat$DomainFacet, levels = domains)
  pdat$Panel <- factor(pdat$Panel, levels = panel_labels[metrics])

  qdat <- list()
  for (domain in domains) {
    for (metric in metrics) {
      for (mod in moderator_order) {
        td <- tests[tests$Domain == domain & tests$Metric == metric & tests$Moderator == mod, , drop = FALSE]
        if (nrow(td) == 0) next
        block_y <- rows$y[rows$Moderator == mod]
        if (length(block_y) == 0) next
        qdat[[length(qdat) + 1]] <- data.frame(
          DomainFacet = domain,
          Metric = metric,
          Panel = panel_labels[metric],
          Moderator = mod,
          y = max(block_y) - 0.35,
          label = paste0(fmt_qm(td$QM[1]), "\n", fmt_p(td$QMp[1])),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  qdat <- bind_rows(qdat)
  if (nrow(qdat) > 0) {
    qdat$DomainFacet <- factor(qdat$DomainFacet, levels = domains)
    qdat$Panel <- factor(qdat$Panel, levels = panel_labels[metrics])
  }

  list(plot = pdat, rows = rows, separators = layout$separators, q = qdat)
}

plot_nc_like_multi_domain <- function(domains, metrics, metric_labels, file_stub) {
  built <- build_plot_data_multi_domain(domains, metrics, metric_labels)
  pdat <- built$plot
  rows <- built$rows
  sep_df <- built$separators
  qdat <- built$q
  if (nrow(pdat) == 0) return(invisible(NULL))

  panel_range <- pdat %>%
    group_by(Panel) %>%
    summarise(
      xmin = min(ci_low, estimate, na.rm = TRUE),
      xmax = max(ci_high, estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      span = pmax(xmax - xmin, 0.2),
      x_qm = xmin - 0.58 * span,
      x_n = xmax + 0.40 * span,
      x_letter = xmax + 0.12 * span
    )

  pdat <- left_join(pdat, panel_range[, c("Panel", "x_n", "x_letter")], by = "Panel")
  if (nrow(qdat) > 0) qdat <- left_join(qdat, panel_range[, c("Panel", "x_qm")], by = "Panel")

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_hline(data = sep_df, aes(yintercept = y), linetype = "dashed",
               color = "grey55", linewidth = 0.45) +
    geom_errorbar(data = pdat, aes(y = y, xmin = ci_low, xmax = ci_high),
                  orientation = "y", width = 0.23, linewidth = 0.8, color = "black") +
    geom_point(data = pdat, aes(x = estimate, y = y, fill = EffectClass),
               shape = 21, size = 3.6, color = "black", stroke = 0.8) +
    geom_text(data = pdat[pdat$RowType != "Overall", , drop = FALSE],
              aes(x = x_letter, y = y, label = letter_label),
              hjust = 0, size = 4.2, color = "black") +
    geom_text(data = pdat, aes(x = x_n, y = y, label = n_label),
              hjust = 0, size = 3.3, color = "black") +
    geom_text(data = qdat, aes(x = x_qm, y = y, label = label),
              hjust = 0, vjust = 1, size = 2.9, color = "black") +
    facet_grid(DomainFacet ~ Panel, scales = "free_x") +
    scale_y_continuous(
      breaks = rows$y,
      labels = rows$RowLabel,
      expand = expansion(add = c(0.6, 0.6))
    ) +
    scale_fill_manual(
      values = c("Positive" = "#4E9BD3", "Negative" = "#F2B84B", "Not significant" = "white")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.52, 0.42))) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 12, face = "bold", lineheight = 1.0),
      strip.text.y = element_text(size = 12, face = "bold", angle = 0),
      axis.text.y = element_text(color = "black", size = 9),
      axis.text.x = element_text(color = "black", size = 8),
      axis.ticks.y = element_line(color = "black"),
      panel.spacing.x = unit(0.7, "lines"),
      panel.spacing.y = unit(0.9, "lines"),
      plot.margin = margin(8, 65, 8, 8)
    )

  height <- max(11, min(30, 0.34 * nrow(rows) * length(domains) + 2.8))
  width <- max(10, 2 * length(metrics))
  pdf_path <- file.path(plot_dir, paste0(file_stub, ".pdf"))
  png_path <- file.path(plot_dir, paste0(file_stub, ".png"))
  ggsave(pdf_path, p, width = width, height = height)
  ggsave(png_path, p, width = width, height = height, dpi = 300)
  invisible(p)
}

plot_nc_like_single_metrics(
  domains = c("Bacterial", "Fungi"),
  metrics = alpha_metrics,
  metric_labels = alpha_metric_labels
)

plot_nc_like_single_metrics(
  domains = c("Bacterial", "Fungi"),
  metrics = beta_metrics,
  metric_labels = beta_metric_labels
)

plot_nc_like_single_metrics(
  domains = "Bacterial",
  metrics = phylum_metrics$Bacterial,
  metric_labels = setNames(phylum_metrics$Bacterial, phylum_metrics$Bacterial)
)

plot_nc_like_single_metrics(
  domains = "Fungi",
  metrics = phylum_metrics$Fungi,
  metric_labels = setNames(phylum_metrics$Fungi, phylum_metrics$Fungi)
)

cat("NC-like categorical forest plots complete.\n")
cat("Output directory:", normalizePath(single_plot_dir, winslash = "/"), "\n")

