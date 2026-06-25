# Continuous moderator analyses and random-forest screening
# Clean GitHub code package for the N deposition soil microbiome synthesis.
# This file combines the R scripts actually used for the manuscript-level analyses.

# Sources combined: 7
# Generated on: 2026-06-25


# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/??????.R
# ================================================================

# MISSING SOURCE FILE: Data_analysis/Remote_working_package/06_R_scripts/??????.R

# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/Meta????????.R
# ================================================================

# MISSING SOURCE FILE: Data_analysis/Remote_working_package/06_R_scripts/Meta????????.R

# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/Meta????????_alpha_beta_combined.R
# ================================================================

# MISSING SOURCE FILE: Data_analysis/Remote_working_package/06_R_scripts/Meta????????_alpha_beta_combined.R

# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/Meta????????_beta.R
# ================================================================

# MISSING SOURCE FILE: Data_analysis/Remote_working_package/06_R_scripts/Meta????????_beta.R

# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/?????????.R
# ================================================================

# MISSING SOURCE FILE: Data_analysis/Remote_working_package/06_R_scripts/?????????.R

# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/??????.R
# ================================================================

# MISSING SOURCE FILE: Data_analysis/Remote_working_package/06_R_scripts/??????.R

# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/plot_continuous_alpha_grid_vi_floor.R
# ================================================================

rm(list = ls())

library(dplyr)
library(ggplot2)
library(metafor)

setwd("E:/BaiduSyncdisk/N_deposition1/Data_analysis")

out_dir <- "outputs"
grid_plot_dir <- file.path(out_dir, "plots", "continuous_meta_regression", "grid_panels")
dir.create(grid_plot_dir, showWarnings = FALSE, recursive = TRUE)

effect_file <- file.path(out_dir, "meta_effects_all.csv")
effects <- read.csv(effect_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!("RR" %in% names(effects))) effects$RR <- effects$yi
if (!("Vi" %in% names(effects))) effects$Vi <- effects$vi
if (!("StudyID" %in% names(effects))) effects$StudyID <- effects$ID

target_moderator <- "N_dose"
target_domains <- c("Bacterial", "Fungi")
target_metrics <- c("richness", "shannon", "simpson")
x_axis_label <- "N addition rate (kg N ha-1 yr-1)"

min_rows_for_model <- 8
min_unique_x <- 4

# Some Simpson-effect variances are extremely tiny and make rma.mv fail.
# Change this value for sensitivity checks, e.g. 1e-8, 1e-9, 1e-10.
vi_floor <- 1e-10

plot_height <- 6.5
plot_width <- 10.0
file_stub <- paste("grid_continuous", target_moderator, "alpha_div_vi_floor", sep = "_")

Ftheme <- theme(
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = NA, fill = NA),
  panel.background = element_blank(),
  plot.title = element_text(size = 12),
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_text(size = 12),
  text = element_text(size = 12),
  strip.background = element_blank(),
  strip.placement = "in"
)

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

valid_effect_rows <- function(dat) {
  dat <- dat[
    !is.na(dat$RR) & !is.na(dat$Vi) & dat$Vi > 0 &
      is.finite(dat$RR) & is.finite(dat$Vi),
    ,
    drop = FALSE
  ]
  dat$Vi_original <- dat$Vi
  dat$Vi_was_floored <- dat$Vi < vi_floor
  dat$Vi <- pmax(dat$Vi, vi_floor)
  dat
}

standardize_with_params <- function(x) {
  center <- mean(x, na.rm = TRUE)
  scale_value <- stats::sd(x, na.rm = TRUE)
  if (is.na(scale_value) || scale_value == 0) scale_value <- 1
  list(value = (x - center) / scale_value, center = center, scale = scale_value)
}

prepare_x <- function(dat, moderator, transform) {
  x_raw <- to_num(dat[[moderator]])

  if (transform == "log") {
    if (any(is.na(x_raw)) || any(x_raw <= 0, na.rm = TRUE)) return(NULL)
    x_trans <- log(x_raw)
  } else {
    x_trans <- x_raw
  }

  z <- standardize_with_params(x_trans)
  dat$x_raw <- x_raw
  dat$x_model <- z$value
  dat$x_center <- z$center
  dat$x_scale <- z$scale
  dat$x_transform <- transform

  dat[
    !is.na(dat$x_raw) & !is.na(dat$x_model) &
      is.finite(dat$x_raw) & is.finite(dat$x_model),
    ,
    drop = FALSE
  ]
}

fit_one_continuous_model <- function(dat, moderator, transform) {
  dat <- valid_effect_rows(dat)
  dat <- prepare_x(dat, moderator, transform)

  if (is.null(dat) ||
      nrow(dat) < min_rows_for_model ||
      length(unique(dat$x_raw)) < min_unique_x) {
    return(NULL)
  }

  dat$StudyID <- as.character(dat$StudyID)

  res <- tryCatch(
    suppressWarnings(
      rma.mv(
        yi = RR,
        V = Vi,
        mods = ~ x_model,
        random = ~ 1 | StudyID,
        data = dat,
        method = "REML"
      )
    ),
    error = function(e) {
      message(
        unique(dat$Domain)[1], " / ", unique(dat$Metric)[1],
        " / ", transform, " failed: ", e$message
      )
      NULL
    }
  )

  if (is.null(res) || length(res$beta) < 2) return(NULL)
  list(model = res, data = dat)
}

summarize_continuous_model <- function(fit, moderator, transform) {
  res <- fit$model
  dat <- fit$data

  data.frame(
    Domain = unique(dat$Domain)[1],
    Metric = unique(dat$Metric)[1],
    Moderator = moderator,
    Transform = transform,
    k = nrow(dat),
    studies = length(unique(dat$StudyID)),
    floored_vi = sum(dat$Vi_was_floored, na.rm = TRUE),
    x_min = min(dat$x_raw, na.rm = TRUE),
    x_max = max(dat$x_raw, na.rm = TRUE),
    intercept = as.numeric(res$beta[1]),
    slope = as.numeric(res$beta[2]),
    slope_p = as.numeric(res$pval[2]),
    AIC = as.numeric(AIC(res)),
    stringsAsFactors = FALSE
  )
}

make_prediction <- function(fit, moderator, transform) {
  res <- fit$model
  dat <- fit$data

  x_seq <- seq(
    min(dat$x_raw, na.rm = TRUE),
    max(dat$x_raw, na.rm = TRUE),
    length.out = 200
  )
  x_trans <- if (transform == "log") log(x_seq) else x_seq
  x_model <- (x_trans - dat$x_center[1]) / dat$x_scale[1]
  pred <- predict(res, newmods = x_model)

  data.frame(
    Domain = unique(dat$Domain)[1],
    Metric = unique(dat$Metric)[1],
    x_raw = x_seq,
    fit = pred$pred,
    ci_low = pred$ci.lb,
    ci_high = pred$ci.ub,
    stringsAsFactors = FALSE
  )
}

target_data <- effects %>%
  filter(Domain %in% target_domains, Metric %in% target_metrics)

analysis_units <- split(
  target_data,
  list(target_data$Domain, target_data$Metric),
  drop = TRUE
)

best_model_summaries <- list()
all_predictions <- list()
best_raw_points <- list()

for (unit_name in names(analysis_units)) {
  unit_dat <- analysis_units[[unit_name]]
  fits <- list()

  for (transform in c("raw", "log")) {
    fit <- fit_one_continuous_model(unit_dat, target_moderator, transform)
    if (!is.null(fit)) fits[[transform]] <- fit
  }

  if (length(fits) == 0) next

  model_table <- bind_rows(lapply(names(fits), function(tr) {
    summarize_continuous_model(fits[[tr]], target_moderator, tr)
  }))

  best_transform <- model_table$Transform[which.min(model_table$AIC)]
  best_fit <- fits[[best_transform]]

  best_model_summaries[[length(best_model_summaries) + 1]] <-
    model_table %>% filter(Transform == best_transform)

  all_predictions[[length(all_predictions) + 1]] <-
    make_prediction(best_fit, target_moderator, best_transform)

  pts <- best_fit$data
  pts$weight <- 1 / pts$Vi
  pts$point_size <- (sqrt(pts$weight) / max(sqrt(pts$weight), na.rm = TRUE)) * 2.5 + 1.2
  best_raw_points[[length(best_raw_points) + 1]] <- pts
}

sdat <- bind_rows(best_model_summaries)
ldat <- bind_rows(all_predictions)
pdat <- bind_rows(best_raw_points)

if (nrow(pdat) > 0 && nrow(ldat) > 0 && nrow(sdat) > 0) {
  pdat$Domain <- factor(pdat$Domain, levels = target_domains)
  pdat$Metric <- factor(pdat$Metric, levels = target_metrics)

  ldat$Domain <- factor(ldat$Domain, levels = target_domains)
  ldat$Metric <- factor(ldat$Metric, levels = target_metrics)

  sdat$Domain <- factor(sdat$Domain, levels = target_domains)
  sdat$Metric <- factor(sdat$Metric, levels = target_metrics)

  sdat$label <- sprintf(
    "%s; p = %.3g\nk = %d, study = %d",
    sdat$Transform,
    sdat$slope_p,
    sdat$k,
    sdat$studies
  )

  label_dat <- sdat %>%
    select(Domain, Metric, label) %>%
    mutate(x_raw = -Inf, RR = Inf)

  ldat <- ldat %>%
    left_join(
      sdat %>% select(Domain, Metric, slope_p),
      by = c("Domain", "Metric")
    ) %>%
    mutate(LineType = ifelse(slope_p < 0.05, "solid", "dashed"))

  p <- ggplot() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.5
    ) +
    geom_point(
      data = pdat,
      aes(x = x_raw, y = RR, size = point_size, color = Domain),
      alpha = 0.35
    ) +
    geom_ribbon(
      data = ldat,
      aes(x = x_raw, ymin = ci_low, ymax = ci_high, fill = Domain),
      alpha = 0.15
    ) +
    geom_line(
      data = ldat,
      aes(x = x_raw, y = fit, color = Domain, linetype = LineType),
      linewidth = 1
    ) +
    geom_text(
      data = label_dat,
      aes(x = x_raw, y = RR, label = label),
      hjust = -0.05,
      vjust = 1.2,
      size = 3.2,
      color = "black",
      lineheight = 0.9
    ) +
    facet_grid(Domain ~ Metric, scales = "free_y") +
    scale_size_identity() +
    scale_linetype_identity() +
    scale_color_manual(values = c(Bacterial = "#4E79A7", Fungi = "#59A14F")) +
    scale_fill_manual(values = c(Bacterial = "#4E79A7", Fungi = "#59A14F")) +
    labs(
      x = x_axis_label,
      y = "Response ratio (RR)"
    ) +
    Ftheme +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )

  ggsave(
    file.path(grid_plot_dir, paste0(file_stub, ".pdf")),
    p,
    width = plot_width,
    height = plot_height
  )
  ggsave(
    file.path(grid_plot_dir, paste0(file_stub, ".png")),
    p,
    width = plot_width,
    height = plot_height,
    dpi = 300
  )

  print(sdat[, c("Domain", "Metric", "Transform", "k", "studies", "floored_vi", "slope_p", "AIC")])
  print(p)
  cat("Output directory:", normalizePath(grid_plot_dir, winslash = "/"), "\n")
} else {
  cat("Target data are insufficient for this continuous meta-regression grid.\n")
}

